/*****************************************************************************\
* (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once
#include "Event/ZipUtils.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/Kernel.h"
#include "Kernel/STLExtensions.h"
#include <boost/type_index.hpp>
#include <numeric>
#include <type_traits>
#include <vector>

namespace Pr {
  namespace details {
    // Tag class to allow optimisation when we know the selection is a dummy.
    struct alwaysTrue {};
    struct alwaysFalse {};

    template <typename T>
    inline auto const typename_v = boost::typeindex::type_id_with_cvr<T>().pretty_name();

    template <typename S>
    void check_for_set_operation( S const& s1, S const& s2, const char* msg_prefix ) {
      if ( s1.m_container.data() != s2.m_container.data() ) {
        throw GaudiException{std::string{msg_prefix} +
                                 "(Selection, Selection) was given views to different underlying storage.",
                             typename_v<S>, StatusCode::SUCCESS};
      }
      assert( std::is_sorted( s1.m_indices.begin(), s1.m_indices.end() ) );
      assert( std::is_sorted( s2.m_indices.begin(), s2.m_indices.end() ) );
    }
  } // namespace details

  /** @class Selection PrSelection.h
   *
   *  Selection<T, IndexSize> represents a selection of some other, contiguous
   *  storage container (typically std::vector<T>, but this is not required).
   *
   *  @tparam T         The selected object type (e.g. Track, Particle, ...). By contruction this is not copied, so
   *                    Selection<T> is particularly useful when T is a large object.
   *  @tparam IndexSize The type used to represent indices into the underlying contiguous storage. Defaults to
   *                    uint16_t, meaning by default you can only select 65536 objects...
   */
  template <typename T, typename IndexSize = uint16_t>
  struct Selection {
    // TODO in optimised build we could reduce sizeof(Selection) by 40% by using a single pointer in container_t and
    //      using an index_vector type that imposes size()==capacity()
    using value_type   = T; // make it easier to write generic code that can handle both containers and Selections
    using container_t  = LHCb::span<const T>;
    using index_vector = typename std::vector<IndexSize>;

    // Custom iterator class for looping through the index vector but dereferencing to values in the actual container
    struct const_iterator {
      using index_iter_type   = typename index_vector::const_iterator;
      using value_type        = T;
      using pointer           = T const*;
      using reference         = T const&;
      using difference_type   = typename index_iter_type::difference_type;
      using iterator_category = std::random_access_iterator_tag;

      container_t     m_container;
      index_iter_type m_iter;

      T const& operator*() const { return m_container[*m_iter]; }

      T const* operator->() const { return &m_container[*m_iter]; }

      const_iterator& operator+=( difference_type n ) {
        m_iter += n;
        return *this;
      }

      friend bool operator==( const_iterator const& lhs, const_iterator const& rhs ) {
        return lhs.m_container.data() == rhs.m_container.data() && lhs.m_iter == rhs.m_iter;
      }

      friend bool operator!=( const_iterator const& lhs, const_iterator const& rhs ) { return !( lhs == rhs ); }

      friend difference_type operator-( const_iterator const& lhs, const_iterator const& rhs ) {
        return lhs.m_iter - rhs.m_iter;
      }

      const_iterator& operator++() {
        ++m_iter;
        return *this;
      }

      const_iterator& operator--() {
        --m_iter;
        return *this;
      }
    };

    Zipping::ZipFamilyNumber m_family{Zipping::generateZipIdentifier()};
    index_vector             m_indices;
    container_t              m_container;

    /** Constructor creating a Selection from a contiguous storage container.
     *
     *  By default all elements of the input span are marked as selected.
     *
     *  @param container       Contiguous storage container that this selection refers to.
     *  @param predicate       Optional predicate applied to elements of the input container.
     *  @param reserveCapacity Optional estimate of the number of elements that will be selected by predicate, which is
     *                         used to initialise the vector of selected indices.
     */
    template <typename Predicate = details::alwaysTrue>
    Selection( container_t container, Predicate&& predicate = {}, int reserveCapacity = -1 ) : m_container{container} {
      if ( !m_container.empty() && ( m_container.size() - 1 > std::numeric_limits<IndexSize>::max() ) ) {
        throw GaudiException{"Index overflow: " + std::to_string( container.size() - 1 ) + " > " +
                                 std::to_string( std::numeric_limits<IndexSize>::max() ),
                             details::typename_v<Selection>, StatusCode::SUCCESS};
      }

      if constexpr ( std::is_same_v<Predicate, details::alwaysTrue> ) {
        m_indices.resize( m_container.size() );
        // container_t is contiguous so we just need {0 .. container.size()-1}
        std::iota( m_indices.begin(), m_indices.end(), 0 );
      } else if constexpr ( std::is_same_v<Predicate, details::alwaysFalse> ) {
        if ( reserveCapacity >= 0 ) { m_indices.reserve( reserveCapacity ); }
      } else {
        m_indices.reserve( reserveCapacity < 0 ? m_container.size() : reserveCapacity );
        auto offset = 0;
        for ( auto const& i : m_container ) {
          if ( std::invoke( predicate, i ) ) { m_indices.push_back( offset ); }
          ++offset;
        }
      }
    }

    /** Constructor creating a Selection from another Selection, optionally applying an extra selection.
     *
     *  By default the same elements of the underlying storage are selected.
     *
     *  @param old_selection   Other Selection object to copy.
     *  @param predicate       Optional predicate to further filter the input selection.
     *  @param reserveCapacity Optional estimate of the number of elements that will be selected by predicate, which is
     *                         used to initialise the vector of selected indices.
     */
    template <typename Predicate>
    Selection( Selection const& old_selection, Predicate&& predicate, int reserveCapacity = -1 )
        : m_container{old_selection.m_container} {
      // It's imposed that old_selection has the same IndexSize as us, so we don't need to do any overflow checking
      // apply an additional selection on the input selection
      m_indices.reserve( reserveCapacity < 0 ? old_selection.size() : reserveCapacity );
      std::copy_if( old_selection.m_indices.begin(), old_selection.m_indices.end(), std::back_inserter( m_indices ),
                    [&]( auto i ) { return std::invoke( predicate, m_container[i] ); } );
    }

    const_iterator begin() const { return {m_container, m_indices.begin()}; }

    const_iterator end() const { return {m_container, m_indices.end()}; }

    std::size_t size() const { return m_indices.size(); }

    bool empty() const { return m_indices.empty(); }

    Zipping::ZipFamilyNumber zipIdentifier() const { return m_family; }

    T const& operator[]( std::size_t i ) const {
      assert( i < size() );
      return m_container[m_indices[i]];
    }

    friend bool operator==( Selection const& lhs, Selection const& rhs ) {
      return lhs.m_container.data() == rhs.m_container.data() && lhs.m_indices == rhs.m_indices;
    }

    friend bool operator!=( Selection const& lhs, Selection const& rhs ) { return !( lhs == rhs ); }

    template <typename Predicate>
    Selection select( Predicate&& p, int reserveCapacity = -1 ) const {
      return Selection{*this, std::forward<Predicate>( p ), reserveCapacity};
    }

    // Set operations -- these throw if two Selections that point to different underlying storage are provided.
    friend Selection set_union( Selection const& s1, Selection const& s2 ) {
      details::check_for_set_operation( s1, s2, "set_union" ); // check s1 and s2 are valid and compatible

      // Shortcuts so we can use .back() below
      if ( s1.empty() ) return s2;
      if ( s2.empty() ) return s1;

      // Create a new Selection object with an empty index vector with an appropriate reserved capacity
      std::size_t est_csize = std::max( s1.m_indices.back(), s2.m_indices.back() ) + 1;
      Selection   ret{s1.m_container, details::alwaysFalse{},
                    static_cast<int>( std::min( est_csize, s1.size() + s2.size() ) )};

      // Use std::set_union to combine s1.m_indices and s2.m_indices into a new index container
      std::set_union( s1.m_indices.begin(), s1.m_indices.end(), s2.m_indices.begin(), s2.m_indices.end(),
                      std::back_inserter( ret.m_indices ) );
      return ret;
    }

    friend Selection set_intersection( Selection const& s1, Selection const& s2 ) {
      details::check_for_set_operation( s1, s2, "set_intersection" ); // check s1 and s2 are valid and compatible

      // Create a new Selection object with an empty index vector with an appropriate reserved capacity
      Selection ret{s1.m_container, details::alwaysFalse{}, static_cast<int>( std::min( s1.size(), s2.size() ) )};

      // Use std::set_intersection to combine s1.m_indices and s2.m_indices into a new index container
      std::set_intersection( s1.m_indices.begin(), s1.m_indices.end(), s2.m_indices.begin(), s2.m_indices.end(),
                             std::back_inserter( ret.m_indices ) );
      return ret;
    }

    friend Selection set_difference( Selection const& s1, Selection const& s2 ) {
      details::check_for_set_operation( s1, s2, "set_difference" ); // check s1 and s2 are valid and compatible

      // Create a new Selection object with an empty index vector with an appropriate reserved capacity
      Selection ret{s1.m_container, details::alwaysFalse{}, static_cast<int>( s1.size() )};

      // Use std::set_difference to combine s1.m_indices and s2.m_indices into a new index container
      std::set_difference( s1.m_indices.begin(), s1.m_indices.end(), s2.m_indices.begin(), s2.m_indices.end(),
                           std::back_inserter( ret.m_indices ) );
      return ret;
    }

    friend bool includes( Selection const& s1, Selection const& s2 ) {
      details::check_for_set_operation( s1, s2, "includes" ); // check s1 and s2 are valid and compatible
      return std::includes( s1.m_indices.begin(), s1.m_indices.end(), s2.m_indices.begin(), s2.m_indices.end() );
    }

    friend Selection set_symmetric_difference( Selection const& s1, Selection const& s2 ) {
      details::check_for_set_operation( s1, s2,
                                        "set_symmetric_difference" ); // check s1 and s2 are valid and compatible

      // Shortcuts so we can use .back() below safely
      if ( s1.empty() ) { return s2; }
      if ( s2.empty() ) { return s1; }

      // Create a new Selection object with an empty index vector with an appropriate reserved capacity
      std::size_t est_csize = std::max( s1.m_indices.back(), s2.m_indices.back() ) + 1;
      Selection   ret{s1.m_container, details::alwaysFalse{},
                    static_cast<int>( std::min( est_csize, s1.size() + s2.size() ) )};

      // Use std::set_difference to combine s1.m_indices and s2.m_indices into a new index container
      std::set_symmetric_difference( s1.m_indices.begin(), s1.m_indices.end(), s2.m_indices.begin(), s2.m_indices.end(),
                                     std::back_inserter( ret.m_indices ) );
      return ret;
    }
  };

  // Template parameter deduction guides
  template <typename T, typename... Args>
  Selection( T, Args... )->Selection<typename T::value_type>;

  namespace details {
    // Helpers for requireSelection<S> below
    template <typename T>
    struct isSelectionHelper : std::false_type {};

    template <typename C, typename I>
    struct isSelectionHelper<Selection<C, I>> : std::true_type {};

    template <typename S>
    inline constexpr bool isSelection_v = details::isSelectionHelper<S>::value;
  } // namespace details

  // Helper to determine whether a type is a Selection
  template <typename S>
  using requireSelection = std::enable_if_t<details::isSelection_v<S>>;

  // Value-by-value comparison, don't care whether lhs and rhs are actually views onto the same container.
  template <typename S, typename = requireSelection<S>>
  bool equal_values( S const& lhs, S const& rhs ) {
    return std::equal( lhs.begin(), lhs.end(), rhs.begin(), rhs.end() );
  }

  template <typename S, typename = requireSelection<S>>
  bool equal_values( S const& lhs, typename S::container_t rhs ) {
    return std::equal( lhs.begin(), lhs.end(), rhs.begin(), rhs.end() );
  }

  template <typename S, typename = requireSelection<S>>
  bool equal_values( typename S::container_t lhs, S const& rhs ) {
    return std::equal( lhs.begin(), lhs.end(), rhs.begin(), rhs.end() );
  }
} // namespace Pr
