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

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/UniqueIDGenerator.h"
#include "Functors/Filter.h"
#include "Functors/with_functors.h"
#include "Kernel/EventLocalUnique.h"
#include "LHCbAlgs/Transformer.h"
#include "SelKernel/TrackZips.h"
#include <vector>

namespace Pr {

  namespace {
    template <typename T>
    struct TypeHelper {};

    // Helper for compile time check if a type is a filtered-zip type
    template <typename>
    struct must_be_wrapped : std::false_type {};

    template <typename... T>
    struct must_be_wrapped<std::tuple<T...>> : std::true_type {};

    template <typename T>
    inline constexpr bool must_be_wrapped_v = must_be_wrapped<T>::value;

    // Given the tag type T, TypeHelper tells us both the return type of the
    // functor (OutputType) and whether this return type is to be converted into
    // a non-owning view (WrapOutput), which means that we can deduce the output
    // type of the convert() function too, just using the tag type T. This is the
    // return type of operator() in the algorithm. We also need to deduce that
    // type with the leading bool stripped off, so that we can write down the
    // name of the functional framework base class that the algorithm inherits
    // from
    template <typename TagType>
    struct OutputHelper {
      using FilteredType               = Functors::filtered_t<TagType>;
      static constexpr bool WrapOutput = must_be_wrapped_v<FilteredType>;

      template <typename KeyValue>
      static auto names() {
        if constexpr ( WrapOutput ) {
          return std::tuple{KeyValue{"Output", ""}, KeyValue{"OutputStorage", ""}};
        } else {
          return KeyValue{"Output", ""};
        }
      }

      static auto convert( EventContext const& evtCtx, FilteredType&& filtered ) {
        if constexpr ( WrapOutput ) {
          // OutputType is some owning type that we don't want to use directly, but
          // which we need to keep alive so that we can use a non-owning view into
          // it. Achieve this by wrapping it up in a unique_ptr -- so the address
          // of the contained object is stable -- and returning that. We use this
          // special make_event_local_unique function to create a unique_ptr that is
          // allocated using the custom, event-local memory pool and deleted with a
          // custom, stateful, deleter.
          auto storage =
              LHCb::make_event_local_unique<FilteredType>( LHCb::getMemResource( evtCtx ), std::move( filtered ) );
          // Make a non-owning view, which will store the address of the object
          // hidden inside the unique_ptr.
          auto view = LHCb::Event::make_zip( std::as_const( *storage.get() ) );
          // The point of the output wrapping is to make sure that the input
          // and output types are the same. Check that's the case.
          static_assert( std::is_same_v<decltype( view ), TagType> );
          return std::tuple{false, std::move( view ), std::move( storage )};
        } else {
          // Similarly if the output wrapping was not deemed necessary this
          // "should" mean that the filtered type already matched the input type.
          // Actually with the v2 zip machinery then `T` maps onto `std::tuple<T>`
          if constexpr ( std::is_same_v<FilteredType, std::tuple<TagType>> ) {
            return std::tuple{false, std::get<0>( std::move( filtered ) )};
          } else {
            return std::tuple{false, std::move( filtered )};
          }
        }
      }

      // this is std::tuple<bool, A[, B]>
      using AlgorithmOutput = std::invoke_result_t<decltype( convert ), EventContext const&, FilteredType>;

      // helper to strip off the bool from the type
      template <typename T>
      struct remove_first_type {};

      template <typename T, typename... Ts>
      struct remove_first_type<std::tuple<T, Ts...>> {
        using type = std::tuple<Ts...>;
      };

      // this is std::tuple<A[, B]>
      using DataTuple = typename remove_first_type<AlgorithmOutput>::type;
    };

    // Just shorthand for below
    template <typename T>
    using FilterTransform = LHCb::Algorithm::MultiTransformerFilter<
        typename OutputHelper<T>::DataTuple( EventContext const&, T const& ),
        Gaudi::Functional::Traits::BaseClass_t<LHCb::DetDesc::AlgorithmWithCondition<>>>;

    template <typename T>
    struct FilterCut {
      constexpr static auto PropertyName = "Cut";
      using Signature                    = Functors::filtered_t<T>( T const& );
    };

    template <typename T>
    struct SOAFilterCut {
      constexpr static auto PropertyName = "Cut";
      using InputType                    = typename TypeHelper<T>::InputType;
      using Signature                    = Functors::filtered_t<InputType>( InputType );
    };
  } // namespace

  /** @class Filter PrFilter.cpp
   *
   *  Filter<T> applies a selection to an input Selection<T> and returns a new Selection<T> object.
   *
   *  @tparam T The selected object type (e.g. Track, Particle, ...). By contruction this is not copied, as the
   *            input/output type Selection<T> is just a view of some other underlying storage.
   */
  template <typename T>
  class Filter final : public with_functors<FilterTransform<T>, FilterCut<T>> {
  public:
    using Base    = with_functors<FilterTransform<T>, FilterCut<T>>;
    using OHelper = OutputHelper<T>;

    Filter( const std::string& name, ISvcLocator* pSvcLocator )
        : Base( name, pSvcLocator, {"Input", ""}, OHelper::template names<typename Base::KeyValue>() ) {}

    // Return type is std::tuple<bool, A[, B]>
    typename OHelper::AlgorithmOutput operator()( EventContext const& evtCtx, T const& in ) const override {
      // Get the functor from the with_functors mixin
      auto const& pred = this->template getFunctor<FilterCut<T>>();

      // Apply the functor to get something we can return
      auto filtered = pred( in );
      // Add an extra conversion step if requested
      // We do this first so the constexpr logic in `convert` can normalise
      // things before we update the statistics and pass/fail flag. For
      // example, one zip implementation maps T -> T and the other maps
      // T -> std::tuple<T> when filtering.
      auto converted = OHelper::convert( evtCtx, std::move( filtered ) );
      // Get the pass/fail flag the control flow sees
      bool& filter_pass = std::get<0>( converted );
      // Get the output view, whose type should match `in`
      auto const& filtered_view = std::get<1>( converted );
      // For use in the control flow: did we select anything?
      filter_pass = !filtered_view.empty();
      // Update the statistics.
      m_cutEff += {filtered_view.size(), in.size()};
      return converted;
    }

  private:
    // Counter for recording cut retention statistics
    mutable Gaudi::Accumulators::BinomialCounter<> m_cutEff{this, "Cut selection efficiency"};
  };
} // namespace Pr
