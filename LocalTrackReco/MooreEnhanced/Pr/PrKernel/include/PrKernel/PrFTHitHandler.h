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
#include "FTDAQ/FTInfo.h"
#include "GaudiKernel/Range.h"
#include "Kernel/MultiIndexedContainer.h"
#include "Kernel/STLExtensions.h"
#include "PrKernel/PrHit.h"

/** @class PrFTHitHandler PrFTHitHandler.h
 *  Handler of Sci-Fit hits, the object is stored in the transient event store and
 *  each algorithm reads the object from there without calling the HitManagers (tools)
 *  @author Renato Quagliani
 *  @author Sebastien Ponce
 */

template <typename Hit>
class PrFTHitHandler final {
public:
  using HitContainer = LHCb::Container::MultiIndexedContainer<Hit, LHCb::Detector::FT::NFTZones>;

  /// Type for range of PrFT Hits within a container
  using HitRange   = typename HitContainer::HitRange;
  using HitIter    = typename HitContainer::Hits::const_iterator;
  using HitRevIter = typename HitContainer::Hits::const_reverse_iterator;

  PrFTHitHandler() = default;

  PrFTHitHandler( typename HitContainer::Hits&& hits, typename HitContainer::Offsets&& offsets )
      : m_hits{std::forward<typename HitContainer::Hits>( hits ),
               std::forward<typename HitContainer::Offsets>( offsets )} {}

  // Constructor with the capacity of the HitContainer
  PrFTHitHandler( int size, LHCb::Allocators::EventLocal<ModPrHit> alloc = {} ) : m_hits( size, alloc ){};

  template <typename I>
  void insert( unsigned int lay, I&& b, I&& e ) {
    m_hits.insert( std::forward<I>( b ), std::forward<I>( e ), lay );
  }

  template <typename... Args>
  void addHitInZone( unsigned int lay, Args&&... args ) {
    m_hits.addHit( std::forward_as_tuple( std::forward<Args>( args )... ), lay );
  }

  void setOffsets() { m_hits.setOffsets(); }

  /// Return the current number of zones
  auto isEmpty( unsigned int lay ) const { return m_hits.empty( lay ); }

  LHCb::span<const Hit> getRange_lowerBound( unsigned int lay, float xMin ) const {
    auto r = hits( lay );
    return LHCb::make_span(
        std::lower_bound( r.begin(), r.end(), xMin, []( const Hit& a, float testval ) { return a.x() < testval; } ),
        r.end() );
  }

  LHCb::span<const Hit> getRange( unsigned int lay, float xMin, float xMax ) const {
    auto r = hits( lay );
    auto first =
        std::lower_bound( r.begin(), r.end(), xMin, []( const Hit& a, float testval ) { return a.x() < testval; } );
    auto last = std::lower_bound( first, r.end(), xMax, []( const Hit& a, float testval ) { return a.x() < testval; } );
    return LHCb::make_span( first, last );
  }

  const HitContainer& hits() const { return m_hits; }

  HitRange hits( const unsigned int layer ) const { return m_hits.range( layer ); }

  Hit& hit( size_t index ) { return m_hits.hit( index ); }

  const Hit& hit( size_t index ) const { return m_hits.hit( index ); }

private:
  HitContainer m_hits;
};
