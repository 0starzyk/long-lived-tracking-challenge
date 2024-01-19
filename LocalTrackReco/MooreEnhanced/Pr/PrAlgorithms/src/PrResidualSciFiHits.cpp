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
// Include files
#include "Event/PrHits.h"
#include "Event/PrLongTracks.h"
#include "FTDAQ/FTInfo.h"
#include "LHCbAlgs/Transformer.h"

#include <array>
#include <bitset>
#include <limits>

//-----------------------------------------------------------------------------
// class : PrResidualSciFiHits
// Store residual SciFiHits after other Algorithms, e.g. PrMatchNN or PrForwardTracking
// the input tracks and SciFiHits are in SOA structure
//
// 2020-04-02 : Peilian Li
// 2020-03-18 : Andre Guenther (make sentinels consistent with PrStoreSciFiHits)
//
//-----------------------------------------------------------------------------

namespace LHCb::Pr::FT {

  class ResidualHits
      : public Algorithm::Transformer<Hits( const EventContext&, const LHCb::Pr::Long::Tracks&, const Hits& )> {
    using Tracks = Long::Tracks;

  public:
    ResidualHits( const std::string& name, ISvcLocator* pSvcLocator );

    Hits operator()( const EventContext&, const Tracks&, const Hits& ) const override;
  };

  // Declaration of the Algorithm Factory
  DECLARE_COMPONENT_WITH_ID( ResidualHits, "PrResidualSciFiHits" )

  //=============================================================================
  // Standard constructor, initializes variables
  //=============================================================================
  ResidualHits::ResidualHits( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"TracksLocation", ""}, KeyValue{"SciFiHitsLocation", PrFTInfo::SciFiHitsLocation}},
                     KeyValue{"SciFiHitsOutput", PrFTInfo::SciFiHitsLocation} ) {}

  //=============================================================================
  // Main execution
  //=============================================================================
  Hits ResidualHits::operator()( const EventContext& evtCtx, const Tracks& tracks, const Hits& fthits ) const {

    Hits hits{LHCb::getMemResource( evtCtx ), fthits.size()};

    if ( tracks.size() == 0 ) {
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Track container '" << inputLocation<Tracks>() << "' is empty" << endmsg;
      // explicit copy required (to generally avoid copies)
      hits.copy_from( fthits );
      return hits;
    }

    std::bitset<Hits::maxNumberOfHits> used{};

    /// mark used SciFi Hits
    for ( const auto& track : tracks.scalar() ) {
      const int nfthits = track.nFTHits().cast();
      for ( int id = 0; id != nfthits; id++ ) {
        const auto idx = track.ft_index( id ).cast();
        used[idx]      = true;
      }
    }
    constexpr auto xu  = LHCb::Detector::FT::xZonesUpper;
    constexpr auto uvu = LHCb::Detector::FT::uvZonesUpper;

    constexpr auto xd       = LHCb::Detector::FT::xZonesLower;
    constexpr auto uvd      = LHCb::Detector::FT::uvZonesLower;
    constexpr auto hitzones = std::array<int, LHCb::Detector::FT::NFTZones>{
        xd[0], uvd[0], uvd[1], xd[1], xd[2], uvd[2], uvd[3], xd[3], xd[4], uvd[4], uvd[5], xd[5],
        xu[0], uvu[0], uvu[1], xu[1], xu[2], uvu[2], uvu[3], xu[3], xu[4], uvu[4], uvu[5], xu[5]};

    for ( auto zone : hitzones ) {
      const auto [zoneBegin, zoneEnd] = fthits.getZoneIndices( zone );
      hits.setZoneIndex( zone, hits.size() );
      for ( auto iHit{zoneBegin}; iHit < zoneEnd; ++iHit ) {
        if ( used[iHit] ) continue;
        hits.appendColumn( fthits.planeCode( iHit ), fthits.x( iHit ), fthits.hotHitInfo( iHit ),
                           fthits.coldHitInfo( iHit ) );
      }
      hits.appendColumn( std::numeric_limits<uint8_t>::max(), 1.e9f, {}, {} );
    }

    hits.setZoneIndex( LHCb::Detector::FT::NFTZones, hits.getZoneIndex( xu[0] ) );
    // so when asking for lastZone, lastZone+2 gives the very end of the container
    hits.setZoneIndex( LHCb::Detector::FT::NFTZones + 1, hits.size() );
    // avoid FPEs
    for ( unsigned i{0}; i < SIMDWrapper::best::types::size; ++i ) {
      hits.appendColumn( std::numeric_limits<uint8_t>::max(), 1.e9f, {}, {} );
    }

    return hits;
  }

} // namespace LHCb::Pr::FT
