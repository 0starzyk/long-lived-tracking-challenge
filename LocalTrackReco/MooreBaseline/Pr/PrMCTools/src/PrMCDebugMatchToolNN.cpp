/*****************************************************************************\
* (c) Copyright 2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "PrKernel/PrDebugTrackingToolBase.h"

#include <string>
#include <vector>

#include "GaudiKernel/DataObjectHandle.h"

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/Track.h"
#include "LHCbMath/SIMDWrapper.h"
#include "Linker/LinkedTo.h"

namespace LHCb::Pr::MatchNN {
  struct PrMCDebugMatchToolNN : public PrDebugTrackingToolBase {

    // inherit standard constructors
    using PrDebugTrackingToolBase::PrDebugTrackingToolBase;

    int check( int veloIndex = -1, int seedIndex = -1, const std::vector<int>& = {} ) const override;

  private:
    DataObjectReadHandle<LHCb::Tracks>      m_veloTracks{this, "VeloTracks", ""};
    DataObjectReadHandle<LHCb::Tracks>      m_seedTracks{this, "SeedTracks", ""};
    DataObjectReadHandle<LHCb::LinksByKey>  m_veloTrackLinks{this, "VeloTrackLinks", ""};
    DataObjectReadHandle<LHCb::LinksByKey>  m_seedTrackLinks{this, "SeedTrackLinks", ""};
    DataObjectReadHandle<LHCb::MCParticles> m_mcparticles{this, "MCParticles", ""};
    DataObjectReadHandle<LHCb::MCProperty>  m_trackInfo{this, "TrackInfo", ""};
  };

  // Declaration of the Tool Factory
  DECLARE_COMPONENT_WITH_ID( PrMCDebugMatchToolNN, "PrMCDebugMatchToolNN" )

  int PrMCDebugMatchToolNN::check( int veloIndex, int seedIndex, const std::vector<int>& ) const {
    assert( veloIndex >= 0 && seedIndex >= 0 );
    std::vector<const LHCb::MCParticle*> velo_mcps{};
    std::vector<const LHCb::MCParticle*> seed_mcps{};

    m_veloTrackLinks.get()->applyToLinks( veloIndex, [&]( auto /*veloIndex*/, auto mcKey, auto /*weight*/ ) {
      velo_mcps.push_back( m_mcparticles.get()->operator()( mcKey ) );
    } );
    m_seedTrackLinks.get()->applyToLinks( seedIndex, [&]( auto /*seedIndex*/, auto mcKey, auto /*weight*/ ) {
      seed_mcps.push_back( m_mcparticles.get()->operator()( mcKey ) );
    } );
    auto       found{0};
    const auto trackInfo = MCTrackInfo{*m_trackInfo.get()};
    if ( !velo_mcps.empty() && !seed_mcps.empty() ) {
      auto veloIter = velo_mcps.begin();
      while ( veloIter != velo_mcps.end() ) {
        assert( *veloIter != nullptr );
        if ( !trackInfo.hasVeloAndT( *veloIter ) ) {
          veloIter = std::next( veloIter );
          continue;
        }
        auto seedIter = seed_mcps.begin();
        while ( seedIter != seed_mcps.end() ) {
          assert( *seedIter != nullptr );
          if ( *veloIter == *seedIter ) {
            if ( 11 == std::abs( ( *veloIter )->particleID().pid() ) ) {
              found = -1;
              break;
            } else {
              found = 1;
              break;
            }
          } else {
            seedIter = std::next( seedIter );
          }
        }
        if ( found )
          break;
        else {
          veloIter = std::next( veloIter );
        }
      }
    }
    return found;
  }

} // namespace LHCb::Pr::MatchNN
