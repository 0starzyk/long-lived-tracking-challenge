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

#include <algorithm>
#include <string>
#include <vector>

#include "GaudiKernel/DataObjectHandle.h"

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/PrHits.h"
#include "Event/Track.h"
#include "Kernel/LHCbID.h"

namespace LHCb::Pr::Forward {
  struct PrMCDebugForwardTool : public PrDebugTrackingToolBase {

    // inherit standard constructors
    using PrDebugTrackingToolBase::PrDebugTrackingToolBase;

    int check( int track_index = -1, int index = -1, const std::vector<int>& = {} ) const override;

  private:
    Gaudi::Property<float> m_matchFrac{this, "MatchFraction", 0.7};

    DataObjectReadHandle<LHCb::Tracks>                m_inputTracks{this, "InputTracks", ""};
    DataObjectReadHandle<LHCb::LinksByKey>            m_inputTrackLinks{this, "InputTrackLinks", ""};
    DataObjectReadHandle<LHCb::LinksByKey>            m_SciFiHitLinks{this, "SciFiHitLinks", ""};
    DataObjectReadHandle<LHCb::MCParticles>           m_mcparticles{this, "MCParticles", ""};
    DataObjectReadHandle<LHCb::Pr::Hits<HitType::FT>> m_SciFiHits{this, "SciFiHits", ""};
    DataObjectReadHandle<LHCb::MCProperty>            m_trackInfo{this, "TrackInfo", ""};
  };

  // Declaration of the Tool Factory
  DECLARE_COMPONENT_WITH_ID( PrMCDebugForwardTool, "PrMCDebugForwardTool" )

  int PrMCDebugForwardTool::check( int track_index, int, const std::vector<int>& scifi_indices ) const {
    assert( track_index >= 0 );

    std::vector<const LHCb::MCParticle*> input_mcps{};
    m_inputTrackLinks.get()->applyToLinks( track_index, [&]( auto /*track_index*/, auto mcKey, auto /*weight*/ ) {
      input_mcps.push_back( m_mcparticles.get()->operator()( mcKey ) );
    } );

    if ( std::none_of( input_mcps.begin(), input_mcps.end(), [trackInfo = MCTrackInfo{*m_trackInfo.get()}]( auto mcp ) {
           return trackInfo.hasVeloAndT( mcp );
         } ) ) {
      return 0;
    }

    std::vector<LHCb::LHCbID> lhcbids;
    lhcbids.reserve( 12 );
    std::transform( scifi_indices.begin(), scifi_indices.end(), std::back_inserter( lhcbids ),
                    [hits = m_SciFiHits.get()]( auto idx ) { return hits->lhcbid( idx ); } );

    std::vector<int> match_counters( input_mcps.size(), 0 );
    for ( size_t i{0}; i < input_mcps.size(); ++i ) {
      for ( auto id : lhcbids ) {
        m_SciFiHitLinks.get()->applyToLinks( id.lhcbID(), [&]( auto /*FTChannelID*/, auto mcKey, auto /*weight*/ ) {
          const auto linked_mcp = m_mcparticles.get()->operator()( mcKey );
          if ( m_mcparticles.get() != linked_mcp->parent() ) {
            throw GaudiException( "SciFiHitLinks do not use the same underlying container as MCParticles!",
                                  this->name(), StatusCode::FAILURE );
          }
          if ( input_mcps[i] == linked_mcp ) ++match_counters[i];
        } );
      }
    }
    const auto best          = std::max_element( match_counters.begin(), match_counters.end() );
    const auto best_mcp      = input_mcps[std::distance( match_counters.begin(), best )];
    const auto matching_frac = static_cast<float>( *best ) / lhcbids.size();
    const auto found         = matching_frac >= m_matchFrac;
    if ( const auto pid = std::abs( best_mcp->particleID().pid() ); pid != 11 ) {
      return found;
    } else {
      return found ? pid : found;
    }
  }

} // namespace LHCb::Pr::Forward
