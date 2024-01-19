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
#include <array>
#include <string>
#include <vector>

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "GaudiKernel/DataObjectHandle.h"
#include "Kernel/LHCbID.h"

namespace LHCb::Pr {
  struct PrMCDebugReconstructibleLong : public PrDebugTrackingToolBase {

    // inherit standard constructors
    using PrDebugTrackingToolBase::PrDebugTrackingToolBase;

    int check( int = -1, int = -1, const std::vector<int>& = {} ) const override;

  private:
    DataObjectReadHandle<LHCb::LinksByKey>  m_hitLinks{this, "HitLinks", ""};
    DataObjectReadHandle<LHCb::MCParticles> m_mcparticles{this, "MCParticles", ""};
  };

  // Declaration of the Tool Factory
  DECLARE_COMPONENT_WITH_ID( PrMCDebugReconstructibleLong, "PrMCDebugReconstructibleLong" )

  int PrMCDebugReconstructibleLong::check( int mcp_index, int, const std::vector<int>& ) const {

    std::vector<LHCb::LHCbID> lhcbids{};
    lhcbids.reserve( 50 );
    m_hitLinks.get()->applyToAllLinks( [mcp_index, mcparticles = m_mcparticles.get(),
                                        &lhcbids]( auto raw_lhcbid, auto mcp_index_key, auto /*weight*/ ) {
      if ( mcp_index == static_cast<int>( mcp_index_key ) ) { lhcbids.emplace_back( raw_lhcbid ); }
    } );
    if ( lhcbids.empty() ) return 0;

    size_t nVeloHits{0};
    auto   stations = std::array{std::pair<bool, bool>{}, std::pair<bool, bool>{}, std::pair<bool, bool>{}};
    for ( auto id : lhcbids ) {
      if ( id.isVP() ) {
        ++nVeloHits;
      } else if ( id.isFT() ) {
        const auto ftid          = id.ftID();
        const auto station       = ftid.globalStationIdx();
        stations[station].first  = stations[station].first || ftid.isX();
        stations[station].second = stations[station].second || !ftid.isX();
      }
    }

    if ( nVeloHits >= 3 &&
         std::all_of( stations.begin(), stations.end(), []( auto pair ) { return pair.first && pair.second; } ) ) {
      return 1;
    }
    return 0;
  }
} // namespace LHCb::Pr
