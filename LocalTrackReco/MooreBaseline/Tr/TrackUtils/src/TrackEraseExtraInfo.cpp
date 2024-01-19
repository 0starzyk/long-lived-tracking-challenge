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

/** @class TrackEraseExtraInfo TrackEraseExtraInfo.h
 *  Algorithm that removes extra info from tracks such that it is not written to DST
 *
 *  Parameters:
 * - InputLocation: Input location for tracks.
 * - ErasableInfo: List of extra infos to erase.
 * - PrintExtraInfo: Print the extra info on the track before erasure.
 *
 *  @author S. Hansmann-Menzemer
 *  @date   20.07.2009
 *
 *  @author Michel De Cian
 *  @date   2015-06-27
 */

#include "Event/Track.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include <string>
#include <vector>

using namespace LHCb;

class TrackEraseExtraInfo final : public GaudiAlgorithm {

public:
  // Constructors and destructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode execute() override;

private:
  DataObjectReadHandle<LHCb::Tracks> m_inputLocation{this, "InputLocation", LHCb::TrackLocation::Default};
  Gaudi::Property<std::vector<int>>  m_erasableInfo{this,
                                                   "ErasableInfo",
                                                   {static_cast<int>( LHCb::Track::AdditionalInfo::PatQuality ),
                                                    static_cast<int>( LHCb::Track::AdditionalInfo::Cand1stQPat ),
                                                    static_cast<int>( LHCb::Track::AdditionalInfo::Cand2ndQPat ),
                                                    static_cast<int>( LHCb::Track::AdditionalInfo::Cand1stChi2Mat ),
                                                    static_cast<int>( LHCb::Track::AdditionalInfo::Cand2ndChi2Mat ),
                                                    static_cast<int>( LHCb::Track::AdditionalInfo::MatchChi2 ),
                                                    static_cast<int>( LHCb::Track::AdditionalInfo::TsaLikelihood ),
                                                    static_cast<int>( LHCb::Track::AdditionalInfo::nPRVeloRZExpect ),
                                                    static_cast<int>( LHCb::Track::AdditionalInfo::nPRVelo3DExpect )}};
  Gaudi::Property<bool>              m_printExtraInfo{this, "PrintExtraInfo", false};
};

DECLARE_COMPONENT( TrackEraseExtraInfo )
//=============================================================================
// Main execution
//=============================================================================
StatusCode TrackEraseExtraInfo::execute() {

  Tracks* inCont = m_inputLocation.getIfExists();
  if ( !inCont ) {
    return Warning( "Input container " + m_inputLocation.objKey() + " does not exist", StatusCode::SUCCESS, 20 );
  }

  // -- Print extra info which is on track
  if ( m_printExtraInfo.value() ) {
    for ( LHCb::Track* track : *inCont ) {
      info() << "ExtraInfo for track: " << track->type() << " : " << track->key() << endmsg;
      const LHCb::Track::ExtraInfo extraInfo = track->extraInfo();
      for ( const auto& ei : extraInfo ) {
        const auto addInfo = static_cast<LHCb::Track::AdditionalInfo>( ei.first );
        info() << " " << addInfo << "=" << ei.second << endmsg;
      }
    }
  }

  for ( auto& track : *inCont ) {
    for ( auto i : m_erasableInfo ) track->eraseInfo( static_cast<LHCb::Track::AdditionalInfo>( i ) );
  }

  return StatusCode::SUCCESS;
}
