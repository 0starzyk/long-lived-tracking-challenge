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
#include "Detector/VP/VPChannelID.h"
#include "Event/Track.h"
#include "TrackGhostClassificationBase.h"

class VeloGhostClassification : public TrackGhostClassificationBase {

public:
  /// constructer
  using TrackGhostClassificationBase::TrackGhostClassificationBase;

private:
  StatusCode specific( LHCbIDs::const_iterator& start, LHCbIDs::const_iterator& stop,
                       LHCb::GhostTrackInfo& tinfo ) const override;
};

DECLARE_COMPONENT( VeloGhostClassification )

using namespace LHCb;

StatusCode VeloGhostClassification::specific( LHCbIDs::const_iterator& start, LHCbIDs::const_iterator& stop,
                                              LHCb::GhostTrackInfo& tinfo ) const {

  LHCbIDs vHits;
  vHits.reserve( 20 );
  for ( LHCbIDs::const_iterator iter = start; iter != stop; ++iter ) {
    if ( iter->detectorType() == LHCbID::channelIDtype::VP ) {
      Detector::VPChannelID vChan = iter->vpID();
      vHits.push_back( vChan );
    }
  } // for iter

  // match the Hits
  LHCb::GhostTrackInfo::LinkPair match = bestPair( vHits );

  if ( match.first == 0 ) { tinfo.setClassification( GhostTrackInfo::Classification::GhostParent ); }

  return StatusCode::SUCCESS;
}
