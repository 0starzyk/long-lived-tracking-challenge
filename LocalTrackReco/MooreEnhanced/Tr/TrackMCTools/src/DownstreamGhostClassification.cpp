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
#include "Event/Track.h"
#include "TrackGhostClassificationBase.h"

class DownstreamGhostClassification : public TrackGhostClassificationBase {

public:
  /// constructer
  using TrackGhostClassificationBase::TrackGhostClassificationBase;

private:
  StatusCode specific( LHCbIDs::const_iterator& start, LHCbIDs::const_iterator& stop,
                       LHCb::GhostTrackInfo& tinfo ) const override;
};

DECLARE_COMPONENT( DownstreamGhostClassification )

using namespace LHCb;

StatusCode DownstreamGhostClassification::specific( LHCbIDs::const_iterator& start, LHCbIDs::const_iterator& stop,
                                                    LHCb::GhostTrackInfo& tinfo ) const {

  // split into velo and T hits
  LHCbIDs::const_iterator iter = start;
  LHCbIDs                 utHits;
  utHits.reserve( 20 );
  LHCbIDs tHits;
  tHits.reserve( 20 );
  for ( ; iter != stop; ++iter ) {
    if ( iter->detectorType() == LHCbID::channelIDtype::UT ) {
      utHits.push_back( *iter );
    } else if ( iter->detectorType() == LHCbID::channelIDtype::FT ) {
      tHits.push_back( *iter );
    }
  } // for iter

  // match the T Hits
  LHCb::GhostTrackInfo::LinkPair tMatch = bestPair( tHits );

  // match the UT Hits
  LHCb::GhostTrackInfo::LinkPair utMatch = bestPair( utHits );

  if ( tMatch.first == 0 || tMatch.second < m_purityCut ) {
    tinfo.setClassification( GhostTrackInfo::Classification::GhostParent );
  }

  if ( isMatched( tMatch ) && isMatched( utMatch ) && tMatch.first != utMatch.first ) {
    tinfo.setClassification( LHCb::GhostTrackInfo::Classification::InconsistentParts );
  }

  return StatusCode::SUCCESS;
}
