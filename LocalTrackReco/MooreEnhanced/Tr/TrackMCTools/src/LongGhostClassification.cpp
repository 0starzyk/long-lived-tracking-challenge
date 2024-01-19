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

class LongGhostClassification : public TrackGhostClassificationBase {

public:
  /// constructor
  using TrackGhostClassificationBase::TrackGhostClassificationBase;

  /**
   *  Check this is a ghost .
   *  specialize for long tracks [check velo and T separately]
   *  @param aTrack to link
   *  @return bool true if a ghost
   */
  using TrackGhostClassificationBase::isGhost;

  bool isGhost( LHCbIDs::const_iterator& start, LHCbIDs::const_iterator& stop ) const override;

private:
  StatusCode specific( LHCbIDs::const_iterator& start, LHCbIDs::const_iterator& stop,
                       LHCb::GhostTrackInfo& tinfo ) const override;
};

DECLARE_COMPONENT( LongGhostClassification )

using namespace LHCb;

StatusCode LongGhostClassification::specific( LHCbIDs::const_iterator& start, LHCbIDs::const_iterator& stop,
                                              LHCb::GhostTrackInfo& tinfo ) const {

  // split into velo and T hits
  LHCbIDs::const_iterator iter = start;
  LHCbIDs                 tHits;
  tHits.reserve( 20 );
  LHCbIDs vHits;
  vHits.reserve( 20 );
  for ( ; iter != stop; ++iter ) {
    if ( iter->detectorType() == LHCbID::channelIDtype::FT ) {
      tHits.push_back( *iter );
    } else if ( iter->detectorType() == LHCbID::channelIDtype::VP ) {
      vHits.push_back( *iter );
    }
  } // for iter

  // match the T Hits
  LHCb::GhostTrackInfo::LinkPair tMatch = bestPair( tHits );

  // match the velo Hits
  LHCb::GhostTrackInfo::LinkPair vMatch = bestPair( vHits );

  if ( ( tMatch.first == 0 || tMatch.second < m_purityCut ) && ( vMatch.first == 0 && vMatch.second < m_purityCut ) ) {
    tinfo.setClassification( GhostTrackInfo::Classification::GhostParent );
  } else if ( ( tMatch.first == 0 || tMatch.second < m_purityCut ) )
    tinfo.setClassification( GhostTrackInfo::Classification::GhostTParent );
  else if ( ( vMatch.first == 0 || vMatch.second < m_purityCut ) )
    tinfo.setClassification( GhostTrackInfo::Classification::GhostVeloParent );

  if ( isMatched( vMatch ) && isMatched( tMatch ) && vMatch.first != tMatch.first ) {
    tinfo.setClassification( LHCb::GhostTrackInfo::Classification::InconsistentParts );
  }

  return StatusCode::SUCCESS;
}

bool LongGhostClassification::isGhost( TrackGhostClassificationBase::LHCbIDs::const_iterator& start,
                                       TrackGhostClassificationBase::LHCbIDs::const_iterator& stop ) const {

  LHCbIDs tHits;
  tHits.reserve( 20 );
  LHCbIDs vHits;
  vHits.reserve( 20 );

  for ( auto iter = start; iter != stop; ++iter ) {
    if ( iter->detectorType() == LHCbID::channelIDtype::FT ) {
      tHits.push_back( *iter );
    } else if ( iter->detectorType() == LHCbID::channelIDtype::VP ) {
      vHits.push_back( *iter );
    }
  } // for iter

  // match the T Hits
  LHCb::GhostTrackInfo::LinkPair tMatch = bestPair( tHits );
  if ( !isReal( tMatch ) ) return true;

  // match the velo Hits
  LHCb::GhostTrackInfo::LinkPair vMatch = bestPair( vHits );
  if ( !isReal( vMatch ) ) return true;

  return vMatch.first != tMatch.first;
}
