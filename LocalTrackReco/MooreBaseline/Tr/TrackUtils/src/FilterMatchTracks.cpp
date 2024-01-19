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
#include "GaudiAlg/GaudiAlgorithm.h"

//-----------------------------------------------------------------------------
// Implementation file for class : FilterMatchTracks
//
// 2010-06-14 : Olivier Callot
//-----------------------------------------------------------------------------

/** @class FilterMatchTracks FilterMatchTracks.h
 *  Fileter Match tracks identical to Forward tracks
 *
 *  @author Olivier Callot
 *  @date   2010-06-14
 */

class FilterMatchTracks : public GaudiAlgorithm {
public:
  /// Standard constructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode execute() override; ///< Algorithm execution

private:
  Gaudi::Property<bool> m_filter{this, "Filter", true};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( FilterMatchTracks )

//=============================================================================
// Main execution
//=============================================================================
StatusCode FilterMatchTracks::execute() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  LHCb::Tracks* forward = get<LHCb::Tracks>( LHCb::TrackLocation::Forward );
  LHCb::Tracks* match   = get<LHCb::Tracks>( LHCb::TrackLocation::Match );

  if ( msgLevel( MSG::DEBUG ) )
    debug() << "Comparing " << forward->size() << " forward to " << match->size() << " match tracks" << endmsg;

  for ( auto itForward = forward->begin(); forward->end() != itForward; ++itForward ) {
    for ( auto itMatch = match->begin(); match->end() != itMatch; ++itMatch ) {
      if ( ( *itForward )->nLHCbIDs() < ( *itMatch )->nLHCbIDs() ) continue;
      LHCb::Track* myMatch = *itMatch;
      unsigned int nCommon = ( *itForward )->nCommonLhcbIDs( *myMatch );
      if ( nCommon == ( *itMatch )->nLHCbIDs() ) {
        if ( msgLevel( MSG::DEBUG ) )
          debug() << "Tracks forward " << ( *itForward )->key() << " is containing match track " << ( *itMatch )->key()
                  << " .";
        if ( m_filter.value() ) {
          match->erase( itMatch );
          itMatch = match->begin() - 1; //== re-initialize the iterator
          if ( msgLevel( MSG::DEBUG ) ) debug() << " Match track removed.";
        }
        if ( msgLevel( MSG::DEBUG ) ) debug() << endmsg;
      }
    }
  }
  if ( msgLevel( MSG::DEBUG ) ) debug() << "After filter, rests " << match->size() << " match tracks to fit" << endmsg;

  return StatusCode::SUCCESS;
}
