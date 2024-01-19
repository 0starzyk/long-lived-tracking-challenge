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
// Implementation file for class : FilterDownstreamTracks
//
// 2010-06-15 : Olivier Callot
//-----------------------------------------------------------------------------

/** @class FilterDownstreamTracks FilterDownstreamTracks.h
 *  Filter Downstream tracks that share the T station part with a Forward
 *
 *  @author Olivier Callot
 *  @date   2010-06-15
 */
class FilterDownstreamTracks : public GaudiAlgorithm {
public:
  /// Standard constructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode execute() override; ///< Algorithm execution

private:
  Gaudi::Property<bool>              m_filter{this, "Filter", true};
  DataObjectReadHandle<LHCb::Tracks> m_forward{LHCb::TrackLocation::Forward, this};
  DataObjectReadHandle<LHCb::Tracks> m_downstream{LHCb::TrackLocation::Downstream, this};
  DataObjectReadHandle<LHCb::Tracks> m_match{LHCb::TrackLocation::Match, this};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( FilterDownstreamTracks )

//=============================================================================
// Main execution
//=============================================================================
StatusCode FilterDownstreamTracks::execute() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  LHCb::Tracks* forward    = m_forward.get();
  LHCb::Tracks* downstream = m_downstream.get();

  if ( msgLevel( MSG::DEBUG ) )
    debug() << "Comparing " << forward->size() << " forward to " << downstream->size() << " downstream tracks"
            << endmsg;

  for ( auto itForward = forward->begin(); forward->end() != itForward; ++itForward ) {
    for ( auto itDown = downstream->begin(); downstream->end() != itDown; ++itDown ) {
      unsigned int nCommon = ( *itForward )->nCommonLhcbIDs( **itDown );
      if ( ( *itDown )->nLHCbIDs() * 0.8 <= nCommon ) {
        if ( msgLevel( MSG::DEBUG ) )
          debug() << "Tracks forward " << ( *itForward )->key() << " shares hits with downstream track "
                  << ( *itDown )->key() << " : nCommon " << nCommon;
        if ( m_filter.value() ) {
          downstream->erase( itDown );
          itDown = downstream->begin() - 1; //== re-initialize the iterator
          if ( msgLevel( MSG::DEBUG ) ) debug() << " Downstream track removed.";
        }
        if ( msgLevel( MSG::DEBUG ) ) debug() << endmsg;
      }
    }
  }

  LHCb::Tracks* match = m_match.get();
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "Comparing " << match->size() << " match to " << downstream->size() << " downstream tracks" << endmsg;

  for ( auto itMatch = match->begin(); match->end() != itMatch; ++itMatch ) {
    for ( auto itDown = downstream->begin(); downstream->end() != itDown; ++itDown ) {
      unsigned int nCommon = ( *itMatch )->nCommonLhcbIDs( **itDown );
      if ( ( *itDown )->nLHCbIDs() * 0.8 <= nCommon ) {
        if ( msgLevel( MSG::DEBUG ) )
          debug() << "Tracks match " << ( *itMatch )->key() << " shares hits with downstream track "
                  << ( *itDown )->key() << " : nCommon " << nCommon;
        if ( m_filter.value() ) {
          downstream->erase( itDown );
          itDown = downstream->begin() - 1; //== re-initialize the iterator
          if ( msgLevel( MSG::DEBUG ) ) debug() << " Downstream track removed.";
        }
        if ( msgLevel( MSG::DEBUG ) ) debug() << endmsg;
      }
    }
  }
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "After filter, rests " << downstream->size() << " downstream tracks to fit" << endmsg;
  return StatusCode::SUCCESS;
}
