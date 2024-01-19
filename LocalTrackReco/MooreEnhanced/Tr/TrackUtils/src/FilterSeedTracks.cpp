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
// Implementation file for class : FilterSeedTracks
//
// 2010-06-15 : Olivier Callot
//-----------------------------------------------------------------------------

/** @class FilterSeedTracks FilterSeedTracks.h
 *  Filter the remaining Seed tracks.
 *
 *  @author Olivier Callot
 *  @date   2010-06-17
 */
class FilterSeedTracks : public GaudiAlgorithm {
public:
  /// Standard constructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode execute() override; ///< Algorithm execution

private:
  Gaudi::Property<bool> m_filter{this, "Filter", true};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( FilterSeedTracks )

//=============================================================================
// Main execution
//=============================================================================
StatusCode FilterSeedTracks::execute() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  LHCb::Tracks* forward = get<LHCb::Tracks>( LHCb::TrackLocation::Forward );
  LHCb::Tracks* seed    = get<LHCb::Tracks>( LHCb::TrackLocation::Seed );

  if ( msgLevel( MSG::DEBUG ) )
    debug() << "Comparing " << forward->size() << " forward to " << seed->size() << " seed tracks" << endmsg;

  for ( auto itForward = forward->begin(); forward->end() != itForward; ++itForward ) {
    for ( auto itSeed = seed->begin(); seed->end() != itSeed; ++itSeed ) {
      unsigned int nCommon = ( *itForward )->nCommonLhcbIDs( **itSeed );
      if ( ( *itSeed )->nLHCbIDs() == nCommon ) {
        if ( msgLevel( MSG::DEBUG ) )
          debug() << "Seed track " << ( *itSeed )->key() << " shares all its " << nCommon << " hits with forward track "
                  << ( *itForward )->key();
        if ( m_filter.value() ) {
          seed->erase( itSeed );
          itSeed = seed->begin() - 1; //== re-initialize the iterator
          if ( msgLevel( MSG::DEBUG ) ) debug() << " Seed track removed.";
        }
        if ( msgLevel( MSG::DEBUG ) ) debug() << endmsg;
      }
    }
  }

  LHCb::Tracks* match = get<LHCb::Tracks>( LHCb::TrackLocation::Match );
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "Comparing " << match->size() << " match to " << seed->size() << " seed tracks" << endmsg;

  for ( auto itMatch = match->begin(); match->end() != itMatch; ++itMatch ) {
    for ( auto itSeed = seed->begin(); seed->end() != itSeed; ++itSeed ) {
      unsigned int nCommon = ( *itMatch )->nCommonLhcbIDs( **itSeed );
      if ( ( *itSeed )->nLHCbIDs() == nCommon ) {
        if ( msgLevel( MSG::DEBUG ) )
          debug() << "Seed track " << ( *itSeed )->key() << " shares all its " << nCommon << " hits with match track "
                  << ( *itMatch )->key();
        if ( m_filter.value() ) {
          seed->erase( itSeed );
          itSeed = seed->begin() - 1; //== re-initialize the iterator
          if ( msgLevel( MSG::DEBUG ) ) debug() << " Seed track removed.";
        }
        if ( msgLevel( MSG::DEBUG ) ) debug() << endmsg;
      }
    }
  }

  LHCb::Tracks* downstream = get<LHCb::Tracks>( LHCb::TrackLocation::Downstream );
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "Comparing " << downstream->size() << " downstream to " << seed->size() << " seed tracks" << endmsg;

  for ( auto itDownstream = downstream->begin(); downstream->end() != itDownstream; ++itDownstream ) {
    for ( auto itSeed = seed->begin(); seed->end() != itSeed; ++itSeed ) {
      unsigned int nCommon = ( *itDownstream )->nCommonLhcbIDs( **itSeed );
      if ( ( *itSeed )->nLHCbIDs() == nCommon ) {
        if ( msgLevel( MSG::DEBUG ) )
          debug() << "Seed track " << ( *itSeed )->key() << " shares all its " << nCommon
                  << " hits with downstream track " << ( *itDownstream )->key();
        if ( m_filter.value() ) {
          seed->erase( itSeed );
          itSeed = seed->begin() - 1; //== re-initialize the iterator
          if ( msgLevel( MSG::DEBUG ) ) debug() << " Seed track removed.";
        }
        if ( msgLevel( MSG::DEBUG ) ) debug() << endmsg;
      }
    }
  }

  if ( msgLevel( MSG::DEBUG ) ) debug() << "After filter, rests " << seed->size() << " seed tracks to fit" << endmsg;
  return StatusCode::SUCCESS;
}
