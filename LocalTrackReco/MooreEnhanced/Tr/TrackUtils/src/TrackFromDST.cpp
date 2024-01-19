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
#include <algorithm>

//-----------------------------------------------------------------------------
// Implementation file for class : TrackFromDST
//
// 2006-09-18 : Eduardo Rodrigues
//-----------------------------------------------------------------------------

/** @class TrackFromDST TrackFromDST.h
 *
 *  Algorithm to classify the tracks given as input according to
 *  their History / pattern recognition algorithms.
 *  Typically, this algorithm takes the tracks from the "best" container
 *  and remakes the original containers that were originally put together
 *  in the best container at the end of the tracking sequence in Brunel.
 *
 *  @author Eduardo Rodrigues
 *  @date   2006-09-18
 */
class TrackFromDST : public GaudiAlgorithm {

public:
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode execute() override; ///< Algorithm execution

private:
  // job options
  // -----------
  // input Track container path
  Gaudi::Property<std::string> m_tracksInContainer{this, "TracksInContainer", LHCb::TrackLocation::Default};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TrackFromDST )

namespace {

  std::pair<std::string, std::vector<LHCb::Track::History>> make( std::string_view                            s,
                                                                  std::initializer_list<LHCb::Track::History> l ) {
    return std::pair{std::string{s}, std::vector<LHCb::Track::History>{l}};
  }

  // Map between track history flags and output containers
  // ------------------------------------------------------------------
  static const auto s_map = std::array{make( LHCb::TrackLocation::Seed, {LHCb::Track::History::PrSeeding} ),
                                       make( LHCb::TrackLocation::Velo, {LHCb::Track::History::PrPixel} ),
                                       make( LHCb::TrackLocation::VeloTT, {LHCb::Track::History::PrVeloUT} ),
                                       make( LHCb::TrackLocation::Forward, {LHCb::Track::History::PrForward} ),
                                       make( LHCb::TrackLocation::Match, {LHCb::Track::History::PrMatch} ),
                                       make( LHCb::TrackLocation::Downstream, {LHCb::Track::History::PrDownstream} )};
} // namespace

//=============================================================================
// Main execution
//=============================================================================
StatusCode TrackFromDST::execute() {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  // Get all the tracks from in input container
  // ------------------------------------------
  const LHCb::Tracks* inTracks = get<LHCb::Tracks>( m_tracksInContainer );
  if ( msgLevel( MSG::DEBUG ) )
    debug() << "# of Tracks in input container \"" << m_tracksInContainer << "\" = " << inTracks->size() << endmsg;
  std::vector<LHCb::Track*> work{inTracks->begin(), inTracks->end()}; // make a copy so we partition them...
  auto                      first = work.begin(), last = work.end();
  for ( const auto& item : s_map ) {
    // partition the tracks according to their History flag
    auto pivot = std::stable_partition( first, last, [&]( const LHCb::Track* t ) {
      auto i = std::find( item.second.begin(), item.second.end(), t->history() );
      return i != item.second.end();
    } );
    // clone them to their corresponding containers
    auto out = std::make_unique<LHCb::Tracks>();
    std::for_each( first, pivot, [&]( const LHCb::Track* t ) { out->add( new LHCb::Track( *t, t->key() ) ); } );
    // put container on the TES
    put( out.release(), item.first );
    if ( msgLevel( MSG::DEBUG ) )
      debug() << "Stored " << std::distance( first, pivot ) << " tracks in " << item.first << endmsg;

    // and go to the next range of tracks...
    first = pivot;
  }
  if ( msgLevel( MSG::DEBUG ) ) {
    std::for_each( first, last,
                   [&]( const LHCb::Track* t ) { debug() << "Invalid track type " << t->history() << endmsg; } );
  }
  return StatusCode::SUCCESS;
}
