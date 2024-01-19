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
#include "Linker/LinkedFrom.h"
#include <string>
#include <vector>

using namespace LHCb;
/** @class TrackCloneCleaner TrackCloneCleaner.h
 *
 *  Clean out clone tracks, using information from the Clone linker table
 *
 *  @author M.Needham
 *  @date   30/05/2006
 */
class TrackCloneCleaner final : public GaudiAlgorithm {

public:
  using GaudiAlgorithm::GaudiAlgorithm;
  StatusCode execute() override;

private:
  Gaudi::Property<std::vector<std::string>> m_inputLocations{
      this, "inputLocations", {TrackLocation::Default}}; ///< Locations of Tracks in TES
  Gaudi::Property<std::string> m_linkerLocation{this, "linkerLocation",
                                                TrackLocation::Default + "Clones"}; ///< Location of Clone linker in TES
  Gaudi::Property<double>      m_cloneCut{this, "CloneCut", 5000};                  ///< Clone cut value
};
namespace {

  /** @class WorkingTrack TrackCloneCleaner.h
   *
   *  Working track object for TrackCloneCleaner algorithm
   *
   *  @author M.Needham
   *  @date   30/05/2006
   */
  struct WorkingTrack final {
    WorkingTrack() = default;
    WorkingTrack( LHCb::Track* _track, const bool _clone = false ) : track( _track ), clone( _clone ) {}
    // Access track Chi^2
    [[nodiscard]] double chi2() const { return track->chi2PerDoF(); }
    /// Access number of LHCbIDs
    [[nodiscard]] double nLHCbIDs() const { return track->nLHCbIDs(); }
    /// Compare to see if its the same track
    bool sameTrack( const LHCb::Track* aTrack ) const { return track == aTrack; }
    /// Return the track type ranking
    [[nodiscard]] int trackTypeRank() const;
    LHCb::Track*      track = nullptr; ///< Pointer to the track object
    bool              clone = false;   ///< Clone flag
    using Vector            = std::vector<WorkingTrack>;
  };

  int WorkingTrack::trackTypeRank() const {
    // CRJ : Should probably make type 'ranking' configurable via options ?
    switch ( track->type() ) {
    case LHCb::Track::Types::Long:
      return 0;
    case LHCb::Track::Types::Downstream:
      return 1;
    case LHCb::Track::Types::Upstream:
      return 2;
    case LHCb::Track::Types::Ttrack:
      return 3;
    case LHCb::Track::Types::Velo:
      return 4;
    case LHCb::Track::Types::VeloBackward:
      return 5;
    case LHCb::Track::Types::Muon:
      return 6;
    default:
      return 999;
    }
  }

} // namespace
DECLARE_COMPONENT( TrackCloneCleaner )

StatusCode TrackCloneCleaner::execute() {

  // Get the clone linker info
  auto links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( m_linkerLocation ) );

  // copy the tracks into a temporary vector
  WorkingTrack::Vector tempTracks;

  // loop and make working tracks
  for ( const auto& loc : m_inputLocations ) {
    // tracks to flag
    auto trackCont = getIfExists<Tracks>( loc );
    if ( !trackCont ) {
      error() << "Container " << loc << " does not exist" << endmsg;
      continue;
    }
    if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "Found " << trackCont->size() << " Tracks at " << loc << endmsg; }
    tempTracks.reserve( tempTracks.size() + trackCont->size() );
    for ( const auto& t : *trackCont ) {
      // only consider tracks with clone info
      if ( !LinkedFrom<LHCb::Track>{links}.range( t ).empty() ) tempTracks.emplace_back( t );
    }
  }

  // sort by type Lowest rank , then highest # of LHCbID, then smallest chi2
  auto order = []( const WorkingTrack& lhs, const WorkingTrack& rhs ) {
    return std::tuple{lhs.trackTypeRank(), rhs.nLHCbIDs(), lhs.chi2()} <
           std::tuple{rhs.trackTypeRank(), lhs.nLHCbIDs(), rhs.chi2()};
  };
  std::sort( tempTracks.begin(), tempTracks.end(), order );

  for ( const auto& track : tempTracks ) {
    if ( msgLevel( MSG::VERBOSE ) ) {
      verbose() << "Trying track key=" << track.track->key() << " " << track.track->history()
                << " chi2=" << track.chi2() << " nMeas=" << track.nLHCbIDs() << endmsg;
    }

    // skips if already tagged as a rejected clone
    if ( track.clone ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Already flagged as a clone. Skipping" << endmsg;
      continue;
    }

    // pick up the clones
    for ( auto const& [cloneTrack, dist] : LinkedFrom<LHCb::Track>{links}.weightedRange( track.track ) ) {
      // double check track is not flag as clone of itself !!
      if ( &cloneTrack == track.track ) {
        Error( "Track flagged as clone of itself !!" ).ignore();
        continue;
      }

      if ( msgLevel( MSG::VERBOSE ) ) {
        verbose() << " -> Clone track key=" << cloneTrack.key() << " " << cloneTrack.history() << " dist=" << dist
                  << endmsg;
      }
      // check clone cut
      if ( dist < m_cloneCut ) {
        auto iter = std::find_if( tempTracks.begin(), tempTracks.end(),
                                  [ptr = &cloneTrack]( const WorkingTrack& t ) { return t.sameTrack( ptr ); } );
        if ( iter == tempTracks.end() ) continue;
        iter->clone = true;
        if ( iter->track->info( LHCb::Track::AdditionalInfo::CloneDist, 1e99 ) > dist ) {
          if ( msgLevel( MSG::VERBOSE ) ) {
            verbose() << "  -> Flagging track " << iter->track << " key=" << iter->track->key() << " "
                      << iter->track->history() << " as a clone" << endmsg;
          }
          iter->track->addInfo( LHCb::Track::AdditionalInfo::CloneDist, dist );
        }
      } // passed cut
    }   // clone track
  }     // tempTracks

  return StatusCode::SUCCESS;
}
