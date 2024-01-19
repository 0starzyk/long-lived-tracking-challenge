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
#include "GaudiKernel/ToolHandle.h"
#include "TrackInterfaces/ITrackCloneFinder.h"
#include "remove_if_partitioned.h"
#include <LHCbMath/BloomFilter.h>
#include <algorithm>

//-----------------------------------------------------------------------------
// Implementation file for class : TrackEventCloneKiller
//
// 2006-03-01 : Eduardo Rodrigues
// Based on the clone killer algorithm of Rutger van der Eijk (2002-06-17)
// 2008-04-15 : Adrian Perieanu
// Update for speed and clone rate
// 2014-12-08 : Manuel Schiller
// update for speed, other performance figures largely unchanged
//-----------------------------------------------------------------------------

/** @class TrackEventCloneKiller TrackEventCloneKiller.h
 *
 *
 *  @author Eduardo Rodrigues
 *  @date   2006-03-01
 *  Based on the clone killer algorithm of Rutger van der Eijk (2002-06-17)
 *
 *  @author Adrian Perieanu
 *  @date   2008-05-05
 *  Update for speed and clone rate
 *
 *  @author Manuel Schiller
 *  @date   2014-12-08
 *  further improve speed, use BloomFilter to recognise non-overlapping tracks
 */
class TrackEventCloneKiller final : public GaudiAlgorithm {
public:
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution

private:
  ToolHandle<ITrackCloneFinder> m_cloneFinder{this, "CloneFinderTool",
                                              "TrackCloneFinder"}; ///< interface to clone finder tool

  // Retrieve the input tracks from specified containers
  void getInputTracks( std::vector<LHCb::Track*>& allTracks );

  // job options
  // -----------
  // input Track container paths
  Gaudi::Property<std::vector<std::string>> m_tracksInContainers{
      this,
      "TracksInContainers",
      {LHCb::TrackLocation::Forward, LHCb::TrackLocation::Match, LHCb::TrackLocation::VeloTT,
       LHCb::TrackLocation::Downstream, LHCb::TrackLocation::Tsa}};
  // output Track container path
  Gaudi::Property<std::string> m_tracksOutContainer{this, "TracksOutContainer", LHCb::TrackLocation::Default};
  // list of track types not to be stored
  Gaudi::Property<std::vector<LHCb::Track::Types>> m_ignoredTrackTypes{
      this, "IgnoredTrackTypes", {}, [=]( auto& ) {
        // sort and sanitise the list of ignored track types
        std::sort( begin( m_ignoredTrackTypes ), end( m_ignoredTrackTypes ) );
        m_ignoredTrackTypes.value().erase( std::unique( begin( m_ignoredTrackTypes ), end( m_ignoredTrackTypes ) ),
                                           end( m_ignoredTrackTypes ) );
      }};
  // do not do a clone compare for tracks from a same container
  Gaudi::Property<bool> m_skipSameContainerTracks{this, "SkipSameContainerTracks", true};
  Gaudi::Property<bool> m_compareInSameContainerForwardUpstream{this, "CompareInSameContainerForwardUpstream", false};

  // In some cases we just want to flag the tracks but not copy them to an output
  // container. The following algorithms can supress clones by checking
  // LHCb::Track::Flags::Clone
  // default is old setting
  Gaudi::Property<bool> m_copyTracks{this, "CopyTracks", true};
};

DECLARE_COMPONENT( TrackEventCloneKiller )

//=============================================================================
// Initialization
//=============================================================================
StatusCode TrackEventCloneKiller::initialize() {
  return GaudiAlgorithm::initialize().andThen( [&] {
    // Print out the user-defined settings
    // -----------------------------------
    info() << endmsg << "============ TrackEventCloneKiller Settings ===========" << endmsg
           << "TracksInContainers                    : " << m_tracksInContainers.value() << endmsg
           << "TrackOutContainer                     : " << m_tracksOutContainer << endmsg
           << "IgnoredTrackTypes                     : " << getProperty( "IgnoredTrackTypes" ).toString() << endmsg
           << "CloneFinderTool                       : " << m_cloneFinder.type() << endmsg
           << "SkipSameContainerTracks               : " << getProperty( "SkipSameContainerTracks" ).toString()
           << endmsg << "CompareInSameContainerForwardUpstream : "
           << getProperty( "CompareInSameContainerForwardUpstream" ).toString() << endmsg
           << "=======================================================" << endmsg << endmsg;
  } );
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode TrackEventCloneKiller::execute() {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  // Put all the input tracks into a temporary vector of pointers
  std::vector<LHCb::Track*> allTracks;
  getInputTracks( allTracks );

  const unsigned nTotal = allTracks.size();
  // weights according to clone prior (less hits higher clone proba)
  // i.e. first long tracks, then upstream/downstream, then Velo and Ttracks
  constexpr auto t2w = [] {
    std::array<int, LHCb::Event::Enum::Track::members_of<LHCb::Event::Enum::Track::Type>().size()> tmp{};
    for ( auto& w : tmp ) w = 2;
    tmp[static_cast<int>( LHCb::Track::Types::Unknown )]    = 0;
    tmp[static_cast<int>( LHCb::Track::Types::Long )]       = 4;
    tmp[static_cast<int>( LHCb::Track::Types::Upstream )]   = 3;
    tmp[static_cast<int>( LHCb::Track::Types::Downstream )] = 3;
    return tmp;
  }();

  // sort such that longer tracks come first (more likely to be non-clones)
  std::stable_sort( std::begin( allTracks ), std::end( allTracks ),
                    [&]( const LHCb::Track* tr1, const LHCb::Track* tr2 ) {
                      const auto w1 = t2w[static_cast<int>( tr1->type() )];
                      const auto w2 = t2w[static_cast<int>( tr2->type() )];
                      if ( w1 > w2 )
                        return true;
                      else if ( w2 > w1 )
                        return false;
                      // in case of a tie based on the track type ordering above, the track
                      // with more LHCbIDs goes first
                      else
                        return tr1->nLHCbIDs() > tr2->nLHCbIDs();
                    } );

  /* strategy:
   *
   * we build up a list of likely non-clone tracks, track by track
   *
   * steps:
   * 1. for each unassessed track, loop over the the list of likely non-clone
   *    tracks, and check if the unassessed track is a clone of the likely
   *    non-clone tracks already present
   *
   * 2. if the unassessed track turns out to be a clone of a track on the
   *    likely non-clone list, discard it
   *
   * 3. else save the track to the likely non-clone list, and proceed to the
   *    next unassessed track
   *
   * the list of likely non-clone tracks is built in place in the allTracks
   * vector by using remove_if_partitioned
   *
   * since the actual check if two tracks are clones is somewhat involved, the
   * isClone predicate below is somewhat involved, trying to avoid the actual
   * test for "cloneness" by various tricks before delegating the hard work to
   * the m_cloneFinder tool
   *
   * for each track put on the likely non-clone list, we save a BloomFilter
   * fingerprinting its hit content, so we can quickly check for
   * non-overlapping tracks
   *
   * the task is further complicated by the fact that we have no strict track
   * quality ordering here, and the m_cloneFinder tool may flag either track as
   * a clone of the other (sorting the tracks as above helps quite a bit in the
   * sense that it's much easier to predict which track is flagged as clone,
   * actually); therefore a final pass over the allTrack array to remove the
   * "out-of-order" clones is required; these should be rare as long as long
   * tracks are fed to the algorithm first, and T (and Velo) tracks are last
   */

  // build a BloomFilter for the hit content of each likely non-clone track
  typedef BloomFilter<LHCb::LHCbID, 64, 154028, 1 << 20> BF;
  std::vector<BF>                                        bloomfilters; // 32 bytes per track, p=14.7%
  bloomfilters.reserve( allTracks.size() );

  // check if a track is a clone of one of the tracks on the likely non-clone
  // list (which is in std:begin(allTracks)...doneend)
  const auto isClone = [&allTracks, &bloomfilters, this]( LHCb::Track*                        tr1,
                                                          decltype( std::begin( allTracks ) ) doneend ) {
    const auto isForward  = []( const LHCb::Track* tr ) { return LHCb::Track::History::PrForward == tr->history(); };
    const auto isUpstream = []( const LHCb::Track* tr ) { return LHCb::Track::Types::Upstream == tr->type(); };

    const auto  parent1( m_skipSameContainerTracks.value() ? tr1->parent() : nullptr );
    const bool  isForwardOrUpstream1 = ( isForward( tr1 ) || isUpstream( tr1 ) );
    const auto& ids1( tr1->lhcbIDs() );
    const BF    bf1( std::begin( ids1 ), std::end( ids1 ) );
    auto        jt = std::begin( bloomfilters );
    for ( auto it = std::begin( allTracks ); doneend != it; ++it, ++jt ) {
      LHCb::Track* tr2( *it );
      // ignore known clones - don't attempt to find clones of clones, unless
      // they're also clones of good tracks (this can only happen for the
      // "out-of-order clones" mentioned above, any they ought to be unlikely,
      // since the list is sorted by decreasing number of LHCbIDs)
      if ( tr2->checkFlag( LHCb::Track::Flags::Clone ) ) continue;
      // skip comparison of tracks from the same container, if so desired,
      // subject to restrictions...
      if ( parent1 == tr2->parent() ) {
        // if m_compareInSameContainerForwardUpstream is true, do compare
        // tracks in the Forward and Upstream containers nevertheless
        const bool isForwardOrUpstream2 = ( isForward( tr2 ) || isUpstream( tr2 ) );
        if ( !m_compareInSameContainerForwardUpstream.value() ||
             !( isForwardOrUpstream1 && isForwardOrUpstream2 && tr1->type() == tr2->type() ) )
          continue;
      }
      // skip tracks that do not overlap in terms of hits - they cannot be
      // clones
      const BF& bf2( *jt );
      if ( ( bf1 & bf2 ).empty() ) continue;
      // flag if tr1 and tr2 are clones of each other
      m_cloneFinder->flagClones( *tr2, *tr1 );
      // if tr1 has just been flagged as a clone, we can stop
      if ( tr1->checkFlag( LHCb::Track::Flags::Clone ) ) return true;
    }
    // tr1 is a good candidate to be kept - save its bloomfilter for later use
    bloomfilters.push_back( bf1 );
    return false;
  };
  // build up the list of likely non-clone tracks, and remove definite clones
  allTracks.erase( remove_if_partitioned( std::begin( allTracks ), std::end( allTracks ), isClone ),
                   std::end( allTracks ) );

  // copy non clone tracks to an output container, if wished
  if ( m_copyTracks.value() ) {
    // remove any remaining "out-of-order clones"
    allTracks.erase(
        std::remove_if( std::begin( allTracks ), std::end( allTracks ),
                        []( const LHCb::Track* tr ) { return tr->checkFlag( LHCb::Track::Flags::Clone ); } ),
        std::end( allTracks ) );
    // count unique (non-clone) tracks
    const unsigned nUnique = allTracks.size();
    // remove ignored track types if so desired
    if ( !m_ignoredTrackTypes.empty() ) {
      allTracks.erase( std::remove_if( std::begin( allTracks ), std::end( allTracks ),
                                       [this]( const LHCb::Track* tr ) {
                                         return std::binary_search( begin( m_ignoredTrackTypes ),
                                                                    end( m_ignoredTrackTypes ), tr->type() );
                                       } ),
                       std::end( allTracks ) );
    }
    const unsigned nFiltered = nUnique - allTracks.size();

    LHCb::Tracks* tracksOutCont = new LHCb::Tracks();
    tracksOutCont->reserve( allTracks.size() );
    put( tracksOutCont, m_tracksOutContainer );
    std::for_each( std::begin( allTracks ), std::end( allTracks ),
                   [tracksOutCont]( const LHCb::Track* tr ) { tracksOutCont->add( new LHCb::Track( *tr ) ); } );

    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << "Stored " << allTracks.size() << " tracks, identified " << ( nTotal - nUnique ) << " clones of which "
              << nFiltered << " were not stored." << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
// Retrieve the input tracks from all the user-specified containers
// Note: are only taken into account Valid and Fitted tracks!
//=============================================================================
void TrackEventCloneKiller::getInputTracks( std::vector<LHCb::Track*>& allTracks ) {
  { // figure out the capacity that allTracks needs to provide
    size_t maxsz = 0;
    for ( const auto& inCont : m_tracksInContainers ) {
      auto inTracks = get<LHCb::Track::Range>( inCont );
      maxsz += inTracks.size();
    }
    allTracks.reserve( maxsz );
  }
  for ( const auto& inCont : m_tracksInContainers ) {
    auto tracks = get<LHCb::Track::Range>( inCont );
    if ( msgLevel( MSG::DEBUG ) ) debug() << "# Tracks in " << inCont << " = " << tracks.size() << endmsg;
    const bool canHaveAncestors = ( "LHCb::TrackLocation::Tsa" != inCont );
    // loop over container
    for ( auto tr : tracks ) {
      // skip invalid tracks up front
      if ( tr->checkFlag( LHCb::Track::Flags::Invalid ) ) continue;
      // label ancestors as clones
      if ( canHaveAncestors ) {
        for ( auto anc : tr->ancestors() ) { anc->setFlag( LHCb::Track::Flags::Clone, true ); }
      }
      // skip clone tracks
      if ( tr->checkFlag( LHCb::Track::Flags::Clone ) ) continue;

      allTracks.push_back( const_cast<LHCb::Track*>( tr ) );
    }
  }
  if ( msgLevel( MSG::DEBUG ) ) debug() << "-> total # of tracks retrieved = " << allTracks.size() << endmsg;
}
