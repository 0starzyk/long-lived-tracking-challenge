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
//-----------------------------------------------------------------------------
/** @class CountVeloTracks CountVeloTracks.h
 *  @file CountVeloTracks.cpp
 *
 *  Tool for counting the distinct VELO tracks in an event
 *
 *  Based on code by Matt Needham
 *
 *  @author David Hutchcroft David.Hutchcroft@cern.ch
 *  @date   21/01/2011
 */
//-----------------------------------------------------------------------------

#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "Kernel/ICountContainedObjects.h"

struct CountVeloTracks : extends<GaudiTool, ICountContainedObjects> {

  /// constructor
  using extends::extends;

  /** Returns number of distinct VELO tracks contributing to the container
   *
   *  @param tracks reference to LHCb::Tracks container
   *
   *  @return number of distinct VELO tracks
   */
  unsigned int nObj( const ObjectContainerBase* cont ) const override;
};

DECLARE_COMPONENT( CountVeloTracks )

namespace {

  /// get a vector of LHCbID for VELO only from a track pointer
  std::vector<LHCb::LHCbID> getVeloIDs( const LHCb::Track& track ) {
    auto vids = track.lhcbIDs();
    vids.erase( std::remove_if( vids.begin(), vids.end(), []( const LHCb::LHCbID& id ) { return !id.isVP(); } ),
                vids.end() );
    return vids;
  }

} // namespace

//-----------------------------------------------------------------------------

unsigned int CountVeloTracks::nObj( const ObjectContainerBase* cont ) const {

  const LHCb::Tracks* tracks = dynamic_cast<const LHCb::Tracks*>( cont );
  if ( !tracks ) {
    error() << "Input is not an LHCb::Tracks container" << endmsg;
    return 0;
  }

  // first is the first VELO LHCbID on the track, second a track pointer
  std::multimap<LHCb::LHCbID, const LHCb::Track*> localCont;
  // somewhere to put tracks with VELO parts

  for ( const auto& track : *tracks ) {
    if ( !track->hasVelo() ) continue; // skip tracks without VELO components

    auto cvids = getVeloIDs( *track );
    // keep track if not a match to an existing track
    auto posRange = localCont.equal_range( cvids.front() );
    // loop over tracks with same first velo hit (if any)
    bool isClonedVelo =
        std::any_of( posRange.first, posRange.second, [&]( const std::pair<LHCb::LHCbID, const LHCb::Track*>& i ) {
          auto tvids = getVeloIDs( *i.second );
          // FIXME/TODO should we not do std::equal over a # of elements corresponding to the shortest of the two
          // ranges?
          return std::equal( cvids.begin(), cvids.end(), tvids.begin() );
        } );

    if ( !isClonedVelo ) { // add this track
      localCont.insert( std::make_pair( cvids.front(), track ) );
    }
  }
  return localCont.size();
}
