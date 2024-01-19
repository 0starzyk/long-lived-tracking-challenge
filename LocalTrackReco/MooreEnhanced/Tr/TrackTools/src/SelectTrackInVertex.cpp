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
#include "Event/RecVertex.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "TrackSelectorBase.h"

//-----------------------------------------------------------------------------
// Implementation file for class : SelectTrackInVertex
//
// 2015-09-18 : Michel De Cian
//
//-----------------------------------------------------------------------------

/** @class SelectTrackInVertex SelectTrackInVertex.h
 *
 *  Selector that returns true if a track is found in a collection of tracks belonging to a vertex (LHCb::RecVertex,
 * e.g. the PV)
 *
 *  Parameters:
 *
 *  - VertexContainer: Container of vertices
 *
 *  @author Michel De Cian
 *  @date   2015-09-18
 */
class SelectTrackInVertex : public extends<TrackSelectorBase, IIncidentListener> {

public:
  /// constructor
  using extends::extends;

  StatusCode initialize() override;

  /** Returns if the given track is selected or not
   *
   *  @param aTrack Reference to the Track to test
   *
   *  @return boolean indicating if the track is selected or not
   *  @retval true  Track is selected
   *  @retval false Track is rejected
   */
  bool accept( const LHCb::Track& aTrack ) const override;

  void handle( const Incident& incident ) override;

private:
  /// Get all tracks belonging to all vertices and put them in a container
  void getTracksFromVertices() const;

  Gaudi::Property<std::string> m_vertexContainerName{this, "VertexContainer", LHCb::RecVertexLocation::Primary};
  mutable bool                 m_newEvent;
  mutable std::vector<const LHCb::Track*> m_tracks;
};

DECLARE_COMPONENT( SelectTrackInVertex )

//=============================================================================
// -- Initialize
//=============================================================================
StatusCode SelectTrackInVertex::initialize() {

  StatusCode sc = TrackSelectorBase::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;                 // error printed already by GaudiAlgorithm

  incSvc()->addListener( this, IncidentType::BeginEvent );
  m_newEvent = true;

  return StatusCode::SUCCESS;
}
//=============================================================================
// -- Check if track is in track collection belonging to a vertex
//=============================================================================
bool SelectTrackInVertex::accept( const LHCb::Track& aTrack ) const {

  if ( m_newEvent ) getTracksFromVertices();

  bool found = std::binary_search( std::begin( m_tracks ), std::end( m_tracks ), &aTrack );
  if ( !found && msgLevel( MSG::DEBUG ) ) {
    debug() << "Track " << &aTrack << " does not belong to a vertex" << endmsg;
  }
  return found;
}
//=============================================================================
// -- Check if new event has occurred. If yes, set flag
// -- Note: The actions of initEvent cannot be executed here,
// -- as this handle runs before decoding the clusters
//=============================================================================
void SelectTrackInVertex::handle( const Incident& incident ) {

  if ( IncidentType::BeginEvent == incident.type() ) m_newEvent = true;
}
//=============================================================================
// -- Get vertices, once per event
//=============================================================================
void SelectTrackInVertex::getTracksFromVertices() const {

  m_tracks.clear();

  const LHCb::RecVertices* vertexContainer = getIfExists<LHCb::RecVertices>( m_vertexContainerName );
  if ( !vertexContainer ) return;

  for ( const LHCb::RecVertex* vert : *vertexContainer ) {
    if ( !vert ) continue;
    const SmartRefVector<LHCb::Track> vTracks = vert->tracks();
    m_tracks.insert( m_tracks.end(), vTracks.begin(), vTracks.end() );
  }
  std::sort( std::begin( m_tracks ), std::end( m_tracks ) );

  m_newEvent = false;
}
