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

/** @class MoveRecVertexTracks MoveRecVertexTracks.h
 *
 *  Copy tracks participating in RecVertices to a new container and
 *  update the vertex' track pointers.
 *
 *  @author Rosen Matev
 *  @date   2016-04-16
 */

#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include <string>

class MoveRecVertexTracks : public GaudiAlgorithm {
public:
  /// Standard constructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution

private:
  DataObjectReadHandle<LHCb::RecVertices> m_vertexLocation{this, "VertexLocation", LHCb::RecVertexLocation::Velo3D};
  DataObjectWriteHandle<LHCb::Tracks>     m_outputLocation{this, "OutputLocation", ""};
};

DECLARE_COMPONENT( MoveRecVertexTracks )

StatusCode MoveRecVertexTracks::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize();
  if ( !sc.isFailure() && m_outputLocation.objKey().empty() ) { return Error( "OutputLocation is not specified" ); }
  return sc;
}

StatusCode MoveRecVertexTracks::execute() {
  auto vertices = m_vertexLocation.getIfExists();
  if ( !vertices ) { return Error( "Container " + m_vertexLocation.objKey() + " does not exist." ); }

  auto newTracks = m_outputLocation.put( std::make_unique<LHCb::Tracks>() );

  LHCb::RecVertex::TrackWithWeightVector newVertexTracks;

  for ( auto vertex : *vertices ) {
    auto nTracks = vertex->tracks().size();

    newVertexTracks.clear();
    newVertexTracks.reserve( nTracks );

    for ( size_t i = 0; i < nTracks; ++i ) {
      auto* newTrack = new LHCb::Track( *vertex->tracks()[i] );
      newTracks->insert( newTrack );
      newVertexTracks.emplace_back( newTrack, vertex->weights()[i] );
    }

    vertex->setTracksWithWeights( newVertexTracks );
  }

  return StatusCode::SUCCESS;
}
