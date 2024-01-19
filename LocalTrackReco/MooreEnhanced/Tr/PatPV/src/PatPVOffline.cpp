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

#include "TrackInterfaces/IPVOfflineTool.h"

#include "DetDesc/DetectorElement.h"
#include "Event/RecVertex.h"
#include "Event/State.h"
#include "Event/Track.h"

#include "GaudiKernel/SystemOfUnits.h"
#include "LHCbAlgs/Transformer.h"

#include <assert.h>

/** @class PatPVOffline PatPVOffline.h
 *
 *
 *  @author Mariusz Witek
 *  @date   2010-10-05
 */
class PatPVOffline : public LHCb::Algorithm::Transformer<std::vector<LHCb::RecVertex>( std::vector<LHCb::Track> const&,
                                                                                       DetectorElement const& )> {
public:
  // Standard constructor
  PatPVOffline( const std::string& name, ISvcLocator* pSvcLocator );

  /// Algorithm execution
  std::vector<LHCb::RecVertex> operator()( std::vector<LHCb::Track> const&, DetectorElement const& ) const override;

private:
  // Tools
  ToolHandle<IPVOfflineTool> m_pvsfit{"PVOfflineTool", this};
};

DECLARE_COMPONENT( PatPVOffline )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PatPVOffline::PatPVOffline( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator,
                   {KeyValue{"InputTracks", LHCb::TrackLocation::Default},
                    KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}},
                   KeyValue( "OutputVertices", LHCb::RecVertexLocation::Primary ) ) {}

//=============================================================================
// Execution
//=============================================================================
std::vector<LHCb::RecVertex> PatPVOffline::operator()( std::vector<LHCb::Track> const& inputTracks,
                                                       DetectorElement const&          geometry ) const {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "Execute" << endmsg;

  std::vector<LHCb::RecVertex> rvts;
  StatusCode                   scfit = m_pvsfit->reconstructMultiPV( inputTracks, rvts, geometry );
  scfit.ignore();
  for ( auto& rvt : rvts ) { rvt.setTechnique( LHCb::RecVertex::RecVertexType::Primary ); }

  // ---> Debug
  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << endmsg;
    debug() << "TES location filled with " << rvts.size() << " PrimVertices" << endmsg;
    int nVtx = 0;
    for ( const auto& vertex : rvts ) {
      debug() << " Vertex " << nVtx << endmsg;
      debug() << " x, y, z: " << vertex.position().x() << " " << vertex.position().y() << " " << vertex.position().z()
              << endmsg;
      debug() << " Errors : " << sqrt( vertex.covMatrix()( 0, 0 ) ) << " " << sqrt( vertex.covMatrix()( 1, 1 ) ) << " "
              << sqrt( vertex.covMatrix()( 2, 2 ) ) << endmsg;
      debug() << " Number of tracks: " << vertex.tracks().size() << endmsg;
      debug() << " Chi2/DoF: " << vertex.chi2() / vertex.nDoF() << endmsg;
    }
    debug() << endmsg;
  }
  return rvts;
}
