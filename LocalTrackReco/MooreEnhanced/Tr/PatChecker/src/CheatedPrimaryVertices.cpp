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
// Include files:
// from Gaudi
#include "GaudiKernel/SystemOfUnits.h"
// from Event
#include "Event/MCTrackInfo.h"
#include "Event/RecVertex.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "Linker/LinkedTo.h"
// from gsl
#include "gsl/gsl_math.h"
// Local
#include "CheatedPrimaryVertices.h"

//-----------------------------------------------------------------------------
// Implementation file for class : CheatedPrimaryVertices
//
// 2008-12-07 : Marcin Kucharczyk
//-----------------------------------------------------------------------------

DECLARE_COMPONENT( CheatedPrimaryVertices )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
CheatedPrimaryVertices::CheatedPrimaryVertices( const std::string& name, ISvcLocator* pSvcLocator )
    : GaudiAlgorithm( name, pSvcLocator ), m_inputTracks( NULL ), m_outputVertices( NULL ) {
  declareProperty( "InputTracksName", m_inputTracksName = "Rec/Track/Best" );
  declareProperty( "OutputVerticesName", m_outputVerticesName = "Hlt/Vertex/PV2D" );
}

//=============================================================================
// Execution
//=============================================================================
StatusCode CheatedPrimaryVertices::execute() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "Execute" << endmsg;

  auto const* track_links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::TrackLocation::Default ) );
  m_inputTracks           = get<LHCb::Tracks>( m_inputTracksName );
  m_outputVertices        = new LHCb::RecVertices();
  put( m_outputVertices, m_outputVerticesName );
  // MC particles assigned to hasVelo tracks
  int                                  nrTracks = 0;
  SmartRefVector<LHCb::Track>          hasVeloTracks;
  std::vector<const LHCb::MCParticle*> mcParts;
  for ( const LHCb::Track* trk : *m_inputTracks ) {
    nrTracks++;
    if ( trk->hasVelo() ) {
      hasVeloTracks.push_back( trk );
      LHCb::MCParticle const* mcPart = LinkedTo<LHCb::MCParticle>{track_links}.range( trk->key() ).try_front();
      if ( mcPart ) { mcParts.push_back( mcPart ); }
    }
  }
  // Visible MC PV's
  LHCb::MCVertices*                  mcVertices = get<LHCb::MCVertices>( LHCb::MCVertexLocation::Default );
  std::vector<const LHCb::MCVertex*> vtcsMCPV;
  int                                nrMCPVs = 0;
  for ( const auto* mcVtx : *mcVertices ) {
    nrMCPVs++;
    const LHCb::MCParticle* motherPart = mcVtx->mother();
    if ( 0 == motherPart ) {
      if ( mcVtx->type() == LHCb::MCVertex::MCVertexType::ppCollision ) {
        int isfromMCPV = 0;
        for ( const LHCb::MCParticle* mcPart : mcParts ) { isfromMCPV = isfromMCPV + fromMCVertex( mcPart, mcVtx ); }
        if ( isfromMCPV > 4 ) { vtcsMCPV.push_back( mcVtx ); }
      }
    }
  }
  // Cheated PV's
  double sigmax = 0.020 * Gaudi::Units::mm;
  double sigmay = 0.020 * Gaudi::Units::mm;
  double sigmaz = 0.080 * Gaudi::Units::mm;
  for ( std::vector<const LHCb::MCVertex*>::iterator itPV = vtcsMCPV.begin(); itPV != vtcsMCPV.end(); itPV++ ) {
    const LHCb::MCVertex* mcVtx     = *itPV;
    LHCb::RecVertex*      cheatedPV = new LHCb::RecVertex();
    cheatedPV->clearTracks();
    Gaudi::XYZPoint cheatedPVposition( mcVtx->position().x(), mcVtx->position().y(), mcVtx->position().z() );
    cheatedPV->setPosition( cheatedPVposition );
    Gaudi::SymMatrix3x3 cheatedPVcov;
    cheatedPVcov( 0, 0 ) = gsl_pow_2( sigmax );
    cheatedPVcov( 1, 1 ) = gsl_pow_2( sigmay );
    cheatedPVcov( 2, 2 ) = gsl_pow_2( sigmaz );
    cheatedPV->setCovMatrix( cheatedPVcov );
    int                         nDoF    = -3;
    double                      chi2Vtx = 0.0;
    SmartRefVector<LHCb::Track> cheatedPVTracks;
    for ( LHCb::Track* vTr : hasVeloTracks ) {
      LHCb::MCParticle const* mcPart = LinkedTo<LHCb::MCParticle>{track_links}.range( vTr->key() ).try_front();
      if ( mcPart && fromMCVertex( mcPart, mcVtx ) ) {
        chi2Vtx += chi2Tr( vTr, cheatedPV );
        nDoF += 2;
        cheatedPVTracks.push_back( vTr );
      }
    }
    cheatedPV->setChi2( chi2Vtx );
    cheatedPV->setNDoF( nDoF );
    cheatedPV->setTracks( cheatedPVTracks );
    cheatedPV->setTechnique( LHCb::RecVertex::RecVertexType::Primary );
    m_outputVertices->insert( cheatedPV );
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
// MC particle from MC visible PV
//=============================================================================
int CheatedPrimaryVertices::fromMCVertex( const LHCb::MCParticle* mcParticle, const LHCb::MCVertex* mcVertex ) {
  int                     isDaugh  = 0;
  const LHCb::MCVertex*   mcVtx    = 0;
  const LHCb::MCParticle* motherMC = mcParticle->mother();
  while ( motherMC ) {
    mcVtx    = motherMC->originVertex();
    motherMC = motherMC->mother();
  }
  if ( mcVertex == mcVtx ) isDaugh = 1;
  return isDaugh;
}

//=============================================================================
// Track chi2
//=============================================================================
double CheatedPrimaryVertices::chi2Tr( const LHCb::Track* pvTr, const LHCb::RecVertex* invt ) {
  // Compute impact parameter
  Gaudi::XYZVector unitVect;
  Gaudi::XYZVector vd0;
  Gaudi::XYZPoint  trkPoint( pvTr->firstState().x(), pvTr->firstState().y(), pvTr->firstState().z() );
  Gaudi::XYZVector diffVect = trkPoint - invt->position();
  unitVect                  = pvTr->firstState().slopes().Unit();
  vd0                       = unitVect.Cross( diffVect.Cross( unitVect ) );
  double d02                = vd0.Mag2();
  // Compute the error on the track impact parameter
  Gaudi::XYZVector               vd0Unit = vd0.Unit();
  ROOT::Math::SVector<double, 2> xyVec;
  xyVec[0]                                 = vd0Unit.x();
  xyVec[1]                                 = vd0Unit.y();
  Gaudi::SymMatrix2x2            covMatrix = pvTr->firstState().covariance().Sub<Gaudi::SymMatrix2x2>( 0, 0 );
  ROOT::Math::SVector<double, 2> product;
  product       = covMatrix * xyVec;
  double err2d0 = xyVec[0] * product[0] + xyVec[1] * product[1];
  double chi2   = d02 / err2d0;
  return chi2;
}
