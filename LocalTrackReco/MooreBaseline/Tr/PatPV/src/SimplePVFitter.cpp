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
#include "Event/State.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "IPVFitter.h"
#include "PVUtils.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

class SimplePVFitter : public extends<GaudiTool, IPVFitter> {

public:
  // Standard constructor
  using extends::extends;
  // Fitting
  StatusCode fitVertex( const Gaudi::XYZPoint& seedPoint, LHCb::span<const LHCb::Track* const> tracks,
                        LHCb::RecVertex& vtx, std::vector<const LHCb::Track*>& tracks2remove,
                        IGeometryInfo const& geometry ) const override;

private:
  Gaudi::Property<int>    m_minTr{this, "MinTracks", 5,
                               [=]( auto& ) {
                                 if ( this->m_minTr < 3 ) { // do not allow less than 3 tracks in a vertex
                                   this->m_minTr = 3;
                                   this->warning() << "MinTracks must be at least 3. Set to 3" << endmsg;
                                 }
                               },
                               "Minimum number of tracks to make a vertex"};
  Gaudi::Property<int>    m_Iterations{this, "Iterations", 30, "Number of iterations for minimisation"};
  Gaudi::Property<double> m_maxChi2{this, "maxChi2", 256.0, "Chi2 of completely wrong tracks"};
  Gaudi::Property<double> m_extrapRCut{this, "LinExtrapolation", 5.0 * Gaudi::Units::mm,
                                       "Radial limit for linear extrapolation"};
  Gaudi::Property<double> m_maxDeltaZ{this, "maxDeltaZ", 0.001 * Gaudi::Units::mm, "Fit convergence condition"};
  Gaudi::Property<double> m_acceptTrack{this, "acceptTrack", 0.00001, "Value of Tukey's weight to accept a track"};
  Gaudi::Property<double> m_trackMaxChi2Remove{this, "trackMaxChi2Remove", 36.,
                                               "Max chi2 tracks to be removed from next PV search"};

  // Extrapolators
  ToolHandle<ITrackExtrapolator> m_linExtrapolator{this, "LinearExtrapolator", "TrackLinearExtrapolator"};
  ToolHandle<ITrackExtrapolator> m_fullExtrapolator{this, "FullExtrapolator", "TrackMasterExtrapolator"};

  // Least square iterative PV fit
  StatusCode fit( LHCb::RecVertex& vtx, std::vector<PVTrack*>& pvTracks, std::vector<const LHCb::Track*>& tracks2remove,
                  IGeometryInfo const& geometry ) const;
  // Add track for PV
  StatusCode addTrackForPV( const LHCb::Track* str, std::vector<PVTrack>& pvTracks, double zseed,
                            IGeometryInfo const& geometry ) const;

  void initVertex( PVTracks& pvTracks, PVVertex& pvVertex, const Gaudi::XYZPoint seedPoint ) const;

  // Prepare hessian matrix and vectorD0
  void prepareVertex( LHCb::RecVertex& vtx, PVTrackPtrs& pvTracks, Gaudi::SymMatrix3x3& hess,
                      ROOT::Math::SVector<double, 3>& d0vec, int iter, IGeometryInfo const& geometry ) const;
  // Extrapolation
  StatusCode trackExtrapolate( PVTrack* pvTrack, const LHCb::RecVertex& vtx, IGeometryInfo const& geometry ) const;
  // Add track to fit
  void addTrack( PVTrack* pTrack, Gaudi::SymMatrix3x3& hessian, ROOT::Math::SVector<double, 3>& vectorD0 ) const;
  // Remove track from fit
  void removeTrack( PVTrack* pTrack, Gaudi::SymMatrix3x3& hessian, ROOT::Math::SVector<double, 3>& vectorD0 ) const;
  // Add subtrack
  void addsubTrack( PVTrack* pTrack, Gaudi::SymMatrix3x3& hessian, ROOT::Math::SVector<double, 3>& vectorD0,
                    double invs ) const;
  // Update matrices and vectors
  StatusCode outVertex( LHCb::RecVertex& vtx, PVTrackPtrs& pvTracks, Gaudi::SymMatrix3x3& hess,
                        ROOT::Math::SVector<double, 3>& d0vec ) const;
  // Set current chi2
  void setChi2( LHCb::RecVertex& vtx, PVTrackPtrs& pvTracks ) const;

  // Get Tukey's weight
  double getTukeyWeight( double trchi2, int iter ) const;
};

DECLARE_COMPONENT( SimplePVFitter )

//=============================================================================
// Least square adaptive fitting method
//=============================================================================
StatusCode SimplePVFitter::fitVertex( const Gaudi::XYZPoint& seedPoint, LHCb::span<const LHCb::Track* const> rTracks,
                                      LHCb::RecVertex& vtx, std::vector<const LHCb::Track*>& tracks2remove,
                                      IGeometryInfo const& geometry ) const {
  PVTracks pvTracks;
  PVVertex pvVertex;
  for ( auto itr = rTracks.begin(); itr != rTracks.end(); itr++ ) {
    const LHCb::Track* track = *itr;
    //    if ( !(track->hasVelo()) ) continue;
    StatusCode sc = addTrackForPV( track, pvTracks, seedPoint.z(), geometry );
    if ( !sc.isSuccess() ) continue;
  }
  initVertex( pvTracks, pvVertex, seedPoint );
  // Initial track cleaning
  auto& trks   = pvVertex.pvTracks;
  auto  itrack = trks.begin();
  while ( itrack != trks.end() ) {
    PVTrack*   pvTrack = *itrack;
    StatusCode sc      = trackExtrapolate( pvTrack, pvVertex.primVtx, geometry );
    //    if((!sc.isSuccess()) || (pvTrack->chi2 >= m_maxChi2)) {
    if ( !sc.isSuccess() ) {
      itrack = trks.erase( itrack );
    } else {
      ++itrack;
    }
  }

  // Check the number of tracks for PV candidate
  if ( (int)pvVertex.pvTracks.size() < m_minTr ) {
    if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "Too few tracks to fit PV" << endmsg; }
    vtx = pvVertex.primVtx;
    vtx.setTechnique( LHCb::RecVertex::RecVertexType::Primary );
    return StatusCode::FAILURE;
  }

  StatusCode scvfit = fit( pvVertex.primVtx, pvVertex.pvTracks, tracks2remove, geometry );
  if ( !scvfit.isSuccess() ) {
    if ( msgLevel( MSG::DEBUG ) ) debug() << "PV fit failed" << endmsg;
  }

  vtx = pvVertex.primVtx;
  vtx.setTechnique( LHCb::RecVertex::RecVertexType::Primary );
  return scvfit;
}

//=============================================================================
// Least square adaptive PV fit
//=============================================================================
StatusCode SimplePVFitter::fit( LHCb::RecVertex& vtx, std::vector<PVTrack*>& pvTracks,
                                std::vector<const LHCb::Track*>& tracks2remove, IGeometryInfo const& geometry ) const {
  if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "Least square adaptive PV fit" << endmsg; }
  // Reset hessian and d0 vector
  Gaudi::SymMatrix3x3            hess;
  ROOT::Math::SVector<double, 3> d0vec;
  for ( int i = 0; i < 3; i++ ) {
    d0vec[i] = 0.0;
    for ( int j = 0; j < 3; j++ ) { hess( i, j ) = 0.0; }
  }
  double zPrevious = 99999.0;
  double zVtx      = 0.0;
  double maxdz     = m_maxDeltaZ;
  bool   converged = false;
  int    nbIter    = 0;
  // Iteration loop. Require at least 3 iterations to reach final weight.
  while ( ( nbIter < 2 ) || ( !converged && nbIter < m_Iterations ) ) {
    if ( msgLevel( MSG::DEBUG ) ) { debug() << "Iteration nr: " << nbIter << endmsg; }
    zVtx = vtx.position().z();
    prepareVertex( vtx, pvTracks, hess, d0vec, nbIter, geometry );

    int ntr = 0;
    for ( PVTrackPtrs::iterator itrack = pvTracks.begin(); itrack != pvTracks.end(); itrack++ ) {
      PVTrack* pvTrack = *itrack;
      if ( pvTrack->weight > 0. ) ntr++;
    }
    if ( ntr < 3 ) {
      if ( msgLevel( MSG::DEBUG ) ) debug() << "# tracks too low. ntr = " << ntr << endmsg;
      break;
    }

    // low mutiplicity, could cycle over different tracks
    if ( nbIter > m_Iterations - 10 && ntr < 8 ) { maxdz = 10. * m_maxDeltaZ; }

    StatusCode sc = outVertex( vtx, pvTracks, hess, d0vec );
    if ( sc != StatusCode::SUCCESS ) { break; }
    zVtx = vtx.position().z();
    if ( fabs( zVtx - zPrevious ) < maxdz ) {
      converged = true;
    } else {
      zPrevious = zVtx;
    }

    ++nbIter;
  }

  if ( nbIter >= m_Iterations ) {
    if ( msgLevel( MSG::DEBUG ) ) debug() << " Reached max # iterations without convergence " << nbIter << endmsg;
  }

  setChi2( vtx, pvTracks );

  // Check the weight of tracks to be accepted for PV candidate
  int outTracks = 0;
  for ( PVTrackPtrs::iterator iFitTr = pvTracks.begin(); iFitTr != pvTracks.end(); iFitTr++ ) {
    PVTrack* fitTrack = *iFitTr;
    double   weight   = fitTrack->weight;
    if ( weight > m_acceptTrack ) { outTracks++; }
  }
  if ( outTracks < m_minTr ) {
    if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "Too few tracks after PV fit" << endmsg; }
    return StatusCode::FAILURE;
  }
  // remove vertices too far from beam axis
  double xva = vtx.position().x();
  double yva = vtx.position().y();
  if ( ( xva * xva + yva * yva ) > 5. * 5. ) { return StatusCode::FAILURE; }

  // fill tracks arounf PV to remove from next PV search
  for ( PVTrackPtrs::iterator itrack = pvTracks.begin(); itrack != pvTracks.end(); itrack++ ) {

    if ( ( *itrack )->chi2 < m_trackMaxChi2Remove ) { tracks2remove.push_back( ( *itrack )->refTrack ); }
  }

  return StatusCode::SUCCESS;
  ;
}

//=============================================================================
// Add track for PV
//=============================================================================
StatusCode SimplePVFitter::addTrackForPV( const LHCb::Track* pvtr, PVTracks& pvTracks, double zseed,
                                          IGeometryInfo const& geometry ) const {
  // Create a new PVTrack to be put on the vecctor
  PVTrack pvtrack;
  pvtrack.isUsed = false;
  pvtrack.chi2   = 0;
  pvtrack.d0sq   = 0;
  pvtrack.err2d0 = 0.0;
  // Keep reference to the original track
  pvtrack.refTrack = pvtr;
  LHCb::State mstate;
  mstate = pvtr->firstState();

  if ( ( mstate.checkLocation( LHCb::State::Location::ClosestToBeam ) != true ) || fabs( mstate.z() ) > 500. ) {
    // extrapolate
    if ( fabs( mstate.qOverP() ) > 0 ) {
      StatusCode sc = m_fullExtrapolator->propagate( mstate, zseed, geometry );
      if ( sc.isFailure() ) return sc;
    } else {
      StatusCode sc = m_linExtrapolator->propagate( mstate, zseed, geometry );
      if ( sc.isFailure() ) return sc;
    }
    // check if track is reasonable
    if ( ( mstate.x() * mstate.x() + mstate.y() * mstate.y() ) > 100. * 100. ) return StatusCode::FAILURE;
  }

  pvtrack.stateG   = mstate;
  pvtrack.unitVect = pvtrack.stateG.slopes().Unit();
  pvTracks.push_back( pvtrack );

  return StatusCode::SUCCESS;
}

//=============================================================================
// initialize vertex
//=============================================================================
void SimplePVFitter::initVertex( PVTracks& pvTracks, PVVertex& pvVtx, const Gaudi::XYZPoint seedPoint ) const {
  if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "initVertex method" << endmsg; }
  double zseed   = seedPoint.z();
  double startzs = zseed - 1000.0 * Gaudi::Units::mm;
  double endzs   = zseed + 1000.0 * Gaudi::Units::mm;
  int    nTracks = 0;
  for ( auto pvTrack = pvTracks.begin(); pvTracks.end() != pvTrack; pvTrack++ ) {
    if ( !pvTrack->isUsed ) {
      if ( ( pvTrack->zClose() >= startzs ) && ( pvTrack->zClose() <= endzs ) ) {
        nTracks++;
        pvVtx.pvTracks.push_back( &( *pvTrack ) );
        pvVtx.primVtx.addToTracks( pvTrack->refTrack );
      }
    }
  }

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << " Collected " << nTracks << " tracks out of " << pvTracks.size() << " for this vtx search" << endmsg;
  }
  pvVtx.primVtx.setPosition( seedPoint );
  pvVtx.primVtx.setNDoF( 2 * nTracks - 3 );
  Gaudi::SymMatrix3x3 hsm;
  for ( int i = 0; i < 3; i++ ) {
    for ( int j = 0; j < 3; j++ ) { hsm( i, j ) = 0.0; }
  }
  pvVtx.primVtx.setCovMatrix( hsm );
  pvVtx.primVtx.setChi2( 0.0 );
}

//=============================================================================
// Prepare vertices
//=============================================================================
void SimplePVFitter::prepareVertex( LHCb::RecVertex& vtx, PVTrackPtrs& pvTracks, Gaudi::SymMatrix3x3& hess,
                                    ROOT::Math::SVector<double, 3>& d0vec, int iter,
                                    IGeometryInfo const& geometry ) const {
  int nbIter = iter;
  // Reset hessian and d0 vector
  for ( int i = 0; i < 3; i++ ) {
    d0vec[i] = 0.0;
    for ( int j = 0; j < 3; j++ ) { hess( i, j ) = 0.0; }
  }
  if ( msgLevel( MSG::VERBOSE ) ) {
    verbose() << "Extrapolate tracks to the vertex at z = " << vtx.position().z() << endmsg;
  }
  // Add whole track vector to hessian and d0 vector
  for ( auto itrack = pvTracks.begin(); itrack != pvTracks.end(); ) {
    PVTrack* pvTrack = *itrack;
    // Extrapolate tracks
    StatusCode sc = trackExtrapolate( pvTrack, vtx, geometry );
    if ( !sc.isSuccess() ) {
      if ( msgLevel( MSG::VERBOSE ) ) {
        verbose() << "Track " << pvTrack->refTrack->key() << " could not be extrapolated to the vertex" << endmsg;
      }
      itrack = pvTracks.erase( itrack );
    } else {
      pvTrack->weight = getTukeyWeight( pvTrack->chi2, nbIter );
      addTrack( pvTrack, hess, d0vec );
      ++itrack;
    }
  }
}

//=============================================================================
// Track extrapolation
//=============================================================================
StatusCode SimplePVFitter::trackExtrapolate( PVTrack* pvTrack, LHCb::RecVertex const& vtx,
                                             IGeometryInfo const& geometry ) const {
  Gaudi::XYZPoint trkPoint( pvTrack->stateG.x(), pvTrack->stateG.y(), pvTrack->stateG.z() );
  // Set the new origin at the current vertex
  Gaudi::XYZVector diffVect = trkPoint - vtx.position();
  // Compute the distance from the vertex to the track
  pvTrack->vd0 = pvTrack->unitVect.Cross( diffVect.Cross( pvTrack->unitVect ) );
  // Compute impact parameter d0
  pvTrack->d0sq = pvTrack->vd0.Mag2();
  // Use linear extrapolator to move the track to the desired position
  LHCb::State stateToMove = pvTrack->stateG;
  StatusCode  sc = m_linExtrapolator->propagate( stateToMove, vtx.position().z() + pvTrack->vd0.z(), geometry );
  if ( !sc.isSuccess() ) {
    if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "Error propagating the track state" << endmsg; }
    return StatusCode::FAILURE;
  }
  // Compute the error on the track impact parameter
  Gaudi::XYZVector               vd0Unit = pvTrack->vd0.unit();
  ROOT::Math::SVector<double, 2> xyVec;
  xyVec[0]                                 = vd0Unit.x();
  xyVec[1]                                 = vd0Unit.y();
  Gaudi::SymMatrix2x2            covMatrix = stateToMove.covariance().Sub<Gaudi::SymMatrix2x2>( 0, 0 );
  ROOT::Math::SVector<double, 2> covMatrix_xyVec_product;
  covMatrix_xyVec_product = covMatrix * xyVec;
  pvTrack->err2d0         = xyVec[0] * covMatrix_xyVec_product[0] + xyVec[1] * covMatrix_xyVec_product[1];
  pvTrack->chi2           = ( pvTrack->d0sq ) / pvTrack->err2d0;
  return StatusCode::SUCCESS;
}

//=============================================================================
// Add track
//=============================================================================
void SimplePVFitter::addTrack( PVTrack* pTrack, Gaudi::SymMatrix3x3& hess,
                               ROOT::Math::SVector<double, 3>& d0vec ) const {
  addsubTrack( pTrack, hess, d0vec, +1.0 );
}

//=============================================================================
// Remove track
//=============================================================================
void SimplePVFitter::removeTrack( PVTrack* pTrack, Gaudi::SymMatrix3x3& hess,
                                  ROOT::Math::SVector<double, 3>& d0vec ) const {
  addsubTrack( pTrack, hess, d0vec, -1.0 );
}

//=============================================================================
// Add subtrack
//=============================================================================
void SimplePVFitter::addsubTrack( PVTrack* pvTrack, Gaudi::SymMatrix3x3& hess, ROOT::Math::SVector<double, 3>& d0vec,
                                  double invs ) const {
  invs *= ( 2.0 / pvTrack->err2d0 ) * pvTrack->weight;
  double unitVectStd[3];
  unitVectStd[0] = pvTrack->unitVect.x();
  unitVectStd[1] = pvTrack->unitVect.y();
  unitVectStd[2] = pvTrack->unitVect.z();
  double vd0Std[3];
  vd0Std[0] = pvTrack->vd0.x();
  vd0Std[1] = pvTrack->vd0.y();
  vd0Std[2] = pvTrack->vd0.z();
  for ( int k = 0; k < 3; k++ ) {
    double vd0 = 0.;
    for ( int l = 0; l < 3; l++ ) {
      double delta_kl = ( k == l ) ? 1.0 : 0.0;
      double val      = delta_kl - unitVectStd[k] * unitVectStd[l];
      vd0 += vd0Std[l] * val * invs;
      if ( l <= k ) { hess( k, l ) += val * invs; }
    }
    d0vec[k] += vd0;
  }
}

//=============================================================================
// Output Vertex
//=============================================================================
StatusCode SimplePVFitter::outVertex( LHCb::RecVertex& vtx, PVTrackPtrs& pvTracks, Gaudi::SymMatrix3x3& hess,
                                      ROOT::Math::SVector<double, 3>& d0vec ) const {
  if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "Computing new vertex position" << endmsg; }
  // Solve the linear equations to get a vertex
  int fail;
  hess.Inverse( fail );
  if ( 0 != fail ) {
    if ( msgLevel( MSG::DEBUG ) ) debug() << "Error inverting hessian matrix" << endmsg;
    return StatusCode::FAILURE;
  } else {
    hess.Invert();
  }
  ROOT::Math::SVector<double, 3> delta;
  for ( int i = 0; i < 3; i++ ) { delta[i] = 0.0; }
  delta = hess * d0vec;
  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Displacement found: " << delta[0] << " , " << delta[1] << " , " << delta[2] << " , " << endmsg;
  }
  Gaudi::XYZPoint newVtx( vtx.position().x() + delta[0], vtx.position().y() + delta[1], vtx.position().z() + delta[2] );
  vtx.setPosition( newVtx );
  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << " New vertex position: " << vtx.position().x() << " , " << vtx.position().y() << " , "
            << vtx.position().z() << endmsg;
  }
  // Set covariance matrix
  hess = hess * 2.0;
  vtx.setCovMatrix( hess );
  // Set tracks
  vtx.clearTracks();
  for ( PVTrackPtrs::iterator itrack = pvTracks.begin(); itrack != pvTracks.end(); itrack++ ) {
    if ( ( *itrack )->weight > m_acceptTrack ) { vtx.addToTracks( ( *itrack )->refTrack, (float)( *itrack )->weight ); }

    if ( msgLevel( MSG::DEBUG ) ) {
      PVTrack* pvtr = ( *itrack );
      debug() << format( "x,y,z: %9.3f %9.3f %9.3f  d0sq,errd0: %9.4f %9.4f  chi2,w: %9.3f %9.5f hasVelo %d ",
                         pvtr->stateG.x(), pvtr->stateG.y(), pvtr->stateG.z(), pvtr->d0sq, sqrt( pvtr->err2d0 ),
                         pvtr->chi2, pvtr->weight, pvtr->refTrack->hasVelo() )
              << endmsg;
    }
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
// Compute PV chi2
//=============================================================================
void SimplePVFitter::setChi2( LHCb::RecVertex& vtx, PVTrackPtrs& pvTracks ) const {
  int    nDoF = -3;
  double chi2 = 0.0;
  for ( PVTrackPtrs::iterator itrack = pvTracks.begin(); itrack != pvTracks.end(); itrack++ ) {
    if ( ( *itrack )->weight > m_acceptTrack ) {
      chi2 += ( *itrack )->chi2 * ( *itrack )->weight;
      nDoF += 2;
    }
  }
  if ( msgLevel( MSG::VERBOSE ) ) {
    verbose() << "Compute chi2 of this vertex: " << chi2 << " for " << nDoF << " DoF ( " << chi2 / nDoF << " / DoF)"
              << endmsg;
  }
  vtx.setChi2( chi2 );
  vtx.setNDoF( nDoF );
}

//=============================================================================
// Get Tukey's weight
//=============================================================================
double SimplePVFitter::getTukeyWeight( double trchi2, int iter ) const {
  if ( iter < 2 ) return 1.;

  double ctrv = 25. - 5. * ( iter - 2 );
  if ( ctrv < 5. ) ctrv = 5.;
  double cT2    = trchi2 / ( ctrv * ctrv );
  double weight = 0.;
  if ( cT2 < 1. ) {
    weight = 1. - cT2;
    weight = weight * weight;
  }
  return weight;
}
