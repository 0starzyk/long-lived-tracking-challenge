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
#include "Event/TrackTypes.h"
#include "Event/TrackUnitsConverters.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "TrackInterfaces/ITrackChi2Calculator.h"

using namespace Gaudi;

//-----------------------------------------------------------------------------
// Implementation file for class : TrackChi2Calculator
//
// 2003-09-18 : Jeroen van Tilburg
//-----------------------------------------------------------------------------

/** @class TrackChi2Calculator TrackChi2Calculator.h "TrackTools/TrackChi2Calculator.h"
 *
 *  This tool can be used to calculate the chi2-distance between
 *  two track states.
 *
 *  @author Jeroen van Tilburg
 *  @date   2003-09-18
 */

class TrackChi2Calculator : public extends<GaudiTool, ITrackChi2Calculator> {
public:
  /// Standard constructor
  using extends::extends;

  /** Calculate the chi2 distance between two track vectors.
   *  The track vectors must be given as (x,y,tx,ty,q/p).
   *  @return StatusCode: Failure if matrix inversion failed
   *  @param  trackVector1 input 1st track HepVector
   *  @param  trackCov1 input covariance matrix corresponding to 1st vector
   *  @param  trackVector2 input 2nd track HepVector
   *  @param  trackCov2 input covariance matrix corresponding to 2nd vector
   *  @param  chi2 output chi2 distance between the two vectors
   */
  StatusCode calculateChi2( const Gaudi::TrackVector& trackVector1, const Gaudi::TrackSymMatrix& trackCov1,
                            const Gaudi::TrackVector& trackVector2, const Gaudi::TrackSymMatrix& trackCov2,
                            double& chi2 ) const override;

private:
  StatusCode calculateChi2( Gaudi::Vector4& trackVector1, Gaudi::Vector4& trackVector2, Gaudi::SymMatrix4x4& trackCov12,
                            double& chi2 ) const;

  /// invert the 5x5 matrix
  StatusCode invertMatrix( Gaudi::TrackSymMatrix& invC ) const;

  /// invert the 4x4 matrix
  StatusCode invertMatrix( Gaudi::SymMatrix4x4& invC ) const;

  // job options
  // -----------
  /// Re-scale the track covariance matrices with this vector
  Gaudi::Property<std::vector<double>> m_scaleVector{this, "ScaleVector"};
  /// Remove Tx from chi2 contribution in case of matching inside the magnet
  /// with straight lines
  Gaudi::Property<bool> m_matchInMagnet{this, "MatchInMagnet", false};
  /// Add momentum information to the matching chi2 criterium
  Gaudi::Property<bool> m_addMomentum{this, "AddMomentum", false};
};

DECLARE_COMPONENT( TrackChi2Calculator )

//=============================================================================
//
//=============================================================================
StatusCode TrackChi2Calculator::calculateChi2( const Gaudi::TrackVector&    trackVector1,
                                               const Gaudi::TrackSymMatrix& trackCov1,
                                               const Gaudi::TrackVector&    trackVector2,
                                               const Gaudi::TrackSymMatrix& trackCov2, double& chi2 ) const {
  if ( !m_addMomentum ) { // then the dimension is 4
    Vector4      vec1       = trackVector1.Sub<Vector4>( 0 );
    Vector4      vec2       = trackVector2.Sub<Vector4>( 0 );
    SymMatrix4x4 trackCov12 = trackCov1.Sub<SymMatrix4x4>( 0, 0 ) + trackCov2.Sub<SymMatrix4x4>( 0, 0 );
    return calculateChi2( vec1, vec2, trackCov12, chi2 );
  }

  // If momentum information is to be used in the matching
  // -----------------------------------------------------
  // initialize chi2
  chi2 = 0.0;

  // copy tracks info
  TrackVector    vec1      = TrackVector( trackVector1 );
  TrackVector    vec2      = TrackVector( trackVector2 );
  TrackSymMatrix trackCinv = trackCov1 + trackCov2;

  // invert the matrix
  StatusCode sc = invertMatrix( trackCinv );
  if ( sc.isFailure() ) return StatusCode::FAILURE;

  // Remove Tx from chi2 in case of matching inside the magnet
  if ( m_matchInMagnet ) {
    vec1[2] = 0.0;
    vec2[2] = 0.0;
  }

  // Re-scale the chi2-contributions in case of error under/over-estimation
  unsigned int scaleVectorSize = m_scaleVector.size();
  if ( scaleVectorSize > 0 ) {
    for ( unsigned int i = 0; i < 5 && i <= scaleVectorSize; ++i ) {
      vec1[i] *= sqrt( fabs( m_scaleVector[i] ) );
      vec2[i] *= sqrt( fabs( m_scaleVector[i] ) );
    }
  }

  // Calculate the chi2 distance between 2 tracks
  chi2 = ROOT::Math::Similarity<double, 5>( vec1 - vec2, trackCinv );

  return StatusCode::SUCCESS;
}

//=============================================================================
//
//=============================================================================
StatusCode TrackChi2Calculator::calculateChi2( Gaudi::Vector4& trackVector1, Gaudi::Vector4& trackVector2,
                                               Gaudi::SymMatrix4x4& trackCov12, double& chi2 ) const {
  // initialize chi2
  chi2 = 0.0;

  // invert the matrix
  StatusCode sc = invertMatrix( trackCov12 );
  if ( sc.isFailure() ) return StatusCode::FAILURE;

  // Remove Tx from chi2 in case of matching inside the magnet
  if ( m_matchInMagnet ) {
    trackVector1[2] = 0.0;
    trackVector2[2] = 0.0;
  }

  // Re-scale the chi2-contributions in case of error under/over-estimation
  int scaleVectorSize = m_scaleVector.size();
  if ( scaleVectorSize > 0 ) {
    for ( int i = 0; i < 4 && i <= scaleVectorSize; ++i ) {
      trackVector1[i] *= sqrt( fabs( m_scaleVector[i] ) );
      trackVector2[i] *= sqrt( fabs( m_scaleVector[i] ) );
    }
  }

  // Calculate the chi2 distance between 2 tracks
  chi2 = ROOT::Math::Similarity<double, 4>( trackVector1 - trackVector2, trackCov12 );

  return StatusCode::SUCCESS;
}

//=============================================================================
//
//=============================================================================
StatusCode TrackChi2Calculator::invertMatrix( Gaudi::TrackSymMatrix& invC ) const {
  // This routine is taken from TrKalmanSmoother.cpp. It rescales
  // the matrix before it actually calls the DSINV wrapper.

  // Invert previous node covariance matrix
  // What follows may seem strange - trust me it works - you
  // are strongly recommended NOT to change it. It turns out that
  // the choice of MeV, mm as units is BAD - the inversion simply fails
  // for numerical reasons. Therefore it is necessary to change back to G3
  // units, invert then go back to G4 units
  // M. Needham 13/6/2000

  // check that the elements are not too large else dsinv will crash
  for ( unsigned int i = 0; i < 4; ++i ) {
    for ( unsigned int j = 0; j < 4; ++j ) {
      if ( invC( i, j ) > 1e20 ) { return Warning( "old covariance errors too big to invert", StatusCode::FAILURE ); }
    }
  }

  // G3 units
  TrackUnitsConverters::convertToG3( invC );

  bool OK = invC.Invert();

  // G4 units
  TrackUnitsConverters::convertToG4( invC );

  return OK ? StatusCode::SUCCESS : StatusCode::FAILURE;
}

//=============================================================================
//
//=============================================================================
StatusCode TrackChi2Calculator::invertMatrix( Gaudi::SymMatrix4x4& invC ) const {
  // This routine is taken from TrKalmanSmoother.cpp. It rescales
  // the matrix before it actually calls the DSINV wrapper.

  // Invert previous node covariance matrix
  // What follows may seem strange - trust me it works - you
  // are strongly recommended NOT to change it. It turns out that
  // the choice of MeV, mm as units is BAD - the inversion simply fails
  // for numerical reasons. Therefore it is necessary to change back to G3
  // units, invert then go back to G4 units
  // M. Needham 13/6/2000

  // check that the elements are not too large else dsinv will crash
  for ( unsigned int i = 0; i < 4; ++i ) {
    for ( unsigned int j = 0; j < 4; ++j ) {
      if ( invC( i, j ) > 1e20 ) { return Warning( "old covariance errors too big to invert", StatusCode::FAILURE ); }
    }
  }

  // G3 units
  TrackUnitsConverters::convertToG3( invC );

  bool OK = invC.Invert();

  // G4 units
  TrackUnitsConverters::convertToG4( invC );

  return OK ? StatusCode::SUCCESS : StatusCode::FAILURE;
}

//=============================================================================
