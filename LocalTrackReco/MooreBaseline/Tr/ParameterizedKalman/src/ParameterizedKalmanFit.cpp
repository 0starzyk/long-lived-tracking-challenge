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
// Include files
#include "DetDesc/DetectorElement.h"
#include "Event/ODIN.h"
#include "Event/State.h"
#include "FTDet/DeFTDetector.h"
#include "Kernel/Trajectory.h"
#include "LHCbMath/Similarity.h"

// local
#include "ParameterizedKalmanFit.h"
#include "SerializeTrack.h"

//########################################################################
//
// Implementation file for class : ParameterizedKalmanFit
//
// 2017-10-26: Simon Stemmle
//
//########################################################################

using namespace ParKalman;

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( ParameterizedKalmanFit )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
ParameterizedKalmanFit::ParameterizedKalmanFit( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator,
                   {KeyValue{"InputName", "Rec/Track/ForwardFast"},
                    KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top},
                    KeyValue{"Magnet", LHCb::Det::Magnet::det_path}, KeyValue{"FT", DeFTDetectorLocation::Default}},
                   KeyValue{"OutputName", "Rec/Track/Best2"} ) {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode ParameterizedKalmanFit::initialize() {
  StatusCode sc = Transformer::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;           // error printed already by GaudiAlgorithm

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Initialize" << endmsg;

  info() << "Use parameters from " << m_ParamFileLocation.value() + "/Mag*" << endmsg;

  m_ParExtrUp   = std::make_unique<KalmanParametrizations>( m_ParamFileLocation, Polarity::Up, m_UseOneParameterSet );
  m_ParExtrDown = std::make_unique<KalmanParametrizations>( m_ParamFileLocation, Polarity::Down, m_UseOneParameterSet );

  // cache information for the smoother of outliers should be removed
  m_do_smoother = m_MaxNoutlier > 0;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
std::vector<ParameterizedKalmanFit::OutputTrack> ParameterizedKalmanFit::
                                                 operator()( std::vector<InputTrack> const& input, DetectorElement const& lhcb, DeMagnet const& magnet,
            DeFT const& deFT ) const {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;
  //============================================================
  //== Main processing: fit Tracks
  //============================================================

  // output tracks
  std::vector<OutputTrack> result;

  if ( input.empty() ) { return result; }

  // since we know the max size let's be safe and allocate a bit more to avoid reallocation
  result.reserve( input.size() );

  // fit status
  StatusCode sc;

  // Check the FT geometry version
  if ( deFT.version() < 62 ) {
    error() << "This version requires FTDet v6.2 or higher" << endmsg;
    throw std::runtime_error{"ParametrizedKalmanFit: This version requires FTDet v6.2 or higher"};
  }

  // struct that contains the intermediate track information
  trackInfo tI;
  tI.m_do_smoother = m_do_smoother;

  // select the respective extrapolator
  tI.m_extr = ( magnet.signedRelativeCurrent() > 0 ? m_ParExtrUp.get() : m_ParExtrDown.get() );

  // Loop over the tracks and fit them
  for ( const auto& inputTrack : input ) {
    // Create a new track
    auto& outputTrack = result.emplace_back();
    outputTrack.setType( OutputTrack::Types::Long );
    outputTrack.addToLhcbIDs( inputTrack.lhcbIDs() );
    // copy the states
    for ( const auto& state : inputTrack.states() ) { outputTrack.addToStates( state ); }

    tI.m_inputTrack = &inputTrack;
    tI.m_track      = &outputTrack;
    assert( tI.m_inputTrack );
    assert( tI.m_track );

    sc = fit( tI, *lhcb.geometry() );

    if ( !sc.isSuccess() ) result.pop_back();
  }

  return result;
}

//=============================================================================
//  Perform the fit
//=============================================================================
StatusCode ParameterizedKalmanFit::fit( trackInfo& tI, IGeometryInfo const& geometry ) const {
  // load hit information
  LoadHitsFromTrackV1( tI, m_measProvider, m_UseUT, m_UseT );

  // current state:
  // current z position
  double lastz = -1.;
  // state vector x: x,y,tx,ty,qop
  Gaudi::Vector5      x;
  Gaudi::SymMatrix5x5 C;

  // best state closest to the beam
  double              zBest = 0;
  Gaudi::Vector5      xBest;
  Gaudi::SymMatrix5x5 CBest;
  // the according best chi2
  double chi2Best = 0;

  // best state at first measurement
  double              zBestFirstMeas = 0;
  Gaudi::Vector5      xBestFirstMeas;
  Gaudi::SymMatrix5x5 CBestFirstMeas;

  // best state at last measurement
  double              zBestLastMeas = 0;
  Gaudi::Vector5      xBestLastMeas;
  Gaudi::SymMatrix5x5 CBestLastMeas;

  //#####################################
  // do a forward and a backward iteration
  //#####################################

  // create the seed state at the first VELO hit
  CreateVeloSeedState( 0, x, C, lastz, tI );

  // reset chi2 for the forward filtering
  tI.m_chi2  = 0;
  tI.m_chi2V = 0;

  // start by updating with the first measurment
  UpdateState( 1, 0, x, C, lastz, tI );

  // forward filtering
  for ( int nhit = 1; nhit < tI.m_NHitsTotal; nhit++ ) {
    if ( !PredictState( 1, nhit, x, C, lastz, tI ) ) return StatusCode::FAILURE;
    UpdateState( 1, nhit, x, C, lastz, tI );
  }

  // we have now a first best state at the laste measurement
  xBestLastMeas = x;
  CBestLastMeas = C;
  zBestLastMeas = lastz;

  // first best momentum estimate
  tI.m_BestMomEst = x[4];
  // take the chi2 from the forward filtering
  chi2Best = tI.m_chi2;

  // reset covariance matrix to represent "no information"
  C( 0, 0 ) = 400;
  C( 0, 1 ) = 0;
  C( 0, 2 ) = 0;
  C( 0, 3 ) = 0;
  C( 0, 4 ) = 0;
  C( 1, 1 ) = 400;
  C( 1, 2 ) = 0;
  C( 1, 3 ) = 0;
  C( 1, 4 ) = 0;
  C( 2, 2 ) = 0.01;
  C( 2, 3 ) = 0;
  C( 2, 4 ) = 0;
  C( 3, 3 ) = 0.01;
  C( 3, 4 ) = 0;
  C( 4, 4 ) = 25 * C( 4, 4 ); // TODO check this

  tI.m_chi2  = 0;
  tI.m_chi2T = 0;

  // backward filtering
  // start with an update using the information of the last hit
  UpdateState( -1, tI.m_NHitsTotal - 1, x, C, lastz, tI );
  for ( int nhit = tI.m_NHitsTotal - 2; nhit >= 0; nhit-- ) {
    if ( !PredictState( -1, nhit, x, C, lastz, tI ) ) return StatusCode::FAILURE;
    UpdateState( -1, nhit, x, C, lastz, tI );
  }

  // we have now a first best state at the first measurement
  xBestFirstMeas = x;
  CBestFirstMeas = C;
  zBestFirstMeas = lastz;
  // and for the extrapolation to the beam line
  xBest = x;
  CBest = C;
  zBest = lastz;

  // set the momentum estimate from the backward filtering (might be better)
  if ( !m_UseForwardMomEstimate ) tI.m_BestMomEst = xBest[4];
  // same for the chi2 (might be better) To test!
  if ( !m_UseForwardChi2Estimate ) chi2Best = tI.m_chi2;

  //##############################################
  // do outlier removal if requested and necessary
  //##############################################
  for ( int i = 0; i < m_MaxNoutlier; i++ ) {
    // start by running the smoother (average forward and backward)
    tI.m_chi2 = 0;
    for ( int nhit = 0; nhit < tI.m_NHitsTotal; nhit++ ) { AverageState( nhit, tI ); }
    // now check if there are outliers and remove the worst one
    if ( DoOutlierRemoval( tI ) ) {
      // do a forward iteration
      // reset chi2 for the forward filtering
      tI.m_chi2  = 0;
      tI.m_chi2V = 0;

      // set the state at the first hit
      lastz = tI.m_StateZPos[0];
      x     = tI.m_StateForwardUpdated[0];
      // reset covariance matrix to represent "no information"
      C( 0, 0 ) = 400;
      C( 0, 1 ) = 0;
      C( 0, 2 ) = 0;
      C( 0, 3 ) = 0;
      C( 0, 4 ) = 0;
      C( 1, 1 ) = 400;
      C( 1, 2 ) = 0;
      C( 1, 3 ) = 0;
      C( 1, 4 ) = 0;
      C( 2, 2 ) = 0.01;
      C( 2, 3 ) = 0;
      C( 2, 4 ) = 0;
      C( 3, 3 ) = 0.01;
      C( 3, 4 ) = 0;
      C( 4, 4 ) = 25 * C( 4, 4 ); // TODO check this

      // start by updating with the first measurment
      UpdateState( 1, 0, x, C, lastz, tI );

      // forward filtering
      for ( int nhit = 1; nhit < tI.m_NHitsTotal; nhit++ ) {
        if ( !PredictState( 1, nhit, x, C, lastz, tI ) ) return StatusCode::FAILURE;
        UpdateState( 1, nhit, x, C, lastz, tI );
      }

      // we have now a better best state at the laste measurement
      xBestLastMeas = x;
      CBestLastMeas = C;
      zBestLastMeas = lastz;

      // better best momentum estimate
      tI.m_BestMomEst = x[4];
      // take the chi2 from the forward filtering
      chi2Best = tI.m_chi2;

      // reset covariance matrix to represent "no information"
      C( 0, 0 ) = 400;
      C( 0, 1 ) = 0;
      C( 0, 2 ) = 0;
      C( 0, 3 ) = 0;
      C( 0, 4 ) = 0;
      C( 1, 1 ) = 400;
      C( 1, 2 ) = 0;
      C( 1, 3 ) = 0;
      C( 1, 4 ) = 0;
      C( 2, 2 ) = 0.01;
      C( 2, 3 ) = 0;
      C( 2, 4 ) = 0;
      C( 3, 3 ) = 0.01;
      C( 3, 4 ) = 0;
      C( 4, 4 ) = 25 * C( 4, 4 ); // TODO check this

      tI.m_chi2  = 0;
      tI.m_chi2T = 0;

      // backward filtering
      // start with an update using the information of the last hit
      UpdateState( -1, tI.m_NHitsTotal - 1, x, C, lastz, tI );
      for ( int nhit = tI.m_NHitsTotal - 2; nhit >= 0; nhit-- ) {
        if ( !PredictState( -1, nhit, x, C, lastz, tI ) ) return StatusCode::FAILURE;
        UpdateState( -1, nhit, x, C, lastz, tI );
      }

      // we have now a better best state at the first measurement
      xBestFirstMeas = x;
      CBestFirstMeas = C;
      zBestFirstMeas = lastz;
      // and for the extrapolation to the beam line
      xBest = x;
      CBest = C;
      zBest = lastz;

      // set the momentum estimate from the backward filtering (might be better)
      if ( !m_UseForwardMomEstimate ) tI.m_BestMomEst = xBest[4];
      // same for the chi2 (might be better) To test!
      if ( !m_UseForwardChi2Estimate ) chi2Best = tI.m_chi2;
    }
    // No outlier found
    else
      break;
  }

  // set best momentum and chi2
  xBest[4]  = tI.m_BestMomEst;
  tI.m_chi2 = chi2Best;
  // do a momentum scaling to account for energy loss:
  xBestFirstMeas[4] = tI.m_BestMomEst * ( 1. + 0.5 * std::fabs( tI.m_BestMomEst ) );
  xBestLastMeas[4]  = tI.m_BestMomEst * ( 1. + 30 * std::fabs( tI.m_BestMomEst ) );

  // extrapolate to the vertex
  ExtrapolateToVertex( xBest, CBest, zBest, geometry, m_extrapolator_toPV );
  // and create a new LHCb::Track
  addInfoToTrack( xBest, CBest, zBest, tI );
  // add some more states
  addStateToTrack( xBestFirstMeas, CBestFirstMeas, zBestFirstMeas, tI, LHCb::State::Location::FirstMeasurement );
  addStateToTrack( xBestLastMeas, CBestLastMeas, zBestLastMeas, tI, LHCb::State::Location::LastMeasurement );

  return StatusCode::SUCCESS;
}

//////////////////////////////////////////
// General method for predicting to a hit
//////////////////////////////////////////
bool ParameterizedKalmanFit::PredictState( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C,
                                           double& lastz, trackInfo& tI ) const {
  // success flag
  bool Succes = true;

  // Choose the appropiate predicting method depending on the detector
  // forward_________________________________________________________
  if ( forward > 0 ) {
    // Predict inside VELO
    if ( nHit < tI.m_NHitsV ) PredictStateV( nHit, x, C, lastz, tI );

    // Predict to first UT layer or directly to T
    else if ( nHit == tI.m_NHitsV ) {
      // To first UT hit
      Succes       = PredictStateVUT( x, C, lastz, tI );
      tI.m_PrevNUT = 0;
      // predict further in UT if there is no hit
      while ( tI.m_HasHitUT[tI.m_PrevNUT] == 0 && Succes && tI.m_PrevNUT < 3 ) {
        tI.m_PrevNUT++;
        PredictStateUT( tI.m_PrevNUT, x, C, lastz, tI );
      }
      // In case there is no UT hit, extrapolate to T
      if ( tI.m_NHitsUT == 0 && Succes ) {
        PredictStateUTT( x, C, lastz, tI );
        tI.m_PrevNT = 0;
        // predict further if there is no hit
        while ( tI.m_HasHitT[tI.m_PrevNT] == 0 ) {
          tI.m_PrevNT++;
          PredictStateT( tI.m_PrevNT, x, C, lastz, tI );
        }
      }
    }
    // Predict inside UT
    else if ( nHit < tI.m_NHitsV + tI.m_NHitsUT ) {
      // predict to next UT layer
      tI.m_PrevNUT++;
      PredictStateUT( tI.m_PrevNUT, x, C, lastz, tI );
      // predict further if there is no hit
      while ( tI.m_HasHitUT[tI.m_PrevNUT] == 0 && tI.m_PrevNUT < 3 ) {
        tI.m_PrevNUT++;
        PredictStateUT( tI.m_PrevNUT, x, C, lastz, tI );
      }
    }

    // Predict from UT to T
    else if ( nHit == tI.m_NHitsV + tI.m_NHitsUT ) {
      // check if we are at last UT station layer
      while ( tI.m_PrevNUT < 3 ) {
        tI.m_PrevNUT++;
        PredictStateUT( tI.m_PrevNUT, x, C, lastz, tI );
      }
      PredictStateUTT( x, C, lastz, tI );
      tI.m_PrevNT = 0;
      // predict further if there is no hit
      while ( tI.m_HasHitT[tI.m_PrevNT] == 0 ) {
        tI.m_PrevNT++;
        PredictStateT( tI.m_PrevNT, x, C, lastz, tI );
      }
    }

    // Predict inside T
    else if ( nHit < tI.m_NHitsTotal ) {
      // predict to next T layer
      tI.m_PrevNT++;
      PredictStateT( tI.m_PrevNT, x, C, lastz, tI );
      // predict further if there is no hit
      while ( tI.m_HasHitT[tI.m_PrevNT] == 0 && tI.m_PrevNT < 11 ) {
        tI.m_PrevNT++;
        PredictStateT( tI.m_PrevNT, x, C, lastz, tI );
      }
    }
  }
  // forward end_____________________________________________________

  // backwards_______________________________________________________
  else {
    // reset prevNT in case there was no forward prediction
    if ( nHit == tI.m_NHitsTotal - 2 ) tI.m_PrevNT = 11;
    // Predict inside T
    if ( nHit >= tI.m_NHitsV + tI.m_NHitsUT ) {
      // predict to next UT layer
      tI.m_PrevNT--;
      PredictStateT( tI.m_PrevNT, x, C, lastz, tI );
      // predict further if there is no hit
      while ( tI.m_HasHitT[tI.m_PrevNT] == 0 && tI.m_PrevNT > 0 ) {
        tI.m_PrevNT--;
        PredictStateT( tI.m_PrevNT, x, C, lastz, tI );
      }
    }

    // Predict to first UT layer or directly to VP
    else if ( nHit == tI.m_NHitsV + tI.m_NHitsUT - 1 ) {
      // To last UT hit
      PredictStateUTT( x, C, lastz, tI );
      tI.m_PrevNUT = 3;
      // predict further if there is no hit
      while ( tI.m_HasHitUT[tI.m_PrevNUT] == 0 && tI.m_PrevNUT > 0 ) {
        tI.m_PrevNUT--;
        PredictStateUT( tI.m_PrevNUT, x, C, lastz, tI );
      }
      // In case there is no UT hit, extrapolate to VP
      if ( tI.m_NHitsUT == 0 ) { Succes &= PredictStateVUT( x, C, lastz, tI ); }
    }

    // Predict inside UT
    else if ( nHit >= tI.m_NHitsV ) {
      // predict to next UT layer
      tI.m_PrevNUT--;
      PredictStateUT( tI.m_PrevNUT, x, C, lastz, tI );
      // predict further if there is no hit
      while ( tI.m_HasHitUT[tI.m_PrevNUT] == 0 && tI.m_PrevNUT > 0 ) {
        tI.m_PrevNUT--;
        PredictStateUT( tI.m_PrevNUT, x, C, lastz, tI );
      }
    }

    // Predict to VP
    else if ( nHit == tI.m_NHitsV - 1 ) {
      Succes &= PredictStateVUT( x, C, lastz, tI );
    }

    // Simple version for the VELO
    else if ( nHit < tI.m_NHitsV - 1 ) {
      PredictStateV( nHit, x, C, lastz, tI );
    }
  }
  // backwards end___________________________________________________

  // Save information for the smoother
  if ( tI.m_do_smoother ) {
    if ( forward > 0 ) {
      tI.m_StateForwardPredicted[nHit] = x;
      tI.m_CovForwardPredicted[nHit]   = C;
    } else {
      tI.m_StateBackwardPredicted[nHit] = x;
      tI.m_CovBackwardPredicted[nHit]   = C;
    }
  }
  return Succes;
}
