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
#include "Event/TrackFitResult.h"
#include "LHCbMath/Similarity.h"
#include "TrackKernel/TrackTraj.h"

#include "ParameterizedKalmanFit_Methods.h"

//########################################################################
//
// Implementation of methods defined in ParameterizedKalmanFit_Methods.h
//
// 2017-10-30: Simon Stemmle
//
//########################################################################

namespace ParKalman {

  ////////////////////////////////////////
  // Load hit information
  ////////////////////////////////////////
  void LoadHits( trackInfo& tI, bool m_UseUT, bool m_UseT ) {
    // reset counting variables
    tI.m_NHitsV     = 0;
    tI.m_NHitsUT    = 0;
    tI.m_NHitsT     = 0;
    tI.m_NHitsTotal = 0;
    tI.m_HasHitUT.fill( 0 );
    tI.m_HasHitT.fill( 0 );
    tI.m_UTHitToUTLayer.fill( 0 );
    tI.m_THitToTLayer.fill( 0 );

    // Load Velo hits directly from TrackHits.
    double lastVeloZ = -1e4;

    for ( const auto& trackHit : tI.m_inputTrack->veloHits ) {
      if ( std::abs( trackHit.beginPoint.z() ) - lastVeloZ < 2 ) continue;
      lastVeloZ = trackHit.beginPoint.z();

      tI.m_XMeasV[tI.m_NHitsV] = trackHit.beginPoint.x();
      tI.m_YMeasV[tI.m_NHitsV] = trackHit.beginPoint.y();
      tI.m_ZMeasV[tI.m_NHitsV] = trackHit.beginPoint.z();

      tI.m_XErrV[tI.m_NHitsV] = trackHit.errorx;
      tI.m_YErrV[tI.m_NHitsV] = trackHit.errory;

      tI.m_NHitsTotal++;
      tI.m_NHitsV++;
    }

    // Load UT hits directly from TrackHits.
    if ( m_UseUT ) {
      for ( const auto& trackHit : tI.m_inputTrack->utHits ) {
        if ( tI.m_HasHitUT[trackHit.layer] == 1 ) { continue; }
        tI.m_measUtLite[trackHit.layer].point     = trackHit.beginPoint;
        tI.m_measUtLite[trackHit.layer].direction = ( trackHit.endPoint - trackHit.beginPoint ).Unit();
        tI.m_measUtLite[trackHit.layer].error     = trackHit.errorx;

        tI.m_NHitsTotal++;
        tI.m_NHitsUT++;
        tI.m_HasHitUT[trackHit.layer] = 1;
      }
    }

    // Load FT hits directly from TrackHits.
    if ( m_UseT ) {
      for ( const auto& trackHit : tI.m_inputTrack->ftHits ) {
        if ( tI.m_HasHitT[trackHit.layer] == 1 ) { continue; }
        tI.m_measFtLite[trackHit.layer].point     = trackHit.beginPoint;
        tI.m_measFtLite[trackHit.layer].direction = ( trackHit.endPoint - trackHit.beginPoint ).Unit();
        tI.m_measFtLite[trackHit.layer].error     = trackHit.errorx;

        tI.m_NHitsTotal++;
        tI.m_NHitsT++;
        tI.m_HasHitT[trackHit.layer] = 1;
      }
    }

    // set ndof
    tI.m_Ndof   = tI.m_NHitsT + tI.m_NHitsUT + tI.m_NHitsV * 2;
    tI.m_NdofT  = tI.m_NHitsT;
    tI.m_NdofUT = tI.m_NHitsUT;
    tI.m_NdofV  = 2 * tI.m_NHitsV;

    // Set the hit status of all hits to 1: Should be used
    tI.m_HitStatus.fill( 1 );
    // create maps for the nHit -> layer matching
    // T stations
    int counter = 0;
    for ( int i = 0; i < 12; i++ ) {
      if ( tI.m_HasHitT[i] ) {
        tI.m_THitToTLayer[counter] = i;
        counter++;
      }
    }
    // UT
    counter = 0;
    for ( int i = 0; i < 4; i++ ) {
      if ( tI.m_HasHitUT[i] ) {
        tI.m_UTHitToUTLayer[counter] = i;
        counter++;
      }
    }
  }

  void LoadHitsFromTrackV1( trackInfo& tI, const ToolHandle<IMeasurementProviderProjector>& m_measProvider,
                            bool m_UseUT, bool m_UseT ) {
    // reset counting variables
    tI.m_NHitsV     = 0;
    tI.m_NHitsUT    = 0;
    tI.m_NHitsT     = 0;
    tI.m_NHitsTotal = 0;
    tI.m_HasHitUT.fill( 0 );
    tI.m_HasHitT.fill( 0 );
    tI.m_UTHitToUTLayer.fill( 0 );
    tI.m_THitToTLayer.fill( 0 );

    // check for double hits in one layer
    double lastVeloZ = -1e4;

    // Load (create) all measurements of the track
    // m_measProvider->load(*(tI.m_track));

    // const std::vector<LHCb::Measurement> measurements =
    // ((LHCb::TrackFitResult*)tI.m_track->fitResult())->measurements();

    LHCb::TrackFitResult::MeasurementContainer measurements;

    auto n2d = std::count_if( tI.m_track->lhcbIDs().begin(), tI.m_track->lhcbIDs().end(),
                              []( LHCb::LHCbID id ) { return id.isVP() || id.isMuon(); } );
    measurements.reserve( tI.m_track->lhcbIDs().size() + n2d ); // note: count VP & Muon, as they have 2 measurments per
                                                                // id

    std::vector<LHCb::LHCbID> ids = tI.m_track->lhcbIDs();

    m_measProvider->addToMeasurements( ids, measurements, LHCb::TrackTraj{*tI.m_track} );

    auto itM    = measurements.begin();
    auto itMEnd = measurements.end();
    while ( !( itM == itMEnd ) ) {
      itM->visit(
          // Load VP hit__________________________________________
          [&]( const LHCb::Measurement::VP& m ) {
            // check if there was a previous hit in this layer
            if ( std::abs( lastVeloZ - m.trajectory.beginPoint().Z() ) < 2 ) {
              ++itM; // VP come in pairs
              return;
            }
            // x measure
            assert( m.projection() == LHCb::Measurement::VP::Projection::X );
            [[maybe_unused]] auto lhcbID_X = itM->lhcbID();
            [[maybe_unused]] auto z_X      = itM->z();

            tI.m_XMeasV[tI.m_NHitsV] = m.trajectory.beginPoint().X();
            // TODO check if beginPoint is the correct point
            tI.m_XErrV[tI.m_NHitsV] = m.errMeasure;

            // z measure
            tI.m_ZMeasV[tI.m_NHitsV] = m.trajectory.beginPoint().Z();

            // y measure
            ++itM;
            itM->visit(
                [&]( const LHCb::Measurement::VP& n ) {
                  assert( n.projection() == LHCb::Measurement::VP::Projection::Y );
                  assert( itM->lhcbID() == lhcbID_X );
                  assert( itM->z() == z_X );
                  tI.m_YMeasV[tI.m_NHitsV] = n.trajectory.beginPoint().Y();
                  tI.m_YErrV[tI.m_NHitsV]  = n.errMeasure;
                },
                []( const auto& ) { throw std::logic_error( "VPX not followed by VPY" ); } );

            // set lhcbID
            tI.m_lhcbIDs[tI.m_NHitsTotal] = ( itM->lhcbID() );

            // increase counters
            tI.m_NHitsTotal++;
            tI.m_NHitsV++;
            lastVeloZ = m.trajectory.beginPoint().Z();
          },
          // Load VP2D hit__________________________________________
          [&]( const LHCb::Measurement::VP2D& m ) {
            // check if there was a previous hit in this layer
            if ( std::abs( lastVeloZ - m.trajectory.Z() ) < 2 ) { return; }

            // x and y measures
            tI.m_XMeasV[tI.m_NHitsV] = m.trajectory.X();
            tI.m_YMeasV[tI.m_NHitsV] = m.trajectory.Y();
            // TODO check if beginPoint is the correct point
            tI.m_XErrV[tI.m_NHitsV] = m.errMeasure[0];
            tI.m_YErrV[tI.m_NHitsV] = m.errMeasure[1];

            // z measure
            tI.m_ZMeasV[tI.m_NHitsV] = m.trajectory.Z();

            // set lhcbID
            tI.m_lhcbIDs[tI.m_NHitsTotal] = ( itM->lhcbID() );

            // increase counters
            tI.m_NHitsTotal++;
            tI.m_NHitsV++;
            lastVeloZ = m.trajectory.Z();
          },
          // Load UT hit__________________________________________
          [&]( const LHCb::Measurement::UT& m ) {
            if ( !m_UseUT ) return;

            int layerNum = ( itM->lhcbID().utID().station() - 1 ) * 2 + itM->lhcbID().utID().layer() - 1;
            // check if there was a previous hit in this layer
            if ( tI.m_HasHitUT[layerNum] == 1 ) { return; }

            tI.m_measUtLite[layerNum].point     = m.trajectory.beginPoint();
            tI.m_measUtLite[layerNum].direction = m.trajectory.direction( 0 );
            tI.m_measUtLite[layerNum].error     = m.errMeasure;

            // set lhcbID
            tI.m_lhcbIDs[tI.m_NHitsTotal] = itM->lhcbID();

            // increase counters
            tI.m_NHitsTotal++;
            tI.m_NHitsUT++;
            tI.m_HasHitUT[layerNum] = 1;
          },
          // Load FT hit__________________________________________
          [&]( const LHCb::Measurement::FT& m ) {
            if ( !m_UseT ) return;

            int layerNum =
                ( to_unsigned( itM->lhcbID().ftID().station() ) - 1 ) * 4 + to_unsigned( itM->lhcbID().ftID().layer() );
            // check if there was a previous hit in this layer
            if ( tI.m_HasHitT[layerNum] == 1 ) { return; }

            tI.m_measFtLite[layerNum].point     = m.trajectory.beginPoint();
            tI.m_measFtLite[layerNum].direction = m.trajectory.direction( 0 );
            tI.m_measFtLite[layerNum].error     = m.errMeasure;

            // set lhcbID
            tI.m_lhcbIDs[tI.m_NHitsTotal] = itM->lhcbID();

            // increase counters
            tI.m_NHitsTotal++;
            tI.m_NHitsT++;
            tI.m_HasHitT[layerNum] = 1;
          },
          []( const auto& ) { return; } );
      ++itM;
    } // end of itM loop
    // set ndof
    tI.m_Ndof   = tI.m_NHitsT + tI.m_NHitsUT + tI.m_NHitsV * 2;
    tI.m_NdofT  = tI.m_NHitsT;
    tI.m_NdofUT = tI.m_NHitsUT;
    tI.m_NdofV  = 2 * tI.m_NHitsV;

    // Set the hit status of all hits to 1: Should be used
    tI.m_HitStatus.fill( 1 );
    // create maps for the nHit -> layer matching
    // T stations
    int counter = 0;
    for ( int i = 0; i < 12; i++ ) {
      if ( tI.m_HasHitT[i] ) {
        tI.m_THitToTLayer[counter] = i;
        counter++;
      }
    }
    // UT
    counter = 0;
    for ( int i = 0; i < 4; i++ ) {
      if ( tI.m_HasHitUT[i] ) {
        tI.m_UTHitToUTLayer[counter] = i;
        counter++;
      }
    }
  }

  ///////////////////////////////////////////
  // Method to create a seed state at the first Velo hit
  ///////////////////////////////////////////
  void CreateVeloSeedState( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // Set hit position as x,y
    x( 0 ) = tI.m_XMeasV[nHit];
    x( 1 ) = tI.m_YMeasV[nHit];
    // Set slope between first and last hit as tx,ty
    x( 2 ) = ( tI.m_XMeasV[tI.m_NHitsV - 1] - tI.m_XMeasV[0] ) / ( tI.m_ZMeasV[tI.m_NHitsV - 1] - tI.m_ZMeasV[0] );
    x( 3 ) = ( tI.m_YMeasV[tI.m_NHitsV - 1] - tI.m_YMeasV[0] ) / ( tI.m_ZMeasV[tI.m_NHitsV - 1] - tI.m_ZMeasV[0] );
    // Set forward algorithm qop estimate as qop
    tI.m_BestMomEst = tI.m_track->firstState().qOverP();
    x( 4 )          = tI.m_BestMomEst;
    // Set hit z position as current z value
    lastz = tI.m_ZMeasV[nHit];

    // Set covariance matrix
    // Large uncertainties and no correlation
    C( 0, 0 ) = 100;
    C( 0, 1 ) = 0;
    C( 0, 2 ) = 0;
    C( 0, 3 ) = 0;
    C( 0, 4 ) = 0;
    C( 1, 1 ) = 100;
    C( 1, 2 ) = 0;
    C( 1, 3 ) = 0;
    C( 1, 4 ) = 0;
    C( 2, 2 ) = 1;
    C( 2, 3 ) = 0;
    C( 2, 4 ) = 0;
    C( 3, 3 ) = 1;
    C( 3, 4 ) = 0;
    C( 4, 4 ) = 0.09 * x( 4 ) * x( 4 );

    // Save information for the smoother
    if ( tI.m_do_smoother ) {
      tI.m_StateForwardPredicted[0] = x;
      tI.m_CovForwardPredicted[0]   = C;
      tI.m_StateZPos[0]             = lastz;
    }
  }

  //////////////////////////////////////////
  // General method for updating at a hit
  //////////////////////////////////////////
  void UpdateState( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // should the hit be used or not
    if ( tI.m_HitStatus[nHit] == 0 ) return;

    // Choose the appropiate updating method depending on the detector
    if ( nHit < tI.m_NHitsV )
      UpdateStateV( forward, nHit, x, C, tI );
    else if ( nHit < tI.m_NHitsV + tI.m_NHitsUT )
      UpdateStateUT( nHit - tI.m_NHitsV, x, C, lastz, tI );
    else
      UpdateStateT( forward, nHit - tI.m_NHitsV - tI.m_NHitsUT, x, C, lastz, tI );

    // Save information for the smoother
    if ( tI.m_do_smoother ) {
      if ( forward > 0 ) {
        tI.m_StateForwardUpdated[nHit] = x;
        tI.m_CovForwardUpdated[nHit]   = C;
      } else {
        tI.m_StateBackwardUpdated[nHit] = x;
        tI.m_CovBackwardUpdated[nHit]   = C;
      }
      tI.m_StateZPos[nHit] = lastz;
    }
  }

  //////////////////////////////////////////
  // Predict inside the VELO
  //////////////////////////////////////////
  void PredictStateV( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // Transporation and noise matrix
    Gaudi::Matrix5x5    F;
    Gaudi::SymMatrix5x5 Q;
    tI.m_extr->ExtrapolateInV( lastz, tI.m_ZMeasV[nHit], x, F, Q );

    // transport covariance matrix
    C = LHCb::Math::Similarity( F, C );

    // Add noise to covariance
    C += Q;

    // Set current z position
    lastz = tI.m_ZMeasV[nHit];
  }

  //////////////////////////////////////////
  // Predict VELO <-> UT
  //////////////////////////////////////////
  bool PredictStateVUT( Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    int forward = lastz < 1000. ? 1 : -1;
    // check prediction
    bool returnVal = true;
    // save old state variables
    Gaudi::Vector5 x_old = x;

    // predicted z position
    double z;
    // noise matrix
    Gaudi::SymMatrix5x5 Q;
    // jacobian matrix
    Gaudi::Matrix5x5 F;

    // Velo to UT prediction
    if ( forward > 0 ) {
      // save reference state for this intermediate extrapolation
      tI.m_RefStateForwardV = x;

      // Transporation and noise matrix
      Gaudi::Matrix5x5    F;
      Gaudi::SymMatrix5x5 Q;

      if ( tI.m_HasHitUT[0] != 0 ) {
        returnVal =
            tI.m_extr->ExtrapolateVUT( lastz, tI.m_measUtLite[0].point, tI.m_measUtLite[0].direction, z, x, F, Q );
      } else {
        z         = tI.m_extr->VUTExtrEndZ();
        returnVal = tI.m_extr->ExtrapolateVUT( lastz, z, x, F, Q );
      }

      // transport covariance matrix
      C = LHCb::Math::Similarity( F, C );

      // Add noise to covariance
      C += Q;

      // save the extrapolation for future usage as reference
      tI.m_RefPropForwardVUT  = F;
      tI.m_RefStateForwardFUT = x;
    }
    // T to Velo prediction (using forward as reference)
    else if ( true ) {
      // z position of last Velo hit
      z = tI.m_ZMeasV[tI.m_NHitsV - 1];

      // propagate deviation from reference (prediction from V to T)
      // propagation matrix from forward filtering
      // x=x_prev_ref+F^{-1}*(x-xPred_ref)
      F = tI.m_RefPropForwardVUT;
      F.InvertFast();

      x = tI.m_RefStateForwardV + F * ( x - tI.m_RefStateForwardFUT );

      // transport covariance matrix
      C = LHCb::Math::Similarity( F, C );

      // add noise
      Gaudi::SymMatrix5x5 Q;
      tI.m_extr->GetNoiseVUTBackw( lastz, z, x_old, Q );

      // Add noise to covariance
      C += Q;
    }

    // Velo to T prediction without reference
    else {
      // z position of last Velo hit
      z = tI.m_ZMeasV[tI.m_NHitsV - 1];

      // Transporation and noise matrix
      Gaudi::Matrix5x5    F;
      Gaudi::SymMatrix5x5 Q;

      returnVal = tI.m_extr->ExtrapolateVUT( lastz, z, x, F, Q );

      // transport covariance matrix
      C = LHCb::Math::Similarity( F, C );

      // Add noise to covariance
      C += Q;
    }

    // set current z position
    lastz = z;

    return returnVal;
  }

  //////////////////////////////////////////
  // Predict UT <-> UT
  //////////////////////////////////////////
  void PredictStateUT( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // new z position
    double z;
    // Transporation and noise matrix
    Gaudi::Matrix5x5    F;
    Gaudi::SymMatrix5x5 Q;

    if ( tI.m_HasHitUT[nHit] != 0 ) {
      tI.m_extr->ExtrapolateInUT( lastz, nHit, tI.m_measUtLite[nHit].point, tI.m_measUtLite[nHit].direction, z, x, F,
                                  Q );
    } else {
      z = lastz;
      tI.m_extr->ExtrapolateInUT( lastz, nHit, z, x, F, Q );
    }

    // transport covariance matrix
    C = LHCb::Math::Similarity( F, C );

    // Add noise to covariance
    C += Q;
    // Reminder: backward T station label (different than for tuples)

    lastz = z;
  }

  //////////////////////////////////////////
  // Predict UT <-> T precise version
  //////////////////////////////////////////
  void PredictStateUTT( Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    int forward = lastz < 5000. ? 1 : -1;

    // cache old state
    Gaudi::Vector5 x_old = x;

    // when going backward: predict to fixed z in T(z=7855)
    if ( forward < 0 ) {
      PredictStateTFT( forward, x, C, lastz, tI );
    }
    // If we are at a different z position: go to the start position of the extrpolation)
    else if ( tI.m_extr->UTTExtrBeginZ() != lastz ) {
      PredictStateUTFUT( forward, x, C, lastz, tI );
    }

    // jacobian matrix
    Gaudi::Matrix5x5 F;

    // extrapolating from last UT hit (if no hit: z is set to a default value)
    // to fixed z in T (z position is a parameters set to the middel of the first layer)
    if ( forward > 0 ) {
      // If we are at a different z position: go to the start position of the extrpolation)
      if ( tI.m_extr->UTTExtrBeginZ() != lastz ) {}

      // Calculate the extrapolation for a refernece state that uses the
      // inital forward momentum estimate
      // cache old state
      Gaudi::Vector5 xref = x;
      // This reference can then also be used for the backward propagation
      // This gives a better momentum estimate
      xref[4] = tI.m_BestMomEst;
      // save reference state for this intermediate extrapolation
      tI.m_RefStateForwardUT = xref;

      // Transporation and noise matrix
      Gaudi::Matrix5x5    F;
      Gaudi::SymMatrix5x5 Q;
      tI.m_extr->ExtrapolateUTT( xref, F, Q );

      // save reference state/jacobian after this intermediate extrapolation
      tI.m_RefStateForwardT  = xref;
      tI.m_RefPropForwardUTT = F;

      // extrapolate the actual state
      // propagate the deviation from reference
      x = tI.m_RefStateForwardT + F * ( x - tI.m_RefStateForwardUT );

      // Set current z position
      lastz = tI.m_extr->UTTExtrEndZ();

      // transport covariance matrix
      C = LHCb::Math::Similarity( F, C );

      // Add noise to covariance
      C += Q;

    }

    //(no parametrization for this -> use reference)
    else {
      // propagate deviation from reference (use forward filtering for this)
      // use inverted jacobian from forward extrapolation TODO use that der_x, der_y=0
      F = tI.m_RefPropForwardUTT;
      F.InvertFast();

      // propagate the eviation from reference
      x = tI.m_RefStateForwardUT + F * ( x_old - tI.m_RefStateForwardT );

      // set current z position
      lastz = tI.m_extr->UTTExtrBeginZ();

      // transport covariance matrix
      C = LHCb::Math::Similarity( F, C );

      // Get noise
      Gaudi::SymMatrix5x5 Q;
      tI.m_extr->GetNoiseUTTBackw( x_old, Q );

      // Add noise to covariance
      C += Q;
    }

    // When going backwards: predict to the last VELO measurement
    if ( forward > 0 ) {
      PredictStateTFT( forward, x, C, lastz, tI );
    }
    // in case of a hit, z might not be exactly the default position:
    else if ( tI.m_HasHitUT[3] == 1 ) {
      PredictStateUTFUT( forward, x, C, lastz, tI );
    }
  }

  //////////////////////////////////////////////
  // Predict UT(fixed z) <-> last UT layer
  //////////////////////////////////////////////
  void PredictStateUTFUT( int forward, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // Transporation and noise matrix
    Gaudi::Matrix5x5 F;

    if ( forward > 0 ) {
      tI.m_extr->ExtrapolateUTFUTDef( lastz, x, F );
    } else {
      tI.m_extr->ExtrapolateUTFUT( lastz, tI.m_measUtLite[3].point, tI.m_measUtLite[3].direction, x, F );
    }

    // transport covariance matrix
    C = LHCb::Math::Similarity( F, C );
  }

  //////////////////////////////////////////
  // Predict T <-> T
  //////////////////////////////////////////
  void PredictStateT( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // new z position
    double z;
    // Transporation and noise matrix
    Gaudi::Matrix5x5    F;
    Gaudi::SymMatrix5x5 Q;
    if ( tI.m_HasHitT[nHit] != 0 )
      tI.m_extr->ExtrapolateInT( lastz, nHit, tI.m_measFtLite[nHit].point, tI.m_measFtLite[nHit].direction, z, x, F,
                                 Q );
    else {
      tI.m_extr->ExtrapolateInT( lastz, nHit, z, x, F, Q );
    }

    // transport covariance matrix
    C = LHCb::Math::Similarity( F, C );

    // Add noise to covariance
    C += Q;

    lastz = z;

    // Check if this is needed
    // lastz = z0+1.0/dydz*(x[1]-y0);//TODO: try to improve the prediction to a non-fix z plane
    //(do iterations) (predict first y,ty at central value and determine then z)
    // more important for the step through the magnet.
    // Maybe predicting first to a vertical plane is usefull
  }

  //////////////////////////////////////////////
  // Predict T(fixed z) <-> first T layer
  //////////////////////////////////////////////
  void PredictStateTFT( int forward, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // new z position
    double z;
    // Transporation and noise matrix
    Gaudi::Matrix5x5    F;
    Gaudi::SymMatrix5x5 Q;

    if ( forward > 0 ) {
      if ( tI.m_HasHitT[0] != 0 ) {
        tI.m_extr->ExtrapolateTFT( lastz, tI.m_measFtLite[0].point, tI.m_measFtLite[0].direction, z, x, F, Q );
      } else {
        tI.m_extr->ExtrapolateTFTDef( lastz, z, x, F, Q );
      }
    } else {
      z = tI.m_extr->UTTExtrEndZ();
      tI.m_extr->ExtrapolateTFT( lastz, z, x, F, Q );
    }

    // transport covariance matrix
    C = LHCb::Math::Similarity( F, C );

    // Add noise to covariance
    C += Q;

    lastz = z;
  }

  /////////////////////////////////////////
  // Update state with velo measurement
  /////////////////////////////////////////
  void UpdateStateV( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, trackInfo& tI ) {
    // set measurement residual
    Gaudi::Vector2 res( tI.m_XMeasV[nHit] - x( 0 ), tI.m_YMeasV[nHit] - x( 1 ) );

    double              CresTmp[3] = {tI.m_XErrV[nHit] * tI.m_XErrV[nHit] + C( 0, 0 ), C( 0, 1 ),
                         tI.m_YErrV[nHit] * tI.m_YErrV[nHit] + C( 1, 1 )};
    Gaudi::SymMatrix2x2 CRes( CresTmp, 3 );

    // Kalman formalism:
    int                 fail;
    Gaudi::SymMatrix2x2 CResInv = CRes.Inverse( fail );

    ROOT::Math::SMatrix<double, 5, 2> K;

    Multiply_S5x5_S2x2( C, CResInv, K );

    x += K * res;

    Gaudi::SymMatrix5x5 KCrKt;
    LHCb::Math::Similarity( K, CRes, KCrKt );

    C -= KCrKt;

    double chi2Tmp = ROOT::Math::Similarity( res, CResInv );
    tI.m_chi2 += chi2Tmp;
    if ( forward > 0 ) tI.m_chi2V += chi2Tmp;
  }

  /////////////////////////////////////////
  // Update state with UT measurement
  /////////////////////////////////////////
  void UpdateStateUT( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // get the layer in which this UT-station is located
    int layer = tI.m_UTHitToUTLayer[nHit];
    // Projection matrix
    // tan(apha)=dx/dy

    auto point = tI.m_measUtLite[layer].point;
    auto dir   = tI.m_measUtLite[layer].direction;

    // project to a coordinate system that is rotated by
    // alpha= atan2(traj.direction(0).X(),traj.direction(0).Y());
    const double         x2y2 = std::sqrt( dir.X() * dir.X() + dir.Y() * dir.Y() );
    const Gaudi::Vector2 H( dir.Y() / x2y2, -dir.X() / x2y2 );

    // Residual
    // take the begin point of trajectory and rotate
    const double res = H[0] * point.X() + H[1] * point.Y() - ( H[0] * x[0] + H[1] * x[1] );
    double       CRes;
    Similarity_1x2_S5x5_2x1( H, C, CRes );
    CRes += tI.m_measUtLite[layer].error * tI.m_measUtLite[layer].error;

    // K = P*H
    Gaudi::Vector5 K;
    Multiply_S5x5_2x1( C, H, K );
    // K * S^-1
    K /= CRes;
    // x += K*res
    x += K * res;

    // K*S*K(T) (is the same as K*H*P)
    Gaudi::SymMatrix5x5 KCResKt;
    SymmetricTensorProduct5( std::sqrt( CRes ) * K, KCResKt );
    // P -= KSK
    C -= KCResKt;

    // calculate chi2
    tI.m_chi2 += res * res / CRes;

    // TODO does it help to reset the z position here?
    lastz = point.Z() + dir.Z() / dir.Y() * ( x[1] - point.Y() );
  }

  /////////////////////////////////////////
  // Update state with T measurement
  /////////////////////////////////////////
  void UpdateStateT( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI ) {
    // get the layer in which this T-station is located
    int layer = tI.m_THitToTLayer[nHit];
    // Projection matrix
    // tan(apha)=dx/dy

    auto point = tI.m_measFtLite[layer].point;
    auto dir   = tI.m_measFtLite[layer].direction;

    // project to a coordinate system that is rotated by
    // alpha= atan2(traj.direction(0).X(),traj.direction(0).Y());
    const double         x2y2 = std::sqrt( dir.X() * dir.X() + dir.Y() * dir.Y() );
    const Gaudi::Vector2 H( dir.Y() / x2y2, -dir.X() / x2y2 );

    // Residual
    // take the begin point of trajectory and rotate
    const double res = H[0] * point.X() + H[1] * point.Y() - ( H[0] * x[0] + H[1] * x[1] );
    double       CRes;
    Similarity_1x2_S5x5_2x1( H, C, CRes );
    CRes += tI.m_measFtLite[layer].error * tI.m_measFtLite[layer].error;

    // K = P*H
    Gaudi::Vector5 K;
    Multiply_S5x5_2x1( C, H, K );
    // K * S^-1
    K /= CRes;
    // x += K*ires
    x += K * res;

    // K*S*K(T) (is the same as K*H*P)
    Gaudi::SymMatrix5x5 KCResKt;
    SymmetricTensorProduct5( std::sqrt( CRes ) * K, KCResKt );
    // P -= KSK
    C -= KCResKt;

    // calculate chi2
    tI.m_chi2 += res * res / CRes;

    if ( forward < 0 ) tI.m_chi2T += res * res / CRes;

    // TODO does it help to reset the z position here?
    lastz = point.Z() + dir.Z() / dir.Y() * ( x[1] - point.Y() );
  }

  ///////////////////////////////////////////////
  //  Smooth/average method
  ///////////////////////////////////////////////
  bool AverageState( int nHit, trackInfo& tI ) {
    // is this hit used?
    if ( tI.m_HitStatus[nHit] == 0 ) return true;
    // calculate the average state
    bool averaged = LHCb::Math::Average( tI.m_StateBackwardUpdated[nHit], tI.m_CovBackwardUpdated[nHit],
                                         tI.m_StateForwardPredicted[nHit], tI.m_CovForwardPredicted[nHit],
                                         tI.m_StateSmoothed[nHit], tI.m_CovSmoothed[nHit] );

    // chi2 calculation
    tI.m_HitChi2[nHit] = 0;
    // For T-stations
    if ( nHit >= tI.m_NHitsV + tI.m_NHitsUT ) {
      // get the layer in which this T-station is located
      int layer = tI.m_THitToTLayer[nHit - tI.m_NHitsV - tI.m_NHitsUT];
      // Projection matrix
      // tan(apha)=dx/dy

      auto point = tI.m_measFtLite[layer].point;
      auto dir   = tI.m_measFtLite[layer].direction;
      auto error = tI.m_measFtLite[layer].error;

      // project to a coordinate system that is rotated by
      // alpha= atan2(dir.X(),dir.Y());
      const double         x2y2 = std::sqrt( dir.X() * dir.X() + dir.Y() * dir.Y() );
      const Gaudi::Vector2 H( dir.Y() / x2y2, -dir.X() / x2y2 );

      // Residual
      // take the begin point of trajectory and rotate
      const double res = H[0] * point.X() + H[1] * point.Y() -
                         ( H[0] * tI.m_StateSmoothed[nHit][0] + H[1] * tI.m_StateSmoothed[nHit][1] );

      double CRes;
      Similarity_1x2_S5x5_2x1( H, tI.m_CovSmoothed[nHit], CRes );
      CRes -= error * error;

      tI.m_HitChi2[nHit] = -res * res / CRes;
    } else if ( nHit >= tI.m_NHitsV ) {
      // get the layer in which this UT-station is located
      int layer = tI.m_UTHitToUTLayer[nHit - tI.m_NHitsV];
      // Projection matrix
      // tan(apha)=dx/dy
      auto point = tI.m_measUtLite[layer].point;
      auto dir   = tI.m_measUtLite[layer].direction;
      auto error = tI.m_measUtLite[layer].error;

      // project to a coordinate system that is rotated by
      // alpha= atan2(dir.X(),dir.Y());
      const double         x2y2 = std::sqrt( dir.X() * dir.X() + dir.Y() * dir.Y() );
      const Gaudi::Vector2 H( dir.Y() / x2y2, -dir.X() / x2y2 );

      // Residual
      // take the begin point of trajectory and rotate
      const double res = H[0] * point.X() + H[1] * point.Y() -
                         ( H[0] * tI.m_StateSmoothed[nHit][0] + H[1] * tI.m_StateSmoothed[nHit][1] );

      double CRes;
      Similarity_1x2_S5x5_2x1( H, tI.m_CovSmoothed[nHit], CRes );
      CRes -= error * error;

      tI.m_HitChi2[nHit] = -res * res / CRes;
    }
    // For pixels in the VELO
    else {
      // set measurement residual
      Gaudi::Vector2 res( tI.m_XMeasV[nHit] - tI.m_StateSmoothed[nHit][0],
                          tI.m_YMeasV[nHit] - tI.m_StateSmoothed[nHit][1] );

      double              CresTmp[3] = {tI.m_XErrV[nHit] * tI.m_XErrV[nHit] - tI.m_CovSmoothed[nHit]( 0, 0 ),
                           -tI.m_CovSmoothed[nHit]( 0, 1 ),
                           tI.m_YErrV[nHit] * tI.m_YErrV[nHit] - tI.m_CovSmoothed[nHit]( 1, 1 )};
      Gaudi::SymMatrix2x2 CRes( CresTmp, 3 );

      int fail;
      tI.m_HitChi2[nHit] = ROOT::Math::Similarity( res, CRes.Inverse( fail ) );
    }
    tI.m_chi2 += tI.m_HitChi2[nHit];
    return averaged;
  }

  ////////////////////////////////////////////////////////////////
  // extrapolate to the vertex using the default extrpolator
  ////////////////////////////////////////////////////////////////
  void ExtrapolateToVertex( Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, IGeometryInfo const& geometry,
                            ToolHandle<ITrackExtrapolator> const& m_extrapolator_toPV ) {
    LHCb::State state( x, C, lastz, LHCb::State::Location::LocationUnknown );

    // determine the z position to extrapolate to
    double z = lastz - ( x[0] * x[2] + x[1] * x[3] ) / ( x[2] * x[2] + x[3] * x[3] );

    // get the extrapolated state
    m_extrapolator_toPV->propagate( state, z, geometry ).ignore();
    x[0]  = state.x();
    x[1]  = state.y();
    x[2]  = state.tx();
    x[3]  = state.ty();
    x[4]  = state.qOverP();
    C     = state.covariance();
    lastz = z;
  }

  ///////////////////////////////////////////////////////////////////////////
  // Create Track
  ///////////////////////////////////////////////////////////////////////////
  void addInfoToTrack( const Gaudi::Vector5& x, const Gaudi::SymMatrix5x5& C, double z, trackInfo& tI ) {
    // delete all states
    tI.m_track->clearStates();
    // add the sate(should be the state at the beam)
    LHCb::State state( x, C, z, LHCb::State::Location::ClosestToBeam );
    // set track chi2 and ndof
    tI.m_track->setChi2PerDoF( tI.m_chi2 / ( tI.m_Ndof - 5 ) );
    tI.m_track->setNDoF( tI.m_Ndof - 5 );

    // set additional information for later usage
    tI.m_track->eraseInfo( OutputTrack::AdditionalInfo::FitVeloChi2 );
    tI.m_track->eraseInfo( OutputTrack::AdditionalInfo::FitVeloNDoF );
    tI.m_track->eraseInfo( OutputTrack::AdditionalInfo::FitTChi2 );
    tI.m_track->eraseInfo( OutputTrack::AdditionalInfo::FitTNDoF );

    tI.m_track->eraseInfo( OutputTrack::AdditionalInfo::FitMatchChi2 );

    tI.m_track->addInfo( OutputTrack::AdditionalInfo::FitVeloChi2, tI.m_chi2V );
    tI.m_track->addInfo( OutputTrack::AdditionalInfo::FitVeloNDoF, tI.m_NdofV - 5 );
    tI.m_track->addInfo( OutputTrack::AdditionalInfo::FitTChi2, tI.m_chi2T );
    tI.m_track->addInfo( OutputTrack::AdditionalInfo::FitTNDoF, tI.m_NdofT - 5 );

    tI.m_track->addInfo( OutputTrack::AdditionalInfo::FitMatchChi2, tI.m_chi2 - tI.m_chi2T - tI.m_chi2V );

    // add lhcb ids. Let all lhcb ids untouched for the moment
    // TrackIDContainer IDContainer;
    // for(int i = 0; i<tI.m_NHitsTotal; i++){
    //  if(m_HitStatus[i]){
    //    IDContainer.push_back(m_lhcbIDs[i]);
    //  }
    //}
    // track->setLhcbIDs(IDContainer);

    tI.m_track->addToStates( state );
  }

  ///////////////////////////////////////////////////////////////////////////
  // Add state to track
  ///////////////////////////////////////////////////////////////////////////
  void addStateToTrack( const Gaudi::Vector5& x, const Gaudi::SymMatrix5x5& C, double z, trackInfo& tI,
                        LHCb::State::Location loc ) {
    // add the state
    LHCb::State state( x, C, z, loc );
    tI.m_track->addToStates( state );
  }

  ///////////////////////////////////////////////////////////////////////////
  // Check if outliers should be removed and remove one of them
  ///////////////////////////////////////////////////////////////////////////
  bool DoOutlierRemoval( trackInfo& tI ) {
    double maxChi2 = 0;
    int    n       = 0;
    for ( int i = 0; i < tI.m_NHitsTotal; i++ ) {
      if ( tI.m_HitChi2[i] * tI.m_HitStatus[i] > maxChi2 ) {
        maxChi2 = tI.m_HitChi2[i];
        n       = i;
      }
    }
    if ( maxChi2 > 9 ) {
      tI.m_HitStatus[n] = 0;
      // update degrees of freedom
      if ( n < tI.m_NHitsV ) {
        tI.m_Ndof  = tI.m_Ndof - 2;
        tI.m_NdofV = tI.m_NdofV - 2;
      } else if ( n < tI.m_NHitsV + tI.m_NHitsUT ) {
        tI.m_Ndof   = tI.m_Ndof - 1;
        tI.m_NdofUT = tI.m_NdofUT - 1;
      } else {
        tI.m_Ndof -= 1;
        tI.m_NdofT -= 1;
      }
      return true;
    } else
      return false;
  }

  void Similarity_1x2_S5x5_2x1( const Gaudi::Vector2& A, const Gaudi::SymMatrix5x5& B, double& R ) {
    R = A( 0 ) * ( A( 0 ) * B( 0, 0 ) + A( 1 ) * B( 0, 1 ) ) + A( 1 ) * ( A( 0 ) * B( 0, 1 ) + A( 1 ) * B( 1, 1 ) );
  }

  void SymmetricTensorProduct5( const Gaudi::Vector5 A, Gaudi::SymMatrix5x5& RM ) {
    double* R = RM.Array();
    R[0]      = A[0] * A[0];
    R[1]      = A[1] * A[0];
    R[2]      = A[1] * A[1];
    R[3]      = A[2] * A[0];
    R[4]      = A[2] * A[1];
    R[5]      = A[2] * A[2];
    R[6]      = A[3] * A[0];
    R[7]      = A[3] * A[1];
    R[8]      = A[3] * A[2];
    R[9]      = A[3] * A[3];
    R[10]     = A[4] * A[0];
    R[11]     = A[4] * A[1];
    R[12]     = A[4] * A[2];
    R[13]     = A[4] * A[3];
    R[14]     = A[4] * A[4];
  }

  void Multiply_S5x5_2x1( const Gaudi::SymMatrix5x5& AM, const Gaudi::Vector2& B, Gaudi::Vector5& R ) {
    const double* A = AM.Array();
    R[0]            = A[0] * B[0] + A[1] * B[1];
    R[1]            = A[1] * B[0] + A[2] * B[1];
    R[2]            = A[3] * B[0] + A[4] * B[1];
    R[3]            = A[6] * B[0] + A[7] * B[1];
    R[4]            = A[10] * B[0] + A[11] * B[1];
  }

  void Multiply_S5x5_S2x2( const Gaudi::SymMatrix5x5& AM, const Gaudi::SymMatrix2x2& BM,
                           ROOT::Math::SMatrix<double, 5, 2>& RM ) {
    const double* A = AM.Array();
    const double* B = BM.Array();
    double*       R = RM.Array();
    R[0]            = A[0] * B[0] + A[1] * B[1];
    R[1]            = A[0] * B[1] + A[1] * B[2];

    R[2] = A[1] * B[0] + A[2] * B[1];
    R[3] = A[1] * B[1] + A[2] * B[2];

    R[4] = A[3] * B[0] + A[4] * B[1];
    R[5] = A[3] * B[1] + A[4] * B[2];

    R[6] = A[6] * B[0] + A[7] * B[1];
    R[7] = A[6] * B[1] + A[7] * B[2];

    R[8] = A[10] * B[0] + A[11] * B[1];
    R[9] = A[10] * B[1] + A[11] * B[2];
  }
} // namespace ParKalman
