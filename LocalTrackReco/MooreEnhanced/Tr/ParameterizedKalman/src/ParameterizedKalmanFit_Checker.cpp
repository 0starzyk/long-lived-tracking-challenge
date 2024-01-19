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

#include "Event/State.h"

#include "Event/MCParticle.h"

#include "Event/ODIN.h"

#include "Event/MCTrackInfo.h"

#include "LHCbMath/Similarity.h"

#include "Kernel/Trajectory.h"

// local
#include "ParameterizedKalmanFit_Checker.h"

using namespace ParKalman;

// position of the individual trees in the tree vector
namespace {
  enum TrPos {
    TrPos_comp      = 0,
    TrPos_crSeed    = 1,
    TrPos_predV     = 2,
    TrPos_predVUT   = 4,
    TrPos_predUT    = 6,
    TrPos_predUTFUT = 14,
    TrPos_predUTTF  = 15,
    TrPos_predUTT   = 17,
    TrPos_predTFT   = 21,
    TrPos_predT     = 23,
    TrPos_upV       = 71,
    TrPos_upLV      = 73,
    TrPos_upFUT     = 74,
    TrPos_upUT      = 75,
    TrPos_upLUT     = 83,
    TrPos_upFT      = 84,
    TrPos_upT       = 85
  };
}

//########################################################################
//
// Implementation file for class : ParameterizedKalmanFit_Checker
//
// 2017-11-02: Simon Stemmle
//
//########################################################################

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( ParameterizedKalmanFit_Checker )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
ParameterizedKalmanFit_Checker::ParameterizedKalmanFit_Checker( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator,
                   {KeyValue{"InputName", "Rec/Track/ForwardFast"},
                    KeyValue{"ODINLocation", LHCb::ODINLocation::Default},
                    KeyValue{"LinkerLocation", Links::location( "Rec/Track/ForwardFast" )},
                    KeyValue{"MCProperty", LHCb::MCPropertyLocation::TrackInfo},
                    KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top},
                    KeyValue{"Magnet", LHCb::Det::Magnet::det_path}},
                   KeyValue{"OutputName", "Rec/Track/ForwardFastFitted_Checker"} ) {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode ParameterizedKalmanFit_Checker::initialize() {
  StatusCode sc = Transformer::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;           // error printed already by GaudiAlgorithm

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Initialize" << endmsg;

  info() << "Use parameters from " << m_ParamFileLocation.value() + "/Mag*" << endmsg;

  m_ParExtrUp   = std::make_unique<KalmanParametrizations>( m_ParamFileLocation, Polarity::Up, m_UseOneParameterSet );
  m_ParExtrDown = std::make_unique<KalmanParametrizations>( m_ParamFileLocation, Polarity::Down, m_UseOneParameterSet );

  // set state to truth at certain step in the kf when we run it for the tuning
  // keep it configurable for debugging
  if ( m_RunForTuning ) {
    m_SetTrueStateAfterUpdate     = true;
    m_SetTrueStateAfterPredict    = true;
    m_SetTrueStateAfterCreateSeed = true;
  }

  // cache information for the smoother of outliers should be removed
  m_do_smoother = m_MaxNoutlier > 0;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
LHCb::Tracks ParameterizedKalmanFit_Checker::operator()( LHCb::Tracks const& input, LHCb::ODIN const& odin,
                                                         LHCb::LinksByKey const& links, LHCb::MCProperty const& mcprop,
                                                         DetectorElement const& lhcb, DeMagnet const& magnet ) const {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;
  //============================================================
  //== Main processing: fit Tracks
  //============================================================

  // output tracks
  LHCb::Tracks result;

  if ( input.empty() ) { return result; }

  // fit status
  StatusCode sc;

  // struct that contains the intermediate track information
  trackInfo tI{links};
  tI.m_do_smoother = m_do_smoother;

  const auto trackInfo = MCTrackInfo{mcprop};

  // select the respective extrapolator
  tI.m_extr = ( magnet.signedRelativeCurrent() > 0 ? m_ParExtrUp.get() : m_ParExtrDown.get() );

  // Create output tuples that contain information for tuning or performance tests
  // Create a new file for every event in order to be threadsafe
  tI.m_TreesFile = new TFile( ( m_TreesFileName.toString() + "_" + std::to_string( odin.runNumber() ) + "_" +
                                std::to_string( odin.eventNumber() ) + ".root" )
                                  .c_str(),
                              "RECREATE" );

  // create the varaibles to be filled
  trackTupleInfo treeVars;

  // create the trees
  std::vector<TTree*> trees;
  tI.m_TreesFile->cd();
  addTrees( trees, &treeVars );

  // do the tuning iterations
  if ( m_RunForTuning ) {
    // Loop over the tracks and fit them
    for ( auto const& trackIn : input ) {
      // Create a new track keeping the same key
      auto track = std::make_unique<LHCb::Track>( *trackIn, trackIn->key() );

      tI.m_track = track.get();

      fitForTuning( tI, trackInfo, &trees, &treeVars, *lhcb.geometry() )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      result.insert( track.release() );
    }
  }
  // do the default fit
  else {
    // Loop over the tracks and fit them
    for ( auto const& trackIn : input ) {
      // Create a new track keeping the same key
      auto track = std::make_unique<LHCb::Track>( *trackIn, trackIn->key() );

      tI.m_track = track.get();

      sc = fit( tI, trackInfo, &trees, &treeVars, *lhcb.geometry() );
      if ( sc.isSuccess() ) result.insert( track.release() );
    }
  }
  tI.m_TreesFile->Write();
  tI.m_TreesFile->Close();

  return result;
}

//=============================================================================
//  Run to extract tuning information
//=============================================================================
StatusCode ParameterizedKalmanFit_Checker::fitForTuning( trackInfo& tI, const MCTrackInfo& trackInfo,
                                                         std::vector<TTree*>* trees, trackTupleInfo* tV,
                                                         IGeometryInfo const& geometry ) const {

  // check if a matching mc particle is found
  tV->m_MC_status = MatchesMC( tI, trackInfo );
  if ( tV->m_MC_status != 1 ) return StatusCode::SUCCESS;

  // load hit information
  LoadHits_Ch( tI, m_measProvider, m_UseUT, m_UseT, tV );
  // use only tracks with hits in every FT, UT layer
  if ( tI.m_NHitsUT != 4 || tI.m_NHitsT != 12 ) return StatusCode::SUCCESS;
  // current state:
  // current z position
  double lastz = -1.;
  // state vector x: x,y,tx,ty,qop
  Gaudi::Vector5 x;
  // covariance
  Gaudi::SymMatrix5x5 C;

  //############################
  // do the actual Kalman filter
  //############################

  // create the seed state at the first VELO hit
  CreateVeloSeedState_Ch( 0, x, C, lastz, tI, trees, tV, geometry );
  // Update with the information of this hit
  UpdateState_Ch( 1, 0, x, C, lastz, tI, trees, tV, geometry );
  // forward filtering
  for ( int nhit = 1; nhit < tI.m_NHitsTotal; nhit++ ) {
    if ( !PredictState_Ch( 1, nhit, x, C, lastz, tI, trees, tV, geometry ) ) return StatusCode::FAILURE;
    UpdateState_Ch( 1, nhit, x, C, lastz, tI, trees, tV, geometry );
  }

  // backward filtering
  // start with an update using the information of the last hit
  UpdateState_Ch( -1, tI.m_NHitsTotal - 1, x, C, lastz, tI, trees, tV, geometry );
  for ( int nhit = tI.m_NHitsTotal - 2; nhit >= 0; nhit-- ) {
    if ( !PredictState_Ch( -1, nhit, x, C, lastz, tI, trees, tV, geometry ) ) return StatusCode::FAILURE;
    UpdateState_Ch( -1, nhit, x, C, lastz, tI, trees, tV, geometry );
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
// Do the default fit
//=============================================================================
StatusCode ParameterizedKalmanFit_Checker::fit( trackInfo& tI, MCTrackInfo const& trackInfo, std::vector<TTree*>* trees,
                                                trackTupleInfo* tV, IGeometryInfo const& geometry ) const {
  // check if the correct matching particle or no particle is found
  tV->m_MC_status = MatchesMC( tI, trackInfo );
  if ( tV->m_MC_status == 2 ) return StatusCode::SUCCESS;

  // load hit information
  LoadHits_Ch( tI, m_measProvider, m_UseUT, m_UseT, tV );

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
  CreateVeloSeedState_Ch( 0, x, C, lastz, tI, trees, tV, geometry );

  // reset chi2 for the forward filtering
  tI.m_chi2  = 0;
  tI.m_chi2V = 0;

  // start by updating with the first measurment
  UpdateState_Ch( 1, 0, x, C, lastz, tI, trees, tV, geometry );

  // forward filtering
  for ( int nhit = 1; nhit < tI.m_NHitsTotal; nhit++ ) {
    if ( !PredictState_Ch( 1, nhit, x, C, lastz, tI, trees, tV, geometry ) ) return StatusCode::FAILURE;
    UpdateState_Ch( 1, nhit, x, C, lastz, tI, trees, tV, geometry );
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
  UpdateState_Ch( -1, tI.m_NHitsTotal - 1, x, C, lastz, tI, trees, tV, geometry );
  for ( int nhit = tI.m_NHitsTotal - 2; nhit >= 0; nhit-- ) {
    if ( !PredictState_Ch( -1, nhit, x, C, lastz, tI, trees, tV, geometry ) ) return StatusCode::FAILURE;
    UpdateState_Ch( -1, nhit, x, C, lastz, tI, trees, tV, geometry );
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
      UpdateState_Ch( 1, 0, x, C, lastz, tI, trees, tV, geometry );

      // forward filtering
      for ( int nhit = 1; nhit < tI.m_NHitsTotal; nhit++ ) {
        if ( !PredictState_Ch( 1, nhit, x, C, lastz, tI, trees, tV, geometry ) ) return StatusCode::FAILURE;
        UpdateState_Ch( 1, nhit, x, C, lastz, tI, trees, tV, geometry );
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
      UpdateState_Ch( -1, tI.m_NHitsTotal - 1, x, C, lastz, tI, trees, tV, geometry );
      for ( int nhit = tI.m_NHitsTotal - 2; nhit >= 0; nhit-- ) {
        if ( !PredictState_Ch( -1, nhit, x, C, lastz, tI, trees, tV, geometry ) ) return StatusCode::FAILURE;
        UpdateState_Ch( -1, nhit, x, C, lastz, tI, trees, tV, geometry );
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

  // Comapre to true state
  FillNtuple( xBest, CBest, zBest, tI, tV, -1000, 0, geometry );
  FillNtuple( xBestFirstMeas, CBestFirstMeas, zBestFirstMeas, tI, tV, -1000, 1, geometry );
  FillNtuple( xBestLastMeas, CBestLastMeas, zBestLastMeas, tI, tV, -1000, 2, geometry );

  ( *trees )[TrPos_comp]->Fill();

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Create trees that should be filled for tuning and perfomance checks
//=============================================================================
void ParameterizedKalmanFit_Checker::addTrees( std::vector<TTree*>& trees, trackTupleInfo* treeVars ) const {
  // Comparing the perfomance________________________________________
  // 0: Best states
  //
  // Seeding_________________________________________________________
  // 1: Create seed state in VELO
  //
  // Predicting______________________________________________________
  // 2-3:    Predict VELO                      <-> VELO
  // 4-5:    Predict last VELO measurement     <-> UT
  // 6-13:  Predict UT                        <-> UT
  // 14-15:  Predict last UT measurement       <-> T (fixed z)
  // 16-19:  Predict last UT measurement       <-> T
  // 20-21:  Predict T (fixed z)               <-> T
  // 22-59:  Predict T                         <-> T
  //
  // Updating________________________________________________________
  // 60-61:  Update inside     VELO
  // 62:     Update last       VELO coming from UT/T
  // 63:     Update first      UT   coming from VELO
  // 64-71:  Update inside     UT
  // 72:     Update last       UT   coming from T
  // 73:     Update first      T    coming from VELO/UT
  // 74-121: Update inside     T

  trees.push_back( new TTree( "compare", "compare" ) );

  trees.push_back( new TTree( "createSeedV_0", "createSeedV_0" ) );

  for ( int i = 0; i < 2; i++ )
    trees.push_back( new TTree( ( "predictStateV_" + std::to_string( i ) ).c_str(),
                                ( "PredictStateV_" + std::to_string( i ) ).c_str() ) );

  for ( int i = 0; i < 2; i++ )
    trees.push_back( new TTree( ( "predictStateVUT_" + std::to_string( i ) ).c_str(),
                                ( "PredictStateVUT_" + std::to_string( i ) ).c_str() ) );

  for ( int i = 0; i < 8; i++ )
    trees.push_back( new TTree( ( "predictStateUT_" + std::to_string( i ) ).c_str(),
                                ( "PredictStateUT_" + std::to_string( i ) ).c_str() ) );

  for ( int i = 0; i < 1; i++ )
    trees.push_back( new TTree( ( "predictStateUTFUT_" + std::to_string( i ) ).c_str(),
                                ( "PredictStateUTFUT_" + std::to_string( i ) ).c_str() ) );

  for ( int i = 0; i < 2; i++ )
    trees.push_back( new TTree( ( "predictStateUTTF_" + std::to_string( i ) ).c_str(),
                                ( "PredictStateUTTF_" + std::to_string( i ) ).c_str() ) );

  for ( int i = 0; i < 4; i++ )
    trees.push_back( new TTree( ( "predictStateUTT_" + std::to_string( i ) ).c_str(),
                                ( "PredictStateUTT_" + std::to_string( i ) ).c_str() ) );

  for ( int i = 0; i < 2; i++ )
    trees.push_back( new TTree( ( "predictStateTFT_" + std::to_string( i ) ).c_str(),
                                ( "PredictStateTFT_" + std::to_string( i ) ).c_str() ) );

  for ( int i = 0; i < 48; i++ )
    trees.push_back( new TTree( ( "predictStateT_" + std::to_string( i ) ).c_str(),
                                ( "PredictStateT_" + std::to_string( i ) ).c_str() ) );

  for ( int i = 0; i < 2; i++ )
    trees.push_back( new TTree( ( "UpdateStateV_" + std::to_string( i ) ).c_str(),
                                ( "UpdateStateV_" + std::to_string( i ) ).c_str() ) );

  trees.push_back( new TTree( "UpdateStateLV_0", "UpdateStateLV_0" ) );

  trees.push_back( new TTree( "UpdateStateFUT_0", "UpdateStateFUT_0" ) );

  for ( int i = 0; i < 8; i++ )
    trees.push_back( new TTree( ( "UpdateStateUT_" + std::to_string( i ) ).c_str(),
                                ( "UpdateStateUT_" + std::to_string( i ) ).c_str() ) );

  trees.push_back( new TTree( "UpdateStateLUT_0", "UpdateStateLUT_0" ) );

  trees.push_back( new TTree( "UpdateStateFT_0", "UpdateStateFT_0" ) );

  for ( int i = 0; i < 48; i++ )
    trees.push_back( new TTree( ( "UpdateStateT_" + std::to_string( i ) ).c_str(),
                                ( "UpdateStateT_" + std::to_string( i ) ).c_str() ) );

  // Set the branches
  trees[0]->Branch( "sF_sigmaxx_T", &( treeVars->m_sF_P[2][0] ), "sF_sigmaxx_T/D" );
  trees[0]->Branch( "sF_sigmayy_T", &( treeVars->m_sF_P[2][2] ), "sF_sigmayy_T/D" );
  trees[0]->Branch( "sF_sigmatxtx_T", &( treeVars->m_sF_P[2][5] ), "sF_sigmatxtx_T/D" );
  trees[0]->Branch( "sF_sigmatyty_T", &( treeVars->m_sF_P[2][9] ), "sF_sigmatyty_T/D" );
  trees[0]->Branch( "sF_sigmaqopqop_T", &( treeVars->m_sF_P[2][14] ), "sF_sigmaqopqop_T/D" );
  trees[0]->Branch( "sF_sigmaxy_T", &( treeVars->m_sF_P[2][1] ), "sF_sigmaxy_T/D" );
  trees[0]->Branch( "sF_sigmaxtx_T", &( treeVars->m_sF_P[2][3] ), "sF_sigmaxtx_T/D" );
  trees[0]->Branch( "sF_sigmaxty_T", &( treeVars->m_sF_P[2][6] ), "sF_sigmaxty_T/D" );
  trees[0]->Branch( "sF_sigmaxqop_T", &( treeVars->m_sF_P[2][10] ), "sF_sigmaxqop_T/D" );
  trees[0]->Branch( "sF_sigmaytx_T", &( treeVars->m_sF_P[2][4] ), "sF_sigmaytx_T/D" );
  trees[0]->Branch( "sF_sigmayty_T", &( treeVars->m_sF_P[2][7] ), "sF_sigmayty_T/D" );
  trees[0]->Branch( "sF_sigmayqop_T", &( treeVars->m_sF_P[2][11] ), "sF_sigmayqop_T/D" );
  trees[0]->Branch( "sF_sigmatxty_T", &( treeVars->m_sF_P[2][8] ), "sF_sigmatxty_T/D" );
  trees[0]->Branch( "sF_sigmatxqop_T", &( treeVars->m_sF_P[2][12] ), "sF_sigmatxqop_T/D" );
  trees[0]->Branch( "sF_sigmatyqop_T", &( treeVars->m_sF_P[2][13] ), "sF_sigmatyqop_T/D" );

  trees[0]->Branch( "sF_x_T", &( treeVars->m_sF_x[2][0] ), "sF_x_T/D" );
  trees[0]->Branch( "sF_y_T", &( treeVars->m_sF_x[2][1] ), "sF_y_T/D" );
  trees[0]->Branch( "sF_tx_T", &( treeVars->m_sF_x[2][2] ), "sF_tx_T/D" );
  trees[0]->Branch( "sF_ty_T", &( treeVars->m_sF_x[2][3] ), "sF_ty_T/D" );
  trees[0]->Branch( "sF_qop_T", &( treeVars->m_sF_x[2][4] ), "sF_qop_T/D" );

  trees[0]->Branch( "sF_z_T", &( treeVars->m_sF_z[2] ), "sF_z_T/D" );

  trees[0]->Branch( "sF_true_x_T", &( treeVars->m_sF_true_x[2][0] ), "sF_true_x_T/D" );
  trees[0]->Branch( "sF_true_y_T", &( treeVars->m_sF_true_x[2][1] ), "sF_true_y_T/D" );
  trees[0]->Branch( "sF_true_tx_T", &( treeVars->m_sF_true_x[2][2] ), "sF_true_tx_T/D" );
  trees[0]->Branch( "sF_true_ty_T", &( treeVars->m_sF_true_x[2][3] ), "sF_true_ty_T/D" );
  trees[0]->Branch( "sF_true_qop_T", &( treeVars->m_sF_true_x[2][4] ), "sF_true_qop_T/D" );

  trees[0]->Branch( "sF_sigmaxx_V", &( treeVars->m_sF_P[0][0] ), "sF_sigmaxx_V/D" );
  trees[0]->Branch( "sF_sigmayy_V", &( treeVars->m_sF_P[0][2] ), "sF_sigmayy_V/D" );
  trees[0]->Branch( "sF_sigmatxtx_V", &( treeVars->m_sF_P[0][5] ), "sF_sigmatxtx_V/D" );
  trees[0]->Branch( "sF_sigmatyty_V", &( treeVars->m_sF_P[0][9] ), "sF_sigmatyty_V/D" );
  trees[0]->Branch( "sF_sigmaqopqop_V", &( treeVars->m_sF_P[0][14] ), "sF_sigmaqopqop_V/D" );
  trees[0]->Branch( "sF_sigmaxy_V", &( treeVars->m_sF_P[0][1] ), "sF_sigmaxy_V/D" );
  trees[0]->Branch( "sF_sigmaxtx_V", &( treeVars->m_sF_P[0][3] ), "sF_sigmaxtx_V/D" );
  trees[0]->Branch( "sF_sigmaxty_V", &( treeVars->m_sF_P[0][6] ), "sF_sigmaxty_V/D" );
  trees[0]->Branch( "sF_sigmaxqop_V", &( treeVars->m_sF_P[0][10] ), "sF_sigmaxqop_V/D" );
  trees[0]->Branch( "sF_sigmaytx_V", &( treeVars->m_sF_P[0][4] ), "sF_sigmaytx_V/D" );
  trees[0]->Branch( "sF_sigmayty_V", &( treeVars->m_sF_P[0][7] ), "sF_sigmayty_V/D" );
  trees[0]->Branch( "sF_sigmayqop_V", &( treeVars->m_sF_P[0][11] ), "sF_sigmayqop_V/D" );
  trees[0]->Branch( "sF_sigmatxty_V", &( treeVars->m_sF_P[0][8] ), "sF_sigmatxty_V/D" );
  trees[0]->Branch( "sF_sigmatxqop_V", &( treeVars->m_sF_P[0][12] ), "sF_sigmatxqop_V/D" );
  trees[0]->Branch( "sF_sigmatyqop_V", &( treeVars->m_sF_P[0][13] ), "sF_sigmatyqop_V/D" );

  trees[0]->Branch( "sF_x_V", &( treeVars->m_sF_x[0][0] ), "sF_x_V/D" );
  trees[0]->Branch( "sF_y_V", &( treeVars->m_sF_x[0][1] ), "sF_y_V/D" );
  trees[0]->Branch( "sF_tx_V", &( treeVars->m_sF_x[0][2] ), "sF_tx_V/D" );
  trees[0]->Branch( "sF_ty_V", &( treeVars->m_sF_x[0][3] ), "sF_ty_V/D" );
  trees[0]->Branch( "sF_qop_V", &( treeVars->m_sF_x[0][4] ), "sF_qop_V/D" );

  trees[0]->Branch( "sF_z_V", &( treeVars->m_sF_z[0] ), "sF_z_V/D" );

  trees[0]->Branch( "sF_true_x_V", &( treeVars->m_sF_true_x[0][0] ), "sF_true_x_V/D" );
  trees[0]->Branch( "sF_true_y_V", &( treeVars->m_sF_true_x[0][1] ), "sF_true_y_V/D" );
  trees[0]->Branch( "sF_true_tx_V", &( treeVars->m_sF_true_x[0][2] ), "sF_true_tx_V/D" );
  trees[0]->Branch( "sF_true_ty_V", &( treeVars->m_sF_true_x[0][3] ), "sF_true_ty_V/D" );
  trees[0]->Branch( "sF_true_qop_V", &( treeVars->m_sF_true_x[0][4] ), "sF_true_qop_V/D" );

  trees[0]->Branch( "sF_sigmaxx", &( treeVars->m_sF_P[1][0] ), "sF_sigmaxx/D" );
  trees[0]->Branch( "sF_sigmayy", &( treeVars->m_sF_P[1][2] ), "sF_sigmayy/D" );
  trees[0]->Branch( "sF_sigmatxtx", &( treeVars->m_sF_P[1][5] ), "sF_sigmatxtx/D" );
  trees[0]->Branch( "sF_sigmatyty", &( treeVars->m_sF_P[1][9] ), "sF_sigmatyty/D" );
  trees[0]->Branch( "sF_sigmaqopqop", &( treeVars->m_sF_P[1][14] ), "sF_sigmaqopqop/D" );
  trees[0]->Branch( "sF_sigmaxy", &( treeVars->m_sF_P[1][1] ), "sF_sigmaxy/D" );
  trees[0]->Branch( "sF_sigmaxtx", &( treeVars->m_sF_P[1][3] ), "sF_sigmaxtx/D" );
  trees[0]->Branch( "sF_sigmaxty", &( treeVars->m_sF_P[1][6] ), "sF_sigmaxty/D" );
  trees[0]->Branch( "sF_sigmaxqop", &( treeVars->m_sF_P[1][10] ), "sF_sigmaxqop/D" );
  trees[0]->Branch( "sF_sigmaytx", &( treeVars->m_sF_P[1][4] ), "sF_sigmaytx/D" );
  trees[0]->Branch( "sF_sigmayty", &( treeVars->m_sF_P[1][7] ), "sF_sigmayty/D" );
  trees[0]->Branch( "sF_sigmayqop", &( treeVars->m_sF_P[1][11] ), "sF_sigmayqop/D" );
  trees[0]->Branch( "sF_sigmatxty", &( treeVars->m_sF_P[1][8] ), "sF_sigmatxty/D" );
  trees[0]->Branch( "sF_sigmatxqop", &( treeVars->m_sF_P[1][12] ), "sF_sigmatxqop/D" );
  trees[0]->Branch( "sF_sigmatyqop", &( treeVars->m_sF_P[1][13] ), "sF_sigmatyqop/D" );

  trees[0]->Branch( "sF_x", &( treeVars->m_sF_x[1][0] ), "sF_x/D" );
  trees[0]->Branch( "sF_y", &( treeVars->m_sF_x[1][1] ), "sF_y/D" );
  trees[0]->Branch( "sF_tx", &( treeVars->m_sF_x[1][2] ), "sF_tx/D" );
  trees[0]->Branch( "sF_ty", &( treeVars->m_sF_x[1][3] ), "sF_ty/D" );
  trees[0]->Branch( "sF_qop", &( treeVars->m_sF_x[1][4] ), "sF_qop/D" );

  trees[0]->Branch( "sF_z", &( treeVars->m_sF_z[1] ), "sF_z/D" );
  trees[0]->Branch( "sF_chi2", &( treeVars->m_sF_chi2 ), "sF_chi2/D" );
  trees[0]->Branch( "sF_chi2_V", &( treeVars->m_sF_chi2_V ), "sF_chi2_V/D" );
  trees[0]->Branch( "sF_chi2_T", &( treeVars->m_sF_chi2_T ), "sF_chi2_T/D" );
  trees[0]->Branch( "sF_ndof", &( treeVars->m_sF_ndof ), "sF_ndof/D" );

  trees[0]->Branch( "sF_true_x", &( treeVars->m_sF_true_x[1][0] ), "sF_true_x/D" );
  trees[0]->Branch( "sF_true_y", &( treeVars->m_sF_true_x[1][1] ), "sF_true_y/D" );
  trees[0]->Branch( "sF_true_tx", &( treeVars->m_sF_true_x[1][2] ), "sF_true_tx/D" );
  trees[0]->Branch( "sF_true_ty", &( treeVars->m_sF_true_x[1][3] ), "sF_true_ty/D" );
  trees[0]->Branch( "sF_true_qop", &( treeVars->m_sF_true_x[1][4] ), "sF_true_qop/D" );

  trees[0]->Branch( "true_qop_PV", &( treeVars->m_true_qop_PV ), "true_qop_PV/D" );

  trees[0]->Branch( "nHitsV", &( treeVars->m_NHitsV ), "nHitsV/I" );
  trees[0]->Branch( "nHitsUT", &( treeVars->m_NHitsUT ), "nHitsUT/I" );
  trees[0]->Branch( "nHitsT", &( treeVars->m_NHitsT ), "nHitsT/I" );
  trees[0]->Branch( "nHitsTotal", &( treeVars->m_NHitsTotal ), "nHitsTotal/I" );

  trees[0]->Branch( "MCstatus", &( treeVars->m_MC_status ), "MCstatus/I" );

  // set branches for the rest of the trees
  // int k=1;
  for ( auto it = std::next( trees.begin() ); it != trees.end(); ++it ) {
    // std::cout << k << " " << (*it)->GetName() << std::endl;
    // k++;
    ( *it )->Branch( "sigmaxx", &( treeVars->m_P[0] ), "sigmaxx/D" );
    ( *it )->Branch( "sigmayy", &( treeVars->m_P[2] ), "sigmayy/D" );
    ( *it )->Branch( "sigmatxtx", &( treeVars->m_P[5] ), "sigmatxtx/D" );
    ( *it )->Branch( "sigmatyty", &( treeVars->m_P[9] ), "sigmatyty/D" );
    ( *it )->Branch( "sigmaqopqop", &( treeVars->m_P[14] ), "sigmaqopqop/D" );
    ( *it )->Branch( "sigmaxy", &( treeVars->m_P[1] ), "sigmaxy/D" );
    ( *it )->Branch( "sigmaxtx", &( treeVars->m_P[3] ), "sigmaxtx/D" );
    ( *it )->Branch( "sigmaxty", &( treeVars->m_P[6] ), "sigmaxty/D" );
    ( *it )->Branch( "sigmaxqop", &( treeVars->m_P[10] ), "sigmaxqop/D" );
    ( *it )->Branch( "sigmaytx", &( treeVars->m_P[4] ), "sigmaytx/D" );
    ( *it )->Branch( "sigmayty", &( treeVars->m_P[7] ), "sigmayty/D" );
    ( *it )->Branch( "sigmayqop", &( treeVars->m_P[11] ), "sigmayqop/D" );
    ( *it )->Branch( "sigmatxty", &( treeVars->m_P[8] ), "sigmatxty/D" );
    ( *it )->Branch( "sigmatxqop", &( treeVars->m_P[12] ), "sigmatxqop/D" );
    ( *it )->Branch( "sigmatyqop", &( treeVars->m_P[13] ), "sigmatyqop/D" );

    ( *it )->Branch( "x", &( treeVars->m_x[0] ), "x/D" );
    ( *it )->Branch( "y", &( treeVars->m_x[1] ), "y/D" );
    ( *it )->Branch( "tx", &( treeVars->m_x[2] ), "tx/D" );
    ( *it )->Branch( "ty", &( treeVars->m_x[3] ), "ty/D" );
    ( *it )->Branch( "qop", &( treeVars->m_x[4] ), "qop/D" );
    ( *it )->Branch( "sigmaxx_extr", &( treeVars->m_P_extr[0] ), "sigmaxx_extr/D" );
    ( *it )->Branch( "sigmayy_extr", &( treeVars->m_P_extr[2] ), "sigmayy_extr/D" );
    ( *it )->Branch( "sigmatxtx_extr", &( treeVars->m_P_extr[5] ), "sigmatxtx_extr/D" );
    ( *it )->Branch( "sigmatyty_extr", &( treeVars->m_P_extr[9] ), "sigmatyty_extr/D" );
    ( *it )->Branch( "sigmaqopqop_extr", &( treeVars->m_P_extr[14] ), "sigmaqopqop_extr/D" );
    ( *it )->Branch( "sigmaxy_extr", &( treeVars->m_P_extr[1] ), "sigmaxy_extr/D" );
    ( *it )->Branch( "sigmaxtx_extr", &( treeVars->m_P_extr[3] ), "sigmaxtx_extr/D" );
    ( *it )->Branch( "sigmaxty_extr", &( treeVars->m_P_extr[6] ), "sigmaxty_extr/D" );
    ( *it )->Branch( "sigmaxqop_extr", &( treeVars->m_P_extr[10] ), "sigmaxqop_extr/D" );
    ( *it )->Branch( "sigmaytx_extr", &( treeVars->m_P_extr[4] ), "sigmaytx_extr/D" );
    ( *it )->Branch( "sigmayty_extr", &( treeVars->m_P_extr[7] ), "sigmayty_extr/D" );
    ( *it )->Branch( "sigmayqop_extr", &( treeVars->m_P_extr[11] ), "sigmayqop_extr/D" );
    ( *it )->Branch( "sigmatxty_extr", &( treeVars->m_P_extr[8] ), "sigmatxty_extr/D" );
    ( *it )->Branch( "sigmatxqop_extr", &( treeVars->m_P_extr[12] ), "sigmatxqop_extr/D" );
    ( *it )->Branch( "sigmatyqop_extr", &( treeVars->m_P_extr[13] ), "sigmatyqop_extr/D" );

    ( *it )->Branch( "x_extr", &( treeVars->m_x_extr[0] ), "x_extr/D" );
    ( *it )->Branch( "y_extr", &( treeVars->m_x_extr[1] ), "y_extr/D" );
    ( *it )->Branch( "tx_extr", &( treeVars->m_x_extr[2] ), "tx_extr/D" );
    ( *it )->Branch( "ty_extr", &( treeVars->m_x_extr[3] ), "ty_extr/D" );
    ( *it )->Branch( "qop_extr", &( treeVars->m_x_extr[4] ), "qop_extr/D" );
    ( *it )->Branch( "x_prev", &( treeVars->m_x_prev[0] ), "x_prev/D" );
    ( *it )->Branch( "y_prev", &( treeVars->m_x_prev[1] ), "y_prev/D" );
    ( *it )->Branch( "tx_prev", &( treeVars->m_x_prev[2] ), "tx_prev/D" );
    ( *it )->Branch( "ty_prev", &( treeVars->m_x_prev[3] ), "ty_prev/D" );
    ( *it )->Branch( "qop_prev", &( treeVars->m_x_prev[4] ), "qop_prev/D" );

    ( *it )->Branch( "z", &( treeVars->m_z ), "z/D" );
    ( *it )->Branch( "z_prev", &( treeVars->m_z_prev ), "z_prev/D" );

    ( *it )->Branch( "true_z", &( treeVars->m_true_z ), "true_z/D" );
    ( *it )->Branch( "true_x", &( treeVars->m_true_x[0] ), "true_x/D" );
    ( *it )->Branch( "true_y", &( treeVars->m_true_x[1] ), "true_y/D" );
    ( *it )->Branch( "true_tx", &( treeVars->m_true_x[2] ), "true_tx/D" );
    ( *it )->Branch( "true_ty", &( treeVars->m_true_x[3] ), "true_ty/D" );
    ( *it )->Branch( "true_qop", &( treeVars->m_true_x[4] ), "true_qop/D" );
    ( *it )->Branch( "true_x_prev", &( treeVars->m_true_x_prev[0] ), "true_x_prev/D" );
    ( *it )->Branch( "true_y_prev", &( treeVars->m_true_x_prev[1] ), "true_y_prev/D" );
    ( *it )->Branch( "true_tx_prev", &( treeVars->m_true_x_prev[2] ), "true_tx_prev/D" );
    ( *it )->Branch( "true_ty_prev", &( treeVars->m_true_x_prev[3] ), "true_ty_prev/D" );
    ( *it )->Branch( "true_qop_prev", &( treeVars->m_true_x_prev[4] ), "true_qop_prev/D" );
    ( *it )->Branch( "true_qop_here", &( treeVars->m_true_qop_here ), "true_qop_here/D" );

    ( *it )->Branch( "dxdy", &( treeVars->m_hit_dxdy ), "dxdy/D" );
    ( *it )->Branch( "dzdy", &( treeVars->m_hit_dzdy ), "dzdy/D" );
    ( *it )->Branch( "z0", &( treeVars->m_hit_z0 ), "z0/D" );
    ( *it )->Branch( "x0", &( treeVars->m_hit_x0 ), "x0/D" );
    ( *it )->Branch( "y0", &( treeVars->m_hit_y0 ), "y0/D" );
    ( *it )->Branch( "x0_err", &( treeVars->m_hit_x0_err ), "x0_err/D" );
    ( *it )->Branch( "y0_err", &( treeVars->m_hit_y0_err ), "y0_err/D" );

    ( *it )->Branch( "nHit", &( treeVars->m_NHit ), "nHit/I" );
    ( *it )->Branch( "nHitsV", &( treeVars->m_NHitsV ), "nHitsV/I" );
    ( *it )->Branch( "nHitsUT", &( treeVars->m_NHitsUT ), "nHitsUT/I" );
    ( *it )->Branch( "nHitsT", &( treeVars->m_NHitsT ), "nHitsT/I" );
    ( *it )->Branch( "nHitsTotal", &( treeVars->m_NHitsTotal ), "nHitsTotal/I" );

    ( *it )->Branch( "MCstatus", &( treeVars->m_MC_status ), "MCstatus/I" );
  }
}
// 1: Create seed state in VELO
//
// Predicting______________________________________________________
// 2-3:    Predict VELO                      <-> VELO
// 4-5:    Predict last VELO measurement     <-> UT
// 6-13:   Predict UT                        <-> UT
// 14:     Predict UT                        <-> UT (fixed z)
// 15-16:  Predict last UT measurement       <-> T (fixed z)
// 17-20:  Predict last UT measurement       <-> T
// 21-22:  Predict T (fixed z)               <-> T
// 23-70:  Predict T                         <-> T
//
// Updating________________________________________________________
// 71-72:  Update inside     VELO
// 73:     Update last       VELO coming from UT/T
// 74:     Update first      UT   coming from VELO
// 75-82:  Update inside     UT
// 83:     Update last       UT   coming from T
// 84:     Update first      T    coming from VELO/UT
// 85-132: Update inside     T

// Further implementations
#include "ParameterizedKalmanFit_Checker_include.icpp"
