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
#pragma once

#include <optional>

#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "GaudiKernel/ToolHandle.h"
#include "Kernel/Trajectory.h"
#include "TrackInterfaces/IMeasurementProviderProjector.h"

#include "Event/MCParticle.h"
#include "Event/Track.h"

#include "TrackInterfaces/IMeasurementProviderProjector.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include "KalmanParametrizations.h"

#include "Associators/Associators.h"

#include <TFile.h>
//########################################################################
//
// Some common content for the different Parameterized Kalman algorithms
//
// 2017-10-30: Simon Stemmle
//
//########################################################################

namespace ParKalman {

  struct FitterHit {
    ROOT::Math::XYZPoint  point;
    ROOT::Math::XYZVector direction;
    double                error;
  };

  ////////////////////////////////////////////////////////////////////////////
  // struct that contains all temporary information needed during a track fit
  ////////////////////////////////////////////////////////////////////////////
  using OutputTrack = LHCb::Event::v1::Track;
  using InputTrack  = LHCb::Event::v2::Track;
  using FTYPE       = double;
  using Point =
      typename ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>;
  struct trackInfo {

    trackInfo() = default;
    trackInfo( const LHCb::LinksByKey& _links ) : links{_links} {}

    // Should additional information be cached to run the smoother afterwards
    bool m_do_smoother = false;
    // file to store output trees (only for this event in order to be threadsafe?)
    TFile* m_TreesFile;
    // pointer to the current track
    const InputTrack* m_inputTrack;
    OutputTrack*      m_track;
    // best momentum estimate of the track
    double m_BestMomEst;
    // total number of hits (no matter if active or not)
    int m_NHitsTotal;
    // number of VELO hits
    int m_NHitsV;
    // number of UT station hits
    int m_NHitsUT;
    // number of T station hits
    int m_NHitsT;
    // x,y position and uncertainty Velo
    std::array<double, 25> m_XMeasV;
    std::array<double, 25> m_YMeasV;
    std::array<double, 25> m_ZMeasV;
    std::array<double, 25> m_XErrV;
    std::array<double, 25> m_YErrV;
    // strip trajectories
    std::array<FitterHit, 4>  m_measUtLite;
    std::array<FitterHit, 12> m_measFtLite;
    std::array<double, 4>     m_MeasErrUT;
    std::array<double, 12>    m_MeasErrT;
    // position of T hits
    std::array<double, 12> m_PosHitT;
    // mask for status of hit: active (use it,1) or inactive (don't use it,0)
    std::array<int, 41> m_HitStatus;
    // array for the chi2/ndf for each hit
    std::array<double, 41> m_HitChi2;
    // variables that track the chi2s
    double m_chi2;
    double m_chi2T;
    double m_chi2V;
    // variables that track the ndofs
    int m_Ndof;
    int m_NdofT;
    int m_NdofUT;
    int m_NdofV;
    // LHCb IDs of initial hits
    std::array<LHCb::LHCbID, 41> m_lhcbIDs;
    // mask for layers in UT stations: is there a hit or not
    std::array<int, 4> m_HasHitUT;
    // matching the hit to the UT layer
    std::array<int, 4> m_UTHitToUTLayer;
    // mask for layers in T stations: is there a hit or not
    std::array<int, 12> m_HasHitT;
    // matching the hit to the TT layer
    std::array<int, 12> m_THitToTLayer;
    // previous hit in UT
    int m_PrevNUT;
    // previous hit in T
    int m_PrevNT;
    // final z positions
    std::array<double, 41> m_StateZPos;
    // forward predicted states
    std::array<Gaudi::Vector5, 41> m_StateForwardPredicted;
    // backward predicted states
    std::array<Gaudi::Vector5, 41> m_StateBackwardPredicted;
    // forward updated states
    std::array<Gaudi::Vector5, 41> m_StateForwardUpdated;
    // backward updated states
    std::array<Gaudi::Vector5, 41> m_StateBackwardUpdated;
    // smoothed/averaged states
    std::array<Gaudi::Vector5, 41> m_StateSmoothed;
    // forward predicted cov
    std::array<Gaudi::SymMatrix5x5, 41> m_CovForwardPredicted;
    // backward predicted cov
    std::array<Gaudi::SymMatrix5x5, 41> m_CovBackwardPredicted;
    // forward updated cov
    std::array<Gaudi::SymMatrix5x5, 41> m_CovForwardUpdated;
    // backward updated cov
    std::array<Gaudi::SymMatrix5x5, 41> m_CovBackwardUpdated;
    // smoothed cov
    std::array<Gaudi::SymMatrix5x5, 41> m_CovSmoothed;
    // propagation matrices
    std::array<Gaudi::Matrix5x5, 41> m_PropForward;
    // reference states/jacobian for the intermediate extrapolations
    Gaudi::Vector5 m_RefStateForwardV;
    Gaudi::Vector5 m_RefStateForwardFUT;
    Gaudi::Vector5 m_RefStateForwardUT;
    Gaudi::Vector5 m_RefStateForwardT;

    Gaudi::Matrix5x5 m_RefPropForwardVUT;
    Gaudi::Matrix5x5 m_RefPropForwardUTT;

    // pointer to the extrapolator that should be used
    const KalmanParametrizations* m_extr;

    // MC linker: only needed for ParameterizedKalmanFit_Checker
    std::optional<InputLinks<ContainedObject, LHCb::MCParticle>> links;
  };

  ////////////////////////////////////////
  // Load hit information
  ////////////////////////////////////////
  void LoadHits( trackInfo& tI, bool m_UseUT, bool m_UseT );

  void LoadHitsFromTrackV1( trackInfo& tI, const ToolHandle<IMeasurementProviderProjector>& measProvider, bool m_UseUT,
                            bool m_UseT );

  ///////////////////////////////////////////
  // Method to create a seed state at the first Velo hit
  ///////////////////////////////////////////
  void CreateVeloSeedState( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  //////////////////////////////////////////
  // General method for updating at a hit
  //////////////////////////////////////////
  void UpdateState( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  //////////////////////////////////////////
  // Predict inside the VELO
  //////////////////////////////////////////
  void PredictStateV( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  //////////////////////////////////////////
  // Predict VELO <-> UT
  //////////////////////////////////////////
  bool PredictStateVUT( Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  //////////////////////////////////////////
  // Predict UT <-> UT
  //////////////////////////////////////////
  void PredictStateUT( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  ///////////////////////////////////////////////////////////////
  // Predict from last UT layer to start of UTTF (or vice versa)
  ///////////////////////////////////////////////////////////////
  void PredictStateUTFUT( int forward, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  //////////////////////////////////////////
  // Predict UT <-> T precise version
  //////////////////////////////////////////
  void PredictStateUTT( Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  //////////////////////////////////////////
  // Predict T <-> T
  //////////////////////////////////////////
  void PredictStateT( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  //////////////////////////////////////////////
  // Predict T(fixed z=7783) <-> first T layer
  //////////////////////////////////////////////
  void PredictStateTFT( int forward, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  /////////////////////////////////////////
  // Update state with velo measurement
  /////////////////////////////////////////
  void UpdateStateV( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, trackInfo& tI );

  /////////////////////////////////////////
  // Update state with UT measurement
  /////////////////////////////////////////
  void UpdateStateUT( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  /////////////////////////////////////////
  // Update state with T measurement
  /////////////////////////////////////////
  void UpdateStateT( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI );

  ///////////////////////////////////////////////
  //  Smothe/average method
  ///////////////////////////////////////////////
  bool AverageState( int nHit, trackInfo& tI );

  ////////////////////////////////////////////////////////////////
  // extrapolate to the vertex using the default extrpolator
  ////////////////////////////////////////////////////////////////
  void ExtrapolateToVertex( Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, IGeometryInfo const& geometry,
                            ToolHandle<ITrackExtrapolator> const& m_extrapolator_toPV );

  ///////////////////////////////////////////////////////////////////////////
  // Create Track
  ///////////////////////////////////////////////////////////////////////////
  void addInfoToTrack( const Gaudi::Vector5& x, const Gaudi::SymMatrix5x5& C, double z, trackInfo& tI );

  ///////////////////////////////////////////////////////////////////////////
  // Add state to track
  ///////////////////////////////////////////////////////////////////////////
  void addStateToTrack( const Gaudi::Vector5& x, const Gaudi::SymMatrix5x5& C, double z, trackInfo& tI,
                        LHCb::State::Location loc );

  ///////////////////////////////////////////////////////////////////////////
  // Check if outliers should be removed and remove one of them
  ///////////////////////////////////////////////////////////////////////////
  bool DoOutlierRemoval( trackInfo& tI );

  // Matrix operations
  void Similarity_1x2_S5x5_2x1( const Gaudi::Vector2& A, const Gaudi::SymMatrix5x5& B, double& R );

  void SymmetricTensorProduct5( const Gaudi::Vector5 A, Gaudi::SymMatrix5x5& RM );

  void Multiply_S5x5_2x1( const Gaudi::SymMatrix5x5& AM, const Gaudi::Vector2& B, Gaudi::Vector5& R );

  void Multiply_S5x5_S2x2( const Gaudi::SymMatrix5x5& AM, const Gaudi::SymMatrix2x2& BM,
                           ROOT::Math::SMatrix<double, 5, 2>& RM );
} // namespace ParKalman
