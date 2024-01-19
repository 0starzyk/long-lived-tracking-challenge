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

#include <memory>

#include "KalmanParametrizations.h"
#include "ParameterizedKalmanFit_Methods.h"

#include "Associators/Associators.h"
#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/ODIN.h"
#include "Event/Track.h"
#include "MCInterfaces/IIdealStateCreator.h"
#include "Magnet/DeMagnet.h"
#include "TrackInterfaces/IMeasurementProviderProjector.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include "GaudiKernel/AnyDataHandle.h"
#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "GaudiKernel/System.h"
#include "LHCbAlgs/Transformer.h"

#include <TFile.h>
#include <TTree.h>

using namespace ParKalman;

/** @class ParameterizedKalmanFit_Checker ParameterizedKalmanFit_Checker.h
 *
 *  Algorithm to perfrom a simplified Kalman filter based on parameterized extrapolations
 *  Version that creates tuples to track performance or to do the parameter tuning afterwards
 *
 *  Parameters:
 *  - UseUTHits:             Use hits in the UT in the fit
 *  - UseTHits:              Use hits in the SciFi in the fit
 *  - MaxNumOutlier:         Maximal number of outliers that should be removed
 *  -
 *  - RunForTuning:          Do Kalman steps using truth information to extract information for
 *                           tuning
 *  - OutputTreesFile:       Output location and file name for the tuning tuples
 *  - UseForwMomEstiamte:    Use the momentum estaimte of an forward iteration
 *  - UseForwChi2Estiamte:   Use the chi2 of an forward estimate
 *  - UseOneParameterSet:    Use the parameters for magnet polarity down also for mag up
 *                           (some signs are reversed)
 *  - ParamFileLocation:     Location of the parameter files:
 *                           Has to contain MagUp/ and MagDown/
 *
 *  @author Simon Stemmle
 *  @date   2017-11-02
 */

struct trackTupleInfo;

class ParameterizedKalmanFit_Checker
    : public LHCb::Algorithm::Transformer<LHCb::Tracks( LHCb::Tracks const&, LHCb::ODIN const&, LHCb::LinksByKey const&,
                                                        LHCb::MCProperty const&, DetectorElement const&,
                                                        DeMagnet const& ),
                                          LHCb::DetDesc::usesConditions<DetectorElement, DeMagnet>> {
public:
  /// Standard constructor
  ParameterizedKalmanFit_Checker( const std::string& name, ISvcLocator* pSvcLocator );

  /// Algorithm initialization
  StatusCode initialize() override;

  /// Algorithm execution
  LHCb::Tracks operator()( LHCb::Tracks const&, LHCb::ODIN const& odin, LHCb::LinksByKey const& links,
                           LHCb::MCProperty const&, DetectorElement const& lhcb,
                           DeMagnet const& magnet ) const override;

protected:
private:
  Gaudi::Property<bool> m_UseUT{this, "UseUTHits", true};

  Gaudi::Property<bool> m_UseT{this, "UseTHits", true};

  Gaudi::Property<int> m_MaxNoutlier{this, "MaxNumOutlier", 1};

  Gaudi::Property<bool> m_UseForwardMomEstimate{this, "UseForwMomEstiamte", true};

  Gaudi::Property<bool> m_UseForwardChi2Estimate{this, "UseForwChi2Estiamte", true};

  Gaudi::Property<bool> m_UseOneParameterSet{this, "UseOneParameterSet", false};

  Gaudi::Property<bool> m_RunForTuning{this, "RunForTuning", false};

  Gaudi::Property<std::string> m_TreesFileName{this, "OutputTreesFile", "ParameterizedKalmanTuning"};

  Gaudi::Property<std::string> m_ParamFileLocation{
      this, "ParamFileLocation", System::getEnv( "PARAMFILESROOT" ) + "/data/ParametrizedKalmanFit/FT6x2"};

  // Control options for tuning,testing
  bool m_SetTrueStateAfterUpdate     = false;
  bool m_SetTrueStateAfterPredict    = false;
  bool m_SetTrueStateAfterCreateSeed = false;

  // cache information for the smoother step
  bool m_do_smoother = true;

  // Parametrizations objects for each polarity
  std::unique_ptr<KalmanParametrizations> m_ParExtrUp;
  std::unique_ptr<KalmanParametrizations> m_ParExtrDown;

  //#####
  // Tools
  //#####
  ToolHandle<IMeasurementProviderProjector> m_measProvider = {this, "MeasProvider",
                                                              "MeasurementProvider/MeasurementProvider"};

  // ideal state creator for tuning and performance checks
  ToolHandle<IIdealStateCreator> m_idealStateCreator = {"IdealStateCreator", this};

  // extrapolators
  // 1.For tuning and performance checks
  ToolHandle<ITrackExtrapolator> m_extrapolator = {"TrackMasterExtrapolator/extr1", this};
  // 2.For the extrapolation to the beam pipe
  ToolHandle<ITrackExtrapolator> m_extrapolator_toPV = {"TrackMasterExtrapolator/extr2", this};

  //#################
  // 1. Level methods
  //#################

  /// Method to run the Kalman filter in order to extract tuning information
  StatusCode fitForTuning( trackInfo& tI, const MCTrackInfo& trackInfo, std::vector<TTree*>* trees, trackTupleInfo* tV,
                           IGeometryInfo const& geometry ) const;

  /// Runs the Kalman filter on a track
  StatusCode fit( trackInfo& tI, MCTrackInfo const& trackInfo, std::vector<TTree*>* trees, trackTupleInfo* tV,
                  IGeometryInfo const& geometry ) const;

  /// Create trees that should be filled for tuning and perfomance checks
  void addTrees( std::vector<TTree*>& trees, trackTupleInfo* treeVars ) const;

  //####################################
  // Main methods for the Kalman filter
  //####################################

  /// Load hit information
  void LoadHits_Ch( trackInfo& tI, const ToolHandle<IMeasurementProviderProjector>& measProvider, bool m_UseUT,
                    bool m_UseT, trackTupleInfo* tV ) const;

  /// Method to create a seed state at the first Velo hit
  void CreateVeloSeedState_Ch( int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI,
                               std::vector<TTree*>* trees, trackTupleInfo* tV, IGeometryInfo const& geometry ) const;

  /// General method for updating at a hit
  void UpdateState_Ch( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI,
                       std::vector<TTree*>* trees, trackTupleInfo* tV, IGeometryInfo const& geometry ) const;

  /// General method for predicting to a hit
  bool PredictState_Ch( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI,
                        std::vector<TTree*>* trees, trackTupleInfo* tV, IGeometryInfo const& geometry ) const;

  /// Fill information for comparing default and this kalman filter
  void FillNtuple( Gaudi::Vector5 x, Gaudi::SymMatrix5x5 C, double z, trackInfo& tI, trackTupleInfo* tV,
                   double position, int pos, IGeometryInfo const& geometry ) const;

  //#######################################
  // Further methods for the Kalman filter
  //#######################################

  /// Check if a MC particle is linked to this track
  int MatchesMC( const trackInfo& tI, const MCTrackInfo& trackInfo ) const;

  /// Get true state at a given z position
  bool TrueState( double zpos, double& trueX, double& trueY, double& truetX, double& truetY, double& trueqop,
                  const trackInfo& tI, IGeometryInfo const& geometry, bool initialQop = true ) const;

  /// Method to set the information of the default extrapolator for the tuning
  void fillInfoForExtrapolation( double z_prev, double z, trackInfo& tI, trackTupleInfo* tV,
                                 IGeometryInfo const& geometry ) const;

  /// Predict UT <-> T precise version
  void PredictStateUTT_Ch( Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz, trackInfo& tI,
                           std::vector<TTree*>* trees, trackTupleInfo* tV, IGeometryInfo const& geometry ) const;

  /// Predict UT(fixed z) <-> UT(last layer)
  void PredictStateUTFUT_Ch( int forward, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& P, double& lastz, trackInfo& tI,
                             std::vector<TTree*>* trees, trackTupleInfo* tV, IGeometryInfo const& geometry ) const;

  /// Predict T(fixed z=7783) <-> first T layer
  void PredictStateTFT_Ch( int forward, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& P, double& lastz, trackInfo& tI,
                           std::vector<TTree*>* trees, trackTupleInfo* tV, IGeometryInfo const& geometry ) const;
};

// struct that contains variables for the tupling
struct trackTupleInfo {
  std::array<double, 5>  m_x;
  std::array<double, 5>  m_xTmp;
  std::array<double, 5>  m_x_prev;
  std::array<double, 5>  m_x_extr;
  std::array<double, 15> m_P;
  std::array<double, 15> m_P_extr;
  double                 m_z;
  double                 m_zTmp;
  double                 m_z_prev;

  std::array<double, 5> m_true_x;
  std::array<double, 5> m_true_xTmp;
  std::array<double, 5> m_true_x_prev;
  double                m_true_z;
  double                m_true_qop_here;
  double                m_true_qop_PV;
  int                   m_MC_status;

  std::array<std::array<double, 5>, 3>  m_sF_x;
  std::array<std::array<double, 15>, 3> m_sF_P;
  std::array<std::array<double, 5>, 3>  m_sF_true_x;
  std::array<double, 3>                 m_sF_z;
  double                                m_sF_chi2;
  double                                m_sF_chi2_V;
  double                                m_sF_chi2_T;
  double                                m_sF_chi2_UT;
  double                                m_sF_ndof;

  double m_hit_dxdy;
  double m_hit_dzdy;
  double m_hit_x0;
  double m_hit_y0;
  double m_hit_z0;
  double m_hit_x0_err;
  double m_hit_y0_err;

  // total number of hits (no matter if active or not)
  int m_NHitsTotal;
  // number of VELO hits
  int m_NHitsV;
  // number of UT station hits
  int m_NHitsUT;
  // number of T station hits
  int m_NHitsT;
  int m_NHit;
};
