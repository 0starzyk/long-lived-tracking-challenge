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

#include "KalmanParametrizations.h"
#include "ParameterizedKalmanFit_Methods.h"

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/Track.h"
#include "Event/Track_v2.h"
#include "FTDet/DeFTDetector.h"
#include "GaudiKernel/AnyDataHandle.h"
#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "GaudiKernel/System.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Transformer.h"
#include "Magnet/DeMagnet.h"
#include "TrackInterfaces/IMeasurementProviderProjector.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include <memory>

using namespace ParKalman;

/** @class ParameterizedKalmanFit ParameterizedKalmanFit.h
 *
 *  Algorithm to perfrom a simplified Kalman filter based on parameterized extrapolations
 *
 *  Parameters:
 *  - UseUTHits:             Use hits in the UT in the fit
 *  - UseTHits:              Use hits in the SciFi in the fit
 *  - MaxNumOutlier:         Maximal number of outliers that should be removed
 *  - UseForwMomEstiamte:    Use the momentum estaimte of an forward iteration
 *  - UseForwChi2Estiamte:   Use the chi2 of an forward estimate
 *  - UseOneParameterSet:    Use the parameters for magnet polarity down also for mag up
 *                           (some signs are reversed)
 *  - ParamFileLocation:     Location of the parameter files:
 *                           Has to contain MagUp/ and MagDown/
 *
 *  @author Simon Stemmle
 *  @date   2017-10-26
 */

class ParameterizedKalmanFit
    : public LHCb::Algorithm::Transformer<std::vector<LHCb::Event::v1::Track>(
                                              std::vector<LHCb::Event::v2::Track> const&, DetectorElement const&,
                                              DeMagnet const&, DeFT const& ),
                                          LHCb::DetDesc::usesConditions<DetectorElement, DeMagnet, DeFT>> {
  using OutputTrack = LHCb::Event::v1::Track;
  using InputTrack  = LHCb::Event::v2::Track;

public:
  /// Standard constructor
  ParameterizedKalmanFit( const std::string& name, ISvcLocator* pSvcLocator );

  /// Algorithm initialization
  StatusCode initialize() override;

  /// Algorithm execution
  std::vector<OutputTrack> operator()( std::vector<InputTrack> const&, DetectorElement const&, DeMagnet const&,
                                       DeFT const& ) const override;

private:
  Gaudi::Property<bool> m_UseUT{this, "UseUTHits", true};

  Gaudi::Property<bool> m_UseT{this, "UseTHits", true};

  Gaudi::Property<int> m_MaxNoutlier{this, "MaxNumOutlier", 1};

  Gaudi::Property<bool> m_UseForwardMomEstimate{this, "UseForwMomEstimate", true};

  Gaudi::Property<bool> m_UseForwardChi2Estimate{this, "UseForwChi2Estimate", true};

  Gaudi::Property<bool> m_UseOneParameterSet{this, "UseOneParameterSet", false};

  Gaudi::Property<std::string> m_ParamFileLocation{
      this, "ParamFileLocation", System::getEnv( "PARAMFILESROOT" ) + "/data/ParametrizedKalmanFit/FT6x2"};

  // cache information for the smoother step
  bool m_do_smoother = true;

  // Parametrizations objects for each polarity
  std::unique_ptr<KalmanParametrizations> m_ParExtrUp;
  std::unique_ptr<KalmanParametrizations> m_ParExtrDown;

  //#####
  // Tools
  //#####

  // extrapolators for the extrapolation to the beam pipe
  ToolHandle<ITrackExtrapolator> m_extrapolator_toPV{this, "Extrapolator", "TrackMasterExtrapolator"};

  ToolHandle<IMeasurementProviderProjector> m_measProvider = {this, "MeasProvider",
                                                              "MeasurementProvider/MeasurementProvider"};

  //#################
  // 1. Main method
  //#################

  /// Runs the Kalman filter on a track
  StatusCode fit( trackInfo& tI, IGeometryInfo const& geometry ) const;

  //##########################################################################################
  // Methods for the Kalman filter that are not implemented in ParameterizedKalmanFit_Methods
  //##########################################################################################

  /// General method for predicting to a hit
  bool PredictState( int forward, int nHit, Gaudi::Vector5& x, Gaudi::SymMatrix5x5& C, double& lastz,
                     trackInfo& tI ) const;
};
