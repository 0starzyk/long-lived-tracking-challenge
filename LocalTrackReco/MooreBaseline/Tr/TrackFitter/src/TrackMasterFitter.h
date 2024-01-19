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
#ifndef TRACKFITTER_TRACKMASTERFITTER_H
#define TRACKFITTER_TRACKMASTERFITTER_H 1

// Include files
// -------------
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/ToolHandle.h"

// interface base class
#include "Event/Track.h"
#include "TrackInterfaces/IMaterialLocator.h"
#include "TrackInterfaces/IMeasurementProviderProjector.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackInterfaces/ITrackFitter.h"
#include "TrackInterfaces/ITrackKalmanFilter.h"

namespace LHCb {
  class FitNode;
  class State;
  class TrackFitResult;
} // namespace LHCb

/** @class TrackMasterFitter TrackMasterFitter.h
 *
 *
 *  @author Jose Angel Hernando Morata, Eduardo Rodrigues
 *  @date   2005-04-15
 *  reusing the previous code
 *  @author Rutger van der Eijk  07-04-1999
 *  @author Matthew Needham
 */

class TrackMasterFitter : public extends<GaudiTool, ITrackFitter> {
public:
  /// Standard constructor
  using extends::extends;

  StatusCode initialize() override;

private:
  StatusCode operator()( LHCb::Track& track, IGeometryInfo const& geometry, const LHCb::Tr::PID& pid ) const override;

  StatusCode operator()( LHCb::span<LHCb::Track> tracks, IGeometryInfo const& geometry,
                         const LHCb::Tr::PID& pid ) const override;

  void reset() override { m_measProvider->reset(); }

  StatusCode fit_r( LHCb::Track& track, std::any& accelCache, IGeometryInfo const& geometry, LHCb::Tr::PID pid ) const;

  //! initialize reference states for initial trajectory
  StatusCode initializeRefStates( LHCb::Track& track, IGeometryInfo const& geometry, LHCb::Tr::PID pid ) const;

  //! determine track state at various z positions
  StatusCode determineStates( LHCb::Track& track ) const;

  //! remove outliers from the node vector
  LHCb::FitNode* outlierRemoved( LHCb::Track& track ) const;

  //! update the reference vector for each measurement before next iteration
  StatusCode updateRefVectors( LHCb::Track& track, const LHCb::Tr::PID pid, bool doUpdateTransport,
                               std::any& accelCache, IGeometryInfo const& geometry ) const;

  //! projectReference state
  StatusCode projectReference( LHCb::Track& track ) const;

  //! Retrieve the number of nodes with a measurement
  unsigned int nNodesWithMeasurement( const LHCb::Track& track ) const;

  //! Create the nodes from the measurements
  StatusCode makeNodes( LHCb::Track& track, const LHCb::Tr::PID pid, std::any& accelCache,
                        IGeometryInfo const& geometry ) const;

  //! Update material corrections stored in nodes
  StatusCode updateMaterialCorrections( LHCb::Track& track, const LHCb::Tr::PID pid, std::any& accelCache,
                                        IGeometryInfo const& geometry ) const;

  //! Update transport matrices stored in nodes
  StatusCode updateTransport( LHCb::Track& track, IGeometryInfo const& geometry ) const;

  const ITrackExtrapolator* extrapolator( LHCb::Track::Types tracktype ) const {
    if ( ( tracktype == LHCb::Track::Types::Velo ) || ( tracktype == LHCb::Track::Types::VeloBackward ) )
      return &( *m_veloExtrapolator );
    return &( *m_extrapolator );
  }

  /// allocate a cache to be used with fit_r
  std::any createCache() const { return m_materialLocator.get()->createCache(); }

private:
  ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackMasterExtrapolator"}; ///< extrapolator
  ToolHandle<ITrackExtrapolator> m_veloExtrapolator{this, "VeloExtrapolator",
                                                    "TrackLinearExtrapolator"}; ///< extrapolator for Velo-only tracks
  ToolHandle<ITrackKalmanFilter> m_trackNodeFitter{
      this, "NodeFitter", "TrackKalmanFilter"}; ///< delegate to actual track fitter (which fits from nodes)
  ToolHandle<IMeasurementProviderProjector> m_measProvider{this, "MeasProvider", "MeasurementProvider"};
  ToolHandle<IMaterialLocator>              m_materialLocator{this, "MaterialLocator", "DetailedMaterialLocator"};

private:
  Gaudi::Property<bool>   m_upstream{this, "FitUpstream", true, "switch between upstream/downstream fit"};
  Gaudi::Property<bool>   m_addDefaultRefNodes{this, "AddDefaultReferenceNodes", true, "add default reference nodes"};
  Gaudi::Property<bool>   m_stateAtBeamLine{this, "StateAtBeamLine", true, "add state closest to the beam-line"};
  Gaudi::Property<int>    m_numFitIter{this, "NumberFitIterations", 10, "number of fit iterations to perform"};
  Gaudi::Property<double> m_chi2Outliers{this, "Chi2Outliers", 9.0, "chi2 of outliers to be removed"};
  Gaudi::Property<int>    m_numOutlierIter{this, "MaxNumberOutliers", 2, "max number of outliers to be removed"};
  Gaudi::Property<bool>   m_useSeedStateErrors{this, "UseSeedStateErrors", false, "use errors of the seed state"};
  Gaudi::Property<bool>   m_useClassicalSmoother{this, "UseClassicalSmoother", false, "Use classical smoother"};
  Gaudi::Property<bool>   m_fillExtraInfo{this, "FillExtraInfo", true, "Fill the extra info"};

  Gaudi::Property<double>              m_errorX{this, "ErrorX", 20.0 * Gaudi::Units::mm, "Seed error on x"};
  Gaudi::Property<double>              m_errorY{this, "ErrorY", 20.0 * Gaudi::Units::mm, "Seed error on y"};
  Gaudi::Property<double>              m_errorTx{this, "ErrorTx", 0.1, "Seed error on slope x"};
  Gaudi::Property<double>              m_errorTy{this, "ErrorTy", 0.1, "Seed error on slope y"};
  Gaudi::Property<std::vector<double>> m_errorQoP{this, "ErrorQoP", {0.0, 0.01}, "Seed error on QoP"};

  Gaudi::Property<bool> m_makeNodes{this, "MakeNodes", false};
  Gaudi::Property<bool> m_makeMeasurements{this, "MakeMeasurements", false};
  Gaudi::Property<bool> m_updateTransport{this, "UpdateTransport", true,
                                          "Update the transport matrices between iterations"};
  Gaudi::Property<int>  m_maxUpdateTransports{this, "MaxUpdateTransports", 10,
                                             "Update transport only n-times during iterations"};
  Gaudi::Property<bool> m_updateMaterial{this, "UpdateMaterial", false,
                                         "Update material corrections between iterations"};
  Gaudi::Property<bool> m_updateReferenceInOutlierIters{
      this, "UpdateReferenceInOutlierIterations", true,
      "Update projection in iterations in which outliers are removed"};
  Gaudi::Property<double> m_minMomentumForELossCorr{this, "MinMomentumELossCorr", 10. * Gaudi::Units::MeV,
                                                    "Minimum momentum used in correction for energy loss"};
  Gaudi::Property<bool>   m_applyMaterialCorrections{this, "ApplyMaterialCorrections", true,
                                                   "Apply material corrections"};
  Gaudi::Property<bool>   m_applyEnergyLossCorrections{this, "ApplyEnergyLossCorr", true,
                                                     "Apply energy loss corrections"};
  Gaudi::Property<double> m_maxDeltaChi2Converged{this, "MaxDeltaChiSqConverged", 0.01,
                                                  "Maximum change in chisquare for converged fit"};

  Gaudi::Property<double> m_scatteringPt{
      this, "TransverseMomentumForScattering", 400. * Gaudi::Units::MeV,
      "transverse momentum used for scattering if track has no good momentum estimate"};
  Gaudi::Property<double> m_scatteringP{this, "MomentumForScattering", -1,
                                        "momentum used for scattering in e.g. magnet off data"};
  Gaudi::Property<double> m_minMomentumForScattering{this, "MinMomentumForScattering", 100. * Gaudi::Units::MeV,
                                                     "Minimum momentum used for scattering"};
  Gaudi::Property<double> m_maxMomentumForScattering{this, "MaxMomentumForScattering", 500. * Gaudi::Units::GeV,
                                                     "Maximum momentum used for scattering"};
  Gaudi::Property<size_t> m_minNumVPHits{this, "MinNumVPHitsForOutlierRemoval", 3, "Minimum number of VP layers"};
  Gaudi::Property<size_t> m_minNumUTHits{this, "MinNumUTHitsForOutlierRemoval", 3, "Minimum number of UT layers"};
  Gaudi::Property<size_t> m_minNumTHits{this, "MinNumTHitsForOutlierRemoval", 6, "Minimum number of T layers"};
  Gaudi::Property<size_t> m_minNumMuonHits{this, "MinNumMuonHitsForOutlierRemoval", 4, "Minimum number of Muon layers"};

  // job options
  std::string m_extrapolatorName;     ///< name of the extrapolator in Gaudi
  std::string m_veloExtrapolatorName; ///< name of the velo-only extrapolator

  //! helper to print a failure comment
  StatusCode failure( const std::string& comment ) const;
  StatusCode failureInfo( const std::string& comment ) const;

  bool                  m_debugLevel;
  Gaudi::Property<bool> m_useFastMaterialApproximation{this, "FastMaterialApproximation", false,
                                                       "Use fast approximation of scattering corrections"};
};

#endif // TRACKFITTER_TRACKKALMANFILTER_H
