/*****************************************************************************\
* (c) Copyright 2000-2021 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Detector/FT/FTChannelID.h"
#include "Detector/FT/FTConstants.h"
#include "Detector/FT/FTUtils.h"
#include "Event/FitNode.h"
#include "Event/Measurement.h"
#include "Event/PrFitNode.h"
#include "Event/PrKalmanFitResult.h"
#include "Event/Track.h"
#include "Event/TrackFitResult.h"
#include "Event/TrackParameters.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Consumer.h"
#include "Map.h"
#include "TrackKernel/TrackFunctors.h"
#include "TrackMonitorTupleBase.h"
#include <Gaudi/Accumulators/Histogram.h>

template <typename TFitResult, typename TNode>
class FTTrackMonitor
    : public LHCb::Algorithm::Consumer<void( const LHCb::Track::Range& tracks, const DetectorElement& ),
                                       LHCb::DetDesc::usesBaseAndConditions<TrackMonitorTupleBase, DetectorElement>> {

public:
  FTTrackMonitor( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"TracksInContainer", LHCb::TrackLocation::Default},
                   KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}} ){};

  StatusCode initialize() override {
    return Consumer::initialize().andThen( [&] {
      if ( m_verboseMode.value() ) {
        for ( int i = 0; i < 256; i++ ) {
          m_histograms_modules.emplace( std::piecewise_construct, std::forward_as_tuple( i ),
                                        std::forward_as_tuple( this, i ) );
        }
      }
    } );
  }

  void operator()( const LHCb::Track::Range& tracks, const DetectorElement& lhcb ) const override;

private:
  void fillContainers( const LHCb::Track& track, IGeometryInfo const& geometry ) const;

  std::array<Gaudi::Property<double>, 3> m_refVec{
      {{this, "ReferenceZT1", 7931.0 * Gaudi::Units::mm, "midpoint of FTStation 1"},
       {this, "ReferenceZT2", 8613.0 * Gaudi::Units::mm, "midpoint of FTStation 2"},
       {this, "ReferenceZT3", 9298.0 * Gaudi::Units::mm, "midpoint of FTStation 3"}}};

  Gaudi::Property<bool> m_expertMode{this, "ExpertMode", false};
  Gaudi::Property<bool> m_verboseMode{this, "VerboseMode", false};

private:
  // Plots for the whole SciFi
  mutable Gaudi::Accumulators::Histogram<1> m_unbiasedFTResidual{
      this, "unbiasedFTResidual", "unbiasedFTResidual", {200, -2. * Gaudi::Units::mm, 2. * Gaudi::Units::mm, "X"}};
  mutable Gaudi::Accumulators::Histogram<1> m_biasedFTResidual{
      this, "biasedFTResidual", "biasedFTResidual", {200, -2. * Gaudi::Units::mm, 2. * Gaudi::Units::mm, "X"}};
  mutable Gaudi::Accumulators::Histogram<2> m_unbiasedFTResidualLayer{
      this,
      "unbiasedResidualLayer",
      "unbiasedResidual per FTLayer",
      {{FTConstants::nLayersTotal, -0.5, FTConstants::nLayersTotal - 0.5, "X"},
       {200, -2. * Gaudi::Units::mm, 2. * Gaudi::Units::mm, "Y"}}};
  mutable Gaudi::Accumulators::Histogram<2> m_biasedFTResidualLayer{
      this,
      "biasedResidualLayer",
      "biasedResidual per FTLayer",
      {{FTConstants::nLayersTotal, -0.5, FTConstants::nLayersTotal - 0.5, "X"},
       {200, -2. * Gaudi::Units::mm, 2. * Gaudi::Units::mm, "Y"}}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_RMSResidualModules{
      this,
      "RMSResidualModules",
      "Mean Residual (rms-unbiased) in each module",
      {FTConstants::nModulesTotal, -0.5, FTConstants::nModulesTotal - 0.5}};
  mutable Gaudi::Accumulators::Histogram<2> m_residualPerModule{
      this,
      "residualPerModule",
      "residual per module",
      {{FTConstants::nModulesTotal, -0.5, FTConstants::nModulesTotal - 0.5, "X"}, {50, -1., 1., "Y"}}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_RMSResidualQuarters{
      this,
      "RMSResidualQuarters",
      "Mean Residual (rms-unbiased) in each quarter",
      {FTConstants::nQuartersTotal, -0.5, FTConstants::nQuartersTotal - 0.5}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_UnbiasedResidualModules{
      this,
      "UnbiasedResidualModules",
      "Unbiased residual in each module",
      {FTConstants::nModulesTotal, -0.5, FTConstants::nModulesTotal - 0.5}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_UnbiasedResidualQuarters{
      this,
      "UnbiasedResidualQuarters",
      "Unbiased residual in each quarter",
      {FTConstants::nQuartersTotal, -0.5, FTConstants::nQuartersTotal - 0.5}};
  mutable Gaudi::Accumulators::Histogram<1> m_globalSiPMIdx{
      this,
      "globalSiPMIdx",
      "globalSiPMIdx of nodes",
      {FTConstants::nSiPMsTotal, 0, FTConstants::nSiPMsTotal, "globalSiPMIdx"}};

  // Plots split by SciFi module
  struct HistoModules {
    mutable Gaudi::Accumulators::Histogram<1> UnbiasedResidualModules_permodule;

    HistoModules( const FTTrackMonitor* owner, unsigned int module_id )
        : UnbiasedResidualModules_permodule{owner,
                                            LHCb::Detector::FTUtils::globalModuleToLocalString( module_id ),
                                            LHCb::Detector::FTUtils::globalModuleToLocalString( module_id ),
                                            {200, -2. * Gaudi::Units::mm, 2. * Gaudi::Units::mm, "X"}} {}
  };

  // Plots split by SciFi station
  struct StationTrackHistogram {
    mutable Gaudi::Accumulators::ProfileHistogram<1> RMSResidualModules;
    mutable Gaudi::Accumulators::ProfileHistogram<1> UnbiasedResidualModules;
    mutable Gaudi::Accumulators::Histogram<1>        x;
    mutable Gaudi::Accumulators::Histogram<1>        y;
    mutable Gaudi::Accumulators::Histogram<1>        Tx;
    mutable Gaudi::Accumulators::Histogram<1>        Ty;
    mutable Gaudi::Accumulators::Histogram<1>        Xdist;
    mutable Gaudi::Accumulators::Histogram<1>        XdistWide;
    mutable Gaudi::Accumulators::Histogram<1>        Ydist;
    mutable Gaudi::Accumulators::Histogram<2>        pos;
    mutable Gaudi::Accumulators::Histogram<2>        slopes;

    StationTrackHistogram( const FTTrackMonitor* owner, std::string const& station )
        : RMSResidualModules{owner,
                             "RMSResidualModules" + station,
                             "Residual (rms-unbiased) in FTStation " + station,
                             {FTConstants::nModulesMax * FTConstants::nQuarters * FTConstants::nLayers, -0.5,
                              ( FTConstants::nModulesMax * FTConstants::nQuarters * FTConstants::nLayers ) - 0.5}}
        , UnbiasedResidualModules{owner,
                                  "UnbiasedResidualModules" + station,
                                  "Unbiased Residual in FTStation " + station,
                                  {FTConstants::nModulesMax * FTConstants::nQuarters * FTConstants::nLayers, -0.5,
                                   ( FTConstants::nModulesMax * FTConstants::nQuarters * FTConstants::nLayers ) - 0.5}}
        , x{owner, "x" + station, "x in FTStation " + station, {200, -3000, 3000}}
        , y{owner, "y" + station, "y in FTStation " + station, {200, -3000, 3000}}
        , Tx{owner, "tx" + station, "tx in FTStation " + station, {200, -0.2, 0.2}}
        , Ty{owner, "ty" + station, "ty in FTStation " + station, {200, -0.2, 0.2}}
        , Xdist{owner,
                "xdist" + station,
                "x dist in FTStation " + station,
                {200, -100 * Gaudi::Units::cm, 100 * Gaudi::Units::cm}}
        , XdistWide{owner,
                    "xdistwide" + station,
                    "x dist in FTStation " + station,
                    {200, -350 * Gaudi::Units::cm, 350 * Gaudi::Units::cm}}
        , Ydist{owner,
                "ydist" + station,
                "y dist in FTStation " + station,
                {200, -100 * Gaudi::Units::cm, 100 * Gaudi::Units::cm}}
        , pos{owner,
              "pos" + station,
              "position in FTStation " + station + ";x [mm]; y [mm]",
              {200, -3000., 3000.},
              {200, -3000., 3000}}
        , slopes{owner,
                 "slopes" + station,
                 "slopes in FTStation " + station + ";t_{x}; t_{y}",
                 {200, -0.2, 0.2},
                 {200, -0.2, 0.2}} {}
  };

  std::map<LHCb::Detector::FTChannelID::StationID, StationTrackHistogram> m_histograms =
      []( const FTTrackMonitor* parent ) {
        std::map<LHCb::Detector::FTChannelID::StationID, StationTrackHistogram> map;
        for ( int i : {1, 2, 3} ) {
          map.emplace( std::piecewise_construct, std::forward_as_tuple( LHCb::Detector::FTChannelID::StationID( i ) ),
                       std::forward_as_tuple( parent, "T" + std::to_string( i ) ) );
        }
        return map;
      }( this );

  std::map<unsigned int, HistoModules> m_histograms_modules;
};

using FTTrMonitor          = FTTrackMonitor<LHCb::TrackFitResult, LHCb::FitNode>;
using FTTrMonitor_PrKalman = FTTrackMonitor<LHCb::PrKalmanFitResult, LHCb::Pr::Tracks::Fit::Node>;
DECLARE_COMPONENT_WITH_ID( FTTrMonitor, "FTTrackMonitor" )
DECLARE_COMPONENT_WITH_ID( FTTrMonitor_PrKalman, "FTTrackMonitor_PrKalman" )

template <typename TFitResult, typename TNode>
void FTTrackMonitor<TFitResult, TNode>::operator()( const LHCb::Track::Range& tracks,
                                                    const DetectorElement&    lhcb ) const {
  auto& geometry = *lhcb.geometry();

  for ( const LHCb::Track* track : tracks ) {
    const auto type = track->type();
    if ( type != LHCb::Track::Types::Downstream && type != LHCb::Track::Types::Long &&
         type != LHCb::Track::Types::Ttrack ) {
      continue;
    }
    if ( track->checkFitStatus( LHCb::Track::FitStatus::Fitted ) ) { fillContainers( *track, geometry ); }
  }
}

template <typename TFitResult, typename TNode>
void FTTrackMonitor<TFitResult, TNode>::fillContainers( const LHCb::Track&   track,
                                                        IGeometryInfo const& geometry ) const {

  const auto fitResult = dynamic_cast<const TFitResult*>( track.fitResult() );
  if ( !fitResult ) {
    throw GaudiException( "Fit result is NULL - track was not fitted or uses wrong fit result type", name(),
                          StatusCode::FAILURE );
  }

  // space to put info about whole tracks
  const double     chi2NDOF    = track.chi2PerDoF();
  const double     probChi2    = track.probChi2();
  const double     ghostProb   = track.ghostProbability();
  const double     phi         = track.phi();
  const double     eta         = track.pseudoRapidity();
  const double     p           = ( track.p() ) / Gaudi::Units::GeV;  // p in GeV
  const double     pt          = ( track.pt() ) / Gaudi::Units::GeV; // p in GeV
  const auto       momentumVec = track.momentum();
  const double     px          = momentumVec.X();
  const double     py          = momentumVec.Y();
  const double     pz          = momentumVec.Z();
  Gaudi::XYZVector slopes      = track.slopes();
  const double     tx          = slopes.X();
  const double     ty          = slopes.Y();

  if ( m_expertMode.value() ) {
    Tuple trackTuple = nTuple( "FTTrackTuple_tracks", "" );

    trackTuple->column( "phi", phi ).ignore();
    trackTuple->column( "eta", eta ).ignore();

    trackTuple->column( "chi2PerDoF", chi2NDOF ).ignore();
    trackTuple->column( "probChi2", probChi2 ).ignore();
    trackTuple->column( "ghostProb", ghostProb ).ignore();

    // TODO: consider including track type

    trackTuple->column( "p", p ).ignore();
    trackTuple->column( "pt", pt ).ignore();
    trackTuple->column( "px", px ).ignore();
    trackTuple->column( "py", py ).ignore();
    trackTuple->column( "pz", pz ).ignore();
    trackTuple->column( "tx", tx ).ignore();
    trackTuple->column( "ty", ty ).ignore();
    trackTuple->write().ignore();
  }

  for ( const auto& node : nodes( *fitResult ) ) {
    if ( !node.hasMeasurement() ) continue;
    if ( node.isHitOnTrack() && node.isFT() ) {
      LHCb::LHCbID lhcbID = id( node );
      assert( lhcbID.isFT() );

      LHCb::Detector::FTChannelID chan                  = lhcbID.ftID();
      unsigned int                station               = to_unsigned( chan.station() );
      unsigned int                uniqueLayer           = chan.globalLayerIdx();
      unsigned int                uniqueQuarter         = chan.globalQuarterIdx();
      unsigned int                uniqueModule          = chan.globalModuleIdx();
      unsigned int                uniqueModuleInStation = chan.localModuleIdx_station();
      const auto                  globalMatIdx  = chan.localMatIdx() + FTConstants::nMats * chan.globalModuleIdx();
      const auto                  globalSiPMIdx = chan.localSiPMIdx() + FTConstants::nSiPM * globalMatIdx;

      ++m_unbiasedFTResidual[node.unbiasedResidual()];
      ++m_biasedFTResidual[node.residual()];
      ++m_globalSiPMIdx[globalSiPMIdx];

      if ( node.errResidual2() > TrackParameters::lowTolerance ) {
        double residualRms = node.residual() * std::sqrt( node.errMeasure2() / node.errResidual2() );
        m_RMSResidualModules[uniqueModule] += residualRms;
        m_RMSResidualQuarters[uniqueQuarter] += residualRms;
      }

      m_UnbiasedResidualModules[uniqueModule] += node.unbiasedResidual();
      m_UnbiasedResidualQuarters[uniqueQuarter] += node.unbiasedResidual();

      ++m_unbiasedFTResidualLayer[{uniqueLayer, node.unbiasedResidual()}];
      ++m_biasedFTResidualLayer[{uniqueLayer, node.residual()}];
      ++m_residualPerModule[{uniqueModule, node.residual()}];

      // plots per module on request
      if ( m_verboseMode.value() ) {
        auto& histos_modules = m_histograms_modules.at( chan.globalModuleIdx() );
        ++histos_modules.UnbiasedResidualModules_permodule[node.unbiasedResidual()];
      }

      // plots per station
      auto& histos = m_histograms.at( chan.station() );

      if ( node.errResidual2() > TrackParameters::lowTolerance ) {
        double residualRms = node.residual() * std::sqrt( node.errMeasure2() / node.errResidual2() );
        histos.RMSResidualModules[uniqueModuleInStation] += residualRms;
      }
      histos.UnbiasedResidualModules[uniqueModuleInStation] += node.residual();

      LHCb::StateVector aState;
      extrapolator()->propagate( track, m_refVec[station - 1], aState, geometry ).ignore();

      ++histos.x[aState.x()];
      ++histos.y[aState.y()];
      ++histos.Tx[aState.tx()];
      ++histos.Ty[aState.ty()];
      ++histos.Xdist[aState.x()];
      ++histos.XdistWide[aState.x()];
      ++histos.Ydist[aState.y()];
      ++histos.pos[{aState.x(), aState.y()}];
      ++histos.slopes[{aState.tx(), aState.ty()}];
      if ( m_expertMode.value() ) {
        // Fill ntuples in expert mode only

        Tuple nodeTuple = nTuple( "FTTrackTuple_nodes", "" );

        // Node properties
        nodeTuple->column( "p", p ).ignore();
        nodeTuple->column( "pt", pt ).ignore();
        nodeTuple->column( "px", px ).ignore();
        nodeTuple->column( "py", py ).ignore();
        nodeTuple->column( "pz", pz ).ignore();
        nodeTuple->column( "tx", tx ).ignore();
        nodeTuple->column( "ty", ty ).ignore();

        // Node errors
        if constexpr ( std::is_same<TNode, LHCb::FitNode>::value ) {
          node.template visit( [&]( auto& n ) {
            for ( int i = 0; i < n.typedim; i++ ) {
              nodeTuple->column( "Error" + std::to_string( i ), n.errMeasure()[i] ).ignore();
            }
          } );
        }
        if constexpr ( std::is_same<TNode, LHCb::Pr::Tracks::Fit::Node>::value ) {
          nodeTuple->column( "Error", node.measurement_error ).ignore();
        }

        // Residuals
        nodeTuple->column( "residual", node.residual() ).ignore();
        nodeTuple->column( "unbiasedResidual", node.unbiasedResidual() ).ignore();
        nodeTuple->column( "errResidual", node.errResidual() ).ignore();
        nodeTuple->column( "errMeasure2", node.errMeasure2() ).ignore();
        nodeTuple->column( "errResidual2", node.errResidual2() ).ignore();

        // SciFi geometry objects
        nodeTuple->column( "station", station ).ignore();
        nodeTuple->column( "globalLayerID", uniqueLayer ).ignore();
        nodeTuple->column( "globalQuarterIndex", uniqueQuarter ).ignore();
        nodeTuple->column( "globalModuleIndex", uniqueModule ).ignore();
        nodeTuple->column( "globalMatID", chan.globalMatID() ).ignore();
        nodeTuple->column( "sipmInMatRaw", chan.sipm() ).ignore();
        nodeTuple->column( "channelInSiPMraw", chan.channel() ).ignore();

        nodeTuple->column( "uniqueModuleInStation", uniqueModuleInStation ).ignore();
        nodeTuple->column( "localMatIdxInQuarter", chan.localMatIdx_quarter() ).ignore();
        nodeTuple->column( "localSiPMIdxInQuarter", chan.localSiPMIdx_quarter() ).ignore();
        nodeTuple->column( "localChannelIdxInModule", chan.localChannelIdx_quarter() ).ignore();

        // Global coordinates for node
        const Gaudi::XYZPoint nodeGlobal( state( node ).position().x(), state( node ).position().y(),
                                          state( node ).position().z() );
        LHCb::State           unbiasedState;
        if constexpr ( std::is_same<TNode, LHCb::FitNode>::value ) {
          unbiasedState = node.template visit( [&]( auto& n ) { return n.unbiasedState( node ); } );
        }
        if constexpr ( std::is_same<TNode, LHCb::Pr::Tracks::Fit::Node>::value ) {
          unbiasedState = node.unbiasedState();
        }
        const Gaudi::XYZPoint node_unbiasedGlobal( unbiasedState.x(), unbiasedState.y(), unbiasedState.z() );

        nodeTuple->column( "node_", nodeGlobal ).ignore();
        nodeTuple->column( "unbiasedNode_", node_unbiasedGlobal ).ignore();

        nodeTuple->write().ignore();
      }
    }
  }
}
