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

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/FitNode.h"
#include "Event/PrFitNode.h"
#include "Event/PrKalmanFitResult.h"
#include "Event/Track.h"
#include "Event/TrackFitResult.h"
#include "Event/VPLightCluster.h"
#include "Gaudi/Accumulators/Histogram.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "LHCbAlgs/Consumer.h"
#include "TrackInterfaces/IDetailedHitExpectation.h"
#include "TrackKernel/TrackFunctors.h"
#include "VPDet/DeVP.h"
#include "fmt/format.h"

//=============================================================================
// Anonymous Namespace
//=============================================================================

namespace {
  enum HitType { VPX = 0, VPY, VP2D, UT, FT, Muon, HitTypeUnknown };

  const std::vector<std::string> HitTypeName{"VPX", "VPY", "VP2D", "UT", "FT", "Muon"};

  template <typename TNode>
  inline HitType hittypemap( const TNode& node ) {
    if ( node.isMuon() )
      return HitType::Muon;
    else if ( node.isFT() )
      return HitType::FT;
    else if ( node.isUT() )
      return HitType::UT;
    else if ( node.isVP() ) {
      if constexpr ( std::is_same<TNode, LHCb::FitNode>::value ) {
        const LHCb::Measurement& meas = node.measurement();
        return meas.visit(
            []( const LHCb::Measurement::VP& vp ) {
              return vp.projection() == LHCb::Measurement::VP::Projection::X ? VPX : VPY;
            },
            []( const LHCb::Measurement::VP2D& ) { return HitType::VP2D; },
            []( ... ) { return HitType::HitTypeUnknown; } );
      }
      if constexpr ( std::is_same<TNode, LHCb::Pr::Tracks::Fit::Node>::value ) {
        return node.measurement_dir[0] == 0 ? HitType::VPY : HitType::VPX;
      }
    } else
      return HitType::HitTypeUnknown;
  }
} // namespace

namespace LHCb::Tr::Monitor {

  namespace {
    template <typename TFitResult>
    int HitsACside( TFitResult fitResult, const DeVP& det ) {
      // -1: left side only, 0: overlap track, +1: right side only
      int  side         = 0;
      bool allhitsleft  = true;
      bool allhitsright = true;

      const auto& nodes_ = nodes( *fitResult );
      for ( const auto& node : nodes_ ) {
        if ( !node.hasMeasurement() ) continue;
        if ( !node.isHitOnTrack() ) continue;
        if ( !node.isVP() ) continue;
        Detector::VPChannelID chan = id( node ).vpID();
        const DeVPSensor&     sens = det.sensor( chan );
        allhitsleft                = allhitsleft && sens.isLeft();
        allhitsright               = allhitsright && sens.isRight();
      }
      if ( allhitsleft ) side = -1;
      if ( allhitsright ) side = 1;
      return side;
    }

  } // namespace

  enum struct OutputMode { Histograms, Tuples };
  template <typename TFitResult, typename TNode, OutputMode outmode>
  struct VPTrackMonitor final
      : LHCb::Algorithm::Consumer<void( LHCb::Track::Range const&, VPLightClusters const&,
                                        IDetailedHitExpectation const&, DeVP const&, DetectorElement const& ),
                                  DetDesc::usesBaseAndConditions<GaudiTupleAlg, DeVP, DetectorElement>> {
    Gaudi::Property<int>  TrackEnum{this, "TrackEnum", 0, "Track type"};
    Gaudi::Property<bool> ResidualDebug{this, "ResidualDebug", 0,
                                        "Swithcing on the option to produce extra plots for residuals"};

    /// Standard constructor
    VPTrackMonitor( const std::string& name, ISvcLocator* pSvcLocator )
        : Consumer{name,
                   pSvcLocator,
                   {KeyValue{"TrackContainer", ""}, KeyValue{"ClusterContainer", VPClusterLocation::Light},
                    KeyValue{"VPHitExpectation", "VPHitExpectation"},
                    KeyValue{"VPDetectorLocation", DeVPLocation::Default},
                    KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}}} {}

    struct Histos {
      // Histograms filled per track
      mutable Gaudi::Accumulators::Histogram<1> m_phi;
      mutable Gaudi::Accumulators::Histogram<1> m_eta;
      mutable Gaudi::Accumulators::Histogram<1> m_theta;
      mutable Gaudi::Accumulators::Histogram<1> m_chi2NDOF;
      mutable Gaudi::Accumulators::Histogram<1> m_probChi2;

      mutable Gaudi::Accumulators::Histogram<1>        m_nClusters;
      mutable Gaudi::Accumulators::Histogram<1>        m_pseudoEfficiency;
      mutable Gaudi::Accumulators::Histogram<1>        m_nVPHits;
      mutable Gaudi::Accumulators::Histogram<1>        m_p;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_theta_nVPHits;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_eta_nVPHits;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_phi_nVPHits;

      // Histograms filled per cluster
      mutable Gaudi::Accumulators::Histogram<1>        m_sensor;
      mutable Gaudi::Accumulators::Histogram<1>        m_module;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_pseudoEfficiency_sens;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_pseudoEfficiency_sens_inner;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_pseudoEfficiency_sens_outer;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_nHit_sens;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_nHit_mod;
      mutable Gaudi::Accumulators::Histogram<2>        m_sensor_res;
      mutable Gaudi::Accumulators::Histogram<2>        m_sensor_resPull;
      mutable Gaudi::Accumulators::Histogram<2>        m_module_res;
      mutable Gaudi::Accumulators::Histogram<2>        m_module_resPull;
      mutable Gaudi::Accumulators::Histogram<2>        m_sensor_resX;
      mutable Gaudi::Accumulators::Histogram<2>        m_sensor_resXPull;
      mutable Gaudi::Accumulators::Histogram<2>        m_module_resX;
      mutable Gaudi::Accumulators::Histogram<2>        m_module_resXPull;
      mutable Gaudi::Accumulators::Histogram<2>        m_sensor_resY;
      mutable Gaudi::Accumulators::Histogram<2>        m_sensor_resYPull;
      mutable Gaudi::Accumulators::Histogram<2>        m_module_resY;
      mutable Gaudi::Accumulators::Histogram<2>        m_module_resYPull;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_res_station;
      mutable Gaudi::Accumulators::ProfileHistogram<2> m_prof_resX_station_tx;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resX_station25_tx;
      mutable Gaudi::Accumulators::ProfileHistogram<2> m_prof_resY_station_ty;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resY_station25_ty;

      mutable Gaudi::Accumulators::Histogram<1> m_nTracksperEvent;

      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resXLocal_Vs_qXLocal;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resYLocal_Vs_qYLocal;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resXGlobal_Vs_qXGlobal;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resYGlobal_Vs_qYGlobal;

      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resXGlobal_Vs_qXLocal;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resYGlobal_Vs_qYLocal;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resXLocal_Vs_qXGlobal;
      mutable Gaudi::Accumulators::ProfileHistogram<1> m_prof_resYLocal_Vs_qYGlobal;

      Histos( const VPTrackMonitor* owner, std::string const& fwdbwd )
          : // Histograms filled per track
          m_phi{owner, fwdbwd + "/phi", "phi", {200, -3.2, 3.2}}
          , m_eta{owner, fwdbwd + "/pseudoRapidity", "pseudoRapidity", {400, -6., 6.}}
          , m_theta{owner, fwdbwd + "/theta", "theta", {200, 0., 700.}}
          , m_chi2NDOF{owner, fwdbwd + "/chi2NDOF", "chi2NDOF", {50, 0., 20.}}
          , m_probChi2{owner, fwdbwd + "/probChi2", "probChi2", {50, -0.1, 1.1}}

          , m_nClusters{owner, fwdbwd + "/Clusters", "#clusters", {71, -0.5, 70.5}}
          , m_pseudoEfficiency{owner, fwdbwd + "/PseudoEfficiency", "Track Pseudoefficiency", {50, -0.1, 1.1}}

          , m_nVPHits{owner, fwdbwd + "/VPHits", "#VP hits", {27, -0.5, 26.5}}
          , m_p{owner, fwdbwd + "/p", "p", {200, 0., 100.}}
          , m_theta_nVPHits{owner, fwdbwd + "/VPHitsPerTheta", "#VP hits vs theta [mrad]", {200, 0., 700.}}
          , m_eta_nVPHits{owner, fwdbwd + "/VPHitsPerEta", "#VP hits vs eta", {200, 0., 6.}}
          , m_phi_nVPHits{owner, fwdbwd + "/VPHitsPerPhi", "#VP hits vs phi", {200, -3.2, 3.2}}

          // Histograms filled per cluster
          , m_sensor{owner, fwdbwd + "/Sensors", "Sensors", {208, -0.5, 207.5}}
          , m_module{owner, fwdbwd + "/Modules", "Modules", {52, -0.5, 51.5}}
          , m_pseudoEfficiency_sens{owner,
                                    fwdbwd + "/Pseudoefficiency_sens",
                                    "Pseudoefficiency of sensors",
                                    {208, -.5, 207.5}}
          , m_pseudoEfficiency_sens_inner{owner,
                                          fwdbwd + "/Pseudoefficiency_sens_inner",
                                          "Pseudoefficiency of inner sensors",
                                          {208, -.5, 207.5}}
          , m_pseudoEfficiency_sens_outer{owner,
                                          fwdbwd + "/Pseudoefficiency_sens_outer",
                                          "Pseudoefficiency of outer sensors",
                                          {208, -.5, 207.5}}
          , m_nHit_sens{owner, fwdbwd + "/N_nhit_sens", "Hit Multiplicity on sensors", {208, -.5, 207.5}}
          , m_nHit_mod{owner, fwdbwd + "/N_nhit_mod", "Hit Multiplicity on modules", {52, -.5, 51.5}}
          , m_sensor_res{owner,
                         fwdbwd + "/BiasedResidualSensor",
                         "Biased Residual per Sensor",
                         {208, -.5, 207.5},
                         {200, -.1, .1}}
          , m_sensor_resPull{owner,
                             fwdbwd + "/BiasedResidualPullSensor",
                             "Biased Residual / Error per Sensor",
                             {208, -.5, 207.5},
                             {200, -5., 5.}}
          , m_module_res{owner,
                         fwdbwd + "/BiasedResidualModule",
                         "Biased Residual per Module",
                         {52, -.5, 51.5},
                         {200, -.1, .1}}
          , m_module_resPull{owner,
                             fwdbwd + "/BiasedResidualPullModule",
                             "Biased Residual / Error per Module",
                             {52, -.5, 51.5},
                             {200, -5., 5.}}
          , m_sensor_resX{owner,
                          fwdbwd + "/BiasedResidualXSensor",
                          "Biased Residual X per Sensor",
                          {208, -.5, 207.5},
                          {200, -.1, .1}}
          , m_sensor_resXPull{owner,
                              fwdbwd + "/BiasedResidualXPullSensor",
                              "Biased Residual X / Error per Sensor",
                              {208, -.5, 207.5},
                              {200, -5., 5.}}
          , m_module_resX{owner,
                          fwdbwd + "/BiasedResidualXModule",
                          "Biased Residual X per Module",
                          {52, -.5, 51.5},
                          {200, -.1, .1}}
          , m_module_resXPull{owner,
                              fwdbwd + "/BiasedResidualXPullModule",
                              "Biased Residual X/ Error per Module",
                              {52, -.5, 51.5},
                              {200, -5., 5.}}
          , m_sensor_resY{owner,
                          fwdbwd + "/BiasedResidualYSensor",
                          "Biased Residual Y per Sensor",
                          {208, -.5, 207.5},
                          {200, -.1, .1}}
          , m_sensor_resYPull{owner,
                              fwdbwd + "/BiasedResidualYPullSensor",
                              "Biased Residual Y / Error per Sensor",
                              {208, -.5, 207.5},
                              {200, -5., 5.}}
          , m_module_resY{owner,
                          fwdbwd + "/BiasedResidualYModule",
                          "Biased Residual Y per Module",
                          {52, -.5, 51.5},
                          {200, -.1, .1}}
          , m_module_resYPull{owner,
                              fwdbwd + "/BiasedResidualYPullModule",
                              "Biased Residual Y / Error per Module",
                              {52, -.5, 51.5},
                              {200, -5., 5.}}
          , m_prof_res_station{owner, fwdbwd + "/ProfResVsStation", "ProfResVsStation", {26, -.5, 25.5}}
          , m_prof_resX_station_tx{owner,
                                   fwdbwd + "/ProfResXVsStationTx",
                                   "ProfResXVsStationTx",
                                   {26, -.5, 25.5},
                                   {100, -0.3, 0.3}}
          , m_prof_resX_station25_tx{owner,
                                     fwdbwd + "/ProfResXVsStation25Tx",
                                     "ProfResXVsStation25Tx",
                                     {100, -0.3, 0.3}}
          , m_prof_resY_station_ty{owner,
                                   fwdbwd + "/ProfResYVsStationTy",
                                   "ProfResYVsStationTy",
                                   {26, -.5, 25.5},
                                   {100, -0.3, 0.3}}
          , m_prof_resY_station25_ty{owner,
                                     fwdbwd + "/ProfResYVsStation25Ty",
                                     "ProfResYVsStation25Ty",
                                     {100, -0.3, 0.3}}

          , m_nTracksperEvent{owner, fwdbwd + "/NtracksperEvent", "nTracks Vs event", {351, -0.5, 350.5}}

          , m_prof_resXLocal_Vs_qXLocal{owner, fwdbwd + "/LocalResX_Vs_LocalQx", "LocalResX Vs LocalQx", {100, -5, 50}}
          , m_prof_resYLocal_Vs_qYLocal{owner, fwdbwd + "/LocalResY_Vs_LocalQy", "LocalResY Vs LocalQy", {100, -5, 20}}
          , m_prof_resXGlobal_Vs_qXGlobal{owner,
                                          fwdbwd + "/GlobalResX_Vs_GlobalQx",
                                          "GlobalResX Vs GlobalQx",
                                          {300, -80, 80}}
          , m_prof_resYGlobal_Vs_qYGlobal{owner,
                                          fwdbwd + "/GlobalResY_Vs_GlobalQy",
                                          "GlobalResY Vs GlobalQy",
                                          {300, -80, 80}}

          , m_prof_resXGlobal_Vs_qXLocal{owner,
                                         fwdbwd + "/GlobalResX_Vs_LocalQx",
                                         "GlobalResX Vs LocalQx",
                                         {100, -5, 50}}
          , m_prof_resYGlobal_Vs_qYLocal{owner,
                                         fwdbwd + "/GlobalResY_Vs_LocalQy",
                                         "GlobalResY Vs LocalQy",
                                         {100, -5, 20}}
          , m_prof_resXLocal_Vs_qXGlobal{owner,
                                         fwdbwd + "/LocalResX_Vs_GlobalQx",
                                         "LocalResX Vs GlobalQx",
                                         {300, -80, 80}}
          , m_prof_resYLocal_Vs_qYGlobal{
                owner, fwdbwd + "/LocalResY_Vs_GlobalQy", "LocalResY Vs GlobalQy", {300, -60, 60}} {}
    };

    std::array<Histos, 8> m_histograms{
        Histos( this, "Forward/TracksAside" ),   Histos( this, "Backward/TracksAside" ),
        Histos( this, "Forward/TracksOverlap" ), Histos( this, "Backward/TracksOverlap" ),
        Histos( this, "Forward/TracksCside" ),   Histos( this, "Backward/TracksCside" ),
        Histos( this, "Forward/TracksAll" ),     Histos( this, "Backward/TracksAll" )};

    struct optionalHistos {
      // Histograms filled per track

      mutable std::array<Gaudi::Accumulators::ProfileHistogram<2>, 208> m_2Dprof_resXLocal_Vs_XYLocal;
      mutable std::array<Gaudi::Accumulators::ProfileHistogram<2>, 52>  m_2Dprof_resXGlobal_Vs_XYGlobal;
      mutable std::array<Gaudi::Accumulators::ProfileHistogram<2>, 208> m_2Dprof_resYLocal_Vs_XYLocal;
      mutable std::array<Gaudi::Accumulators::ProfileHistogram<2>, 52>  m_2Dprof_resYGlobal_Vs_XYGlobal;

      mutable std::array<Gaudi::Accumulators::Histogram<1>, 52> m_res1D_module;
      mutable std::array<Gaudi::Accumulators::Histogram<1>, 52> m_resX1D_module;
      mutable std::array<Gaudi::Accumulators::Histogram<1>, 52> m_resY1D_module;

      template <std::size_t... IDXs>
      static std::array<Gaudi::Accumulators::ProfileHistogram<2>, sizeof...( IDXs )>
      profhisto2DArrayBuilder( const VPTrackMonitor* owner, std::string const& fwdbwd, const std::string& name,
                               const std::string& title, std::tuple<unsigned, double, double> xbins,
                               std::tuple<unsigned, double, double> ybins, std::index_sequence<IDXs...> ) {
        return {{{owner,
                  fwdbwd + name + std::to_string( IDXs ),
                  title + std::to_string( IDXs ),
                  {std::get<0>( xbins ), std::get<1>( xbins ), std::get<2>( xbins )},
                  {std::get<0>( ybins ), std::get<1>( ybins ), std::get<2>( ybins )}}...}};
      }

      template <std::size_t... IDXs>
      static std::array<Gaudi::Accumulators::Histogram<1>, sizeof...( IDXs )>
      histo1DArrayBuilder( const VPTrackMonitor* owner, std::string const& fwdbwd, const std::string& name,
                           const std::string& title, std::tuple<unsigned, double, double> xbins,
                           std::index_sequence<IDXs...> ) {
        return {{{owner,
                  fwdbwd + name + std::to_string( IDXs ),
                  title + std::to_string( IDXs ),
                  {std::get<0>( xbins ), std::get<1>( xbins ), std::get<2>( xbins )}}...}};
      }

      optionalHistos( const VPTrackMonitor* owner, std::string const& fwdbwd )
          : // Histograms filled per track

          m_2Dprof_resXLocal_Vs_XYLocal{profhisto2DArrayBuilder( owner, fwdbwd, "/LocalResX_Vs_LocalXY_Sens",
                                                                 "LocalResX Vs LocalXY", {100, -5., 50.},
                                                                 {100, -5., 20.}, std::make_index_sequence<208>() )}
          , m_2Dprof_resXGlobal_Vs_XYGlobal{profhisto2DArrayBuilder( owner, fwdbwd, "/GlobalResX_Vs_GlobalXY_Mod",
                                                                     "GlobalResX Vs GlobalXY", {300, -80., 80.},
                                                                     {300, -80., 80.}, std::make_index_sequence<52>() )}
          , m_2Dprof_resYLocal_Vs_XYLocal{profhisto2DArrayBuilder( owner, fwdbwd, "/LocalResY_Vs_LocalXY_Sens",
                                                                   "LocalResYVs LocalXY", {100, -5., 50.},
                                                                   {100, -5., 20.}, std::make_index_sequence<208>() )}
          , m_2Dprof_resYGlobal_Vs_XYGlobal{profhisto2DArrayBuilder( owner, fwdbwd, "/GlobalResY_Vs_GlobalXY_Mod",
                                                                     "GlobalResY Vs GlobalXY", {300, -80., 80.},
                                                                     {300, -80., 80.}, std::make_index_sequence<52>() )}
          , m_res1D_module{histo1DArrayBuilder( owner, fwdbwd, "/BiasedResidual1DMod", "Biased Residual Module ",
                                                {200, -0.1, 0.1}, std::make_index_sequence<52>() )}
          , m_resX1D_module{histo1DArrayBuilder( owner, fwdbwd, "/BiasedResidualX1DMod", "Biased X Residual Module ",
                                                 {200, -0.1, 0.1}, std::make_index_sequence<52>() )}
          , m_resY1D_module{histo1DArrayBuilder( owner, fwdbwd, "/BiasedResidualY1DMod", "Biased Y Residual Module ",
                                                 {200, -0.1, 0.1}, std::make_index_sequence<52>() )} {}
    };

    std::array<std::unique_ptr<optionalHistos>, 8> m_optionalhistograms;

    StatusCode initialize() override {
      return Consumer::initialize().andThen( [&] {
        if ( ResidualDebug ) {

          m_optionalhistograms[0] = std::make_unique<optionalHistos>( this, "Forward/TracksAside" );
          m_optionalhistograms[1] = std::make_unique<optionalHistos>( this, "Backward/TracksAside" );
          m_optionalhistograms[2] = std::make_unique<optionalHistos>( this, "Forward/TracksOverlap" );
          m_optionalhistograms[3] = std::make_unique<optionalHistos>( this, "Backward/TracksOverlap" );
          m_optionalhistograms[4] = std::make_unique<optionalHistos>( this, "Forward/TracksCside" );
          m_optionalhistograms[5] = std::make_unique<optionalHistos>( this, "Backward/TracksCside" );
          m_optionalhistograms[6] = std::make_unique<optionalHistos>( this, "Forward/TracksAll" );
          m_optionalhistograms[7] = std::make_unique<optionalHistos>( this, "Backward/TracksAll" );
        }
      } );
    }

    ///< Algorithm execution
    void operator()( const LHCb::Track::Range& tracks, const VPLightClusters& clusters,
                     IDetailedHitExpectation const& expectTool, const DeVP& det,
                     const DetectorElement& lhcb ) const override {

      const Gaudi::XYZPoint origin( 0., 0., 0. );

      size_t tracknumber = 0;

      auto get_cluster_ptr = [&clusters]( const auto& node ) {
        auto i = std::find_if( begin( clusters ), end( clusters ),
                               [id = id( node ).vpID()]( const auto& c ) { return c.channelID() == id; } );
        return i != end( clusters ) ? &*i : nullptr;
      };

      std::array<std::array<std::size_t, 208>, 8> N_exp{{{0}, {0}}};
      std::array<std::array<std::size_t, 208>, 8> N_rec{{{0}, {0}}};
      std::array<std::array<std::size_t, 52>, 8>  N_rec_mod{{{0}, {0}}};

      std::vector<int> nTracks( 8, 0 );

      auto& geometry = *lhcb.geometry();

      // start of track-loop
      for ( const LHCb::Track* track : tracks ) {
        const bool bwd  = track->isVeloBackward();
        const auto type = track->type();

        if constexpr ( outmode == OutputMode::Histograms ) {
          // When making histograms (Velo monitoring), skip tracks which have no hits in the VELO
          if ( !track->hasVelo() ) { continue; }
          if ( TrackEnum == 3 ) {
            if ( type != LHCb::Track::Types::Long ) { continue; }
          }
        } else {
          // When making tuples (alignment), skip tracks which are not VELO tracks and not long tracks and not backward
          // tracks
          if ( type != LHCb::Track::Types::Velo && type != LHCb::Track::Types::Long && !bwd ) { continue; }
        }

        const bool fitted = track->checkFitStatus( LHCb::Track::FitStatus::Fitted );
        if ( !fitted ) continue;

        const double     chi2NDOF    = track->chi2PerDoF();
        const double     probChi2    = track->probChi2();
        const double     ghostProb   = track->ghostProbability();
        const double     phi         = track->phi();
        const double     eta         = bwd ? ( -1.0 * track->pseudoRapidity() ) : track->pseudoRapidity();
        const double     p           = ( track->p() ) / 1000.0;                    // p in GeV
        const double     pt          = ( track->pt() ) / 1000.0;                   // p in GeV
        const double     theta       = 2. * std::atan( std::exp( -eta ) ) * 1000.; // theta in mrad
        const auto       momentumVec = track->momentum();
        const double     px          = momentumVec.X();
        const double     py          = momentumVec.Y();
        const double     pz          = momentumVec.Z();
        Gaudi::XYZVector slopes      = track->slopes();
        const double     tx          = slopes.X();
        const double     ty          = slopes.Y();
        const int        charge      = track->charge();

        const std::vector<LHCb::LHCbID>& ids = track->lhcbIDs();
        const auto                       nVPHits =
            std::count_if( ids.begin(), ids.end(), []( const LHCb::LHCbID& id ) { return id.isVP(); } );

        const auto   info          = expectTool.detailedExpectation( *track, geometry );
        const double nExpectedHits = info.size();
        const double nFoundHits =
            std::count_if( info.begin(), info.end(), []( const auto& element ) { return element.found; } );
        const double PseudoEfficiency = nFoundHits / double( nExpectedHits );

        const auto fitResult = dynamic_cast<const TFitResult*>( track->fitResult() );
        if ( !fitResult ) {
          Warning( "Empty fit result. Track might be using wrong fit result type.", StatusCode::SUCCESS ).ignore();
          continue;
        }

        auto         hitsACside = HitsACside( fitResult, det );
        const auto&  nodes_     = nodes( *fitResult );
        unsigned int nClusters  = nodes_.size();
        size_t       nodenumber = 0;

        for ( int hitsSide = 0; hitsSide <= 1; hitsSide++ ) {
          int sides = 2; // 2 for All; -1 for Aside; 0 for Overlap; 1 for C side
          if ( hitsSide == 1 ) { sides = hitsACside; }

          nTracks[2 * ( sides + 1 ) + int( bwd )]++;

          for ( const auto& element : info ) {
            const LHCb::LHCbID          lhcbid = element.id;
            const Detector::VPChannelID chanid( lhcbid.vpID() );
            const unsigned              sensorNumber = static_cast<unsigned>( chanid.sensor() );

            N_exp[2 * ( sides + 1 ) + int( bwd )][sensorNumber]++;
            N_rec[2 * ( sides + 1 ) + int( bwd )][sensorNumber] += element.found;
            N_rec_mod[2 * ( sides + 1 ) + int( bwd )][element.global_id] += element.found;
          }
        }

        auto velonodenumber = std::count_if( nodes_.begin(), nodes_.end(), []( const TNode& node ) {
          return node.hasMeasurement() && node.isHitOnTrack() && node.isVP();
        } );

        for ( const auto& node : nodes_ ) { // start cluster loop
          if ( !node.hasMeasurement() ) continue;
          if ( !node.isHitOnTrack() ) continue;
          // Skip non-VP measurements.
          if ( !node.isVP() ) continue;

          auto iclus  = get_cluster_ptr( node );
          auto isXRes = hittypemap( node ) == VPX;
          auto isYRes = hittypemap( node ) == VPY;

          // Get the channel.
          Detector::VPChannelID chan = id( node ).vpID();
          // Get the sensor.
          const DeVPSensor& sens = det.sensor( chan );

          const auto            corner = sens.localToGlobal( origin );
          const Gaudi::XYZPoint cluGlobal( iclus->x(), iclus->y(), iclus->z() );
          const Gaudi::XYZPoint nodeGlobal( state( node ).position().x(), state( node ).position().y(),
                                            state( node ).position().z() );

          LHCb::State unbiasedState;
          if constexpr ( std::is_same<TNode, LHCb::FitNode>::value ) {
            unbiasedState = node.template visit( [&]( auto& n ) { return n.unbiasedState( node ); } );
          }
          if constexpr ( std::is_same<TNode, LHCb::Pr::Tracks::Fit::Node>::value ) {
            unbiasedState = node.unbiasedState();
          }
          const Gaudi::XYZPoint node_unbiasedGlobal( unbiasedState.x(), unbiasedState.y(), unbiasedState.z() );
          const auto            node_unbiasedLocal = sens.globalToLocal( node_unbiasedGlobal );
          const int             sensor             = to_unsigned( chan.sensor() );
          const int             module             = sens.module();
          const int             station            = sens.station();

          if constexpr ( outmode == OutputMode::Histograms ) {

            for ( int hitsSide = 0; hitsSide <= 1; hitsSide++ ) {
              int sides = 2; // 2 for All; -1 for Aside; 0 for Overlap; 1 for C side
              if ( hitsSide == 1 ) { sides = hitsACside; }

              auto& histos = m_histograms[2 * ( sides + 1 ) + int( bwd )];

              ++histos.m_sensor[sensor];
              ++histos.m_module[module];

              ++histos.m_sensor_res[{sensor, node.residual()}];
              ++histos.m_sensor_resPull[{sensor, node.residual() / node.errResidual()}];

              ++histos.m_module_res[{module, node.residual()}];
              ++histos.m_module_resPull[{module, node.residual() / node.errResidual()}];

              if ( isXRes ) {

                ++histos.m_sensor_resX[{sensor, node.residual()}];
                ++histos.m_sensor_resXPull[{sensor, node.residual() / node.errResidual()}];

                ++histos.m_module_resX[{module, node.residual()}];
                ++histos.m_module_resXPull[{module, node.residual() / node.errResidual()}];
              } else {

                ++histos.m_sensor_resY[{sensor, node.residual()}];
                ++histos.m_sensor_resYPull[{sensor, node.residual() / node.errResidual()}];

                ++histos.m_module_resY[{module, node.residual()}];
                ++histos.m_module_resYPull[{module, node.residual() / node.errResidual()}];
              }

              const Gaudi::XYZPoint q( state( node ).position().x(), state( node ).position().y(),
                                       state( node ).position().z() );

              const auto ql = sens.globalToLocal( q );

              if ( isXRes ) {
                histos.m_prof_resXLocal_Vs_qXLocal[ql.x()] += node.residual();
                histos.m_prof_resXGlobal_Vs_qXGlobal[q.x()] += ( nodeGlobal.x() - cluGlobal.x() );

                histos.m_prof_resXLocal_Vs_qXGlobal[q.x()] += node.residual();
                histos.m_prof_resXGlobal_Vs_qXLocal[ql.x()] += ( nodeGlobal.x() - cluGlobal.x() );

              } else {
                histos.m_prof_resYLocal_Vs_qYLocal[ql.y()] += node.residual();
                histos.m_prof_resYGlobal_Vs_qYGlobal[q.y()] += ( nodeGlobal.y() - cluGlobal.y() );

                histos.m_prof_resYLocal_Vs_qYGlobal[q.y()] += node.residual();
                histos.m_prof_resYGlobal_Vs_qYLocal[ql.y()] += ( nodeGlobal.y() - cluGlobal.y() );
              }

              if ( ResidualDebug ) {
                auto& optionalhistos = *m_optionalhistograms[2 * ( sides + 1 ) + int( bwd )].get();
                ++optionalhistos.m_res1D_module[module][node.residual()];
                if ( isXRes ) {
                  ++optionalhistos.m_resX1D_module[module][node.residual()];
                  optionalhistos.m_2Dprof_resXLocal_Vs_XYLocal[sensor][{ql.x(), ql.y()}] += node.residual();
                  optionalhistos.m_2Dprof_resXGlobal_Vs_XYGlobal[module][{q.x(), q.y()}] +=
                      ( nodeGlobal.x() - cluGlobal.x() );
                } else {
                  ++optionalhistos.m_resY1D_module[module][node.residual()];
                  optionalhistos.m_2Dprof_resYLocal_Vs_XYLocal[sensor][{ql.x(), ql.y()}] += node.residual();
                  optionalhistos.m_2Dprof_resYGlobal_Vs_XYGlobal[module][{q.x(), q.y()}] +=
                      ( nodeGlobal.y() - cluGlobal.y() );
                }
              }

              histos.m_prof_res_station[station] += node.residual();
              if ( isXRes ) {
                histos.m_prof_resX_station_tx[{station, tx}] += node.residual();
                if ( station == 25 ) histos.m_prof_resX_station25_tx[tx] += node.residual();
              } else {
                histos.m_prof_resY_station_ty[{station, ty}] += node.residual();
                if ( station == 25 ) histos.m_prof_resY_station25_ty[ty] += node.residual();
              }
            }
          } else {
            Tuple nodeTuple = nTuple( "VPTrackMonitor_nodes", "" );

            nodeTuple->column( "residual", node.residual() ).ignore();
            nodeTuple->column( "unbiasedResidual", node.unbiasedResidual() ).ignore();
            nodeTuple->column( "errResidual", node.errResidual() ).ignore();

            // GLOBAL:
            nodeTuple->column( "clus", cluGlobal ).ignore();
            nodeTuple->column( "node_", nodeGlobal ).ignore();
            nodeTuple->column( "bwd", bwd ).ignore();
            nodeTuple->column( "unbiasedNode_", node_unbiasedGlobal ).ignore();

            // LOCAL:
            nodeTuple->column( "clusLocal_", sens.globalToLocal( cluGlobal ) ).ignore();
            nodeTuple->column( "nodeLocal_", sens.globalToLocal( nodeGlobal ) ).ignore();

            nodeTuple->column( "unbiasedNodeLocal_", node_unbiasedLocal ).ignore();

            nodeTuple->column( "sensEdgeX", corner.x() ).ignore();
            nodeTuple->column( "sensEdgeY", corner.y() ).ignore();
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

            nodeTuple->column( "module", module ).ignore();
            nodeTuple->column( "station", station ).ignore();
            nodeTuple->column( "sensor", sensor ).ignore();
            nodeTuple->column( "isRight", sens.isRight() ).ignore();
            nodeTuple->column( "isLeft", sens.isLeft() ).ignore();
            nodeTuple->column( "nodenumber", static_cast<unsigned long long>( nodenumber ) ).ignore();
            nodeTuple->column( "velonodenumber", static_cast<unsigned long long>( velonodenumber ) ).ignore();
            nodeTuple->column( "eta", eta ).ignore();
            nodeTuple->column( "phi", phi ).ignore();
            nodeTuple->column( "p", p ).ignore();
            nodeTuple->column( "pt", pt ).ignore();
            nodeTuple->column( "px", px ).ignore();
            nodeTuple->column( "py", py ).ignore();
            nodeTuple->column( "pz", pz ).ignore();
            nodeTuple->column( "tx", tx ).ignore();
            nodeTuple->column( "ty", ty ).ignore();
            nodeTuple->column( "probChi2", probChi2 ).ignore();
            nodeTuple->column( "chi2NDOF", chi2NDOF ).ignore();
            nodeTuple->column( "tracknumber", static_cast<unsigned long long>( tracknumber ) ).ignore();
            nodeTuple->column( "tracktype", static_cast<int>( type ) ).ignore();
            nodeTuple->column( "charge", charge ).ignore();
            nodeTuple->column( "isXRes", isXRes ).ignore();
            nodeTuple->column( "isYRes", isYRes ).ignore();
            nodeTuple->column( "hitsACside", hitsACside ).ignore();

            nodeTuple->write().ignore();
          }
          nodenumber++;
        } // end cluster loop

        if constexpr ( outmode == OutputMode::Histograms ) {

          for ( int hitsSide = 0; hitsSide <= 1; hitsSide++ ) {
            int sides = 2; // 2 for All; -1 for Aside; 0 for Overlap; 1 for C side
            if ( hitsSide == 1 ) { sides = hitsACside; }

            auto& histos = m_histograms[2 * ( sides + 1 ) + int( bwd )];

            ++histos.m_phi[phi];
            ++histos.m_eta[eta];

            ++histos.m_theta[theta];
            ++histos.m_pseudoEfficiency[PseudoEfficiency];
            ++histos.m_nClusters[nClusters];
            ++histos.m_nVPHits[nVPHits];

            histos.m_theta_nVPHits[theta] += nVPHits;
            histos.m_eta_nVPHits[eta] += nVPHits;
            histos.m_phi_nVPHits[phi] += nVPHits;

            ++histos.m_p[p];
            ++histos.m_chi2NDOF[chi2NDOF];
            ++histos.m_probChi2[probChi2];
          }
        } else {
          Tuple trackTuple = nTuple( "VPTrackMonitor_tracks", "" );

          trackTuple->column( "phi", phi ).ignore();
          trackTuple->column( "eta", eta ).ignore();

          trackTuple->column( "chi2PerDoF", chi2NDOF ).ignore();
          trackTuple->column( "probChi2", probChi2 ).ignore();
          trackTuple->column( "ghostProb", ghostProb ).ignore();

          trackTuple->column( "nClusters", nClusters ).ignore();
          trackTuple->column( "TrackType", static_cast<int>( type ) ).ignore();
          trackTuple->column( "bwd", bwd ).ignore();

          trackTuple->column( "p", p ).ignore();
          trackTuple->column( "pt", pt ).ignore();
          trackTuple->column( "px", px ).ignore();
          trackTuple->column( "py", py ).ignore();
          trackTuple->column( "pz", pz ).ignore();
          trackTuple->column( "tx", tx ).ignore();
          trackTuple->column( "ty", ty ).ignore();
          trackTuple->column( "tracknumber", static_cast<unsigned long long>( tracknumber ) ).ignore();
          trackTuple->column( "hitsACside", hitsACside ).ignore();
          trackTuple->write().ignore();
        }
        tracknumber++;
      } // end of for loop over tracks

      if constexpr ( outmode == OutputMode::Histograms ) {
        for ( int j = 0; j < 8; j++ ) { // loop for forward and backward tracks

          auto& histos = m_histograms[j];

          ++histos.m_nTracksperEvent[nTracks[j]];

          constexpr unsigned int sensmin = 0, sensmax = 207;
          for ( unsigned int i = sensmin; i <= sensmax; i++ ) {
            if ( N_exp[j][i] == 0 ) continue;
            auto& hist_sens_inner_outer =
                ( ( i % 2 ) == 0 ) ? histos.m_pseudoEfficiency_sens_inner : histos.m_pseudoEfficiency_sens_outer;
            hist_sens_inner_outer[i] += ( N_rec[j][i] / double( N_exp[j][i] ) );
            if ( N_rec[j][i] > 0 ) { histos.m_nHit_sens[i] += N_rec[j][i]; }
            histos.m_pseudoEfficiency_sens[i] += ( N_rec[j][i] / double( N_exp[j][i] ) );
          }

          constexpr unsigned int modmin = 0, modmax = 51;
          for ( unsigned int i = modmin; i <= modmax; i++ ) {
            if ( N_rec_mod[j][i] > 0 ) { histos.m_nHit_mod[i] += N_rec_mod[j][i]; }
          }
        } // end of for loop for different types of tracks
      }   // end of if loop of outmode
    }
  };
  using VPMonitorTuples     = VPTrackMonitor<LHCb::TrackFitResult, LHCb::FitNode, OutputMode::Tuples>;
  using VPMonitorHistograms = VPTrackMonitor<LHCb::TrackFitResult, LHCb::FitNode, OutputMode::Histograms>;
  using VPMonitorTuples_PrKalman =
      VPTrackMonitor<LHCb::PrKalmanFitResult, LHCb::Pr::Tracks::Fit::Node, OutputMode::Tuples>;
  using VPMonitorHistograms_PrKalman =
      VPTrackMonitor<LHCb::PrKalmanFitResult, LHCb::Pr::Tracks::Fit::Node, OutputMode::Histograms>;
  DECLARE_COMPONENT_WITH_ID( VPMonitorTuples, "VPTrackMonitorNT" )
  DECLARE_COMPONENT_WITH_ID( VPMonitorHistograms, "VPTrackMonitor" )
  DECLARE_COMPONENT_WITH_ID( VPMonitorTuples_PrKalman, "VPTrackMonitorNT_PrKalman" )
  DECLARE_COMPONENT_WITH_ID( VPMonitorHistograms_PrKalman, "VPTrackMonitor_PrKalman" )
} // namespace LHCb::Tr::Monitor
