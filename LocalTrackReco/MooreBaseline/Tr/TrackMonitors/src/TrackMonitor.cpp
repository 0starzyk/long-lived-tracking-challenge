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
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/FitNode.h"
#include "Event/Measurement.h"
#include "Event/PrFitNode.h"
#include "Event/PrKalmanFitResult.h"
#include "Event/RecVertex.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "Event/TrackFitResult.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/ToStream.h"
#include "Kernel/HitPattern.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Consumer.h"
#include "Map.h"
#include "TrackInterfaces/IHitExpectation.h"
#include "TrackKernel/TrackFunctors.h"
#include <Gaudi/Accumulators/Histogram.h>
#include <array>
#include <cstddef>
#include <map>
#include <mutex>
#include <string>
#include <tuple>
#include <vector>

using namespace LHCb;
using namespace Gaudi;

/** @class TrackMonitor TrackMonitor.h "TrackCheckers/TrackMonitor"
 *
 * Class for track monitoring
 *  @author M. Needham.
 *  @date   6-5-2007
 */

//=============================================================================
// Anonymous Namespace with some useful things
//=============================================================================

namespace {
  enum HitType { VPX = 0, VPY, VP2D, UT, FT, Muon, NHitTypes };

  const std::vector<std::string> HitTypeName{"VPX", "VPY", "VP2D", "UT", "FT", "Muon"};

  // TODO: not sure if these values still make sense for Run3?
  constexpr auto HitTypeMaxRes = std::array{0.1, 0.1, 0.1, 0.5, 1.0, 10.0};
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
            []( const LHCb::Measurement::VP2D& ) { return HitType::VP2D; }, []( ... ) { return HitType::NHitTypes; } );
      }
      if constexpr ( std::is_same<TNode, LHCb::Pr::Tracks::Fit::Node>::value ) {
        return node.measurement_dir[0] == 0 ? HitType::VPY : HitType::VPX;
      }
    } else
      return HitType::NHitTypes;
  }
} // namespace

//=============================================================================
// TrackMonitor Class Definition
//=============================================================================

class TrackMonitor : public LHCb::Algorithm::Consumer<void( const LHCb::Track::Range& )> {
public:
  TrackMonitor( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator, {KeyValue{"TracksInContainer", LHCb::TrackLocation::Default}} ) {}

  StatusCode initialize() override;

  void operator()( LHCb::Track::Range const& ) const override;

private:
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_unknownFitType{this, "Unknown fit result type."};

  //=============================================================================
  // Struct Definition (TYPE + NAME)
  //=============================================================================

  struct HitHistograms {
    Gaudi::Accumulators::Histogram<1> Residual;
    Gaudi::Accumulators::Histogram<1> ResidualPull;
    HitHistograms( const TrackMonitor* owner, std::string const& prefix, std::string const& name, float residualmax )
        : Residual{owner, prefix + name + "Residual", name + "residual (rms-unbiased)", {50, -residualmax, residualmax}}
        , ResidualPull{owner, prefix + name + "residualPull", name + " residual pull", {50, -5, 5}} {}
  };

  //=============================================================================
  // Struct Definition (TYPE)
  //=============================================================================
  struct TypeTrackHistogram {
    Gaudi::Accumulators::Histogram<1>                     probChi2;
    Gaudi::Accumulators::Histogram<1>                     chi2PerDoF;
    Gaudi::Accumulators::Histogram<1>                     ghostProbability;
    Gaudi::Accumulators::Histogram<1>                     nLHCbIDs;
    Gaudi::Accumulators::Histogram<1>                     pseudoRapidity;
    Gaudi::Accumulators::Histogram<1>                     phi;
    Gaudi::Accumulators::Histogram<1>                     nDoF;
    Gaudi::Accumulators::Histogram<1>                     flag;
    Gaudi::Accumulators::Histogram<1>                     typeHistory;
    Gaudi::Accumulators::Histogram<1>                     fitStatus;
    Gaudi::Accumulators::Histogram<1>                     nMeasurements;
    Gaudi::Accumulators::Histogram<1>                     firststateX;
    Gaudi::Accumulators::Histogram<1>                     firststateY;
    Gaudi::Accumulators::Histogram<1>                     firststateZ;
    Gaudi::Accumulators::Histogram<1>                     firststateTX;
    Gaudi::Accumulators::Histogram<1>                     firststateTY;
    Gaudi::Accumulators::Histogram<1>                     firststateQoverP;
    Gaudi::Accumulators::Histogram<1>                     Z;
    Gaudi::Accumulators::Histogram<1>                     P;
    Gaudi::Accumulators::Histogram<1>                     PT;
    Gaudi::Accumulators::Histogram<1>                     nIter;
    Gaudi::Accumulators::Histogram<1>                     pScatter;
    Gaudi::Accumulators::Histogram<1>                     nMuonHits;
    Gaudi::Accumulators::Histogram<1>                     nUTHits;
    Gaudi::Accumulators::Histogram<1>                     nVPHits;
    Gaudi::Accumulators::Histogram<1>                     nFTHits;
    Gaudi::Accumulators::Histogram<1>                     NumOutliers;
    Gaudi::Accumulators::Histogram<1>                     nVeloALayers;
    Gaudi::Accumulators::Histogram<1>                     nVeloCLayers;
    Gaudi::Accumulators::Histogram<1>                     nVeloOverlapLayers;
    Gaudi::Accumulators::Histogram<1>                     nVeloHoles;
    Gaudi::Accumulators::Histogram<1>                     nFTHoles;
    Gaudi::Accumulators::Histogram<1>                     nFTLayers;
    Gaudi::Accumulators::Histogram<1>                     nUTLayers;
    Gaudi::Accumulators::Histogram<1>                     HitVeloALayers;
    Gaudi::Accumulators::Histogram<1>                     HitVeloCLayers;
    Gaudi::Accumulators::Histogram<1>                     HitUTLayers;
    Gaudi::Accumulators::Histogram<1>                     HitFTLayers;
    Gaudi::Accumulators::Histogram<1>                     TmpChi2PerDoF;
    Gaudi::Accumulators::Histogram<1>                     TmpChi2DownStream;
    Gaudi::Accumulators::Histogram<1>                     TmpChi2Match;
    Gaudi::Accumulators::Histogram<1>                     TmpChi2Muon;
    Gaudi::Accumulators::Histogram<1>                     FirstMStateCov0;
    Gaudi::Accumulators::Histogram<1>                     FirstMStateCov1;
    Gaudi::Accumulators::Histogram<1>                     FirstMStateCov2;
    Gaudi::Accumulators::Histogram<1>                     FirstMStateCov3;
    Gaudi::Accumulators::Histogram<1>                     FirstMStateCov4;
    Gaudi::Accumulators::Histogram<1>                     LastMStateCov0;
    Gaudi::Accumulators::Histogram<1>                     LastMStateCov1;
    Gaudi::Accumulators::Histogram<1>                     LastMStateCov2;
    Gaudi::Accumulators::Histogram<1>                     LastMStateCov3;
    Gaudi::Accumulators::Histogram<1>                     LastMStateCov4;
    Gaudi::Accumulators::Histogram<1>                     Multiplicity;
    Gaudi::Accumulators::StatCounter<>                    MultiplicityCounter;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbVeloVsMom;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbVeloVsPhi;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbDownstreamVsMom;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbDownstreamVsPhi;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbMatchVsMom;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbMatchVsPhi;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbVsMom;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbVsEta;
    Gaudi::Accumulators::ProfileHistogram<1>              Chi2ProbVsPhi;
    std::array<std::unique_ptr<HitHistograms>, NHitTypes> hithistograms;

    TypeTrackHistogram( const TrackMonitor* owner, std::string const& type, std::string const& prefix = "" )
        : probChi2{owner, prefix + "probChi2", "probChi2", {50, 0, 1}}
        , chi2PerDoF{owner, prefix + "chi2_per_ndof", "chi2/ndof", {50, 0, 8}}
        , ghostProbability{owner, prefix + "ghostProb", "ghostProb", {50, 0, 1.0}}
        , nLHCbIDs{owner, prefix + "nLHCBIDs", "#nLHCbIDs", {61, -0.5, 60.5}}
        , pseudoRapidity{owner, prefix + "eta", "eta", {50, 0.95, 6.05}}
        , phi{owner, prefix + "phi", "phi", {50, -M_PI, M_PI}}
        , nDoF{owner, prefix + "ndof", "ndof", {51, -0.5, 50.5}}
        , flag{owner, prefix + "flag", "flag", {256, -0.5, 255.5}}
        , typeHistory{owner,
                      prefix + "history",
                      "history",
                      {int( LHCb::Event::Enum::Track::History::Last ), -0.5,
                       int( LHCb::Event::Enum::Track::History::Last ) - 0.5}}
        , fitStatus{owner, prefix + "fitStatus", "fit status", {6, -0.5, 5.5}}
        , nMeasurements{owner, prefix + "nMeasurements", "#nMeas", {61, -0.5, 60.}}
        , firststateX{owner,
                      prefix + "x_firststate",
                      "x of first state",
                      {50, -100 * Gaudi::Units::mm, 100 * Gaudi::Units::mm}}
        , firststateY{owner,
                      prefix + "y_firststate",
                      "y of first state",
                      {50, -100 * Gaudi::Units::mm, 100 * Gaudi::Units::mm}}
        , firststateZ{owner,
                      prefix + "z_firststate",
                      "z of first state",
                      {50, -200 * Gaudi::Units::mm, 200 * Gaudi::Units::mm}}
        , firststateTX{owner, prefix + "tx_firststate", "tx of first state", {50, -1.0, 1.0}}
        , firststateTY{owner, prefix + "ty_firststate", "ty of first state", {50, -1.0, 1.0}}
        , firststateQoverP{owner, prefix + "qop_firststate", "q/p of first state", {50, -0.001, 0.001}}
        , Z{owner, prefix + "z_closest_tozaxis", "z closest to z-axis", {50, -2000, 2000}}
        , P{owner, prefix + "p", "momentum", {100, 0, 100000 * Gaudi::Units::MeV}}
        , PT{owner, prefix + "pt", "pt", {100, 0, 10000 * Gaudi::Units::MeV}}
        , nIter{owner, prefix + "numiter", "number of fit iterations", {11, -0.5, 10.5}}
        , pScatter{owner,
                   prefix + "pscatter",
                   "momentum used for material corrections",
                   {50, 0, 100000 * Gaudi::Units::MeV}}
        , nMuonHits{owner, prefix + "nMuonHits", "# Muon hits", {21, -0.5, 20.5}}
        , nUTHits{owner, prefix + "nUTHits", "# UT hits", {11, -0.5, 10.5}}
        , nVPHits{owner, prefix + "nVPHits", "# VP hits", {27, -0.5, 26.5}}
        , nFTHits{owner, prefix + "nFTHits", "# FT hits", {16, -0.5, 15.5}}
        , NumOutliers{owner, prefix + "noutliers", "#outliers", {11, -0.5, 10.5}}
        , nVeloALayers{owner, prefix + "nVeloAHitLayers", "# Velo-A hit layers", {26, -0.5, 25.5}}
        , nVeloCLayers{owner, prefix + "nVeloCHitLayers", "# Velo-C hit layers", {26, -0.5, 25.5}}
        , nVeloOverlapLayers{owner, prefix + "nVeloOverlapLayers", "# Velo overlap layers", {26, -0.5, 25.5}}
        , nVeloHoles{owner, prefix + "nVeloHoles", "# Velo holes", {11, -0.5, 10.5}}
        , nFTHoles{owner, prefix + "nFTHoles", "# FT holes", {11, -0.5, 10.5}}
        , nFTLayers{owner, prefix + "nFT", "# FT layers", {11, -0.5, 10.5}}
        , nUTLayers{owner, prefix + "nUT", "# UT layers", {11, -0.5, 10.5}}
        , HitVeloALayers{owner, prefix + "HitVeloALayers", "Hit Velo-A layers", {26, -0.5, 25.5}}
        , HitVeloCLayers{owner, prefix + "HitVeloCLayers", "Hit Velo-C layers", {26, -0.5, 25.5}}
        , HitUTLayers{owner, prefix + "HitUTLayers", "Hit UT layers", {4, -0.5, 3.5}}
        , HitFTLayers{owner, prefix + "HitFTLayers", "Hit FT layers", {12, -0.5, 11.5}}
        , TmpChi2PerDoF{owner, prefix + "chi2PerDofVelo", "chi/dof for velo segment", {50, 0, 5.}}
        , TmpChi2DownStream{owner, prefix + "chi2PerDofDownstream", "chi/dof for T(Muon) segment", {50, 0, 5.}}
        , TmpChi2Match{owner, prefix + "chi2PerDofMatch", "chi/dof upstream-downstream match", {50, 0, 5.}}
        , TmpChi2Muon{owner, prefix + "chi2PerDofMuon", "chi/dof for muon segment", {50, 0, 5.}}
        , FirstMStateCov0{owner, prefix + "xerrorAtFirst", "10log(x error) at first measurement", {50, -3, 2}}
        , FirstMStateCov1{owner, prefix + "yerrorAtFirst", "10log(y error) at first measurement", {50, -3, 2}}
        , FirstMStateCov2{owner, prefix + "zerrorAtFirst", "10log(z error) at first measurement", {50, -7, 0}}
        , FirstMStateCov3{owner, prefix + "tyerrorAtFirst", "10log(ty error) at first measurement", {50, -7, 0}}
        , FirstMStateCov4{owner, prefix + "qoperrorAtFirst", "10log(qop error) at first measurement", {50, -8, 0}}
        , LastMStateCov0{owner, prefix + "xerrorAtLast", "10log(x error) at last measurement", {50, -3, 2}}
        , LastMStateCov1{owner, prefix + "yerrorAtLast", "10log(y error) at last measurement", {50, -3, 2}}
        , LastMStateCov2{owner, prefix + "txerrorAtLast", "10log(tx error) at last measurement", {50, -7, 0}}
        , LastMStateCov3{owner, prefix + "tyerrorAtLast", "10log(ty error) at last measurement", {50, -7, 0}}
        , LastMStateCov4{owner, prefix + "qoperrorAtLast", "10log(qop error) at last measurement", {50, -8, 0}}
        , Multiplicity{owner, prefix + "multiplicity", type + " multiplicity", {101, -1., 201.}}
        , MultiplicityCounter{owner, "#" + type}
        , Chi2ProbVeloVsMom{owner,
                            prefix + "chi2ProbVeloVsMom",
                            "chi2 prob for velo segment versus momentum",
                            {50, 0, 100000 * Gaudi::Units::MeV}}
        , Chi2ProbVeloVsPhi{owner,
                            prefix + "chi2ProbVeloVsPhi",
                            "chi2 prob for velo segment versus phi",
                            {50, -M_PI, M_PI}}
        , Chi2ProbDownstreamVsMom{owner,
                                  prefix + "chi2ProbDownstreamVsMom",
                                  "chi2 prob for T(muon) segment versus momentum",
                                  {50, 0, 100000 * Gaudi::Units::MeV}}
        , Chi2ProbDownstreamVsPhi{owner,
                                  prefix + "chi2ProbDownstreamVsPhi",
                                  "chi2 prob for T(muon) segment versus phi",
                                  {50, -M_PI, M_PI}}
        , Chi2ProbMatchVsMom{owner,
                             prefix + "chi2ProbMatchVsMom",
                             "chi2 prob upstream-downstream match versus momentum",
                             {50, 0, 100000 * Gaudi::Units::MeV}}
        , Chi2ProbMatchVsPhi{owner,
                             prefix + "chi2ProbMatchVsPhi",
                             "chi2 prob upstream-downstream match versus phi",
                             {50, -M_PI, M_PI}}
        , Chi2ProbVsMom{owner,
                        prefix + "chi2ProbVsMom",
                        "chi2 prob versus momentum",
                        {50, 0, 100000 * Gaudi::Units::MeV}}
        , Chi2ProbVsEta{owner, prefix + "chi2ProbVsEta", "chi2 prob versus eta", {30, 2, 5}}
        , Chi2ProbVsPhi{owner, prefix + "chi2ProbVsPhi", "chi2 prob versus phi", {50, -M_PI, M_PI}} {}
  };
  //=============================================================================
  // Struct Definition END
  //=============================================================================

  //=============================================================================
  // Struct Definition END
  //=============================================================================

  template <typename FitResultType>
  void fillFitResultHistograms( const Track&, const FitResultType&, TypeTrackHistogram& ) const;

  void initializeHistogramMap();

  // GAUDI PROPERTIES
  Gaudi::Property<bool> m_splitByType{this, "SplitByType", true};
  using Type = LHCb::Event::Enum::Track::Type;
  Gaudi::Property<std::vector<Type>> m_typesToMonitor{
      this,
      "typesToMonitor",
      {Type::Velo, Type::VeloBackward, Type::Long, Type::Upstream, Type::Downstream, Type::Ttrack, Type::Muon}};

  // HISTOGRAM DEFINITION
  std::vector<std::unique_ptr<TypeTrackHistogram>>   m_histograms;
  std::array<TypeTrackHistogram*, int( Type::Last )> m_histogrammap = {};

  mutable Gaudi::Accumulators::Histogram<1> m_trackNumbers{this, "nTracks", "# tracks", {50, 0, 500}};
  mutable Gaudi::Accumulators::Histogram<1> m_trackMultiplicity{
      this, "TrackMultiplicityFine", "# tracks", {200, 0.0, 2000.}};
  mutable Gaudi::Accumulators::Histogram<1> m_type{this, "trackType", "track type", {11, -0.5, 10.5}};
  mutable Gaudi::Accumulators::Histogram<1> m_history{this, "history", "track history", {24, -0.5, 23.5}};

  mutable Gaudi::Accumulators::StatCounter<> m_NTracks{this, "#Tracks"};
};

// using TrackM          = TrackMonitor<LHCb::TrackFitResult, LHCb::FitNode>;
// using TrackM_PrKalman = TrackMonitor<LHCb::PrKalmanFitResult, LHCb::Pr::Tracks::Fit::Node>;
DECLARE_COMPONENT_WITH_ID( TrackMonitor, "TrackMonitor" )
DECLARE_COMPONENT_WITH_ID( TrackMonitor, "TrackMonitor_PrKalman" )

//=============================================================================
// Implementation of class functions
//=============================================================================

StatusCode TrackMonitor::initialize() {
  return Consumer::initialize().andThen( [&] { initializeHistogramMap(); } );
}

void TrackMonitor::initializeHistogramMap() {
  // range for residuals for different hittypes
  const std::vector<std::tuple<std::string, double, double>> range{
      {{"8", 0., 100.}, {"9", 0., 50.}, {"10", 0., 100.}, {"11", 0., 50.}, {"12", 0., 100.}}};

  if ( m_splitByType ) {
    // make seperate histogrammer for each requested type
    for ( const auto& t : m_typesToMonitor ) {
      const std::string type{toString( t )};
      const std::string prefix{type + "/"};
      auto              histogrammer = std::make_unique<TypeTrackHistogram>( this, type, prefix );
      for ( auto hittype = 0; hittype < HitType::NHitTypes; ++hittype ) {
        const auto hitname = HitTypeName[hittype];
        histogrammer->hithistograms[hittype] =
            std::make_unique<HitHistograms>( this, prefix, hitname, HitTypeMaxRes[hittype] );
      }
      m_histograms.push_back( std::move( histogrammer ) );
      m_histogrammap[int( t )] = &( *m_histograms.back() );
    }
  } else {
    // make one histogrammer and add all requested track types
    auto histogrammer = std::make_unique<TypeTrackHistogram>( this, "All", "" );
    for ( auto hittype = 0; hittype < HitType::NHitTypes; ++hittype ) {
      const auto hitname = HitTypeName[hittype];
      histogrammer->hithistograms[hittype] =
          std::make_unique<HitHistograms>( this, "", hitname, HitTypeMaxRes[hittype] );
    }
    m_histograms.push_back( std::move( histogrammer ) );
    for ( const auto& t : m_typesToMonitor ) m_histogrammap[int( t )] = &( *m_histograms.back() );
  }
}

void TrackMonitor::operator()( LHCb::Track::Range const& tracks ) const {
  ++m_trackNumbers[tracks.size()];
  ++m_trackMultiplicity[tracks.size()];
  m_NTracks += tracks.size();

  std::array<unsigned int, int( Type::Last )> multiplicityMap = {};
  for ( const LHCb::Track* track : tracks ) {

    ++m_type[static_cast<int>( track->type() )];
    ++m_history[static_cast<int>( track->history() )];

    auto& histos = m_histogrammap[int( track->type() )];
    if ( !histos ) continue;
    multiplicityMap[int( track->type() )] += 1;

    ++histos->probChi2[track->probChi2()];
    ++histos->chi2PerDoF[track->chi2PerDoF()];
    ++histos->ghostProbability[track->ghostProbability()];
    ++histos->nLHCbIDs[track->nLHCbIDs()];
    ++histos->pseudoRapidity[track->pseudoRapidity()];
    ++histos->phi[track->phi()];
    ++histos->nDoF[track->nDoF()];
    ++histos->flag[static_cast<int>( track->flag() )];
    ++histos->typeHistory[static_cast<int>( track->history() )];
    ++histos->fitStatus[static_cast<int>( track->fitStatus() )];

    const LHCb::State& firststate = track->firstState();
    ++histos->firststateX[firststate.x()];
    ++histos->firststateY[firststate.y()];
    ++histos->firststateZ[firststate.z()];
    ++histos->firststateTX[firststate.tx()];
    ++histos->firststateTY[firststate.ty()];
    ++histos->firststateQoverP[firststate.qOverP()];

    if ( firststate.tx() != 0 || firststate.ty() != 0 ) {
      const TrackVector& vec = firststate.stateVector();

      double z = firststate.z();
      z -= ( vec[0] * vec[2] + vec[1] * vec[3] ) / ( vec[2] * vec[2] + vec[3] * vec[3] );
      ++histos->Z[z];
    }

    if ( firststate.qOverP() != 0 ) {
      ++histos->P[track->p()];
      ++histos->PT[track->pt()];
    }

    const std::vector<LHCb::LHCbID>& ids = track->lhcbIDs();
    const auto nUTHits = std::count_if( ids.begin(), ids.end(), []( const LHCb::LHCbID& id ) { return id.isUT(); } );
    const auto nVPHits = std::count_if( ids.begin(), ids.end(), []( const LHCb::LHCbID& id ) { return id.isVP(); } );
    const auto nFTHits = std::count_if( ids.begin(), ids.end(), []( const LHCb::LHCbID& id ) { return id.isFT(); } );
    const auto nMuonHits =
        std::count_if( ids.begin(), ids.end(), []( const LHCb::LHCbID& id ) { return id.isMuon(); } );

    ++histos->nMuonHits[nMuonHits];
    ++histos->nUTHits[nUTHits];
    ++histos->nVPHits[nVPHits];
    ++histos->nFTHits[nFTHits];

    const LHCb::HitPattern hitpattern{track->lhcbIDs()};
    ++histos->nVeloALayers[hitpattern.numVeloA()];
    ++histos->nVeloCLayers[hitpattern.numVeloC()];
    ++histos->nVeloHoles[hitpattern.numVeloHoles()];
    ++histos->nUTLayers[hitpattern.numUT()];
    ++histos->nFTLayers[hitpattern.numFT()];
    ++histos->nFTHoles[hitpattern.numFTHoles()];

    for ( size_t ilay = 0; ilay < hitpattern.veloA().size(); ++ilay )
      if ( hitpattern.veloA().test( ilay ) ) ++histos->HitVeloALayers[ilay];
    for ( size_t ilay = 0; ilay < hitpattern.veloC().size(); ++ilay )
      if ( hitpattern.veloC().test( ilay ) ) ++histos->HitVeloCLayers[ilay];
    for ( size_t ilay = 0; ilay < hitpattern.ut().size(); ++ilay )
      if ( hitpattern.ut().test( ilay ) ) ++histos->HitUTLayers[ilay];
    for ( size_t ilay = 0; ilay < hitpattern.ft().size(); ++ilay )
      if ( hitpattern.ft().test( ilay ) ) ++histos->HitFTLayers[ilay];

    // construct the Velo left-right overlap pattern: this should move to LHCb::HitPattern
    ++histos->nVeloOverlapLayers[hitpattern.numVeloStationsOverlap()];

    // Loop over the range of nodes, the proper function to retrieve it from fit result
    // will be found through ADL. Some functions acting on fit nodes will be called in the same way.
    // For PrKalmanFilter monitoring see PrKalmanFitResult.h and PrFitNode.h
    // For TrackMasterFitter monitoring see TrackFitResult.h and FitNode.h
    if ( track->fitResult() ) {
      auto prFitResult = dynamic_cast<const LHCb::PrKalmanFitResult*>( track->fitResult() );
      if ( prFitResult )
        fillFitResultHistograms( *track, *prFitResult, *histos );
      else {
        auto masterFitResult = dynamic_cast<const LHCb::TrackFitResult*>( track->fitResult() );
        if ( masterFitResult )
          fillFitResultHistograms( *track, *masterFitResult, *histos );
        else {
          ++m_unknownFitType;
        }
      }
    }
  }

  for ( int type = 0; type < int( Type::Last ); ++type ) {
    auto histos = m_histogrammap[type];
    if ( histos ) {
      const auto count = multiplicityMap[type];
      ++histos->Multiplicity[count];
      if ( count > 0 ) // this was probably an unintentional bug, but I'll leave it, to make it easier to test
        histos->MultiplicityCounter += count;
    }
  }
}

template <typename FitResultType>
void TrackMonitor::fillFitResultHistograms( const Track& track, const FitResultType& fitResult,
                                            TrackMonitor::TypeTrackHistogram& histos ) const {
  ++histos.nMeasurements[fitResult.nActiveMeasurements()];
  ++histos.nIter[fitResult.nIter()];
  ++histos.pScatter[fitResult.pScatter()];

  size_t  numoutliers( 0 );
  HitType mtype = VPX; // initialize to avoid compiler warning

  for ( const auto& node : nodes( fitResult ) ) {
    // discard extremely small fraction of hits with zero error
    // on residual. (e.g. a downstream track with only one
    // active TT hit)
    if ( node.isHitOnTrack() && node.errResidual2() > TrackParameters::lowTolerance &&
         ( mtype = hittypemap( node ) ) != NHitTypes ) {
      // factor for unbiasing the rms (not the mean!)
      double f           = std::sqrt( node.errMeasure2() / node.errResidual2() );
      auto&  namedhistos = histos.hithistograms[mtype];
      ++namedhistos->Residual[f * node.residual()];
      ++namedhistos->ResidualPull[node.residual() / node.errResidual()];

      // these should be expert plots because 2D
      // if ( ( mtype == VPX || mtype == VPY ) && fullDetail() ) {
      // 	// calculate R in the local frame
      // 	Gaudi::XYZPoint globalPoint = state( node ).position();
      // 	Gaudi::XYZPoint localPoint  = vpdet.sensor( id( node ).vpID() ).toLocal( globalPoint );
      // 	double          r           = localPoint.Rho();

      // 	// factor to calculate residual in detector plane
      // 	double cosalpha( 1.0 );

      // 	Gaudi::XYZVector localUnitPocaVector = vpdet.sensor( id( node ).vpID() ).toLocal( pocaVector( node ) );
      // 	cosalpha                             = localUnitPocaVector.Rho() / localUnitPocaVector.R();

      // 	++namedhistos->R[{r, node.residual() * f / cosalpha}];
      // }
    } else if ( node.isOutlier() ) {
      ++numoutliers;
    }
  }

  ++histos.NumOutliers[numoutliers];

  const double mom = track.p();
  const double phi = track.phi();

  if ( auto tmp = fitResult.chi2Velo(); tmp.nDoF() > 0 ) ++histos.TmpChi2PerDoF[tmp.chi2PerDoF()];
  if ( auto tmp = fitResult.chi2Downstream(); tmp.nDoF() > 0 ) ++histos.TmpChi2DownStream[tmp.chi2PerDoF()];
  if ( auto tmp = fitResult.chi2Match(); tmp.nDoF() > 0 ) ++histos.TmpChi2Match[tmp.chi2PerDoF()];
  if ( auto tmp = fitResult.chi2Muon(); tmp.nDoF() > 0 ) ++histos.TmpChi2Muon[tmp.chi2PerDoF()];

  if ( fitResult.chi2Velo().nDoF() > 0 ) {
    auto prob = fitResult.chi2Velo().prob();
    histos.Chi2ProbVeloVsMom[mom] += prob;
    histos.Chi2ProbVeloVsPhi[phi] += prob;
  }
  if ( fitResult.chi2Downstream().nDoF() > 0 ) {
    const LHCb::State* Tstate = track.stateAt( LHCb::State::Location::AtT );
    const double       phiT   = Tstate ? std::atan2( Tstate->y(), Tstate->x() ) : phi;
    auto               prob   = fitResult.chi2Downstream().prob();
    histos.Chi2ProbDownstreamVsMom[mom] += prob;
    histos.Chi2ProbDownstreamVsPhi[phiT] += prob;
  }
  if ( fitResult.chi2Match().nDoF() > 0 ) {
    auto prob = fitResult.chi2Match().prob();
    histos.Chi2ProbMatchVsMom[mom] += prob;
    histos.Chi2ProbMatchVsPhi[phi] += prob;
  }
  if ( fitResult.chi2().nDoF() > 0 ) {
    auto prob = fitResult.chi2().prob();
    histos.Chi2ProbVsMom[mom] += prob;
    histos.Chi2ProbVsEta[track.pseudoRapidity()] += prob;
    histos.Chi2ProbVsPhi[phi] += prob;
  }

  // expert checks
  static const double halfOverLog10 = 0.5 / std::log( 10.0 );
  // find first and last node with measurement
  // First locate the first and last node that actually have information
  const typename FitResultType::NodeType *firstMNode( nullptr ), *lastMNode( nullptr );
  for ( const auto& node : nodes( fitResult ) ) {
    if ( node.isHitOnTrack() ) {
      if ( !firstMNode || node.z() < firstMNode->z() ) firstMNode = &node;
      if ( !lastMNode || node.z() > lastMNode->z() ) lastMNode = &node;
    }
  }
  if ( firstMNode ) {
    ++histos.FirstMStateCov0[log( state( *firstMNode ).covariance()( 0, 0 ) ) * halfOverLog10];
    ++histos.FirstMStateCov1[log( state( *firstMNode ).covariance()( 1, 1 ) ) * halfOverLog10];
    ++histos.FirstMStateCov2[log( state( *firstMNode ).covariance()( 2, 2 ) ) * halfOverLog10];
    ++histos.FirstMStateCov3[log( state( *firstMNode ).covariance()( 3, 3 ) ) * halfOverLog10];
    ++histos.FirstMStateCov4[log( state( *firstMNode ).covariance()( 4, 4 ) ) * halfOverLog10];
  }
  if ( lastMNode ) {
    ++histos.LastMStateCov0[log( state( *lastMNode ).covariance()( 0, 0 ) ) * halfOverLog10];
    ++histos.LastMStateCov1[log( state( *lastMNode ).covariance()( 1, 1 ) ) * halfOverLog10];
    ++histos.LastMStateCov2[log( state( *lastMNode ).covariance()( 2, 2 ) ) * halfOverLog10];
    ++histos.LastMStateCov3[log( state( *lastMNode ).covariance()( 3, 3 ) ) * halfOverLog10];
    ++histos.LastMStateCov4[log( state( *lastMNode ).covariance()( 4, 4 ) ) * halfOverLog10];
  }
}
