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
#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "Event/TwoProngVertex.h"
#include "TrackInterfaces/ITrackVertexer.h"
#include "TrackKernel/TrackStateVertex.h"

#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Consumer.h"
#include <Gaudi/Accumulators/Histogram.h>

#include <algorithm>
#include <mutex>

namespace {
  std::vector<const LHCb::State*> firstStates( LHCb::span<const LHCb::Track* const> tracks ) {
    std::vector<const LHCb::State*> states;
    states.reserve( tracks.size() );
    for ( const auto& track : tracks ) { states.push_back( &track->firstState() ); }
    return states;
  }
} // namespace

class TrackVertexMonitor
    : public LHCb::Algorithm::Consumer<void( LHCb::RecVertex::Range const&, LHCb::Track::Range const&,
                                             DetectorElement const& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DetectorElement>> {
public:
  /** Standard construtor */
  TrackVertexMonitor( const std::string& name, ISvcLocator* pSvcLocator );

  /** Algorithm execute */
  void operator()( LHCb::RecVertex::Range const& pvcontainer, LHCb::Track::Range const& alltracks,
                   DetectorElement const& lhcb ) const override;

private:
  Gaudi::Property<double> m_ipmax{this, "MaxIP", 0.5 * Gaudi::Units::mm};
  Gaudi::Property<double> m_ipmaxprof{this, "MaxIPProfile", 0.1 * Gaudi::Units::mm};
  Gaudi::Property<double> m_dzmax{this, "MaxDz", 5 * Gaudi::Units::mm};
  Gaudi::Property<double> m_xpvmax{this, "MaxXPV", 2 * Gaudi::Units::mm};
  Gaudi::Property<double> m_ypvmax{this, "MaxYPV", 1 * Gaudi::Units::mm};
  Gaudi::Property<double> m_zpvmin{this, "MinZPV", -20 * Gaudi::Units::cm};
  Gaudi::Property<double> m_zpvmax{this, "MaxZPV", 20 * Gaudi::Units::cm};
  Gaudi::Property<double> m_zpvmin_wide{this, "MinZPV_Wide", -150 * Gaudi::Units::cm, "Wide z window for PV plot"};
  Gaudi::Property<double> m_zpvmax_wide{this, "MaxZPV_Wide", 150 * Gaudi::Units::cm, "Wide z window for PV plot"};
  Gaudi::Property<double> m_maxLongTrackChisqPerDof{this, "MaxLongTrackChisqPerDof", 5};
  Gaudi::Property<double> m_minLongTrackMomentum{this, "MinLongTrackMomentum", 5};
  Gaudi::Property<unsigned int> m_nprbins{this, "NumProfileBins", 20};
  Gaudi::Property<unsigned int> m_ntracksPV{this, "NumTracksPV", 2};

  ToolHandle<ITrackVertexer> m_vertexer{"TrackVertexer"};

  mutable Gaudi::Accumulators::Histogram<1> m_numTracksPerPV{
      this, "NumTracksPerPV", "NumTracksPerPV", {50, -0.5, 99.5}};
  mutable Gaudi::Accumulators::Histogram<1> m_numLongTracksPerPV{
      this, "NumLongTracksPerPV", "NumLong", {50, -0.5, 99.5}};
  mutable Gaudi::Accumulators::Histogram<1> m_numBackTracksPerPV{
      this, "NumBackTracksPerPV", "NumBackTracksPerPV", {50, -0.5, 99.5}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvChisquarePerDof{
      this, "PV chisquare per dof", "PV chisquare per dof", {150, 0., 3.}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvXPosition{
      this, "PV x position", "PV x position", {200, -m_xpvmax, m_xpvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvYPosition{
      this, "PV y position", "PV y position", {200, -m_ypvmax, m_ypvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvZPosition{
      this, "PV z position", "PV z position", {200, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvZPositionWide{
      this, "PV z position (wide)", "PV z position (wide)", {200, m_zpvmin_wide, m_zpvmax_wide}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLongChisquarePerDof{
      this, "PV long chisquare per dof", "PV long chisquare per dof", {50, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftX{this, "PV left x", "PV left x", {200, -m_xpvmax, m_xpvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftY{this, "PV left y", "PV left y", {200, -m_ypvmax, m_ypvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftZ{this, "PV left z", "PV left z", {200, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvRightX{this, "PV right x", "PV right x", {200, -m_xpvmax, m_xpvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvRightY{this, "PV right y", "PV right y", {200, -m_ypvmax, m_ypvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvRightZ{this, "PV right z", "PV right z", {200, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftRightDeltaX{
      this, "PV left-right delta x", "PV left-right delta x", {50, -0.1, 0.1}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftRightDeltaY{
      this, "PV left-right delta y", "PV left-right delta y", {50, -0.1, 0.1}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftRightDeltaZ{
      this, "PV left-right delta z", "PV left-right delta z", {50, -1, 1}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftRightDeltaXPull{
      this, "PV left-right delta x pull", "PV left-right delta x pull", {50, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftRightDeltaYPull{
      this, "PV left-right delta y pull", "PV left-right delta y pull", {50, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftRightDeltaZPull{
      this, "PV left-right delta z pull", "PV left-right delta z pull", {50, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvLeftChisquarePerDof{
      this, "PV left chisquare per dof", "PV left chisquare per dof", {50, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvRightChisquarePerDof{
      this, "PV right chisquare per dof", "PV right chisquare per dof", {50, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvForwardBackwardDeltaX{
      this, "PV forward-backward delta x", "PV forward-backward delta x", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvForwardBackwardDeltaY{
      this, "PV forward-backward delta y", "PV forward-backward delta y", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvForwardBackwardDeltaZ{
      this, "PV forward-backward delta z", "PV forward-backward delta z", {50, -m_dzmax, m_dzmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvForwardBackwardDeltaXPull{
      this, "PV forward-backward delta x pull", "PV forward-backward delta x pull", {50, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvForwardBackwardDeltaYPull{
      this, "PV forward-backward delta y pull", "PV forward-backward delta y pull", {50, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvForwardBackwardDeltaZPull{
      this, "PV forward-backward delta z pull", "PV forward-backward delta z pull", {50, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvForwardChisquareDof{
      this, "PV forward chisquare per dof", "PV forward chisquare per dof", {50, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1> m_pvBackwardChisquareDof{
      this, "PV backward chisquare per dof", "PV backward chisquare per dof", {50, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1> m_trackIPX{
      this, "track IP X", "track IP X (biased)", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_trackIPY{
      this, "track IP Y", "track IP Y (biased)", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_trackTransverseIP{
      this, "fast track transverse IP", "fast track transverse IP", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_trackLongitudinalIP{
      this, "fast track longitudinal IP", "fast track longitudinal IP", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_fastTrackIPX{
      this, "fast track IP X", "fast track IP X", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_fastTrackIPY{
      this, "fast track IP Y", "fast track IP Y", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngMass{
      this, "twoprong mass (GeV)", "twoprong mass (GeV)", {50, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngMomentum{
      this, "twoprong momentum (GeV)", "twoprong momentum (GeV)", {50, 0, 200}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngDoca{
      this, "twoprong momentum (GeV)", "twoprong momentum (GeV)", {50, 0, 200}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngDocaPull{
      this, "twoprong doca pull", "twoprong doca pull", {50, -m_ipmax, m_ipmax}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngDecayLength{
      this, "twoprong decaylength", "twoprong decaylength", {50, -2, 2}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngDecayLengthSig{
      this, "twoprong decaylength significance", "twoprong decaylength significance", {50, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngCTau{this, "twoprong ctau", "twoprong ctau", {50, -0.1, 0.1}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngProperLifetime{
      this, "twoprong proper lifetime (ps)", "twoprong proper lifetime (ps)", {50, -0.2, 0.2}};
  mutable Gaudi::Accumulators::Histogram<1> m_twoProngIPChi2PerDof{
      this, "twoprong IP chi2 per dof", "twoprong IP chi2 per dof", {50, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1> m_numPrimaryVertices{
      this, "NumPrimaryVertices", "NumPrimaryVertices", {11, -0.5, 10.5}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_pvYvsZ{
      this, "PV y versus z", "PV y versus z", {m_nprbins, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_pvXvsZ{
      this, "PV x versus z", "PV x versus z", {m_nprbins, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_pvLeftRightDeltaYvsZ{
      this, "PV left-right delta y versus z", "PV left-right delta y versus z", {m_nprbins, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_pvLeftRightDeltaXvsZ{
      this, "PV left-right delta x versus z", "PV left-right delta x versus z", {m_nprbins, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_pvForwardBackwardDeltaYvsZ{this,
                                                                                "PV forward-backward delta y versus z",
                                                                                "PV forward-backward delta y versus z",
                                                                                {m_nprbins, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_pvForwardBackwardDeltaXvsZ{this,
                                                                                "PV forward-backward delta x versus z",
                                                                                "PV forward-backward delta x versus z",
                                                                                {m_nprbins, m_zpvmin, m_zpvmax}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_trackIPXvsPhi{
      this, "track IP X vs phi", "track IP X vs phi (biased)", {m_nprbins, -Gaudi::Units::pi, Gaudi::Units::pi}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_trackIPXvsEta{
      this, "track IP X vs eta", "track IP X vs eta (biased)", {m_nprbins, 2.0, 5.0}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_trackIPYvsPhi{
      this, "track IP Y vs phi", "track IP Y vs phi (biased)", {m_nprbins, -Gaudi::Units::pi, Gaudi::Units::pi}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_trackIPYvsEta{
      this, "track IP Y vs eta", "track IP Y vs eta (biased)", {m_nprbins, 2.0, 5.0}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_trackTransverseIPvsPhi{
      this,
      "fast track transverse IP vs phi",
      "fast track transverse IP vs phi",
      {m_nprbins, -Gaudi::Units::pi, Gaudi::Units::pi}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_trackTransverseIPvsEta{
      this, "fast track transverse IP vs eta", "fast track transverse IP vs eta", {m_nprbins, 2.0, 5.0}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_trackLongitudinalIPvsPhi{
      this,
      "fast track longitudinal IP vs phi",
      "fast track longitudinal IP vs phi",
      {m_nprbins, -Gaudi::Units::pi, Gaudi::Units::pi}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_trackLongitudinalIPvsEta{
      this, "fast track longitudinal IP vs eta", "fast track longitudinal IP vs eta", {m_nprbins, 2.0, 5.0}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_fastTrackIPXvsPhi{
      this, "fast track IP X vs phi", "fast track IP X vs phi", {m_nprbins, -Gaudi::Units::pi, Gaudi::Units::pi}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_fastTrackIPXvsEta{
      this, "fast track IP X vs eta", "fast track IP X vs eta", {m_nprbins, 2.0, 5.0}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_fastTrackIPYvsPhi{
      this, "fast track IP Y vs phi", "fast track IP Y vs phi", {m_nprbins, -Gaudi::Units::pi, Gaudi::Units::pi}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_fastTrackIPYvsEta{
      this, "fast track IP Y vs eta", "fast track IP Y vs eta", {m_nprbins, 2.0, 5.0}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_twoProngDocavsEta{
      this, "twoprong doca vs eta", "twoprong doca vs eta", {m_nprbins, 2.0, 5.0}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_twoProngDocavsPhi{
      this, "twoprong doca vs phi", "twoprong doca vs phi", {m_nprbins, -Gaudi::Units::pi, Gaudi::Units::pi}};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TrackVertexMonitor )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TrackVertexMonitor::TrackVertexMonitor( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer{name,
               pSvcLocator,
               {KeyValue{"PVContainer", LHCb::RecVertexLocation::Primary},
                KeyValue{"TrackContainer", LHCb::TrackLocation::Default},
                KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}}} {}

//=============================================================================
// Initialization
//=============================================================================

namespace {
  template <class TrackContainer, class Predicate>
  std::vector<const LHCb::Track*> myselect( const TrackContainer& tracks, Predicate&& selector ) {
    std::vector<const LHCb::Track*> rc;
    std::copy_if( tracks.begin(), tracks.end(), std::back_inserter( rc ), std::forward<Predicate>( selector ) );
    return rc;
  }

  std::vector<const LHCb::Track*> myconvert( const SmartRefVector<LHCb::Track>& tracks ) {
    return myselect( tracks, []( const LHCb::Track* t ) { return t != nullptr; } );
  }

  auto TrackTypePredicate = []( LHCb::Track::Types atype ) {
    return [=]( const LHCb::Track* track ) { return track->type() == atype; };
  };

  auto TrackBackwardPredicate = []() { return [=]( const LHCb::Track* track ) { return track->isVeloBackward(); }; };

  auto TrackForwardPredicate = []() { return [=]( const LHCb::Track* track ) { return !track->isVeloBackward(); }; };

  auto TrackVeloSidePredicate = []( int asign ) {
    // +1: left side only, 0: overlap track, -1: right side only
    // for asign > 0 select left-side tracks only, for a < 0 select right-side tracks, reject overlap tracks
    return [=]( const LHCb::Track* track ) {
      int  side         = 0;
      bool allhitsleft  = true;
      bool allhitsright = true;

      const std::vector<LHCb::LHCbID>& track_ids = track->lhcbIDs();
      for ( const auto& track_id : track_ids ) {
        if ( !track_id.isVP() ) continue;

        auto vp_id = track_id.vpID();
        // 0 should be right, 1 left side
        bool sidepos = vp_id.sidepos();

        allhitsleft  = allhitsleft && sidepos;
        allhitsright = allhitsright && !sidepos;
      }
      if ( allhitsleft ) side = +1;
      if ( allhitsright ) side = -1;
      return side * asign > 0;
    };
  };

} // namespace

void TrackVertexMonitor::operator()( LHCb::RecVertex::Range const& pvcontainer, LHCb::Track::Range const& alltracks,
                                     DetectorElement const& lhcb ) const {
  const auto isLong     = TrackTypePredicate( LHCb::Track::Types::Long );
  const auto isBackward = TrackBackwardPredicate();
  const auto isForward  = TrackForwardPredicate();

  // lists needed
  // - primary vertices
  // - all tracks
  // - long tracks
  // - backward tracks
  // for now I'll just create the track lists from the Best container

  // number of primary vertices
  ++m_numPrimaryVertices[pvcontainer.size()];

  for ( const LHCb::RecVertex* pv : pvcontainer ) {
    auto tracks         = myconvert( pv->tracks() );
    auto forwardtracks  = myselect( tracks, isForward );
    auto backwardtracks = myselect( tracks, isBackward );
    auto longtracks     = myselect( tracks, isLong );

    // number of tracks per primary vertex
    ++m_numTracksPerPV[tracks.size()];
    // number of long tracks per primary vertex
    ++m_numLongTracksPerPV[longtracks.size()];
    // number of backward tracks per primary vertex
    ++m_numBackTracksPerPV[backwardtracks.size()];
    // chisquare
    ++m_pvChisquarePerDof[pv->chi2()];
    // position with crap hack for vertices at exactly 0
    if ( std::abs( pv->position().x() ) > 0.00001 && std::abs( pv->position().y() ) > 0.00001 ) {
      // info() << "pvx " << pv->position().x() << endmsg;
      ++m_pvXPosition[pv->position().x()];
      ++m_pvYPosition[pv->position().y()];
      ++m_pvZPosition[pv->position().z()];
      ++m_pvZPositionWide[pv->position().z()];
    }

    if ( std::abs( pv->position().y() ) < m_ypvmax ) m_pvYvsZ[pv->position().z()] += pv->position().y();
    if ( std::abs( pv->position().x() ) < m_xpvmax ) m_pvXvsZ[pv->position().z()] += pv->position().x();

    // refit the primary vertex with only the long tracks
    if ( longtracks.size() >= 2 ) {
      auto longvertex = m_vertexer->fit( firstStates( longtracks ), *lhcb.geometry() );
      if ( longvertex ) ++m_pvLongChisquarePerDof[longvertex->chi2() / longvertex->nDoF()];
    }

    // now split the primary vertex in left and right tracks
    auto lefttracks  = myselect( tracks, TrackVeloSidePredicate( +1 ) );
    auto righttracks = myselect( tracks, TrackVeloSidePredicate( -1 ) );
    if ( lefttracks.size() >= m_ntracksPV && righttracks.size() >= m_ntracksPV ) {
      // fit two vertices
      auto leftvertex = m_vertexer->fit( firstStates( lefttracks ), *lhcb.geometry() );

      if ( leftvertex ) {
        ++m_pvLeftX[leftvertex->position().x()];
        ++m_pvLeftY[leftvertex->position().y()];
        ++m_pvLeftZ[leftvertex->position().z()];
        /* PK-R3C undefined
        if ( m_leftSensor ) {
          plot( -( m_leftSensor->globalToVeloHalfBox( leftvertex->position() ) ).x(), "PV left-Left half x",
                -m_xpvmax / 4, m_xpvmax / 4 );
          plot( -( m_leftSensor->globalToVeloHalfBox( leftvertex->position() ) ).y(), "PV left-Left half y",
                -m_ypvmax / 2, m_ypvmax / 2 );
        }
        */
      }
      auto rightvertex = m_vertexer->fit( firstStates( righttracks ), *lhcb.geometry() );
      if ( rightvertex ) {
        ++m_pvRightX[rightvertex->position().x()];
        ++m_pvRightY[rightvertex->position().y()];
        ++m_pvRightZ[rightvertex->position().z()];
        /* PK-R3C
        if ( m_rightSensor ) {
          plot( -( m_rightSensor->globalToVeloHalfBox( rightvertex->position() ) ).x(), "PV right-Right half x",
                -m_xpvmax / 4, m_xpvmax / 4 );
          plot( -( m_rightSensor->globalToVeloHalfBox( rightvertex->position() ) ).y(), "PV right-Right half y",
                -m_ypvmax / 2, m_ypvmax / 2 );
        }
        */
      }
      if ( leftvertex && rightvertex ) {
        // draw the difference
        Gaudi::XYZVector dx = leftvertex->position() - rightvertex->position();

        ++m_pvLeftRightDeltaX[dx.x()];
        ++m_pvLeftRightDeltaY[dx.y()];
        ++m_pvLeftRightDeltaZ[dx.z()];
        if ( std::abs( dx.y() ) < m_ipmaxprof ) m_pvLeftRightDeltaYvsZ[pv->position().z()] += dx.y();
        if ( std::abs( dx.x() ) < m_ipmaxprof ) m_pvLeftRightDeltaXvsZ[pv->position().z()] += dx.x();

        // draw the pull of the difference
        Gaudi::SymMatrix3x3 cov = leftvertex->covMatrix() + rightvertex->covMatrix();

        // cov(0,0)
        if ( cov( 0, 0 ) > 1e-10 ) {
          ++m_pvLeftRightDeltaXPull[dx.x() / std::sqrt( cov( 0, 0 ) )];
        } else {
          Info( "cov(0,0) too small", StatusCode::SUCCESS, 10 ).ignore();
        }
        // cov(1,1)
        if ( cov( 1, 1 ) > 1e-10 ) {
          ++m_pvLeftRightDeltaYPull[dx.y() / std::sqrt( cov( 1, 1 ) )];
        } else {
          Info( "cov(1,1) too small", StatusCode::SUCCESS, 10 ).ignore();
        }
        // cov(2,2)
        if ( cov( 2, 2 ) > 1e-10 ) {
          ++m_pvLeftRightDeltaZPull[dx.z() / std::sqrt( cov( 2, 2 ) )];
        } else {
          Info( "cov(2,2) too small", StatusCode::SUCCESS, 10 ).ignore();
        }

        // draw the chisquares
        if ( leftvertex->nDoF() > 0 ) {
          ++m_pvLeftChisquarePerDof[leftvertex->chi2() / leftvertex->nDoF()];
        } else {
          Info( "left ndof = 0", StatusCode::SUCCESS, 10 ).ignore();
        }
        if ( rightvertex->nDoF() > 0 ) {
          ++m_pvRightChisquarePerDof[rightvertex->chi2() / rightvertex->nDoF()];
        } else {
          Info( "right ndof = 0", StatusCode::SUCCESS, 10 ).ignore();
        }
      }
    }

    if ( forwardtracks.size() >= 2 && backwardtracks.size() >= 2 ) {
      // fit two vertices
      auto forwardvertex  = m_vertexer->fit( firstStates( forwardtracks ), *lhcb.geometry() );
      auto backwardvertex = m_vertexer->fit( firstStates( backwardtracks ), *lhcb.geometry() );
      if ( forwardvertex && backwardvertex ) {
        Gaudi::XYZVector dx = forwardvertex->position() - backwardvertex->position();

        // draw the difference
        ++m_pvForwardBackwardDeltaX[dx.x()];
        ++m_pvForwardBackwardDeltaY[dx.y()];
        ++m_pvForwardBackwardDeltaZ[dx.z()];
        if ( std::abs( dx.y() ) < m_ipmaxprof ) m_pvForwardBackwardDeltaYvsZ[pv->position().z()] += dx.y();
        if ( std::abs( dx.x() ) < m_ipmaxprof ) m_pvForwardBackwardDeltaXvsZ[pv->position().z()] += dx.x();

        // draw the pull of the difference
        Gaudi::SymMatrix3x3 cov = forwardvertex->covMatrix() + backwardvertex->covMatrix();
        // cov(0,0)
        if ( cov( 0, 0 ) > 1e-10 ) {
          ++m_pvForwardBackwardDeltaXPull[dx.x() / std::sqrt( cov( 0, 0 ) )];
        } else {
          Info( "cov(0,0) too small", StatusCode::SUCCESS, 10 ).ignore();
        }
        // cov(1,1)
        if ( cov( 1, 1 ) > 1e-10 ) {
          ++m_pvForwardBackwardDeltaYPull[dx.y() / std::sqrt( cov( 1, 1 ) )];
        } else {
          Info( "cov(1,1) too small", StatusCode::SUCCESS, 10 ).ignore();
        }
        // cov(2,2)
        if ( cov( 2, 2 ) > 1e-10 ) {
          ++m_pvForwardBackwardDeltaZPull[dx.z() / std::sqrt( cov( 2, 2 ) )];
        } else {
          Info( "cov(2,2) too small", StatusCode::SUCCESS, 10 ).ignore();
        }
        // draw the chisquares
        if ( forwardvertex->nDoF() > 0 ) {
          ++m_pvForwardChisquareDof[forwardvertex->chi2()];
        } else {
          Info( "forward ndof = 0", StatusCode::SUCCESS, 10 ).ignore();
        }
        if ( backwardvertex->nDoF() > 0 ) {
          ++m_pvBackwardChisquareDof[backwardvertex->chi2()];
        } else {
          Info( "backward ndof = 0", StatusCode::SUCCESS, 10 ).ignore();
        }
      }
    }

    // for events with a single vertex, do something with IP of
    // highest momentum track, as function of phi and eta.
    if ( pvcontainer.size() == 1 && tracks.size() >= 10 ) {

      // now get all good long tracks from the best container:
      auto goodlongtracks = myselect( alltracks, [&]( const LHCb::Track* tr ) {
        return isLong( tr ) && tr->chi2PerDoF() < m_maxLongTrackChisqPerDof && tr->p() > m_minLongTrackMomentum;
      } );

      for ( const LHCb::Track* tr : goodlongtracks ) {
        const LHCb::State& firststate = tr->firstState();
        double             dz         = pv->position().z() - firststate.z();
        double             dx         = firststate.x() + dz * firststate.tx() - pv->position().x();
        double             dy         = firststate.y() + dz * firststate.ty() - pv->position().y();
        Gaudi::XYZVector   p3         = firststate.momentum();
        ++m_trackIPX[dx];
        ++m_trackIPY[dy];
        // apply a cut for the profiles
        if ( std::abs( dx ) < m_ipmaxprof && std::abs( dy ) < m_ipmaxprof ) {
          double phi = p3.phi();
          double eta = p3.eta();
          m_trackIPXvsEta[eta] += dx;
          m_trackIPXvsPhi[phi] += dx;
          m_trackIPYvsEta[eta] += dy;
          m_trackIPYvsPhi[phi] += dy;
        }
      }

      if ( goodlongtracks.size() >= 2 ) {

        std::sort( goodlongtracks.begin(), goodlongtracks.end(), []( const LHCb::Track* lhs, const LHCb::Track* rhs ) {
          return lhs->firstState().pt() < rhs->firstState().pt();
        } );

        const LHCb::Track* firsttrack = goodlongtracks.back();
        goodlongtracks.pop_back();

        // now pick a 2nd track that makes the highest possible invariant mass with this one
        double             highestmass2( 0 );
        const LHCb::Track* secondtrack = nullptr;
        Gaudi::XYZVector   firstp3     = firsttrack->firstState().momentum();
        for ( const auto& t : goodlongtracks ) {
          Gaudi::XYZVector p3    = t->firstState().momentum();
          double           mass2 = p3.r() * firstp3.r() - p3.Dot( firstp3 );
          if ( secondtrack == 0 || highestmass2 < mass2 ) {
            highestmass2 = mass2;
            secondtrack  = t;
          }
        }

        // recompute the vertex without these tracks
        auto newend = tracks.end();
        newend      = std::remove( tracks.begin(), newend, firsttrack );
        newend      = std::remove( tracks.begin(), newend, secondtrack );
        tracks.erase( newend, tracks.end() );
        auto restvertex = m_vertexer->fit( firstStates( tracks ), *lhcb.geometry() );
        if ( restvertex && firsttrack->nStates() != 0 ) {
          const LHCb::State& firststate = firsttrack->firstState();
          double             dz         = restvertex->position().z() - firststate.z();
          double             dx         = firststate.x() + dz * firststate.tx() - restvertex->position().x();
          double             dy         = firststate.y() + dz * firststate.ty() - restvertex->position().y();
          double             nt = std::sqrt( firststate.tx() * firststate.tx() + firststate.ty() * firststate.ty() );
          // transverse and longitudinal impact parameter
          double           iptrans = ( dx * firststate.ty() - dy * firststate.tx() ) / nt;
          double           iplong  = ( dx * firststate.tx() + dy * firststate.ty() ) / nt;
          Gaudi::XYZVector p3      = firststate.momentum();
          double           phi     = p3.phi();
          double           eta     = p3.eta();

          ++m_trackTransverseIP[iptrans];
          ++m_trackLongitudinalIP[iplong];
          ++m_fastTrackIPX[dx];
          ++m_fastTrackIPY[dy];
          // apply a cut for the profiles
          if ( std::abs( iptrans ) < m_ipmaxprof && std::abs( iplong ) < m_ipmaxprof ) {
            m_trackTransverseIPvsEta[eta] += iptrans;
            m_trackTransverseIPvsPhi[phi] += iptrans;
            m_trackLongitudinalIPvsEta[eta] += iplong;
            m_trackLongitudinalIPvsPhi[phi] += iplong;
          }
          if ( std::abs( dx ) < m_ipmaxprof && std::abs( dy ) < m_ipmaxprof ) {
            m_fastTrackIPXvsEta[eta] += dx;
            m_fastTrackIPXvsPhi[phi] += dx;
            m_fastTrackIPYvsEta[eta] += dy;
            m_fastTrackIPYvsPhi[phi] += dy;
          }

          // The two-track cuts we only make for relatively heavy objects
          double mass = std::sqrt( highestmass2 );
          ++m_twoProngMass[mass / Gaudi::Units::GeV];
          if ( mass > 1 * Gaudi::Units::GeV ) {
            // compute doca of two tracks
            Gaudi::XYZVector dx3  = firsttrack->firstState().position() - secondtrack->firstState().position();
            Gaudi::XYZVector n3   = firsttrack->firstState().slopes().Cross( secondtrack->firstState().slopes() );
            double           doca = dx3.Dot( n3 ) / n3.R();
            ++m_twoProngDoca[doca];
            if ( std::abs( doca ) < 200 ) {
              m_twoProngDocavsEta[firstp3.eta()] += doca;
              m_twoProngDocavsPhi[firstp3.phi()] += doca;
            }
            // the easiest way to compute the pull is with a vertex fit
            auto twoprong = m_vertexer->fit( firsttrack->firstState(), secondtrack->firstState(), *lhcb.geometry() );
            if ( twoprong ) {
              double pc = twoprong->p3().R();
              ++m_twoProngMomentum[pc / Gaudi::Units::GeV];
              ++m_twoProngDocaPull[std::sqrt( twoprong->chi2() ) * ( doca > 0 ? 1 : -1 )];
              double chi2, decaylength, decaylengtherr;
              m_vertexer->computeDecayLength( *twoprong, *restvertex, chi2, decaylength, decaylengtherr );
              ++m_twoProngDecayLength[decaylength];
              ++m_twoProngDecayLengthSig[decaylength / decaylengtherr];
              ++m_twoProngIPChi2PerDof[chi2 / 2];
              ++m_twoProngCTau[decaylength * mass / pc];
              ++m_twoProngIPChi2PerDof[decaylength * mass / ( pc * Gaudi::Units::c_light * Gaudi::Units::picosecond )];
            }
          }
        }
      }
    }
  }
}
