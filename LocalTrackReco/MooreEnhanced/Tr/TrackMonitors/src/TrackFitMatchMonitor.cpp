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
#include "Event/FitNode.h"
#include "Event/PrFitNode.h"
#include "Event/PrKalmanFitResult.h"
#include "Event/Track.h"
#include "Event/TrackFitResult.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/ToolHandle.h"
#include "Kernel/HitPattern.h"
#include "LHCbAlgs/Consumer.h"
#include "TrackKernel/TrackFunctors.h"
#include <Gaudi/Accumulators/Histogram.h>
#include <string>
#include <type_traits>

namespace {
  template <typename... Ts, typename I, typename F,
            typename = std::enable_if_t<std::conjunction_v<std::has_virtual_destructor<I>, std::is_base_of<I, Ts>...>>>
  bool dispatch( I* i, F f ) {
    return ( [&]( auto* t ) -> bool {
      if ( t ) f( t );
      return t;
    }( dynamic_cast<Ts*>( i ) ) || ... );
  }
} // namespace

namespace {
  // the end of T1 station of SciFi
  constexpr float ZAfterT1 = 8062 * Gaudi::Units::mm;
  // approximate z position where the kick distribution peaks
  // corresponds roughly to the center of the magnet
  constexpr float ZMagnetKickMean = 5200 * Gaudi::Units::mm;
} // namespace

class TrackFitMatchMonitor : public LHCb::Algorithm::Consumer<void( const LHCb::Track::Range& tracks )> {
public:
  enum ConstrainMethod { All = 0, QOverP = 1, Projective = 2 };

  /** Standard construtor */
  TrackFitMatchMonitor( const std::string& name, ISvcLocator* pSvcLocator );

  /** Algorithm execute */
  void operator()( const LHCb::Track::Range& tracks ) const override;

private:
  // Trackers that have plots produced (may need to change in order to get the relevant Trackers from some centralized
  // variable rather than hard-coding them)
  enum struct Trackers { VeloUT = 0, TUT, VeloT };

  template <typename TNode>
  void plotDelta( Trackers tracker, const TNode& node, bool upstream ) const;
  void makePlots( Trackers tracker, Gaudi::TrackVector deltac, Gaudi::TrackVector deltacpull,
                  Gaudi::TrackVector state ) const;
  void computeConstrainedDelta( Gaudi::TrackVector delta, Gaudi::TrackSymMatrix cov, Gaudi::TrackVector& deltac,
                                Gaudi::TrackVector& deltacpull ) const;
  template <typename TFitResult>
  void fill( const TFitResult& fr, const LHCb::Track& track ) const;

private:
  Gaudi::Property<int>                                  m_constrainMethod{this, "ConstrainMethod", Projective};
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_nullFitResult{this, "Fit result is NULL"};
  mutable Gaudi::Accumulators::ProfileHistogram<1>      m_curvatureRatioTToLongPr{
      this, "curvatureRatioTToLongVsQoP", "curvature ratio T to Long versus q/p", {40, -0.4, 0.4}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_curvatureRatioVeloUTToLongPr{
      this, "curvatureRatioVeloUTToLongVsQoP", "curvature ratio Velo-UT to Long versus q/p", {40, -0.4, 0.4}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_curvatureRatioTToLongVsTxPos{
      this, "curvatureRatioTToLongVsTx", "curvature ratio T to Long versus tx for pos", {40, -0.25, 0.25}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_curvatureRatioVeloUTToLongVsTxPos{
      this, "curvatureRatioVeloUTToLongVsTx", "curvature ratio Velo-UT to Long versus tx for pos", {40, -0.25, 0.25}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_curvatureRatioTToLongVsTxNeg{
      this, "curvatureRatioTToLongVsTx", "curvature ratio T to Long versus tx for neg", {40, -0.25, 0.25}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_curvatureRatioVeloUTToLongVsTxNeg{
      this,
      "curvatureRatioVeloUTToLongVsTxNeg",
      "curvature ratio Velo-UT to Long versus tx for neg",
      {40, -0.25, 0.25}};

  mutable Gaudi::Accumulators::Histogram<1> m_curvatureRatioTToLongH1{
      this, "curvatureRatioTToLong", "curvature ratio T to Long", {40, 0, 2}};
  mutable Gaudi::Accumulators::Histogram<1> m_curvatureRatioVeloUTToLongH1{
      this, "curvatureRatioVeloUTToLong", "curvature ratio Velo-UT to Long", {40, 0, 2}};
  mutable Gaudi::Accumulators::Histogram<1> m_curvatureRatioTToLongPullH1{
      this, "curvatureRatioTToLongPull", "curvature ratio T to Long pull", {40, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_curvatureRatioVeloUTToLongPullH1{
      this, "curvatureRatioVeloUTToLongPull", "curvature ratio Velo-UT to Long pull", {40, -5, 5}};
  mutable Gaudi::Accumulators::Histogram<1> m_kickZH1{this, "kickZ", "Z position of magnet kick", {40, 4900, 5400}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_kickZVsXPr{
      this, "kickZVsXPr", "Z position of magnet kick versus x", {40, -1500, 1500}};

  //=============================================================================
  // Define struct to store delta plot histograms
  //=============================================================================
  struct TrackingDeltaPlots {
    Gaudi::Accumulators::Histogram<1>        dx_dtx_0_hist;
    Gaudi::Accumulators::Histogram<1>        dy_dtx_0_hist;
    Gaudi::Accumulators::Histogram<1>        dtx_dx_0_hist;
    Gaudi::Accumulators::Histogram<1>        dty_dy_0_hist;
    Gaudi::Accumulators::Histogram<1>        dx_pull_hist;
    Gaudi::Accumulators::Histogram<1>        dy_pull_hist;
    Gaudi::Accumulators::Histogram<1>        dtx_pull_hist;
    Gaudi::Accumulators::Histogram<1>        dty_pull_hist;
    Gaudi::Accumulators::ProfileHistogram<1> dx_tx_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dx_ty_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dy_tx_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dy_ty_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dx_pull_tx_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dx_pull_ty_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dx_pull_qop_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dx_qop_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dy_pull_tx_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dy_pull_ty_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dtx_qop_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dtx_pull_qop_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dtx_pull_tx_profile;
    Gaudi::Accumulators::ProfileHistogram<1> dty_pull_ty_profile;
    TrackingDeltaPlots( const TrackFitMatchMonitor* owner, std::string const& Tracker )
        : dx_dtx_0_hist{owner, Tracker + "/dx for dtx==0", Tracker + " dx for dtx==0", {100, -20, 20}}
        , dy_dtx_0_hist{owner, Tracker + "/dy for dty==0", Tracker + " dy for dty==0", {100, -20, 20}}
        , dtx_dx_0_hist{owner, Tracker + "/dtx for dx==0", Tracker + " dtx for dx==0", {100, -0.010, 0.010}}
        , dty_dy_0_hist{owner, Tracker + "/dty for dy==0", Tracker + " dty for dy==0", {100, -0.010, 0.010}}
        , dx_pull_hist{owner, Tracker + "/dx pull", Tracker + " dx pull", {100, -10, 10}}
        , dy_pull_hist{owner, Tracker + "/dy pull", Tracker + " dy pull", {100, -10, 10}}
        , dtx_pull_hist{owner, Tracker + "/dtx pull", Tracker + " dtx pull", {100, -10, 10}}
        , dty_pull_hist{owner, Tracker + "/dty pull", Tracker + " dty pull", {100, -10, 10}}
        , dx_tx_profile{owner, Tracker + "/dx vs tx", Tracker + " dx vs tx", {100, -0.25, 0.25}}
        , dx_ty_profile{owner, Tracker + "/dx vs ty", Tracker + " dx vs ty", {100, -0.25, 0.25}}
        , dy_tx_profile{owner, Tracker + "/dy vs tx", Tracker + " dy vs tx", {100, -0.25, 0.25}}
        , dy_ty_profile{owner, Tracker + "/dy vs ty", Tracker + " dy vs ty", {100, -0.25, 0.25}}
        , dx_pull_tx_profile{owner, Tracker + "/dx pull vs tx", Tracker + " dx pull vs tx", {100, -0.25, 0.25}}
        , dx_pull_ty_profile{owner, Tracker + "/dx pull vs ty", Tracker + " dx pull vs ty", {100, -0.25, 0.25}}
        , dx_pull_qop_profile{owner, Tracker + "/dx pull vs qop", Tracker + " dx pull vs qop", {40, -0.2, 0.2}}
        , dx_qop_profile{owner, Tracker + "/dx vs qop", Tracker + " dx vs qop", {40, -0.2, 0.2}}
        , dy_pull_tx_profile{owner, Tracker + "/dy pull vs tx", Tracker + " dy pull vs tx", {20, -0.25, 0.25}}
        , dy_pull_ty_profile{owner, Tracker + "/dy pull vs ty", Tracker + " dy pull vs ty", {20, -0.25, 0.25}}
        , dtx_qop_profile{owner, Tracker + "/dtx vs qop", Tracker + " dtx vs qop", {40, -0.2, 0.2}}
        , dtx_pull_qop_profile{owner, Tracker + "/dtx pull vs qop", Tracker + " dtx pull vs qop", {40, -0.2, 0.2}}
        , dtx_pull_tx_profile{owner, Tracker + "/dtx pull vs tx", Tracker + " dtx pull vs tx", {20, -0.25, 0.25}}
        , dty_pull_ty_profile{owner, Tracker + "/dty pull vs ty", Tracker + " dty pull vs ty", {20, -0.25, 0.25}} {}
  };
  // map to associate tracking sub-detector to the relevant struct of histograms
  std::array<std::unique_ptr<TrackingDeltaPlots>, 3> m_histograms = [&] {
    auto create = [&]( Trackers t ) { return std::make_unique<TrackingDeltaPlots>( this, toString( t ) ); };
    return std::array{create( Trackers::VeloUT ), create( Trackers::TUT ), create( Trackers::VeloT )};
  }();

  /** Friend function for converting enum values to strings **/
  std::string toString( Trackers t ) const {
    switch ( t ) {
    case Trackers::VeloUT:
      return "Velo-UT";
    case Trackers::TUT:
      return "T-UT";
    case Trackers::VeloT:
      return "Velo-T";
    }
    throw GaudiException( "Unknown Tracker enum value", __func__, StatusCode::FAILURE );
  }
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT_WITH_ID( TrackFitMatchMonitor, "TrackFitMatchMonitor" )
DECLARE_COMPONENT_WITH_ID( TrackFitMatchMonitor, "TrackFitMatchMonitor_PrKalman" )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TrackFitMatchMonitor::TrackFitMatchMonitor( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer{name, pSvcLocator, {"TrackContainer", LHCb::TrackLocation::Default}} {}

//=============================================================================

void TrackFitMatchMonitor::makePlots( Trackers tracker, Gaudi::TrackVector deltac, Gaudi::TrackVector deltacpull,
                                      Gaudi::TrackVector state ) const {
  // get the pointer to the relevant histograms to fill
  auto& histos = m_histograms.at( int( tracker ) );
  // these titles are only right if you choose 'Projective'
  ++histos->dx_dtx_0_hist[deltac( 0 )];
  ++histos->dy_dtx_0_hist[deltac( 1 )];
  ++histos->dtx_dx_0_hist[deltac( 2 )];
  ++histos->dty_dy_0_hist[deltac( 3 )];
  ++histos->dx_pull_hist[deltacpull( 0 )];
  ++histos->dy_pull_hist[deltacpull( 1 )];
  ++histos->dtx_pull_hist[deltacpull( 2 )];
  ++histos->dty_pull_hist[deltacpull( 3 )];

  if ( std::abs( deltacpull( 0 ) ) < 5 ) {
    histos->dx_pull_tx_profile[state[2]] += deltacpull( 0 );
    histos->dx_tx_profile[state[2]] += deltac( 0 );
    histos->dx_ty_profile[state[3]] += deltac( 0 );
    histos->dy_tx_profile[state[2]] += deltac( 1 );
    histos->dy_ty_profile[state[3]] += deltac( 1 );
    histos->dx_pull_ty_profile[state[3]] += deltacpull( 0 );
    histos->dx_pull_qop_profile[state[4] * Gaudi::Units::GeV] += deltacpull( 0 );
    histos->dx_qop_profile[state[4] * Gaudi::Units::GeV] += deltac( 0 );
  }
  if ( std::abs( deltacpull( 1 ) ) < 5 ) {
    histos->dy_pull_tx_profile[state[2]] += deltacpull( 1 );
    histos->dy_pull_ty_profile[state[3]] += deltacpull( 1 );
  }
  if ( std::abs( deltacpull( 2 ) ) < 5 ) {
    histos->dtx_qop_profile[state[4] * Gaudi::Units::GeV] += deltac( 2 );
    histos->dtx_pull_qop_profile[state[4] * Gaudi::Units::GeV] += deltacpull( 2 );
  }
  if ( std::abs( deltacpull( 2 ) ) < 5 ) {
    histos->dtx_pull_tx_profile[state[2]] += deltacpull( 2 );
    histos->dty_pull_ty_profile[state[3]] += deltacpull( 2 );
  }
}

template <typename TNode>
void TrackFitMatchMonitor::plotDelta( Trackers tracker, const TNode& node, bool upstream ) const {
  // It's the only place where the predicted state is requested.
  const LHCb::State& stateUp   = upstream ? filteredStateForward( node ) : predictedStateForward( node );
  const LHCb::State& stateDown = upstream ? predictedStateBackward( node ) : filteredStateBackward( node );

  // compute the difference
  Gaudi::TrackVector    delta = stateUp.stateVector() - stateDown.stateVector();
  Gaudi::TrackSymMatrix cov   = stateUp.covariance() + stateDown.covariance();

  // now, if we just look at the difference, then the problem is that
  // there are very large correlations: e.g. when you step through the
  // magnet, you don't know the momentum yet. so, if the momentum is
  // off, then so are x and tx. To solve this problem we compute a
  // 'constrained' difference: compute the difference constraining the
  // difference in other parameters to zero. There are two modes of operation:
  // a) constrain all 'other' variables (so for 'dx' constraint 'dy=dty=dtx=dqop=0')
  // b) constrain only qop
  // The results are different, mainly because there is something
  // wrong in the fitted tracks already: We find for MC tracks that
  // things don't match very well.

  Gaudi::TrackVector deltac, deltacpull;
  computeConstrainedDelta( delta, cov, deltac, deltacpull );

  makePlots( tracker, deltac, deltacpull, state( node ).stateVector() );
}

void TrackFitMatchMonitor::computeConstrainedDelta( Gaudi::TrackVector delta, Gaudi::TrackSymMatrix cov,
                                                    Gaudi::TrackVector& deltac, Gaudi::TrackVector& deltacpull ) const {
  if ( m_constrainMethod == All ) {
    // now remove the contribution from the difference in the 'other' variables
    for ( size_t irow = 0; irow < 5; ++irow ) {
      // remove this row from delta and cov
      Gaudi::Vector4      subdelta, subcor;
      Gaudi::SymMatrix4x4 subcov;
      size_t              krow( 0 );
      for ( size_t jrow = 0; jrow < 5; ++jrow )
        if ( jrow != irow ) {
          subdelta( krow ) = delta( jrow );
          subcor( krow )   = cov( irow, jrow );
          size_t kcol( 0 );
          for ( size_t jcol = 0; jcol <= jrow; ++jcol )
            if ( jcol != irow ) {
              subcov( krow, kcol ) = cov( jrow, jcol );
              ++kcol;
            }
          ++krow;
        }

      // now invert the subcov
      subcov.Invert();
      // compute delta and its covariance
      Gaudi::Vector4 tmp = subcov * subcor;
      deltac( irow )     = delta( irow ) - ROOT::Math::Dot( subdelta, tmp );
      double covc        = cov( irow, irow ) - ROOT::Math::Dot( subcor, tmp );
      if ( covc > 0 )
        deltacpull( irow ) = deltac( irow ) / std::sqrt( covc );
      else
        warning() << "problem with covc: " << irow << " " << covc << " " << cov( irow, irow ) << " "
                  << ROOT::Math::Dot( subcor, tmp ) << endmsg;
    }
  } else {
    int map[4] = {2, 3, 0, 1};
    for ( size_t irow = 0; irow < 4; ++irow ) {
      int ref            = m_constrainMethod == QOverP ? 4 : map[irow];
      deltac( irow )     = delta( irow ) - cov( irow, ref ) / cov( ref, ref ) * delta( ref );
      double covc        = cov( irow, irow ) - cov( irow, ref ) / cov( ref, ref ) * cov( ref, irow );
      deltacpull( irow ) = deltac( irow ) / std::sqrt( covc );
    }
  }
}

void TrackFitMatchMonitor::operator()( const LHCb::Track::Range& tracks ) const {
  for ( const LHCb::Track* track : tracks ) {
    if ( !track->fitResult() ) {
      ++m_nullFitResult;
      continue;
    }
    // dispatch based on trackfitresult type
    dispatch<const LHCb::PrKalmanFitResult, const LHCb::TrackFitResult>(
        track->fitResult(), [&]( const auto* fr ) { fill( *fr, *track ); } );
  } // loop over tracks
}

template <typename TFitResult>
void TrackFitMatchMonitor::fill( const TFitResult& fitResult, const LHCb::Track& track ) const {
  if ( nodes( fitResult ).size() <= 0 ) return;

  const typename TFitResult::NodeType *lastVelo( nullptr ), *firstUT( nullptr ), *lastUT( nullptr ), *firstT( nullptr );
  // The appropriate function node( TFitResult* ) will be found using ADL.
  // For TrackFitResult or KalmanFitResult it returns Range of LHCb::FitNodes, see TrackFitResult.h and FitNode.h
  // For PrKalmanFitResult it returns std::span of PrFitNodes, see PrKalmanFitResult.h and PrFitNode.h
  // The Range object for LHCb::FitNode is written in such a way that we can make loops over nodes look the same for
  // both types. This makes it easier to template the algorithm. Note that the methods of both fit nodes with same
  // names work internally in a different way as the two objects are different.
  for ( const auto& node : nodes( fitResult ) ) {
    if ( node.isVP() ) {
      if ( !lastVelo || lastVelo->z() < node.z() ) lastVelo = &node;
    } else if ( node.isUT() ) {
      if ( !firstUT || firstUT->z() > node.z() ) firstUT = &node;
      if ( !lastUT || lastUT->z() < node.z() ) lastUT = &node;
    } else if ( node.isFT() ) {
      if ( !firstT || firstT->z() > node.z() ) firstT = &node;
    }
  }
  if ( lastVelo ) {
    if ( lastUT ) {
      plotDelta( Trackers::VeloUT, *firstUT, true );
      if ( firstT ) plotDelta( Trackers::TUT, *firstT, true );
    } else if ( firstT ) {
      plotDelta( Trackers::VeloT, *firstT, true );
    }
  }

  // inspired by the problems we see in the field. see also UT field study
  LHCb::HitPattern hitpattern{track.lhcbIDs()};
  const bool       hasT    = hitpattern.numFT() > 0;
  const bool       hasVelo = hitpattern.numVelo() > 0;
  const bool       hasUT   = hitpattern.numUT() > 0;

  if ( hasT && hasVelo && std::abs( track.firstState().qOverP() ) > 0 ) {
    // first make sure that we have hits in all 3 T stations
    auto hitsInStation = std::array{0, 0, 0};
    for ( const auto& node : nodes( fitResult ) ) {
      if ( node.isHitOnTrack() ) {
        const LHCb::LHCbID lhcbid = id( node );
        if ( lhcbid.isFT() ) hitsInStation[to_unsigned( lhcbid.ftID().station() ) - 1] += 1;
      }
    }

    if ( hitsInStation[0] >= 3 && hitsInStation[1] >= 3 && hitsInStation[2] >= 3 ) {
      // first get the 3 measurements of the curvature with error
      // nodes are sorted in decreasing z. find the nodes around the magnet
      // note that when using LHCb::FitNode the node before might be EndRich1
      // while that is not possible with PrKalman - here the node before UTHit
      // is VPHit
      const typename TFitResult::NodeType *nodeAfter( nullptr ), *nodeBefore( nullptr ), *firstNodeAfterT1( nullptr );
      for ( const auto& node : nodes( fitResult ) ) {
        // reject reference nodes and outliers
        // for LHCb::FitNode this affects the distribution of zkick!
        if ( !node.isHitOnTrack() ) continue;
        if ( node.z() > ZMagnetKickMean ) {
          if ( !nodeAfter || nodeAfter->z() > node.z() ) nodeAfter = &node;
        } else {
          if ( !nodeBefore || nodeBefore->z() < node.z() ) nodeBefore = &node;
        }
        if ( node.z() > ZAfterT1 )
          if ( !firstNodeAfterT1 || firstNodeAfterT1->z() > node.z() ) firstNodeAfterT1 = &node;
      }

      if ( nodeBefore && nodeAfter && firstNodeAfterT1 ) {

        const auto qop = state( *nodeBefore ).qOverP();
        const auto tx  = state( *nodeBefore ).tx();

        // extract the 'upstream' filtered state of T segment
        const auto& first    = nodes( fitResult ).front();
        const auto& last     = nodes( fitResult ).back();
        const bool  upstream = first.z() > last.z();

        const LHCb::State stateT = upstream ? filteredStateForward( *nodeAfter ) : filteredStateBackward( *nodeAfter );

        const auto qopT    = stateT.qOverP();
        const auto qoperrT = std::sqrt( stateT.covariance()( 4, 4 ) );
        ++m_curvatureRatioTToLongH1[qopT / qop];
        ++m_curvatureRatioTToLongPullH1[( qopT - qop ) / qoperrT];
        if ( std::abs( qopT / qop - 1 ) < 1 ) {
          m_curvatureRatioTToLongPr[qop * Gaudi::Units::GeV] += qopT / qop;
          if ( qop > 0 )
            m_curvatureRatioTToLongVsTxPos[tx] += qopT / qop;
          else
            m_curvatureRatioTToLongVsTxNeg[tx] += qopT / qop;
          m_curvatureRatioTToLongPr[qop * Gaudi::Units::GeV] += qopT / qop;
        }

        // extract the 'downstream' filtered state of Velo-UT segment
        const LHCb::State stateVeloUT =
            upstream ? filteredStateBackward( *nodeBefore ) : filteredStateForward( *nodeBefore );

        if ( hasUT ) {
          const auto qopVeloUT    = stateVeloUT.qOverP();
          const auto qoperrVeloUT = std::sqrt( stateVeloUT.covariance()( 4, 4 ) );

          ++m_curvatureRatioVeloUTToLongH1[qopVeloUT / qop];
          ++m_curvatureRatioVeloUTToLongPullH1[( qopVeloUT - qop ) / qoperrVeloUT];

          if ( std::abs( qopVeloUT / qop - 1 ) < 1 ) {
            m_curvatureRatioVeloUTToLongPr[qop * Gaudi::Units::GeV] += qopVeloUT / qop;
            if ( qop > 0 )
              m_curvatureRatioVeloUTToLongVsTxPos[tx] += qopVeloUT / qop;
            else
              m_curvatureRatioVeloUTToLongVsTxNeg[tx] += qopVeloUT / qop;
          }
        }

        // compute the (x,z) point of the intersection of the 2 segments for linear propagation
        // FIXME: it must be better to take a fixed z position in T.
        if ( 1 / std::abs( qop ) > 5 * Gaudi::Units::GeV ) {
          const auto xT  = stateT.x();
          const auto txT = stateT.tx();
          const auto zT  = stateT.z();

          const auto xVeloUT  = stateVeloUT.x();
          const auto txVeloUT = stateVeloUT.tx();
          const auto zVeloUT  = stateVeloUT.z();

          const double zkick = ( zVeloUT * txVeloUT - xVeloUT + xT - zT * txT ) / ( txVeloUT - txT );
          const double xkick = xT + ( zkick - zT ) * txT;
          // double xkickprime = stateVeloUT.x() + (zkick - stateVeloUT.z()) * stateVeloUT.tx() ;
          ++m_kickZH1[zkick];
          // when rejecting reference nodes zkick distribution peaks around 5200mm
          // and has a longer tail on left side,
          // for profile plot remove outliers
          if ( ( ZMagnetKickMean - 200 ) < zkick && zkick < ( ZMagnetKickMean + 200 ) ) m_kickZVsXPr[xkick] += zkick;
        }
      }
    }
  }
}
