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
#include "Event/PrFitNode.h"
#include "Event/PrKalmanFitResult.h"
#include "Event/Track.h"
#include "Event/VPLightCluster.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "Kernel/HitPattern.h"
#include "LHCbAlgs/Consumer.h"
#include "TrackKernel/TrackFunctors.h"
#include "VPDet/DeVP.h"
#include "fmt/format.h"
#include <Gaudi/Accumulators/Histogram.h>

namespace LHCb::Tr::Monitor {

  namespace {

    template <typename TNode>
    inline const Gaudi::TrackProjectionMatrix1D& projectionMatrix( const TNode& node ) {
      // const Gaudi::TrackProjectionMatrix1D& Hk;
      if constexpr ( std::is_same<TNode, LHCb::Pr::Tracks::Fit::Node>::value ) { return node.projectionMatrix(); }
      if constexpr ( std::is_same<TNode, LHCb::FitNode>::value ) {
        return node.template visit_r<const Gaudi::TrackProjectionMatrix1D&>(
            []( const LHCb::FitNode::DimInfos<LHCb::Enum::nDim::Type::one>& d )
                -> const Gaudi::TrackProjectionMatrix1D& { return d.projectionMatrix(); },
            []( ... ) -> const Gaudi::TrackProjectionMatrix1D& {
              throw std::logic_error( "Alignment not implemented for 2D measurements" );
            } );
      }
    }

    template <typename TNode>
    inline LHCb::State& state( const TNode& node ) {
      // const Gaudi::TrackProjectionMatrix1D& Hk;
      if constexpr ( std::is_same<TNode, LHCb::Pr::Tracks::Fit::Node>::value ) {
        return LHCb::Pr::Tracks::Fit::state( node );
      }
      if constexpr ( std::is_same<TNode, LHCb::FitNode>::value ) { return node.state(); }
    }

    constexpr auto nStationsInVelo = VP::NModules / 2; // 2 modules per station
    namespace GA                   = Gaudi::Accumulators;
    const auto modulelimits        = GA::Axis<double>{VP::NModules, -0.5, VP::NModules - 0.5, "module"};
    const auto stationlimits       = GA::Axis<double>{nStationsInVelo, -0.5, nStationsInVelo - 0.5, "station"};
    const auto reslimits           = GA::Axis<double>{100, -0.2, 0.2, "residual [mm]"};

    // Helper class to store an {x,y} pair of nodes
    template <typename TNode>
    struct VPXYNode {
      const TNode*                xnode{0};
      const TNode*                ynode{0};
      LHCb::Detector::VPChannelID vpid;
      double                      residualX{0};
      double                      residualY{0};
    };

  } // namespace

  struct TrackVPOverlapMonitor final : LHCb::Algorithm::Consumer<void( const LHCb::Track::Range& )> {
    /// Standard constructor
    TrackVPOverlapMonitor( const std::string& name, ISvcLocator* pSvcLocator )
        : Consumer{name, pSvcLocator, {KeyValue{"TrackContainer", TrackLocation::Default}}} {}

  private:
    // define histograms which don't need string formatting
    mutable Gaudi::Accumulators::Histogram<2> x_y_stationoverlap_histo{
        this, "x vs y station overlap", "x vs y station overlap", {100, -80.0, 80.0, "X"}, {100, -80.0, 80.0, "Y"}};
    mutable Gaudi::Accumulators::Histogram<2> x_y_tileoverlap_histo{
        this, "x vs y tile overlap", "x vs y tile overlap", {100, -80.0, 80.0, "X"}, {100, -80.0, 80.0, "Y"}};
    mutable Gaudi::Accumulators::Histogram<1> module_histo{this, "module", "module", {52, 0.0, 52.0}};
    mutable Gaudi::Accumulators::Histogram<1> station_histo{this, "station", "station", {27, 0.0, 27.0}};

    mutable Gaudi::Accumulators::Histogram<2> m_overlap_residual_x_CLI_NLO{
        this, "x overlap residual CLI-NLO", "x overlap residual CLI-NLO", modulelimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_overlap_residual_y_CLI_NLO{
        this, "y overlap residual CLI-NLO", "y overlap residual CLI-NLO", modulelimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_overlap_residual_x_CLI_NSI{
        this, "x overlap residual CLI-NSI", "x overlap residual CLI-NSI", modulelimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_overlap_residual_y_CLI_NSI{
        this, "y overlap residual CLI-NSI", "y overlap residual CLI-NSI", modulelimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_overlap_residual_x_CSO_NSI{
        this, "x overlap residual CSO-NSI", "x overlap residual CSO-NSI", modulelimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_overlap_residual_y_CSO_NSI{
        this, "y overlap residual CSO-NSI", "y overlap residual CSO-NSI", modulelimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_ACoverlap_residual_x_vs_station{
        this, "x A-C overlap residual", "x A/C overlap residual", stationlimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_ACoverlap_residual_y_vs_station{
        this, "y A-C overlap residual", "y A/C overlap residual", stationlimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_Aresidual_x_vs_station{
        this, "x residual A side vs station", "x residual A side vs station", stationlimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_Aresidual_y_vs_station{
        this, "y residual A side vs station", "y residual A side vs station", stationlimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_Cresidual_x_vs_station{
        this, "x residual C side vs station", "x residual C side vs station", stationlimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_Cresidual_y_vs_station{
        this, "y residual C side vs station", "y residual C side vs station", stationlimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_residual_x_vs_module{this, "x residual by module",
                                                                     "x residual by module", modulelimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<2> m_residual_y_vs_module{this, "y residual by module",
                                                                     "y residual by module", modulelimits, reslimits};
    mutable Gaudi::Accumulators::Histogram<1> m_breakpoint_delta_x_aside{
        this, "breakpoint delta-x A-side", "breakpoint delta-x A-side", {50, -0.5, 0.5}};
    mutable Gaudi::Accumulators::Histogram<1> m_breakpoint_delta_x_cside{
        this, "breakpoint delta-x C-side", "breakpoint delta-x C-side", {50, -0.5, 0.5}};
    mutable Gaudi::Accumulators::Histogram<1> m_breakpoint_delta_x{
        this, "breakpoint delta-x", "breakpoint delta-x", {50, -0.5, 0.5}};
    mutable Gaudi::Accumulators::Histogram<1> m_breakpoint_delta_y{
        this, "breakpoint delta-y", "breakpoint delta-y", {50, -0.5, 0.5}};
    mutable Gaudi::Accumulators::Histogram<1> m_breakpoint_delta_tx{
        this, "breakpoint delta-tx", "breakpoint delta-tx", {50, -0.01, 0.01}};
    mutable Gaudi::Accumulators::Histogram<1> m_breakpoint_delta_ty{
        this, "breakpoint delta-ty", "breakpoint delta-ty", {50, -0.01, 0.01}};
    mutable Gaudi::Accumulators::ProfileHistogram<1> m_breakpoint_delta_x_vs_z{
        this, "breakpoint delta-x vs z", "breakpoint delta-x vs z", {50, -300., 700.}};
    mutable Gaudi::Accumulators::ProfileHistogram<1> m_breakpoint_delta_y_vs_z{
        this, "breakpoint delta-y vs z", "breakpoint delta-y vs z", {50, -300., 700.}};
    mutable Gaudi::Accumulators::Histogram<1> m_breakpoint_delta_x_pull{
        this, "breakpoint delta-x pull", "breakpoint delta-x pull", {50, -5, 5}};

  public:
    ///< Algorithm execution
    void operator()( const LHCb::Track::Range& tracks ) const override {

      for ( const Track* track : tracks ) { // start of track-loop
        if ( track->hasVelo() && track->checkFitStatus( Track::FitStatus::Fitted ) ) {
          // dispatch based on trackfitresult type
          const auto prfr = dynamic_cast<const LHCb::PrKalmanFitResult*>( track->fitResult() );
          if ( prfr )
            fill( *track, *prfr );
          else {
            const auto fr = dynamic_cast<const LHCb::TrackFitResult*>( track->fitResult() );
            if ( fr ) fill( *track, *fr );
          }
        }
      }
    }

    template <typename TFitResult>
    void fill( const LHCb::Track& track, const TFitResult& fr ) const {
      using TNode = typename TFitResult::NodeType;

      // As it turns out, it doesn't
      // really make sense to look for 'overlap' modules with the
      // same station number: Very often we skip a couple of
      // stations when we switch sides.
      //
      // So, let's do something else: We record the residuals for
      // the 'prediction from tracks that switch sides (like a
      // breakpoint analysis). See also TrackFitMatchMonitor.
      LHCb::HitPattern hitpattern{track.lhcbIDs()};
      const bool       hasacoverlap = hitpattern.numVeloA() >= 2 && hitpattern.numVeloC() >= 2;
      if ( hasacoverlap ) {
        // info() << "Overlap track: " << std::endl << hitpattern << endmsg ;
        const TNode* inode{0};
        for ( const auto& node : nodes( fr ) ) {
          const TNode* jnode = &( node );
          if ( jnode->isVP() ) {
            if ( inode ) {
              auto iside = id( *inode ).vpID().sidepos();
              auto jside = id( *jnode ).vpID().sidepos();
              if ( iside != jside ) { // breakpoint found
                // perhaps it is easier to use the filtered states and monitor the difference in the middle.
                auto       istate = filteredStateForward( *inode );
                auto       jstate = filteredStateBackward( *jnode );
                const auto zmid   = 0.5 * ( istate.z() + jstate.z() );
                istate.linearTransportTo( zmid );
                jstate.linearTransportTo( zmid );
                const Gaudi::TrackVector delta = jstate.stateVector() - istate.stateVector();
                const auto               cov   = jstate.covariance() + istate.covariance();
                if ( iside == 0 ) {
                  ++m_breakpoint_delta_x_aside[delta[0]];
                } else {
                  ++m_breakpoint_delta_x_cside[delta[0]];
                }
                const int sign = iside == 0 ? 1 : -1;
                ++m_breakpoint_delta_x[sign * delta[0]];
                ++m_breakpoint_delta_y[sign * delta[1]];
                ++m_breakpoint_delta_tx[sign * delta[2]];
                ++m_breakpoint_delta_ty[sign * delta[3]];
                m_breakpoint_delta_x_vs_z[zmid] += sign * delta[0];
                m_breakpoint_delta_y_vs_z[zmid] += sign * delta[1];
                ++m_breakpoint_delta_x_pull[sign * delta[0] / std::sqrt( cov( 0, 0 ) )];
              }
            }
            inode = jnode;
          }
        }
      }

      // fill the list of 2D nodes
      using VPXYNodeType = VPXYNode<TNode>;
      std::vector<VPXYNodeType> xynodes;
      xynodes.reserve( nodes( fr ).size() / 2 );
      VPXYNodeType currentnode;
      for ( const auto& node : nodes( fr ) ) {
        if ( node.isVP() ) {
          const auto vpid = id( node ).vpID();
          if ( vpid != currentnode.vpid ) {
            currentnode.vpid  = vpid;
            currentnode.xnode = &node;
          } else {
            currentnode.ynode = &node;
            xynodes.push_back( currentnode );
          }
        }
      }

      if ( xynodes.size() > 1 ) {
        // compute the x and y residual in the global frame.
        // we do not actually need to know which of the original hits was x or y:
        // we exploit the derivative to global x and y (which are in the projection matrix) to rotate the residual
        // to the right frame.
        //
        // res_x = a * node1.residual + b * node2.residual
        //  with a and b such that d res_x / dx = 1 and d res_x / dy = 0, e.g.
        //     1 = d res_x / dx = a * node1.Hx + b * node2.Hx
        //     0 = d res_x / dy = a * node1.Hy + b * node2.Hy
        // -->  b = - a * node1.Hy / node2.Hy
        // --> 1 = a * node1.Hx  - a *  node2.Hx * node1.Hy / node2.Hy
        // --> a = 1 / ( node1.Hx - node2.Hx * node1.Hy / node2.Hy ) = node2.Hy / ( node2.Hy * node1.Hx - node2.Hx *
        // node1.Hy )
        // --> b = - node1.Hy / ( node2.Hy * node1.Hx - node2.Hx * node1.Hy )
        // -->
        //    res_x = 1/( node2.Hy * node1.Hx - node2.Hx * node1.Hy ) * (  node2.Hy * node1.residual - node1.Hy *
        //    node2.residual )
        //
        // In the same fashion,
        //
        //    res_y = 1/( node2.Hx * node1.Hy - node2.Hy * node1.Hx ) * ( node2.Hx * node1.residual - node1.Hx *
        //    node2.residual )
        //          = 1/( node2.Hy * node1.Hx - node2.Hx * node1.Hy ) * ( node1.Hx * node2.residual - node2.Hx *
        //          node1.residual )
        //
        for ( auto& xynode : xynodes ) {
          const auto& node1 = xynode.xnode;
          const auto& node2 = xynode.ynode;
          const auto  res1  = node1->residual();
          const auto  res2  = node2->residual();
          const auto  H1    = projectionMatrix( *node1 )[0];
          const auto  H2    = projectionMatrix( *node2 )[0];
          const auto  mu    = 1. / ( H1[0] * H2[1] - H1[1] * H2[0] );
          xynode.residualX  = mu * ( H2[1] * res1 - H1[1] * res2 );
          xynode.residualY  = mu * ( H1[0] * res2 - H2[0] * res1 );
          ++m_residual_x_vs_module[{xynode.vpid.module(), xynode.residualX}];
          ++m_residual_y_vs_module[{xynode.vpid.module(), xynode.residualY}];
        }

        // Now that we have efficiently computed the residuals in the global frame, we can compute the overlap
        // residuals
        auto jnode = xynodes.begin();
        auto inode = jnode;
        for ( ++jnode; jnode < xynodes.end(); ++jnode, ++inode ) {
          // station overlaps
          if ( jnode->vpid.station() == inode->vpid.station() ) {
            const auto iside = inode->vpid.sidepos();
            const auto jside = jnode->vpid.sidepos();
            if ( iside != jside ) {
              auto nodeA = *inode;
              auto nodeC = *jnode;
              if ( iside != 0 ) std::swap( nodeA, nodeC );
              ++station_histo[nodeA.vpid.station()];
              const auto globalpos = state( *( nodeA.xnode ) ).position();
              ++x_y_stationoverlap_histo[{globalpos.x(), globalpos.y()}];
              ++m_ACoverlap_residual_x_vs_station[{nodeA.vpid.station(), nodeA.residualX - nodeC.residualX}];
              ++m_ACoverlap_residual_y_vs_station[{nodeA.vpid.station(), nodeA.residualY - nodeC.residualY}];
              ++m_Aresidual_x_vs_station[{nodeA.vpid.station(), nodeA.residualX}];
              ++m_Aresidual_y_vs_station[{nodeA.vpid.station(), nodeA.residualY}];
              ++m_Cresidual_x_vs_station[{nodeA.vpid.station(), nodeC.residualX}];
              ++m_Cresidual_y_vs_station[{nodeA.vpid.station(), nodeC.residualY}];
            } else {
              // Collect nodes in different sensors ont he same module. These are always consecutive.
              // Ladder map: { L0: CLI, L1: NLO, L2, NSI, L3: CSO }
              // Expected overlaps: CLI-NLO, CSO-NSI, CLI-NSI; Not allowed: CSO-CLI, NLO-NSI and CSO-NLO
              // const char* tilenames[4] = { "CLI","NLO","NSI","CSO" } ;
              auto node1 = *inode;
              auto node2 = *jnode;
              if ( node2.vpid.sensor() < node1.vpid.sensor() ) std::swap( node1, node2 );
              const int sensor1 = int( node1.vpid.sensor() ) % 4;
              const int sensor2 = int( node2.vpid.sensor() ) % 4;
              const int module  = node1.vpid.module();
              ++module_histo[module];
              const auto globalpos = state( *( node1.xnode ) ).position();
              ++x_y_tileoverlap_histo[{globalpos.x(), globalpos.y()}];

              if ( sensor1 == 0 && sensor2 == 1 ) { // CLI-NLO
                ++m_overlap_residual_x_CLI_NLO[{module, node1.residualX - node2.residualX}];
                ++m_overlap_residual_y_CLI_NLO[{module, node1.residualY - node2.residualY}];
              } else if ( sensor1 == 0 && sensor2 == 2 ) { // CLI-NSI
                ++m_overlap_residual_x_CLI_NSI[{module, node1.residualX - node2.residualX}];
                ++m_overlap_residual_y_CLI_NSI[{module, node1.residualY - node2.residualY}];
              } else if ( sensor1 == 2 && sensor2 == 3 ) { // NSI-CSO
                ++m_overlap_residual_x_CSO_NSI[{module, node1.residualX - node2.residualX}];
                ++m_overlap_residual_y_CSO_NSI[{module, node1.residualY - node2.residualY}];
              } else {
                warning() << "How can these sensors overlap?: " << sensor1 << " " << sensor2 << endmsg;
              }
            }
          }
        }
      }
    }
  };

  DECLARE_COMPONENT_WITH_ID( TrackVPOverlapMonitor, "TrackVPOverlapMonitor" )
  DECLARE_COMPONENT_WITH_ID( TrackVPOverlapMonitor, "TrackVPOverlapMonitor_PrKalman" )
} // namespace LHCb::Tr::Monitor
