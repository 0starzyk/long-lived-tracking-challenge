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
#include "AIDA/IHistogram1D.h"
#include "AIDA/IProfile1D.h"
#include "CaloDet/DeCalorimeter.h"
#include "DetDesc/DetectorElement.h"
#include "Event/CaloCluster.h"
#include "Event/CaloPosition.h"
#include "Event/FitNode.h"
#include "Event/State.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "LHCbAlgs/Consumer.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackKernel/TrackFunctors.h"

#include "GaudiAlg/GaudiHistoAlg.h"

#include <optional>

namespace {
  using GeometricZ = double;
} // namespace

namespace LHCb {

  class TrackCaloMatchMonitor
      : public Algorithm::Consumer<void( CaloClusters const&, Track::Range const&, DetectorElement const&,
                                         GeometricZ const& ),
                                   DetDesc::usesBaseAndConditions<GaudiHistoAlg, DetectorElement, GeometricZ>> {
  public:
    TrackCaloMatchMonitor( const std::string& name, ISvcLocator* pSvcLocator )
        : Consumer( name, pSvcLocator,
                    {{"CaloClustersLocation", CaloClusterLocation::Ecal},
                     {"TrackLocation", TrackLocation::Default},
                     {"StandardGeometryTop", standard_geometry_top},
                     {"GeometricZLocation", name + "GeometricZLocation"}} ) {}
    StatusCode initialize() override;
    void       operator()( CaloClusters const&, Track::Range const&, DetectorElement const&,
                     GeometricZ const& ) const override;

  private:
    Gaudi::Property<bool>          m_requireTHits{this, "RequireTHits", true}; // this you need for cosmics (MIPs)
    Gaudi::Property<double>        m_clusterZCorrection{this, "ClusterZCorrection", 0};
    ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackMasterExtrapolator"};

    IHistogram1D* m_dxASideH1[3]{nullptr, nullptr, nullptr};
    IHistogram1D* m_dyASideH1[3]{nullptr, nullptr, nullptr};
    IHistogram1D* m_dxCSideH1[3]{nullptr, nullptr, nullptr};
    IHistogram1D* m_dyCSideH1[3]{nullptr, nullptr, nullptr};
    IHistogram1D* m_dyVeloASideH1[3]{nullptr, nullptr, nullptr};
    IHistogram1D* m_dyVeloCSideH1[3]{nullptr, nullptr, nullptr};
    IHistogram1D* m_zH1[3]{nullptr, nullptr, nullptr};
    IHistogram1D* m_eOverPH1[3]{nullptr, nullptr, nullptr};
    IProfile1D*   m_dyVsYPr;
    IProfile1D*   m_dyVsXPr;
    IProfile1D*   m_dyVsTyPr;
    IProfile1D*   m_dyVsTxPr;
    IProfile1D*   m_dxVsYPr;
    IProfile1D*   m_dxVsXPr;
    IProfile1D*   m_dxVsTyPr;
    IProfile1D*   m_dxVsTxPr;
  };

  // Declaration of the Algorithm Factory
  DECLARE_COMPONENT_WITH_ID( TrackCaloMatchMonitor, "TrackEcalMatchMonitor" )

} // namespace LHCb

StatusCode LHCb::TrackCaloMatchMonitor::initialize() {
  return GaudiHistoAlg::initialize().andThen( [&] {
    addConditionDerivation( {DeCalorimeterLocation::Ecal}, inputLocation<GeometricZ>(),
                            [this]( DeCalorimeter const& caloDet ) {
                              // This should not be called twice, as we would mess up the histos
                              // In other terms, this code can only run with all events sharing the
                              // same condition for caloDet. Let's fail if we are called more than once
                              if ( m_zH1[0] ) {
                                throw GaudiException( "This code does not support multiple versions of a condition",
                                                      "TrackCaloMatchMonitor::initialize", StatusCode::FAILURE );
                              }
                              auto geoZ             = caloDet.toGlobal( Gaudi::XYZPoint() ).z() + caloDet.zOffset();
                              char systitle[3][128] = {"outer", "middle", "inner"};
                              char histitle[512];
                              for ( int i = 0; i < 3; ++i ) {
                                sprintf( histitle, "zMatch (%s)", systitle[i] );
                                m_zH1[i] = book1D( histitle, geoZ - 400, geoZ + 400 );
                              }
                              return geoZ;
                            } );
    setHistoTopDir( "Track/" );
    char systitle[3][128] = {"outer", "middle", "inner"};
    char histitle[512];
    for ( int i = 0; i < 3; ++i ) {
      sprintf( histitle, "xEcal - xTRK (%s A-side)", systitle[i] );
      m_dxASideH1[i] = book1D( histitle, -200, 200 );
      sprintf( histitle, "yxEcal - yTRK (%s A-side)", systitle[i] );
      m_dyASideH1[i] = book1D( histitle, -200, 200 );
      sprintf( histitle, "xxEcal - xTRK (%s C-side)", systitle[i] );
      m_dxCSideH1[i] = book1D( histitle, -200, 200 );
      sprintf( histitle, "yxEcal - yTRK (%s C-side)", systitle[i] );
      m_dyCSideH1[i] = book1D( histitle, -200, 200 );
      sprintf( histitle, "yxEcal - yVELO (%s A-side)", systitle[i] );
      m_dyVeloASideH1[i] = book1D( histitle, -200, 200 );
      sprintf( histitle, "yxEcal - yVELO (%s C-side)", systitle[i] );
      m_dyVeloCSideH1[i] = book1D( histitle, -200, 200 );
      sprintf( histitle, "E over P (%s)", systitle[i] );
      m_eOverPH1[i] = book1D( histitle, -2, 2 );
    }

    m_dxVsXPr  = bookProfile1D( "dxVsX", "dx versus x", -3500, 3500 );
    m_dxVsYPr  = bookProfile1D( "dxVsY", "dx versus y", -3500, 3500 );
    m_dxVsTxPr = bookProfile1D( "dxVsTx", "dx versus tx", -0.6, 0.6 );
    m_dxVsTyPr = bookProfile1D( "dxVsTy", "dx versus ty", -0.3, 0.3 );

    m_dyVsXPr  = bookProfile1D( "dyVsX", "dy versus x", -3500, 3500 );
    m_dyVsYPr  = bookProfile1D( "dyVsY", "dy versus y", -3500, 3500 );
    m_dyVsTxPr = bookProfile1D( "dyVsTx", "dy versus tx", -0.6, 0.6 );
    m_dyVsTyPr = bookProfile1D( "dyVsTy", "dy versus ty", -0.3, 0.3 );

    return StatusCode::SUCCESS;
  } );
}

namespace {

  struct MyCaloPosition {
    MyCaloPosition( const LHCb::Detector::Calo::CellID& _cell, const LHCb::CaloPosition& _pos )
        : cell( _cell ), pos( _pos ) {}
    LHCb::Detector::Calo::CellID cell;
    LHCb::CaloPosition           pos;
  };

  const LHCb::State* unbiasedVeloState( const LHCb::Track& track ) {
    auto nodes_ = nodes( track );
    if ( !track.hasVelo() || nodes_.empty() ) return nullptr;
    // find the most-downstream node still within velo
    const LHCb::FitNode* node =
        std::accumulate( nodes_.begin(), nodes_.end(), static_cast<const LHCb::FitNode*>( nullptr ),
                         []( const LHCb::FitNode* n, const LHCb::FitNode* i ) {
                           if ( i->z() >= StateParameters::ZEndVelo ) return n;
                           return ( !n || ( i->z() > n->z() ) ) ? i : n;
                         } );

    const bool upstream =
        ( ( *nodes_.begin() )->z() > ( *nodes_.rbegin() )->z() ); // FIXME: when we use range::span, we can use 'front'
                                                                  // and 'back' again...
    if ( node ) {
      // TrackMasterFit
      return upstream ? &node->filteredStateBackward() : &node->filteredStateForward();
    }
    return nullptr;
  }
} // namespace

void LHCb::TrackCaloMatchMonitor::operator()( CaloClusters const& clusters, Track::Range const& trackcontainer,
                                              DetectorElement const& lhcb, GeometricZ const& geometricZ ) const {
  auto& geometry = *lhcb.geometry();
  // make a list of calo positions, depending on the subsystem use clusters or cells
  std::vector<MyCaloPosition> calopositions;

  calopositions.reserve( clusters.size() );
  for ( const CaloCluster* cluster : clusters ) { calopositions.emplace_back( cluster->seed(), cluster->position() ); }

  for ( MyCaloPosition& cluster : calopositions ) cluster.pos.setZ( geometricZ );

  for ( const Track* track : trackcontainer ) {
    if ( !m_requireTHits || track->hasT() ) {
      const State& closest = closestState( *track, geometricZ );
      StateVector  state   = {closest.stateVector(), closest.z()};
      m_extrapolator->propagate( state, geometricZ, geometry ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      const State* fullvelostate = unbiasedVeloState( *track );

      std::optional<StateVector> velostate;
      if ( fullvelostate ) {
        velostate = StateVector( fullvelostate->stateVector(), fullvelostate->z() );
        velostate->setQOverP( state.qOverP() );
        m_extrapolator->propagate( *velostate, geometricZ, geometry )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      }

      for ( const MyCaloPosition& cluster : calopositions ) {
        // state = &(closestState(*track,pos.z())) ;
        double dz     = cluster.pos.z() + m_clusterZCorrection.value() - state.z();
        double xtrack = state.x() + state.tx() * dz;
        double ytrack = state.y() + state.ty() * dz;
        double dx     = cluster.pos.x() - xtrack;
        double dy     = cluster.pos.y() - ytrack;
        if ( std::abs( dy ) < 200 && std::abs( dx ) < 200 ) {
          if ( xtrack > 0 ) {
            m_dxASideH1[cluster.cell.area()]->fill( dx );
            m_dyASideH1[cluster.cell.area()]->fill( dy );
          } else {
            m_dxCSideH1[cluster.cell.area()]->fill( dx );
            m_dyCSideH1[cluster.cell.area()]->fill( dy );
          }
          m_eOverPH1[cluster.cell.area()]->fill( cluster.pos.e() * state.qOverP() );
          // compute the z-coordinate for which sqrt(dx^2+dy^2) is minimal
          double tx  = state.tx();
          double ty  = state.ty();
          double ddz = ( tx * dx + ty * dy ) / ( tx * tx + ty * ty );
          m_zH1[cluster.cell.area()]->fill( cluster.pos.z() + ddz );
          if ( std::abs( dy ) < 200 && std::abs( dx ) < 100 ) {
            m_dxVsXPr->fill( xtrack, dx );
            m_dxVsYPr->fill( ytrack, dx );
            m_dxVsTxPr->fill( state.tx(), dx );
            m_dxVsTyPr->fill( state.ty(), dx );
            m_dyVsXPr->fill( xtrack, dy );
            m_dyVsYPr->fill( ytrack, dy );
            m_dyVsTxPr->fill( state.tx(), dy );
            m_dyVsTyPr->fill( state.ty(), dy );
          }
        }
        if ( velostate && std::abs( dx ) < 200 ) {
          ytrack = velostate->y() + velostate->ty() * dz;
          dy     = cluster.pos.y() - ytrack;
          if ( cluster.pos.x() > 0 ) {
            m_dyVeloASideH1[cluster.cell.area()]->fill( dy );
          } else {
            m_dyVeloCSideH1[cluster.cell.area()]->fill( dy );
          }
        }
      }
    }
  }
}
