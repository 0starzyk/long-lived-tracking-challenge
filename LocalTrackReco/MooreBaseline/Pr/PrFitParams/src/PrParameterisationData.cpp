/*****************************************************************************\
* (c) Copyright 2022 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "Detector/FT/FTConstants.h"
#include "Magnet/DeMagnet.h"

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Event/StateVector.h"
#include "FTDet/DeFTDetector.h"
#include "LHCbAlgs/Consumer.h"

#include "PrFitParams/IPrFitTool.h"
#include "PrKernel/IPrDebugTrackingTool.h"
#include "PrKernel/PrFTZoneHandler.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

/**
 * This class can be used to produce an NTuple containing magnetic field dependend quantities
 * that need to be parameterised for the pattern recognition.
 *
 * See https://gitlab.cern.ch/gunther/prforwardtracking-parametrisation-tuner/-/tree/master for
 * more information.
 *
 */

namespace {
  using FTZoneCache::ZoneCache;
  using namespace LHCb::Detector::FT;
} // namespace

namespace LHCb::Pr {
  class PrParameterisationData
      : public Algorithm::Consumer<void( const MCParticles&, const MCHits&, const MCHits&, const DetectorElement&,
                                         const ZoneCache&, const DeMagnet& ),
                                   DetDesc::usesConditions<DetectorElement, ZoneCache, DeMagnet>> {
  public:
    PrParameterisationData( const std::string& name, ISvcLocator* pSvcLocator )
        : Consumer( name, pSvcLocator,
                    {KeyValue{"MCParticles", ""}, KeyValue{"MCVPHits", ""}, KeyValue{"MCFTHits", ""},
                     KeyValue{"StandardGeometryTop", standard_geometry_top},
#ifdef USE_DD4HEP
                     KeyValue{"FTZoneCache", "/world:AlgorithmSpecific-" + name + "-FTZoneCache"},
#else
                     KeyValue{"FTZoneCache", "AlgorithmSpecific-" + name + "-FTZoneCache"},
#endif
                     KeyValue{"Magnet", Det::Magnet::det_path}} ){};

    using base_class_t = LHCb::Algorithm::Consumer<void( const MCParticles&, const MCHits&, const MCHits&,
                                                         const DetectorElement&, const ZoneCache&, const DeMagnet& ),
                                                   DetDesc::usesConditions<DetectorElement, ZoneCache, DeMagnet>>;

    StatusCode initialize() override {
      auto sc = base_class_t::initialize();
      if ( sc.isFailure() ) return sc;
      addConditionDerivation<ZoneCache( const DeFT& )>( {DeFTDetectorLocation::Default},
                                                        this->template inputLocation<ZoneCache>() );
      return sc;
    }

    void operator()( const MCParticles&, const MCHits&, const MCHits&, const DetectorElement&, const ZoneCache&,
                     const DeMagnet& ) const override;

  private:
    Gaudi::Property<double> m_zEndVelo{this, "zEndVelo", 770. * Gaudi::Units::mm};
    Gaudi::Property<double> m_zRef{this, "zRef", 8520. * Gaudi::Units::mm};
    Gaudi::Property<double> m_zOutOfMagField{this, "zOutOfMagField", 10. * Gaudi::Units::m};

    ToolHandle<IPrDebugTrackingTool> m_debugTool{this, "DebugTool", "PrMCDebugReconstructibleLong"};
    ToolHandle<IPrFitTool>           m_fitTool{this, "FitTool", "PrFitTool"};
    ToolHandle<ITrackExtrapolator>   m_extrap{this, "Extrapolator", "TrackRungeKuttaExtrapolator"};
  };

  DECLARE_COMPONENT_WITH_ID( PrParameterisationData, "PrParameterisationData" )

  void PrParameterisationData::operator()( const MCParticles& mc_particles, const MCHits& mc_vphits,
                                           const MCHits& mc_fthits, const DetectorElement& lhcb,
                                           const ZoneCache& ft_cache, const DeMagnet& magnet ) const {

    for ( auto mcp : mc_particles ) {
      const auto mcp_key = mcp->index();
      // check if long reconstructible (3 velo hits + one x and one uv hit in every SciFi station)
      if ( !m_debugTool->check( mcp_key ) ) continue;

      // now copy all MC velo hits that belong to this MC particle (MCHit and positions)
      std::vector<MCHit*>          velo_hits{};
      std::vector<Gaudi::XYZPoint> velo_hit_positions{};
      velo_hits.reserve( 50 );
      velo_hit_positions.reserve( 50 );
      std::copy_if( mc_vphits.begin(), mc_vphits.end(), std::back_inserter( velo_hits ),
                    [&mcp]( auto hit ) { return mcp == hit->mcParticle(); } );
      std::transform( velo_hits.begin(), velo_hits.end(), std::back_inserter( velo_hit_positions ),
                      []( auto hit ) { return hit->midPoint(); } );

      // now copy all MC scifi hits that belong to this MC particle (MCHit and positions)
      std::vector<MCHit*> scifi_hits{};
      scifi_hits.reserve( 16 );
      std::vector<Gaudi::XYZPoint> scifi_hit_positions{};
      scifi_hit_positions.reserve( 16 );
      std::copy_if( mc_fthits.begin(), mc_fthits.end(), std::back_inserter( scifi_hits ),
                    [&mcp]( auto hit ) { return mcp == hit->mcParticle(); } );
      std::transform( scifi_hits.begin(), scifi_hits.end(), std::back_inserter( scifi_hit_positions ),
                      []( auto hit ) { return hit->midPoint(); } );

      // to not depend on the velo tracking, the velo hits are simply fitted with a straight line without errors
      const auto velo_x_result = m_fitTool->fitLine( velo_hit_positions, IPrFitTool::XY::X, m_zEndVelo );
      const auto velo_y_result = m_fitTool->fitLine( velo_hit_positions, IPrFitTool::XY::Y, m_zEndVelo );
      assert( velo_x_result && velo_y_result );
      const auto [x0, tx] = velo_x_result.value();
      const auto [y0, ty] = velo_y_result.value();

      // A track is defined as x = AX + BX*dz + CX*dz*dz + DX*dz*dz*dz
      //                       y = AY + BY*dz + CY*dz*dz
      const auto trFitX = m_fitTool->fitCubic( scifi_hit_positions, IPrFitTool::XY::X, m_zRef );
      assert( trFitX );
      // these fit parameters also contain scattering, later we will fit also with "true" extrapolated positions
      const auto [AX, BX, CX, DX] = trFitX.value();
      const auto trFitY           = m_fitTool->fitParabola( scifi_hit_positions, IPrFitTool::XY::Y, m_zRef );
      assert( trFitY );
      const auto [AY, BY, CY] = trFitY.value();

      const auto [chi2_x, chi2_y, chi2_comb] = [&, AX = AX, BX = BX, CX = CX, DX = DX, AY = AY, BY = BY, CY = CY,
                                                x_err = 0.1 * Gaudi::Units::mm, y_err = 1. * Gaudi::Units::mm] {
        auto chi2_x{0.};
        auto chi2_y{0.};
        auto chi2_comb{0.};
        for ( auto pos : scifi_hit_positions ) {
          const auto dz = pos.z() - m_zRef;
          const auto dx = pos.x() - ( AX + dz * ( BX + dz * ( CX + dz * DX ) ) );
          const auto dy = pos.y() - ( AY + dz * ( BY + dz * CY ) );
          chi2_x += dx * dx / ( x_err * x_err );
          chi2_y += dy * dy / ( y_err * y_err );
          chi2_comb += dx * dx / ( x_err * x_err ) + dy * dy / ( y_err * y_err );
        }
        return std::tuple{chi2_x, chi2_y, chi2_comb};
      }();
      const int ndof_x = scifi_hit_positions.size() - 4;
      const int ndof_y = scifi_hit_positions.size() - 3;

      // let's get some info about the particle we are dealing with
      const auto momentum            = mcp->p();
      const auto transverse_momentum = mcp->pt();
      const auto qop                 = ( mcp->particleID().threeCharge() / 3 ) / momentum;
      const auto pid                 = std::abs( mcp->particleID().pid() );

      // propagate velo state to reference plane
      auto track_state = StateVector{{x0, y0, tx, ty, qop}, m_zEndVelo};
      if ( m_extrap->propagate( track_state, m_zRef, *lhcb.geometry() ).isFailure() ) { continue; }
      // fix the state values at the reference plane
      const auto x_ref  = track_state.x();
      const auto y_ref  = track_state.y();
      const auto z_ref  = m_zRef;
      const auto tx_ref = track_state.tx();
      const auto ty_ref = track_state.ty();
      // propagate to behind SciFi stations outside of magnetic field
      if ( m_extrap->propagate( track_state, m_zOutOfMagField, *lhcb.geometry() ).isFailure() ) { continue; }
      const auto x_out      = track_state.x();
      const auto y_out      = track_state.y();
      const auto z_out      = m_zOutOfMagField;
      const auto tx_out     = track_state.tx();
      const auto ty_out     = track_state.ty();
      const auto z_mag_x    = ( x_out - x0 - tx_out * m_zOutOfMagField + tx * m_zEndVelo ) / ( tx - tx_out );
      const auto z_mag_y    = ( y_out - y0 - ty_out * m_zOutOfMagField + ty * m_zEndVelo ) / ( ty - ty_out );
      const auto dSlope_out = tx_out - tx;
      // now prepare stuff for tupling
      std::vector<double> scifi_hit_x_positions{};
      scifi_hit_x_positions.reserve( 16 );
      std::vector<double> scifi_hit_y_positions{};
      scifi_hit_y_positions.reserve( 16 );
      std::vector<double> scifi_hit_z_positions{};
      scifi_hit_z_positions.reserve( 16 );
      std::transform( scifi_hit_positions.begin(), scifi_hit_positions.end(),
                      std::back_inserter( scifi_hit_x_positions ), []( auto pos ) { return pos.x(); } );
      std::transform( scifi_hit_positions.begin(), scifi_hit_positions.end(),
                      std::back_inserter( scifi_hit_y_positions ), []( auto pos ) { return pos.y(); } );
      std::transform( scifi_hit_positions.begin(), scifi_hit_positions.end(),
                      std::back_inserter( scifi_hit_z_positions ), []( auto pos ) { return pos.z(); } );
      const auto                                     magScaleFactor = magnet.signedRelativeCurrent();
      std::vector<IPrDebugTrackingTool::VariableDef> variables{{"x", x0},
                                                               {"y", y0},
                                                               {"z", m_zEndVelo},
                                                               {"tx", tx},
                                                               {"ty", ty},
                                                               {"p", momentum},
                                                               {"pt", transverse_momentum},
                                                               {"qop", qop},
                                                               {"x_ref", x_ref},
                                                               {"y_ref", y_ref},
                                                               {"z_ref", z_ref},
                                                               {"tx_ref", tx_ref},
                                                               {"ty_ref", ty_ref},
                                                               {"x_out", x_out},
                                                               {"y_out", y_out},
                                                               {"z_out", z_out},
                                                               {"tx_out", tx_out},
                                                               {"ty_out", ty_out},
                                                               {"z_mag_x", z_mag_x},
                                                               {"z_mag_y", z_mag_y},
                                                               {"dSlope_out", dSlope_out},
                                                               {"AX", AX},
                                                               {"BX", BX},
                                                               {"CX", CX},
                                                               {"DX", DX},
                                                               {"AY", AY},
                                                               {"BY", BY},
                                                               {"CY", CY},
                                                               {"chi2_x", chi2_x},
                                                               {"chi2_y", chi2_y},
                                                               {"chi2_comb", chi2_comb},
                                                               {"ndof_x", ndof_x},
                                                               {"ndof_y", ndof_y},
                                                               {"pid", pid},
                                                               {"scifi_hit_x", scifi_hit_x_positions},
                                                               {"scifi_hit_y", scifi_hit_y_positions},
                                                               {"scifi_hit_z", scifi_hit_z_positions},
                                                               {"signed_rel_current", magScaleFactor}};

      // propagate to each scifi layer
      std::array<std::pair<std::string, double>, 5 * NFTLayers> layer_states;
      std::vector<Gaudi::XYZPoint>                              extrapolated_positions( NFTLayers );
      if ( [&] {
             // the zone z positions are taken according to the half the track passes through
             const auto upper_half_track = y_ref > 0;
             for ( auto iLayer{0}; iLayer < NFTLayers; ++iLayer ) {
               const auto iZone = 2 * iLayer + upper_half_track;
               const auto zZone = ft_cache.zoneZPos[iZone];
               if ( m_extrap->propagate( track_state, zZone, *lhcb.geometry() ).isFailure() ) { return true; }
               layer_states[iLayer * 5]       = {"x_l" + std::to_string( iLayer ), track_state.x()};
               layer_states[iLayer * 5 + 1]   = {"y_l" + std::to_string( iLayer ), track_state.y()};
               layer_states[iLayer * 5 + 2]   = {"z_l" + std::to_string( iLayer ), track_state.z()};
               layer_states[iLayer * 5 + 3]   = {"tx_l" + std::to_string( iLayer ), track_state.tx()};
               layer_states[iLayer * 5 + 4]   = {"ty_l" + std::to_string( iLayer ), track_state.ty()};
               extrapolated_positions[iLayer] = track_state.position();
             }
             return false;
           }() ) {
        continue;
      }
      variables.insert( variables.end(), layer_states.begin(), layer_states.end() );

      // these fits now use "true" positions without scattering, which is very useful for parameterisations
      const auto trFitX_extrapolated = m_fitTool->fitCubic( extrapolated_positions, IPrFitTool::XY::X, m_zRef );
      assert( trFitX_extrapolated );
      const auto [AX_ex, BX_ex, CX_ex, DX_ex] = trFitX_extrapolated.value();
      const auto trFitY_extrapolated = m_fitTool->fitParabola( extrapolated_positions, IPrFitTool::XY::Y, m_zRef );
      assert( trFitY_extrapolated );
      const auto [AY_ex, BY_ex, CY_ex] = trFitY_extrapolated.value();
      std::array<std::pair<std::string, double>, 7> coef_from_extrapolation{{{"AX_ex", AX_ex},
                                                                             {"BX_ex", BX_ex},
                                                                             {"CX_ex", CX_ex},
                                                                             {"DX_ex", DX_ex},
                                                                             {"AY_ex", AY_ex},
                                                                             {"BY_ex", BY_ex},
                                                                             {"CY_ex", CY_ex}}};
      variables.insert( variables.end(), coef_from_extrapolation.begin(), coef_from_extrapolation.end() );

      m_debugTool->storeData( variables );
    }
  }

} // namespace LHCb::Pr
