/*****************************************************************************\
 * * (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
 * *                                                                             *
 * * This software is distributed under the terms of the GNU General Public      *
 * * Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
 * *                                                                             *
 * * In applying this licence, CERN does not waive the privileges and immunities *
 * * granted to it by virtue of its status as an Intergovernmental Organization  *
 * * or submit itself to any jurisdiction.                                       *
 * \*****************************************************************************/

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "DetDesc/IConditionDerivationMgr.h"
#include "Event/LinksByKey.h"
#include "Event/MCHit.h"
#include "Event/PrHits.h"
#include "Event/VPFullCluster.h"
#include "Gaudi/Accumulators/Histogram.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Consumer.h"
#include "VPDet/DeVP.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrVPHitsChecker
//
// 2020-04: valeriia lukashenko
//-----------------------------------------------------------------------------
// @class PrVPHitsChecker PrVPHitsChecker.h
// 	    * class for plotting redisuals from PrVeloHits to MCHits
//
// 	    * isGlobal - if true use global coordinates, if false use local coordinates
// 	    * performStudy  - if true plot residuals, residulas vs VPCluster size, residuals per cluster size, residuals
// vs track angle. 	   		      if false just plot residuals
// 	    * min is the left side of the residuals range
// 	    * max is the right side of the residuals range
// 	    * nbins is number of the bins for residuals
//          * To be used in the Dashboard
//          * example: Moore/Hlt/RecoCond/options/tracking_developments/mc_hit_resolution_monitor.py
//
// *  @author valeriia lukashenko
// *  @date   2020-04-23
// *  reviewed for Dashboard 2021-05-04

namespace VPHitsTag = LHCb::Pr::VP::VPHitsTag;
namespace {
  using MCHits      = LHCb::MCHits;
  using MCParticles = LHCb::MCParticles;
  using MCHit       = LHCb::MCHit;
  using Hits        = LHCb::Pr::VP::Hits;
  using simd        = SIMDWrapper::scalar::types;
  using VPFullClus  = LHCb::VPFullCluster;
} // namespace
class PrVPHitsChecker final
    : public LHCb::Algorithm::Consumer<void( const Hits&, const std::vector<VPFullClus>&, const LHCb::LinksByKey&,
                                             const MCHits&, const MCParticles&, const DeVP& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DeVP>> {
public:
  PrVPHitsChecker( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer{name,
                 pSvcLocator,
                 {KeyValue{"VPHitsLocation", "Raw/VP/Hits"},
                  KeyValue{"VPFullClusterLocation", LHCb::VPFullClusterLocation::Default},
                  KeyValue{"VPHits2MCHitLinksLocation", "Link/Pr/LHCbID"},
                  KeyValue{"MCHitsLocation", LHCb::MCHitLocation::VP},
                  KeyValue{"MCParticleLocation", LHCb::MCParticleLocation::Default},
                  KeyValue{"DeVP", LHCb::Det::VP::det_path}}} {}

  void operator()( const Hits&, const std::vector<VPFullClus>&, const LHCb::LinksByKey&, const MCHits&,
                   const MCParticles&, const DeVP& ) const override;

private:
  Gaudi::Property<bool>         m_isGlobal{this, "isGlobal", true, "Boolean for choosing global coordinate system"};
  Gaudi::Property<bool>         m_performStudy{this, "performStudy", false,
                                       "Boolean for performing studies with cluster size and diraction of mcparticle"};
  Gaudi::Property<double>       m_min{this, "min_bin", -0.5, "Double for left boundary of residuals range"};
  Gaudi::Property<double>       m_max{this, "max_bin", 0.5, "Double for right boundary of residuals range"};
  Gaudi::Property<unsigned int> m_num_bins{this, "num_bins", 100, "Int for number of bins"};
  const std::string             geom_name = m_isGlobal ? "" : "_local";

  mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR>   m_nomchit{this, "No MCHit is found in the plane of VPHit"};
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_nomcparticle{this, "No MCParticle is found for MCHit"};
  mutable Gaudi::Accumulators::Histogram<1>             m_res_x{
      this, "residuals" + geom_name + "_X", "residuals" + geom_name + "_X", {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_y{
      this, "residuals" + geom_name + "_Y", "residuals" + geom_name + "_Y", {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z{
      this, "residuals" + geom_name + "_Z", "residuals" + geom_name + "_Z", {m_num_bins * 10, m_min, m_max}};

  mutable Gaudi::Accumulators::Histogram<2> res_vs_p_x{this,
                                                       "residuals" + geom_name + "_X_VS_slope_dxdz",
                                                       "residuals" + geom_name + "_X_VS_slope_dxdz",
                                                       {{m_num_bins, m_min, m_max}, {m_num_bins, -1., 1.}}};
  mutable Gaudi::Accumulators::Histogram<2> res_vs_p_y{this,
                                                       "residuals" + geom_name + "_Y_VS_slope_dydz",
                                                       "residuals" + geom_name + "_Y_VS_slope_dydz",
                                                       {{m_num_bins, m_min, m_max}, {m_num_bins, -1., 1.}}};
  mutable Gaudi::Accumulators::Histogram<1> cluster{
      this, "cluster_size", "cluster_size; cluster size; counts", {11, -0.5, 10.5}};
  mutable Gaudi::Accumulators::Histogram<2>        res_x_vs_cluster_size{this,
                                                                  "residuals" + geom_name + "_X_VS_cluster_size",
                                                                  "residuals" + geom_name + "_X_VS_cluster_size",
                                                                  {{m_num_bins, m_min, m_max}, {11, -0.5, 10.5}}};
  mutable Gaudi::Accumulators::Histogram<2>        res_y_vs_cluster_size{this,
                                                                  "residuals" + geom_name + "_Y_VS_cluster_size",
                                                                  "residuals" + geom_name + "_Y_VS_cluster_size",
                                                                  {{m_num_bins, m_min, m_max}, {11, -0.5, 10.5}}};
  mutable Gaudi::Accumulators::Histogram<2>        res_z_vs_cluster_size{this,
                                                                  "residuals" + geom_name + "_Z_VS_cluster_size",
                                                                  "residuals" + geom_name + "_Z_VS_cluster_size",
                                                                  {{m_num_bins, m_min, m_max}, {11, -0.5, 10.5}}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> tprof_res_x_vs_slope_dxdz{
      this,
      "profile_res" + geom_name + "_X_VS_slope_dxdz",
      "profile_res" + geom_name + "_X_VS_slope_dxdz",
      {m_num_bins, m_min, m_max}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> tprof_res_y_vs_slope_dydz{
      this,
      "profile_res" + geom_name + "_Y_VS_slope_dydz",
      "profile_res" + geom_name + "_Y_VS_slope_dydz",
      {m_num_bins, m_min, m_max}};
};

DECLARE_COMPONENT( PrVPHitsChecker )
void PrVPHitsChecker::operator()( const Hits& VPHits, const std::vector<VPFullClus>& clusters,
                                  const LHCb::LinksByKey& links, const MCHits& mchits, const MCParticles&,
                                  const DeVP& det ) const {
  if ( !m_performStudy ) { // in global coordinates

    for ( const auto& vphit : VPHits.scalar() ) {

      // get vp channel id
      const int unwrapped_num = vphit.get<VPHitsTag::ChannelId>().cast();

      // get vp hit position
      const auto            vec = vphit.get<VPHitsTag::pos>().vec3();
      const Gaudi::XYZPoint VPhit( vec.x.cast(), vec.y.cast(), vec.z.cast() );

      // get MCParticles that have the same VPChannelID
      LHCb::MCHit const* MChit{nullptr};
      float              max_weight{0};

      links.applyToLinks( unwrapped_num,
                          [&max_weight, &MChit, &mchits]( unsigned int, unsigned int mcHitKey, float weight ) {
                            if ( weight > max_weight ) MChit = mchits[mcHitKey];
                          } );

      if ( !MChit ) continue;

      const auto residual = VPhit - ( MChit->midPoint() ); // calcualte the residual
      ++m_res_x[residual.X()];
      ++m_res_y[residual.Y()];
      ++m_res_z[residual.Z()];

    } // loop over hits
  }   // if not perform study

  if ( m_performStudy ) {

    std::map<int, double> bin_edges;
    for ( int i = 0; i < 100; ++i ) { bin_edges[i] = -1. + i * 2. / 100; }

    for ( auto& ch : clusters ) {
      const unsigned int    channelID = ch.channelID();
      const Gaudi::XYZPoint VPCluster( ch.x(), ch.y(), ch.z() );
      const auto            id   = LHCb::LHCbID{LHCb::Detector::VPChannelID( channelID )};
      const auto&           sens = det.sensor( id.vpID() );

      LHCb::MCHit const* MChit{nullptr};
      float              max_weight{0};
      links.applyToLinks( channelID,
                          [&max_weight, &MChit, &mchits]( unsigned int, unsigned int mcHitKey, float weight ) {
                            if ( weight > max_weight ) MChit = mchits[mcHitKey];
                          } );
      if ( !MChit ) continue;
      const auto MCparticle   = MChit->mcParticle();
      const auto mcparticle_p = MCparticle->p();
      const auto momentum     = MCparticle->momentum();

      if ( mcparticle_p < 2000 ) continue;

      Gaudi::XYZVector residual = m_isGlobal
                                      ? VPCluster - ( MChit->midPoint() )
                                      : sens.globalToLocal( VPCluster ) - sens.globalToLocal( MChit->midPoint() );

      const auto pixels = ch.pixels();
      ++cluster[pixels.size()];

      ++m_res_x[residual.X()];
      ++m_res_y[residual.Y()];
      ++m_res_z[residual.Z()];

      ++res_x_vs_cluster_size[{residual.X(), pixels.size()}];
      ++res_y_vs_cluster_size[{residual.Y(), pixels.size()}];
      ++res_z_vs_cluster_size[{residual.Z(), pixels.size()}];

      ++res_vs_p_x[{residual.X(), momentum.X() / momentum.Z()}];
      ++res_vs_p_y[{residual.Y(), momentum.Y() / momentum.Z()}];

      tprof_res_x_vs_slope_dxdz[momentum.X() / momentum.Z()] += residual.X();
      tprof_res_y_vs_slope_dydz[momentum.Y() / momentum.Z()] += residual.Y();

    } // loop over clusters
  }   // if perform study
}
