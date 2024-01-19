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
#include "Event/MCParticle.h"
#include "Event/PrHits.h"
#include "Event/UTCluster.h"
#include "Gaudi/Accumulators/Histogram.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Consumer.h"
#include "UTDAQ/UTInfo.h"
#include "UTDet/DeUTDetector.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrUTHitsChecker
//
// 2021-05: valeriia lukashenko
//-----------------------------------------------------------------------------
// @class PrUTHitsChecker PrUTHitsChecker.h
// 	    * class for plotting redisuals from UTHits to MCHits
//
// 	    * isGlobal - if true use global coordinates, if false use local coordinates
// 	    * performStudy  - if true plot residuals, residulas in different planes,
// residuals vs track angle. 	   		      if false just plot residuals
// 	    * min is the left side of the residuals range
// 	    * max is the right side of the residuals range
// 	    * nbins is number of the bins for residuals
//          * to be used in the Dashboard
//          * example: Moore/Hlt/RecoCond/options/tracking_developments/mc_hit_resolution_monitor.py
//
//
// *  @author valeriia lukashenko
// *  @date   2021-05-04
//
namespace {
  using MCHits      = LHCb::MCHits;
  using MCHit       = LHCb::MCHit;
  using MCParticles = LHCb::MCParticles;
  using UTHits      = LHCb::Pr::UT::Hits;
  using simd        = SIMDWrapper::scalar::types;
} // namespace
class PrUTHitsChecker final
    : public LHCb::Algorithm::Consumer<void( const UTHits&, const LHCb::LinksByKey&, const MCHits&, const MCParticles&,
                                             const DeUTDetector& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DeUTDetector>> {
public:
  PrUTHitsChecker( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"UTHitsLocation", UTInfo::HitLocation},
                   KeyValue{"UTHits2MCHitLinksLocation", LHCb::UTClusterLocation::UTClusters + "2MCHits"},
                   KeyValue{"MCHitsLocation", LHCb::MCHitLocation::UT},
                   KeyValue{"MCParticleLocation", LHCb::MCParticleLocation::Default},
                   KeyValue{"DeUT", DeUTDetLocation::location()}} ) {}

  void operator()( const UTHits&, const LHCb::LinksByKey&, const MCHits&, const MCParticles&,
                   const DeUTDetector& ) const override;

private:
  Gaudi::Property<bool>         m_isGlobal{this, "isGlobal", true, "Boolean for choosing global coordinate system"};
  Gaudi::Property<bool>         m_performStudy{this, "performStudy", false,
                                       "Boolean for performing studies with cluster size and diraction of mcparticle"};
  Gaudi::Property<double>       m_min{this, "m_min_bin", -1., "Double for left boundary of residuals range"};
  Gaudi::Property<double>       m_max{this, "m_max_bin", 1., "Double for right boundary of residuals range"};
  Gaudi::Property<unsigned int> m_num_bins{this, "m_num_bins", 100, "Int for number of bins"};
  const std::string             geom_name = m_isGlobal ? "" : "_local";

  mutable Gaudi::Accumulators::Histogram<1> m_res_x{
      this, "residuals" + geom_name + "_X", "residuals" + geom_name + "_X", {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_zoom{
      this, "residuals" + geom_name + "_X_zoom", "residuals" + geom_name + "_X_zoom", {m_num_bins * 10, -0.5, 0.5}};

  mutable Gaudi::Accumulators::Histogram<1> m_res_x_0{this,
                                                      "residuals" + geom_name + "_X_plane_0",
                                                      "residuals" + geom_name + "_X_plane_0",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_1{this,
                                                      "residuals" + geom_name + "_X_plane_1",
                                                      "residuals" + geom_name + "_X_plane_1",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_2{this,
                                                      "residuals" + geom_name + "_X_plane_2",
                                                      "residuals" + geom_name + "_X_plane_2",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_3{this,
                                                      "residuals" + geom_name + "_X_plane_3",
                                                      "residuals" + geom_name + "_X_plane_3",
                                                      {m_num_bins * 10, m_min, m_max}};

  mutable Gaudi::Accumulators::Histogram<1> m_res_z{
      this, "residuals" + geom_name + "_Z", "residuals" + geom_name + "_Z", {m_num_bins * 10, m_min * 20, m_max * 20}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_0{this,
                                                      "residuals" + geom_name + "_Z_plane_0",
                                                      "residuals" + geom_name + "_Z_plane_0",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_1{this,
                                                      "residuals" + geom_name + "_Z_plane_1",
                                                      "residuals" + geom_name + "_Z_plane_1",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_2{this,
                                                      "residuals" + geom_name + "_Z_plane_2",
                                                      "residuals" + geom_name + "_Z_plane_2",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_3{this,
                                                      "residuals" + geom_name + "_Z_plane_3",
                                                      "residuals" + geom_name + "_Z_plane_3",
                                                      {m_num_bins * 10, m_min, m_max}};

  mutable Gaudi::Accumulators::Histogram<2> m_res_vs_p_x{this,
                                                         "residuals" + geom_name + "_X_VS_slope_dxdz",
                                                         "residuals" + geom_name + "_X_VS_slope_dxdz",
                                                         {{m_num_bins, m_min, m_max}, {m_num_bins, -1., 1.}}};
  mutable Gaudi::Accumulators::Histogram<1> m_cluster{
      this, "cluster_size", "cluster_size; cluster size; counts", {11, -0.5, 10.5}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_tprof_res_x_vs_slope_dxdz{
      this,
      "profile_res" + geom_name + "_X_VS_slope_dxdz",
      "profile_res" + geom_name + "_X_VS_slope_dxdz",
      {m_num_bins, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<2> m_res_x_vs_cluster_size{this,
                                                                    "residuals" + geom_name + "_X_VS_cluster_size",
                                                                    "residuals" + geom_name + "_X_VS_cluster_size",
                                                                    {{m_num_bins, m_min, m_max}, {11, -0.5, 10.5}}};
  mutable Gaudi::Accumulators::Histogram<2> m_res_z_vs_cluster_size{this,
                                                                    "residuals" + geom_name + "_Z_VS_cluster_size",
                                                                    "residuals" + geom_name + "_Z_VS_cluster_size",
                                                                    {{m_num_bins, m_min, m_max}, {11, -0.5, 10.5}}};
};

DECLARE_COMPONENT( PrUTHitsChecker )

void PrUTHitsChecker::operator()( const UTHits& UThits, const LHCb::LinksByKey& links, const MCHits& mchits,
                                  const MCParticles&, const DeUTDetector& DeUT ) const {

  std::map<const unsigned int, std::vector<LHCb::MCHit const*>> mcHitForId;
  links.applyToAllLinks( [&mcHitForId, &mchits]( unsigned int id, unsigned int mcHitKey, float ) {
    mcHitForId[id].emplace_back( mchits[mcHitKey] );
  } );

  // Map for MCHits and keys
  if ( !m_performStudy ) {
    const int fullChanIdx =
        static_cast<int>( UTInfo::DetectorNumbers::Layers ) * static_cast<int>( UTInfo::DetectorNumbers::Stations ) *
        static_cast<int>( UTInfo::DetectorNumbers::Regions ) * static_cast<int>( UTInfo::DetectorNumbers::Sectors );
    // new SOA for UTHits: hit.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>
    for ( int fullchan = 0; fullchan < fullChanIdx; fullchan++ ) {
      const auto indexs = UThits.indices( fullchan );

      for ( int i = indexs.first; i != indexs.second; i++ ) {
        const auto        hit           = UThits.scalar()[i];
        const simd::int_v simd_chid     = hit.get<LHCb::Pr::UT::UTHitsTag::channelID>();
        const int         unwrapped_num = simd_chid.cast();

        for ( auto& mcHit : ( *mcHitForId.find( unwrapped_num ) ).second ) {
          const auto hit_x = hit.get<LHCb::Pr::UT::UTHitsTag::xAtYEq0>().cast() +
                             mcHit->midPoint().Y() * hit.get<LHCb::Pr::UT::UTHitsTag::dxDy>().cast();
          Gaudi::XYZPoint UThit = Gaudi::XYZPoint( hit_x, 0, hit.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>().cast() );

          auto residual = UThit - ( mcHit->midPoint() );
          ++m_res_x[residual.X()];
          ++m_res_z[residual.Z()];
        } // while mcHit
      }   // hit loop
    }     // channel id loop
  }       // if not  performStudy
  if ( m_performStudy ) {

    const int fullChanIdx =
        static_cast<int>( UTInfo::DetectorNumbers::Layers ) * static_cast<int>( UTInfo::DetectorNumbers::Stations ) *
        static_cast<int>( UTInfo::DetectorNumbers::Regions ) * static_cast<int>( UTInfo::DetectorNumbers::Sectors );

    for ( int fullchan = 0; fullchan < fullChanIdx; fullchan++ ) {
      const auto indexs = UThits.indices( fullchan );

      for ( int i = indexs.first; i != indexs.second; i++ ) {
        const auto        hit           = UThits.scalar()[i];
        const simd::int_v simd_chid     = hit.get<LHCb::Pr::UT::UTHitsTag::channelID>().cast();
        const int         unwrapped_num = simd_chid.cast();
        auto              sector        = DeUT.findSector( LHCb::Detector::UT::ChannelID( unwrapped_num ) );
        auto              result        = ( mcHitForId.find( unwrapped_num ) );
        if ( result == mcHitForId.end() ) continue;
        for ( auto& mcHit : ( *result ).second ) {
          const auto hit_x = hit.get<LHCb::Pr::UT::UTHitsTag::xAtYEq0>().cast() +
                             mcHit->midPoint().Y() * hit.get<LHCb::Pr::UT::UTHitsTag::dxDy>().cast();
          Gaudi::XYZPoint UThit = Gaudi::XYZPoint( hit_x, 0, hit.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>().cast() );
          LHCb::Detector::UT::ChannelID utid = LHCb::Detector::UT::ChannelID( (unsigned int)unwrapped_num );
          if ( !sector )
            error() << "No sector is found for "
                    << "(" << hit.get<LHCb::Pr::UT::UTHitsTag::xAtYEq0>().cast() << ", 0, "
                    << hit.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>().cast() << ")" << std::endl;
          const auto UThit_local = sector->toLocal( UThit );
          const auto mcHit_local = sector->toLocal( mcHit->midPoint() );

          const int        plane    = 2 * ( utid.station() - 1 ) + ( utid.layer() - 1 ) % 2;
          Gaudi::XYZVector residual = m_isGlobal ? UThit - ( mcHit->midPoint() ) : UThit_local - ( mcHit_local );
          ++m_res_x[residual.X()];
          ++m_res_z[residual.Z()];

          ++m_res_x_zoom[residual.X()];

          switch ( plane ) {
          case 0:
            ++m_res_x_0[residual.X()];
            ++m_res_z_0[residual.Z()];
            return;
          case 1:
            ++m_res_x_1[residual.X()];
            ++m_res_z_1[residual.Z()];
            return;
          case 2:
            ++m_res_x_2[residual.X()];
            ++m_res_z_2[residual.Z()];
            return;
          case 3:
            ++m_res_x_3[residual.X()];
            ++m_res_z_3[residual.Z()];
            return;
          }
          auto       MCparticle = mcHit->mcParticle();
          const auto momentum   = MCparticle->momentum();
          if ( momentum.Z() != 0 ) {
            m_tprof_res_x_vs_slope_dxdz[momentum.X() / momentum.Z()] += residual.X();
            ++m_res_vs_p_x[{residual.X(), momentum.X() / momentum.Z()}];

          } else {
            m_tprof_res_x_vs_slope_dxdz[-9999.] += -9999.;
            ++m_res_vs_p_x[{-9999., -9999.}];
          }
        } // loop mcHits
      }   // loop indices
    }     // loop over channels
  }       // if performStudy
}
