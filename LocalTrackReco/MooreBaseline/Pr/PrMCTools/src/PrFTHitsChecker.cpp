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
#include "Event/FTLiteCluster.h"
#include "Event/LinksByKey.h"
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "FTDet/DeFTDetector.h"
#include "Gaudi/Accumulators/Histogram.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Consumer.h"
#include "PrKernel/PrFTHitHandler.h"
//-----------------------------------------------------------------------------
// Implementation file for class : PrFTHitsChecker
//
// 2021-05: valeriia lukashenko
//-----------------------------------------------------------------------------
// @class PrFTHitsChecker PrFTHitsChecker.h
// 	    * class for plotting redisuals from PrFTHits to MCHits
//
// 	    * isGlobal - if true use global coordinates, if false use local coordinates
// 	    * performStudy  - if true plot residuals, residulas for different FTLiteCluster sizes,
// residuals vs track angle. 	   		      if false just plot residuals
// 	    * min is the left side of the residuals range
// 	    * max is the right side of the residuals range
// 	    * nbins is number of the bins for residuals
//          * to be used in the Dashboard
//          * example: Moore/Hlt/RecoCond/options/tracking_developments/mc_hit_resolution_monitor.py
//
// *  @author valeriia lukashenko
// *  @date   2021-05-04
//

namespace {
  using MCHits         = LHCb::MCHits;
  using MCHit          = LHCb::MCHit;
  using MCParticles    = LHCb::MCParticles;
  using FTHits         = PrFTHitHandler<PrHit>;
  using FTLiteClusters = LHCb::FTLiteCluster::FTLiteClusters;
} // namespace

class PrFTHitsChecker final
    : public LHCb::Algorithm::Consumer<void( const FTHits&, const FTLiteClusters&, const LHCb::LinksByKey&,
                                             const MCHits&, const MCParticles&, const DeFT& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DeFT>> {
public:
  PrFTHitsChecker( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"FTHitsLocation", "FT/FTHits"},
                   KeyValue{"FTLiteClusterLocation", LHCb::FTLiteClusterLocation::Default},
                   KeyValue{"FTHits2MCHitLinksLocation", LHCb::FTLiteClusterLocation::Default + "2MCHits"},
                   KeyValue{"MCHitsLocation", LHCb::MCHitLocation::FT},
                   KeyValue{"MCParticleLocation", LHCb::MCParticleLocation::Default},
                   KeyValue{"DeFT", DeFTDetectorLocation::Default}} ) {}

  void operator()( const FTHits&, const FTLiteClusters&, const LHCb::LinksByKey&, const MCHits&, const MCParticles&,
                   const DeFT& ) const override;

private:
  Gaudi::Property<bool>   m_isGlobal{this, "isGlobal", true, "Boolean for choosing global coordinate system"};
  Gaudi::Property<bool>   m_performStudy{this, "performStudy", false,
                                       "Boolean for performing studies with cluster size and diraction of mcparticle"};
  Gaudi::Property<double> m_min{this, "m_min_bin", -1., "Double for left boundary of residuals range"};
  Gaudi::Property<double> m_max{this, "m_max_bin", 1., "Double for right boundary of residuals range"};

  using ErrorCounter = Gaudi::Accumulators::MsgCounter<MSG::ERROR>;
  mutable ErrorCounter          m_noClus{this, "No cluster found"};
  mutable ErrorCounter          m_noMat{this, "No mat found"};
  mutable ErrorCounter          m_noMCP{this, "No MCParticle found"};
  Gaudi::Property<unsigned int> m_num_bins{this, "m_num_bins", 100, "Int for number of bins"};
  const std::string             geom_name = m_isGlobal ? "" : "_local";

  mutable Gaudi::Accumulators::Histogram<1> m_res_x{
      this, "residuals" + geom_name + "_X", "residuals" + geom_name + "_X", {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_0{this,
                                                      "residuals" + geom_name + "_X_cluster_size_bigger_than_8",
                                                      "residuals" + geom_name + "_X_cluster_size_bigger_than_8",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_1{this,
                                                      "residuals" + geom_name + "_X_cluster_size_between_1_4",
                                                      "residuals" + geom_name + "_X_cluster_size_between_1_4",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_5{this,
                                                      "residuals" + geom_name + "_X_cluster_size_5",
                                                      "residuals" + geom_name + "_X_cluster_size_5",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_6{this,
                                                      "residuals" + geom_name + "_X_cluster_size_6",
                                                      "residuals" + geom_name + "_X_cluster_size_6",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_7{this,
                                                      "residuals" + geom_name + "_X_cluster_size_7",
                                                      "residuals" + geom_name + "_X_cluster_size_7",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_x_8{this,
                                                      "residuals" + geom_name + "_X_cluster_size_8",
                                                      "residuals" + geom_name + "_X_cluster_size_8",
                                                      {m_num_bins * 10, m_min, m_max}};

  mutable Gaudi::Accumulators::Histogram<1> m_res_z{
      this, "residuals" + geom_name + "_Z", "residuals" + geom_name + "_Z", {m_num_bins * 10, m_min * 20, m_max * 20}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_0{this,
                                                      "residuals" + geom_name + "_Z_cluster_size_bigger_than_8",
                                                      "residuals" + geom_name + "_Z_cluster_size_bigger_than_8",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_1{this,
                                                      "residuals" + geom_name + "_Z_cluster_size_between_1_4",
                                                      "residuals" + geom_name + "_Z_cluster_size_between_1_4",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_5{this,
                                                      "residuals" + geom_name + "_Z_cluster_size_5",
                                                      "residuals" + geom_name + "_Z_cluster_size_5",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_6{this,
                                                      "residuals" + geom_name + "_Z_cluster_size_6",
                                                      "residuals" + geom_name + "_Z_cluster_size_6",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_7{this,
                                                      "residuals" + geom_name + "_Z_cluster_size_7",
                                                      "residuals" + geom_name + "_Z_cluster_size_7",
                                                      {m_num_bins * 10, m_min, m_max}};
  mutable Gaudi::Accumulators::Histogram<1> m_res_z_8{this,
                                                      "residuals" + geom_name + "_Z_cluster_size_8",
                                                      "residuals" + geom_name + "_Z_cluster_size_8",
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

DECLARE_COMPONENT( PrFTHitsChecker )

void PrFTHitsChecker::operator()( const FTHits& FThits, const FTLiteClusters& clusters, const LHCb::LinksByKey& links,
                                  const MCHits& mchits, const MCParticles&, const DeFT& DeFT ) const {

  std::map<const unsigned int, std::vector<LHCb::MCHit const*>> mcHitForId;
  links.applyToAllLinks( [&mcHitForId, &mchits]( unsigned int id, unsigned int mcHitKey, float ) {
    mcHitForId[id].emplace_back( mchits[mcHitKey] );
  } );

  if ( !m_performStudy ) {
    for ( int i = 0; i != LHCb::Detector::FT::nbZones(); ++i ) {

      for ( const auto hit : FThits.hits( i ) ) {

        for ( auto& mcHit : ( *mcHitForId.find( hit.id().ftID() ) ).second ) {

          const Gaudi::XYZPoint FThit =
              Gaudi::XYZPoint( hit.x( mcHit->midPoint().Y() ), 0, hit.z( mcHit->midPoint().Y() ) );
          const Gaudi::XYZVector residual = FThit - ( mcHit->midPoint() );
          ++m_res_x[residual.X()];
          ++m_res_z[residual.Z()];
        } // loop mcHits
      }   // loop over hits
    }     // loop over zones
  }       // bool performStudy

  if ( m_performStudy ) {
    for ( int i = 0; i != LHCb::Detector::FT::nbZones(); ++i ) {
      for ( const auto hit : FThits.hits( i ) ) {
        const auto& r = clusters.range();
        auto        it =
            find_if( r.begin(), r.end(), [id = hit.id().ftID()]( const auto& h ) { return h.channelID() == id; } );
        if ( it == r.end() ) {
          debug() << "No cluster is found for " << hit.id().ftID() << endmsg;
          ++m_noClus;
          continue;
        }

        const auto cluster_size = it->pseudoSize();

        auto mat = DeFT.findMat( hit.id().ftID() );

        if ( !mat ) {
          debug() << "No mat is found for hit (" << hit.x() << ", 0, " << hit.z() << ")" << endmsg;
          ++m_noMat;
          continue;
        }

        auto result = mcHitForId.find( hit.id().ftID() );
        if ( result == mcHitForId.end() ) continue;

        for ( auto& mcHit : ( *result ).second ) {
          auto             mcHit_local = mat->toLocal( mcHit->midPoint() ); // FIXME
          const auto       FThit = Gaudi::XYZPoint( hit.x( mcHit->midPoint().Y() ), 0, hit.z( mcHit->midPoint().Y() ) );
          auto             FThit_local = mat->toLocal( FThit ); // FIXME
          Gaudi::XYZVector residual    = m_isGlobal ? FThit - ( mcHit->midPoint() ) : FThit_local - ( mcHit_local );

          ++m_cluster[cluster_size];

          ++m_res_x[residual.X()];
          ++m_res_z[residual.Z()];

          const auto MCparticle = mcHit->mcParticle();
          if ( !MCparticle ) {
            debug() << "No mcParticle is found for the mcHit" << endmsg;
            ++m_noMCP;
            continue;
          }
          const auto& momentum = MCparticle->momentum();
          if ( momentum.Z() != 0. ) {
            m_tprof_res_x_vs_slope_dxdz[momentum.X() / momentum.Z()] += residual.X();
            ++m_res_vs_p_x[{residual.X(), momentum.X() / momentum.Z()}];
          } else {
            m_tprof_res_x_vs_slope_dxdz[-9999.] += -9999.;
            ++m_res_vs_p_x[{-9999., -9999.}];
          }

          ++m_res_x_vs_cluster_size[{residual.X(), cluster_size}];
          ++m_res_z_vs_cluster_size[{residual.Z(), cluster_size}];

          if ( cluster_size == 4 ) {
            ++m_res_x_1[residual.X()];
            ++m_res_z_1[residual.Z()];
          } else if ( cluster_size == 5 ) {
            ++m_res_x_5[residual.X()];
            ++m_res_z_5[residual.Z()];
          } else if ( cluster_size == 6 ) {
            ++m_res_x_6[residual.X()];
            ++m_res_z_6[residual.Z()];
          } else if ( cluster_size == 7 ) {
            ++m_res_x_7[residual.X()];
            ++m_res_z_7[residual.Z()];
          } else if ( cluster_size == 8 ) {
            ++m_res_x_8[residual.X()];
            ++m_res_z_8[residual.Z()];
          } else if ( cluster_size == 0 ) {
            ++m_res_x_0[residual.X()];
            ++m_res_z_0[residual.Z()];
          }
        } // while MChit
      }   // loop over hits
    }     // loop over xzones
  }       // bool perfromStudy
} // void
