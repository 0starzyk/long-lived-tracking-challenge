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
#include "DetDesc/Condition.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "LHCbAlgs/Transformer.h"

#include "Kernel/LHCbID.h"

#include "FTDAQ/FTInfo.h"
#include "SciFiTrackForwardingHits.h"

#include <array>
#include <limits>
#include <string>

#include "Event/FTLiteCluster.h"
#include "FTDet/DeFTDetector.h"

#include <boost/numeric/conversion/cast.hpp>

/** @class SciFiTrackForwardingStoreHit SciFiTrackForwardingStoreHit.cpp
 *
 *  \brief Transforms FTLiteClusters into the input format needed by the SciFiTrackForwarding
 */

namespace {
  using FTLiteClusters = LHCb::FTLiteCluster::FTLiteClusters;

  struct MatsCache {
    /**
     * partial SoA cache for mats, reserve enough (here 4096 which is more than enough)
     * space for all mats ( all mats should be less than 2 * 8 mats * 12 modules * 12 layers)
     */
    std::array<float, LHCb::Detector::FT::maxNumberMats> m_mats_x0;
    std::array<float, LHCb::Detector::FT::maxNumberMats> m_mats_dx;

    float m_dieGap;
    float m_sipmPitch;
    float m_uBegin;
    float m_halfChannelPitch;

    MatsCache() = default;

    MatsCache( const DeFT& ftDet ) {
      auto const first_mat =
          ftDet.firstMat(); // FIXME ftDet.stations()[0]->layers()[0]->quarters()[0]->modules()[0]->mats()[0];

// This parameters are constant accross all mats:
#ifdef USE_DD4HEP
      m_dieGap           = first_mat.dieGap();
      m_sipmPitch        = first_mat.sipmPitch();
      m_uBegin           = first_mat.uBegin();
      m_halfChannelPitch = first_mat.halfChannelPitch();
#else
      m_dieGap           = first_mat->dieGap();
      m_sipmPitch        = first_mat->sipmPitch();
      m_uBegin           = first_mat->uBegin();
      m_halfChannelPitch = first_mat->halfChannelPitch();
#endif
      // FIXME
      auto func = [this]( const DeFTMat& mat ) {
        assert( this->m_dieGap == mat.dieGap() && "Unexpected difference in dieGap" );
        assert( this->m_sipmPitch == mat.sipmPitch() && "Unexpected difference in sipmPitch" );
        assert( this->m_uBegin == mat.uBegin() && "Unexpected difference in uBegin" );
        assert( this->m_halfChannelPitch == mat.halfChannelPitch() && "Unexpected difference in halfChannelPitch" );
        auto index = mat.elementID().globalMatID();

        auto dxdy        = mat.dxdy();
        auto mirrorPoint = mat.mirrorPoint();
        auto ddx         = mat.ddx();

        this->m_mats_x0[index] = mirrorPoint.x() - mirrorPoint.y() * dxdy;
        this->m_mats_dx[index] = ddx.x() - ddx.y() * dxdy;
      };
      ftDet.applyToAllMats( func );
    }
  };
} // namespace

class SciFiTrackForwardingStoreHit
    : public LHCb::Algorithm::Transformer<SciFiTrackForwardingHits( EventContext const&, FTLiteClusters const&,
                                                                    MatsCache const& ),
                                          LHCb::DetDesc::usesConditions<MatsCache>> {
public:
  SciFiTrackForwardingStoreHit( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"HitsLocation", LHCb::FTLiteClusterLocation::Default},
                      // KeyValue{"MatsCache", "AlgorithmSpecific-" + name + "-MatsCache"}},
                      KeyValue{"MatsCache", "FT:" + name + "-MatsCache"}},
                     KeyValue{"Output", "Rec/SciFiTrackForwarding/Hits"} ) {}

  StatusCode initialize() override {
    return Transformer::initialize().andThen( [&] {
      addConditionDerivation<MatsCache( const DeFT& )>( {DeFTDetectorLocation::Default}, inputLocation<MatsCache>() );

      // TODO: this should be ~80 micron; get this from a tool?
      // std::array<float, 9> clusRes = {0.05f, 0.08f, 0.11f, 0.14f, 0.17f, 0.20f, 0.23f, 0.26f, 0.29f};
      // if ( m_updatedRes ) { clusRes = {0.5, 0.1, 0.1, 0.1, 0.080, 0.163, 0.38, 0.48, 0.58}; }
      // std::transform( clusRes.begin(), clusRes.end(), m_invClusResolution.begin(),
      //                []( const float& c ) { return 1.f / c; } );
    } );
  }

  SciFiTrackForwardingHits operator()( EventContext const&, FTLiteClusters const&, MatsCache const& ) const override;

  // new resolutions for testing
  // Gaudi::Property<bool> m_updatedRes{this, "UpdatedResolutions", false};
  /// Cached resolution
  // std::array<float, 9> m_invClusResolution;
};

DECLARE_COMPONENT( SciFiTrackForwardingStoreHit )

SciFiTrackForwardingHits SciFiTrackForwardingStoreHit::
                         operator()( EventContext const& evtCtx, FTLiteClusters const& clusters, MatsCache const& cache ) const {

  int                      size = clusters.size() + 2 * 24; // N clusters + 2 guards / zone
  SciFiTrackForwardingHits tmp{LHCb::getMemResource( evtCtx )};

  auto& hitvec = tmp.hits;
  auto& IDvec  = tmp.IDs;

  hitvec.reserve( size );
  IDvec.reserve( size );

  // TODO: Verify that the hits are sorted as expected
  /*assert( hitHandler.hits().is_sorted( []( const auto& lhs, const auto& rhs ) { return lhs.x() < rhs.x(); } ) &&
          "FT hits must be properly sorted for the pattern recognition "
          "Lower by X for each zone" );*/

  for ( unsigned i{0}; i < tmp.zonerange.size(); ++i ) {
    size = hitvec.size();
    hitvec.emplace_back( std::numeric_limits<float>::lowest() );
    IDvec.emplace_back( 0 );

    for ( int quarter = 0; quarter < 2; quarter++ ) {
      int iQuarter = i * 2 + quarter;
      for ( auto const& clus : clusters.range( iQuarter ) ) {
        LHCb::Detector::FTChannelID id    = clus.channelID();
        auto                        index = id.globalMatID();

        float uFromChannel = cache.m_uBegin + ( 2 * id.channel() + 1 + clus.fractionBit() ) * cache.m_halfChannelPitch;
        uFromChannel += id.die() * cache.m_dieGap;
        uFromChannel += id.sipm() * cache.m_sipmPitch;

        float x0 = cache.m_mats_x0[index] + cache.m_mats_dx[index] * uFromChannel;

        hitvec.emplace_back( x0 );
        IDvec.emplace_back( LHCb::LHCbID( id ).lhcbID() );
      }
    }

    // padding of one simd length
    for ( int idx{0}; idx < 8; ++idx ) {
      hitvec.emplace_back( std::numeric_limits<float>::max() );
      IDvec.emplace_back( 0 );
    }

    tmp.zonerange[i] = {size, hitvec.size() - size};
  }

  return tmp;
}
