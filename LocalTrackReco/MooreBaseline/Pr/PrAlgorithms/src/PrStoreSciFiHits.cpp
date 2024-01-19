
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
#include "Detector/FT/FTChannelID.h"
#include "Event/FTLiteCluster.h"
#include "Event/PrHits.h"
#include "FTDAQ/FTInfo.h"
#include "FTDet/DeFTDetector.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Transformer.h"
#include "LHCbMath/bit_cast.h"
#include "PrKernel/FTMatsCache.h"
#include "PrKernel/PrFTZoneHandler.h"

#include <array>
#include <limits>
#include <memory>
#include <string>

#include "GaudiKernel/StdArrayAsProperty.h"

/** @class StoreHits PrStoreSciFiHits.cpp
 *
 *  @brief Transforms FTLiteClusters into the input format needed by the PrForwardTracking
 */

namespace LHCb::Pr::FT {
  using FTLiteClusters = FTLiteCluster::FTLiteClusters;
  using MatsCache      = FTMatsCache::MatsCache;

  // TODO: get this from a tool?
  constexpr auto invClusRes2 = [] {
    auto tmp = std::array{0.05f, 0.08f, 0.11f, 0.14f, 0.17f, 0.20f, 0.23f, 0.26f, 0.29f};
    for ( std::size_t i{0}; i < tmp.size(); ++i ) {
      const auto inv = 1.f / tmp[i];
      tmp[i]         = inv * inv;
    }
    return tmp;
  }();

  class StoreHits : public Algorithm::Transformer<Hits( const EventContext&, const FTLiteClusters&, const MatsCache& ),
                                                  Algorithm::Traits::usesConditions<MatsCache>> {
  public:
    StoreHits( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer( name, pSvcLocator,
                       {KeyValue{"HitsLocation", FTLiteClusterLocation::Default},
                        KeyValue{"FTMatsCache", FTZoneCache::MatLocation + name}},
                       KeyValue{"Output", PrFTInfo::SciFiHitsLocation} ) {}

    StatusCode initialize() override {
      return Transformer::initialize().andThen( [&] {
        addConditionDerivation<MatsCache( const DeFT& )>( {DeFTDetectorLocation::Default}, inputLocation<MatsCache>() );
      } );
    }

    Hits operator()( const EventContext&, const FTLiteClusters&, const MatsCache& ) const override;

  private:
    // Calib
    Gaudi::Property<bool> m_applyMatContractionCalibration{this, "ApplyMatContractionCalibration",
                                                           false}; // TODO: change to true when calibration is available
    Gaudi::Property<std::array<bool, FTConstants::nLayersTotal>> m_layerMasks{this, "LayerMasks", {}};
    // Counters
    using SC = Gaudi::Accumulators::StatCounter<>;
    using SCbuf =
        Gaudi::Accumulators::Buffer<Gaudi::Accumulators::StatAccumulator, Gaudi::Accumulators::atomicity::full, double>;
    mutable SC                                               m_cntTotalHits{this, "Total number of hits"};
    mutable std::array<SC, LHCb::Detector::FT::nLayersTotal> m_cntHitsPerLayer{
        SC{this, "Hits in T1X1"}, SC{this, "Hits in T1U"}, SC{this, "Hits in T1V"}, SC{this, "Hits in T1X2"},
        SC{this, "Hits in T2X1"}, SC{this, "Hits in T2U"}, SC{this, "Hits in T2V"}, SC{this, "Hits in T2X2"},
        SC{this, "Hits in T3X1"}, SC{this, "Hits in T3U"}, SC{this, "Hits in T3V"}, SC{this, "Hits in T3X2"}};
    mutable std::array<SC, LHCb::Detector::FT::nLayersTotal> m_cntXPerLayer{
        SC{this, "Average X in T1X1"}, SC{this, "Average X in T1U"},  SC{this, "Average X in T1V"},
        SC{this, "Average X in T1X2"}, SC{this, "Average X in T2X1"}, SC{this, "Average X in T2U"},
        SC{this, "Average X in T2V"},  SC{this, "Average X in T2X2"}, SC{this, "Average X in T3X1"},
        SC{this, "Average X in T3U"},  SC{this, "Average X in T3V"},  SC{this, "Average X in T3X2"}};
  };

  DECLARE_COMPONENT_WITH_ID( StoreHits, "PrStoreSciFiHits" )

  Hits StoreHits::operator()( const EventContext& evtCtx, FTLiteClusters const& clusters,
                              MatsCache const& cache ) const {

    Hits hits{LHCb::getMemResource( evtCtx )};

    constexpr auto xu  = LHCb::Detector::FT::xZonesUpper;
    constexpr auto uvu = LHCb::Detector::FT::uvZonesUpper;

    constexpr auto xd  = LHCb::Detector::FT::xZonesLower;
    constexpr auto uvd = LHCb::Detector::FT::uvZonesLower;

    // hits are stored in the same order as the layers are in z, i.e. x-u-v-x
    // hit zone order here: 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23
    constexpr auto hitzones = std::array<int, LHCb::Detector::FT::NFTZones>{
        xd[0], uvd[0], uvd[1], xd[1], xd[2], uvd[2], uvd[3], xd[3], xd[4], uvd[4], uvd[5], xd[5],
        xu[0], uvu[0], uvu[1], xu[1], xu[2], uvu[2], uvu[3], xu[3], xu[4], uvu[4], uvu[5], xu[5]};

    auto                                                bufTotalHits    = m_cntTotalHits.buffer();
    std::array<SCbuf, LHCb::Detector::FT::nLayersTotal> bufHitsPerLayer = {
        m_cntHitsPerLayer[0].buffer(), m_cntHitsPerLayer[1].buffer(),  m_cntHitsPerLayer[2].buffer(),
        m_cntHitsPerLayer[3].buffer(), m_cntHitsPerLayer[4].buffer(),  m_cntHitsPerLayer[5].buffer(),
        m_cntHitsPerLayer[6].buffer(), m_cntHitsPerLayer[7].buffer(),  m_cntHitsPerLayer[8].buffer(),
        m_cntHitsPerLayer[9].buffer(), m_cntHitsPerLayer[10].buffer(), m_cntHitsPerLayer[11].buffer()};
    std::array<SCbuf, LHCb::Detector::FT::nLayersTotal> bufXPerLayer = {
        m_cntXPerLayer[0].buffer(), m_cntXPerLayer[1].buffer(),  m_cntXPerLayer[2].buffer(),
        m_cntXPerLayer[3].buffer(), m_cntXPerLayer[4].buffer(),  m_cntXPerLayer[5].buffer(),
        m_cntXPerLayer[6].buffer(), m_cntXPerLayer[7].buffer(),  m_cntXPerLayer[8].buffer(),
        m_cntXPerLayer[9].buffer(), m_cntXPerLayer[10].buffer(), m_cntXPerLayer[11].buffer()};
    bufTotalHits += clusters.size();

    for ( auto i : hitzones ) {

      hits.setZoneIndex( i, hits.size() );

      if ( m_layerMasks[(unsigned int)( i / 2 )] ) {
        hits.appendColumn( std::numeric_limits<uint8_t>::max(), 1.e9f, {}, {} );
        continue;
      }

      for ( auto quarter{0}; quarter < 2; ++quarter ) {
        const auto iQuarter = bit_cast<unsigned>( i * 2 + quarter );
        const auto info =
            static_cast<uint8_t>( ( iQuarter >> 1 ) | ( ( ( iQuarter << 4 ) ^ ( iQuarter << 5 ) ^ 128u ) & 128u ) );
        bufHitsPerLayer[iQuarter / 4] += clusters.range( iQuarter ).size();
        for ( const auto& clus : clusters.range( iQuarter ) ) {
          const auto id = clus.channelID();

          const auto hfChPitch    = ( 2 * id.channel() + 1 + clus.fractionBit() ) * cache.halfChannelPitch;
          const auto dieGap       = id.die() * cache.dieGap;
          const auto sipmPitch    = id.sipm() * cache.sipmPitch;
          const auto uFromChannel = cache.uBegin + hfChPitch + dieGap + sipmPitch;
          const auto index        = id.globalMatID();
          const auto endPoint     = cache.mirrorPoint[index] + cache.ddx[index] * uFromChannel;

          const auto dxdy    = cache.dxdy[index];
          const auto dzdy    = cache.dzdy[index];
          auto       yMin    = endPoint.y();
          const auto x0_orig = endPoint.x() - dxdy * yMin;
          auto       x0      = x0_orig; // to be updated by conditions if available
          // const auto x0   = endPoint.x() - dxdy * yMin;
          bufXPerLayer[iQuarter / 4] += x0_orig;
          const auto z0   = endPoint.z() - dzdy * yMin;
          auto       yMax = yMin + cache.globaldy[index];
          if ( id.isBottom() ) std::swap( yMin, yMax );
          assert( clus.pseudoSize() < 9 && "Pseudosize of cluster is > 8. Out of range." );

          const auto matContractionParameterVector = cache.matContractionParameterVector[index];
          const bool matContractionAvailable       = matContractionParameterVector.size() != 0;
          if ( m_applyMatContractionCalibration ) {
            if ( !matContractionAvailable ) {
              warning() << "Requested to correct FT mats for temperature distortions but conditions not available. "
                           "These corrections cannot be applied. Check your conditions tags if you intend to apply "
                           "these corrections or set property ApplyMatContractionCalibration to false if you do not. "
                           "Continuing without applying corrections."
                        << endmsg;
            } else {
              const auto matContraction = matContractionParameterVector.at( ( id.sipm() * 128 ) + id.channel() );
              const auto x0Calibration =
                  ( cache.ddx[index] * matContraction ).x() - dxdy * ( cache.ddx[index] * matContraction ).y();
              const auto calibratedx0 = x0_orig + x0Calibration;
              x0                      = calibratedx0;
              debug() << "x0 before contraction correction: " << x0_orig
                      << ", x0 after contraction correction: " << calibratedx0 << endmsg;
            }
          }
          hits.appendColumn( ( info & 63u ) >> 1, x0, {invClusRes2[clus.pseudoSize()], dzdy, dxdy, z0},
                             {id, yMin, yMax} );
        }
      }
      // add a large number at the end of x hits for each zone to stop binary search before end
      hits.appendColumn( std::numeric_limits<uint8_t>::max(), 1.e9f, {}, {} );
    }
    /** the start index is given by zonesIndexes[zone] but the end index is given by
     * zoneIndexes[zone+2] because upper and lower zones are adjacent in the container
     * so when asking for start and end of zone 22, 24 gives the beginning of the first
     * upper zone which is the end of zone 22
     */
    hits.setZoneIndex( LHCb::Detector::FT::NFTZones, hits.getZoneIndex( xu[0] ) );
    // so when asking for lastZone, lastZone+2 gives the very end of the container
    hits.setZoneIndex( LHCb::Detector::FT::NFTZones + 1, hits.size() );
    // avoid FPEs
    for ( unsigned i{0}; i < SIMDWrapper::best::types::size; ++i ) {
      hits.appendColumn( std::numeric_limits<uint8_t>::max(), 1.e9f, {}, {} );
    }
    return hits;
  }
} // namespace LHCb::Pr::FT
