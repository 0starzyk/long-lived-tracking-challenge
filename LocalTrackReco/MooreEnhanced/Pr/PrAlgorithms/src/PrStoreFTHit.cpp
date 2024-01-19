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
#include "Detector/FT/FTConstants.h"
#include "Event/FTLiteCluster.h"
#include "FTDet/DeFTDetector.h"
#include "LHCbAlgs/Transformer.h"
#include "PrKernel/FTMatsCache.h"
#include "PrKernel/PrFTHitHandler.h"
#include "range/v3/view/transform.hpp"
#include <algorithm>
#include <array>
#include <string>

//-----------------------------------------------------------------------------
// Implementation file for class : PrStoreFTHit
// This algorithms is expected to store the hits in TES together with Geometry
// information ( should be fixed ?) . Geometry information are encoded in the
// PrFTHiHandler::m_zones private variable
//
// 2016-07-07 : Renato Quagliani
//-----------------------------------------------------------------------------
namespace {
  using FTLiteClusters = LHCb::FTLiteCluster::FTLiteClusters;
  using MatsCache      = FTMatsCache::MatsCache;
} // namespace

class PrStoreFTHit
    : public LHCb::Algorithm::Transformer<PrFTHitHandler<PrHit>( const FTLiteClusters&, const MatsCache& ),
                                          LHCb::DetDesc::usesConditions<MatsCache>> {

public:
  //=============================================================================
  // Standard constructor, initializes variables
  //=============================================================================
  PrStoreFTHit( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"InputLocation", LHCb::FTLiteClusterLocation::Default},
                      KeyValue{"FTMatsCache", FTMatsCache::Location + name}},
                     KeyValue{"FTHitsLocation", PrFTInfo::FTHitsLocation} ) {}

  StatusCode initialize() override {
    return Transformer::initialize().andThen( [&] {
      addConditionDerivation<MatsCache( const DeFT& )>( {DeFTDetectorLocation::Default}, inputLocation<MatsCache>() );
      std::array<float, 9> clusRes = {0.05, 0.08, 0.11, 0.14, 0.17, 0.20, 0.23, 0.26, 0.29};
      if ( m_updatedRes ) { clusRes = {0.5, 0.1, 0.1, 0.1, 0.080, 0.163, 0.38, 0.48, 0.58}; }
      std::transform( clusRes.begin(), clusRes.end(), m_invClusResolution.begin(), []( auto c ) { return 1.f / c; } );
    } );
  }

  PrFTHitHandler<PrHit> operator()( const FTLiteClusters& clusters, const MatsCache& cache ) const override;

private:
  // Calib
  Gaudi::Property<bool> m_applyMatContractionCalibration{this, "ApplyMatContractionCalibration",
                                                         false}; // TODO: change to true when calibration is available
  // Type of resolution
  Gaudi::Property<bool> m_updatedRes{this, "UpdatedResolutions", false};

  // Cached resolution
  std::array<float, 9> m_invClusResolution;
};
DECLARE_COMPONENT( PrStoreFTHit )
//=============================================================================
// Main execution
//=============================================================================
PrFTHitHandler<PrHit> PrStoreFTHit::operator()( const FTLiteClusters& clusters, const MatsCache& cache ) const {
  // create a hitHandler to be returned and stored in the TES
  PrFTHitHandler<PrHit> hitHandler( clusters.size() );
  if ( msgLevel( MSG::DEBUG ) ) debug() << "Retrieved " << clusters.size() << " clusters" << endmsg;

  for ( uint16_t iQuarter = 0; iQuarter < LHCb::Detector::FT::nQuartersTotal; ++iQuarter ) {
    uint info = ( iQuarter >> 1 ) | ( ( ( iQuarter << 4 ) ^ ( iQuarter << 5 ) ^ 128u ) & 128u ); // FIXME

    auto r = ranges::views::transform( clusters.range( iQuarter ), [&]( LHCb::FTLiteCluster clus ) -> PrHit {
      const auto id           = clus.channelID();
      const auto index        = id.globalMatID();
      const auto dxdy         = cache.dxdy[index];
      const auto dzdy         = cache.dzdy[index];
      const auto globaldy     = cache.globaldy[index];
      auto       uFromChannel = cache.uBegin + ( 2 * id.channel() + 1 + clus.fractionBit() ) * cache.halfChannelPitch;
      if ( id.die() ) uFromChannel += cache.dieGap;
      uFromChannel += id.sipm() * cache.sipmPitch;
      const auto endPoint = cache.mirrorPoint[index] + cache.ddx[index] * uFromChannel;
      const auto x0_orig  = endPoint.x() - dxdy * endPoint.y();
      auto       x0       = x0_orig; // to be updated by conditions if available

      const auto matContractionParameterVector = cache.matContractionParameterVector[index];
      const bool matContractionAvailable       = matContractionParameterVector.size() != 0;
      if ( m_applyMatContractionCalibration ) {
        if ( !matContractionAvailable ) {
          warning() << "Requested to correct FT mats for temperature distortions but conditions not available. These "
                       "corrections cannot be applied. Check your conditions tags if you intend to apply these "
                       "corrections or set property ApplyMatContractionCalibration to false if you do not. Continuing "
                       "without applying corrections."
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

      const auto z0   = endPoint.z() - dzdy * endPoint.y();
      auto       yMin = endPoint.y();
      auto       yMax = yMin + globaldy;
      if ( id.isBottom() ) std::swap( yMin, yMax );
      assert( clus.pseudoSize() < 9 && "Pseudosize of cluster is > 8. Out of range." );
      const auto werrX = m_invClusResolution[clus.pseudoSize()];
      return {LHCb::LHCbID( id ), x0, z0, dxdy, dzdy, yMin, yMax, werrX, werrX * werrX, info};
    } );

    hitHandler.insert( ( iQuarter >> 1 ), r.begin(), r.end() );
  }

  // Verify that the hits are sorted as expected
  assert( hitHandler.hits().is_sorted( []( const auto& lhs, const auto& rhs ) { return lhs.x() < rhs.x(); } ) &&
          "FT hits must be properly sorted for the pattern recognition "
          "Lower by X for each zone" );

  if ( msgLevel( MSG::VERBOSE ) ) {
    const auto& container = hitHandler.hits();
    const auto& offsets   = container.offsets();
    for ( size_t i = 0; i < container.size(); ++i ) {
      unsigned id = 2 * container.hit( i ).planeCode() + container.hit( i ).zone();
      verbose() << std::setw( 6 ) << std::right << i << " [" << offsets[id].first << ";" << offsets[id].second << "] "
                << std::setw( 6 ) << std::right << id << std::setw( 10 ) << std::right << container.hit( i ).x()
                << endmsg;
    }
  }
  return hitHandler;
}
