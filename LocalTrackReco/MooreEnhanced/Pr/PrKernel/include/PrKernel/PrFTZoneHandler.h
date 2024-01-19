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
#pragma once

// Include files
#include "FTDAQ/FTInfo.h"
#include "FTDet/DeFTDetector.h"
#include "Kernel/DetectorSegment.h"
#include "PrKernel/PrHitZone.h"
#include <cassert>
#include <stdexcept>

namespace FTZoneCache {

  const std::string Location    = "AlgorithmSpecific-FTZoneCache";
  const std::string MatLocation = "AlgorithmSpecific-FTMatsCache";

  /** @class PrFTZoneHandler PrFTZoneHandler.h
   * Handlers of zones, the object is stored in the detector event store as a condition
   * and each algorithms reads the object from there without calling the HitManagers (tools)
   * @author Renato Quagliani
   * @author Sebastien Ponce
   */

  class PrFTZoneHandler final {

  public:
    PrFTZoneHandler( DeFT const& ftDet ) {
#ifdef USE_DD4HEP
      auto func = [this]( const DeFTLayer& layer ) {
        const auto id = layer.layerIdx();
        // fixme
        DetectorSegment seg( 0, layer.globalZ(), layer.dxdy(), layer.dzdy(), 0., 0. );
        const auto      xmax = 0.5f * layer.sizeX();
        const auto      ymax = 0.5f * layer.sizeY();
        // The setGeometry defines the z at y=0, the dxDy and the dzDy, as well as the isX properties of the zone.
        // This is important, since these are used in the following.
        // They are set once for each zone in this method.
        const bool ok =
            ( this->MakeZone( 2 * id + 1, seg, -xmax, xmax, -25.f, ymax ) && // Small overlap (25 mm) for stereo layers
              this->MakeZone( 2 * id, seg, -xmax, xmax, -ymax, 25.f ) );     // Small overlap (25 mm) for stereo layers
        if ( !ok ) { throw std::runtime_error( "Failed to create DeFT Zones for ID = " + std::to_string( id ) ); }
      };
      ftDet.applyToAllLayers( func );
#else
      for ( auto station : ftDet.stations() ) {
        for ( auto layer : station->layers() ) {
          const auto      id = 4 * ( station->stationID() - 1 ) + layer->layerID();
          DetectorSegment seg( 0, layer->globalZ(), layer->dxdy(), layer->dzdy(), 0., 0. );
          const auto      xmax = 0.5f * layer->sizeX();
          const auto      ymax = 0.5f * layer->sizeY();
          // The setGeometry defines the z at y=0, the dxDy and the dzDy, as well as the isX properties of the zone.
          // This is important, since these are used in the following.
          // They are set once for each zone in this method.
          const bool ok =
              ( this->MakeZone( 2 * id + 1, seg, -xmax, xmax, -25.f, ymax ) && // Small overlap (25 mm) for stereo
                                                                               // layers
                this->MakeZone( 2 * id, seg, -xmax, xmax, -ymax, 25.f ) ); // Small overlap (25 mm) for stereo layers
          if ( !ok ) { throw std::runtime_error( "Failed to create DeFT Zones for ID = " + std::to_string( id ) ); }
        }
      }
#endif
    }
    /// Standard constructor
    PrFTZoneHandler() = default;

    bool MakeZone( unsigned int n, DetectorSegment& seg, float xMin, float xMax, float yMin, float yMax ) {
      if ( n < m_zones.size() ) {
        m_zones[n].setZone( n, seg, xMin, xMax, yMin, yMax );
        return true;
      } else {
        return false;
      }
    }
    const PrHitZone& zone( unsigned int n ) const {
      if ( n >= m_zones.size() ) {
        throw std::runtime_error( "Zone index " + std::to_string( n ) + " is out-of-range" );
      }
      return m_zones[n];
    }

    template <PrHitZone::Side SIDE>
    static int getXZone( int layer ) {
      if constexpr ( SIDE == PrHitZone::Side::Upper ) {
        return LHCb::Detector::FT::xZonesUpper[layer];
      } else {
        return LHCb::Detector::FT::xZonesLower[layer];
      }
    }

    template <PrHitZone::Side SIDE>
    static int getUVZone( int layer ) {
      if constexpr ( SIDE == PrHitZone::Side::Upper ) {
        return LHCb::Detector::FT::uvZonesUpper[layer];
      } else {
        return LHCb::Detector::FT::uvZonesLower[layer];
      }
    }

    template <PrHitZone::Side SIDE>
    static int getTriangleZone( int layer ) {
      if constexpr ( SIDE == PrHitZone::Side::Upper ) {
        return LHCb::Detector::FT::uvZonesLower[layer];
      } else {
        return LHCb::Detector::FT::uvZonesUpper[layer];
      }
    }

  private:
    // plain vector with indexing needed!
    std::array<PrHitZone, LHCb::Detector::FT::Numbers::NFTZones> m_zones;
  };

  /**
   * @class ZoneCache PrFTZoneHandler.h
   * @brief Caches derived conditions, e.g. positions of SciFi layers
   * @author André Günther
   */
  struct ZoneCache final {
    PrFTZoneHandler                                 handler;
    PrHitZone                                       lowerLastZone{}, upperLastZone{};
    std::array<float, LHCb::Detector::FT::NFTZones> zoneDxDy{std::numeric_limits<float>::signaling_NaN()};
    std::array<float, LHCb::Detector::FT::NFTZones> zoneDzDy{std::numeric_limits<float>::signaling_NaN()};
    std::array<float, LHCb::Detector::FT::NFTZones> zoneZPos{std::numeric_limits<float>::signaling_NaN()};
    ZoneCache( const DeFT& ftDet ) : handler( ftDet ) {
      for ( int i{0}; i < LHCb::Detector::FT::NFTLayers; ++i ) {
        zoneZPos[2 * i]     = handler.zone( 2 * i ).z();
        zoneZPos[2 * i + 1] = handler.zone( 2 * i + 1 ).z();
        zoneDxDy[2 * i]     = handler.zone( 2 * i ).dxDy();
        zoneDxDy[2 * i + 1] = handler.zone( 2 * i + 1 ).dxDy();
        zoneDzDy[2 * i]     = handler.zone( 2 * i ).dzDy();
        zoneDzDy[2 * i + 1] = handler.zone( 2 * i + 1 ).dzDy();
      }
      lowerLastZone = handler.zone( LHCb::Detector::FT::NFTZones - 2 );
      upperLastZone = handler.zone( LHCb::Detector::FT::NFTZones - 1 );
    }
    // for debugging
    friend std::ostream& operator<<( std::ostream& os, const ZoneCache& c ) {
      os << "ZoneCache:" << std::endl;
      os << "zoneDxDy = [ ";
      for ( auto z : c.zoneDxDy ) { os << z << " "; }
      os << "]" << std::endl;
      os << "zoneDzDy = [ ";
      for ( auto z : c.zoneDzDy ) { os << z << " "; }
      os << "]" << std::endl;
      os << "zoneZPos = [ ";
      for ( auto z : c.zoneZPos ) { os << z << " "; }
      os << "]" << std::endl;
      os << c.lowerLastZone << c.upperLastZone;
      return os;
    }
  };

} // namespace FTZoneCache
