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
#include "Event/PrHits.h"
#include "Gaudi/Accumulators/Histogram.h"
#include "LHCbAlgs/Consumer.h"
#include <array>
#include <string>

namespace LHCb::Pr::FT {
  using namespace LHCb::Detector::FT;
  class SciFiHitsMonitor : public Algorithm::Consumer<void( const Hits& )> {
  public:
    SciFiHitsMonitor( const std::string& name, ISvcLocator* pSvcLocator )
        : Consumer( name, pSvcLocator, {KeyValue{"SciFiHits", ""}} ) {}

    void operator()( const Hits& hits ) const override {
      for ( auto zone{0u}; zone < NFTZones; ++zone ) {
        const auto [start, end] = hits.getZoneIndices( zone );
        auto xPosHistoBuffer    = m_hitsXPosHistos.at( zone / 2u ).buffer();
        for ( auto iHit = start; iHit < end; ++iHit ) { ++xPosHistoBuffer[hits.x( iHit )]; }
      }
    };

  private:
    using Histo = Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, float>;
    mutable std::array<Histo, NFTLayers> m_hitsXPosHistos{{
        {this, "T1X1", "SciFi hits x positions T1X1", {2000, -3000, 3000}},
        {this, "T1U", "SciFi hits x positions T1U", {2000, -3000, 3000}},
        {this, "T1V", "SciFi hits x positions T1V", {2000, -3000, 3000}},
        {this, "T1X2", "SciFi hits x positions T1X2", {2000, -3000, 3000}},
        {this, "T2X1", "SciFi hits x positions T2X1", {2000, -3000, 3000}},
        {this, "T2U", "SciFi hits x positions T2U", {2000, -3000, 3000}},
        {this, "T2V", "SciFi hits x positions T2V", {2000, -3000, 3000}},
        {this, "T2X2", "SciFi hits x positions T2X2", {2000, -3000, 3000}},
        {this, "T3X1", "SciFi hits x positions T3X1", {2000, -3200, 3200}},
        {this, "T3U", "SciFi hits x positions T3U", {2000, -3200, 3200}},
        {this, "T3V", "SciFi hits x positions T3V", {2000, -3200, 3200}},
        {this, "T3X2", "SciFi hits x positions T3X2", {2000, -3200, 3200}},
    }};
  };

  DECLARE_COMPONENT_WITH_ID( SciFiHitsMonitor, "PrSciFiHitsMonitor" )

} // namespace LHCb::Pr::FT
