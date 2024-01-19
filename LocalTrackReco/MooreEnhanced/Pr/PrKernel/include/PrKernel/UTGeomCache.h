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
#include <Kernel/IUTReadoutTool.h>
#include <Kernel/STLExtensions.h>
#include <PrKernel/UTHitHandler.h>
#include <UTDAQ/UTInfo.h>
#include <UTDet/DeUTDetector.h>
#include <cassert>
#include <cstdint>
#include <vector>

namespace LHCb::Pr::UT {

  struct UTGeomCache {
    static constexpr int NSectorPerBoard = 6;
    static constexpr int NSectorAllBoard = 240 * NSectorPerBoard;
    /// Faster access to sectors
#ifdef USE_DD4HEP
    std::array<DeUTSector, NSectorAllBoard> sectors{};
#else
    std::array<DeUTSector const*, NSectorAllBoard> sectors{};
#endif

    struct FullChan {
      unsigned int idx{0};
      unsigned int chanID{0};
    };
    std::array<FullChan, NSectorAllBoard> fullchan;

    UTGeomCache() = default;
    UTGeomCache( const DeUTDetector& utDet, const IUTReadoutTool& ro, const IUTReadoutTool::ReadoutInfo& roInfo ) {
      sectors.fill( nullptr );
      for ( unsigned int srcId = 0; srcId < ro.nBoard( &roInfo ); srcId++ ) {
        for ( const auto& [idx, sector] :
              range::enumerate( ro.findBySourceID( srcId, &roInfo ),
                                ro.TELLNumberToSourceID( srcId + 1, &roInfo ) * UTGeomCache::NSectorPerBoard ) ) {
          assert( idx < NSectorAllBoard );
#ifdef USE_DD4HEP
          if ( sectors[idx].isValid() ) throw std::runtime_error( "UTGeomCache: duplicate sector???" );
          sectors[idx] = utDet.getSector( sector );
#else
          if ( sectors[idx] ) throw std::runtime_error( "UTGeomCache: duplicate sector???" );
          sectors[idx] = &utDet.getSector( sector );
#endif
          fullchan[idx] = {static_cast<unsigned int>( ::UT::HitHandler::HitsInUT::idx(
                               sector.station(), sector.layer(), sector.detRegion(), sector.sector() ) ),
                           static_cast<unsigned int>( sector )};
        }
      }
    }
  };
} // namespace LHCb::Pr::UT
