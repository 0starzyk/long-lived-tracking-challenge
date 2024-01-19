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
// Include files
#include "DetDesc/IConditionDerivationMgr.h"
#include "Event/ODIN.h"
#include "Event/PrLongTracks.h"
#include "Gaudi/Accumulators.h"
#include "GaudiKernel/IRegistry.h"
#include "LHCbAlgs/Transformer.h"
#include "PrKernel/PrHit.h"
#include "PrKernel/UTHit.h"
#include "PrKernel/UTHitHandler.h"
#include "UTDAQ/UTInfo.h"
#include "UTDet/DeUTDetector.h"
#include <Vc/Vc>
#include <vector>

#include "boost/container/small_vector.hpp"
#include "boost/container/static_vector.hpp"
#include <memory>

namespace LHCb::Pr {

  //-----------------------------------------------------------------------------
  // class : PrResidualUTHits
  // Store residual UTHits after other Algorithms, e.g. PrMatchNN/PrForward used
  //
  // 2020-04-21 : Peilian Li
  //
  //-----------------------------------------------------------------------------

  class PrResidualUTHits
      : public Algorithm::Transformer<UT::HitHandler( Long::Tracks const&, UT::HitHandler const&, DeUTDetector const& ),
                                      DetDesc::usesConditions<DeUTDetector>> {
  public:
    PrResidualUTHits( const std::string& name, ISvcLocator* pSvcLocator );
    UT::HitHandler operator()( Long::Tracks const&, UT::HitHandler const&, DeUTDetector const& ) const override;
  };

  // Declaration of the Algorithm Factory
  DECLARE_COMPONENT_WITH_ID( LHCb::Pr::PrResidualUTHits, "PrResidualUTHits" )

} // namespace LHCb::Pr

LHCb::Pr::PrResidualUTHits::PrResidualUTHits( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator,
                   {{"TracksLocation", ""}, {"UTHitsLocation", ""}, {"DeUTLocation", DeUTDetLocation::location()}},
                   {"UTHitsOutput", ""} ) {}

UT::HitHandler LHCb::Pr::PrResidualUTHits::operator()( Long::Tracks const& tracks, UT::HitHandler const& uthithandler,
                                                       DeUTDetector const& utDet ) const {
  UT::HitHandler tmp{};

  if ( tracks.size() == 0 ) {
    if ( msgLevel( MSG::DEBUG ) )
      debug() << "Track container '" << inputLocation<Long::Tracks>() << "' is empty" << endmsg;
    return uthithandler;
  }

  std::vector<long unsigned int> usedUTHits{};
  usedUTHits.reserve( uthithandler.nbHits() );

  for ( const auto& track : tracks.scalar() ) {
    const auto ids = track.lhcbIDs();
    for ( auto id : ids ) {
      if ( !( id.isUT() ) ) continue;
      usedUTHits.emplace_back( id.utID().channelID() );
    }
  }

  for ( int iStation = 1; iStation <= static_cast<int>( UTInfo::DetectorNumbers::Stations ); ++iStation ) {
    for ( int iLayer = 1; iLayer <= static_cast<int>( UTInfo::DetectorNumbers::Layers ); ++iLayer ) {
      for ( int iRegion = 1; iRegion <= static_cast<int>( UTInfo::DetectorNumbers::Regions ); ++iRegion ) {
        for ( int iSector = 1; iSector <= static_cast<int>( UTInfo::DetectorNumbers::Sectors ); ++iSector ) {
          for ( auto& uthit : uthithandler.hits( iStation, iLayer, iRegion, iSector ) ) {
            bool used = std::any_of( usedUTHits.begin(), usedUTHits.end(),
                                     [utid = uthit.chanID().channelID()]( const auto& id ) { return utid == id; } );

            if ( used ) continue;
            const unsigned int fullChanIdx = UT::HitHandler::HitsInUT::idx( iStation, iLayer, iRegion, iSector );
            auto&              aSector     = utDet.getSector( uthit.chanID() );
            tmp.emplace_back( aSector, fullChanIdx, uthit.strip(), uthit.fracStrip(), uthit.chanID(), uthit.size(),
                              uthit.highThreshold() );
          }
        }
      }
    }
  }
  return tmp;
}
