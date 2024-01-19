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
#include "Event/PrHits.h"
#include "Event/RawBank.h"
#include "Kernel/IUTReadoutTool.h"
#include "Kernel/UTDAQBoard.h"
#include "Kernel/UTDecoder.h"
#include "Kernel/UTTell1Board.h"
#include "LHCbAlgs/Transformer.h"
#include "PrKernel/UTGeomCache.h"
#include "PrKernel/UTHitHandler.h"
#include "UTDet/DeUTDetector.h"
#include "fmt/format.h"
#include <cassert>
#include <memory>
#include <numeric>

namespace LHCb::Pr::UT {

  template <typename HANDLER>
  using Transformer =
      LHCb::Algorithm::Transformer<HANDLER( const EventContext&, const RawBank::View&, const UTGeomCache& ),
                                   LHCb::DetDesc::usesConditions<UTGeomCache>>;

  template <typename HANDLER>
  class StoreHit : public Transformer<HANDLER> {
  public:
    using KeyValue = typename Transformer<HANDLER>::KeyValue;
    using Transformer<HANDLER>::inputLocation;

    StoreHit( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer<HANDLER>( name, pSvcLocator,
                                {KeyValue{"RawBanks", "DAQ/RawBanks/UT"},
                                 KeyValue{"GeomCache", "AlgorithmSpecific-" + name + "-UTGeomCache"}},
                                KeyValue{"UTHitsLocation", UTInfo::HitLocation} ) {}

    StatusCode initialize() override {
      return Transformer<HANDLER>::initialize().andThen( [&] {
        // TODO : alignment need the updateSvc for detector ( UT experts needed )
        this->addConditionDerivation( {DeUTDetLocation::location(), m_readoutTool->getReadoutInfoKey()},
                                      this->template inputLocation<UTGeomCache>(),
                                      [this]( const DeUTDetector& utDet, IUTReadoutTool::ReadoutInfo const& roInfo ) {
                                        return UTGeomCache{utDet, *m_readoutTool, roInfo};
                                      } );
      } );
    }

    HANDLER operator()( const EventContext& evtCtx, const LHCb::RawBank::View& tBanks,
                        const UTGeomCache& cache ) const override {
      HANDLER hitHandler{Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx )};
      hitHandler.reserve( 10000 );

      try {

        for ( const auto& bank : tBanks ) {
          // make local decoder
          if ( bank->size() == 0 ) continue;

          auto decode = [&hitHandler, geomOffset = bank->sourceID() * UTGeomCache::NSectorPerBoard,
                         &cache]( auto decoder_range ) {
            for ( const auto& aWord : decoder_range ) {
              const std::size_t geomIdx = geomOffset + ( aWord.channelID() / 512 );
              assert( geomIdx < cache.sectors.size() );
              assert( geomIdx < cache.fullchan.size() );

              auto        aSector  = cache.sectors[geomIdx];
              const auto& fullChan = cache.fullchan[geomIdx];

              // FIXME: move functionality into channelID -- why is this _different_ from UTChennelID::strip() ??
              const auto strip = ( aWord.channelID() & 0x1ff ) + 1;

#ifdef USE_DD4HEP
              hitHandler.emplace_back( aSector, fullChan.idx, strip, aWord.fracStripBits(),
#else
              hitHandler.emplace_back( *aSector, fullChan.idx, strip, aWord.fracStripBits(),
#endif
                                       Detector::UT::ChannelID{fullChan.chanID + strip}, aWord.pseudoSizeBits(),
                                       aWord.hasHighThreshold() );
            }
          };
          switch ( UTDAQ::version{bank->version()} ) {
          case UTDAQ::version::v5:
            if ( !isCluster )
              decode( UTDecoder<UTDAQ::version::v5>{*bank}.posRange() );
            else
              decode( UTDecoder<UTDAQ::version::v5>{*bank}.posAdcRange( isAdcMax, stripMax ) );
            break;
          case UTDAQ::version::v4:
            decode( UTDecoder<UTDAQ::version::v4>{*bank}.posRange() );
            break;
          default:
            throw std::runtime_error{"unknown version of the RawBank"}; /* OOPS: unknown format */
          };
        }

        m_nBanks += tBanks.size();

        if constexpr ( std::is_same_v<HANDLER, ::LHCb::Pr::UT::Hits> ) hitHandler.addPadding();

      } catch ( std::runtime_error& ) { // FIXME: temporary work around for MC produced in dec 2020 - may 2021 -- for
                                        // FEST only!
        hitHandler.clear();
        ++m_bad_data;
      }

      return hitHandler;
    }

  private:
    //---Properties
    Gaudi::Property<bool>                               isCluster{this, "isCluster", true};
    Gaudi::Property<bool>                               isAdcMax{this, "isAdcMax", false};
    Gaudi::Property<unsigned int>                       stripMax{this, "stripMax", 4};
    mutable Gaudi::Accumulators::SummingCounter<>       m_nBanks{this, "#banks"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_bad_data{this, "Decoding Error -- dropping all UT hits"};
    ToolHandle<IUTReadoutTool>                          m_readoutTool{this, "ReadoutTool", "UTReadoutTool"};
  };
  // Declaration of the Algorithm Factory
  DECLARE_COMPONENT_WITH_ID( StoreHit<::UT::HitHandler>, "PrStoreUTHit" ) // scalar hits
  DECLARE_COMPONENT_WITH_ID( StoreHit<Hits>, "PrStorePrUTHits" )          // SoA hits

} // namespace LHCb::Pr::UT
