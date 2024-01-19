/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "Event/PrHits.h"
#include "Event/PrSeedTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Event/ZipUtils.h"
#include "FTDAQ/FTInfo.h"
#include "GaudiKernel/StdArrayAsProperty.h"
#include "LHCbAlgs/ScalarTransformer.h"
#include "LHCbAlgs/Transformer.h"
#include "LHCbMath/MatVec.h"
#include "LHCbMath/bit_cast.h"
#include <vector>

/**
 * Converter between LHCb::Pr::Fitted::Forward::Tracks ( SoA PoD ) and vector<Track_v2>
 * @author Sascha Stahl
 */
namespace LHCb::Converters::Track::PrSeeding {
  using Tag = LHCb::Pr::Seeding::Tag;
  using LHCb::Pr::FT::Hits;
  class fromTrackv2PrSeedingTracks
      : public Algorithm::Transformer<LHCb::Pr::Seeding::Tracks( const EventContext& evtCtx,
                                                                 const std::vector<Event::v2::Track>&, const Hits& )> {

  public:
    fromTrackv2PrSeedingTracks( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer( name, pSvcLocator,
                       {KeyValue{"InputTracks", ""}, KeyValue{"InputSciFiHits", PrFTInfo::SciFiHitsLocation}},
                       KeyValue{"OutputTracks", ""} ) {}

    LHCb::Pr::Seeding::Tracks operator()( const EventContext& evtCtx, const std::vector<Event::v2::Track>& inputTracks,
                                          const Hits& fthits ) const override {

      LHCb::Pr::Seeding::Tracks outputTracks{Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx )};
      outputTracks.reserve( inputTracks.size() );

      for ( unsigned int t = 0; t < inputTracks.size(); t++ ) {
        const auto& inTrack  = inputTracks[t];
        auto        outTrack = outputTracks.emplace_back<SIMDWrapper::InstructionSet::Scalar>();

        const auto qop = inTrack.firstState().qOverP();

        outTrack.field<Tag::Chi2PerDoF>().set( static_cast<float>( inTrack.chi2PerDoF() ) );

        // -- copy LHCbIDs
        outTrack.field<Tag::FTHits>().resize( inTrack.lhcbIDs().size() );
        int i = 0;
        for ( auto id : inTrack.lhcbIDs() ) {
          if ( i == LHCb::Pr::Seeding::Tracks::MaxFTHits ) {
            ++m_exceed_max_hits;
            break;
          }
          const auto lhcbid = LHCb::Event::lhcbid_v<SIMDWrapper::scalar::types>( id );
          outTrack.field<Tag::FTHits>()[i].template field<Tag::LHCbID>().set( lhcbid );
          for ( unsigned int ihit = 0; ihit != fthits.size(); ihit++ ) {
            if ( id == fthits.lhcbid( ihit ) ) {
              outTrack.field<Tag::FTHits>()[i].template field<Tag::Index>().set( bit_cast<int>( ihit ) );
              break;
            }
          }
          ++i;
        }

        // -- copy states
        int istate = 0;
        for ( const auto& state : inTrack.states() ) {
          if ( istate == LHCb::Pr::Seeding::Tracks::NumSeedStates ) {
            Transformer::error() << "Reached maximum number of states in LHCb::Pr::Seeding::Tracks "
                                 << LHCb::Pr::Seeding::Tracks::NumSeedStates << "No more states will be added"
                                 << endmsg;
            break;
          }
          auto outState = outTrack.field<Tag::States>( istate );
          outState.setPosition( state.x(), state.y(), state.z() );
          outState.setDirection( state.tx(), state.ty() );
          outState.setQOverP( qop );
          ++istate;
        }
        if ( istate < static_cast<int>( LHCb::Pr::Seeding::Tracks::NumSeedStates ) ) {
          Transformer::error() << "Missing states in LHCb::Pr::Seeding::Tracks got only" << istate << endmsg;
        }
      }
      m_nbTracksCounter += outputTracks.size();
      return outputTracks;
    }

  private:
    mutable Gaudi::Accumulators::SummingCounter<>       m_nbTracksCounter{this, "Nb of converted Tracks"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_exceed_max_hits{
        this, "Reached maximum number of hits in LHCb::Pr::Seeding::Tracks 12. No more hits will be added"};
  }; // namespace

  DECLARE_COMPONENT( fromTrackv2PrSeedingTracks )

} // namespace LHCb::Converters::Track::PrSeeding
