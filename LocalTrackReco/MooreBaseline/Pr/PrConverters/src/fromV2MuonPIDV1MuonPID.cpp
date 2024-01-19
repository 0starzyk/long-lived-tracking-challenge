/*****************************************************************************\
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

// Gaudi
#include "LHCbAlgs/Transformer.h"

// LHCb
#include "Event/MuonPID.h"
#include "Event/MuonPIDs_v2.h"
#include "Event/Track.h"

/**
 * Converter from MuonPID SoA PoD to LHCb::MuonPIDs (KeyedContainer)
 *
 * @author Ricardo Vazquez Gomez (UB)
 *
 * Based on https://gitlab.cern.ch/lhcb/Rec/blob/master/Pr/PrConverters/src/fromPrVeloUTTrack.cpp
 * from Michel De Cian
 */

namespace LHCb::Converters::Muon {

  class fromV2MuonPIDV1MuonPID : public Algorithm::Transformer<LHCb::MuonPIDs(
                                     const Event::v2::Muon::PIDs&, const LHCb::Tracks&, const LHCb::Tracks& )> {

  public:
    fromV2MuonPIDV1MuonPID( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer( name, pSvcLocator,
                       {KeyValue{"InputMuonPIDs", LHCb::MuonPIDLocation::Default},
                        KeyValue{"InputTracks", LHCb::TrackLocation::Default}, KeyValue{"InputMuonTracks", ""}},
                       KeyValue{"OutputMuonPIDs", "Rec/Muon/MuonPID"} ) {}

    LHCb::MuonPIDs operator()( const Event::v2::Muon::PIDs& muonPIDs, const LHCb::Tracks& tracks,
                               const LHCb::Tracks& muontracks ) const override {
      LHCb::MuonPIDs out;
      out.reserve( muonPIDs.size() );

      m_nbMuonPIDsCounter += muonPIDs.size();

      if ( !tracks.empty() && !muonPIDs.empty() ) {
        auto offsetTrack = std::find_if( tracks.begin(), tracks.end(),
                                         [&]( const auto& tk ) { return tk->type() == m_track_type.value(); } );
        auto offsetIndex = ( *offsetTrack )->index();

        for ( auto const& muonPID : muonPIDs.scalar() ) {
          auto newMuonPID = new LHCb::MuonPID();
          newMuonPID->setChi2Corr( muonPID.Chi2Corr().cast() );
          newMuonPID->setIsMuon( muonPID.IsMuon().cast() );
          newMuonPID->setIsMuonTight( muonPID.IsMuonTight().cast() );
          newMuonPID->setInAcceptance( muonPID.InAcceptance().cast() );
          newMuonPID->setPreSelMomentum( muonPID.PreSelMomentum().cast() );
          newMuonPID->setMuonLLMu( muonPID.LLMu().cast() );
          newMuonPID->setMuonLLBg( muonPID.LLBg().cast() );
          newMuonPID->setMuonMVA2( muonPID.CatBoost().cast() );
          auto index = muonPID.indices();
          // get the correct track
          auto it = tracks.begin() + index.cast() + offsetIndex;
          newMuonPID->setIDTrack( ( *it ) );
          auto muonit = muontracks.begin() + index.cast();
          newMuonPID->setMuonTrack( *muonit );
          out.insert( newMuonPID );
        }
      }
      return out;
    };

  private:
    Gaudi::Property<Event::Enum::Track::Type> m_track_type{
        this, "RestrictToType", Event::Enum::Track::Type::Unknown,
        "If set, filter the input tracks and only write those of the given type. Otherwise the full set is processed"};
    mutable Gaudi::Accumulators::SummingCounter<> m_nbMuonPIDsCounter{this, "Nb of Produced MuonPIDs"};
  };

  DECLARE_COMPONENT_WITH_ID( fromV2MuonPIDV1MuonPID, "fromV2MuonPIDV1MuonPID" )
} // namespace LHCb::Converters::Muon
