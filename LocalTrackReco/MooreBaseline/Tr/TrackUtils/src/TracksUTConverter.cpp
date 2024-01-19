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
#include "Event/PrUpstreamTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Kernel/VPConstants.h"
#include "LHCbAlgs/Transformer.h"
#include <vector>

/**
 * Converter between TracksUT SoA PoD and vector<Track_v2>
 *
 * @author Arthur Hennequin (CERN, LIP6)
 *
 * Based on https://gitlab.cern.ch/lhcb/Rec/blob/master/Pr/PrConverters/src/fromPrVeloUTTrack.cpp
 * from Michel De Cian
 */

class TracksUTConverter : public LHCb::Algorithm::Transformer<std::vector<LHCb::Event::v2::Track>(
                              const std::vector<LHCb::Event::v2::Track>&, const LHCb::Pr::Upstream::Tracks& )> {
  using Track  = LHCb::Event::v2::Track;
  using Tracks = LHCb::Pr::Upstream::Tracks;

public:
  TracksUTConverter( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"TracksVPLocation", "Rec/Track/v2/Velo"}, KeyValue{"TracksUTLocation", "Rec/Track/UT"}},
                     KeyValue{"OutputTracksLocation", "Rec/Track/v2/UT"} ) {}

  std::vector<Track> operator()( const std::vector<Track>& tracksVP, const Tracks& tracksUT ) const override {
    std::vector<Track> out;
    out.reserve( tracksUT.size() );

    m_nbTracksCounter += tracksUT.size();
    for ( auto const& track : tracksUT.scalar() ) {
      auto  trackVP  = tracksVP[track.trackVP().cast()];
      auto& newTrack = out.emplace_back( trackVP );

      // set q/p in all of the existing states
      for ( auto& state : newTrack.states() ) state.setQOverP( track.qOverP().cast() );

      // Add LHCbIds
      newTrack.setLhcbIDs( track.lhcbIDs(), LHCb::Tag::Unordered );

      // As we don't need the state in the UT, it is not added in PrVeloUT
      // and can't be added here.
      newTrack.setType( Track::Type::Upstream );
      newTrack.setHistory( Track::History::PrVeloUT ); // Track::History::PatVeloTT
      newTrack.addToAncestors( trackVP );
      newTrack.setPatRecStatus( Track::PatRecStatus::PatRecIDs );
    }

    return out;
  };

private:
  mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
};

DECLARE_COMPONENT( TracksUTConverter )
