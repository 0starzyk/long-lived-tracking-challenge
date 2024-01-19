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
#include "Detector/UT/ChannelID.h"
#include "Detector/VP/VPChannelID.h"
#include "Event/PrLongTracks.h"
#include "Event/PrVeloTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "GaudiKernel/StdArrayAsProperty.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Transformer.h"
#include <vector>

/**
 * Converter between TracksFT SoA PoD and vector<Track_v2>
 *
 * @author Arthur Hennequin (CERN, LIP6)
 */

class TracksMatchConverter : public LHCb::Algorithm::Transformer<std::vector<LHCb::Event::v2::Track>(
                                 const std::vector<LHCb::Event::v2::Track>&, const std::vector<LHCb::Event::v2::Track>&,
                                 const LHCb::Pr::Long::Tracks& )> {
  using Track  = LHCb::Event::v2::Track;
  using Tracks = LHCb::Pr::Long::Tracks;
  // From PrGeometryTool in PrAlgorithms

public:
  TracksMatchConverter( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"TracksSeedLocation", "Rec/Track/Seed"},
                      KeyValue{"TracksVeloLocation", "Rec/Track/Velo"},
                      KeyValue{"TracksMatchLocation", "Rec/Track/MatchSOA"}},
                     KeyValue{"OutputTracksLocation", "Rec/Track/Match"} ) {}

  Gaudi::Property<std::array<float, 5>> m_covarianceValues{this, "covarianceValues", {4.0, 400.0, 4.e-6, 1.e-4, 0.1}};

  std::vector<Track> operator()( const std::vector<Track>& tracksSeed, const std::vector<Track>& tracksVelo,
                                 const Tracks& tracksMatch ) const override {
    std::vector<Track> out;
    out.reserve( tracksMatch.size() );
    m_nbTracksCounter += tracksMatch.size();

    for ( auto const& track : tracksMatch.scalar() ) {
      auto& trackSeed = tracksSeed[track.trackSeed().cast()];
      auto& trackVelo = tracksVelo[track.trackVP().cast()];
      auto& newTrack  = out.emplace_back( trackVelo );
      newTrack.addToAncestors( trackSeed );
      newTrack.addToAncestors( trackVelo );

      for ( auto& state : trackSeed.states() ) { newTrack.addToStates( state ); }

      // set q/p in all of the existing states
      auto const qop     = track.qOverP().cast();
      auto const errQop2 = m_covarianceValues[4] * qop * qop;

      for ( auto& state : newTrack.states() ) {
        state.setQOverP( qop );
        state.setErrQOverP2( errQop2 );
      }

      /// add lhcbIDs
      newTrack.setLhcbIDs( track.lhcbIDs(), LHCb::Tag::Unordered );

      newTrack.setType( Track::Type::Long );
      newTrack.setHistory( Track::History::PrMatch );
      newTrack.setPatRecStatus( Track::PatRecStatus::PatRecIDs );
    }

    return out;
  };

private:
  mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
};

DECLARE_COMPONENT( TracksMatchConverter )
