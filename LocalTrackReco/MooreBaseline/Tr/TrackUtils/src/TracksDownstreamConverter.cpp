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

#include "Event/PrDownstreamTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Kernel/VPConstants.h"
#include "LHCbAlgs/Transformer.h"
#include <vector>

using DownTag = LHCb::Pr::Downstream::Tag;

/**
 * Converter between Pr::Downstream::Tracks and vector<Track_v2>
 */
class TracksDownstreamConverter
    : public LHCb::Algorithm::Transformer<std::vector<LHCb::Event::v2::Track>(
          const std::vector<LHCb::Event::v2::Track>&, const LHCb::Pr::Downstream::Tracks& )> {
  using Track  = LHCb::Event::v2::Track;
  using Tracks = LHCb::Pr::Downstream::Tracks;

public:
  TracksDownstreamConverter( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"TracksFTLocation", "Rec/Track/v2/Seed"},
                      KeyValue{"TracksDownstreamLocation", "Rec/Track/Downstream"}},
                     KeyValue{"OutputTracksLocation", "Rec/Track/v2/Downstream"} ) {}

  std::vector<Track> operator()( const std::vector<Track>& tracksFT, const Tracks& tracksUT ) const override {
    std::vector<Track> out;
    out.reserve( tracksUT.size() );

    m_nbTracksCounter += tracksUT.size();
    for ( auto const& track : tracksUT.scalar() ) {
      auto  trackSeed = tracksFT[track.get<DownTag::trackSeed>().cast()];
      auto& newTrack  = out.emplace_back( trackSeed );

      const auto QOverP    = track.qOverP().cast();
      const auto errQOverP = m_stateErrorP * QOverP;

      // Adjust q/p and its uncertainty
      for ( auto& state : newTrack.states() ) {
        state.covariance()( 4, 4 ) = errQOverP * errQOverP;
        state.setQOverP( QOverP );
      }

      // Create a state at zUTa
      auto        s = track.get<DownTag::State>();
      LHCb::State ttState;
      ttState.setState( s.x().cast(), s.y().cast(), s.z().cast(), s.tx().cast(), s.ty().cast(), QOverP );

      Gaudi::TrackSymMatrix cov;
      cov( 0, 0 ) = m_stateErrorX2;
      cov( 1, 1 ) = m_stateErrorY2;
      cov( 2, 2 ) = m_stateErrorTX2;
      cov( 3, 3 ) = m_stateErrorTY2;
      cov( 4, 4 ) = errQOverP * errQOverP;

      ttState.setCovariance( cov );
      newTrack.addToStates( ttState );

      // Add LHCbIds
      for ( auto const lhcbid : track.lhcbIDs() ) { newTrack.addToLhcbIDs( lhcbid ); }

      newTrack.setType( Track::Type::Downstream );
      newTrack.setHistory( Track::History::PrDownstream );
      newTrack.setPatRecStatus( Track::PatRecStatus::PatRecIDs );
      newTrack.addToAncestors( trackSeed );
    }

    return out;
  };

private:
  mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};

  // - StateErrorX2: Error^2 on x-position (for making Track)
  Gaudi::Property<double> m_stateErrorX2{this, "StateErrorX2", 4.0};
  // - StateErrorY2: Error^2 on y-position (for making Track)
  Gaudi::Property<double> m_stateErrorY2{this, "StateErrorY2", 400.};
  // - StateErrorTX2: Error^2 on tx-slope (for making Track)
  Gaudi::Property<double> m_stateErrorTX2{this, "StateErrorTX2", 6.e-5};
  // - StateErrorTY2: Error^2 on ty-slope (for making Track)
  Gaudi::Property<double> m_stateErrorTY2{this, "StateErrorTY2", 1.e-4};
  // - StateErrorP:  Error^2 on momentum (for making Track)
  Gaudi::Property<double> m_stateErrorP{this, "StateErrorP", 0.15};
};

DECLARE_COMPONENT( TracksDownstreamConverter )
