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
#include "Event/PrVeloTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Event/VPLightCluster.h"
#include "Kernel/VPConstants.h"
#include "LHCbAlgs/Transformer.h"
#include <vector>

/**
 * Converter between TracksVP SoA PoD and vector<Track_v2>
 *
 * @author Arthur Hennequin (CERN, LIP6)
 */
namespace {
  using Track  = LHCb::Event::v2::Track;
  using Tracks = LHCb::Pr::Velo::Tracks;
  using TracksProxy =
      typename Tracks::template proxy_type<SIMDWrapper::Scalar, LHCb::Pr::ProxyBehaviour::Contiguous, Tracks const>;

  using dType = SIMDWrapper::scalar::types;
  using I     = dType::int_v;
  using F     = dType::float_v;

  void SetFlagsAndPt( LHCb::Event::v2::Track& outtrack, float ptVelo ) {
    outtrack.setHistory( LHCb::Event::v2::Track::History::PrPixel );
    outtrack.setPatRecStatus( LHCb::Event::v2::Track::PatRecStatus::PatRecIDs );
    const int firstRow = outtrack.lhcbIDs()[0].channelID();
    const int charge   = ( firstRow % 2 == 0 ? -1 : 1 );
    for ( auto& aState : outtrack.states() ) {
      const float tx1    = aState.tx();
      const float ty1    = aState.ty();
      const float slope2 = std::max( tx1 * tx1 + ty1 * ty1, 1.e-20f );
      const float qop    = charge * std::sqrt( slope2 ) / ( ptVelo * std::sqrt( 1.f + slope2 ) );
      aState.setQOverP( qop );
      aState.setErrQOverP2( 1e-6 );
    }
  }

  LHCb::State getState( const TracksProxy& track, int index ) {
    LHCb::State           state;
    LHCb::StateVector     s;
    Gaudi::TrackSymMatrix c;
    // Add state closest to beam
    Vec3<F> pos  = track.StatePos( index );
    Vec3<F> dir  = track.StateDir( index );
    Vec3<F> covX = track.StateCovX( index );
    Vec3<F> covY = track.StateCovY( index );
    s.setX( pos.x.cast() );
    s.setY( pos.y.cast() );
    s.setZ( pos.z.cast() );
    s.setTx( dir.x.cast() );
    s.setTy( dir.y.cast() );
    s.setQOverP( 0. );
    c( 0, 0 ) = covX.x.cast();
    c( 2, 0 ) = covX.y.cast();
    c( 2, 2 ) = covX.z.cast();
    c( 1, 1 ) = covY.x.cast();
    c( 3, 1 ) = covY.y.cast();
    c( 3, 3 ) = covY.z.cast();
    c( 4, 4 ) = 1.f;
    state.setState( s );
    state.setCovariance( c );
    return state;
  }
} // namespace

class TracksVPConverter
    : public LHCb::Algorithm::Transformer<std::vector<LHCb::Event::v2::Track>( const LHCb::Pr::Velo::Tracks& )> {

  Gaudi::Property<float> m_ptVelo{this, "ptVelo", 400 * Gaudi::Units::MeV, "Default pT for Velo tracks"};

public:
  TracksVPConverter( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator, KeyValue{"TracksLocation", "Rec/Track/Velo"},
                     KeyValue{"OutputTracksLocation", "Rec/Track/v2/Velo"} ) {}

  std::vector<LHCb::Event::v2::Track> operator()( const Tracks& tracks ) const override {
    std::vector<LHCb::Event::v2::Track> out;
    out.reserve( tracks.size() );

    m_nbTracksCounter += tracks.size();

    for ( auto const& track : tracks.scalar() ) {
      auto& newTrack = out.emplace_back();

      newTrack.setLhcbIDs( track.lhcbIDs(), LHCb::Tag::Unordered );
      newTrack.states().reserve( 2 );
      auto state_beam = getState( track, 0 );
      state_beam.setLocation( LHCb::State::Location::ClosestToBeam );
      newTrack.addToStates( state_beam );
      auto state_endvelo = getState( track, 1 );
      state_endvelo.setLocation( LHCb::State::Location::EndVelo );
      newTrack.addToStates( state_endvelo );
      newTrack.setType( LHCb::Event::v2::Track::Type::Velo );
      SetFlagsAndPt( newTrack, m_ptVelo );
    }

    return out;
  };

private:
  mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
};

class TracksVPMergerConverter : public LHCb::Algorithm::Transformer<std::vector<LHCb::Event::v2::Track>(
                                    const LHCb::Pr::Velo::Tracks&, const LHCb::Pr::Velo::Tracks& )> {

  Gaudi::Property<float> m_ptVelo{this, "ptVelo", 400 * Gaudi::Units::MeV, "Default pT for Velo tracks"};

public:
  TracksVPMergerConverter( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator, {KeyValue{"TracksForwardLocation", ""}, KeyValue{"TracksBackwardLocation", ""}},
                     KeyValue{"OutputTracksLocation", ""} ) {}

  std::vector<LHCb::Event::v2::Track> operator()( const Tracks& fwd_tracks, const Tracks& bwd_tracks ) const override {
    std::vector<LHCb::Event::v2::Track> out;
    out.reserve( fwd_tracks.size() + bwd_tracks.size() );

    m_nbTracksCounter += fwd_tracks.size() + bwd_tracks.size();

    for ( auto const& fwdtrack : fwd_tracks.scalar() ) {
      auto& newTrack = out.emplace_back();

      newTrack.setLhcbIDs( fwdtrack.lhcbIDs(), LHCb::Tag::Unordered );

      newTrack.states().reserve( 2 );
      auto state_beam = getState( fwdtrack, 0 );
      state_beam.setLocation( LHCb::State::Location::ClosestToBeam );
      newTrack.addToStates( state_beam );
      auto state_endvelo = getState( fwdtrack, 1 );
      state_endvelo.setLocation( LHCb::State::Location::EndVelo );
      newTrack.addToStates( state_endvelo );
      newTrack.setType( LHCb::Event::v2::Track::Type::Velo );
      SetFlagsAndPt( newTrack, m_ptVelo );
    }

    for ( auto const& bwdtrack : bwd_tracks.scalar() ) {
      auto& newTrack = out.emplace_back();
      newTrack.setLhcbIDs( bwdtrack.lhcbIDs(), LHCb::Tag::Unordered );
      newTrack.states().reserve( 1 );
      auto state_beam = getState( bwdtrack, 0 );
      state_beam.setLocation( LHCb::State::Location::ClosestToBeam );
      newTrack.addToStates( state_beam );
      newTrack.setType( LHCb::Event::v2::Track::Type::VeloBackward );
      SetFlagsAndPt( newTrack, m_ptVelo );
    }

    return out;
  };

private:
  mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
};

DECLARE_COMPONENT_WITH_ID( TracksVPConverter, "TracksVPConverter" )
DECLARE_COMPONENT_WITH_ID( TracksVPMergerConverter, "TracksVPMergerConverter" )
