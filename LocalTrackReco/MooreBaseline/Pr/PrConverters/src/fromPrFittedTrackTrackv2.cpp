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
#include "Event/PrFittedForwardTracks.h"
#include "Event/PrLongTracks.h"
#include "Event/PrUpstreamTracks.h"
#include "Event/PrVeloTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "GaudiKernel/StdArrayAsProperty.h"
#include "LHCbAlgs/Transformer.h"
#include "SelKernel/TrackZips.h"
#include <vector>

/**
 * Converter between LHCb::Pr::Fitted::Forward::Tracks ( SoA PoD ) and vector<Track_v2>
 *
 *
 * @author Sascha Stahl
 */

namespace {
  using dType = SIMDWrapper::scalar::types;
  using I     = dType::int_v;
  using F     = dType::float_v;
  constexpr std::array<float, 5> default_covarianceValues{4.0, 400.0, 4.e-6, 1.e-4, 0.1};

  template <typename ZippedFittedTrackType>
  auto& get_fitted_tracks( const ZippedFittedTrackType& fitted_tracks ) {
    return fitted_tracks.template get<LHCb::Pr::Fitted::Forward::Tracks>();
  }

  template <>
  auto& get_fitted_tracks( const LHCb::Pr::Fitted::Forward::Tracks& fitted_tracks ) {
    return fitted_tracks;
  }

  LHCb::State get_state( Vec3<F> const& pos, Vec3<F> const& dir, F const qop, Vec3<F> const& covX, Vec3<F> const& covY,
                         F const qopError ) {
    LHCb::State state;
    state.setState( pos.x.cast(), pos.y.cast(), pos.z.cast(), dir.x.cast(), dir.y.cast(), qop.cast() );
    state.covariance()( 0, 0 ) = covX.x.cast();
    state.covariance()( 0, 2 ) = covX.y.cast();
    state.covariance()( 2, 2 ) = covX.z.cast();
    state.covariance()( 1, 1 ) = covY.x.cast();
    state.covariance()( 1, 3 ) = covY.y.cast();
    state.covariance()( 3, 3 ) = covY.z.cast();
    state.covariance()( 4, 4 ) = qopError.cast();
    return state;
  }

  std::vector<LHCb::Event::v2::Track> convert_tracks( LHCb::Pr::Long::Tracks const&            forward_tracks,
                                                      LHCb::Pr::Fitted::Forward::Tracks const& fitted_tracks,
                                                      std::array<float, 5> const               covarianceValues ) {
    std::vector<LHCb::Event::v2::Track> out;
    out.reserve( fitted_tracks.size() );

    auto const fwdzipped    = forward_tracks.scalar();
    auto const fitfwdzipped = fitted_tracks.scalar();
    for ( auto const& fitfwd : fitfwdzipped ) {
      auto  forward_track_index = fitfwd.trackSeed().cast();
      auto& newTrack            = out.emplace_back();
      // set track flags
      newTrack.setType( LHCb::Event::v2::Track::Type::Long );
      newTrack.setHistory( LHCb::Event::v2::Track::History::PrForward );
      newTrack.setPatRecStatus( LHCb::Event::v2::Track::PatRecStatus::PatRecIDs );
      newTrack.setFitStatus( LHCb::Event::v2::Track::FitStatus::Fitted );
      // get momentum
      F qop      = fitfwd.qOverP();
      F qopError = covarianceValues[4] * qop * qop;

      // closest to beam state
      LHCb::State closesttobeam_state = get_state( fitfwd.closestToBeamStatePos(), fitfwd.closestToBeamStateDir(), qop,
                                                   fitfwd.covX(), fitfwd.covY(), qopError );
      closesttobeam_state.setLocation( LHCb::State::Location::ClosestToBeam );

      const Vec3<F> covX{covarianceValues[0], 0.f, covarianceValues[2]};
      const Vec3<F> covY{covarianceValues[1], 0.f, covarianceValues[3]};

      // scifi state
      LHCb::State scifi_state = get_state( fwdzipped[forward_track_index].StatePos( 1 ),
                                           fwdzipped[forward_track_index].StateDir( 1 ), qop, covX, covY, qopError );
      scifi_state.setLocation( LHCb::State::Location::AtT );

      // add states
      newTrack.addToStates( closesttobeam_state );
      newTrack.addToStates( scifi_state );

      // set chi2 / chi2ndof
      newTrack.setChi2PerDoF( LHCb::Event::v2::Track::Chi2PerDoF{fitfwd.chi2().cast(), fitfwd.chi2nDoF().cast()} );

      // If we rely on pointers internally stored in the classes we can take it from fitted tracks
      auto lhcbids = fwdzipped[forward_track_index].lhcbIDs();
      newTrack.addToLhcbIDs( lhcbids, LHCb::Tag::Unordered_tag{} );
    }
    return out;
  }

} // namespace

namespace LHCb::Converters::Track::v2 {
  template <typename FittedTrackType>
  class fromPrFittedForwardTrack
      : public Algorithm::Transformer<std::vector<Event::v2::Track>( const FittedTrackType& )> {

  public:
    using base_class = Algorithm::Transformer<std::vector<Event::v2::Track>( const FittedTrackType& )>;
    using KeyValue   = typename base_class::KeyValue;

    fromPrFittedForwardTrack( const std::string& name, ISvcLocator* pSvcLocator )
        : base_class( name, pSvcLocator, {KeyValue{"FittedTracks", ""}}, KeyValue{"OutputTracks", ""} ) {}
    Gaudi::Property<std::array<float, 5>> m_covarianceValues{this, "covarianceValues", default_covarianceValues};

    std::vector<Event::v2::Track> operator()( const FittedTrackType& fitted_tracks_like ) const override {

      auto const& fitted_tracks  = get_fitted_tracks( fitted_tracks_like );
      auto const* forward_tracks = fitted_tracks.getForwardAncestors();
      if ( forward_tracks == nullptr ) {
        base_class::error()
            << "Forward tracks container of Fitted forward tracks does not exist. Conversion to track v2 will not work."
            << endmsg;
        return std::vector<Event::v2::Track>{};
      }
      std::vector<Event::v2::Track> out = convert_tracks( *forward_tracks, fitted_tracks, m_covarianceValues );
      m_nbTracksCounter += out.size();
      return out;
    }

  private:
    mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
  }; // namespace

  DECLARE_COMPONENT_WITH_ID( fromPrFittedForwardTrack<Pr::Fitted::Forward::Tracks>,
                             "LHCb__Converters__Track__v2__fromPrFittedForwardTrack" )
  DECLARE_COMPONENT_WITH_ID( fromPrFittedForwardTrack<Pr::Fitted::Forward::TracksWithMuonID>,
                             "LHCb__Converters__Track__v2__fromPrFittedForwardTrackWithMuonID" )
  DECLARE_COMPONENT_WITH_ID( fromPrFittedForwardTrack<Pr::Fitted::Forward::TracksWithPVs>,
                             "LHCb__Converters__Track__v2__fromPrFittedForwardTrackWithPVs" )

} // namespace LHCb::Converters::Track::v2
