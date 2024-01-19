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

#include "Event/PrFittedForwardTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "LHCbAlgs/Transformer.h"
#include <vector>

/**
 * Converter between TracksFit SoA PoD and vector<Track_v2>
 *
 * @author Arthur Hennequin (CERN, LIP6)
 */

class TracksFitConverter : public LHCb::Algorithm::Transformer<std::vector<LHCb::Event::v2::Track>(
                               const std::vector<LHCb::Event::v2::Track>&, const LHCb::Pr::Fitted::Forward::Tracks& )> {
  using Track  = LHCb::Event::v2::Track;
  using Tracks = LHCb::Pr::Fitted::Forward::Tracks;

public:
  TracksFitConverter( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"TracksFTLocation", "Rec/Track/v2/FT"}, KeyValue{"TracksFitLocation", "Rec/Track/Fit"}},
                     KeyValue{"OutputTracksLocation", "Rec/Track/v2/Fit"} ) {}

  StatusCode initialize() override {
    StatusCode sc = Transformer::initialize();
    if ( sc.isFailure() ) return sc;
    return StatusCode::SUCCESS;
  };

  std::vector<Track> operator()( const std::vector<Track>& tracksFT, const Tracks& tracksFit ) const override {
    std::vector<Track> out;
    out.reserve( tracksFit.size() );

    m_nbTracksCounter += tracksFit.size();

    using dType = SIMDWrapper::scalar::types;
    using F     = dType::float_v;

    for ( auto const& track : tracksFit.scalar() ) {
      auto  trackFT  = tracksFT[track.trackSeed().cast()];
      auto& newTrack = out.emplace_back( trackFT );

      // set q/p in all of the existing states
      F qop = track.qOverP();
      for ( auto& state : newTrack.states() ) state.setQOverP( qop.cast() );

      // update closest to beam state
      const Vec3<F> pos = track.closestToBeamStatePos();
      const Vec3<F> dir = track.closestToBeamStateDir();

      auto beamState = newTrack.stateAt( LHCb::State::Location::ClosestToBeam );
      beamState->setState( pos.x.cast(), pos.y.cast(), pos.z.cast(), dir.x.cast(), dir.y.cast(), qop.cast() );

      // update cov
      const Vec3<F> covX = track.covX();
      const Vec3<F> covY = track.covY();

      beamState->covariance()( 0, 0 ) = covX.x.cast();
      beamState->covariance()( 0, 2 ) = covX.y.cast();
      beamState->covariance()( 2, 2 ) = covX.z.cast();
      beamState->covariance()( 1, 1 ) = covY.x.cast();
      beamState->covariance()( 1, 3 ) = covY.y.cast();
      beamState->covariance()( 3, 3 ) = covY.z.cast();

      // set chi2 / chi2ndof
      newTrack.setChi2PerDoF( LHCb::Event::v2::Track::Chi2PerDoF{track.chi2().cast(), track.chi2nDoF().cast()} );

      // set history
      newTrack.addToAncestors( trackFT );
    }

    return out;
  };

private:
  mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
};

DECLARE_COMPONENT( TracksFitConverter )
