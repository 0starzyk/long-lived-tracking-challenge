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

#include <vector>

// Gaudi
#include "LHCbAlgs/Transformer.h"

// LHCb
#include "Event/StateParameters.h"
#include "Event/Track.h"

#include "Event/PrFittedForwardTracks.h"
#include "Event/PrHits.h"
#include "Event/PrLongTracks.h"
#include "Event/PrVeloTracks.h"

#include "VeloKalmanHelpers.h"

/**
 * Velo only Kalman fit
 *
 * @author Arthur Hennequin (CERN, LIP6)
 */
namespace LHCb::Pr::Velo {

  using TrackTag = LHCb::Pr::Fitted::Forward::Tag;
  using VP::Hits;

  class Kalman
      : public LHCb::Algorithm::Transformer<Fitted::Forward::Tracks(
            const EventContext&, LHCb::UniqueIDGenerator const&, const Hits&, const Tracks&, const Long::Tracks& )> {
    using TracksVP  = Tracks;
    using TracksFT  = Long::Tracks;
    using TracksFit = Fitted::Forward::Tracks;
    using simd      = SIMDWrapper::best::types;
    using I         = simd::int_v;
    using F         = simd::float_v;

  public:
    Kalman( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer( name, pSvcLocator,
                       {KeyValue{"UniqueIDGenerator", LHCb::UniqueIDGeneratorLocation::Default},
                        KeyValue{"HitsLocation", "Raw/VP/Hits"}, KeyValue{"TracksVPLocation", "Rec/Track/VP"},
                        KeyValue{"TracksFTLocation", "Rec/Track/FT"}},
                       KeyValue{"OutputTracksLocation", "Rec/Track/Fit"} ) {}

    StatusCode initialize() override {
      StatusCode sc = Transformer::initialize();
      if ( sc.isFailure() ) return sc;
      return StatusCode::SUCCESS;
    };

    TracksFit operator()( const EventContext& evtCtx, LHCb::UniqueIDGenerator const& unique_id_gen, const Hits& hits,
                          const TracksVP& tracksVP, const TracksFT& tracksFT ) const override {
      // Forward tracks and its fit are zipable as there is a one to one correspondence.
      TracksFit out{&tracksFT, unique_id_gen, tracksFT.history(), tracksFT.zipIdentifier(),
                    LHCb::getMemResource( evtCtx )};
      out.reserve( tracksFT.size() );
      m_nbTracksCounter += tracksFT.size();

      for ( auto const& track : tracksFT.simd() ) {
        auto       loop_mask = track.loop_mask();
        auto const idxVP     = track.trackVP();
        auto const qop       = track.qOverP();

        auto [stateInfo, chi2, nDof] = fitBackwardWithMomentum( loop_mask, tracksVP, idxVP, qop, hits, 0 );

        auto outTrack = out.emplace_back( tracksFT.size() );

        outTrack.field<TrackTag::trackSeed>().set( track.indices() );
        outTrack.field<TrackTag::State>().setQOverP( qop );
        outTrack.field<TrackTag::Chi2>().set( chi2 / F( nDof ) );
        outTrack.field<TrackTag::Chi2nDoF>().set( nDof );
        outTrack.field<TrackTag::State>().setPosition( stateInfo.pos().x, stateInfo.pos().y, stateInfo.pos().z );
        outTrack.field<TrackTag::State>().setDirection( stateInfo.dir().x, stateInfo.dir().y );

        outTrack.field<TrackTag::StateCovX>( 0 ).set( stateInfo.covX().x );
        outTrack.field<TrackTag::StateCovX>( 1 ).set( stateInfo.covX().y );
        outTrack.field<TrackTag::StateCovX>( 2 ).set( stateInfo.covX().z );

        outTrack.field<TrackTag::StateCovY>( 0 ).set( stateInfo.covY().x );
        outTrack.field<TrackTag::StateCovY>( 1 ).set( stateInfo.covY().y );
        outTrack.field<TrackTag::StateCovY>( 2 ).set( stateInfo.covY().z );
      }

      return out;
    };

  private:
    mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
  };
} // namespace LHCb::Pr::Velo

DECLARE_COMPONENT_WITH_ID( LHCb::Pr::Velo::Kalman, "VeloKalman" )
