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

#include "TrackInterfaces/IPrFitterTool.h"

#include "Event/PartialChiSquareds.h"
#include "Event/PrDownstreamTracks.h"
#include "Event/PrKalmanFitResult.h"
#include "Event/PrLongTracks.h"
#include "Event/PrSciFiHits.h"
#include "Event/PrSeedTracks.h"
#include "Event/PrVeloHits.h"
#include "Event/PrVeloTracks.h"
#include "Event/Track_v1.h"
#include "Event/Track_v3.h"

#include "LHCbAlgs/LHCbAlgsHelpers.h"
#include "PrKalmanFilter/KF.h"

#include "Gaudi/Accumulators.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Transformer.h"

namespace {

  using namespace LHCb::Pr::Tracks::Fit;
  using TrackV1 = LHCb::Event::v1::Track;
} // namespace

template <typename InputTrackType, typename OutputTrackType>
class KalmanFilterToolExample
    : public LHCb::Algorithm::MultiTransformer<OutputTrackType(
          InputTrackType const&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
          const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&, const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&,
          const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>&, const IPrFitterTool& )> {
public:
  using base_t = LHCb::Algorithm::MultiTransformer<OutputTrackType(
      const InputTrackType&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&, const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&, const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>&,
      const IPrFitterTool& )>;

  using base_t::info;
  KalmanFilterToolExample( const std::string& name, ISvcLocator* pSvcLocator )
      : base_t( name, pSvcLocator,
                {typename base_t::KeyValue{"Input", ""}, typename base_t::KeyValue{"HitsVP", ""},
                 typename base_t::KeyValue{"HitsUT", ""}, typename base_t::KeyValue{"HitsFT", ""},
                 typename base_t::KeyValue{"HitsMuon", ""}, typename base_t::KeyValue{"TrackFitter", ""}},
                LHCb::Algorithm::IOHelper<InputTrackType, OutputTrackType, base_t>::OutputNames(
                    get_out_names<OutputTrackType>() ) ) {}

  OutputTrackType operator()( InputTrackType const&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                              const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                              const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&,
                              const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>&, const IPrFitterTool& ) const override;
};

template <typename InputTrackType, typename OutputTrackType>
class KalmanFilterToolExample_noUT : public LHCb::Algorithm::MultiTransformer<OutputTrackType(
                                         InputTrackType const&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>&, const IPrFitterTool& )> {
public:
  using base_t = LHCb::Algorithm::MultiTransformer<OutputTrackType(
      const InputTrackType&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&, const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&,
      const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>&, const IPrFitterTool& )>;

  using base_t::info;
  KalmanFilterToolExample_noUT( const std::string& name, ISvcLocator* pSvcLocator )
      : base_t( name, pSvcLocator,
                {typename base_t::KeyValue{"Input", ""}, typename base_t::KeyValue{"HitsVP", ""},
                 typename base_t::KeyValue{"HitsFT", ""}, typename base_t::KeyValue{"HitsMuon", ""},
                 typename base_t::KeyValue{"TrackFitter", ""}},
                LHCb::Algorithm::IOHelper<InputTrackType, OutputTrackType, base_t>::OutputNames(
                    get_out_names<OutputTrackType>() ) ) {}

  OutputTrackType operator()( InputTrackType const&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                              const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&,
                              const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>&, const IPrFitterTool& ) const override;
};

template <typename InputTrackType, typename OutputTrackType>
OutputTrackType KalmanFilterToolExample<InputTrackType, OutputTrackType>::
                operator()( InputTrackType const& tracks, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>& hits_vp,
            const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>& hits_ut, const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& hits_ft,
            const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>& hits_muon, const IPrFitterTool& fit ) const {
  if constexpr ( std::is_same_v<OutputTrackType, V1Output> ) {
    if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Velo::Tracks> ) {
      // by construction same as PrKalmanFilter_Velo
      return fit( tracks, hits_vp );
    } else if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Downstream::Tracks> ) {
      // by construction same as PrKalmanFilter_Downstream
      return fit( tracks, hits_ut, hits_ft );
    } else if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Long::Tracks> ) {
      // by construction same as PrKalmanFilter
      return fit( tracks, hits_vp, hits_ut, hits_ft );
    } else if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Seeding::Tracks> ) {
      // by construction same as PrKalmanFilter_Seed
      return fit( tracks, hits_ft );
    } else if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Upstream::Tracks> ) {
      // by construction same as PrKalmanFilter_Upstream
      return fit( tracks, hits_vp, hits_ut );
    } else {
      // v1 tracks as input!
      // let's first fit the whole container
      if ( tracks.empty() ) return V1Output{};
      const auto* single_track = *tracks.begin();
      auto        tracks_v1    = [&] {
        if ( single_track->checkType( TrackV1::Types::Velo ) ) {
          return fit( tracks, hits_vp );
        } else if ( single_track->checkType( TrackV1::Types::Downstream ) ) {
          return fit( tracks, hits_ut, hits_ft );
        } else if ( single_track->checkType( TrackV1::Types::Ttrack ) ) {
          return fit( tracks, hits_ft );
        } else if ( single_track->checkType( TrackV1::Types::Upstream ) ) {
          return fit( tracks, hits_vp, hits_ut );
        } else if ( single_track->checkType( TrackV1::Types::LongMuon ) ) {
          return fit( tracks, hits_vp, hits_ut, hits_ft, hits_muon );
        } else {
          // long tracks ... but in principle this operator is as well suitable for velo, downstream and upstream, cool!
          return fit( tracks, hits_vp, hits_ut, hits_ft );
        }
      }();

      // we can also fit single v1 tracks, nice!
      auto track_v1 = [&] {
        if ( single_track->checkType( TrackV1::Types::Velo ) ) {
          return fit( *single_track, hits_vp );
        } else if ( single_track->checkType( TrackV1::Types::Downstream ) ) {
          return fit( *single_track, hits_ut, hits_ft );
        } else if ( single_track->checkType( TrackV1::Types::Ttrack ) ) {
          return fit( *single_track, hits_ft );
        } else if ( single_track->checkType( TrackV1::Types::Upstream ) ) {
          return fit( *single_track, hits_vp, hits_ut );
        } else if ( single_track->checkType( TrackV1::Types::LongMuon ) ) {
          return fit( *single_track, hits_vp, hits_ut, hits_ft, hits_muon );
        } else {
          // long tracks ... but in principle this operator is as well suitable for velo, downstream and upstream, cool!
          return fit( *single_track, hits_vp, hits_ut, hits_ft );
        }
      }();
      // let's add the track to the container, why not? one reason not to might be that track is not valid :)
      if ( !track_v1.checkFlag( TrackV1::Flags::Invalid ) )
        std::get<TracksV1>( tracks_v1 ).add( new TrackV1( std::move( track_v1 ) ) );
      return tracks_v1;
    }
  } else if constexpr ( std::is_same_v<OutputTrackType, V3FullOutput> ) {

    static_assert( std::is_same_v<InputTrackType, LHCb::Pr::Velo::Tracks> ||
                   std::is_same_v<InputTrackType, LHCb::Pr::Downstream::Tracks> ||
                   std::is_same_v<InputTrackType, LHCb::Pr::Long::Tracks> ||
                   std::is_same_v<InputTrackType, LHCb::Pr::Seeding::Tracks> ||
                   std::is_same_v<InputTrackType, LHCb::Pr::Upstream::Tracks> );

    if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Velo::Tracks> ) {
      // by construction same as PrKalmanFilter_Velo
      return fit.fitted_v3_tracks( tracks, hits_vp );
    } else if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Downstream::Tracks> ) {
      // by construction same as PrKalmanFilter_Downstream
      return fit.fitted_v3_tracks( tracks, hits_ut, hits_ft );
    } else if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Long::Tracks> ) {
      // by construction same as PrKalmanFilter
      return fit.fitted_v3_tracks( tracks, hits_vp, hits_ut, hits_ft );
    } else if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Seeding::Tracks> ) {
      // by construction same as PrKalmanFilter_Seed
      return fit.fitted_v3_tracks( tracks, hits_ft );
    } else if constexpr ( std::is_same_v<InputTrackType, LHCb::Pr::Upstream::Tracks> ) {
      // by construction same as PrKalmanFilter_Upstream
      return fit.fitted_v3_tracks( tracks, hits_vp, hits_ut );
    }
  }
}

template <typename InputTrackType, typename OutputTrackType>
OutputTrackType KalmanFilterToolExample_noUT<InputTrackType, OutputTrackType>::
                operator()( InputTrackType const& tracks, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>& hits_vp,
            const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&   hits_ft,
            const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>& hits_muon, const IPrFitterTool& fit ) const {
  if constexpr ( std::is_same_v<OutputTrackType, V1Output> ) {
    // v1 tracks as input!
    // let's first fit the whole container
    if ( tracks.empty() ) return V1Output{};
    const auto* single_track = *tracks.begin();
    auto        tracks_v1    = [&] {
      if ( single_track->checkType( TrackV1::Types::Long ) ) {
        return fit( tracks, hits_vp, hits_ft );
      } else if ( single_track->checkType( TrackV1::Types::LongMuon ) ) {
        return fit( tracks, hits_vp, hits_ft, hits_muon );
      } else {
        return V1Output{};
      }
    }();
    return tracks_v1;
  }
}

using PrKalmanFilterToolExampleAlgo            = KalmanFilterToolExample<LHCb::Pr::Long::Tracks, V1Output>;
using PrKalmanFilterToolExampleAlgo_Downstream = KalmanFilterToolExample<LHCb::Pr::Downstream::Tracks, V1Output>;
using PrKalmanFilterToolExampleAlgo_Seed       = KalmanFilterToolExample<LHCb::Pr::Seeding::Tracks, V1Output>;
using PrKalmanFilterToolExampleAlgo_Velo       = KalmanFilterToolExample<LHCb::Pr::Velo::Tracks, V1Output>;
using PrKalmanFilterToolExampleAlgo_Upstream   = KalmanFilterToolExample<LHCb::Pr::Upstream::Tracks, V1Output>;
using PrKalmanFilterToolExampleAlgo_V1V1       = KalmanFilterToolExample<LHCb::Event::v1::Tracks, V1Output>;

using PrKalmanFilterToolExampleAlgo_V3            = KalmanFilterToolExample<LHCb::Pr::Long::Tracks, V3FullOutput>;
using PrKalmanFilterToolExampleAlgo_Downstream_V3 = KalmanFilterToolExample<LHCb::Pr::Downstream::Tracks, V3FullOutput>;
using PrKalmanFilterToolExampleAlgo_Seed_V3       = KalmanFilterToolExample<LHCb::Pr::Seeding::Tracks, V3FullOutput>;
using PrKalmanFilterToolExampleAlgo_Velo_V3       = KalmanFilterToolExample<LHCb::Pr::Velo::Tracks, V3FullOutput>;
using PrKalmanFilterToolExampleAlgo_Upstream_V3   = KalmanFilterToolExample<LHCb::Pr::Upstream::Tracks, V3FullOutput>;

using PrKalmanFilterToolExampleAlgo_V1V1_noUT = KalmanFilterToolExample_noUT<LHCb::Event::v1::Tracks, V1Output>;

DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo, "PrKalmanFilterToolExampleAlgo" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_Downstream, "PrKalmanFilterToolExampleAlgo_Downstream" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_Seed, "PrKalmanFilterToolExampleAlgo_Seed" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_Velo, "PrKalmanFilterToolExampleAlgo_Velo" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_Upstream, "PrKalmanFilterToolExampleAlgo_Upstream" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_V1V1, "PrKalmanFilterToolExampleAlgo_V1V1" )

DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_V3, "PrKalmanFilterToolExampleAlgo_V3" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_Downstream_V3, "PrKalmanFilterToolExampleAlgo_Downstream_V3" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_Seed_V3, "PrKalmanFilterToolExampleAlgo_Seed_V3" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_Velo_V3, "PrKalmanFilterToolExampleAlgo_Velo_V3" )
DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_Upstream_V3, "PrKalmanFilterToolExampleAlgo_Upstream_V3" )

DECLARE_COMPONENT_WITH_ID( PrKalmanFilterToolExampleAlgo_V1V1_noUT, "PrKalmanFilterToolExampleAlgo_V1V1_noUT" )
