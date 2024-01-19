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

#include "PrKalmanFilter/KF.h"
// Rec
#include "Event/ParametrisedScatters.h"
#include "Event/PrFitNode.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
// LHCb
#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/ChiSquare.h"
#include "Event/PartialChiSquareds.h"
#include "Event/PrDownstreamTracks.h"
#include "Event/PrLongTracks.h"
#include "Event/PrSeedTracks.h"
#include "Event/PrTracksTag.h"
#include "Event/PrVeloTracks.h"
#include "Event/SIMDEventTypes.h"
#include "Event/SOACollection.h"
#include "Event/SOAUtils.h"
#include "Event/StateParameters.h"
#include "Event/StateVector.h"
#include "Event/Track.h"
#include "Event/TrackTypes.h"
#include "Event/Track_v3.h"
#include "Event/UniqueIDGenerator.h"
#include "Kernel/STLExtensions.h"
#include "LHCbAlgs/LHCbAlgsHelpers.h"
#include "LHCbAlgs/Transformer.h"
#include "LHCbMath/LHCbMath.h"
#include "LHCbMath/SIMDWrapper.h"
#include "LHCbMath/Similarity.h"
// Gaudi
#include "Gaudi/Accumulators.h"
#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/ToolHandle.h"
// std
#include <algorithm>
#include <string>
#include <type_traits>

namespace LHCb::Pr {

  namespace {
    using namespace LHCb::Pr::Tracks::Fit;

    template <HitType T>
    std::string getKey();
    template <>
    std::string getKey<HitType::VP>() {
      return "HitsVP";
    }
    template <>
    std::string getKey<HitType::UT>() {
      return "HitsUT";
    }
    template <>
    std::string getKey<HitType::FT>() {
      return "HitsFT";
    }
  } // namespace

  template <typename InputTrackType, typename OutputTrackType, HitType... Types>
  class KalmanFilter
      : public LHCb::Algorithm::MultiTransformer<OutputTrackType( InputTrackType const&, const Hits<Types>&...,
                                                                  DetectorElement const&,
                                                                  const LHCb::UniqueIDGenerator& ),
                                                 LHCb::DetDesc::usesConditions<DetectorElement>> {
  public:
    using base_t =
        LHCb::Algorithm::MultiTransformer<OutputTrackType( InputTrackType const&, const Hits<Types>&...,
                                                           DetectorElement const&, const LHCb::UniqueIDGenerator& ),
                                          LHCb::DetDesc::usesConditions<DetectorElement>>;

    using base_t::info;
    using KeyValue = typename base_t::KeyValue;

    KalmanFilter( std::string const& name, ISvcLocator* pSvcLocator )
        : base_t( name, pSvcLocator,
                  {KeyValue( "Input", "" ), KeyValue( getKey<Types>(), "" )...,
                   KeyValue( "StandardGeometryTop", LHCb::standard_geometry_top ),
                   KeyValue( "InputUniqueIDGenerator", LHCb::UniqueIDGeneratorLocation::Default )},
                  LHCb::Algorithm::IOHelper<InputTrackType, OutputTrackType, base_t>::OutputNames(
                      get_out_names<OutputTrackType>() ) ) {}

    using tag        = LHCb::Pr::tag_type_t<InputTrackType>;
    using proxy_type = decltype( *( std::declval<InputTrackType const>().scalar().begin() ) );

    OutputTrackType operator()( InputTrackType const& tracks, const Hits<Types>&..., DetectorElement const& lhcb,
                                LHCb::UniqueIDGenerator const& unique_id_gen ) const override;

  private:
    Gaudi::Property<bool>   m_classic_smoothing_post{this, "ClassicSmoothing", true,
                                                   "Run classical smoother as post processing step for alignment."};
    Gaudi::Property<double> m_errorX{this, "ErrorX", 20.0 * Gaudi::Units::mm, "Seed error on x"};
    Gaudi::Property<double> m_errorY{this, "ErrorY", 20.0 * Gaudi::Units::mm, "Seed error on y"};
    Gaudi::Property<double> m_errorTx{this, "ErrorTx", 0.1, "Seed error on slope x"};
    Gaudi::Property<double> m_errorTy{this, "ErrorTy", 0.1, "Seed error on slope y"};
    Gaudi::Property<double> m_errorQoP{this, "ErrorQoP", 0.01, "Seed error on QoP"};
    Gaudi::Property<bool>   m_fill_fitresult{this, "FillFitResult", false, "Fill PrKalmanFitResult"};
    Gaudi::Property<double> m_maxchi2perdof{this, "MaxChi2", 9999999, "Maximum Chi2 per DoF"};
    Gaudi::Property<double> m_min_chi2_for_outlier{this, "MinChi2Outlier", 9,
                                                   "Minimum Chi2 of a node to be considered for outlier removal"};
    Gaudi::Property<double> m_maxchi2perdof_pre_outlier{this, "MaxChi2PreOutlierRemoval", 9999999,
                                                        "Maximum Chi2 per DoF before outlier removal"};
    Gaudi::Property<int>    m_max_fit_iter{this, "MaxFitIterations", 10, "max number of fit iterations to perform"};
    Gaudi::Property<int> m_max_outlier_iter{this, "MaxOutlierIterations", 2, "max number of fit iterations to perform"};
    Gaudi::Property<double> m_ptVelo{this, "VeloTrackPT", 400, "PT to use when fitting VELO tracks"};
    Gaudi::Property<size_t> m_minNumVPLayers{this, "MinNumVPHitsForOutlierRemoval", 3, "Minimum number of VP layers"};
    Gaudi::Property<size_t> m_minNumUTLayers{this, "MinNumUTHitsForOutlierRemoval", 3, "Minimum number of UT layers"};
    Gaudi::Property<size_t> m_minNumFTLayers{this, "MinNumFTHitsForOutlierRemoval", 6, "Minimum number of FT layers"};
    Gaudi::Property<size_t> m_minNumMuonLayers{this, "MinNumMuonHitsForOutlierRemoval", 4,
                                               "Minimum number of Muon layers"};

    ToolHandle<ITrackExtrapolator> m_extrap{this, "ReferenceExtrapolator", "TrackMasterExtrapolator"};

    mutable Gaudi::Accumulators::SummingCounter<unsigned int>   m_counter_states_failed{this, "Add states failed"};
    mutable Gaudi::Accumulators::Counter<>                      m_counter_cut{this, "chi2 cut"};
    mutable Gaudi::Accumulators::AveragingCounter<unsigned int> m_counter_iterations{this, "nIterations"};
    mutable Gaudi::Accumulators::AveragingCounter<unsigned int> m_counter_outlier_iterations{this,
                                                                                             "nOutlierIterations"};
    mutable Gaudi::Accumulators::SummingCounter<unsigned int>   m_counter_tracks_in{this, "nTracksInput"};
    mutable Gaudi::Accumulators::SummingCounter<unsigned int>   m_counter_tracks_out{this, "nTracksOutput"};
    mutable Gaudi::Accumulators::Counter<>                      m_counter_pre_outlier_cut{this, "Pre outlier chi2 cut"};
    mutable Gaudi::Accumulators::SummingCounter<unsigned int>   m_counter_transport_failed{this, "Transport failed"};
  };

  template <typename InputTrackType, typename OutputTrackType, HitType... Types>
  OutputTrackType KalmanFilter<InputTrackType, OutputTrackType, Types...>::
                  operator()( InputTrackType const& tracks, const Hits<Types>&... hits, DetectorElement const& lhcb,
              LHCb::UniqueIDGenerator const& unique_id_gen ) const {

    m_counter_tracks_in += tracks.size();
    const auto* extrap = m_extrap.get();

    auto outlier_iter_buffer{m_counter_outlier_iterations.buffer()};
    auto iter_buffer{m_counter_iterations.buffer()};
    auto transport_failed_buffer{m_counter_transport_failed.buffer()};
    auto states_failed_buffer{m_counter_states_failed.buffer()};
    auto pre_outlier_chi2_cut_buffer{m_counter_pre_outlier_cut.buffer()};
    auto chi2_cut_buffer{m_counter_cut.buffer()};

    const auto fit_config =
        KF::FitConfiguration{m_errorX * m_errorX,    m_errorY * m_errorY,     m_errorTx * m_errorTx,
                             m_errorTy * m_errorTy,  m_errorQoP * m_errorQoP, m_maxchi2perdof_pre_outlier,
                             m_min_chi2_for_outlier, m_maxchi2perdof,         m_max_outlier_iter,
                             m_max_fit_iter,         m_minNumVPLayers,        m_minNumUTLayers,
                             m_minNumFTLayers,       m_minNumMuonLayers};

    bool backward = [&] {
      if constexpr ( isVelo<InputTrackType> )
        return tracks.backward();
      else
        return false;
    }();

    auto output = initialize_output<InputTrackType, OutputTrackType>( tracks.size(), unique_id_gen, backward );

    std::vector<Node> fitnodes;
    fitnodes.reserve( 50 );

    auto const input_tracks = tracks.scalar();

    for ( auto const track : input_tracks ) {
      mydebug( "intput track", track.offset() );

      sort_hits_prepare_fitnodes<InputTrackType>( fitnodes, track, std::forward_as_tuple( hits... ) );

      auto const ctb_state_vec = get_upstream_seed( tracks, track, m_ptVelo );
      // To get the downstream seed for velo backward tracks
      // we use ctb state to get tx, ty and qop
      // and the most downstream node to get x, y, z.
      // Note that for velo backward tracks fitnodes
      // are ordered with increasing z.
      auto const downstream_state_vec = get_downstream_seed( tracks, track, m_ptVelo, ctb_state_vec, fitnodes.front() );

      const auto scatteringMomentum =
          std::clamp( std::abs( 1. / downstream_state_vec.qOverP() ), scatter_min, scatter_max );
      if ( init_nodes( fitnodes, scatteringMomentum, ctb_state_vec, downstream_state_vec, *lhcb.geometry(), *extrap,
                       isLong<InputTrackType> || isDownstream<InputTrackType> )
               .isFailure() ) {
        continue;
      }

      auto [chi2, success, nIter] = KF::iterate_fit( fitnodes, fit_config, *lhcb.geometry(), *extrap, iter_buffer,
                                                     pre_outlier_chi2_cut_buffer, transport_failed_buffer );
      if ( !success ) { continue; }

      if ( !KF::remove_outliers( fitnodes, fit_config, chi2, outlier_iter_buffer, chi2_cut_buffer ) ) { continue; }

      if constexpr ( isV1Tracks<OutputTrackType> ) {
        auto new_track = make_output_track( tracks, track, fitnodes, chi2, nIter, scatteringMomentum, m_fill_fitresult,
                                            m_classic_smoothing_post, *lhcb.geometry(), *extrap );
        if ( !new_track ) {
          ++states_failed_buffer;
          continue;
        }
        std::get<TracksV1>( output ).add( new_track.release() );
      } else if constexpr ( isV3Tracks<OutputTrackType> ) {
        auto sc = add_output_v3_track( std::get<LHCb::Event::v3::Tracks>( output ), tracks, track, fitnodes, chi2,
                                       scatteringMomentum, *lhcb.geometry(), *extrap, unique_id_gen );
        if ( sc.isFailure() ) {
          ++states_failed_buffer;
          continue;
        }
      } else if constexpr ( isV3TracksExtra<OutputTrackType> || isV3TracksFull<OutputTrackType> ) {
        auto& new_tracks        = std::get<LHCb::Event::v3::Tracks>( output );
        auto& new_partial_chi2s = std::get<LHCb::Event::v3::Track::PartialChiSquareds>( output );
        auto  sc = add_output_v3_track( new_tracks, tracks, track, fitnodes, chi2, scatteringMomentum, *lhcb.geometry(),
                                       *extrap, unique_id_gen );
        if ( sc.isFailure() ) {
          ++states_failed_buffer;
          continue;
        }
        add_output_v3_partial_chi2( new_partial_chi2s, fitnodes, chi2 );
        if constexpr ( isV3TracksFull<OutputTrackType> ) {
          auto& new_fit_results = std::get<std::vector<LHCb::PrKalmanFitResult>>( output );
          new_fit_results.emplace_back(
              make_fit_result( fitnodes, nIter, scatteringMomentum, m_classic_smoothing_post ) );
        }
      }

    } // loop over tracks

    if constexpr ( isV1Tracks<OutputTrackType> ) {
      m_counter_tracks_out += std::get<TracksV1>( output ).size();
      return output;
    } else if constexpr ( isV3Tracks<OutputTrackType> ) {
      m_counter_tracks_out += std::get<LHCb::Event::v3::Tracks>( output ).size();
      return output;
    } else if constexpr ( isV3TracksExtra<OutputTrackType> ) {
      auto& [new_tracks, new_partial_chi2s] = output;
      assert( new_tracks.size() == new_partial_chi2s.size() );
      m_counter_tracks_out += new_tracks.size();
      return output;
    } else if constexpr ( isV3TracksFull<OutputTrackType> ) {
      auto& [new_tracks, new_partial_chi2s, new_fit_results] = output;
      assert( new_tracks.size() == new_partial_chi2s.size() );
      assert( new_tracks.size() == new_fit_results.size() );
      m_counter_tracks_out += new_tracks.size();
      return output;
    }
  }

  using PrKalmanFilter_V1            = KalmanFilter<Long::Tracks, V1Output, HitType::VP, HitType::UT, HitType::FT>;
  using PrKalmanFilter_noUT_V1       = KalmanFilter<Long::Tracks, V1Output, HitType::VP, HitType::FT>;
  using PrKalmanFilter_Downstream_V1 = KalmanFilter<Downstream::Tracks, V1Output, HitType::UT, HitType::FT>;
  using PrKalmanFilter_Seed_V1       = KalmanFilter<Seeding::Tracks, V1Output, HitType::FT>;
  using PrKalmanFilter_Velo_V1       = KalmanFilter<Velo::Tracks, V1Output, HitType::VP>;
  using PrKalmanFilter_Upstream_V1   = KalmanFilter<Upstream::Tracks, V1Output, HitType::VP, HitType::UT>;

  using PrKalmanFilter_V3            = KalmanFilter<Long::Tracks, V3Output, HitType::VP, HitType::UT, HitType::FT>;
  using PrKalmanFilter_noUT_V3       = KalmanFilter<Long::Tracks, V3Output, HitType::VP, HitType::FT>;
  using PrKalmanFilter_Downstream_V3 = KalmanFilter<Downstream::Tracks, V3Output, HitType::UT, HitType::FT>;
  using PrKalmanFilter_Seed_V3       = KalmanFilter<Seeding::Tracks, V3Output, HitType::FT>;
  using PrKalmanFilter_Velo_V3       = KalmanFilter<Velo::Tracks, V3Output, HitType::VP>;
  using PrKalmanFilter_Upstream_V3   = KalmanFilter<Upstream::Tracks, V3Output, HitType::VP, HitType::UT>;

  using PrKalmanFilter_V3Extra      = KalmanFilter<Long::Tracks, V3ExtraOutput, HitType::VP, HitType::UT, HitType::FT>;
  using PrKalmanFilter_noUT_V3Extra = KalmanFilter<Long::Tracks, V3ExtraOutput, HitType::VP, HitType::FT>;
  using PrKalmanFilter_Downstream_V3Extra = KalmanFilter<Downstream::Tracks, V3ExtraOutput, HitType::UT, HitType::FT>;
  using PrKalmanFilter_Seed_V3Extra       = KalmanFilter<Seeding::Tracks, V3ExtraOutput, HitType::FT>;
  using PrKalmanFilter_Velo_V3Extra       = KalmanFilter<Velo::Tracks, V3ExtraOutput, HitType::VP>;
  using PrKalmanFilter_Upstream_V3Extra   = KalmanFilter<Upstream::Tracks, V3ExtraOutput, HitType::VP, HitType::UT>;

  using PrKalmanFilter_V3Full      = KalmanFilter<Long::Tracks, V3FullOutput, HitType::VP, HitType::UT, HitType::FT>;
  using PrKalmanFilter_noUT_V3Full = KalmanFilter<Long::Tracks, V3FullOutput, HitType::VP, HitType::FT>;
  using PrKalmanFilter_Downstream_V3Full = KalmanFilter<Downstream::Tracks, V3FullOutput, HitType::UT, HitType::FT>;
  using PrKalmanFilter_Seed_V3Full       = KalmanFilter<Seeding::Tracks, V3FullOutput, HitType::FT>;
  using PrKalmanFilter_Velo_V3Full       = KalmanFilter<Velo::Tracks, V3FullOutput, HitType::VP>;
  using PrKalmanFilter_Upstream_V3Full   = KalmanFilter<Upstream::Tracks, V3FullOutput, HitType::VP, HitType::UT>;

  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_V1, "PrKalmanFilter" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_noUT_V1, "PrKalmanFilter_noUT" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Downstream_V1, "PrKalmanFilter_Downstream" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Seed_V1, "PrKalmanFilter_Seed" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Velo_V1, "PrKalmanFilter_Velo" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Upstream_V1, "PrKalmanFilter_Upstream" )

  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_V3, "PrKalmanFilter_V3" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_noUT_V3, "PrKalmanFilter_noUT_V3" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Downstream_V3, "PrKalmanFilter_Downstream_V3" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Seed_V3, "PrKalmanFilter_Seed_V3" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Velo_V3, "PrKalmanFilter_Velo_V3" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Upstream_V3, "PrKalmanFilter_Upstream_V3" )

  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_V3Extra, "PrKalmanFilter_V3Extra" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_noUT_V3Extra, "PrKalmanFilter_noUT_V3Extra" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Downstream_V3Extra, "PrKalmanFilter_Downstream_V3Extra" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Seed_V3Extra, "PrKalmanFilter_Seed_V3Extra" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Velo_V3Extra, "PrKalmanFilter_Velo_V3Extra" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Upstream_V3Extra, "PrKalmanFilter_Upstream_V3Extra" )

  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_V3Full, "PrKalmanFilter_V3Full" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_noUT_V3Full, "PrKalmanFilter_noUT_V3Full" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Downstream_V3Full, "PrKalmanFilter_Downstream_V3Full" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Seed_V3Full, "PrKalmanFilter_Seed_V3Full" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Velo_V3Full, "PrKalmanFilter_Velo_V3Full" )
  DECLARE_COMPONENT_WITH_ID( PrKalmanFilter_Upstream_V3Full, "PrKalmanFilter_Upstream_V3Full" )
} // namespace LHCb::Pr
