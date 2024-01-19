/*****************************************************************************\
* (c) Copyright 2020 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once
// Rec
#include "Event/ParametrisedScatters.h"
#include "Event/PrFitNode.h"
#include "Event/PrKalmanFitResult.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
// LHCb
#include "Event/ChiSquare.h"
#include "Event/PrDownstreamTracks.h"
#include "Event/PrHits.h"
#include "Event/PrLongTracks.h"
#include "Event/PrSeedTracks.h"
#include "Event/PrVeloTracks.h"
#include "Event/TrackTypes.h"
#include "Kernel/HitPattern.h"
#include "Kernel/TransformedRange.h"
#include "LHCbMath/Similarity.h"

#include "Event/PartialChiSquareds.h"
#include "Event/Track_v3.h"

// Gaudi
#include "GaudiKernel/Kernel.h"
// std
#include <boost/container/static_vector.hpp>
#include <cstddef>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

namespace LHCb::Pr::Tracks::Fit {

  constexpr inline auto scatter_min = 100. * Gaudi::Units::MeV;
  constexpr inline auto scatter_max = 500000. * Gaudi::Units::MeV;

  constexpr inline auto ZBegVelo = -400. * Gaudi::Units::mm;

  // in PR weight = 12.0 / ( pitch * pitch );
  // Measurement Provider error: pitch * sqrt( eValue * eValue + m_APE * m_APE );
  // evalue is .22, .14, or .25 based on pseudosize.
  // ???APE is currently set to zero but is a property that could change?
  template <typename F>
  auto ut_err( F weight ) {
    return ( 1. / sqrt( weight ) ).cast();
  }
  // FIXME Pr gets weight from an array via clustersize
  // MeasurementProvider does 0.04 + 0.01 * itH->pseudoSize()
  // we don't know the pseudosize anymore at this stage though :(
  inline auto           ft_err( double weight ) { return .5 / std::sqrt( weight ); }
  constexpr inline auto vp_err = 0.0125;

  template <HitType Type, HitType... Types>
  inline constexpr bool has_hit_type = ( is_equal<Type, Types> || ... );

  template <typename T>
  inline constexpr bool isLong = std::is_same_v<T, Long::Tracks>;

  template <typename T>
  inline constexpr bool isDownstream = std::is_same_v<T, Downstream::Tracks>;

  template <typename T>
  inline constexpr bool isUpstream = std::is_same_v<T, Upstream::Tracks>;

  template <typename T>
  inline constexpr bool isSeed = std::is_same_v<T, Seeding::Tracks>;

  template <typename T>
  inline constexpr bool isVelo = std::is_same_v<T, Velo::Tracks>;

  using dType = SIMDWrapper::scalar::types;
  using I     = dType::int_v;
  using F     = dType::float_v;

  using TracksV1      = LHCb::Event::v1::Tracks;
  using V1Output      = std::tuple<LHCb::Event::v1::Tracks>;
  using V3Output      = std::tuple<LHCb::Event::v3::Tracks>;
  using V3ExtraOutput = std::tuple<LHCb::Event::v3::Tracks, LHCb::Event::v3::Track::PartialChiSquareds>;
  using V3FullOutput  = std::tuple<LHCb::Event::v3::Tracks, LHCb::Event::v3::Track::PartialChiSquareds,
                                  std::vector<LHCb::PrKalmanFitResult>>;

  template <typename T>
  auto get_out_names() {
    if constexpr ( std::is_same_v<T, V1Output> )
      return std::array<std::string, 1>{"OutputTracks"};
    else if constexpr ( std::is_same_v<T, V3Output> )
      return std::array<std::string, 1>{"OutputTracks"};
    else if constexpr ( std::is_same_v<T, V3ExtraOutput> )
      return std::array<std::string, 2>{"OutputTracks", "OutputPartialChi2s"};
    else if constexpr ( std::is_same_v<T, V3FullOutput> )
      return std::array<std::string, 3>{"OutputTracks", "OutputPartialChi2s", "OutputTrackFitResults"};
    else
      throw GaudiException( "Output type is not supported", "PrKalmanFilter", StatusCode::FAILURE );
  }

  template <typename T>
  inline constexpr bool isV1Tracks = std::is_same_v<T, V1Output>;

  template <typename T>
  inline constexpr bool isV3Tracks = std::is_same_v<T, V3Output>;

  template <typename T>
  inline constexpr bool isV3TracksExtra = std::is_same_v<T, V3ExtraOutput>;

  template <typename T>
  inline constexpr bool isV3TracksFull = std::is_same_v<T, V3FullOutput>;

  template <typename T>
  constexpr inline auto v3_track_type() {
    if ( isLong<T> )
      return LHCb::Event::v3::TrackType::Long;
    else if ( isDownstream<T> )
      return LHCb::Event::v3::TrackType::Downstream;
    else if ( isVelo<T> )
      return LHCb::Event::v3::TrackType::Velo;
    else if ( isSeed<T> )
      return LHCb::Event::v3::TrackType::Ttrack;
    else if ( isUpstream<T> )
      return LHCb::Event::v3::TrackType::Upstream;
  };

  template <typename T>
  inline constexpr bool hasFT = isLong<T> || isDownstream<T> || isSeed<T>;

  template <typename T>
  inline constexpr bool hasUT = isLong<T> || isDownstream<T> || isUpstream<T>;

  template <typename T>
  inline constexpr bool hasVP = isVelo<T> || isLong<T> || isUpstream<T>;

  template <typename InputTrackType>
  using proxy_type = decltype( *( std::declval<InputTrackType const>().scalar().begin() ) );

  template <bool enable = false, typename... Args>
  inline void mydebug( Args&&... args ) {
    if constexpr ( enable ) { ( ( std::cout << std::forward<Args>( args ) << " " ), ... ) << std::endl; }
  }

  template <typename T>
  void apply_delta_E( T& state, double dE ) {
    // apply a momentum loss/gain to a state
    if ( std::abs( state.qOverP() ) > LHCb::Math::lowTolerance ) {
      auto charge = state.qOverP() > 0 ? 1. : -1.;
      auto momnew = std::max( 10., std::abs( 1 / state.qOverP() ) + dE );
      if ( std::abs( momnew ) > 10 ) state.setQOverP( charge / momnew );
    }
  }

  template <typename InputTrackType, HitType... Types>
  void sort_hits_prepare_fitnodes( std::vector<Node>& fitnodes, const proxy_type<InputTrackType>& track,
                                   std::tuple<const Hits<Types>&...> containers ) {

    static_assert( TracksInfo::MaxFTHits >= TracksInfo::MaxUTHits );
    static_assert( isLong<InputTrackType> || isVelo<InputTrackType> || isDownstream<InputTrackType> ||
                       isUpstream<InputTrackType> || isSeed<InputTrackType>,
                   "Only Velo, Downstream, Long or Seed tracks are supported currently." );

    boost::container::static_vector<int, TracksInfo::MaxFTHits> sort_order{};
    fitnodes.clear();
    // adding ft hits to fitnodes
    if constexpr ( has_hit_type<HitType::FT, Types...> ) {
      static_assert( isLong<InputTrackType> || isDownstream<InputTrackType> || isSeed<InputTrackType> );

      const auto& hits_ft = get<HitType::FT>( containers );
      sort_order.resize( track.nFTHits().cast() );

      // get the indices into the hit container for all our ft hits
      std::generate( sort_order.begin(), sort_order.end(),
                     [&, cnt = 0]() mutable { return track.ft_index( cnt++ ).cast(); } );

      // sort the indices according to the hit z-coordinate
      std::sort( sort_order.begin(), sort_order.end(),
                 [&]( int const i, int const j ) { return hits_ft.z( i ) > hits_ft.z( j ); } );

      for ( auto const idx : sort_order ) {
        // X and Z are defined at y=0
        const auto x      = hits_ft.x( idx );
        const auto z      = hits_ft.z( idx );
        const auto w      = hits_ft.w( idx );
        const auto dxdy   = hits_ft.dxDy( idx );
        const auto dzdy   = hits_ft.dzDy( idx );
        const auto lhcbid = hits_ft.lhcbid( idx );
        fitnodes.emplace_back( x, 0, z, dxdy, 1, dzdy, ft_err( w ), Node::Type::FTHit, lhcbid );
      }
    }

    // adding ut hits to fitnodes
    if constexpr ( has_hit_type<HitType::UT, Types...> ) {
      static_assert( isLong<InputTrackType> || isDownstream<InputTrackType> || isUpstream<InputTrackType> );

      const auto& hits_ut = get<HitType::UT>( containers );
      sort_order.resize( track.nUTHits().cast() );

      // get the indices into the hit container for all our ut hits
      std::generate( sort_order.begin(), sort_order.end(),
                     [&, cnt = 0]() mutable { return track.ut_index( cnt++ ).cast(); } );

      const auto hits_vec = hits_ut.scalar();
      // sort the indices according to the hit z-coordinate
      std::sort( sort_order.begin(), sort_order.end(), [&hits_vec]( int const i, int const j ) {
        return hits_vec[i].template get<UT::UTHitsTag::zAtYEq0>().cast() >
               hits_vec[j].template get<UT::UTHitsTag::zAtYEq0>().cast();
      } );

      for ( auto const idx : sort_order ) {

        const auto x0     = hits_vec[idx].template get<UT::UTHitsTag::xAtYEq0>();
        const auto z0     = hits_vec[idx].template get<UT::UTHitsTag::zAtYEq0>();
        const auto w      = hits_vec[idx].template get<UT::UTHitsTag::weight>();
        const auto dxdy   = hits_vec[idx].template get<UT::UTHitsTag::dxDy>();
        const auto id     = hits_vec[idx].template get<UT::UTHitsTag::channelID>();
        const auto ut_id  = LHCb::Detector::UT::ChannelID{bit_cast<unsigned>( id.cast() )};
        const auto lhcbid = LHCb::LHCbID{ut_id};
        fitnodes.emplace_back( x0.cast(), 0, z0.cast(), dxdy.cast(), 1, 0, ut_err( w ), Node::Type::UTHit, lhcbid );
      }
    }

    // adding vp nodes if we have a track that has vp hits
    if constexpr ( has_hit_type<HitType::VP, Types...> ) {
      static_assert( isLong<InputTrackType> || isVelo<InputTrackType> || isUpstream<InputTrackType> );

      const auto& hits_vp = get<HitType::VP>( containers );
      const auto  n_hits  = track.nVPHits().cast();
      // note velo is by default ordered from higher to lower z-coordinate
      for ( int j{0}; j < n_hits; ++j ) {
        const auto idx    = track.vp_index( j ).cast();
        const auto vphit  = hits_vp.scalar();
        const auto pos    = vphit[idx].template get<VP::VPHitsTag::pos>().vec3();
        const auto id     = vphit[idx].template get<VP::VPHitsTag::ChannelId>();
        const auto vp_id  = LHCb::Detector::VPChannelID{bit_cast<unsigned>( id.cast() )};
        const auto lhcbid = LHCb::LHCbID{vp_id};

        fitnodes.emplace_back( pos.x.cast(), pos.y.cast(), pos.z.cast(), 0, 1, 0, vp_err, Node::Type::VPHit, lhcbid );
        fitnodes.emplace_back( pos.x.cast(), pos.y.cast(), pos.z.cast(), 1, 0, 0, vp_err, Node::Type::VPHit, lhcbid );
      }
      // For backward velo tracks we should reverse the order of nodes,
      // then the rest of fitter logic can stay the same.
      // Note that the order of each VP node pair does not matter,
      // see https://gitlab.cern.ch/lhcb/Rec/-/merge_requests/2734
      // for resolution plots.
      if constexpr ( isVelo<InputTrackType> )
        if ( track.backward() ) std::reverse( fitnodes.begin(), fitnodes.end() );
    }
  }

  template <typename InputTrackType>
  auto get_upstream_seed( const InputTrackType& tracks, const proxy_type<InputTrackType>& track, float ptVelo ) {

    if constexpr ( isLong<InputTrackType> ) {

      auto const vp_tracks = tracks.getVeloAncestors()->scalar();
      auto const qop       = track.qOverP().cast();
      auto const vp_tr     = track.trackVP().cast();
      auto const vp_track  = vp_tracks[vp_tr];
      auto const pos       = vp_track.StatePos( 0 );
      auto const dir       = vp_track.StateDir( 0 );

      return LHCb::StateVector{{pos.x.cast(), pos.y.cast(), dir.x.cast(), dir.y.cast(), qop}, pos.z.cast()};

    } else if constexpr ( isDownstream<InputTrackType> ) {

      auto const qop = track.qOverP().cast();
      auto const s   = track.template get<LHCb::Pr::Downstream::Tag::State>();

      return LHCb::StateVector{{s.x().cast(), s.y().cast(), s.tx().cast(), s.ty().cast(), qop}, s.z().cast()};

    } else if constexpr ( isSeed<InputTrackType> ) {

      const auto qop = track.qOverP().cast();
      // get first state (this is also where qop is defined)
      const auto s = track.template get<Seeding::Tag::States>( 0 );
      return LHCb::StateVector{{s.x().cast(), s.y().cast(), s.tx().cast(), s.ty().cast(), qop}, s.z().cast()};

    } else if constexpr ( isVelo<InputTrackType> ) {

      using tag         = LHCb::Pr::Velo::Tag;
      auto const x      = track.template get<tag::States>( 0 ).x().cast();
      auto const y      = track.template get<tag::States>( 0 ).y().cast();
      auto const z      = track.template get<tag::States>( 0 ).z().cast();
      auto const tx     = track.template get<tag::States>( 0 ).tx().cast();
      auto const ty     = track.template get<tag::States>( 0 ).ty().cast();
      auto const slope2 = std::max( tx * tx + ty * ty, float{1e-20} );
      auto const charge = 1.f;
      auto const qop    = charge * sqrt( slope2 ) / ( ptVelo * sqrt( 1. + slope2 ) );

      return LHCb::StateVector{{x, y, tx, ty, qop}, z};

    } else if constexpr ( isUpstream<InputTrackType> ) {

      auto const vp_tracks = tracks.getVeloAncestors()->scalar();
      auto const qop       = track.qOverP().cast();
      auto const vp_tr     = track.trackVP().cast();
      auto const vp_track  = vp_tracks[vp_tr];
      auto const pos       = vp_track.StatePos( 0 );
      auto const dir       = vp_track.StateDir( 0 );

      return LHCb::StateVector{{pos.x.cast(), pos.y.cast(), dir.x.cast(), dir.y.cast(), qop}, pos.z.cast()};
    }
  }

  template <typename InputTrackType>
  auto get_downstream_seed( const InputTrackType& tracks, const proxy_type<InputTrackType>& track, float ptVelo,
                            const LHCb::StateVector& ctb_state, const Node& downstream_node ) {

    if constexpr ( isLong<InputTrackType> ) {
      using tag = LHCb::Pr::Long::Tag;

      auto const qop = track.qOverP().cast();
      auto const x   = track.template get<tag::States>( 1 ).x().cast();
      auto const y   = track.template get<tag::States>( 1 ).y().cast();
      auto const z   = track.template get<tag::States>( 1 ).z().cast();
      auto const tx  = track.template get<tag::States>( 1 ).tx().cast();
      auto const ty  = track.template get<tag::States>( 1 ).ty().cast();

      return LHCb::StateVector{{x, y, tx, ty, qop}, z};

    } else if constexpr ( isDownstream<InputTrackType> ) {
      auto const seed_tracks = tracks.getFTAncestors()->scalar();
      auto const ft_tr       = track.template get<Downstream::Tag::trackSeed>().cast();
      auto const qop         = track.qOverP().cast();
      auto const ft_track    = seed_tracks[ft_tr];
      auto const s           = ft_track.template get<Seeding::Tag::States>( TracksInfo::NumSeedStates - 1 );

      return LHCb::StateVector{{s.x().cast(), s.y().cast(), s.tx().cast(), s.ty().cast(), qop}, s.z().cast()};

    } else if constexpr ( isSeed<InputTrackType> ) {
      const auto qop = track.qOverP().cast();
      // get last state
      const auto s = track.template get<Seeding::Tag::States>( TracksInfo::NumSeedStates - 1 );
      return LHCb::StateVector{{s.x().cast(), s.y().cast(), s.tx().cast(), s.ty().cast(), qop}, s.z().cast()};

    } else if constexpr ( isVelo<InputTrackType> ) {

      if ( !track.backward() ) {
        using tag         = LHCb::Pr::Velo::Tag;
        auto const x      = track.template get<tag::States>( 1 ).x().cast();
        auto const y      = track.template get<tag::States>( 1 ).y().cast();
        auto const z      = track.template get<tag::States>( 1 ).z().cast();
        auto const tx     = track.template get<tag::States>( 1 ).tx().cast();
        auto const ty     = track.template get<tag::States>( 1 ).ty().cast();
        auto const slope2 = std::max( tx * tx + ty * ty, float{1e-20} );
        auto const charge = 1.f;
        auto const qop    = charge * sqrt( slope2 ) / ( ptVelo * sqrt( 1. + slope2 ) );

        return LHCb::StateVector{{x, y, tx, ty, qop}, z};
      } else {
        auto const x   = downstream_node.measurement_pos[0];
        auto const y   = downstream_node.measurement_pos[1];
        auto const z   = downstream_node.measurement_pos[2];
        auto const tx  = ctb_state.tx();
        auto const ty  = ctb_state.ty();
        auto const qop = ctb_state.qOverP();

        return LHCb::StateVector{{x, y, tx, ty, qop}, z};
      }

    } else if constexpr ( isUpstream<InputTrackType> ) {
      auto const qop = track.qOverP().cast();
      auto const s   = track.template get<LHCb::Pr::Upstream::Tag::State>();

      return LHCb::StateVector{{s.x().cast(), s.y().cast(), s.tx().cast(), s.ty().cast(), qop}, s.z().cast()};
    }
  }

  inline StatusCode init_nodes( LHCb::span<Node> fitnodes, double scatteringMomentum,
                                const LHCb::StateVector& ctb_state_vec, const LHCb::StateVector& last_measurement,
                                const IGeometryInfo& geometry, const ITrackExtrapolator& extrap,
                                bool track_through_magnet ) {

    mydebug( "Start of init_nodes" );
    // FIXME: this could use some rethinking
    // Long tracks are currently initialized from CTB and AtT states but
    // we would probably do better using CTB, something in UT and something in the middle of SciFi
    // This however requires more logic than the simple extrapolate to first and last hit and then go step
    // by step that is done below
    // FIXME: furhtermore this logic isn't really ideal for downstream where the upstream seed is
    // in the middel of the UT, thus we extrapolate that to the first UT it and then go hit by hit
    // instead of going lef and right from the ut state.
    Gaudi::TrackMatrix F{ROOT::Math::SMatrixIdentity()};

    // we need want to init the first node based on scifi state, otherwise ctb state
    StatusCode ret;
    {
      auto start_vec = last_measurement;
      ret            = extrap.propagate( start_vec, fitnodes.front().z(), geometry );
      if ( ret.isFailure() ) { return ret; }
      // set reference state
      fitnodes.front().lhcb_ref_vector.parameters() = start_vec.parameters();
    }

    {
      auto end_vec = ctb_state_vec;
      ret          = extrap.propagate( end_vec, fitnodes.back().z(), geometry );
      if ( ret.isFailure() ) { return ret; }
      fitnodes.back().lhcb_ref_vector.parameters() = end_vec.parameters();
    }

    // if we are only processing velo tracks, we don't need to do
    // the step of first setting ref vectors for upstream of magnet nodes

    auto const is_on_other_side_of_mag = [track_through_magnet]( auto const& node ) {
      return ( track_through_magnet && node->z() < 4000 );
    };

    // init all nodes on the velo side of magnet
    // if we aren't fitting velo tracks (see lambda above)
    for ( auto node = std::next( fitnodes.rbegin() ); is_on_other_side_of_mag( node ); ++node ) {
      auto state_vec = std::prev( node )->lhcb_ref_vector;
      ret            = extrap.propagate( state_vec, node->z(), geometry );
      if ( ret.isFailure() ) { return ret; }
      node->lhcb_ref_vector.parameters() = state_vec.parameters();
    }

    auto lhcb_ref_vector{std::cref( fitnodes.begin()->lhcb_ref_vector )};
    // now we are going to go from scifi -> velo
    // and also create the transport matrices and vetors
    for ( auto node = std::next( fitnodes.begin() ); node != fitnodes.end(); ++node ) {
      // wanted explicit copy
      auto       state_vec = lhcb_ref_vector.get();
      auto const z         = node->z();
      ret                  = extrap.propagate( state_vec, z, geometry, &F );
      if ( ret.isFailure() ) { return ret; }

      if ( not is_on_other_side_of_mag( node ) ) node->lhcb_ref_vector.parameters() = state_vec.parameters();

      std::tie( node->noise_matrix, node->delta_energy ) =
          TrackFit::param_scatter_impl::computeNoiseAndDeltaE( *( node - 1 ), *node, scatteringMomentum );

      apply_delta_E( state_vec, node->delta_energy );

      Gaudi::TrackVector transportvec = state_vec.parameters() - F * lhcb_ref_vector.get().parameters();
      node->set_transport( F, transportvec );
      // update the reference
      lhcb_ref_vector = std::cref( node->lhcb_ref_vector );
    }
    mydebug( "End of init_nodes" );
    if ( ret.isSuccess() )
      for ( auto& node : fitnodes ) { node.project_reference(); }
    return ret;
  }

  namespace KF {
    constexpr auto nan = std::numeric_limits<double>::signaling_NaN();
    struct FitConfiguration {
      double varianceX{nan};
      double varianceY{nan};
      double varianceTx{nan};
      double varianceTy{nan};
      double varianceQoP{nan};
      double maxchi2perdof_pre_outlier{nan};
      double min_chi2_for_outlier{nan};
      double maxchi2perdof{nan};
      int    max_outlier_iter{0};
      int    max_fit_iter{0};
      size_t minNumVPLayers{3};
      size_t minNumUTLayers{3};
      size_t minNumFTLayers{6};
      size_t minNumMuonLayers{4};
    };

    inline void pre_fit_init( LHCb::span<Node> fitnodes, const FitConfiguration& fit_config ) {
      auto& seed_cov   = fitnodes.back().predicted_state_cov[Node::backward];
      seed_cov( 0, 0 ) = fit_config.varianceX;
      seed_cov( 1, 1 ) = fit_config.varianceY;
      seed_cov( 2, 2 ) = fit_config.varianceTx;
      seed_cov( 3, 3 ) = fit_config.varianceTy;
      seed_cov( 4, 4 ) = fit_config.varianceQoP;

      fitnodes.front().predicted_state_cov[Node::forward] = fitnodes.back().predicted_state_cov[Node::backward];
    }

    template <int direction = Node::forward>
    void predict( Node const& prevnode, Node& node ) {
      mydebug( "predict of node at z= ", node.z() );

      if constexpr ( direction == Node::forward ) {
        auto const& F = node.transport_matrix;
        mydebug( "transport Matrix F\n", F );
        mydebug( "transport Vector\n", node.transport_vector );
        node.predicted_state_vec[direction] = F * prevnode.filtered_state_vec[direction] + node.transport_vector;
        LHCb::Math::Similarity( F, prevnode.filtered_state_cov[direction], node.predicted_state_cov[direction] );
        node.predicted_state_cov[direction] += node.noise_matrix;
      } else {
        auto const& F_inv = prevnode.transport_matrix_inverse;
        mydebug( "inverse prevnode transport Matrix F\n", F_inv );
        mydebug( "prevnode transport Vector\n", prevnode.transport_vector );
        node.predicted_state_vec[direction] =
            F_inv * ( prevnode.filtered_state_vec[direction] - prevnode.transport_vector );

        LHCb::Math::Similarity( F_inv, prevnode.filtered_state_cov[direction] + prevnode.noise_matrix,
                                node.predicted_state_cov[direction] );
      }

      mydebug( "predict statevec:  ", node.predicted_state_vec[direction] );
      mydebug( "predict statecov:\n", node.predicted_state_cov[direction] );
    }

    template <int direction = Node::forward>
    LHCb::ChiSquare filter( Node& node ) {
      mydebug( "filter of node at z= ", node.z() );

      node.filtered_state_vec[direction] = node.predicted_state_vec[direction];
      node.filtered_state_cov[direction] = node.predicted_state_cov[direction];

      if ( node.m_is_outlier ) {
        node.delta_chi2[direction] = LHCb::ChiSquare{};
        return node.delta_chi2[direction];
      }

      // EKF res = residual + H * (x_ref - x)
      mydebug( "RefVector:  ", node.ref_vector() );
      mydebug( "residual:  ", node.ref_residual );
      mydebug( "error:  ", node.measurement_error );
      mydebug( "projection:  ", node.projection );
      double chi2 =
          LHCb::Math::Filter( node.filtered_state_vec[direction], node.filtered_state_cov[direction], node.ref_vector(),
                              node.projection, node.ref_residual, node.measurement_error * node.measurement_error );

      node.delta_chi2[direction] = LHCb::ChiSquare( chi2, 1 );

      mydebug( "filter statevec:  ", node.filtered_state_vec[direction] );
      mydebug( "filter statecov:\n", node.filtered_state_cov[direction] );
      mydebug( "delta chi2: ", node.delta_chi2[direction] );
      return node.delta_chi2[direction];
    }

    inline void average_node( Node& node ) {
      // Reminder: This function doesn't care if a node is an outlier
      // If it is, the filtering step was skipped for that node. Thus we end up
      // simply averaging two predictions without taking into account the measurement
      mydebug( "average of node at z= ", node.z() );

      int const filter_idx =
          ( node.predicted_state_cov[Node::backward]( 0, 0 ) > node.predicted_state_cov[Node::forward]( 0, 0 ) )
              ? Node::backward
              : Node::forward;
      int const predict_idx = ( filter_idx + 1 ) % 2;

      LHCb::Math::Average( node.filtered_state_vec[filter_idx], node.filtered_state_cov[filter_idx],
                           node.predicted_state_vec[predict_idx], node.predicted_state_cov[predict_idx],
                           node.final_state_vec, node.final_state_cov );

      mydebug( "average filter_idx:", filter_idx, "predict_idx:", predict_idx );
      mydebug( "average statevec:  ", node.final_state_vec );
      mydebug( "average statecov:\n", node.final_state_cov, "\n" );
    }

    template <int direction = Node::forward>
    LHCb::ChiSquare predict_and_filter( Node const& prevnode, Node& node ) {
      predict<direction>( prevnode, node );
      return filter<direction>( node );
    }

    inline LHCb::ChiSquare fit( LHCb::span<Node> fitnodes ) {
      // not checking that first_hit_node_fwd != end or that I can increment the iterator safely.
      // But it should be safe since we know our track is made of a couple of hits
      // and I can't think of a scenario where this could crash.

      fitnodes.front().predicted_state_vec[Node::forward] = fitnodes.front().ref_vector();
      auto chi2_fwd = filter<Node::forward>( fitnodes.front() ) + LHCb::ChiSquare{0, -5};
      for ( size_t i{1}; i < fitnodes.size(); ++i ) {
        chi2_fwd += predict_and_filter<Node::forward>( fitnodes[i - 1], fitnodes[i] );
      }

      // start from forward filter
      fitnodes.back().predicted_state_vec[Node::backward] = fitnodes.back().filtered_state_vec[Node::forward];
      // start independent from ref_vector
      // fitnodes.back().predicted_state_vec[FitNode::backward] = fitnodes.back().ref_vector();
      auto chi2_bkwd = filter<Node::backward>( fitnodes.back() ) + LHCb::ChiSquare{0, -5};
      for ( int i( fitnodes.size() - 2 ); i >= 0; --i ) {
        chi2_bkwd += predict_and_filter<Node::backward>( fitnodes[i + 1], fitnodes[i] );
      }

      // use the smaller chi2
      mydebug( chi2_fwd, chi2_bkwd );
      return chi2_fwd.chi2() < chi2_bkwd.chi2() ? chi2_fwd : chi2_bkwd;
    }

    inline void smooth_and_update_ref_vector( LHCb::span<Node> fitnodes ) {
      fitnodes.front().ref_vector() = fitnodes.front().filtered_state_vec[Node::backward];
      fitnodes.back().ref_vector()  = fitnodes.back().filtered_state_vec[Node::forward];

      std::for_each( fitnodes.begin() + 1, fitnodes.end() - 1, []( Node& node ) {
        average_node( node );
        node.ref_vector() = node.final_state_vec;
      } );
    }

    inline void smooth_and_update_ref_vector_outlier( LHCb::span<Node> fitnodes ) {

      auto forward_iter  = fitnodes.begin();
      auto backward_iter = fitnodes.rbegin();

      // the first hit is never averaged and following ones only if the previous ones aren't all outliers
      do {
        forward_iter->ref_vector()    = forward_iter->filtered_state_vec[Node::backward];
        forward_iter->final_state_cov = forward_iter->filtered_state_cov[Node::backward];
        ++forward_iter;
      } while ( forward_iter->m_is_outlier );

      // same as above from the other end
      do {
        backward_iter->ref_vector()    = backward_iter->filtered_state_vec[Node::forward];
        backward_iter->final_state_cov = backward_iter->filtered_state_cov[Node::forward];
        ++backward_iter;
      } while ( backward_iter->m_is_outlier );

      std::for_each( fitnodes.begin() + std::distance( fitnodes.begin(), forward_iter ),
                     fitnodes.end() - std::distance( fitnodes.rbegin(), backward_iter ), []( Node& node ) {
                       average_node( node );
                       node.ref_vector() = node.final_state_vec;
                     } );
    }

    inline void classical_smoother_iteration_for_alignment( LHCb::span<Node>               fitnodes,
                                                            LHCb::span<Gaudi::TrackMatrix> gain_matrices ) {
      // classical smoother does NOT rely on any information from the backward pass.
      // Usually one would do one forward iteration and then run the classical smoother
      // thus any use of backward information in the below would we a bug!

      // WARNING
      // Assumes, fit procedure has already run and outlier removal is over.
      // should be seen as a posprocessing step for use of fitnodes in alignment
      // FIXME we can probably assert some of this

      // start at the last node in forward direction
      int idx = fitnodes.size() - 1;
      // we start smoothing at the first node that has upstream info
      do {
        auto& node           = fitnodes[idx];
        node.final_state_vec = node.filtered_state_vec[Node::forward];
        node.final_state_cov = node.filtered_state_cov[Node::forward];
        // FIXME should we set ref vector and update ref-residual?
        // I think so right? best info for alignment would be the updated info?
        node.ref_vector() = node.filtered_state_vec[Node::forward];
        // updates H and the ref residual based on the new ref_vector
        node.project_reference();
        //
        --idx;
      } while ( idx > 0 && fitnodes[idx].m_is_outlier );

      // all remaining nodes need to be smoothed
      for ( ; idx >= 0; --idx ) {

        auto& node = fitnodes[idx];
        // next node in forward direction
        auto const& nextnode = fitnodes[idx + 1];

        // Get the filtered result from this node
        node.final_state_vec = node.filtered_state_vec[Node::forward];
        node.final_state_cov = node.filtered_state_cov[Node::forward];

        auto const& F = nextnode.transport_matrix;
        auto const& Q = nextnode.noise_matrix;
        auto&       A = gain_matrices[idx];

        bool nonZeroNoise = ( Q( 2, 2 ) + Q( 3, 3 ) + Q( 4, 4 ) ) > 0;
        if ( nonZeroNoise ) {

          // invert the covariance matrix
          auto inv_nextnode_cov = nextnode.predicted_state_cov[Node::forward];
          if ( !inv_nextnode_cov.InvertChol() ) {
            // FIXME
            throw std::runtime_error( "inversion in classical smoother failed" );
          }

          // calculate gain matrix A. we can make this quicker by epxloiting that F is empty
          A                            = node.final_state_cov * Transpose( F ) * inv_nextnode_cov;
          Gaudi::TrackSymMatrix FCFinv = LHCb::Math::Similarity( F, node.final_state_cov ); // is also nextNodeC - Q
          FCFinv.InvertChol();
          Gaudi::TrackSymMatrix sum = nextnode.final_state_cov + Q + LHCb::Math::Similarity( Q, FCFinv );
          LHCb::Math::Similarity( A, sum, node.final_state_cov );
        } else {
          // if there is no noise, the gain matrix is just the inverse of
          // the transport matrix
          A = F;
          A.Invert();
          // the update of the covariance matrix becomes a lot simpler
          LHCb::Math::Similarity( A, nextnode.final_state_cov, node.final_state_cov );
        }
        node.final_state_vec += A * ( nextnode.final_state_vec - nextnode.predicted_state_vec[Node::forward] );
        // FIXME should we set ref vector and update ref-residual?
        // I think so right? best info for alignment would be the updated info?
        node.ref_vector() = node.final_state_vec;
        // updates H and the ref residual based on the new ref_vector
        node.project_reference();
      }
    }

    inline Gaudi::Matrix1x6 alignmentDerivatives( const Node& node, const Gaudi::XYZPoint& pivot ) {
      Gaudi::Vector3 const r_pos{node.lhcb_ref_vector.x(), node.lhcb_ref_vector.y(), node.z()};
      Gaudi::Vector3 const r_dir{node.lhcb_ref_vector.tx(), node.lhcb_ref_vector.ty(), 1.0};

      Gaudi::Vector3 dist = r_pos - node.measurement_pos;
      Gaudi::Vector3 d0   = node.measurement_dir;
      Gaudi::Vector3 d1   = r_dir;
      Gaudi::Vector3 c1{0, 0, 0};
      const auto     mat    = std::array{Dot( d0, d0 ), -Dot( d1, d0 ), Dot( d1, d1 )};
      auto           mu     = std::array{-Dot( dist, d0 ), Dot( d1, dist )};
      const auto     decomp = ROOT::Math::CholeskyDecomp<double, 2>{mat.data()};

      if ( !decomp.Solve( mu ) ) throw std::runtime_error( "singular matrix" );

      dist += mu[0] * d0 - mu[1] * d1;
      // Set up the vector onto which we project everything. This should
      // actually be parallel to dist.
      auto             unit_poca = Cross( d0, d1 ).Unit();
      Gaudi::Matrix1x3 dual( unit_poca.Array(), 3 );

      ROOT::Math::SMatrix<double, 3, 5> deriv{};
      deriv( 0, 0 ) = deriv( 1, 1 ) = 1.;
      deriv( 0, 2 ) = deriv( 1, 3 ) = mu[1];

      // Calculate the derivative of the poca on measTraj to alignment parameters.
      // Only non-zero elements:
      Gaudi::Vector3                    tmp = node.measurement_pos - node.measurement_dir * mu[0];
      Gaudi::XYZPoint                   meas_poca{tmp[0], tmp[1], tmp[2]};
      Gaudi::XYZVector                  arm = meas_poca - pivot;
      ROOT::Math::SMatrix<double, 3, 6> dPosdAlpha;
      // Derivative to translation
      dPosdAlpha( 0, 0 ) = dPosdAlpha( 1, 1 ) = dPosdAlpha( 2, 2 ) = 1;
      // Derivative to rotation around x-axis
      dPosdAlpha( 1, 3 ) = -arm.z();
      dPosdAlpha( 2, 3 ) = arm.y();
      // Derivative to rotation around y-axis
      dPosdAlpha( 0, 4 ) = arm.z();
      dPosdAlpha( 2, 4 ) = -arm.x();
      // Derivative to rotation around z-axis
      dPosdAlpha( 0, 5 ) = -arm.y();
      dPosdAlpha( 1, 5 ) = arm.x();

      return dual * dPosdAlpha;
    }

    inline StatusCode update_transport_matrix( LHCb::span<Node> fitnodes, IGeometryInfo const& geometry,
                                               const ITrackExtrapolator& extrap ) {
      mydebug( "Start of update_transport" );
      Gaudi::TrackMatrix F{ROOT::Math::SMatrixIdentity()};
      auto               lhcb_ref_vector{std::cref( fitnodes.front().lhcb_ref_vector )};

      StatusCode ret;
      for ( auto node = std::next( fitnodes.begin() ); node != fitnodes.end(); ++node ) {
        // wanted explicit copy
        auto       state_vec = lhcb_ref_vector.get();
        auto const z         = node->z();
        ret                  = extrap.propagate( state_vec, z, geometry, &F );
        if ( ret.isFailure() ) { return ret; }

        apply_delta_E( state_vec, node->delta_energy );

        Gaudi::TrackVector transportvec = state_vec.parameters() - F * lhcb_ref_vector.get().parameters();
        node->set_transport( F, transportvec );
        // update the reference
        lhcb_ref_vector = std::cref( node->lhcb_ref_vector );
      }
      mydebug( "End of update_transport" );
      return ret;
    }
    template <typename Buffer, typename Buffer2>
    auto iterate_fit( LHCb::span<Node> fitnodes, const FitConfiguration& fit_config, const IGeometryInfo& geo,
                      const ITrackExtrapolator& extrap, Buffer& iter_buffer, Buffer2& pre_outlier_chi2_cut_buffer,
                      Buffer& transport_failed_buffer ) {

      pre_fit_init( fitnodes, fit_config );
      // first fit execution
      auto       prev_chi2 = fit( fitnodes );
      auto const tolerance = prev_chi2.nDoF() * 0.01;
      // iterate until max iterations or chi2 change below tolerance
      // iter starts at 2 to align max_fit_iter property with TrackMasterFitter
      int iter{2};
      for ( bool has_converged{false}; iter <= fit_config.max_fit_iter && !has_converged; ++iter ) {
        mydebug( "Iteration:", iter );
        smooth_and_update_ref_vector( fitnodes );
        for ( auto& node : fitnodes ) node.project_reference();
        if ( update_transport_matrix( fitnodes, geo, extrap ).isFailure() ) {
          ++transport_failed_buffer;
          iter_buffer += ( iter - 1 );
          return std::tuple{prev_chi2, false, iter - 1};
        }
        auto const chi2  = fit( fitnodes );
        auto const dchi2 = prev_chi2.chi2() - chi2.chi2();
        prev_chi2        = chi2;
        has_converged    = std::abs( dchi2 ) < tolerance;
      }
      iter_buffer += ( iter - 1 );
      if ( prev_chi2.chi2PerDoF() >= fit_config.maxchi2perdof_pre_outlier ) {
        ++pre_outlier_chi2_cut_buffer;
        return std::tuple{prev_chi2, false, iter - 1};
      }
      return std::tuple{prev_chi2, true, iter - 1};
    }

    /**
     * @brief OUTLIER REMOVAL
     * @noteFIXME this is currently more of an ad-hoc implementation
     *Check if it's faster to first check based on each nodes delta_chi2
     *similar to what is done in TMF.
     *Also given that I now need to smooth preety much every time I also fit
     *it should be checked if It's faster to fuse the backwards and smoothing
     *operations
     */
    template <typename Buffer, typename Buffer2>
    auto remove_outliers( LHCb::span<Node> fitnodes, const FitConfiguration& fit_config, LHCb::ChiSquare& prev_chi2,
                          Buffer& outlier_iter_buffer, Buffer2& chi2_cut_buffer ) {

      const auto min_hits = std::array<size_t, 4>{fit_config.minNumVPLayers, fit_config.minNumUTLayers,
                                                  fit_config.minNumFTLayers, fit_config.minNumMuonLayers};
      // VP=0, UT=1, FT=2, Muon=3
      LHCb::HitPattern pattern{LHCb::TransformedRange{fitnodes, []( const auto& node ) { return node.lhcbID; },
                                                      []( const auto& node ) { return node.isHitOnTrack(); }}};

      int outlier_iter{0};
      for ( ; outlier_iter < fit_config.max_outlier_iter; ++outlier_iter ) {
        smooth_and_update_ref_vector_outlier( fitnodes );

        const auto num_hits = std::array{pattern.numVelo(), pattern.numUT(), pattern.numFT(), pattern.numMuon()};

        Node* outlier     = nullptr;
        auto  worstChi2   = fit_config.min_chi2_for_outlier;
        auto  outlier_pos = fitnodes.begin();
        auto  it          = fitnodes.begin();
        for ( auto& node : fitnodes ) {
          if ( node.m_is_outlier ) {
            ++it;
            continue;
          }
          node.project_reference();
          auto type_idx = static_cast<int>( node.type() ) - 1;
          // not allowed to remove this one so skip the chi2 calc
          if ( num_hits[type_idx] <= min_hits[type_idx] ) {
            ++it;
            continue;
          }

          double V    = node.measurement_error * node.measurement_error;
          double HCH  = LHCb::Math::Similarity( node.projection, node.final_state_cov )[0][0];
          double R    = V - HCH;
          double chi2 = ( node.ref_residual * node.ref_residual ) / R;

          if ( chi2 > worstChi2 ) {
            worstChi2   = chi2;
            outlier     = &node;
            outlier_pos = it;
          }
          ++it;
        }

        // if I don't find an outlier, we are done.
        if ( !outlier ) break;
        if ( outlier->m_type == Node::Type::VPHit ) {

          auto index = std::distance( outlier_pos, fitnodes.end() );
          // our partner is +1 or -1 depending on if we are 1st or 2nd part of a vp-hit
          auto& partner        = *( outlier_pos + ( index % 2 == 0 ? 1 : -1 ) );
          partner.m_is_outlier = true;
          pattern.remove( partner.lhcbID );
        }
        outlier->m_is_outlier = true;
        pattern.remove( outlier->lhcbID );
        prev_chi2 = fit( fitnodes );
      }
      outlier_iter_buffer += outlier_iter;

      if ( prev_chi2.chi2PerDoF() >= fit_config.maxchi2perdof ) {
        ++chi2_cut_buffer;
        return false;
      }
      return true;
    }

  } // namespace KF

  inline double calc_ctb_z( LHCb::State const& state ) {
    auto const& vec = state.stateVector();
    auto        z   = state.z();
    // check on division by zero (track parallel to beam line!)
    if ( vec[2] != 0 || vec[3] != 0 ) {
      z -= ( vec[0] * vec[2] + vec[1] * vec[3] ) / ( vec[2] * vec[2] + vec[3] * vec[3] );
    }
    return z;
  }

  inline auto calc_extra_info( const std::vector<Node>& fitnodes ) {

    auto velo_chi2     = LHCb::ChiSquare{0, -5};
    auto down_chi2     = LHCb::ChiSquare{0, -5};
    auto upstream_chi2 = LHCb::ChiSquare{};
    auto NUTOutliers   = 0;

    for ( auto const& node : fitnodes ) {
      switch ( node.type() ) {
      case Node::Type::VPHit:
        velo_chi2 += node.delta_chi2[Node::backward];
        break;
      case Node::Type::UTHit:
        upstream_chi2 += node.delta_chi2[Node::backward];
        if ( node.m_is_outlier ) ++NUTOutliers;
        break;
      case Node::Type::FTHit:
        down_chi2 += node.delta_chi2[Node::forward];
        break;
      default:
        break;
      }
    }
    upstream_chi2 += velo_chi2;

    return std::tuple{velo_chi2, upstream_chi2, down_chi2, NUTOutliers};
  }

  inline auto make_fit_result( std::vector<Node>& fitnodes, int nIter, double scatteringMomentum,
                               bool classic_smoothing_post ) {

    auto kfr = LHCb::PrKalmanFitResult{};

    if ( classic_smoothing_post ) {
      // run a classical smoothing iteration for the alignment
      // and store the gain matrices from that step
      kfr.gain_matrices.resize( fitnodes.size() );
      KF::classical_smoother_iteration_for_alignment( fitnodes, kfr.gain_matrices );
    }

    // transfer our nodes into the FitResult
    kfr.fitnodes            = std::move( fitnodes );
    kfr.scattering_momentum = scatteringMomentum;
    kfr.number_of_iter      = nIter;
    // need to re-reserve since I just moved out of a container that I use every loop
    // that should make my "moved from" vector usable again.
    fitnodes.reserve( 50 );
    return kfr;
  }

  template <typename TrackProxy>
  void update_state( TrackProxy& newtrack, const LHCb::State& state ) {

    if constexpr ( std::is_same_v<TrackProxy, LHCb::Event::v1::Track> ) {
      auto& statevec = newtrack.states();
      statevec.emplace_back( new LHCb::State{std::move( state )} );
    } else {
      using SL       = LHCb::Event::v3::Tracks::StateLocation;
      const auto loc = [&] {
        switch ( state.location() ) {
        case State::ClosestToBeam:
          return SL::ClosestToBeam;
        case State::FirstMeasurement:
          return SL::FirstMeasurement;
        case State::LastMeasurement:
          return SL::LastMeasurement;
        case State::BegRich1:
          return SL::BegRich1;
        case State::EndRich1:
          return SL::EndRich1;
        case State::BegRich2:
          return SL::BegRich2;
        case State::EndRich2:
          return SL::EndRich2;
        default:
          throw GaudiException( "Invalid State Location", "PrKalmanFilter", StatusCode::FAILURE );
        }
      }();

      // positions, slopes, q/p
      namespace tag = LHCb::Event::v3::Tag;
      newtrack.template field<tag::States>( loc ).setPosition( state.x(), state.y(), state.z() );
      newtrack.template field<tag::States>( loc ).setDirection( state.tx(), state.ty() );
      newtrack.template field<tag::States>( loc ).setQOverP( state.qOverP() );

      // covariance
      auto const& cov = state.covariance();
      newtrack.template field<tag::StateCovs>( loc ).set(
          cov( 0, 0 ), cov( 0, 1 ), cov( 0, 2 ), cov( 0, 3 ), cov( 0, 4 ), cov( 1, 1 ), cov( 1, 2 ), cov( 1, 3 ),
          cov( 1, 4 ), cov( 2, 2 ), cov( 2, 3 ), cov( 2, 4 ), cov( 3, 3 ), cov( 3, 4 ), cov( 4, 4 ) );
    }
  }

  template <typename InputTrackType, typename TrackProxy>
  StatusCode add_fitted_states( TrackProxy& newtrack, LHCb::span<Node> fitnodes, double const scatteringMomentum,
                                IGeometryInfo const& geometry, const ITrackExtrapolator& extrap ) {
    if constexpr ( std::is_same_v<TrackProxy, LHCb::Event::v1::Track> )
      assert( newtrack.states().size() == 0 && "this should be a new track without states on it" );

    Gaudi::TrackMatrix transMat              = ROOT::Math::SMatrixIdentity();
    auto const         extrapolate_and_noise = [&extrap, scatteringMomentum, &transMat, &geometry](
                                           LHCb::State& state, double z_target, auto node_type, auto prev_node_type ) {
      auto&      vec   = state.stateVector();
      auto const z_old = state.z();
      return extrap.propagate( vec, state.z(), z_target, &transMat, geometry ).andThen( [&] {
        state.setZ( z_target );
        // transport
        state.covariance()   = LHCb::Math::Similarity( transMat, state.covariance() );
        auto const [cov, dE] = TrackFit::param_scatter_impl::computeNoiseAndDeltaE(
            node_type, vec[0], vec[1], state.z(), vec[2], vec[3], prev_node_type, z_old, scatteringMomentum );

        state.covariance() += cov;
        apply_delta_E( state, dE );
      } );
    };

    auto const& first_hit = fitnodes.back();
    auto const& last_hit  = fitnodes.front();
    auto const  fwd       = Node::forward;
    auto const  bkwd      = Node::backward;

    StatusCode ret = StatusCode::SUCCESS;

    // for long, velo and upstream we add CTB first measurement here
    if constexpr ( isVelo<InputTrackType> || isLong<InputTrackType> || isUpstream<InputTrackType> ) {
      auto state = LHCb::State{first_hit.filtered_state_vec[fwd], first_hit.filtered_state_cov[fwd], first_hit.z(),
                               LHCb::State::ClosestToBeam};

      // don't go outside the Velo volume
      // if CTB state is outside Velo: define CTB state to be the same as FirstMeasurement
      auto is_ctb_first = true;
      if ( const auto ctb_z = calc_ctb_z( state ); ctb_z < StateParameters::ZEndVelo && ctb_z > ZBegVelo ) {
        ret = extrapolate_and_noise( state, ctb_z, TrackFit::param_scatter_impl::NodeType::ClosestToBeam,
                                     TrackFit::param_scatter_impl::NodeType::VPHit );
        if ( ret.isFailure() ) { return ret; }

        is_ctb_first = ctb_z < first_hit.z();
        if constexpr ( std::is_same_v<TrackProxy, LHCb::Event::v1::Track> ) {
          // no effect for v3::Tracks
          if ( newtrack.isVeloBackward() ) is_ctb_first = !is_ctb_first;
        }
      }

      auto first_meas = LHCb::State{first_hit.filtered_state_vec[fwd], first_hit.filtered_state_cov[fwd], first_hit.z(),
                                    LHCb::State::Location::FirstMeasurement};

      if ( is_ctb_first ) {
        update_state( newtrack, state );
        update_state( newtrack, first_meas );
      } else {
        update_state( newtrack, first_meas );
        update_state( newtrack, state );
      }
    }

    if constexpr ( isLong<InputTrackType> || isUpstream<InputTrackType> ) {
      auto const [upstream_hit, downstream_hit] = [&]() {
        auto const up = std::find_if( fitnodes.begin(), fitnodes.end(),
                                      []( auto const& node ) { return node.z() < StateParameters::ZBegRich1; } );
        assert( up != fitnodes.begin() );
        return std::pair{up, std::prev( up )};
      }();

      // TODO speedup by not using the big FitNode classes here
      std::vector<Node> rich1_nodes;
      rich1_nodes.reserve( 4 );
      rich1_nodes.push_back( *downstream_hit );
      rich1_nodes.emplace_back( 0, 0, StateParameters::ZEndRich1, 0, 0, 0, 0, Node::Type::EndRich1, LHCbID{} );
      rich1_nodes.emplace_back( 0, 0, StateParameters::ZBegRich1, 0, 0, 0, 0, Node::Type::BegRich1, LHCbID{} );
      rich1_nodes.push_back( *upstream_hit );

      auto               lhcb_ref_vector{std::cref( rich1_nodes.begin()->lhcb_ref_vector )};
      Gaudi::TrackMatrix TransMat{ROOT::Math::SMatrixIdentity()};

      for ( auto node = std::next( rich1_nodes.begin() ); node != rich1_nodes.end(); ++node ) {
        // wanted explicit copy
        auto       state_vec = lhcb_ref_vector.get();
        auto const z         = node->z();
        ret                  = extrap.propagate( state_vec, z, geometry, &TransMat );
        if ( ret.isFailure() ) { return ret; }

        node->lhcb_ref_vector.parameters() = state_vec.parameters();

        std::tie( node->noise_matrix, node->delta_energy ) =
            TrackFit::param_scatter_impl::computeNoiseAndDeltaE( *( node - 1 ), *node, scatteringMomentum );

        apply_delta_E( state_vec, node->delta_energy );

        Gaudi::TrackVector transportvec = state_vec.parameters() - TransMat * lhcb_ref_vector.get().parameters();
        node->set_transport( TransMat, transportvec );
        // update the reference
        lhcb_ref_vector = std::cref( node->lhcb_ref_vector );
      }

      KF::predict<Node::forward>( rich1_nodes[0], rich1_nodes[1] );

      rich1_nodes[1].filtered_state_vec[Node::forward] = rich1_nodes[1].predicted_state_vec[Node::forward];
      rich1_nodes[1].filtered_state_cov[Node::forward] = rich1_nodes[1].predicted_state_cov[Node::forward];

      KF::predict<Node::forward>( rich1_nodes[1], rich1_nodes[2] );
      KF::predict<Node::backward>( rich1_nodes[3], rich1_nodes[2] );

      rich1_nodes[2].filtered_state_vec = rich1_nodes[2].predicted_state_vec;
      rich1_nodes[2].filtered_state_cov = rich1_nodes[2].predicted_state_cov;

      KF::predict<Node::backward>( rich1_nodes[2], rich1_nodes[1] );

      rich1_nodes[1].filtered_state_vec[Node::backward] = rich1_nodes[1].predicted_state_vec[Node::backward];
      rich1_nodes[1].filtered_state_cov[Node::backward] = rich1_nodes[1].predicted_state_cov[Node::backward];
      KF::average_node( rich1_nodes[1] );
      KF::average_node( rich1_nodes[2] );

      update_state( newtrack, LHCb::State{rich1_nodes[2].final_state_vec, rich1_nodes[2].final_state_cov,
                                          rich1_nodes[2].z(), LHCb::State::Location::BegRich1} );
      update_state( newtrack, LHCb::State{rich1_nodes[1].final_state_vec, rich1_nodes[1].final_state_cov,
                                          rich1_nodes[1].z(), LHCb::State::Location::EndRich1} );

    } else if constexpr ( isSeed<InputTrackType> ) {

      update_state( newtrack, LHCb::State{first_hit.filtered_state_vec[fwd], first_hit.filtered_state_cov[fwd],
                                          first_hit.z(), LHCb::State::Location::FirstMeasurement} );

    } else if constexpr ( isDownstream<InputTrackType> ) {
      // for downstream tracks, first measurement is after RICH1
      // thus create state here but first add the RICH states
      auto first_meas = LHCb::State{first_hit.filtered_state_vec[fwd], first_hit.filtered_state_cov[fwd], first_hit.z(),
                                    LHCb::State::Location::FirstMeasurement};

      // RICH1 states are much simpler because they are only an extrapolation
      auto end_rich1 = LHCb::State{first_meas};
      ret            = extrapolate_and_noise( end_rich1, StateParameters::ZEndRich1,
                                   TrackFit::param_scatter_impl::NodeType::EndRich1,
                                   TrackFit::param_scatter_impl::nodetype( first_hit ) );
      if ( ret.isFailure() ) { return ret; }
      end_rich1.setLocation( LHCb::State::Location::EndRich1 );

      // to keep the right order we first insert the BegRICH1 state
      auto beg_rich1 = LHCb::State{end_rich1};
      ret            = extrapolate_and_noise( beg_rich1, StateParameters::ZBegRich1,
                                   TrackFit::param_scatter_impl::NodeType::BegRich1,
                                   TrackFit::param_scatter_impl::NodeType::EndRich1 );
      if ( ret.isFailure() ) { return ret; }
      beg_rich1.setLocation( LHCb::State::Location::BegRich1 );
      update_state( newtrack, beg_rich1 );
      // then the EndRICH1 state
      update_state( newtrack, end_rich1 );
      // now add first measurement
      update_state( newtrack, first_meas );
    }

    auto last_meas = LHCb::State{last_hit.filtered_state_vec[bkwd], last_hit.filtered_state_cov[bkwd], last_hit.z(),
                                 LHCb::State::Location::LastMeasurement};

    auto last_hit_type = TrackFit::param_scatter_impl::nodetype( last_hit );

    if ( last_hit.isMuon() ) { // for LongMuon tracks change last_meas to be last SciFi hit
      const auto it = std::find_if( fitnodes.begin(), fitnodes.end(), []( auto const& node ) {
        return node.type() == LHCb::Pr::Tracks::Fit::Node::Type::FTHit;
      } );

      last_meas = LHCb::State{( *it ).filtered_state_vec[bkwd], ( *it ).filtered_state_cov[bkwd], ( *it ).z(),
                              LHCb::State::Location::AtT};

      last_hit_type = TrackFit::param_scatter_impl::NodeType::FTHit;
    }

    update_state( newtrack, last_meas );

    // Rich2 only for downstream and long and seed
    if constexpr ( isLong<InputTrackType> || isDownstream<InputTrackType> || isSeed<InputTrackType> ) {
      // begRich2 is copy of last measurement
      // now extrapolate
      auto beg_rich2 = LHCb::State{last_meas};
      ret            = extrapolate_and_noise( beg_rich2, StateParameters::ZBegRich2,
                                   TrackFit::param_scatter_impl::NodeType::BegRich2, last_hit_type );
      if ( ret.isFailure() ) { return ret; }
      beg_rich2.setLocation( LHCb::State::Location::BegRich2 );

      auto end_rich2 = LHCb::State{beg_rich2};
      ret            = extrapolate_and_noise( end_rich2, StateParameters::ZEndRich2,
                                   TrackFit::param_scatter_impl::NodeType::EndRich2,
                                   TrackFit::param_scatter_impl::NodeType::BegRich2 );
      if ( ret.isFailure() ) { return ret; }
      end_rich2.setLocation( LHCb::State::Location::EndRich2 );

      update_state( newtrack, beg_rich2 );
      update_state( newtrack, end_rich2 );

      if ( last_hit.isMuon() ) {
        update_state( newtrack, LHCb::State{last_hit.filtered_state_vec[bkwd], last_hit.filtered_state_cov[bkwd],
                                            last_hit.z(), LHCb::State::Location::LastMeasurement} );
      }
    }
    return ret;
  }

  template <typename InputTrackType>
  auto make_output_track( const InputTrackType& tracks, const proxy_type<InputTrackType>& track,
                          std::vector<Node>& fitnodes, LHCb::ChiSquare prev_chi2, int nIter, double scatteringMomentum,
                          bool fill_fitresult, bool classic_smoothing_post, const IGeometryInfo& geo,
                          const ITrackExtrapolator& extrap ) {
    //
    // start creating our output track
    //
    using out_track      = LHCb::Event::v1::Track;
    using Type           = LHCb::Event::Enum::Track::Type;
    using PatRecStatus   = LHCb::Event::Enum::Track::PatRecStatus;
    using History        = LHCb::Event::Enum::Track::History;
    using AdditionalInfo = LHCb::Event::Enum::Track::AdditionalInfo;
    using FitStatus      = LHCb::Event::Enum::Track::FitStatus;

    auto new_track = [&] {
      if constexpr ( isLong<InputTrackType> ) {

        return std::make_unique<out_track>( tracks.history(), Type::Long, PatRecStatus::PatRecIDs );

      } else if constexpr ( isVelo<InputTrackType> ) {

        if ( track.backward() )
          return std::make_unique<out_track>( History::PrPixel, Type::VeloBackward, PatRecStatus::PatRecIDs );
        else
          return std::make_unique<out_track>( History::PrPixel, Type::Velo, PatRecStatus::PatRecIDs );

      } else if constexpr ( isSeed<InputTrackType> ) {

        return std::make_unique<out_track>( History::PrSeeding, Type::Ttrack, PatRecStatus::PatRecIDs );

      } else if constexpr ( isDownstream<InputTrackType> ) {

        return std::make_unique<out_track>( History::PrDownstream, Type::Downstream, PatRecStatus::PatRecIDs );

      } else if constexpr ( isUpstream<InputTrackType> ) {

        return std::make_unique<out_track>( History::PrVeloUT, Type::Upstream, PatRecStatus::PatRecIDs );
      }
    }();

    // set the LHCbIds of our new track
    new_track->setLhcbIDs( track.lhcbIDs() );

    new_track->setFitStatus( FitStatus::Fitted );
    new_track->setNDoF( prev_chi2.nDoF() );
    new_track->setChi2PerDoF( prev_chi2.chi2() / prev_chi2.nDoF() );

    const auto [velo_chi2, upstream_chi2, down_chi2, n_ut_outliers] = calc_extra_info( fitnodes );

    if constexpr ( isLong<InputTrackType> ) {

      new_track->addInfo( AdditionalInfo::FitTChi2, down_chi2.chi2() );
      new_track->addInfo( AdditionalInfo::FitTNDoF, down_chi2.nDoF() );
      new_track->addInfo( AdditionalInfo::FitVeloChi2, velo_chi2.chi2() );
      new_track->addInfo( AdditionalInfo::FitVeloNDoF, velo_chi2.nDoF() );
      new_track->addInfo( AdditionalInfo::FitMatchChi2, prev_chi2.chi2() - upstream_chi2.chi2() - down_chi2.chi2() );
      new_track->addInfo( AdditionalInfo::NUTOutliers, n_ut_outliers );
    } else if constexpr ( isVelo<InputTrackType> ) {

      new_track->addInfo( AdditionalInfo::FitVeloChi2, velo_chi2.chi2() );
      new_track->addInfo( AdditionalInfo::FitVeloNDoF, velo_chi2.nDoF() );

    } else if constexpr ( isDownstream<InputTrackType> || isSeed<InputTrackType> ) {

      new_track->addInfo( AdditionalInfo::FitTChi2, down_chi2.chi2() );
      new_track->addInfo( AdditionalInfo::FitTNDoF, down_chi2.nDoF() );
      new_track->addInfo( AdditionalInfo::NUTOutliers, n_ut_outliers );

    } else if constexpr ( isUpstream<InputTrackType> ) {

      new_track->addInfo( AdditionalInfo::FitVeloChi2, velo_chi2.chi2() );
      new_track->addInfo( AdditionalInfo::FitVeloNDoF, velo_chi2.nDoF() );
      new_track->addInfo( AdditionalInfo::NUTOutliers, n_ut_outliers );
    }

    if ( add_fitted_states<InputTrackType>( *new_track, fitnodes, scatteringMomentum, geo, extrap ).isFailure() ) {
      return std::unique_ptr<out_track>{nullptr};
    }

    if ( fill_fitresult ) {
      auto kfr = make_fit_result( fitnodes, nIter, scatteringMomentum, classic_smoothing_post );
      new_track->setFitResult( new LHCb::PrKalmanFitResult{std::move( kfr )} );
    }

    return new_track;
  }

  template <typename InputTrackType>
  StatusCode add_output_v3_track( LHCb::Event::v3::Tracks& new_tracks, const InputTrackType& tracks,
                                  const proxy_type<InputTrackType>& track, std::vector<Node>& fitnodes,
                                  LHCb::ChiSquare prev_chi2, double scatteringMomentum, const IGeometryInfo& geo,
                                  const ITrackExtrapolator& extrap, LHCb::UniqueIDGenerator const& unique_id_gen ) {
    //
    // start creating our output track
    //

    static_assert( isLong<InputTrackType> || isVelo<InputTrackType> || isDownstream<InputTrackType> ||
                       isSeed<InputTrackType> || isUpstream<InputTrackType>,
                   "Only Velo, Downstream, Long, Seed or Upstream tracks are supported currently." );

    auto new_track = new_tracks.emplace_back<SIMDWrapper::InstructionSet::Scalar>();

    namespace tag = LHCb::Event::v3::Tag;
    using int_v   = decltype( new_track.template field<tag::UniqueID>().get() );
    new_track.field<tag::UniqueID>().set( unique_id_gen.generate<int_v>().value() );

    const auto history = [&] {
      if constexpr ( isLong<InputTrackType> )
        return tracks.history();
      else if constexpr ( isDownstream<InputTrackType> )
        return LHCb::Event::Enum::Track::History::PrDownstream;
      else if constexpr ( isVelo<InputTrackType> )
        return LHCb::Event::Enum::Track::History::PrPixel;
      else if constexpr ( isSeed<InputTrackType> )
        return LHCb::Event::Enum::Track::History::PrSeeding;
      else if constexpr ( isUpstream<InputTrackType> )
        return LHCb::Event::Enum::Track::History::PrVeloUT;
    }();
    new_track.field<tag::history>().set( history );
    new_track.field<tag::nDoF>().set( prev_chi2.nDoF() );
    new_track.field<tag::Chi2>().set( prev_chi2.chi2() );

    if constexpr ( hasVP<InputTrackType> ) {
      if constexpr ( isVelo<InputTrackType> ) {
        new_track.field<tag::trackVP>().set( track.indices() );
      } else {
        new_track.field<tag::trackVP>().set( track.trackVP() );
      }
      auto n_vp_hits = track.nVPHits();
      new_track.template field<tag::VPHits>().resize( n_vp_hits );
      for ( int i{0}; i < n_vp_hits; ++i ) {
        auto id = track.vp_lhcbID( i );
        new_track.template field<tag::VPHits>()[i].template field<tag::LHCbID>().set( id );
      }
    } else {
      new_track.field<tag::trackVP>().set( I{-1} );
    }

    if constexpr ( hasUT<InputTrackType> ) {
      auto n_ut_hits = track.nUTHits();
      new_track.template field<tag::UTHits>().resize( n_ut_hits );
      for ( int i{0}; i < n_ut_hits; ++i ) {
        auto id = track.ut_lhcbID( i );
        new_track.template field<tag::UTHits>()[i].template field<tag::LHCbID>().set( id );
      }
      if constexpr ( isLong<InputTrackType> )
        new_track.field<tag::trackUT>().set( track.trackUT() );
      else
        new_track.field<tag::trackUT>().set( I{-1} );
    } else {
      new_track.field<tag::trackUT>().set( I{-1} );
    }

    if constexpr ( hasFT<InputTrackType> ) {
      auto n_ft_hits = track.nFTHits();
      new_track.template field<tag::FTHits>().resize( n_ft_hits );
      for ( int i{0}; i < n_ft_hits; ++i ) {
        auto id = track.ft_lhcbID( i );
        new_track.template field<tag::FTHits>()[i].template field<tag::LHCbID>().set( id );
      }
      if constexpr ( isDownstream<InputTrackType> )
        new_track.field<tag::trackSeed>().set( track.seed_track_index() );
      else if constexpr ( isLong<InputTrackType> )
        new_track.field<tag::trackSeed>().set( track.trackSeed() );
      else if constexpr ( isSeed<InputTrackType> )
        new_track.field<tag::trackSeed>().set( track.indices() );
    } else {
      new_track.field<tag::trackSeed>().set( I{-1} );
    }

    StatusCode ret = add_fitted_states<InputTrackType>( new_track, fitnodes, scatteringMomentum, geo, extrap );
    return ret;
  }

  inline void add_output_v3_partial_chi2( LHCb::Event::v3::Track::PartialChiSquareds& out_partial_chi2s,
                                          std::vector<Node>& fitnodes, LHCb::ChiSquare prev_chi2 ) {

    auto new_partial_chi2 = out_partial_chi2s.emplace_back<SIMDWrapper::InstructionSet::Scalar>();

    const auto [velo_chi2, upstream_chi2, down_chi2, NUTOutliers] = calc_extra_info( fitnodes );

    namespace tag = LHCb::Event::v3::Track::PartialChiSquaredsTag;
    new_partial_chi2.field<tag::FitVeloChi2>().set( velo_chi2.chi2() );
    new_partial_chi2.field<tag::FitVeloNDoF>().set( velo_chi2.nDoF() );
    new_partial_chi2.field<tag::FitTChi2>().set( down_chi2.chi2() );
    new_partial_chi2.field<tag::FitTNDoF>().set( down_chi2.nDoF() );
    new_partial_chi2.field<tag::FitMatchChi2>().set( prev_chi2.chi2() - upstream_chi2.chi2() - down_chi2.chi2() );
    new_partial_chi2.field<tag::NUTOutliers>().set( NUTOutliers );
  }

  template <typename InputTrackType, typename OutputTrackType>
  inline auto initialize_output( int size, LHCb::UniqueIDGenerator const& unique_id_gen, bool backward = false ) {
    if constexpr ( isV1Tracks<OutputTrackType> ) {
      auto out_tracks = LHCb::Event::v1::Tracks{};
      out_tracks.reserve( size );
      return std::make_tuple( std::move( out_tracks ) );
    } else if constexpr ( isV3Tracks<OutputTrackType> ) {
      auto out_tracks = LHCb::Event::v3::Tracks{v3_track_type<InputTrackType>(), backward, unique_id_gen};
      out_tracks.reserve( size );
      return std::make_tuple( std::move( out_tracks ) );
    } else if constexpr ( isV3TracksExtra<OutputTrackType> || isV3TracksFull<OutputTrackType> ) {
      auto out_tracks = LHCb::Event::v3::Tracks{v3_track_type<InputTrackType>(), backward, unique_id_gen};
      out_tracks.reserve( size );
      auto out_partial_chi2s = LHCb::Event::v3::Track::PartialChiSquareds{out_tracks.zipIdentifier()};
      out_partial_chi2s.reserve( size );
      if constexpr ( isV3TracksExtra<OutputTrackType> ) {
        return std::make_tuple( std::move( out_tracks ), std::move( out_partial_chi2s ) );
      } else if constexpr ( isV3TracksFull<OutputTrackType> ) {
        std::vector<LHCb::PrKalmanFitResult> fit_results;
        fit_results.reserve( size );
        return std::make_tuple( std::move( out_tracks ), std::move( out_partial_chi2s ), std::move( fit_results ) );
      }
    } else
      throw GaudiException( "Output type is not supported", "PrKalmanFilter", StatusCode::FAILURE );
  }
} // namespace LHCb::Pr::Tracks::Fit
