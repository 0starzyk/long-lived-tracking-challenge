/*****************************************************************************\
* (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/PrFitNode.h"
#include "Event/UniqueIDGenerator.h"
#include "GaudiAlg/FunctionalTool.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IBinder.h"
#include "Kernel/LHCbID.h"
#include "Kernel/STLExtensions.h"
#include "LHCbMath/bit_cast.h"
#include "PrKalmanFilter/KF.h"
#include "TrackInterfaces/IPrFitterTool.h"

#include <boost/range/adaptor/reversed.hpp>
#include <vector>

namespace LHCb::Pr {
  namespace {
    using TracksV1           = LHCb::Event::v1::Tracks;
    using TracksV3           = LHCb::Event::v3::Tracks;
    using PartialChiSquareds = LHCb::Event::v3::Track::PartialChiSquareds;
    using TrackFitResults    = std::vector<LHCb::PrKalmanFitResult>;
    using TrackV1            = LHCb::Event::v1::Track;

    using namespace LHCb::Pr::Tracks::Fit;

    template <typename... T>
    struct always_false : std::false_type {};

    TrackV1 invalid_track( TrackV1&& track ) {
      track.setFlag( TrackV1::Flags::Invalid, true );
      return std::move( track );
    }

    /**
     * @brief Create fitnodes from lhbids on a v1 track
     *
     * @tparam Types hit types for containers used in the pattern recognition
     * @param fitnodes vector that will hold the fitnode objects
     * @param track v1 track
     * @param containers tuple of the hit containers needed for the track type
     * @note The current implementation only works for Velo, Downstream, Upstream Long tracks but can be extended
     * by making this function undestand how to create nodes for other hit types.
     */
    template <HitType... Types>
    void get_fitnodes_from_lhcbids( std::vector<Node>& fitnodes, const TrackV1& track,
                                    std::tuple<const Hits<Types>&...> containers ) {

      if constexpr ( has_hit_type<HitType::VP, Types...> && !has_hit_type<HitType::FT, Types...> &&
                     !has_hit_type<HitType::UT, Types...> && !has_hit_type<HitType::Muon, Types...> ) {
        assert( track.checkType( TrackV1::Types::Velo ) || track.isVeloBackward() );
      } else if constexpr ( !has_hit_type<HitType::VP, Types...> && has_hit_type<HitType::FT, Types...> &&
                            has_hit_type<HitType::UT, Types...> && !has_hit_type<HitType::Muon, Types...> ) {
        assert( track.checkType( TrackV1::Types::Downstream ) );
      } else if constexpr ( has_hit_type<HitType::VP, Types...> && !has_hit_type<HitType::FT, Types...> &&
                            has_hit_type<HitType::UT, Types...> && !has_hit_type<HitType::Muon, Types...> ) {
        assert( track.checkType( TrackV1::Types::Upstream ) );
      } else if constexpr ( !has_hit_type<HitType::VP, Types...> && has_hit_type<HitType::FT, Types...> &&
                            !has_hit_type<HitType::UT, Types...> && !has_hit_type<HitType::Muon, Types...> ) {
        assert( track.checkType( TrackV1::Types::Ttrack ) );
      } else if constexpr ( has_hit_type<HitType::VP, Types...> && has_hit_type<HitType::FT, Types...> &&
                            has_hit_type<HitType::UT, Types...> && !has_hit_type<HitType::Muon, Types...> ) {
        assert( track.checkType( TrackV1::Types::Long ) || track.checkType( TrackV1::Types::Velo ) ||
                track.isVeloBackward() || track.checkType( TrackV1::Types::Downstream ) ||
                track.checkType( TrackV1::Types::Upstream ) );
      } else if constexpr ( has_hit_type<HitType::VP, Types...> && has_hit_type<HitType::FT, Types...> &&
                            has_hit_type<HitType::UT, Types...> && has_hit_type<HitType::Muon, Types...> ) {
        assert( track.checkType( TrackV1::Types::LongMuon ) );
      } else if constexpr ( has_hit_type<HitType::VP, Types...> && has_hit_type<HitType::FT, Types...> &&
                            !has_hit_type<HitType::UT, Types...> && !has_hit_type<HitType::Muon, Types...> ) {
        assert( track.checkType( TrackV1::Types::Long ) );
      } else if constexpr ( has_hit_type<HitType::VP, Types...> && has_hit_type<HitType::FT, Types...> &&
                            !has_hit_type<HitType::UT, Types...> && has_hit_type<HitType::Muon, Types...> ) {
        assert( track.checkType( TrackV1::Types::LongMuon ) );
      } else {
        static_assert( always_false<std::integral_constant<HitType, Types>...>::value,
                       "Hit Containers do not match any supported track type" );
      }

      fitnodes.clear();
      // FIXME: if UT -> lhcbids would be ordered correctly like they are in the detector this would give the nodes
      // in the right order
      // lhcbids should be VP -> UT -> FT and NOT VP -> FT -> UT ! will be fixed with LHCb!3226
      for ( auto lhcbid : boost::adaptors::reverse( track.lhcbIDs() ) ) {
        assert( lhcbid.isVP() || lhcbid.isUT() || lhcbid.isFT() || lhcbid.isMuon() );
        if constexpr ( has_hit_type<HitType::VP, Types...> ) {
          if ( lhcbid.isVP() ) {
            const auto& hits_vp = get<HitType::VP>( containers );
            const auto  vp_hits = hits_vp.scalar();
            const auto  vp_hit  = std::find_if( vp_hits.begin(), vp_hits.end(), [=]( const auto& proxy ) {
              const auto vp_id =
                  LHCb::Detector::VPChannelID{bit_cast<unsigned>( proxy.template get<VP::VPHitsTag::ChannelId>() )};
              return lhcbid == LHCb::LHCbID{vp_id};
            } );
            assert( vp_hit != vp_hits.end() );

            const auto pos = ( *vp_hit ).template get<VP::VPHitsTag::pos>().vec3();
            fitnodes.emplace_back( pos.x.cast(), pos.y.cast(), pos.z.cast(), 0, 1, 0, vp_err, Node::Type::VPHit,
                                   lhcbid );
            fitnodes.emplace_back( pos.x.cast(), pos.y.cast(), pos.z.cast(), 1, 0, 0, vp_err, Node::Type::VPHit,
                                   lhcbid );
            continue;
          }
        }
        if constexpr ( has_hit_type<HitType::FT, Types...> ) {
          if ( lhcbid.isFT() ) {
            const auto& hits_ft = get<HitType::FT>( containers );
            const auto  idx     = hits_ft.find_index_of( lhcbid );
            const auto  x       = hits_ft.x( idx );
            const auto  z       = hits_ft.z( idx );
            const auto  w       = hits_ft.w( idx );
            const auto  dxdy    = hits_ft.dxDy( idx );
            const auto  dzdy    = hits_ft.dzDy( idx );

            fitnodes.emplace_back( x, 0, z, dxdy, 1, dzdy, ft_err( w ), Node::Type::FTHit, lhcbid );
            continue;
          }
        }
        if constexpr ( has_hit_type<HitType::UT, Types...> ) {
          if ( lhcbid.isUT() ) {
            const auto& hits_ut = get<HitType::UT>( containers );
            const auto  ut_hits = hits_ut.scalar();
            const auto  ut_hit  = std::find_if( ut_hits.begin(), ut_hits.end(), [=]( const auto& proxy ) {
              const auto ut_id =
                  LHCb::Detector::UT::ChannelID{bit_cast<unsigned>( proxy.template get<UT::UTHitsTag::channelID>() )};
              return lhcbid == LHCb::LHCbID{ut_id};
            } );
            assert( ut_hit != ut_hits.end() );

            const auto x0   = ( *ut_hit ).template get<UT::UTHitsTag::xAtYEq0>();
            const auto z0   = ( *ut_hit ).template get<UT::UTHitsTag::zAtYEq0>();
            const auto w    = ( *ut_hit ).template get<UT::UTHitsTag::weight>();
            const auto dxdy = ( *ut_hit ).template get<UT::UTHitsTag::dxDy>();

            fitnodes.emplace_back( x0.cast(), 0, z0.cast(), dxdy.cast(), 1, 0, ut_err( w ), Node::Type::UTHit, lhcbid );
            continue;
          }
        }
        if constexpr ( has_hit_type<HitType::Muon, Types...> ) {
          if ( lhcbid.isMuon() ) {
            const auto& hits_muon = get<HitType::Muon>( containers );

            const auto tileid          = Detector::Muon::TileID{lhcbid.muonID()};
            const auto station         = tileid.station();
            const auto hits_in_station = hits_muon.hits( station );

            const auto it       = std::find_if( hits_in_station.begin(), hits_in_station.end(),
                                          [=]( auto hit ) { return tileid == hit.tile(); } );
            const auto muon_hit = *it;
            const auto x        = muon_hit.x();
            const auto z        = muon_hit.z();
            const auto y        = muon_hit.y();
            const auto dx       = muon_hit.dx();
            const auto dy       = muon_hit.dy();

            fitnodes.emplace_back( x, y, z, 0, 1, 0, 2 * dx * LHCb::Math::inv_sqrt_12, Node::Type::MuonHit,
                                   lhcbid ); // node in Y direction
            fitnodes.emplace_back( x, y, z, 1, 0, 0, 2 * dy * LHCb::Math::inv_sqrt_12, Node::Type::MuonHit,
                                   lhcbid ); // node in X direction
            continue;
          }
        }
      }
      assert( !fitnodes.empty() );
      // FIXME: with reasonably ordered lhcbid's this would be unnecessary! remove with LHCb!3226
      std::sort( fitnodes.begin(), fitnodes.end(), [&]( const auto& i, const auto& j ) { return i.z() > j.z(); } );
      assert( std::is_sorted( fitnodes.begin(), fitnodes.end(),
                              [&]( const auto& i, const auto& j ) { return i.z() > j.z(); } ) );
    }

    /**
     * @brief makes the v1 output track
     *
     * @param track v1 track
     * @param fitnodes vector of fitnodes
     * @param prev_chi2 fit chi2
     * @param scatteringMomentum
     * @param fill_fitresult
     * @param classic_smoothing_post
     * @param geo lhcb
     * @param extrap track extrapolator
     * @return TrackV1 if this fails the Invalid flag is set
     * @note The current implementation only knows Velo, Downstream and Long tracks but can easily be extended for other
     * track types.
     *
     */
    TrackV1 make_output_track( const TrackV1& track, std::vector<Node>& fitnodes, LHCb::ChiSquare prev_chi2, int nIter,
                               double scatteringMomentum, bool fill_fitresult, bool classic_smoothing_post,
                               const IGeometryInfo& geo, const ITrackExtrapolator& extrap ) {

      // start creating our output track
      auto new_track = TrackV1{track.history(), track.type(), track.patRecStatus()};

      // set the LHCbIds of our new track
      new_track.setSortedLhcbIDs( track.lhcbIDs() );

      new_track.setFitStatus( TrackV1::FitStatus::Fitted );
      new_track.setNDoF( prev_chi2.nDoF() );
      new_track.setChi2PerDoF( prev_chi2.chi2() / prev_chi2.nDoF() );

      auto velo_chi2 =
          LHCb::ChiSquare{0, -5}; // initialized with -5 because of the number of free param in the track fit
      auto down_chi2     = LHCb::ChiSquare{0, -5};
      auto upstream_chi2 = LHCb::ChiSquare{}; // not initialized because always used combined with one of the above
      auto muon_chi2     = LHCb::ChiSquare{};

      for ( auto const& node : fitnodes ) {
        switch ( node.type() ) {
        case Node::Type::VPHit:
          velo_chi2 += node.delta_chi2[Node::backward];
          break;
        case Node::Type::UTHit:
          upstream_chi2 += node.delta_chi2[Node::backward];
          break;
        case Node::Type::FTHit:
          down_chi2 += node.delta_chi2[Node::forward];
          break;
        case Node::Type::MuonHit:
          muon_chi2 += node.delta_chi2[Node::forward];
          break;
        default:
          break;
        }
      }

      switch ( new_track.type() ) {
      case TrackV1::Types::Long:
        upstream_chi2 += velo_chi2;
        new_track.addInfo( TrackV1::AdditionalInfo::FitTChi2, down_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitTNDoF, down_chi2.nDoF() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitVeloChi2, velo_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitVeloNDoF, velo_chi2.nDoF() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitMatchChi2,
                           prev_chi2.chi2() - upstream_chi2.chi2() - down_chi2.chi2() );
        if ( add_fitted_states<Long::Tracks>( new_track, fitnodes, scatteringMomentum, geo, extrap ).isFailure() ) {
          return invalid_track( std::move( new_track ) );
        }
        break;
      case TrackV1::Types::LongMuon:
        upstream_chi2 += velo_chi2;
        new_track.addInfo( TrackV1::AdditionalInfo::FitTChi2, down_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitTNDoF, down_chi2.nDoF() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitVeloChi2, velo_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitVeloNDoF, velo_chi2.nDoF() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitMuonChi2, muon_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitMatchChi2,
                           prev_chi2.chi2() - upstream_chi2.chi2() - down_chi2.chi2() - muon_chi2.chi2() );
        if ( add_fitted_states<Long::Tracks>( new_track, fitnodes, scatteringMomentum, geo, extrap ).isFailure() ) {
          return invalid_track( std::move( new_track ) );
        }
        break;
      case TrackV1::Types::Velo:
      case TrackV1::Types::VeloBackward:
        new_track.addInfo( TrackV1::AdditionalInfo::FitVeloChi2, velo_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitVeloNDoF, velo_chi2.nDoF() );
        if ( add_fitted_states<Velo::Tracks>( new_track, fitnodes, scatteringMomentum, geo, extrap ).isFailure() ) {
          return invalid_track( std::move( new_track ) );
        }
        break;
      case TrackV1::Types::Downstream:
        new_track.addInfo( TrackV1::AdditionalInfo::FitTChi2, down_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitTNDoF, down_chi2.nDoF() );
        if ( add_fitted_states<Downstream::Tracks>( new_track, fitnodes, scatteringMomentum, geo, extrap )
                 .isFailure() ) {
          return invalid_track( std::move( new_track ) );
        }
        break;
      case TrackV1::Types::Upstream:
        upstream_chi2 += velo_chi2;
        new_track.addInfo( TrackV1::AdditionalInfo::FitVeloChi2, velo_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitVeloNDoF, velo_chi2.nDoF() );
        if ( add_fitted_states<Upstream::Tracks>( new_track, fitnodes, scatteringMomentum, geo, extrap ).isFailure() ) {
          return invalid_track( std::move( new_track ) );
        }
        break;
      case TrackV1::Types::Ttrack:
        new_track.addInfo( TrackV1::AdditionalInfo::FitTChi2, down_chi2.chi2() );
        new_track.addInfo( TrackV1::AdditionalInfo::FitTNDoF, down_chi2.nDoF() );
        if ( add_fitted_states<Seeding::Tracks>( new_track, fitnodes, scatteringMomentum, geo, extrap ).isFailure() ) {
          return invalid_track( std::move( new_track ) );
        }
      default:
        break;
      }

      if ( fill_fitresult ) {
        auto* kfr = new LHCb::PrKalmanFitResult{};
        new_track.setFitResult( kfr );

        if ( classic_smoothing_post ) {
          // run a classical smoothing iteration for the alignment
          // and store the gain matrices from that step
          kfr->gain_matrices.resize( fitnodes.size() );
          KF::classical_smoother_iteration_for_alignment( fitnodes, kfr->gain_matrices );
        }

        // transfer our nodes into the FitResult
        kfr->fitnodes            = std::move( fitnodes );
        kfr->scattering_momentum = scatteringMomentum;
        kfr->number_of_iter      = nIter;
        // need to re-reserve since I just moved out of a container that I use every loop
        // that should make my "moved from" vector usable again.
        fitnodes.reserve( 50 );
      }
      return new_track;
    }

  } // namespace

  class KalmanFilterTool
      : public Gaudi::Functional::ToolBinder<
            Gaudi::Interface::Bind::Box<IPrFitterTool>( const DetectorElement&, const ITrackExtrapolator&,
                                                        const LHCb::UniqueIDGenerator& ),
            LHCb::DetDesc::usesBaseAndConditions<LHCb::DetDesc::ConditionAccessorHolder<FixTESPath<AlgTool>>,
                                                 DetectorElement>> {
  private:
    class BoundInstance final : public Gaudi::Interface::Bind::Stub<IPrFitterTool> {
      const KalmanFilterTool*        m_parent;
      const DetectorElement&         m_geometry;
      const ITrackExtrapolator&      m_extrap;
      const LHCb::UniqueIDGenerator& m_unique_id_gen;

    public:
      BoundInstance( const KalmanFilterTool* parent, const DetectorElement& lhcb, const ITrackExtrapolator& extrap,
                     const LHCb::UniqueIDGenerator& unique_id_gen )
          : m_parent{parent}, m_geometry{lhcb}, m_extrap{extrap}, m_unique_id_gen{unique_id_gen} {}

      // these operators implement the interface!
      V1Output operator()( const Long::Tracks& tracks, const Hits<HitType::VP>& hits_vp,
                           const Hits<HitType::UT>& hits_ut, const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_pr_tracks<V1Output>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_vp, hits_ut,
                                                  hits_ft );
      }

      V1Output operator()( const Downstream::Tracks& tracks, const Hits<HitType::UT>& hits_ut,
                           const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_pr_tracks<V1Output>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_ut, hits_ft );
      }

      V1Output operator()( const Seeding::Tracks& tracks, const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_pr_tracks<V1Output>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_ft );
      }

      V1Output operator()( const Velo::Tracks& tracks, const Hits<HitType::VP>& hits_vp ) const override {
        return m_parent->fit_pr_tracks<V1Output>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_vp );
      }

      V1Output operator()( const Upstream::Tracks& tracks, const Hits<HitType::VP>& hits_vp,
                           const Hits<HitType::UT>& hits_ut ) const override {
        return m_parent->fit_pr_tracks<V1Output>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_vp, hits_ut );
      }

      V1Output operator()( const TracksV1& tracks, const Hits<HitType::VP>& hits_vp, const Hits<HitType::UT>& hits_ut,
                           const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_v1_tracks( tracks, m_geometry, m_extrap, hits_vp, hits_ut, hits_ft );
      }

      V1Output operator()( const TracksV1& tracks, const Hits<HitType::VP>& hits_vp,
                           const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_v1_tracks( tracks, m_geometry, m_extrap, hits_vp, hits_ft );
      }

      V1Output operator()( const TracksV1& tracks, const Hits<HitType::VP>& hits_vp, const Hits<HitType::UT>& hits_ut,
                           const Hits<HitType::FT>& hits_ft, const Hits<HitType::Muon>& hits_muon ) const override {
        return m_parent->fit_v1_tracks( tracks, m_geometry, m_extrap, hits_vp, hits_ut, hits_ft, hits_muon );
      }

      V1Output operator()( const TracksV1& tracks, const Hits<HitType::VP>& hits_vp, const Hits<HitType::FT>& hits_ft,
                           const Hits<HitType::Muon>& hits_muon ) const override {
        return m_parent->fit_v1_tracks( tracks, m_geometry, m_extrap, hits_vp, hits_ft, hits_muon );
      }

      V1Output operator()( const TracksV1& tracks, const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_v1_tracks( tracks, m_geometry, m_extrap, hits_ft );
      }

      V1Output operator()( const TracksV1& tracks, const Hits<HitType::VP>& hits_vp ) const override {
        return m_parent->fit_v1_tracks( tracks, m_geometry, m_extrap, hits_vp );
      }

      V1Output operator()( const TracksV1& tracks, const Hits<HitType::UT>& hits_ut,
                           const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_v1_tracks( tracks, m_geometry, m_extrap, hits_ut, hits_ft );
      }

      V1Output operator()( const TracksV1& tracks, const Hits<HitType::VP>& hits_vp,
                           const Hits<HitType::UT>& hits_ut ) const override {
        return m_parent->fit_v1_tracks( tracks, m_geometry, m_extrap, hits_vp, hits_ut );
      }

      V3FullOutput fitted_v3_tracks( const Long::Tracks& tracks, const Hits<HitType::VP>& hits_vp,
                                     const Hits<HitType::UT>& hits_ut,
                                     const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_pr_tracks<V3FullOutput>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_vp, hits_ut,
                                                      hits_ft );
      }

      V3FullOutput fitted_v3_tracks( const Downstream::Tracks& tracks, const Hits<HitType::UT>& hits_ut,
                                     const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_pr_tracks<V3FullOutput>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_ut, hits_ft );
      }

      V3FullOutput fitted_v3_tracks( const Seeding::Tracks& tracks, const Hits<HitType::FT>& hits_ft ) const override {
        return m_parent->fit_pr_tracks<V3FullOutput>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_ft );
      }

      V3FullOutput fitted_v3_tracks( const Velo::Tracks& tracks, const Hits<HitType::VP>& hits_vp ) const override {
        return m_parent->fit_pr_tracks<V3FullOutput>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_vp );
      }

      V3FullOutput fitted_v3_tracks( const Upstream::Tracks& tracks, const Hits<HitType::VP>& hits_vp,
                                     const Hits<HitType::UT>& hits_ut ) const override {
        return m_parent->fit_pr_tracks<V3FullOutput>( tracks, m_geometry, m_extrap, m_unique_id_gen, hits_vp, hits_ut );
      }

      TrackV1 operator()( const TrackV1& track, const Hits<HitType::VP>& hits_vp, const Hits<HitType::UT>& hits_ut,
                          const Hits<HitType::FT>& hits_ft ) const override {
        std::vector<Node> fitnodes;
        fitnodes.reserve( 50 );
        return m_parent->fit_v1_track( track, fitnodes, m_geometry, m_extrap, hits_vp, hits_ut, hits_ft );
      }

      TrackV1 operator()( const TrackV1& track, const Hits<HitType::VP>& hits_vp, const Hits<HitType::UT>& hits_ut,
                          const Hits<HitType::FT>& hits_ft, const Hits<HitType::Muon>& hits_muon ) const override {
        std::vector<Node> fitnodes;
        fitnodes.reserve( 50 );
        return m_parent->fit_v1_track( track, fitnodes, m_geometry, m_extrap, hits_vp, hits_ut, hits_ft, hits_muon );
      }

      TrackV1 operator()( const TrackV1& track, const Hits<HitType::FT>& hits_ft ) const override {
        std::vector<Node> fitnodes;
        fitnodes.reserve( 50 );
        return m_parent->fit_v1_track( track, fitnodes, m_geometry, m_extrap, hits_ft );
      }

      TrackV1 operator()( const TrackV1& track, const Hits<HitType::VP>& hits_vp ) const override {
        std::vector<Node> fitnodes;
        fitnodes.reserve( 50 );
        return m_parent->fit_v1_track( track, fitnodes, m_geometry, m_extrap, hits_vp );
      }

      TrackV1 operator()( const TrackV1& track, const Hits<HitType::UT>& hits_ut,
                          const Hits<HitType::FT>& hits_ft ) const override {
        std::vector<Node> fitnodes;
        fitnodes.reserve( 50 );
        return m_parent->fit_v1_track( track, fitnodes, m_geometry, m_extrap, hits_ut, hits_ft );
      }

      TrackV1 operator()( const TrackV1& track, const Hits<HitType::VP>& hits_vp,
                          const Hits<HitType::UT>& hits_ut ) const override {
        std::vector<Node> fitnodes;
        fitnodes.reserve( 50 );
        return m_parent->fit_v1_track( track, fitnodes, m_geometry, m_extrap, hits_vp, hits_ut );
      }
    };

  public:
    KalmanFilterTool( std::string type, std::string name, const IInterface* parent )
        : ToolBinder{std::move( type ),
                     std::move( name ),
                     parent,
                     {{"StandardGeometryTop", LHCb::standard_geometry_top},
                      {"ReferenceExtrapolator", "TrackMasterExtrapolator"},
                      {"UniqueIDGenerator", LHCb::UniqueIDGeneratorLocation::Default}},
                     construct<BoundInstance>( this )} {}

    // these functions are called by BoundInstance::operator()
    template <typename OutputTrackType, typename InputTrackType, HitType... Types>
    OutputTrackType fit_pr_tracks( const InputTrackType&, const DetectorElement&, const ITrackExtrapolator& extrap,
                                   const LHCb::UniqueIDGenerator& unique_id_gen, const Hits<Types>&... ) const;

    template <HitType... Types>
    TrackV1 fit_v1_track( const TrackV1&, std::vector<Node>&, const DetectorElement&, const ITrackExtrapolator& extrap,
                          const Hits<Types>&... ) const;

    template <HitType... Types>
    TracksV1 fit_v1_tracks( const TracksV1&, const DetectorElement&, const ITrackExtrapolator& extrap,
                            const Hits<Types>&... ) const;

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
    Gaudi::Property<double> m_maxchi2perdof_pre_outlier{this, "MaxChi2PreOutlierRemoval", 9999999,
                                                        "Maximum Chi2 per DoF before outlier removal"};
    Gaudi::Property<double> m_min_chi2_for_outlier{this, "MinChi2Outlier", 9,
                                                   "Minimum Chi2 of a node to be considered for outlier removal"};
    Gaudi::Property<int>    m_max_fit_iter{this, "MaxFitIterations", 10, "max number of fit iterations to perform"};
    Gaudi::Property<int> m_max_outlier_iter{this, "MaxOutlierIterations", 2, "max number of fit iterations to perform"};
    Gaudi::Property<double> m_ptVelo{this, "VeloTrackPT", 400, "PT to use when fitting VELO tracks"};
    Gaudi::Property<size_t> m_minNumVPLayers{this, "MinNumVPHitsForOutlierRemoval", 3, "Minimum number of VP layers"};
    Gaudi::Property<size_t> m_minNumUTLayers{this, "MinNumUTHitsForOutlierRemoval", 3, "Minimum number of UT layers"};
    Gaudi::Property<size_t> m_minNumFTLayers{this, "MinNumFTHitsForOutlierRemoval", 6, "Minimum number of FT layers"};
    Gaudi::Property<size_t> m_minNumMuonLayers{this, "MinNumMuonHitsForOutlierRemoval", 4,
                                               "Minimum number of Muon layers"};

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

  DECLARE_COMPONENT_WITH_ID( KalmanFilterTool, "KalmanFilterTool" )

  /**
   * @brief fit a single v1 track using the PrKalmanfilter code
   *
   * @tparam Types types of hits for the containers
   * @param track v1 track
   * @param fitnodes vector of fitnodes
   * @param lhcb our lovely detector geometry
   * @param hits needed hit containers for track type as variadic arg
   * @return TrackV1 if fit fails the Invalid flag is set
   */
  template <HitType... Types>
  TrackV1 KalmanFilterTool::fit_v1_track( const TrackV1& track, std::vector<Node>& fitnodes,
                                          const DetectorElement& lhcb, const ITrackExtrapolator& extrap,
                                          const Hits<Types>&... hits ) const {
    const auto fit_config =
        KF::FitConfiguration{m_errorX * m_errorX,    m_errorY * m_errorY,     m_errorTx * m_errorTx,
                             m_errorTy * m_errorTy,  m_errorQoP * m_errorQoP, m_maxchi2perdof_pre_outlier,
                             m_min_chi2_for_outlier, m_maxchi2perdof,         m_max_outlier_iter,
                             m_max_fit_iter,         m_minNumVPLayers,        m_minNumUTLayers,
                             m_minNumFTLayers,       m_minNumMuonLayers};

    get_fitnodes_from_lhcbids( fitnodes, track, std::forward_as_tuple( hits... ) );

    if ( !track.hasStateAt( LHCb::State::Location::FirstMeasurement ) )
      throw GaudiException{"Track doesn't have a FirstMeasurement state", "fit_v1_track", StatusCode::FAILURE};
    if ( !track.hasStateAt( LHCb::State::Location::LastMeasurement ) )
      throw GaudiException{"Track doesn't have a LastMeasurement state", "fit_v1_track", StatusCode::FAILURE};
    const auto& first = *track.stateAt( LHCb::State::Location::FirstMeasurement );
    const auto& last  = *track.stateAt( LHCb::State::Location::LastMeasurement );
    const auto  first_state_vec =
        LHCb::StateVector{{first.x(), first.y(), first.tx(), first.ty(), first.qOverP()}, first.z()};
    const auto last_state_vec = LHCb::StateVector{{last.x(), last.y(), last.tx(), last.ty(), last.qOverP()}, last.z()};

    const auto scatteringMomentum = std::clamp( track.p(), scatter_min, scatter_max );
    if ( auto tr_thru_mag = track.checkType( TrackV1::Types::Long ) || track.checkType( TrackV1::Types::Downstream ) ||
                            track.checkType( TrackV1::Types::LongMuon );
         init_nodes( fitnodes, scatteringMomentum, first_state_vec, last_state_vec, *lhcb.geometry(), extrap,
                     tr_thru_mag )
             .isFailure() ) {
      return invalid_track( TrackV1{} );
    }

    auto [chi2, success, nIter] = KF::iterate_fit( fitnodes, fit_config, *lhcb.geometry(), extrap, m_counter_iterations,
                                                   m_counter_pre_outlier_cut, m_counter_transport_failed );
    if ( !success ) { return TrackV1{}; }

    if ( !KF::remove_outliers( fitnodes, fit_config, chi2, m_counter_outlier_iterations, m_counter_cut ) ) {
      return invalid_track( TrackV1{} );
    }
    const auto new_track = make_output_track( track, fitnodes, chi2, nIter, scatteringMomentum, m_fill_fitresult,
                                              m_classic_smoothing_post, *lhcb.geometry(), extrap );
    if ( new_track.checkFlag( TrackV1::Flags::Invalid ) ) { ++m_counter_states_failed; }
    return new_track;
  }

  /**
   * @brief fit tracks from a PrTracksXXX container just like PrKalmanFilter does
   *
   * @tparam InputTrackType e.g. PrLongTracks
   * @tparam Types types of hits for containers
   * @param tracks e.g. tracks created by PrForwardTracking
   * @param lhcb your favourite detector geometry
   * @param hits hit containers for hits that are present on the track you want to fit
   * @return TracksV1 That's a KeyedContainer of pointers to v1 Tracks /-.-/
   * @note the underlying code is the same as for the algorithm PrKalmanFilter and thus gives the
   * same results
   */
  template <typename OutputTrackType, typename InputTrackType, HitType... Types>
  OutputTrackType KalmanFilterTool::fit_pr_tracks( const InputTrackType& tracks, const DetectorElement& lhcb,
                                                   const ITrackExtrapolator& extrap,
                                                   const UniqueIDGenerator&  unique_id_gen,
                                                   const Hits<Types>&... hits ) const {
    static_assert( isV1Tracks<OutputTrackType> || isV3TracksFull<OutputTrackType> );

    m_counter_tracks_in += tracks.size();

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
                             m_max_fit_iter};

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
      if ( init_nodes( fitnodes, scatteringMomentum, ctb_state_vec, downstream_state_vec, *lhcb.geometry(), extrap,
                       isLong<InputTrackType> || isDownstream<InputTrackType> )
               .isFailure() ) {
        continue;
      }

      auto [chi2, success, nIter] = KF::iterate_fit( fitnodes, fit_config, *lhcb.geometry(), extrap, iter_buffer,
                                                     pre_outlier_chi2_cut_buffer, transport_failed_buffer );
      if ( !success ) { continue; }

      if ( !KF::remove_outliers( fitnodes, fit_config, chi2, outlier_iter_buffer, chi2_cut_buffer ) ) { continue; }

      if constexpr ( isV1Tracks<OutputTrackType> ) {
        auto new_track = make_output_track( tracks, track, fitnodes, chi2, nIter, scatteringMomentum, m_fill_fitresult,
                                            m_classic_smoothing_post, *lhcb.geometry(), extrap );
        if ( !new_track ) {
          ++states_failed_buffer;
          continue;
        }
        std::get<TracksV1>( output ).add( new_track.release() );
      } else if constexpr ( isV3Tracks<OutputTrackType> ) {
        auto sc = add_output_v3_track( output, tracks, track, fitnodes, chi2, scatteringMomentum, *lhcb.geometry(),
                                       extrap, unique_id_gen );
        if ( sc.isFailure() ) {
          ++states_failed_buffer;
          continue;
        }
      } else if constexpr ( isV3TracksExtra<OutputTrackType> || isV3TracksFull<OutputTrackType> ) {
        auto& new_tracks        = std::get<TracksV3>( output );
        auto& new_partial_chi2s = std::get<PartialChiSquareds>( output );
        auto  sc = add_output_v3_track( new_tracks, tracks, track, fitnodes, chi2, scatteringMomentum, *lhcb.geometry(),
                                       extrap, unique_id_gen );
        if ( sc.isFailure() ) {
          ++states_failed_buffer;
          continue;
        }
        add_output_v3_partial_chi2( new_partial_chi2s, fitnodes, chi2 );
        if constexpr ( isV3TracksFull<OutputTrackType> ) {
          auto& new_fit_results = std::get<TrackFitResults>( output );
          new_fit_results.emplace_back(
              make_fit_result( fitnodes, nIter, scatteringMomentum, m_classic_smoothing_post ) );
        }
      }

    } // loop over tracks

    if constexpr ( isV1Tracks<OutputTrackType> ) {
      m_counter_tracks_out += std::get<TracksV1>( output ).size();
      return output;
    } else if constexpr ( isV3Tracks<OutputTrackType> ) {
      m_counter_tracks_out += std::get<TracksV3>( output ).size();
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

  /**
   * @brief fit a v1 tracks using the PrKalmanfilter code
   *
   * @tparam Types
   * @param tracks v1 track
   * @param fitnodes vector of fitnodes
   * @param lhcb our lovely detector geometry
   * @param hits needed hit containers for track type as variadic arg
   * @return TracksV1
   */
  template <HitType... Types>
  TracksV1 KalmanFilterTool::fit_v1_tracks( const TracksV1& tracks, const DetectorElement& lhcb,
                                            const ITrackExtrapolator& extrap, const Hits<Types>&... hits ) const {
    m_counter_tracks_in += tracks.size();
    TracksV1 output;
    output.reserve( tracks.size() );
    std::vector<Node> fitnodes;
    fitnodes.reserve( 50 );
    for ( const auto* track : tracks ) {
      const auto new_track = fit_v1_track( *track, fitnodes, lhcb, extrap, hits... );
      if ( new_track.checkFlag( TrackV1::Flags::Invalid ) ) continue;
      output.add( new TrackV1( std::move( new_track ) ) );
    }
    m_counter_tracks_out += output.size();
    return output;
  }

} // namespace LHCb::Pr
