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
#include "Event/GhostProbability.h"
#include "Event/PrDownstreamTracks.h"
#include "Event/PrLongTracks.h"
#include "Event/PrSeedTracks.h"
#include "Event/Track.h"
#include "Event/Track_v1.h"
#include "Event/Track_v3.h"
#include "Gaudi/Accumulators.h"
#include "Gaudi/PluginServiceV2.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Transformer.h"
#include "LHCbMath/SIMDWrapper.h"
#include "TrackKernel/TrackCloneData.h"
#include <algorithm>
#include <functional>
#include <stdexcept>

namespace {

  struct SimpleTrackData : public LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSliceWithLHCbIDs<> {
    SimpleTrackData( std::vector<LHCb::LHCbID> const& trackids, float qp ) : m_qOverP( qp ) {
      for ( auto it = trackids.cbegin(); trackids.cend() != it; ++it ) {
        switch ( it->detectorType() ) {
        case LHCb::LHCbID::channelIDtype::VP: {
          auto jt = it + 1;
          for ( ; trackids.cend() != jt && ( LHCb::LHCbID::channelIDtype::VP == jt->detectorType() ); ++jt ) {}
          m_ids[HitType::VP] = LHCbIDs( it, jt );
          this->m_bloomfilters[HitType::VP].insert( m_ids[HitType::VP] );
          it = --jt;
        } break;
        case LHCb::LHCbID::channelIDtype::UT: {
          auto jt = it + 1;
          for ( ; trackids.cend() != jt && ( LHCb::LHCbID::channelIDtype::UT == jt->detectorType() ); ++jt ) {}
          m_ids[HitType::UT] = LHCbIDs( it, jt );
          this->m_bloomfilters[HitType::UT].insert( m_ids[HitType::UT] );
          it = --jt;
        } break;
        case LHCb::LHCbID::channelIDtype::FT: {
          auto jt = it + 1;
          for ( ; trackids.cend() != jt && ( LHCb::LHCbID::channelIDtype::FT == jt->detectorType() ); ++jt ) {}
          m_ids[HitType::T] = LHCbIDs( it, jt );
          this->m_bloomfilters[HitType::T].insert( m_ids[HitType::T] );
          it = --jt;
        } break;
        default:
          break;
        }
      }
    }

    SimpleTrackData( LHCb::Event::Track const* tr ) : SimpleTrackData( tr->lhcbIDs(), tr->firstState().qOverP() ) {}

    float qOverP() const { return m_qOverP; }
    float m_qOverP{0};
  };

} // namespace

/** @brief Kills clones of fitted tracks wrt to reference container.
 */
template <typename InTracksType, typename RefTracksType>
class PrCloneKiller final
    : public LHCb::Algorithm::Transformer<InTracksType( InTracksType const&, RefTracksType const& )> {
public:
  using base_class_t = LHCb::Algorithm::Transformer<InTracksType( InTracksType const&, RefTracksType const& )>;
  using TrackType    = LHCb::Track::Types;

  PrCloneKiller( const std::string& name, ISvcLocator* pSvcLocator )
      : base_class_t( name, pSvcLocator,
                      {typename base_class_t::KeyValue{"TracksInContainer", ""},
                       typename base_class_t::KeyValue{"TracksRefContainer", ""}},
                      typename base_class_t::KeyValue{"TracksOutContainer", ""} ) {}

  InTracksType operator()( const InTracksType& inTracks, const RefTracksType& refTracks ) const override {

    auto const input_tracks = inTracks.scalar();

    if ( refTracks.empty() ) {
      m_counter_input += input_tracks.size();
      m_counter_selected += input_tracks.size();
      return input_tracks.filter( []( auto const& ) { return true; } );
    }

    std::vector<SimpleTrackData> refdatapool;
    // -- This is only needed in case we have v3::Tracks, but it has to be declared outside the 'if' clause:
    // -- SimpleTrackData contains iterators to the LHCbIDs, and if it is declared inside the clause, they become
    // -- invalid.
    std::vector<std::vector<LHCb::LHCbID>> idsCollection;

    if constexpr ( std::is_same_v<RefTracksType, LHCb::Tracks> ) {
      refdatapool = std::vector<SimpleTrackData>{refTracks.begin(), refTracks.end()};
    } else {
      // -- This is for  v3::Tracks and zips with it
      for ( const auto& track : refTracks.scalar() ) { idsCollection.push_back( std::move( track.lhcbIDs() ) ); }
      int i = 0;
      for ( const auto& id_vec : idsCollection ) {
        auto const t = refTracks.scalar()[i];
        refdatapool.emplace_back( id_vec, t.qOverP( t.defaultState() ).cast() );
        ++i;
      }
    }

    // those are the only ones needed for a light sequence so let's keep it simple
    constexpr bool islong       = std::is_same_v<InTracksType, LHCb::Pr::Long::Tracks>;
    constexpr bool isdownstream = std::is_same_v<InTracksType, LHCb::Pr::Downstream::Tracks>;
    constexpr bool isupstream   = std::is_same_v<InTracksType, LHCb::Pr::Upstream::Tracks>;
    constexpr bool isT          = std::is_same_v<InTracksType, LHCb::Pr::Seeding::Tracks>;
    constexpr bool isTrackV3    = std::is_same_v<InTracksType, LHCb::Event::v3::Tracks>;
    static_assert( islong || isdownstream || isupstream || isT || isTrackV3 );

    auto const in_type = [&]() {
      if constexpr ( isTrackV3 ) {
        return inTracks.type();
      } else if constexpr ( islong ) {
        return TrackType::Long;
      } else if constexpr ( isdownstream ) {
        return TrackType::Downstream;
      } else if constexpr ( isupstream ) {
        return TrackType::Upstream;
      } else if constexpr ( isT ) {
        return TrackType::Ttrack;
      } else {
        throw GaudiException( "Only Long, Downstream, Upstream and T tracks are supported.", "PrCloneKiller",
                              StatusCode::FAILURE );
      }
    }();

    auto ref_type = TrackType::Long;
    if constexpr ( std::is_same_v<RefTracksType, LHCb::Tracks> ) {
      ref_type = ( *refTracks.begin() )->type();
    } else {
      ref_type = refTracks.scalar()[0].type();
    }

    // this algorithm is implemented with the assumption that all tracks
    // in the ref container are of the same type
    // inTracks can only be one type as it's a PrTracks container
    // For v3 tracks, so far they can only be of one type
    if constexpr ( std::is_same_v<RefTracksType, LHCb::Tracks> ) {
      assert( std::all_of( refTracks.begin(), refTracks.end(),
                           [ref_type]( auto const* tr ) { return tr->type() == ref_type; } ) );
    }

    auto const areClones = getCloneLambda( in_type, ref_type );

    std::vector<LHCb::LHCbID> id_vec;
    id_vec.reserve( LHCb::Pr::Long::Tracks::MaxLHCbIDs );

    auto const qopLambda = []( auto const& track ) {
      if constexpr ( isTrackV3 ) {
        return track.qOverP( track.defaultState() ).cast();
      } else {
        return track.qOverP().cast();
      }
    };

    auto const pred = [&]( auto const& track ) {
      id_vec.clear();
      auto const qop = qopLambda( track );
      id_vec         = track.lhcbIDs();
      std::sort( id_vec.begin(), id_vec.end() );

      auto const t = SimpleTrackData{id_vec, qop};

      return std::none_of( refdatapool.begin(), refdatapool.end(),
                           [&]( SimpleTrackData const& t2 ) { return areClones( t, t2 ); } );
    };

    auto output = input_tracks.filter( pred );
    m_counter_input += input_tracks.size();
    m_counter_selected += output.size();
    return output;
  }

  std::function<bool( SimpleTrackData const&, SimpleTrackData const& )> getCloneLambda( TrackType a,
                                                                                        TrackType b ) const {

    constexpr auto track_combi = []( TrackType a, TrackType b ) {
      return static_cast<int>( a ) + 256 * static_cast<int>( b );
    };

    switch ( track_combi( a, b ) ) {
    case track_combi( TrackType::Long, TrackType::Long ):
      return [&]( SimpleTrackData const& it, SimpleTrackData const& jt ) {
        const auto dqop = it.qOverP() - jt.qOverP();
        return it.overlapFraction( jt, SimpleTrackData::T ) > m_maxOverlapFracT &&
               ( std::abs( dqop ) < m_minLongLongDeltaQoP ||
                 it.overlapFraction( jt, SimpleTrackData::VP ) > m_maxOverlapFracVelo );
      };

    case track_combi( TrackType::Long, TrackType::Downstream ):
    case track_combi( TrackType::Downstream, TrackType::Long ):
      return [&]( SimpleTrackData const& it, SimpleTrackData const& jt ) {
        const auto dqop = it.qOverP() - jt.qOverP();
        return it.overlapFraction( jt, SimpleTrackData::T ) > m_maxOverlapFracT &&
               ( std::abs( dqop ) < m_minLongDownstreamDeltaQoP ||
                 it.overlapFraction( jt, SimpleTrackData::UT ) > m_maxOverlapFracUT );
      };

    case track_combi( TrackType::Downstream, TrackType::Downstream ):
      // it seems that there are no down stream tracks that share T hits ...
      return [&]( SimpleTrackData const& it, SimpleTrackData const& jt ) {
        return it.overlapFraction( jt, SimpleTrackData::T ) > m_maxOverlapFracT &&
               it.overlapFraction( jt, SimpleTrackData::UT ) > m_maxOverlapFracUT;
      };

    case track_combi( TrackType::Long, TrackType::Upstream ):
    case track_combi( TrackType::Upstream, TrackType::Long ):
    case track_combi( TrackType::Upstream, TrackType::Upstream ):
      return [&]( SimpleTrackData const& it, SimpleTrackData const& jt ) {
        return it.overlapFraction( jt, SimpleTrackData::VP ) > m_maxOverlapFracVelo &&
               it.overlapFraction( jt, SimpleTrackData::UT ) > m_maxOverlapFracUT;
      };

    case track_combi( TrackType::Long, TrackType::Velo ):
      return [&]( SimpleTrackData const& it, SimpleTrackData const& jt ) {
        return it.overlapFraction( jt, SimpleTrackData::VP ) > m_maxOverlapFracVelo;
      };

    case track_combi( TrackType::Long, TrackType::Ttrack ):
    case track_combi( TrackType::Ttrack, TrackType::Long ):
    case track_combi( TrackType::Downstream, TrackType::Ttrack ):
      return [&]( SimpleTrackData const& it, SimpleTrackData const& jt ) {
        return it.overlapFraction( jt, SimpleTrackData::T ) > m_maxOverlapFracT;
      };

    case track_combi( TrackType::Upstream, TrackType::Downstream ):
    case track_combi( TrackType::Downstream, TrackType::Upstream ):
      return [&]( SimpleTrackData const&, SimpleTrackData const& ) { return false; };

    default:
      return [&]( SimpleTrackData const&, SimpleTrackData const& ) {
        throw std::runtime_error( "PrCloneKiller unhandeld switch case" );
        return false;
      };
    }
  }

private:
  Gaudi::Property<double> m_maxOverlapFracVelo{this, "MaxOverlapFracVelo", 0.5};
  Gaudi::Property<double> m_maxOverlapFracT{this, "MaxOverlapFracT", 0.5};
  Gaudi::Property<double> m_maxOverlapFracUT{this, "MaxOverlapFracUT", 0.35, "essentially: max 1 common hit"};
  Gaudi::Property<double> m_minLongLongDeltaQoP{this, "MinLongLongDeltaQoP", -1};
  Gaudi::Property<double> m_minLongDownstreamDeltaQoP{this, "MinLongDownstreamDeltaQoP", 5e-6};
  Gaudi::Property<bool>   m_internalKilling{this, "InternalKilling", false};

  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_counter_input{this, "nTracksInput"};
  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_counter_selected{this, "nTracksSelected"};
};

using PrCloneKillerLongV1 = PrCloneKiller<LHCb::Pr::Long::Tracks, LHCb::Tracks>;
using PrCloneKillerDownV1 = PrCloneKiller<LHCb::Pr::Downstream::Tracks, LHCb::Tracks>;
using PrCloneKillerUpV1   = PrCloneKiller<LHCb::Pr::Upstream::Tracks, LHCb::Tracks>;
using PrCloneKillerSeedV1 = PrCloneKiller<LHCb::Pr::Seeding::Tracks, LHCb::Tracks>;
using PrCloneKillerLongV3 = PrCloneKiller<LHCb::Pr::Long::Tracks, LHCb::Event::v3::Tracks>;
using PrCloneKillerDownV3 = PrCloneKiller<LHCb::Pr::Downstream::Tracks, LHCb::Event::v3::Tracks>;
using PrCloneKillerUpV3   = PrCloneKiller<LHCb::Pr::Upstream::Tracks, LHCb::Event::v3::Tracks>;
using PrCloneKillerSeedV3 = PrCloneKiller<LHCb::Pr::Seeding::Tracks, LHCb::Event::v3::Tracks>;

DECLARE_COMPONENT_WITH_ID( PrCloneKillerLongV1, "PrCloneKillerLong" )
DECLARE_COMPONENT_WITH_ID( PrCloneKillerDownV1, "PrCloneKillerDown" )
DECLARE_COMPONENT_WITH_ID( PrCloneKillerUpV1, "PrCloneKillerUp" )
DECLARE_COMPONENT_WITH_ID( PrCloneKillerSeedV1, "PrCloneKillerSeed" )
DECLARE_COMPONENT_WITH_ID( PrCloneKillerLongV3, "PrCloneKillerLongV3" )
DECLARE_COMPONENT_WITH_ID( PrCloneKillerDownV3, "PrCloneKillerDownV3" )
DECLARE_COMPONENT_WITH_ID( PrCloneKillerUpV3, "PrCloneKillerUpV3" )
DECLARE_COMPONENT_WITH_ID( PrCloneKillerSeedV3, "PrCloneKillerSeedV3" )
