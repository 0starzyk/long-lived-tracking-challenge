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
#include "Event/KalmanFitResult.h"
#include "Event/Track.h"
#include "Event/TrackFunctor.h"
#include "Gaudi/Accumulators.h"
#include "GaudiKernel/ToolHandle.h"
#include "Kernel/HitPattern.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Transformer.h"
#include "TrackInterfaces/ITrackCloneFinder.h"
#include "TrackKernel/TrackCloneData.h"
#include "TrackKernel/TrackFunctors.h"
#include "range/v3/version.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <exception>
#include <functional>
#include <memory>
#include <type_traits>
#include <unordered_map>
#if DEBUGHISTOGRAMS
#  include "GaudiAlg/GaudiHistoAlg.h"
using TrackCloneKillerBase = GaudiHistoAlg;
#else
using TrackCloneKillerBase = GaudiAlgorithm;
#endif

namespace {

  /// structure to save some data for each track
  class TrackData : public LHCb::TrackCloneData<false> {
  private:
    bool                   m_isAccepted{false};
    double                 m_qOverP{0};
    LHCb::Track::FitStatus m_fitStatus{LHCb::Track::FitStatus::Unknown};
    enum { Clone = 1 };

  public:
    /// constructor
    TrackData( LHCb::Track* tr )
        : TrackCloneData<false>{tr}, m_qOverP( track().firstState().qOverP() ), m_fitStatus( track().fitStatus() ) {}
    /// return q/p (or what it was at construction time)
    double                 qOverP() const { return m_qOverP; }
    LHCb::Track::FitStatus previousStatus() const { return m_fitStatus; }
    bool                   cloneFlag() const { return userFlags() & Clone; }
    void                   setCloneFlag() { setUserFlags( userFlags() | Clone ); }
    bool                   isAccepted() const { return m_isAccepted; }
    void                   setAccepted( const bool isAccepted ) { m_isAccepted = isAccepted; }

    std::vector<std::reference_wrapper<TrackData>> clones;
  };

} // namespace

/** @brief Kills clones of fitted tracks wrt to reference container and inside input container (optionally).
 *         Optionally not fitted tracks can be used as well in both input and reference containers.
 *
 *
 * @author Andrii Usachov
 * - initial release, largely copied from TrackBestTrackCreator
 */
class TrackCloneKiller final
    : public LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Tracks&, const LHCb::Tracks& ),
                                          Gaudi::Functional::Traits::BaseClass_t<TrackCloneKillerBase>> {
public:
  /// Standard constructor

  using base_class_t = LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Tracks&, const LHCb::Tracks& ),
                                                    Gaudi::Functional::Traits::BaseClass_t<TrackCloneKillerBase>>;
  //   using base_class_t::addConditionDerivation;
  using base_class_t::debug;
  using base_class_t::error;
  using base_class_t::info;
  using base_class_t::inputLocation;
  using base_class_t::msgLevel;

  TrackCloneKiller( const std::string& name, ISvcLocator* pSvcLocator );

  virtual StatusCode initialize() override; ///< Algorithm initialization
  LHCb::Tracks       operator()( const LHCb::Tracks& inTracks, const LHCb::Tracks& refTracks ) const override;

private:
  Gaudi::Property<double> m_maxOverlapFracVelo{this, "MaxOverlapFracVelo", 0.5};
  Gaudi::Property<double> m_maxOverlapFracT{this, "MaxOverlapFracT", 0.5};
  Gaudi::Property<double> m_maxOverlapFracUT{this, "MaxOverlapFracUT", 0.35, "essentially: max 1 common hit"};
  Gaudi::Property<double> m_minLongLongDeltaQoP{this, "MinLongLongDeltaQoP", -1};
  Gaudi::Property<double> m_minLongDownstreamDeltaQoP{this, "MinLongDownstreamDeltaQoP", 5e-6};
  Gaudi::Property<bool>   m_keepUnFitted{this, "KeepUnFitted", false, "Keep unfitted tracks"};
  Gaudi::Property<bool>   m_useUnFittedRef{this, "UseUnFittedRef", false, "Use unfitted tracks from reference"};
  Gaudi::Property<bool>   m_skipSameContainerTracks{this, "SkipSameContainerTracks", true};

protected:
  std::vector<TrackData> fillDataPool( const LHCb::Tracks& inTracks, bool keepUnFitted = false ) const;

  /// are tracks clones in Velo
  bool veloClones( const TrackData&, const TrackData& ) const;
  /// are tracks clones in Velo
  bool veloOrClones( const TrackData&, const TrackData& ) const;

  /// are tracks clones in T
  bool TClones( const TrackData&, const TrackData& ) const;
  /// are tracks clones in UT
  bool UTClones( const TrackData&, const TrackData& ) const;

  /// check if tracks pointed to by their TrackData objects are clones
  bool areClones( const TrackData& it, const TrackData& jt ) const;

  /// mapping between original track and the index of its copy
  using CopyMapEntry = std::pair<const LHCb::Track*, size_t>;
  using CopyMap      = std::vector<CopyMapEntry>;
};

DECLARE_COMPONENT( TrackCloneKiller )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TrackCloneKiller::TrackCloneKiller( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator, {KeyValue{"TracksInContainer", {}}, KeyValue{"TracksRefContainer", {}}},
                   {"TracksOutContainer", LHCb::TrackLocation::Default} ) {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode TrackCloneKiller::initialize() {
  return Transformer::initialize().andThen( [&] {
    // Print out the user-defined settings
    if ( msgLevel( MSG::DEBUG ) )
      debug() << endmsg << "============ TrackCloneKiller Settings ===========" << endmsg
              << "TracksInContainer : " << getProperty( "TracksInContainer" ).toString() << endmsg
              << "TrackOutContainer  : " << getProperty( "TracksOutContainer" ).toString() << endmsg
              << "=======================================================" << endmsg << endmsg;
  } );
}

std::vector<TrackData> TrackCloneKiller::fillDataPool( const LHCb::Tracks& inTracks, bool keepUnFitted ) const {
  // create pool for TrackData objects for all input tracks
  std::vector<TrackData> trackdatapool;

  // reserve enough space so we don't have to reallocate
  trackdatapool.reserve( inTracks.size() );

  // generate the TrackData objects for the input tracks, initialising the
  // States for use in the Kalman filter on the way
  for ( auto& oldtr : inTracks ) {
    // pre-initialise (if required)
    if ( !keepUnFitted ) {
      const bool fitted = ( oldtr->fitStatus() == LHCb::Track::FitStatus::Fitted ||
                            oldtr->fitStatus() == LHCb::Track::FitStatus::FitFailed );
      if ( !fitted || oldtr->fitStatus() == LHCb::Track::FitStatus::FitFailed ) continue;
    }

    // keep a record where this track came from
    trackdatapool.emplace_back( oldtr );
  }
  return trackdatapool;
}

//=============================================================================
// Main execution
//=============================================================================
LHCb::Tracks TrackCloneKiller::operator()( const LHCb::Tracks& inTracks, const LHCb::Tracks& refTracks ) const {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  // take a vector of "references" which is much easier to sort (because less
  // data is moved around)
  std::vector<TrackData>                         trackdatapool = fillDataPool( inTracks, m_keepUnFitted.value() );
  std::vector<std::reference_wrapper<TrackData>> alltracks( trackdatapool.begin(), trackdatapool.end() );

  // sort them by quality
  auto qualitySort = []( const TrackData& t1, const TrackData& t2 ) { return t1 < t2; };
  std::stable_sort( alltracks.begin(), alltracks.end(), qualitySort );

  std::vector<TrackData>                         refdatapool = fillDataPool( refTracks, m_useUnFittedRef.value() );
  std::vector<std::reference_wrapper<TrackData>> ref_tracks( refdatapool.begin(), refdatapool.end() );

  // Prepare TrackData reference containers
  std::vector<std::reference_wrapper<TrackData>> successful_tracks;

  // Helper function to verify if t is a clone

  // Conditions for being a clone:
  // * cloneFlag is true
  // * It's a clone of any track in v
  auto isClone = [&]( TrackData& t, const std::vector<std::reference_wrapper<TrackData>>& v ) {
    const auto firstClone =
        std::find_if( v.begin(), v.end(), [&]( const TrackData& t2 ) { return areClones( t, t2 ); } );
    if ( firstClone != v.end() ) return true;
    return false;
  };

  // Sequential treatment
  if ( m_skipSameContainerTracks.value() )
    std::for_each( alltracks.begin(), alltracks.end(), [&]( TrackData& t ) {
      if ( !t.cloneFlag() && !isClone( t, ref_tracks ) ) { successful_tracks.emplace_back( t ); }
    } );
  else
    std::for_each( alltracks.begin(), alltracks.end(), [&]( TrackData& t ) {
      if ( !t.cloneFlag() && !isClone( t, ref_tracks ) && !isClone( t, successful_tracks ) ) {
        successful_tracks.emplace_back( t );
      }
    } );

  // create output container, and put selected tracks there
  LHCb::Tracks tracksOutCont;
  tracksOutCont.reserve( successful_tracks.size() );

  // insert selected tracks
  for ( TrackData& tr : successful_tracks ) { tracksOutCont.add( std::move( tr ).trackptr() ); }

  if ( msgLevel( MSG::DEBUG ) ) {
    size_t nTracks = inTracks.size();
    size_t nClones = nTracks - successful_tracks.size();
    debug() << "Selected " << successful_tracks.size() << " out of " << nTracks << " tracks. Rejected  clones"
            << nClones << endmsg;
  }

  return tracksOutCont;
}

bool TrackCloneKiller::veloOrClones( const TrackData& lhs, const TrackData& rhs ) const {
  const double f = lhs.overlapFraction( rhs, TrackData::VP );
#ifdef DEBUGHISTOGRAMS
  if ( f > 0 ) plot1D( fR, "veloOverlapFractionH1", 0, 1 );
#else
  return f > m_maxOverlapFracVelo;
#endif
}

bool TrackCloneKiller::TClones( const TrackData& lhs, const TrackData& rhs ) const {
  const double f = lhs.overlapFraction( rhs, TrackData::T );
#ifdef DEBUGHISTOGRAMS
  if ( f > 0 ) plot1D( f, "TOverlapFractionH1", 0, 1 );
#endif
  return f > m_maxOverlapFracT;
}

bool TrackCloneKiller::UTClones( const TrackData& lhs, const TrackData& rhs ) const {
  const double f = lhs.overlapFraction( rhs, TrackData::UT );
#ifdef DEBUGHISTOGRAMS
  if ( f > 0 ) plot1D( f, "UTOverlapFractionH1", 0, 1 );
#endif
  return f > m_maxOverlapFracUT;
}

constexpr int to_index( LHCb::Track::Types i, LHCb::Track::Types j ) {
  constexpr int offset = 256;
  return static_cast<int>( i ) + offset * static_cast<int>( j );
}

bool TrackCloneKiller::areClones( const TrackData& it, const TrackData& jt ) const {
  const LHCb::Track &itrack( it.track() ), &jtrack( jt.track() );
  const double       dqop = it.qOverP() - jt.qOverP();
  switch ( to_index( itrack.type(), jtrack.type() ) ) {
  case to_index( LHCb::Track::Types::Long, LHCb::Track::Types::Long ):
#ifdef DEBUGHISTOGRAMS
    if ( TClones( it, jt ) && veloOrClones( it, jt ) ) {
      plot( dqop, "LLDqopClones", -1e-5, 1e-5 );
    } else if ( TClones( it, jt ) ) {
      plot( dqop, "LLDqopTClones", -1e-5, 1e-5 );
    } else if ( veloOrClones( it, jt ) ) {
      plot( dqop, "LLDqopVeloOrClones", -1e-5, 1e-5 );
    }
#endif
    return TClones( it, jt ) && ( std::abs( dqop ) < m_minLongLongDeltaQoP || veloOrClones( it, jt ) );
  case to_index( LHCb::Track::Types::Long, LHCb::Track::Types::Downstream ):
  case to_index( LHCb::Track::Types::Downstream, LHCb::Track::Types::Long ):
#ifdef DEBUGHISTOGRAMS
    if ( TClones( it, jt ) ) {
      plot( dqop, "DLDqop", -2e-5, 2e-5 );
      if ( UTClones( it, jt ) ) plot( dqop, "DLDqopUTClones", -2e-5, 2e-5 );
    }
#endif
    return TClones( it, jt ) && ( std::abs( dqop ) < m_minLongDownstreamDeltaQoP || UTClones( it, jt ) );
  case to_index( LHCb::Track::Types::Downstream, LHCb::Track::Types::Downstream ):
    // it seems that there are no down stream tracks that share T hits ...
#ifdef DEBUGHISTOGRAMS
    if ( TClones( it, jt ) ) { plot( dqop, "DDDqop", -1e-4, 1e-4 ); }
#endif
    return TClones( it, jt ) && UTClones( it, jt );
  case to_index( LHCb::Track::Types::Long, LHCb::Track::Types::Upstream ):
  case to_index( LHCb::Track::Types::Upstream, LHCb::Track::Types::Long ):
  case to_index( LHCb::Track::Types::Upstream, LHCb::Track::Types::Upstream ):
    return veloOrClones( it, jt ) && UTClones( it, jt );
  case to_index( LHCb::Track::Types::Long, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::Long ):
  case to_index( LHCb::Track::Types::Upstream, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::Upstream ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Long, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::Long ):
  case to_index( LHCb::Track::Types::Upstream, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::Upstream ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::Velo ):
    return veloOrClones( it, jt );
  case to_index( LHCb::Track::Types::Long, LHCb::Track::Types::Ttrack ):
  case to_index( LHCb::Track::Types::Ttrack, LHCb::Track::Types::Long ):
  case to_index( LHCb::Track::Types::Downstream, LHCb::Track::Types::Ttrack ):
  case to_index( LHCb::Track::Types::Ttrack, LHCb::Track::Types::Downstream ):
  case to_index( LHCb::Track::Types::Ttrack, LHCb::Track::Types::Ttrack ):
    return TClones( it, jt );
  case to_index( LHCb::Track::Types::Ttrack, LHCb::Track::Types::Upstream ):
  case to_index( LHCb::Track::Types::Upstream, LHCb::Track::Types::Ttrack ):
  case to_index( LHCb::Track::Types::Ttrack, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::Ttrack ):
  case to_index( LHCb::Track::Types::Downstream, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::Downstream ):
  case to_index( LHCb::Track::Types::Ttrack, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::Ttrack ):
  case to_index( LHCb::Track::Types::Downstream, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::Downstream ):
  case to_index( LHCb::Track::Types::Downstream, LHCb::Track::Types::Upstream ):
  case to_index( LHCb::Track::Types::Upstream, LHCb::Track::Types::Downstream ):
    break;
  default:
    error() << "Don't know how to handle combi: " << itrack.type() << " " << jtrack.type() << endmsg;
  }
  return false;
}
