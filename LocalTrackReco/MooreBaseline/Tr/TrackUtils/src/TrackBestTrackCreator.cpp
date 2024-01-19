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
#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/ChiSquare.h"
#include "Event/KalmanFitResult.h"
#include "Event/Track.h"
#include "Event/TrackFunctor.h"
#include "Kernel/HitPattern.h"
#include "Kernel/LHCbID.h"
#include "Kernel/STLExtensions.h"
#include "LHCbAlgs/MergingTransformer.h"
#include "TrackInterfaces/IGhostProbability.h"
#include "TrackInterfaces/ITrackFitter.h"
#include "TrackKernel/TrackCloneData.h"
#include "TrackKernel/TrackFunctors.h"

#include "Gaudi/Accumulators.h"
#include "GaudiKernel/ToolHandle.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

#if DEBUGHISTOGRAMS
#  include "GaudiAlg/GaudiHistoAlg.h"
using TrackBestTrackCreatorBase = LHCb::DetDesc::ConditionAccessorHolder<GaudiHistos<FixTESPath<Gaudi::Algorithm>>>;
#else
using TrackBestTrackCreatorBase = LHCb::DetDesc::ConditionAccessorHolder<FixTESPath<Gaudi::Algorithm>>;
#endif

namespace {

  /// structure to save some data for each track
  /// first parameter(false) indicates that the TrackData
  // does NOT own the track
  class TrackData : public LHCb::TrackCloneData<false> {
  private:
    bool                   m_isAccepted{false};
    double                 m_qOverP{0};
    LHCb::Track::FitStatus m_fitStatus{LHCb::Track::FitStatus::Unknown};
    enum { Clone = 1 };

  public:
    /// constructor
    TrackData( LHCb::Track& tr )
        : TrackCloneData<false>( &tr ), m_qOverP( track().firstState().qOverP() ), m_fitStatus( track().fitStatus() ) {}
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

namespace LHCb {

  /** @brief kill clones among tracks in input container and Kalman-fit survivors
   *
   * @author Wouter Hulsbergen
   * - initial release
   *
   * @author Manuel Schiller
   * @date 2014-11-18
   * - simplify, C++11 cleanups, use BloomFilter based TrackCloneData
   * @date 2014-12-09
   * - add code to use ancestor information on tracks to deal with obvious clones
   */
  class TrackBestTrackCreator final
      : public Algorithm::MergingTransformer<Tracks( Gaudi::Functional::vector_of_const_<Track::Range> const& ),
                                             Gaudi::Functional::Traits::BaseClass_t<TrackBestTrackCreatorBase>> {
  public:
    TrackBestTrackCreator( const std::string& name, ISvcLocator* pSvcLocator );

    virtual StatusCode initialize() override; ///< Algorithm initialization
    Tracks             operator()( Gaudi::Functional::vector_of_const_<Track::Range> const& ) const override;
    StatusCode         stop() override {
      m_fitter->reset();
      return StatusCode::SUCCESS;
    }

    template <typename CounterType>
    class CounterSet {
    public:
      template <typename OWNER>
      CounterSet( OWNER* owner, const std::string& tag ) {
        using Type           = Track::Types;
        constexpr auto ilast = static_cast<int>( Type::Last );
        // setup counters for each type of forward track
        for ( auto itype{0}; itype < ilast; ++itype ) {
          m_counters.emplace_back( owner, toString( static_cast<Type>( itype ) ) + "." + tag );
        }
      }
      CounterType& operator()( LHCb::Track::Types type ) { return m_counters[static_cast<int>( type )]; }

    private:
      // the actual counters
      std::deque<CounterType> m_counters;
    };

  private:
    ToolHandle<ITrackFitter>      m_fitter{this, "Fitter", "TrackMasterFitter"};
    ToolHandle<IGhostProbability> m_ghostTool{this, "GhostIdTool", "UpgradeGhostId"};

    DetDesc::ConditionAccessor<DetectorElement> m_lhcb{this, "StandardGeometryTop", standard_geometry_top};

    Gaudi::Property<double> m_maxOverlapFracVelo{this, "MaxOverlapFracVelo", 0.5};
    Gaudi::Property<double> m_maxOverlapFracT{this, "MaxOverlapFracT", 0.5};
    Gaudi::Property<double> m_maxOverlapFracUT{this, "MaxOverlapFracUT", 0.35, "essentially: max 1 common hit"};
    Gaudi::Property<double> m_minLongLongDeltaQoP{this, "MinLongLongDeltaQoP", -1};
    Gaudi::Property<double> m_minLongDownstreamDeltaQoP{this, "MinLongDownstreamDeltaQoP", 5e-6};
    Gaudi::Property<double> m_maxChi2DoF{this, "MaxChi2DoF", 3};
    Gaudi::Property<double> m_maxChi2DoFVelo{this, "MaxChi2DoFVelo", 999};
    Gaudi::Property<double> m_maxChi2DoFT{this, "MaxChi2DoFT", 999};
    Gaudi::Property<double> m_maxChi2DoFMatchAndUT{this, "MaxChi2DoFMatchUT", 999};
    Gaudi::Property<double> m_maxGhostProb{this, "MaxGhostProb", 99999};
    Gaudi::Property<bool>   m_fitTracks{this, "FitTracks", true, "fit the tracks using the Kalman filter"};
    Gaudi::Property<bool>   m_addGhostProb{this, "AddGhostProb", false, "Add the Ghost Probability to a track"};
    Gaudi::Property<bool>   m_useAncestorInfo{this, "UseAncestorInfo", true,
                                            "use ancestor information to identify obvious clones"};
    Gaudi::Property<bool>   m_doNotRefit{this, "DoNotRefit", false, "Do not refit already fitted tracks"};

    // A merging transformer does not easily allow this to be passed in the operator(), so we just do it manually and
    // pass it as a member variable.
    LHCb::DetDesc::ConditionAccessor<LHCb::Detector::DeLHCb> m_det{this, "DeLHCb", LHCb::standard_geometry_top};

    mutable Gaudi::Accumulators::Counter<>         m_fittedBeforeCounter{this, "FittedBefore"};
    mutable Gaudi::Accumulators::Counter<>         m_fitFailedBeforeCounter{this, "FitFailedBefore"};
    mutable Gaudi::Accumulators::BinomialCounter<> m_badInputCounter{this, "BadInput"};
    mutable Gaudi::Accumulators::BinomialCounter<> m_fitFailedCounter{this, "FitFailed"};

    mutable CounterSet<Gaudi::Accumulators::AveragingCounter<double>> m_chisqProbSumCounters{this, "chisqProbSum"};
    mutable CounterSet<Gaudi::Accumulators::AveragingCounter<double>> m_ghostProbCounters{this, "ghostProbability"};
    mutable CounterSet<Gaudi::Accumulators::AveragingCounter<unsigned int>> m_numOutliersCounters{this, "numOutliers"};
    mutable CounterSet<Gaudi::Accumulators::AveragingCounter<unsigned int>> m_fitFailedCounters{this, "FitFailed"};
    mutable CounterSet<Gaudi::Accumulators::BinomialCounter<>>              m_badChisqCounters{this, "badChisq"};
    mutable CounterSet<Gaudi::Accumulators::BinomialCounter<>>              m_flipChargeCounters{this, "flipCharge"};

  protected:
    /// are tracks clones in Velo
    bool veloClones( const TrackData&, const TrackData& ) const;
    /// are tracks clones in Velo
    bool veloOrClones( const TrackData&, const TrackData& ) const;
    /// are tracks clones in T
    bool TClones( const TrackData&, const TrackData& ) const;
    /// are tracks clones in UT
    bool UTClones( const TrackData&, const TrackData& ) const;

    int  fitAndSelect( TrackData&, IGeometryInfo const& geometry ) const;
    void fitAndUpdateCounters( TrackData&, IGeometryInfo const& geometry ) const;

    /// check if tracks pointed to by their TrackData objects are clones
    bool areClones( const TrackData& it, const TrackData& jt ) const;
  };

  DECLARE_COMPONENT_WITH_ID( TrackBestTrackCreator, "TrackBestTrackCreator" )

} // namespace LHCb

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
LHCb::TrackBestTrackCreator::TrackBestTrackCreator( const std::string& name, ISvcLocator* pSvcLocator )
    : MergingTransformer( name, pSvcLocator,
                          // NOTE that the algorithm behaviour is weakly dependent on the ordering of the input
                          // containers, because if the track 'quality' is the same the ordering comes from the input
                          // container ordering. The quality is related to the track type and number of LHCbIDs. See
                          // LHCBPS-1757 for more details. One should be careful to set this order deterministically.
                          // NOTE that optimising the order may improve algorithm performance -- to be checked
                          {"TracksInContainers",
                           {TrackLocation::Forward, TrackLocation::Match, TrackLocation::VeloTT,
                            TrackLocation::Downstream, TrackLocation::Tsa, TrackLocation::Velo}},
                          {"TracksOutContainer", TrackLocation::Default} ) {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode LHCb::TrackBestTrackCreator::initialize() {
  return MergingTransformer::initialize().andThen( [&] {
    if ( !m_fitTracks ) m_fitter.disable();
    if ( !m_addGhostProb ) m_ghostTool.disable();

    // Print out the user-defined settings
    if ( msgLevel( MSG::DEBUG ) )
      debug() << endmsg << "============ TrackBestTrackCreator Settings ===========" << endmsg
              << "TracksInContainers : " << getProperty( "TracksInContainers" ).toString() << endmsg
              << "TrackOutContainer  : " << getProperty( "TracksOutContainer" ).toString() << endmsg
              << "=======================================================" << endmsg << endmsg;
  } );
}

//=============================================================================
// Main execution
//=============================================================================
LHCb::Tracks LHCb::TrackBestTrackCreator::
             operator()( Gaudi::Functional::vector_of_const_<LHCb::Track::Range> const& ranges ) const {

  // This ifdef only exists because we use ConditionAccessors by hand
  // This should not be done, the Conditions should be declared at the algorithm level
  // and be used transparently via the functional framework. But MergingTransformer
  // does not allow that for now.
  // On top, we should get the context from the operator() arguments. But this
  // is also missing for the MergingTransformer. FIXME
#ifdef USE_DD4HEP
  auto lhcb = m_lhcb.get( getConditionContext( Gaudi::Hive::currentContext() ) );
#else
  auto& lhcb = m_lhcb.get( getConditionContext( Gaudi::Hive::currentContext() ) );
#endif

  auto& geometry = *lhcb.geometry();

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;
  int nGhosts = 0;

  // create pool for TrackData objects for all input tracks
  std::vector<Track>     tracksToFit;
  std::vector<TrackData> trackdatapool;
  // get total number of tracks
  size_t nTracks = std::accumulate( std::begin( ranges ), std::end( ranges ), size_t( 0 ),
                                    []( size_t sz, const Track::Range& range ) { return sz + range.size(); } );
  trackdatapool.reserve( nTracks );
  tracksToFit.reserve( nTracks );

  // generate the TrackData objects for the input tracks
  nTracks = 0;
  for ( auto& range : ranges ) {
    for ( auto& oldtr : range ) {

      const bool fitted = m_doNotRefit.value() && ( oldtr->fitStatus() == Track::FitStatus::Fitted ||
                                                    oldtr->fitStatus() == Track::FitStatus::FitFailed );
      if ( fitted && oldtr->fitStatus() == Track::FitStatus::FitFailed ) continue;

      // clone track
      auto& tr = tracksToFit.emplace_back( *oldtr );
      ++nTracks;
      // keep a record where this track came from
      tr.addToAncestors( oldtr );
      trackdatapool.emplace_back( tr );
    }
  }

  // take a vector of "references" which is much easier to sort (because less
  // data is moved around)
  std::vector<std::reference_wrapper<TrackData>> alltracks( trackdatapool.begin(), trackdatapool.end() );

  // sort them by quality
  auto qualitySort = []( const TrackData& t1, const TrackData& t2 ) { return t1 < t2; };
  std::stable_sort( alltracks.begin(), alltracks.end(), qualitySort );

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

  std::for_each( alltracks.begin(), alltracks.end(), [&]( TrackData& t ) {
    if ( !t.cloneFlag() && !isClone( t, successful_tracks ) ) {
      nGhosts += fitAndSelect( t, geometry );
      if ( t.isAccepted() ) { successful_tracks.emplace_back( t ); }
    }
  } );

  // create output container, and put selected tracks there
  Tracks tracksOutCont;
  tracksOutCont.reserve( successful_tracks.size() );
  for ( TrackData& tr : successful_tracks ) { tracksOutCont.add( new Track( std::move( tr.track() ) ) ); }

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Selected " << successful_tracks.size() << " out of " << nTracks << " tracks. Rejected " << nGhosts
            << endmsg;
  }

  return tracksOutCont;
}

int LHCb::TrackBestTrackCreator::fitAndSelect( TrackData& td, IGeometryInfo const& geometry ) const {

  fitAndUpdateCounters( td, geometry );

  Track& track = td.track();
  td.setAccepted( track.fitStatus() == Track::FitStatus::Fitted );

  // Impose tighter restrictions to composite tracks
  if ( td.isAccepted() && ( track.type() == Track::Types::Long || track.type() == Track::Types::Upstream ||
                            track.type() == Track::Types::Downstream ) ) {

    const ChiSquare& chi2 = {track.chi2(), track.nDoF()};

    // compatibility for tracks fitted with simple_fitter
    // since they don't have a fitresult
    auto const* fit = fitResult( track );

    const ChiSquare chi2T = ( fit ? fit->chi2Downstream()
                                  : ChiSquare{track.info( Track::AdditionalInfo::FitTChi2, 0 ),
                                              static_cast<int>( track.info( Track::AdditionalInfo::FitTNDoF, 0 ) )} );

    const ChiSquare chi2Velo =
        ( fit ? fit->chi2Velo()
              : ChiSquare{track.info( Track::AdditionalInfo::FitVeloChi2, 0 ),
                          static_cast<int>( track.info( Track::AdditionalInfo::FitVeloNDoF, 0 ) )} );

    // note: this includes UT hit contribution
    const ChiSquare chi2MatchAndUT = chi2 - chi2T - chi2Velo;

    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << "Track " << chi2.chi2() << " " << chi2.nDoF() << " ("
              << ( chi2.chi2PerDoF() < m_maxChi2DoF ? "1" : "0" ) << "), " << chi2T.chi2() << " " << chi2T.nDoF()
              << " (" << ( chi2T.chi2PerDoF() < m_maxChi2DoFT ? "1" : "0" ) << "), " << chi2Velo.chi2() << " "
              << chi2Velo.nDoF() << " (" << ( chi2Velo.chi2PerDoF() < m_maxChi2DoFVelo ? "1" : "0" ) << "), "
              << chi2MatchAndUT.chi2() << " " << chi2MatchAndUT.nDoF() << " ("
              << ( chi2MatchAndUT.chi2PerDoF() < m_maxChi2DoFMatchAndUT ? "1" : "0" ) << "), "
              << track.ghostProbability() << " (" << ( track.ghostProbability() < m_maxGhostProb ? "1" : "0" ) << ") "
              << endmsg;
    }

    td.setAccepted( chi2.chi2PerDoF() < m_maxChi2DoF && chi2T.chi2PerDoF() < m_maxChi2DoFT &&
                    chi2Velo.chi2PerDoF() < m_maxChi2DoFVelo && chi2MatchAndUT.chi2PerDoF() < m_maxChi2DoFMatchAndUT &&
                    ( !m_addGhostProb || track.ghostProbability() < m_maxGhostProb ) );
  }

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Track " << ( td.isAccepted() ? " accepted" : " not accepted" ) << endmsg;
  }
  return !td.isAccepted();
}

/**
 * @brief      Invokes the fitter and sets some counters.
 */
void LHCb::TrackBestTrackCreator::fitAndUpdateCounters( TrackData& td, IGeometryInfo const& geometry ) const {

  Track& track = td.track();

  // Conditions for fitting a track
  if ( m_fitTracks.value() &&
       ( !m_doNotRefit ||
         ( track.fitStatus() != Track::FitStatus::Fitted && track.fitStatus() != Track::FitStatus::FitFailed ) ) &&
       ( track.nStates() != 0 && !track.checkFlag( Track::Flags::Invalid ) ) ) {
    ( *m_fitter )( track, geometry, Tr::PID::Pion() ).ignore();
  }

  bool badinput  = false;
  bool fitfailed = false;

  if ( m_doNotRefit &&
       ( td.previousStatus() == Track::FitStatus::Fitted || td.previousStatus() == Track::FitStatus::FitFailed ) ) {
    if ( td.previousStatus() == Track::FitStatus::Fitted ) {
      ++m_fittedBeforeCounter;
      if ( m_addGhostProb ) m_ghostTool->execute( track ).ignore();
    } else {
      /// fit failed before
      /// This should always be 0 as this type is filtered out when initializing the tracks
      ++m_fitFailedBeforeCounter;
      track.setFlag( Track::Flags::Invalid, true );
    }
  } else {
    if ( track.nStates() == 0 || track.checkFlag( Track::Flags::Invalid ) ) {
      track.setFlag( Track::Flags::Invalid, true );
      badinput = true;
    } else {
      // Note: The track has already been fitted
      if ( track.fitStatus() == Track::FitStatus::Fitted ) {
        if ( m_addGhostProb ) m_ghostTool->execute( track ).ignore();
        // Update counters
        if ( track.nDoF() > 0 ) {
          double chisqProb = track.probChi2();
          m_chisqProbSumCounters( track.type() ) += chisqProb;
          m_badChisqCounters( track.type() ) += bool( chisqProb < 0.01 );
        }
        m_flipChargeCounters( track.type() ) += bool( td.qOverP() * track.firstState().qOverP() < 0 );
        m_numOutliersCounters( track.type() ) += nMeasurementsRemoved( track );
        m_ghostProbCounters( track.type() ) += track.ghostProbability();
      } else {
        track.setFlag( Track::Flags::Invalid, true );
        fitfailed = true;
        m_fitFailedCounters( track.type() ) += int( fitfailed );
      }
    }
  }
  m_badInputCounter += badinput;
  m_fitFailedCounter += fitfailed;
}

bool LHCb::TrackBestTrackCreator::veloOrClones( const TrackData& lhs, const TrackData& rhs ) const {
  const double f = lhs.overlapFraction( rhs, TrackData::VP );
#ifdef DEBUGHISTOGRAMS
  if ( f > 0 ) plot1D( f, "veloOverlapFractionH1", 0, 1 );
#endif
  return ( f > m_maxOverlapFracVelo );
}

bool LHCb::TrackBestTrackCreator::TClones( const TrackData& lhs, const TrackData& rhs ) const {
  const double f = lhs.overlapFraction( rhs, TrackData::T );
#ifdef DEBUGHISTOGRAMS
  if ( f > 0 ) plot1D( f, "TOverlapFractionH1", 0, 1 );
#endif
  return f > m_maxOverlapFracT;
}

bool LHCb::TrackBestTrackCreator::UTClones( const TrackData& lhs, const TrackData& rhs ) const {
  const double f = lhs.overlapFraction( rhs, TrackData::UT );
#ifdef DEBUGHISTOGRAMS
  if ( f > 0 ) plot1D( f, "UTOverlapFractionH1", 0, 1 );
#endif
  return f > m_maxOverlapFracUT;
}

constexpr int to_index( LHCb::Track::Types i, LHCb::Track::Types j ) {
  return static_cast<int>( i ) + 256 * static_cast<int>( j );
}

bool LHCb::TrackBestTrackCreator::areClones( const TrackData& it, const TrackData& jt ) const {
  const Track &itrack( it.track() ), &jtrack( jt.track() );
  const double dqop = it.qOverP() - jt.qOverP();
  switch ( to_index( itrack.type(), jtrack.type() ) ) {
  case to_index( Track::Types::Long, Track::Types::Long ):
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
  case to_index( Track::Types::Long, Track::Types::Downstream ):
  case to_index( Track::Types::Downstream, Track::Types::Long ):
#ifdef DEBUGHISTOGRAMS
    if ( TClones( it, jt ) ) {
      plot( dqop, "DLDqop", -2e-5, 2e-5 );
      if ( UTClones( it, jt ) ) plot( dqop, "DLDqopUTClones", -2e-5, 2e-5 );
    }
#endif
    return TClones( it, jt ) && ( std::abs( dqop ) < m_minLongDownstreamDeltaQoP || UTClones( it, jt ) );
  case to_index( Track::Types::Downstream, Track::Types::Downstream ):
    // it seems that there are no down stream tracks that share T hits ...
#ifdef DEBUGHISTOGRAMS
    if ( TClones( it, jt ) ) { plot( dqop, "DDDqop", -1e-4, 1e-4 ); }
#endif
    return TClones( it, jt ) && UTClones( it, jt );
  case to_index( Track::Types::Long, Track::Types::Upstream ):
  case to_index( Track::Types::Upstream, Track::Types::Long ):
  case to_index( Track::Types::Upstream, Track::Types::Upstream ):
    return veloOrClones( it, jt ) && UTClones( it, jt );
  case to_index( LHCb::Track::Types::Long, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::Long ):
  case to_index( LHCb::Track::Types::Upstream, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::Upstream ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Velo, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::Velo ):
  case to_index( LHCb::Track::Types::Long, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::Long ):
  case to_index( LHCb::Track::Types::Upstream, LHCb::Track::Types::VeloBackward ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::Upstream ):
  case to_index( LHCb::Track::Types::VeloBackward, LHCb::Track::Types::VeloBackward ):
    return veloOrClones( it, jt );
  case to_index( Track::Types::Long, Track::Types::Ttrack ):
  case to_index( Track::Types::Ttrack, Track::Types::Long ):
  case to_index( Track::Types::Downstream, Track::Types::Ttrack ):
  case to_index( Track::Types::Ttrack, Track::Types::Downstream ):
  case to_index( Track::Types::Ttrack, Track::Types::Ttrack ):
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
