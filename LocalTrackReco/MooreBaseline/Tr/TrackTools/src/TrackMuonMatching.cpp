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

#include "TrackInterfaces/ITrackChi2Calculator.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include "MuonDet/DeMuonDetector.h"
#include "MuonDet/MuonNamespace.h"

#include "Event/FitNode.h"
#include "Event/State.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Event/TrackParameters.h"

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "DetDesc/IDetectorElement.h"

#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Transformer.h"

#include <functional>
#include <string>
#include <vector>

// helper using to design the base detectorElementClass, as it's different
// for DDDB and DD4hep
#ifdef USE_DD4HEP
using GenericDetElem = LHCb::Detector::DeIOV;
#else
using GenericDetElem = ::DetectorElement;
#endif

/**
 *  (based on OTMuonMatching by Jan Amoraal)
 *  @author Stefania Vecchi
 *  @date   2010-06-04
 */
class TrackMuonMatching
    : public LHCb::Algorithm::Transformer<
          LHCb::Tracks( const LHCb::Tracks&, const LHCb::Tracks&, const GenericDetElem&, const DeMuonDetector& ),
          LHCb::DetDesc::usesBaseAndConditions<GaudiTupleAlg, GenericDetElem, DeMuonDetector>> {

public:
  TrackMuonMatching( const std::string& name, ISvcLocator* pSvcLocator );

  LHCb::Tracks operator()( const LHCb::Tracks& longTracks, const LHCb::Tracks& muonTracks, const GenericDetElem& geo,
                           const DeMuonDetector& ) const override;

private:
  StatusCode matchChi2( LHCb::State& longTrack, LHCb::State& mTrack, const double& atZ, double& chi2,
                        GenericDetElem const& detelem ) const;
  StatusCode longTmuonExtrap( LHCb::State* lState, const double& atZ, GenericDetElem const& detelem ) const;

  auto createMatchedTrack( LHCb::Track& longt, LHCb::Track& muont ) const;

  Gaudi::Property<double> m_matchAtZ{this, "MatchAtZ", 12500 * Gaudi::Units::mm};
  Gaudi::Property<bool>   m_matchAtFirstMuonHit{this, "MatchAtFirstMuonHit", false};
  Gaudi::Property<double> m_matchChi2Cut{this, "MatchChi2Cut", 100.0};
  Gaudi::Property<bool>   m_allCombinations{this, "AllCombinations", true};
  Gaudi::Property<bool>   m_returnLongMuon{this, "ReturnLongMuon", true};

  ToolHandle<ITrackExtrapolator>   m_extrapolator{this, "Extrapolator", "TrackLinearExtrapolator"};
  ToolHandle<ITrackChi2Calculator> m_chi2Calculator{this, "Chi2Calculator", "TrackChi2Calculator"};

  mutable Gaudi::Accumulators::Counter<> m_countMatchedTracks{this, "nMatchedTracks"};
};

using namespace LHCb;

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TrackMuonMatching )

auto TrackMuonMatching::createMatchedTrack( LHCb::Track& longTrack, LHCb::Track& muonTrack ) const {
  auto matchedTrack = std::make_unique<LHCb::Track>();
  matchedTrack->copy( longTrack );
  // Remove long track's LastMeasurement state
  // or convert it in an AtT if there isn't one yet
  if ( matchedTrack->hasStateAt( LHCb::State::Location::AtT ) ) {
    matchedTrack->removeFromStates( matchedTrack->stateAt( State::LastMeasurement ) );
  } else {
    assert( matchedTrack->stateAt( State::LastMeasurement )->z() >= StateParameters::ZBegT &&
            matchedTrack->stateAt( State::LastMeasurement )->z() <= StateParameters::ZEndT );
    matchedTrack->stateAt( State::LastMeasurement )->setLocation( LHCb::State::Location::AtT );
  }
  // Add LastMeasurement from muon track
  // and set the momentum of state AtT from long track
  if ( muonTrack.hasStateAt( LHCb::State::Location::LastMeasurement ) ) {
    matchedTrack->addToStates( *( muonTrack.stateAt( State::LastMeasurement ) ) );
    matchedTrack->stateAt( State::LastMeasurement )->setQOverP( matchedTrack->stateAt( State::AtT )->qOverP() );
  }
  /// Add muon ids to copied T track
  for ( LHCbID id : muonTrack.lhcbIDs() ) matchedTrack->addToLhcbIDs( id );
  // FIXME returning LongMuon type should be standard case
  if ( m_returnLongMuon.value() ) matchedTrack->setType( LHCb::Track::Types::LongMuon );

  return matchedTrack;
}

TrackMuonMatching::TrackMuonMatching( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer(
          name, pSvcLocator,
          {KeyValue{"LongTracksLocation", TrackLocation::Default}, KeyValue{"MuonTracksLocation", TrackLocation::Muon},
           KeyValue{"StandardGeometryTop", "/dd/Structure/LHCb"}, KeyValue{"DeMuonDetector", DeMuonLocation::Default}},
          KeyValue{"TracksOutputLocation", "Rec/Track/Best/TMuon"} ) {}

Tracks TrackMuonMatching::operator()( const Tracks& longTracks, const Tracks& muonTracks, const GenericDetElem& geo,
                                      const DeMuonDetector& ) const {
  Tracks matchedTracks;

  bool flaglongT = false;
  if ( longTracks.size() != 0 ) flaglongT = true;

  std::vector<int> mStation;
  mStation.reserve( muonTracks.size() * longTracks.size() );
  std::vector<int> mRegion;
  mRegion.reserve( muonTracks.size() * longTracks.size() );
  std::vector<double> matchChi2s;
  matchChi2s.reserve( muonTracks.size() * longTracks.size() );
  std::vector<State> matchedLStates;
  matchedLStates.reserve( muonTracks.size() * longTracks.size() );
  std::vector<State> matchedMStates;
  matchedMStates.reserve( muonTracks.size() * longTracks.size() );
  std::vector<State*> bestlState;
  matchedMStates.reserve( muonTracks.size() );
  std::vector<double> bestchi2;
  bestchi2.reserve( muonTracks.size() );

  int    i = -1;
  double z = m_matchAtZ;

  /// Now match this T track to muon tracks

  for ( Tracks::const_iterator m = muonTracks.begin(), mEnd = muonTracks.end(); m != mEnd; ++m ) {
    debug() << " MuonTrack chi2 " << ( *m )->chi2PerDoF() << endmsg;
    if ( ( *m )->chi2PerDoF() > 5. ) continue;

    State* lState = 0;
    State* mState = 0;
    i++;
    bestlState.push_back( NULL );
    double minchi2 = 10000;

    if ( flaglongT ) {
      if ( m_matchAtFirstMuonHit.value() ) {
        if ( ( *m )->hasStateAt( State::Muon ) )
          mState = ( *m )->stateAt( State::Muon );
        else if ( ( *m )->hasStateAt( State::FirstMeasurement ) )
          mState = ( *m )->stateAt( State::FirstMeasurement );

        z = mState->z();
        verbose() << "Found muon state. Going to extrapolate to this state with z = " << z << endmsg;
      }
      /// Matched Muon-T track
      auto best_matchedTrack = std::make_unique<LHCb::Track>();
      for ( Tracks::const_iterator t = longTracks.begin(), tEnd = longTracks.end(); t != tEnd; ++t ) {
        if ( !( *t )->hasT() ) continue;
        if ( ( *t )->chi2PerDoF() > 5. ) continue;
        if ( ( *t )->ghostProbability() > 0.7 ) continue;
        if ( !( *t )->checkType( LHCb::Track::Types::Long ) ) continue;

        double chi2 = -9999.0;
        /// Get the longTrack state closest to this z
        lState = &( *t )->closestState( z );
        /// Get the Muon state closest to this z
        mState = &( *m )->closestState( z );
        //  Calculate mach chi2
        StatusCode sc = matchChi2( *lState, *mState, z, chi2, geo );

        if ( sc.isSuccess() && chi2 > -1.0 && chi2 < m_matchChi2Cut ) {
          debug() << "chi2 Matching is " << chi2 << endmsg;

          matchChi2s.push_back( chi2 );
          lState->setLocation( State::Muon ); // Muon state added to long track here
          mState->setLocation( State::Muon );
          matchedLStates.push_back( *lState );
          matchedMStates.push_back( *mState );
          mStation.push_back( ( *m )->lhcbIDs().front().muonID().station() );
          mRegion.push_back( ( *m )->lhcbIDs().front().muonID().region() );

          if ( m_allCombinations.value() ) { matchedTracks.insert( createMatchedTrack( *( *t ), *( *m ) ).release() ); }

          if ( chi2 < minchi2 ) {
            bestlState[i] = lState;
            minchi2       = chi2;
            if ( !m_allCombinations.value() ) { best_matchedTrack = createMatchedTrack( *( *t ), *( *m ) ); }
          }
        } else {
          debug() << "matching failed " << chi2 << endmsg;
          sc = StatusCode::SUCCESS;
        }
      }
      if ( !m_allCombinations.value() && minchi2 < m_matchChi2Cut.value() ) {
        matchedTracks.insert( best_matchedTrack.release() );
        ++m_countMatchedTracks;
      }

      bestchi2.push_back( minchi2 );
    }
  }
  return matchedTracks;
}

#ifdef USE_DD4HEP

StatusCode TrackMuonMatching::matchChi2( LHCb::State&, LHCb::State&, const double&, double&,
                                         GenericDetElem const& ) const {
  throw GaudiException( "matchChi2 not implemented for DD4hep", "TrackMuonMatching", StatusCode::FAILURE );
}
StatusCode TrackMuonMatching::longTmuonExtrap( LHCb::State*, const double&, GenericDetElem const& ) const {
  throw GaudiException( "longTmuonExtrap not implemented for DD4hep", "TrackMuonMatching", StatusCode::FAILURE );
}

#else

StatusCode TrackMuonMatching::matchChi2( LHCb::State& lState, LHCb::State& mState, const double& atZ, double& chi2,
                                         GenericDetElem const& detelem ) const {
  /// Extrapolate states
  StatusCode sc = m_extrapolator->propagate( lState, atZ, *detelem.geometry() );
  if ( !sc.isSuccess() ) { return Warning( "Could not propagate longTrack state", StatusCode::FAILURE, 5 ); }

  sc = m_extrapolator->propagate( mState, atZ, *detelem.geometry() );
  if ( !sc.isSuccess() ) { return Warning( "Could not propagate Muon state", StatusCode::FAILURE, 5 ); }

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Extrapolated longTrack state to z = " << atZ << " is " << lState << endmsg
            << "Extrapolated Muon state to z = " << atZ << " is " << mState << endmsg;
  }

  /// Now calculate the match chi2
  sc = m_chi2Calculator->calculateChi2( lState.stateVector(), lState.covariance(), mState.stateVector(),
                                        mState.covariance(), chi2 );
  if ( !sc.isSuccess() ) Error( "Could not invert matrices", StatusCode::FAILURE ).ignore();

  return sc;
}

StatusCode TrackMuonMatching::longTmuonExtrap( LHCb::State* lState, const double& atZ,
                                               GenericDetElem const& detelem ) const {
  /// Extrapolate states
  StatusCode sc = m_extrapolator->propagate( ( *lState ), atZ, *detelem.geometry() );
  if ( !sc.isSuccess() ) Warning( "Could not propagate longTrack state", StatusCode::FAILURE, 5 ).ignore();
  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Extrapolated longTrack state to z = " << atZ << " is " << ( *lState ) << endmsg;
  }
  return sc;
}

#endif
