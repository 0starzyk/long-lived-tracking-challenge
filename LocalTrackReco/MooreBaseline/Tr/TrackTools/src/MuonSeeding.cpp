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
#include "Event/RecVertex.h"
#include "Event/RecVertex_v2.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "Event/VPCluster.h"
#include "GaudiKernel/Chrono.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Transformer.h"
#include "StandaloneMuonTrack.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackInterfaces/ITrackFitter.h"

//-----------------------------------------------------------------------------
// Implementation file for class : MuonSeeding
//
// 2010-09-14 : Michel De Cian
//-----------------------------------------------------------------------------

/** @class MuonSeeding MuonSeeding.h
 *
 * \brief  Make a MuonSeeding: Get muon standalone tracks
 *
 * Parameters:
 * - InputMuonTracks: The location the input tracks read from.
 * - OutputMuonTracks: The location the tracks should be written to.
 * Properties
 * - Extrapolator: Name for the track extrapolator.
 * - Fittter: Name of fitter used for the track fit.
 *
 *  @author Michel De Cian
 *  @date   2010-09-20
 */

class MuonSeeding : public LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Tracks&, DetectorElement const& ),
                                                        LHCb::DetDesc::usesConditions<DetectorElement>> {

public:
  using base_t = LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Tracks&, DetectorElement const& ),
                                              LHCb::DetDesc::usesConditions<DetectorElement>>;
  /// Standard constructor
  MuonSeeding( const std::string& name, ISvcLocator* pSvcLocator )
      : base_t( name, pSvcLocator,
                {KeyValue{"InputMuonTracks", "Rec/Track/Muon"},
                 KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}},
                KeyValue{"OutputMuonTracks", "Rec/Track/MuonSeed"} ) {}

  LHCb::Tracks operator()( const LHCb::Tracks&    tracks,
                           DetectorElement const& lhcb ) const override; ///< Algorithm execution

private:
  // -- Methods
  StatusCode iterateToPV( LHCb::Track* track, LHCb::State& muonState, LHCb::State& veloState,
                          const Gaudi::XYZPoint& PVPos, double qOverP, IGeometryInfo const& lhcb ) const;

  // -- Properties
  Gaudi::Property<bool> m_fitTracks{this, "FitTracks", true};

  // -- Tools
  ToolHandle<ITrackFitter>       m_trackFitter{this, "Fitter", "TrackMasterFitter/Fitter"};
  ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackMasterExtrapolator/Extrapolator"};

  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_countMuCandidates{this, "nMuonTrackCandidates"};
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING>     m_failed_iteration{this, "Failed iteration to PV!"};
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING>     m_failed_trackfit{this, "Failed track fit!"};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( MuonSeeding )

//=============================================================================
// Main execution
//=============================================================================
LHCb::Tracks MuonSeeding::operator()( const LHCb::Tracks& muonTracks, DetectorElement const& lhcb ) const {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  // -- This is where we assume the track came from, Do not use the real PV
  // for the standalone Muon, which won't help too much due to bad resolution
  // and the hits search is (0,0,0) oriented
  Gaudi::XYZPoint PVPos( 0., 0., 0. );

  LHCb::Tracks outputTracks{};
  outputTracks.reserve( muonTracks.size() );

  // -- Loop over all Muon Tracks
  for ( const auto& muontrack : muonTracks ) {

    if ( m_fitTracks ) {
      // Try and improve the covariance information
      auto sc = ( *m_trackFitter )( *muontrack, *lhcb.geometry(), LHCb::Tr::PID::Muon() );
      if ( sc.isFailure() ) { ++m_failed_trackfit; }
    }
    // -- Change q/p until it points to the origin (adapted from Wouter)
    LHCb::State veloState, muonState;
    auto sc = iterateToPV( muontrack, muonState, veloState, PVPos, muontrack->firstState().qOverP(), *lhcb.geometry() );
    if ( sc.isFailure() ) { ++m_failed_iteration; }

    // -- Set Pattern Reco status and track type, finally fit the track
    muontrack->setPatRecStatus( LHCb::Track::PatRecStatus::PatRecIDs );
    muontrack->setType( LHCb::Track::Types::Muon );

    if ( !m_fitTracks ) {
      muontrack->clearStates(); // remove the state recMomentum created
    }

    muontrack->addToStates( veloState );
    muontrack->addToStates( muonState );

    outputTracks.add( muontrack );
  }

  m_countMuCandidates += outputTracks.size();

  return outputTracks;
}

//=============================================================================
//  Change the q/p till the track points to the PV (stolen from Wouter)
//=============================================================================
StatusCode MuonSeeding::iterateToPV( LHCb::Track* track, LHCb::State& muonState, LHCb::State& veloState,
                                     const Gaudi::XYZPoint& PVPos, double qOverP,
                                     IGeometryInfo const& geometry ) const {
  muonState = track->closestState( 15000 );

  muonState.setQOverP( qOverP );

  // Set the y slope based on the target position at ~the origin
  muonState.setTy( ( muonState.y() - PVPos.y() ) / ( muonState.z() - PVPos.z() ) );

  // Set the uncertainty on ty to just come from the y uncertainty from the muon stations
  auto cov    = muonState.covariance();
  cov( 3, 3 ) = muonState.ty() * muonState.ty() * ( cov( 1, 1 ) / ( muonState.y() - PVPos.y() ) );
  muonState.setCovariance( cov );

  // -- Now call the extrapolator and iterate until we have the desired accuracy.
  const double tolerance     = 0.5; // [mm]
  const int    maxIterations = 10;

  auto veloStateVec = muonState.stateVector();
  for ( int i = 0; i < maxIterations; ++i ) {
    auto sc =
        m_extrapolator->propagate( veloStateVec, muonState.z(), PVPos.z(), nullptr, geometry, LHCb::Tr::PID::Muon() );
    if ( sc.isFailure() ) { return StatusCode::FAILURE; }
    const auto deltaX = -( muonState.x() - PVPos.x() );
    if ( std::abs( deltaX ) > tolerance ) break;
  }

  veloState.setState( veloStateVec );
  veloState.setLocation( LHCb::State::Location::ClosestToBeam );
  muonState.setLocation( LHCb::State::Location::Muon );

  return m_extrapolator->propagate( veloState.stateVector(), veloState.z(), PVPos.z(), nullptr, geometry,
                                    LHCb::Tr::PID::Muon() );
}
