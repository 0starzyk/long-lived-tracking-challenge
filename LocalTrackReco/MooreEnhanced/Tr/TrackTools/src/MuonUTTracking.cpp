/*****************************************************************************\
* (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
// Include files
#include "Event/PrTracksTag.h"
#include "Event/Track.h"
#include "Event/Track_v1.h"
#include "LHCbAlgs/Transformer.h"

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/PrimaryVertices.h"
#include "Event/State.h"
#include "Event/StateParameters.h"

#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackInterfaces/ITrackFitter.h"
#include "boost/container/static_vector.hpp"

#include "PrKernel/IPrAddUTHitsTool.h"
#include "PrKernel/PrMutUTHits.h"

// Declaration of the Algorithm Factory

//-----------------------------------------------------------------------------
// Implementation file for class : MuonUTTrackinging
/*
 * \brief  Make a MuonUTTracking: Get muon standalone track, add UT hits, refit
 **
 ** Parameters:
 ** - Extrapolator: Name for the track extrapolator.
 ** - Fitter: Name for the TrackMasterFitter.
 ** - MinNUTHits: Minimal number of UT hits that need to be added to save the track.
 ** - InputMuonTracks: The location the input standalone muon tracks read from.
 ** - RecVertices: The location of the Primary vertices.
 ** - OutputTracks: The location the tracks should be written to.
 *
 */
//-----------------------------------------------------------------------------
namespace {
  using Vertices = LHCb::Event::PV::PrimaryVertexContainer;
}

class MuonUTTracking final : public LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Tracks&, const Vertices&,
                                                                               const IPrAddUTHitsTool& utHitAddingTool,
                                                                               DetectorElement const& ),
                                                                 LHCb::DetDesc::usesConditions<DetectorElement>> {

public:
  MuonUTTracking( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"InputMuonTracks", ""}, KeyValue{"RecVertices", ""},
                      KeyValue{"AddUTHitsTool", "PrAddUTHitsTool"},
                      KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}},
                     KeyValue{"OutputTracks", ""} ) {}

  LHCb::Tracks operator()( const LHCb::Tracks&, const Vertices&, const IPrAddUTHitsTool&,
                           DetectorElement const& ) const override;

private:
  Gaudi::XYZPoint fillPVs( const Vertices& ) const;
  LHCb::State     iterateToPV( LHCb::Track*, Gaudi::XYZPoint&, double, IGeometryInfo const& ) const;

  //---Properties
  Gaudi::Property<double>       m_tolerance{this, "Tolerance", 0.5};
  Gaudi::Property<unsigned int> m_maxIterations{this, "MaxIteration", 10};
  Gaudi::Property<unsigned int> m_minNumberUTHits{this, "MinNUTHits", 2};

  //--- Tools
private:
  ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackMasterExtrapolator"}; ///< extrapolator
  ToolHandle<ITrackFitter>       m_trackFitter{this, "Fitter", "TrackMasterFitter"};

  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING>     m_failed_propagation{this,
                                                                             "Could not propagate state to VELO!"};
  mutable Gaudi::Accumulators::Counter<>                    m_countEvents{this, "nEvents"};
  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_countMuonTracks{this, "nMuonTracks"};
  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_countMuonUTTracks{this, "nMuonTracks with UT hits added"};
  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_fitMuonUTTracks{this, "nMuonUTTracks passing fit"};
};

DECLARE_COMPONENT( MuonUTTracking )

//=============================================================================
// Main execution
//=============================================================================
LHCb::Tracks MuonUTTracking::operator()( const LHCb::Tracks& muontracks, const Vertices& pvs,
                                         const IPrAddUTHitsTool& utHitAddingTool, DetectorElement const& lhcb ) const {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  // -- Fill PVs
  auto PVPos = fillPVs( pvs );
  ++m_countEvents;

  LHCb::Tracks outputTracks;
  outputTracks.reserve( muontracks.size() );
  // -- Collect the UT hits
  boost::container::static_vector<LHCb::LHCbID, LHCb::Pr::TracksInfo::MaxUTHits> utHits;

  m_countMuonTracks += muontracks.size();
  // -- Loop over all Muon Tracks
  for ( auto& track : muontracks ) {
    // -- Change Q/p until it points to the PV (stolen from Wouter)
    auto muonState = iterateToPV( track, PVPos, track->firstState().qOverP(),
                                  *lhcb.geometry() ); // -- This is the function that iterates

    auto veloStateVec = muonState.stateVector();
    auto sc           = m_extrapolator->propagate( veloStateVec, muonState.z(), PVPos.z(), nullptr, *lhcb.geometry() );
    if ( !sc ) {
      ++m_failed_propagation;
      continue;
    }

    LHCb::State veloState( veloStateVec, PVPos.z(), LHCb::State::Vertex );
    /// -- Add UT hits
    utHits.clear();
    utHitAddingTool.getUTHits( veloState, utHits );
    if ( msgLevel( MSG::DEBUG ) ) debug() << "Found " << utHits.size() << " UT hits to add" << endmsg;

    // -- Skip if not enough UT hits were found
    if ( utHits.size() < m_minNumberUTHits ) {
      if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Not enough hits in UT found" << endmsg;
      continue;
    }
    // -- Add the UT hits
    for ( const auto hit : utHits ) {
      const auto lhcbid = LHCb::LHCbID( hit );
      track->addToLhcbIDs( lhcbid );
    }

    // -- Set Pattern Reco status and track type, finally fit the track
    track->setPatRecStatus( LHCb::Track::PatRecStatus::PatRecIDs );
    track->setType( LHCb::Track::Types::MuonUT );
    track->setChi2PerDoF( 1000 ); // -- Try to set a dummy value to satisfy the tupletool
    track->setNDoF( 1 );

    track->clearStates();
    track->addToStates( veloState );
    track->addToStates( muonState );
    m_countMuonUTTracks += 1;
    sc = ( *m_trackFitter )( *track, *lhcb.geometry(), LHCb::Tr::PID::Muon() );
    if ( !sc && msgLevel( MSG::DEBUG ) ) debug() << "Fit failed" << endmsg;
    if ( sc ) outputTracks.insert( track );
  }

  if ( msgLevel( MSG::DEBUG ) ) debug() << "Filling number of tracks " << outputTracks.size() << endmsg;

  m_fitMuonUTTracks += outputTracks.size();
  return outputTracks;
}

//=============================================================================
//  Fill PVs
//=============================================================================
Gaudi::XYZPoint MuonUTTracking::fillPVs( const Vertices& pvs ) const {

  // Return the position of the PV with the largest multiplicity
  auto it = std::max_element( pvs.begin(), pvs.end(),
                              []( const auto& lhs, const auto& rhs ) { return lhs.nDoF() < rhs.nDoF(); } );
  return it != pvs.end() ? it->position() : Gaudi::XYZPoint{0.0, 0.0, 0.0};
}

//=============================================================================
//  Change the q/p till the track points to the PV (stolen from Wouter)
//=============================================================================
LHCb::State MuonUTTracking::iterateToPV( LHCb::Track* track, Gaudi::XYZPoint& PVPos, double qOverP,
                                         IGeometryInfo const& geometry ) const {

  LHCb::State muonState = track->closestState( 15270. );
  muonState.setQOverP( qOverP );

  // -- Now call the extrapolator and iterate until we have the desired accuracy.
  LHCb::State dummyState = muonState;

  for ( unsigned int i = 0; i < m_maxIterations.value(); ++i ) {
    dummyState = muonState;
    m_extrapolator
        ->propagate( dummyState.stateVector(), dummyState.z(), PVPos.z(), nullptr, geometry, LHCb::Tr::PID::Muon() )
        .ignore();
    const auto deltaX = -( dummyState.x() - PVPos.x() );

    if ( std::abs( deltaX ) > m_tolerance.value() ) break;
  }

  return muonState;
}
