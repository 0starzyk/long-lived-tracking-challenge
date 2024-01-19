/*****************************************************************************\
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

// local
#include "MeasurementProvider.h"

// from TrackFitEvent
#include "Event/FitNode.h"

using namespace LHCb;

//-----------------------------------------------------------------------------
// Implementation file for class : MeasurementProvider
//
// 2005-04-14 : Jose Angel Hernando Morata
//-----------------------------------------------------------------------------

DECLARE_COMPONENT( MeasurementProvider )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MeasurementProvider::MeasurementProvider( const std::string& type, const std::string& name, const IInterface* parent )
    : extends( type, name, parent ) {
  declareProperty( "VPProvider", m_vpProvider );
  declareProperty( "UTProvider", m_utProvider );
  declareProperty( "FTProvider", m_ftProvider );
  declareProperty( "MuonProvider", m_muonProvider );
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode MeasurementProvider::initialize() {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "MeasurementProvider::initialize()" << endmsg;
  StatusCode sc = extends::initialize();
  if ( sc.isFailure() ) return sc; // error already reported by base class

  if ( !m_ignoreVP.value() ) {
    sc = m_vpProvider.retrieve();
    if ( sc.isFailure() ) return sc;
    m_providers.emplace_back( LHCb::Measurement::Type::VP, &( *m_vpProvider ) );
  } else {
    m_vpProvider.disable();
  }

  if ( !m_ignoreUT.value() ) {
    sc = m_utProvider.retrieve();
    if ( sc.isFailure() ) return sc;
    m_providers.emplace_back( LHCb::Measurement::Type::UT, &( *m_utProvider ) );
  } else {
    m_utProvider.disable();
  }

  if ( !m_ignoreFT.value() ) {
    sc = m_ftProvider.retrieve();
    if ( sc.isFailure() ) return sc;
    m_providers.emplace_back( LHCb::Measurement::Type::FT, &( *m_ftProvider ) );
  } else {
    m_ftProvider.disable();
  }

  if ( !m_ignoreMuon.value() ) {
    sc = m_muonProvider.retrieve();
    if ( sc.isFailure() ) return sc;
    m_providers.emplace_back( LHCb::Measurement::Type::Muon, &( *m_muonProvider ) );
  } else {
    m_muonProvider.disable();
  }
  return sc;
}

StatusCode MeasurementProvider::finalize() {
  StatusCode sc;
  if ( msgLevel( MSG::DEBUG ) ) debug() << "In MeasurementProvider::finalize. Releasing tool handles." << endmsg;
  // make sure to release all toolhandles
  if ( !m_ignoreVP.value() ) {
    sc = m_vpProvider.release();
    if ( sc.isFailure() ) return sc;
  }
  if ( !m_ignoreUT.value() ) {
    sc = m_utProvider.release();
    if ( sc.isFailure() ) return sc;
  }
  if ( !m_ignoreFT.value() ) {
    sc = m_ftProvider.release();
    if ( sc.isFailure() ) return sc;
  }
  if ( !m_ignoreMuon.value() ) {
    sc = m_muonProvider.release();
    if ( sc.isFailure() ) return sc;
  }
  sc = GaudiTool::finalize();
  return sc;
}

//=============================================================================
// Load all the Measurements from the list of LHCbIDs on the input Track
//=============================================================================
StatusCode MeasurementProvider::load( Track& track ) const {
  // make sure we have a place to store the result
  if ( !track.fitResult() ) track.setFitResult( new LHCb::TrackFitResult() );

  auto* fit = fitResult( track );

  std::vector<LHCbID> newids;
  newids.reserve( track.lhcbIDs().size() );
  for ( const auto& id : track.lhcbIDs() ) {
    // First look if the Measurement corresponding to this LHCbID
    // is already in the Track, i.e. whether it has already been loaded!
    if ( fit->measurement( id ) ) {
      Warning( "Found measurements already loaded on track!", StatusCode::SUCCESS, 0 ).ignore();
      if ( msgLevel( MSG::DEBUG ) || msgLevel( MSG::VERBOSE ) )
        debug() << "Measurement had already been loaded for the LHCbID"
                << " channelID, detectorType = " << id.channelID() << " , " << id.detectorType()
                << "  -> Measurement loading skipped for this LHCbID!" << endmsg;
    } else
      newids.push_back( id );
  }

  // create all measurements for selected IDs
  LHCb::TrackFitResult::MeasurementContainer newmeasurements;
  auto n2d = std::count_if( newids.begin(), newids.end(), []( LHCb::LHCbID id ) { return id.isVP() || id.isMuon(); } );
  newmeasurements.reserve( newids.size() + n2d ); // note: count VP & Muon, as they have 2 measurments/id

  // provide  a reference trajectory to addToMeasurements...
  addToMeasurements( newids, newmeasurements, LHCb::TrackTraj{track} );

  // add the measurements to the track
  fit->setMeasurements( std::move( newmeasurements ) );

  // Update the status flag of the Track
  track.setPatRecStatus( Track::PatRecStatus::PatRecMeas );

  return StatusCode::SUCCESS;
}

namespace {
  LHCb::Measurement::Type measurementtype( LHCb::LHCbID id ) {
    switch ( id.detectorType() ) {
    case LHCb::LHCbID::channelIDtype::VP:
      return LHCb::Measurement::Type::VP;
    case LHCb::LHCbID::channelIDtype::UT:
      return LHCb::Measurement::Type::UT;
    case LHCb::LHCbID::channelIDtype::FT:
      return LHCb::Measurement::Type::FT;
    case LHCb::LHCbID::channelIDtype::Muon:
      return LHCb::Measurement::Type::Muon;
    default:
      return LHCb::Measurement::Type::Unknown;
    }
  }
} // namespace

//-----------------------------------------------------------------------------
/// Create a list of measurements from a list of LHCbIDs
//-----------------------------------------------------------------------------
void MeasurementProvider::addToMeasurements( LHCb::span<LHCb::LHCbID> ids, std::vector<LHCb::Measurement>& measurements,
                                             const LHCb::ZTrajectory<double>& tracktraj ) const {
  // dispatch the measurements according to their type.
  for ( auto&& [type, provider] : m_providers ) {
    assert( provider != nullptr );
    auto pivot = std::stable_partition( ids.begin(), ids.end(),
                                        [&, type = type]( const auto id ) { return measurementtype( id ) == type; } );
    auto n     = std::distance( ids.begin(), pivot );
    provider->addToMeasurements( ids.first( n ), measurements, tracktraj );
    ids = ids.subspan( n );
  }
  assert( ids.empty() );
}

//-----------------------------------------------------------------------------
/// Project the state vector in this fitnode and update projection matrix and reference residual
//-----------------------------------------------------------------------------
StatusCode MeasurementProvider::projectReference( LHCb::FitNode& node ) const {
  StatusCode sc = StatusCode::FAILURE;
  if ( node.hasMeasurement() ) {
    auto mtype = measurementtype( node.measurement().lhcbID() );
    try {
      for ( auto&& [type, provider] : m_providers )
        if ( type == mtype ) {
          sc = provider->projectReference( node );
          break;
        }
    } catch ( StatusCode scr ) { sc = scr; }
  }
  return sc;
}

//-----------------------------------------------------------------------------
/// reset internal state, if any
//-----------------------------------------------------------------------------
void MeasurementProvider::reset() {
  for ( auto&& [type, provider] : m_providers ) {
    assert( provider != nullptr );
    provider->reset();
  }
}
