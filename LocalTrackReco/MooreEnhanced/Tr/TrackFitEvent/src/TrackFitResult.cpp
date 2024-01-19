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
#include "Event/TrackFitResult.h"
#include "Event/TrackFunctor.h"
#include <algorithm>

using namespace LHCb;

//=============================================================================
// Copy constructor. Hidden: use clone method instead
//=============================================================================
TrackFitResult::TrackFitResult( const TrackFitResult& rhs ) {
  m_nIter          = rhs.m_nIter;
  m_pScatter       = rhs.m_pScatter;
  m_chi2           = rhs.m_chi2;
  m_chi2Velo       = rhs.m_chi2Velo;
  m_chi2Upstream   = rhs.m_chi2Upstream;
  m_chi2Long       = rhs.m_chi2Long;
  m_chi2Muon       = rhs.m_chi2Muon;
  m_chi2Downstream = rhs.m_chi2Downstream;
  m_measurements   = rhs.m_measurements;

  // copy the nodes. be sure to remap the measurement.
  m_nodes.reserve( rhs.m_nodes.size() );
  for ( const auto& node : rhs.m_nodes ) {
    m_nodes.push_back( new LHCb::FitNode( *node ) );
    if ( node->hasMeasurement() ) {
      auto it = std::find_if( rhs.m_measurements.begin(), rhs.m_measurements.end(),
                              [&]( const LHCb::Measurement& m ) { return &m == &( node->measurement() ); } );
      if ( it == rhs.m_measurements.end() ) {
        throw GaudiException( "TrackFitResult::copy: found a node pointing to a measurement not on track!", __func__,
                              StatusCode::FAILURE );
      }
      m_nodes.back()->setMeasurement( m_measurements[it - rhs.m_measurements.begin()] );
    }
  }
}

//=============================================================================
// Destructor
//=============================================================================
TrackFitResult::~TrackFitResult() { std::for_each( m_nodes.begin(), m_nodes.end(), TrackFunctor::deleteObject() ); }

//=============================================================================
// Clone the track
//=============================================================================
std::unique_ptr<LHCb::ITrackFitResult> TrackFitResult::clone() const {
  return std::unique_ptr<TrackFitResult>( new TrackFitResult{*this} );
}

//=============================================================================
// Add a list of measurement to the list associated to the Track. This takes ownership.
//=============================================================================
void TrackFitResult::setMeasurements( std::vector<LHCb::Measurement>&& measurements ) {
  if ( !m_measurements.empty() &&
       std::any_of( m_nodes.begin(), m_nodes.end(), []( const auto* n ) { return n->hasMeasurement(); } ) ) {
    throw GaudiException( "attempt to remove measurements which are still referenced", __func__, StatusCode::FAILURE );
  }
  m_measurements = std::move( measurements );
}

//=============================================================================
// Return the Measurement on the Track corresponding to the input LHCbID
//=============================================================================
const Measurement* TrackFitResult::measurement( const LHCbID& value ) const {
  auto it = std::find_if( m_measurements.begin(), m_measurements.end(),
                          [&]( const Measurement& m ) { return m.lhcbID() == value; } );
  return it != m_measurements.end() ? &( *it ) : nullptr;
}

//=============================================================================
// Remove all measurements from the track
//=============================================================================
void TrackFitResult::clearMeasurements() {
  // remove all nodes first
  clearNodes();
  // now remove the measurements
  m_measurements.clear();
}

//=============================================================================
// Remove all nodes from the track
//=============================================================================
void TrackFitResult::clearNodes() {
  std::for_each( m_nodes.begin(), m_nodes.end(), TrackFunctor::deleteObject() );
  m_nodes.clear();
}

//=============================================================================
// reset the track
//=============================================================================
void TrackFitResult::reset() {
  m_nIter = 0;
  clearMeasurements();
}

//=============================================================================
// Count the number of outliers
//=============================================================================

unsigned int LHCb::TrackFitResult::nOutliers() const {
  return std::count_if( nodes().begin(), nodes().end(),
                        []( const FitNode* node ) { return node->type() == LHCb::FitNode::Type::Outlier; } );
}

//=============================================================================
// Count the number of outliers
//=============================================================================

unsigned int LHCb::TrackFitResult::nActiveMeasurements() const { return m_measurements.size() - nOutliers(); }
