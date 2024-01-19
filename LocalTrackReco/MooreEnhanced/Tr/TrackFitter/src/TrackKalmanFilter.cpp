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
// Include files
// -------------
// from TrackEvent
#include "Event/FitNode.h"
#include "Event/KalmanFitResult.h"
#include "GaudiKernel/reverse.h"

// local
#include "TrackKalmanFilter.h"
#include <boost/foreach.hpp>

//-----------------------------------------------------------------------------
// Implementation file for class : TrackKalmanFilter
//
//  Original Author: Rutger van der Eijk
//  Created: 07-04-1999
//  Adapted: 21-03-2002  Olivier Callot
//  Adapted: 15-04-2005  Jose A. Hernando, Eduardo Rodrigues
//  Adapted: 20-10-2006  J. van Tilburg
//-----------------------------------------------------------------------------

DECLARE_COMPONENT( TrackKalmanFilter )

//=========================================================================
// Fit the nodes
//=========================================================================
StatusCode TrackKalmanFilter::fit( LHCb::Track& track ) const {
  // if( msgLevel( MSG::DEBUG ) ) debug() << "TrackKalmanFilter::fit" << endmsg ;

  StatusCode sc{StatusCode::SUCCESS};

  // The seed covariance comes from the KalmanFitResult
  LHCb::KalmanFitResult* kalfit = dynamic_cast<LHCb::KalmanFitResult*>( track.fitResult() );
  if ( !kalfit ) return Warning( "No kalfit on track", StatusCode::FAILURE, 0 );

  const auto& nodes = kalfit->nodes();
  if ( nodes.empty() ) return Warning( "Fit failure: track has no nodes", StatusCode::FAILURE, 0 );

  // This is set up with the aim to trigger the cache such that there
  // will be no nested calls. That makes it easier to profile the
  // fit. Eventually, we could do without all of this, since
  // everything is anyway computed on demand.
  for ( auto& node : nodes ) {
    node->predictedState( LHCb::FitNode::Forward );
    node->filteredState( LHCb::FitNode::Forward );
  }
  if ( m_forceBiDirectionalFit ) {
    for ( auto& node : reverse( nodes ) ) {
      node->predictedState( LHCb::FitNode::Backward );
      node->filteredState( LHCb::FitNode::Backward );
    }
  }
  // force the smoothing.
  if ( m_forceSmooth ) {
    for ( auto& node : nodes ) node->state();
  }

  if ( kalfit->inError() ) { sc = Warning( kalfit->getError(), StatusCode::FAILURE, 0 ); }

  // This is the only thing we need KalmanFilter still for: set the
  // total chisquare of the track. This could also be done from
  // TrackMasterFitter.

  // Count the number of active track parameters. For now, just look at the momentum.
  size_t npar = m_DoF;
  if ( npar == 5u ) {
    auto i = std::find_if( nodes.rbegin(), nodes.rend(),
                           []( const auto* node ) { return node->type() == LHCb::FitNode::Type::HitOnTrack; } );
    if ( i != nodes.rend() ) {
      const LHCb::State& laststate = ( *i )->filteredState( LHCb::FitNode::Forward );
      const double       threshold = 0.1;
      npar = ( laststate.covariance()( 4, 4 ) / kalfit->seedCovariance()( 4, 4 ) < threshold ? 5 : 4 );
      dynamic_cast<LHCb::KalmanFitResult*>( track.fitResult() )->setNTrackParameters( npar );
    }
  }

  LHCb::ChiSquare chisq = ( *nodes.rbegin() )->totalChi2( LHCb::FitNode::Forward );
  int             ndof  = chisq.nDoF() - ( npar - kalfit->nTrackParameters() );
  if ( m_forceBiDirectionalFit ) {
    // FIXME: why don't we take the maximum, like in KalmanFitResult?
    LHCb::ChiSquare chisqbw = ( *nodes.begin() )->totalChi2( LHCb::FitNode::Backward );
    if ( chisqbw.chi2() < chisq.chi2() ) chisq = chisqbw;
  }
  track.setChi2AndDoF( chisq.chi2(), ndof );

  if ( msgLevel( MSG::DEBUG ) ) printStates( track );
  return sc;
}

void TrackKalmanFilter::printErrMeasures( LHCb::Track& track ) const {
  LHCb::KalmanFitResult* kalfit = dynamic_cast<LHCb::KalmanFitResult*>( track.fitResult() );
  if ( kalfit ) {
    for ( const auto* inode : kalfit->nodes() ) {
      std::string errMeasure_s = inode->visit_r<std::string>( [&]( const auto& n ) {
        const auto& em = n.errMeasure();
        return std::accumulate( em.begin(), em.end(), std::string{},
                                []( auto s, auto i ) { return s + std::to_string( i ) + "  "; } );
      } );
      debug() << "FitNode errMeasure " << errMeasure_s << endmsg;
    }
  }
}

void TrackKalmanFilter::printStates( LHCb::Track& track ) const {
  LHCb::KalmanFitResult* kalfit = dynamic_cast<LHCb::KalmanFitResult*>( track.fitResult() );
  if ( kalfit ) {
    for ( const auto* inode : kalfit->nodes() ) { debug() << *inode << endmsg; }
  }
}
