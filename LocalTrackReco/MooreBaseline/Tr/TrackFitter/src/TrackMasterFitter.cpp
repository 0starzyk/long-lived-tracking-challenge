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
// from Gaudi
#include "GaudiKernel/System.h"
#include "GaudiKernel/SystemOfUnits.h"

// from LHCb
#include "Kernel/HitPattern.h"
#include "Kernel/TransformedRange.h"
#include "LHCbMath/LHCbMath.h"

// from TrackInterfaces
#include "TrackInterfaces/IMaterialLocator.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

// from TrackEvent
#include "Event/FitNode.h"
#include "Event/KalmanFitResult.h"
#include "Event/Measurement.h"
#include "Event/State.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Event/TrackFunctor.h"
#include "TrackKernel/TrackFunctors.h"
#include "TrackKernel/TrackTraj.h"

// local
#include "TrackMasterFitter.h"

#include <cmath>

// stl
#include <algorithm>

#include "boost/container/static_vector.hpp"

#include "Event/ParametrisedScatters.h"

using namespace Gaudi;
using namespace LHCb;

namespace {

  //! Add info from fitter as extrainfo to track
  void fillExtraInfo( Track& track ) {
    // TODO: migrate clients to use KalmanFitResult directly. Then
    // remove the extrainfo field.

    // Clean up the track info
    track.eraseInfo( Track::AdditionalInfo::FitVeloChi2 );
    track.eraseInfo( Track::AdditionalInfo::FitVeloNDoF );
    track.eraseInfo( Track::AdditionalInfo::FitTChi2 );
    track.eraseInfo( Track::AdditionalInfo::FitTNDoF );
    track.eraseInfo( Track::AdditionalInfo::FitMatchChi2 );
    track.eraseInfo( Track::AdditionalInfo::NUTOutliers );

    auto kalfit = static_cast<const LHCb::KalmanFitResult*>( track.fitResult() );

    if ( track.hasT() ) {
      track.addInfo( Track::AdditionalInfo::FitTChi2, kalfit->chi2Downstream().chi2() );
      track.addInfo( Track::AdditionalInfo::FitTNDoF, kalfit->chi2Downstream().nDoF() );
    }

    if ( track.hasVelo() ) {
      track.addInfo( Track::AdditionalInfo::FitVeloChi2, kalfit->chi2Velo().chi2() );
      track.addInfo( Track::AdditionalInfo::FitVeloNDoF, kalfit->chi2Velo().nDoF() );
      if ( track.hasT() ) track.addInfo( Track::AdditionalInfo::FitMatchChi2, kalfit->chi2Match().chi2() );
    }

    if ( track.hasUT() ) {
      int n_ut_outliers =
          kalfit->nMeasurements<LHCb::Measurement::UT>() - kalfit->nActiveMeasurements<LHCb::Measurement::UT>();
      track.addInfo( Track::AdditionalInfo::NUTOutliers, n_ut_outliers );
    }
  }

  //=========================================================================
  // Determine the z-position of the closest approach to the beam line
  //  by linear extrapolation.
  //=========================================================================
  double closestToBeamLine( const State& state ) {
    const TrackVector& vec = state.stateVector();
    auto               z   = state.z();
    // check on division by zero (track parallel to beam line!)
    if ( vec[2] != 0 || vec[3] != 0 ) {
      z -= ( vec[0] * vec[2] + vec[1] * vec[3] ) / ( vec[2] * vec[2] + vec[3] * vec[3] );
    }
    // don't go outside the sensible volume
    return std::min( std::max( z, -100 * Gaudi::Units::cm ), StateParameters::ZBegRich2 );
  }

  const auto FirstLess = []( const auto& lhs, const auto& rhs ) { return lhs.first < rhs.first; };

  const auto DecreasingChi2 = []( const auto& lhs, const auto& rhs ) { return lhs.second > rhs.second; };

  const auto hasMeasurement = []( const LHCb::FitNode* node ) { return node->hasMeasurement(); };

} // namespace
//-----------------------------------------------------------------------------
// Implementation file for class : TrackMasterFitter
//
//  Original Author: Rutger van der Eijk
//  Created: 07-04-1999
//  Adapted: 21-03-2002  Olivier Callot
//  Adapted: 15-04-2005  Jose A. Hernando, Eduardo Rodrigues
//-----------------------------------------------------------------------------

DECLARE_COMPONENT( TrackMasterFitter )

//=========================================================================
// Initialize
//=========================================================================
StatusCode TrackMasterFitter::initialize() {
  StatusCode sc = base_class::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;

  m_debugLevel = msgLevel( MSG::DEBUG ) || msgLevel( MSG::VERBOSE );

  if ( m_debugLevel )
    debug() << " " << endmsg << "============ TrackMasterFitter Settings ===========" << endmsg
            << ( ( m_upstream.value() ) ? " Upstream fit" : " Downstream fit" ) << endmsg
            << " Number of fit iterations: " << m_numFitIter << endmsg << " Max " << m_numOutlierIter
            << " outliers removed with outliers"
            << " at chi2 > " << m_chi2Outliers << endmsg << " State z positions at: " << endmsg
            << ( ( m_stateAtBeamLine.value() ) ? " beam line," : "" ) << " first/last measurement"
            << ( m_addDefaultRefNodes.value() ? ", default reference positions" : "" ) << endmsg
            << "==================================================" << endmsg;

  return StatusCode::SUCCESS;
}

//=========================================================================
// Helper to print a failure comment
//=========================================================================
StatusCode TrackMasterFitter::failure( const std::string& comment ) const {
  if ( m_debugLevel ) debug() << "TrackMasterFitter failure: " + comment << endmsg;
  return Warning( comment, StatusCode::FAILURE, 1 );
}
//=========================================================================
// Helper to print a failure comment (with info instead of warning)
//=========================================================================
StatusCode TrackMasterFitter::failureInfo( const std::string& comment ) const {
  if ( m_debugLevel ) debug() << "TrackMasterFitter failure: " + comment << endmsg;
  return StatusCode::FAILURE;
}

//=========================================================================
// Fit a set of tracks
//=========================================================================
StatusCode TrackMasterFitter::operator()( LHCb::span<LHCb::Track> tracks, IGeometryInfo const& geometry,
                                          const LHCb::Tr::PID& pid ) const {
  auto       cache = createCache();
  StatusCode sc( StatusCode::SUCCESS );
  for ( auto& track : tracks ) {
    StatusCode sc_temp = fit_r( track, cache, geometry, pid );
    if ( sc_temp.isFailure() ) { sc = StatusCode::FAILURE; }
  }
  return sc;
}

//=========================================================================
// Fit the track
//=========================================================================
StatusCode TrackMasterFitter::operator()( LHCb::Track& track, IGeometryInfo const& geometry,
                                          const LHCb::Tr::PID& pid ) const {
  auto cache = createCache();
  return fit_r( track, cache, geometry, pid );
}

//=========================================================================
// Fit the track
//=========================================================================
StatusCode TrackMasterFitter::fit_r( Track& track, std::any& accelCache, IGeometryInfo const& geometry,
                                     const LHCb::Tr::PID pid ) const {
  if ( m_debugLevel ) debug() << " TrackMasterFitter::fit" << endmsg;

  // any track that doesnt make it to the end is failed
  track.setFitStatus( Track::FitStatus::FitFailed );

  // create the KalmanFitResult if it doesn't exist yet
  LHCb::KalmanFitResult* kalfitresult = dynamic_cast<LHCb::KalmanFitResult*>( track.fitResult() );
  if ( !kalfitresult ) {
    auto* fr     = fitResult( track );
    kalfitresult = ( fr ? new LHCb::KalmanFitResult( *fr ) : new LHCb::KalmanFitResult() );
    if ( m_useClassicalSmoother.value() ) kalfitresult->setBiDirectionnalSmoother( false );
    track.setFitResult( kalfitresult );
  }

  // Make the nodes from the measurements
  StatusCode  sc;
  const auto& nodes = kalfitresult->nodes();
  if ( nodes.empty() || m_makeNodes.value() ) {
    sc = makeNodes( track, pid, accelCache, geometry );
    if ( sc.isFailure() ) {
      kalfitresult->clearNodes();
      return failureInfo( "unable to make nodes from the measurements" );
    }
  } else {
    sc = updateRefVectors( track, pid, m_updateTransport.value(), accelCache, geometry );
    if ( sc.isFailure() ) return failureInfo( "unable to update the ref vectors" );
  }

  // create a covariance matrix to seed the Kalman fit
  TrackSymMatrix seedCov; // Set off-diagonal elements to zero
  if ( m_useSeedStateErrors.value() ) {
    State state0 = track.firstState();
    auto  z1     = nodes.front()->z();
    extrapolator( track.type() )
        ->propagate( state0, z1, geometry )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    seedCov = state0.covariance();
    if ( m_debugLevel )
      debug() << " state0 at z " << z1 << " vector " << state0.stateVector() << "\n"
              << " covariance " << state0.covariance() << endmsg;
  } else {
    seedCov( 0, 0 ) = m_errorX * m_errorX;
    seedCov( 1, 1 ) = m_errorY * m_errorY;
    seedCov( 2, 2 ) = m_errorTx * m_errorTx;
    seedCov( 3, 3 ) = m_errorTy * m_errorTy;
    seedCov( 4, 4 ) = std::pow( m_errorQoP[0] * track.firstState().qOverP(), 2 ) + std::pow( m_errorQoP[1], 2 );
  }
  kalfitresult->setSeedCovariance( seedCov );
  nodes.front()->setSeedCovariance( seedCov );
  nodes.back()->setSeedCovariance( seedCov );

  if ( m_debugLevel )
    debug() << "SeedState: z = " << nodes.front()->z() << " stateVector = " << nodes.front()->refVector()
            << " covariance  = " << seedCov << endmsg;

  // Iterate the track fit for linearisation. Be careful with chi2
  // convergence here: The first iteration might not be using OT
  // drifttimes in which case the chi2 can actually go up in the 2nd
  // iteration.
  int  iter      = 1;
  bool converged = false;
  for ( ; iter <= m_numFitIter && !converged; ++iter ) {
    if ( m_debugLevel ) debug() << "Iteration # " << iter << endmsg;

    // update reference trajectories with smoothed states
    // TODO: combine this with the projection of the residuals which now resides in TrackKalmanFilter
    if ( iter > 1 ) {
      sc = updateRefVectors( track, pid, m_updateTransport.value() && iter <= m_maxUpdateTransports + 1, accelCache,
                             geometry );
      if ( sc.isFailure() ) return failureInfo( "problem updating ref vectors" );
    }

    auto prevchi2 = track.chi2();
    sc            = m_trackNodeFitter->fit( track );

    if ( sc.isFailure() ) return failureInfo( std::string( "unable to fit the track. " ) + kalfitresult->getError() );

    if ( m_debugLevel )
      debug() << "chi2 =  " << track.chi2() << " ref state = (" << nodes.back()->refVector()
              << ") at z= " << nodes.back()->z() << endmsg;
    auto dchi2 = prevchi2 - track.chi2();
    // require at least 3 iterations, because of the OT prefit.
    converged = iter > 1 && std::abs( dchi2 ) < m_maxDeltaChi2Converged * track.nDoF();
  }
  kalfitresult->setNIter( iter - 1 );

  // Outlier removal iterations
  iter                   = kalfitresult->nOutliers();
  LHCb::FitNode* outlier = nullptr;
  while ( iter < m_numOutlierIter && track.nDoF() > 1 && ( outlier = outlierRemoved( track ) ) ) {
    if ( m_debugLevel ) debug() << "Outlier iteration # " << iter << endmsg;

    // update reference trajectories with smoothed states
    if ( m_updateReferenceInOutlierIters.value() ) {
      sc = updateRefVectors( track, pid, m_updateTransport.value(), accelCache, geometry );
      if ( sc.isFailure() ) return failureInfo( "problem updating ref vectors" );
    }

    // catch the outlier partner hit if the outlier is VP
    LHCb::FitNode* outlier_partner = nullptr;
    if ( outlier->measurement().lhcbID().isVP() && outlier->hasMeasurement() ) {
      auto prev_node = outlier->prevNode( FitNode::Direction::Backward );
      auto next_node = outlier->prevNode( FitNode::Direction::Forward );
      auto outlierID = outlier->measurement().lhcbID().channelID();
      // prevNode returns const, a workaround without unconst:
      if ( prev_node && prev_node->hasMeasurement() && prev_node->measurement().lhcbID().channelID() == outlierID ) {
        outlier_partner = *( std::find( nodes.begin(), nodes.end(), prev_node ) );
      }
      if ( next_node && next_node->hasMeasurement() && next_node->measurement().lhcbID().channelID() == outlierID ) {
        outlier_partner = *( std::find( nodes.begin(), nodes.end(), next_node ) );
      }
    }

    // deactivate the outlier
    outlier->deactivateMeasurement();
    if ( outlier_partner ) { outlier_partner->deactivateMeasurement(); }

    // Call the track fit
    sc = m_trackNodeFitter->fit( track );

    if ( sc.isFailure() ) {
      // Should only be used if the track belongs to a container and therefore has a key!
      // std::ostringstream mess;
      // mess << "unable to fit the track # " << track.key();
      // return failure(mess.str());
      return failureInfo( "unable to fit the track" );
    }

    if ( m_debugLevel )
      debug() << "chi2 =  " << track.chi2() << " ref state = (" << nodes.back()->refVector()
              << ") at z= " << nodes.back()->z() << endmsg;
    ++iter;
  }

  // TODO: if outliers are removed, we add one more iteration to update the reference

  // determine the track states at user defined z positions
  sc = determineStates( track );
  if ( sc.isFailure() ) return failureInfo( "failed in determining states" );

  if ( m_debugLevel && !track.states().empty() ) debug() << "first state = " << track.firstState() << endmsg;

  // fill extra info
  if ( m_fillExtraInfo.value() ) fillExtraInfo( track );

  // declare the track successful if there are no errors
  if ( kalfitresult->inError() ) {
    sc = Warning( kalfitresult->getError(), StatusCode::FAILURE, 0 );
  } else {
    track.setFitStatus( Track::FitStatus::Fitted );
  }

  return sc;
}

//=============================================================================
//
//=============================================================================
StatusCode TrackMasterFitter::determineStates( Track& track ) const {
  StatusCode sc = StatusCode::SUCCESS;

  // clean the non-fitted states in the track!
  track.clearStates();

  const auto& nodes_ = nodes( track );

  // Add the state at the first and last measurement position
  // -----------------------------------------------
  auto  inode                = std::find_if( nodes_.cbegin(), nodes_.cend(), hasMeasurement );
  auto* firstMeasurementNode = ( inode != nodes_.cend() ? *inode : nullptr );
  auto  jnode                = std::find_if( nodes_.crbegin(), nodes_.crend(), hasMeasurement );
  auto* lastMeasurementNode  = ( jnode != nodes_.crend() ? *jnode : nullptr );

  bool upstream = nodes_.front()->z() > nodes_.back()->z();
  bool reversed = ( upstream && !track.isVeloBackward() ) || ( !upstream && track.isVeloBackward() );

  // This state is not filtered for a forward only fit.
  if ( m_addDefaultRefNodes.value() ) {
    State firststate = firstMeasurementNode->state();
    firststate.setLocation( reversed ? State::Location::LastMeasurement : State::Location::FirstMeasurement );
    track.addToStates( firststate );
  }
  // This state is always filtered
  State laststate = lastMeasurementNode->state();
  laststate.setLocation( reversed ? State::Location::FirstMeasurement : State::Location::LastMeasurement );
  track.addToStates( laststate );

  // Add the states at the reference positions
  // ------------------------------------------
  for ( const auto& node : nodes_ )
    if ( node->type() == LHCb::FitNode::Type::Reference ) track.addToStates( node->state() );

  if ( m_debugLevel ) {
    debug() << "Track " << track.key() << " has " << track.nStates() << " states after fit\n  at z = ";
    for ( const auto& s : track.states() ) debug() << s->z() << ", ";
    debug() << nMeasurements( track ) << " measurements used for the fit (out of " << track.nLHCbIDs() << ")."
            << endmsg;
  }
  return sc;
}

//=========================================================================
//
//=========================================================================
namespace {
  enum HitType { VP, UT, T, Muon };
  int hittypemap( const LHCb::Measurement& m ) {
    return m.visit( []( const auto& arg ) {
      using arg_t = std::decay_t<decltype( arg )>;
      if constexpr ( std::is_base_of_v<LHCb::Measurement::VP, arg_t> ) {
        return HitType::VP;
      } else if constexpr ( std::is_base_of_v<LHCb::Measurement::UT, arg_t> ) {
        return HitType::UT;
      } else if constexpr ( std::is_base_of_v<LHCb::Measurement::FT, arg_t> ) {
        return HitType::T;
      } else if constexpr ( std::is_base_of_v<LHCb::Measurement::Muon, arg_t> ) {
        return HitType::Muon;
      } else if constexpr ( std::is_base_of_v<LHCb::Measurement::VP2D, arg_t> ) {
        return HitType::VP;
      }
    } );
  }
} // namespace

LHCb::FitNode* TrackMasterFitter::outlierRemoved( Track& track ) const {
  // return true if outlier chi2 cut < 0
  LHCb::FitNode* outlier = nullptr;
  if ( m_chi2Outliers < 0.0 ) return outlier;

  // Count the number of hit layers of each type
  auto        fr    = dynamic_cast<const LHCb::KalmanFitResult*>( track.fitResult() );
  const auto& nodes = fr->nodes();
  // Note: this only works if VP==0,UT==1,T==2,Muon==3!
  const size_t     minNumPlanes[4] = {m_minNumVPHits, m_minNumUTHits, m_minNumTHits, m_minNumMuonHits};
  LHCb::HitPattern pattern{
      LHCb::TransformedRange{nodes, []( const auto& node ) { return node->measurement().lhcbID(); },
                             []( const auto& node ) { return node->type() == LHCb::FitNode::Type::HitOnTrack; }}};
  const size_t numPlanes[4] = {pattern.numVelo(), pattern.numUT(), pattern.numFT(), pattern.numMuon()};

  // loop over the nodes and find the one with the highest chi2 >
  // m_chi2Outliers, provided there is enough hits of this type left.

  // Computing the chi2 will trigger the smoothing. Especially in the
  // trigger, where we don't update the reference, this makes a big
  // difference. One way to save some time is to first test a number
  // of which we know that it is either equal to or bigger than the
  // chi2 contribution of the hit, namely the chi2 sum of the match
  // and the hit at this node. We then sort the hits in this number,
  // and only compute chi2 if it can be bigger than the current worst
  // one.
  LHCb::ChiSquare totalchi2   = ( *nodes.begin() )->totalChi2( LHCb::FitNode::Backward );
  LHCb::ChiSquare totalchi2_2 = ( *nodes.rbegin() )->totalChi2( LHCb::FitNode::Forward );
  if ( totalchi2_2.chi2() > totalchi2.chi2() ) totalchi2 = totalchi2_2;

  if ( m_debugLevel ) { debug() << "removeWorstOutlier: total chi2 " << totalchi2.chi2() << endmsg; }

  using NodeWithChi2 = std::pair<const FitNode*, double>;
  std::vector<NodeWithChi2> nodesWithChi2UL;
  nodesWithChi2UL.reserve( nodes.size() );
  int numtried( 0 ), numcalled( 0 );
  for ( const auto& node : nodes ) {
    if ( node->hasMeasurement() && node->type() == LHCb::FitNode::Type::HitOnTrack ) {
      int hittype = hittypemap( node->measurement() );
      if ( numPlanes[hittype] > minNumPlanes[hittype] ) {
        ++numtried;
        auto chi2MatchAndHit = totalchi2;
        for ( FitNode::Direction dir : {FitNode::Direction::Forward, FitNode::Direction::Backward} ) {
          if ( const auto* tmpnode = node->prevNode( dir ); tmpnode ) chi2MatchAndHit -= tmpnode->totalChi2( dir );
        }
        // correct DoF for the track parameters still included in the top approach
        // chi2MatchAndHit -= LHCb::ChiSquare(0,5); // OK but ignores edge hit cases
        int ndof        = node->visit_r<int>( []( auto& f ) { return f.typedim; } );
        chi2MatchAndHit = LHCb::ChiSquare( chi2MatchAndHit.chi2(), ndof );

        if ( m_debugLevel ) {
          const auto* prevnode = node->prevNode( FitNode::Direction::Forward );
          const auto* nextnode = node->prevNode( FitNode::Direction::Backward );

          debug() << "node LHCbID " << node->measurement().lhcbID().channelID() << " chi2Contribution "
                  << chi2MatchAndHit.chi2() << " [ndf=" << chi2MatchAndHit.nDoF() << "]"
                  << " (" << totalchi2.chi2() << " [ndf=" << totalchi2.nDoF() << "] - "
                  << ( prevnode ? prevnode->totalChi2( FitNode::Direction::Forward ).chi2() : 0 )
                  << " [ndf=" << ( prevnode ? prevnode->totalChi2( FitNode::Direction::Forward ).nDoF() : 0 ) << "] - "
                  << ( nextnode ? nextnode->totalChi2( FitNode::Direction::Backward ).chi2() : 0 )
                  << " [ndf=" << ( nextnode ? nextnode->totalChi2( FitNode::Direction::Backward ).nDoF() : 0 ) << "]"
                  << ")" << endmsg;
        }
        if ( chi2MatchAndHit.chi2PerDoF() > m_chi2Outliers ) {
          nodesWithChi2UL.emplace_back( node, chi2MatchAndHit.chi2PerDoF() );
        }
      }
    }
  }

  // now sort them
  auto worstChi2 = m_chi2Outliers;
  std::sort( nodesWithChi2UL.begin(), nodesWithChi2UL.end(), DecreasingChi2 );

  if ( m_debugLevel ) {
    debug() << "All bad measurements, ordered: " << endmsg;
    for ( const auto& nodeUL : nodesWithChi2UL ) {
      auto& node = nodeUL.first;
      debug() << "Measurement of type " << node->measurement().type() << " LHCbID "
              << node->measurement().lhcbID().channelID() << " at z " << node->z() << " with chi2Contribution/dof "
              << nodeUL.second << " and chi2/dof " << node->chi2().chi2PerDoF() << endmsg;
    }
  }

  // -- if the first is smaller than worstChi2, and it is sorted in decreasing order
  // -- the 'if' will never be true and we can return here (M. De Cian)
  if ( !nodesWithChi2UL.empty() && nodesWithChi2UL.front().second < worstChi2 ) return outlier;

  for ( const auto& node : nodesWithChi2UL ) {
    if ( node.second > worstChi2 ) {
      const auto chi2 = node.first->chi2().chi2PerDoF();
      ++numcalled;
      if ( chi2 > worstChi2 ) {
        worstChi2 = chi2;
        outlier   = const_cast<LHCb::FitNode*>( node.first );
      }
    }
  }
  if ( m_debugLevel ) {
    if ( outlier ) {
      debug() << "Measurement of type " << outlier->measurement().type() << " LHCbID "
              << outlier->measurement().lhcbID().channelID() << " at z " << outlier->z() << " with chi2/dof "
              << outlier->chi2().chi2PerDoF() << " removed." << endmsg;
    } else {
      debug() << "Outlier not found" << endmsg;
    }
  }

  return outlier;
}

//=========================================================================
// Update the measurements before a refit
//=========================================================================
StatusCode TrackMasterFitter::updateRefVectors( Track& track, const LHCb::Tr::PID pid, bool doUpdateTransport,
                                                std::any& accelCache, IGeometryInfo const& geometry ) const {
  if ( m_debugLevel ) debug() << "TrackMasterFitter::updateRefVectors" << endmsg;

  StatusCode sc = StatusCode::SUCCESS;
  for ( auto& node : nodes( track ) ) { node->setRefVector( node->state().stateVector() ); }
  if ( m_debugLevel ) debug() << "Ref vector for first state: " << nodes( track ).front()->refVector() << endmsg;

  // update the projections. need to be done every time ref is
  // updated. we can move this code here at some point.
  sc = projectReference( track );
  if ( sc.isFailure() ) return failureInfo( "problem projecting reference" );

  // update the material using the new ref vectors
  if ( m_applyMaterialCorrections.value() && m_updateMaterial.value() ) {
    sc = updateMaterialCorrections( track, pid, accelCache, geometry );
    if ( sc.isFailure() ) return failureInfo( "problem updating material" );
  }

  // update the transport using the new ref vectors
  if ( doUpdateTransport ) {
    sc = updateTransport( track, geometry );
    if ( sc.isFailure() ) return failureInfo( "problem updating transport" );
  }
  return sc;
}

//=========================================================================
// Update the measurements before a refit
//=========================================================================
StatusCode TrackMasterFitter::projectReference( LHCb::Track& track ) const {
  StatusCode sc = StatusCode::SUCCESS;
  for ( LHCb::FitNode* node : nodes( track ) ) {
    if ( !node->refIsSet() ) {
      sc = Warning( "Node without reference", StatusCode::FAILURE, 0 );
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Node without reference" << endmsg;
      break;
    }
    if ( !node->hasMeasurement() ) continue;
    // if the reference is not set, issue an error

    sc = m_measProvider->projectReference( *node );

    if ( sc.isFailure() ) {
      Warning( "unable to project statevector", sc, 0 ).ignore();
      if ( msgLevel( MSG::DEBUG ) ) debug() << "unable to project this statevector: " << node->refVector() << endmsg;
      break;
    }
  }
  return sc;
}

//=============================================================================
// Create the nodes from the measurements
//=============================================================================
StatusCode TrackMasterFitter::makeNodes( Track& track, const LHCb::Tr::PID pid, std::any& accelCache,
                                         IGeometryInfo const& geometry ) const {
  // Clear the nodes
  LHCb::KalmanFitResult& fitresult = static_cast<LHCb::KalmanFitResult&>( *( track.fitResult() ) );
  fitresult.clearNodes();

  // Clear the measurements if asked for
  if ( m_makeMeasurements.value() ) fitresult.clearMeasurements();

  // make sure the track has sufficient reference states
  if ( m_debugLevel ) debug() << "Track before making nodes: " << track << endmsg;
  StatusCode sc = initializeRefStates( track, geometry, pid );
  if ( sc.isFailure() ) return Warning( "Problems setting reference info", StatusCode::FAILURE, 1 );

  // Check if it is needed to populate the track with measurements
  if ( track.checkPatRecStatus( Track::PatRecStatus::PatRecIDs ) || fitresult.measurements().empty() ) {
    StatusCode sc1 = m_measProvider->load( track );
    if ( sc1.isFailure() ) return Error( "Unable to load measurements!", StatusCode::FAILURE );
    track.setPatRecStatus( Track::PatRecStatus::PatRecMeas );
    if ( m_debugLevel )
      debug() << "# LHCbIDs, Measurements = " << track.nLHCbIDs() << ", " << fitresult.nMeasurements() << endmsg;
  }

  // check that there are sufficient measurements. in fact, one is
  // enough for the fit not to fail
  if ( fitresult.measurements().empty() ) return Warning( "No measurements on track", StatusCode::FAILURE, 0 );

  // Create the nodes for the measurements.
  const auto&                    measures = fitresult.measurements();
  TrackFitResult::NodeContainer& nodes    = fitresult.nodes();
  nodes.reserve( measures.size() + 6 );
  std::transform( measures.rbegin(), measures.rend(), std::back_inserter( nodes ), []( const Measurement& m ) {
    // return new FitNode( m );
    return m.visit(
        [&]( LHCb::Measurement::VP2D const& ) -> FitNode* { return new FitNode( m, LHCb::FitNode::FitNode2DType{} ); },
        [&]( auto const& ) -> FitNode* { return new FitNode( m, LHCb::FitNode::FitNode1DType{} ); } );
  } );

  // Add reference nodes depending on track type
  // Note: all Muon track types are required to have EndVelo and closestToBeamLine states
  if ( m_addDefaultRefNodes.value() ) {
    if ( ( track.hasVelo() || track.hasMuon() ) && !track.isVeloBackward() )
      nodes.push_back( new FitNode( StateParameters::ZEndVelo, State::Location::EndVelo ) );
    if ( track.hasUT() || track.hasMuon() ) {
      nodes.push_back( new FitNode( StateParameters::ZBegRich1, State::Location::BegRich1 ) );
      nodes.push_back( new FitNode( StateParameters::ZEndRich1, State::Location::EndRich1 ) );
    }
    if ( track.hasT() || track.hasMuon() ) {
      nodes.push_back( new FitNode( StateParameters::ZBegT, State::Location::AtT ) );
      nodes.push_back( new FitNode( StateParameters::ZBegRich2, State::Location::BegRich2 ) );
      nodes.push_back( new FitNode( StateParameters::ZEndRich2, State::Location::EndRich2 ) );
    }
  }

  // At a node for the position at the beamline
  if ( m_stateAtBeamLine.value() && ( track.hasVelo() || track.hasMuon() ) ) {
    const LHCb::State& refstate = *( track.isVeloBackward() ? track.states().back() : track.states().front() );
    nodes.push_back( new FitNode( closestToBeamLine( refstate ), State::Location::ClosestToBeam ) );
  }

  // Sort the nodes in z
  bool backward = track.isVeloBackward();
  bool upstream = ( m_upstream.value() && !backward ) || ( !m_upstream.value() && backward );
  if ( upstream ) {
    std::stable_sort( nodes.begin(), nodes.end(), TrackFunctor::decreasingByZ() );
  } else {
    std::stable_sort( nodes.begin(), nodes.end(), TrackFunctor::increasingByZ() );
  }

  // Set the reference using a TrackTraj

  LHCb::TrackTraj tracktraj( track.states(), LHCb::Tag::State::AssumeSorted );
  std::for_each( nodes.begin(), nodes.end(),
                 [&]( LHCb::FitNode* node ) { node->setRefVector( tracktraj.stateVector( node->z() ) ); } );

  // set links between nodes
  fitresult.establishNodeLinks();

  // update the projections. need to be done every time ref is updated
  sc = projectReference( track );

  // add all the noise, if required
  if ( m_applyMaterialCorrections.value() && sc.isSuccess() ) {
    sc = updateMaterialCorrections( track, pid, accelCache, geometry );
  }

  // create the transport of the reference (after the noise, because it uses energy loss)
  return sc.isSuccess() ? updateTransport( track, geometry ) : sc;
}

//=========================================================================

StatusCode TrackMasterFitter::updateMaterialCorrections( LHCb::Track& track, const LHCb::Tr::PID pid,
                                                         std::any& accelCache, IGeometryInfo const& geometry ) const {
  if ( m_debugLevel ) debug() << "TrackMasterFitter::updateMaterialCorrections" << endmsg;

  // the noise in each node is the noise in the propagation between
  // the previous node and this node.

  // first collect all volumes on the track. The advantages of collecting them all at once
  // is that it is much faster. (Less call to TransportSvc.)
  LHCb::KalmanFitResult* fit   = static_cast<LHCb::KalmanFitResult*>( track.fitResult() );
  const auto&            nodes = fit->nodes();

  if ( nodes.size() > 1 ) {
    // only apply energyloss correction for tracks that traverse magnet
    bool applyenergyloss = m_applyEnergyLossCorrections.value() &&
                           ( track.hasMuon() || ( track.hasT() && ( track.hasVelo() || track.hasUT() ) ) );

    // if only velo, or magnet off, use a fixed momentum based on pt.
    auto scatteringMomentum = track.firstState().p();
    if ( m_scatteringPt > 0 && !track.hasT() && !track.hasUT() && !track.hasMuon() ) {
      auto tx            = track.firstState().tx();
      auto ty            = track.firstState().ty();
      auto slope2        = tx * tx + ty * ty;
      auto tanth         = std::max( std::sqrt( slope2 / ( 1 + slope2 ) ), 1e-4 );
      scatteringMomentum = m_scatteringPt / tanth;
    }
    // set some limits for the momentum used for scattering
    scatteringMomentum = std::min( scatteringMomentum, m_maxMomentumForScattering.value() );
    scatteringMomentum = std::max( scatteringMomentum, m_minMomentumForScattering.value() );

    // if m_scatteringP is set, use it
    if ( m_scatteringP > 0 ) scatteringMomentum = m_scatteringP;
    fit->setPScatter( scatteringMomentum );

    if ( m_debugLevel ) debug() << "scattering momentum: " << scatteringMomentum << endmsg;

    // this is farily tricky now: we want to use TracjTraj, but we
    // cannot create it directly from the nodes, because that would
    // trigger the filter!

    if ( m_useFastMaterialApproximation ) {
      // still missing the energyloss
      for ( size_t i{1}; i < nodes.size(); ++i ) {
        const auto [Q, deltaE] =
            TrackFit::param_scatter_impl::computeNoiseAndDeltaE( *nodes[i - 1], *nodes[i], scatteringMomentum );
        nodes[i]->setNoiseMatrix( Q );
        nodes[i]->setDeltaEnergy( deltaE );
      }

    } else {

      // note that this is the same traj as used in setting the first
      // ref. cannot we just save it somewhere?
      LHCb::TrackTraj tracktraj( track.states(), LHCb::Tag::State::AssumeSorted );
      auto            zmin = nodes.front()->z();
      auto            zmax = nodes.back()->z();
      if ( zmin > zmax ) std::swap( zmin, zmax );
      tracktraj.setRange( zmin, zmax );

      fit->intersections() = m_materialLocator->intersect( tracktraj, accelCache, geometry );
      // shall can we free up the space?
      fit->intersections().shrink_to_fit();

      // now we need to redistribute the result between the nodes. the first node cannot have any noise.
      auto inode   = nodes.begin();
      auto zorigin = ( *inode )->z();
      for ( ++inode; inode != nodes.end(); ++inode ) {
        FitNode*    node    = dynamic_cast<FitNode*>( *inode );
        auto        ztarget = node->z();
        LHCb::State state( node->refVector() );
        state.covariance() = Gaudi::TrackSymMatrix();
        state.setQOverP( 1 / scatteringMomentum );
        m_materialLocator->applyMaterialCorrections( state, fit->intersections(), zorigin, pid, true, applyenergyloss );
        auto deltaE = 1 / state.qOverP() - scatteringMomentum;

        node->setNoiseMatrix( state.covariance() );
        node->setDeltaEnergy( deltaE );
        /*
        auto Qalt = computeNoiseMatrix(*node,scatteringMomentum)  ;
        if( state.covariance()(2,2)>0 &&
            std::abs( Qalt(2,2) / state.covariance()(2,2) -1 ) > 1 ) {
          std::cout << "odd alternative noise matrix: " << state.z() << std::endl
                    << computeNoiseMatrix(*node,scatteringMomentum) << std::endl
                    << state.covariance() << std::endl ;
        }
        */

        zorigin = ztarget;
      }
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode TrackMasterFitter::updateTransport( LHCb::Track& track, IGeometryInfo const& geometry ) const {
  if ( m_debugLevel ) debug() << "TrackMasterFitter::updateTransport" << endmsg;

  // FIXME: This is probably the right place to set the number of
  // track parameters in the fitresult: simply check that the
  // propagation does not depend on momentum
  bool hasMomentum = false;

  // sets the propagation between the previous node and this. node that the reference
  // of the previous node is used.
  StatusCode sc = StatusCode::SUCCESS;

  LHCb::KalmanFitResult* fitresult = static_cast<LHCb::KalmanFitResult*>( track.fitResult() );
  const auto&            nodes     = fitresult->nodes();
  if ( nodes.size() > 1 ) {
    const ITrackExtrapolator* extrap    = extrapolator( track.type() );
    auto                      inode     = nodes.begin();
    const LHCb::StateVector*  refvector = &( ( *inode )->refVector() );
    TrackMatrix               F         = ROOT::Math::SMatrixIdentity();
    // Reset the first node
    ( *inode )->resetFilterStatus();
    for ( ++inode; inode != nodes.end() && sc.isSuccess(); ++inode ) {

      FitNode*          node        = *inode;
      auto              z           = node->z();
      LHCb::StateVector statevector = *refvector;

      if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "Propagating node " << **inode << endmsg; }

      StatusCode thissc = extrap->propagate( statevector, z, geometry, &F );
      if ( thissc.isFailure() ) {
        if ( m_debugLevel )
          debug() << "unable to propagate reference vector from z=" << refvector->z() << " to " << z
                  << "; track type = " << track.type() << ": vec = " << refvector->parameters() << endmsg;
        sc = thissc;
      }

      // correct for energy loss
      auto dE = node->deltaEnergy();
      if ( std::abs( statevector.qOverP() ) > LHCb::Math::lowTolerance ) {
        auto charge = statevector.qOverP() > 0 ? 1. : -1.;
        auto momnew = std::max( m_minMomentumForELossCorr.value(), std::abs( 1 / statevector.qOverP() ) + dE );
        if ( std::abs( momnew ) > m_minMomentumForELossCorr.value() ) statevector.setQOverP( charge / momnew );
      }

      // calculate the 'transport vector' (need to replace that)
      Gaudi::TrackVector tranportvec = statevector.parameters() - F * refvector->parameters();
      node->setTransportMatrix( F );
      node->setTransportVector( tranportvec );

      // update the reference
      refvector = &( ( *inode )->refVector() );

      // test dtx/dqop to see if the momentum affects this track.
      if ( std::abs( F( 2, 4 ) ) != 0 ) hasMomentum = true;

      if ( msgLevel( MSG::VERBOSE ) ) { verbose() << "Calculated node " << **inode << endmsg; }
    }
  }

  if ( m_useSeedStateErrors.value() ) {
    // we need to do this until we can properly deal with the seed state
    fitresult->setNTrackParameters( 0 );
  } else {
    fitresult->setNTrackParameters( hasMomentum ? 5 : 4 );
  }

  if ( m_debugLevel ) debug() << "End of TrackMasterFitter::updateTransport" << endmsg;
  return sc;
}

StatusCode TrackMasterFitter::initializeRefStates( LHCb::Track& track, IGeometryInfo const& geometry,
                                                   LHCb::Tr::PID pid ) const {
  if ( m_debugLevel ) debug() << "TrackMasterFitter::initializeRefStates" << endmsg;

  // given existing states on the track, this tool adds states at fixed
  // z-positions along the track. if a track state already exists
  // sufficiently close to the desired state, it will not add the
  // state.
  StatusCode sc = StatusCode::SUCCESS;

  // first fix the momentum of states on the track. need to make sure this works for Velo-TT as well.
  if ( track.states().empty() ) { return Error( "Track has no state! Can not fit.", StatusCode::FAILURE ); }
  // first need to make sure all states already on track have
  // reasonable momentum. still needs to check that this works for
  // velo-TT
  const LHCb::State* stateAtT = track.stateAt( LHCb::State::Location::AtT );
  const LHCb::State& refstate =
      stateAtT ? *stateAtT : *( track.isVeloBackward() ? track.states().front() : track.states().back() );
  for ( auto* state : track.states() ) const_cast<LHCb::State*>( state )->setQOverP( refstate.qOverP() );

  // collect the z-positions where we want the states
  boost::container::static_vector<double, 4> zpositions;
  if ( track.hasT() || track.hasMuon() ) {
    zpositions.push_back( StateParameters::ZBegT );
    zpositions.push_back( StateParameters::ZEndT );
  }
  if ( track.hasUT() || track.hasMuon() || ( track.hasT() && track.hasVelo() ) )
    zpositions.push_back( StateParameters::ZEndTT );
  if ( track.hasVelo() || track.hasMuon() ) zpositions.push_back( StateParameters::ZEndVelo );

  // the following container is going to hold pairs of 'desired'
  // z-positionds and actual states. the reason for the gymnastics
  // is that we always want to propagate from the closest availlable
  // state, but then recursively. this will make the parabolic
  // approximation reasonably accurate.
  typedef std::pair<double, const LHCb::State*> ZPosWithState;
  std::vector<ZPosWithState>                    states;
  states.reserve( track.states().size() );
  // we first add the states we already have
  const auto& tstates = track.states();
  std::transform( tstates.begin(), tstates.end(), std::back_inserter( states ), []( const LHCb::State* s ) {
    return std::pair{s->z(), s};
  } );

  // now add the other z-positions, provided nothing close exists
  const double maxDistance = 50 * Gaudi::Units::cm;
  for ( auto z : zpositions ) {
    bool not_found = std::none_of( states.begin(), states.end(),
                                   [&]( const ZPosWithState& s ) { return std::abs( z - s.first ) < maxDistance; } );
    if ( not_found ) states.emplace_back( z, nullptr );
  }
  std::sort( states.begin(), states.end(), FirstLess );

  // create the states in between
  const ITrackExtrapolator*   extrap = extrapolator( track.type() );
  LHCb::Track::StateContainer newstates;
  for ( auto it = states.begin(); it != states.end() && sc.isSuccess(); ++it ) {
    if ( it->second ) continue;

    // find the nearest existing state to it
    auto best = states.end();
    for ( auto jt = states.begin(); jt != states.end(); ++jt )
      if ( it != jt && jt->second &&
           ( best == states.end() || std::abs( jt->first - it->first ) < std::abs( best->first - it->first ) ) )
        best = jt;

    assert( best != states.end() );

    // move from that state to this iterator, using the extrapolator and filling all states in between.
    int               direction = best > it ? -1 : +1;
    LHCb::StateVector statevec( best->second->stateVector(), best->second->z() );
    for ( auto jt = best + direction; jt != it + direction && sc.isSuccess(); jt += direction ) {
      StatusCode thissc = extrap->propagate( statevec, jt->first, geometry, 0, pid );
      if ( !thissc.isSuccess() ) {
        Warning( "initializeRefStates() fails in propagating state", StatusCode::FAILURE ).ignore();
        if ( m_debugLevel ) debug() << "Problem propagating state: " << statevec << " to z= " << jt->first << endmsg;
        sc = thissc;
      } else {
        newstates.push_back( new LHCb::State( statevec ) );
        jt->second = newstates.back();
      }
    }
  }
  // finally, insert the new states to the track.
  if ( sc.isSuccess() ) {
    /// Exceptional unodered states: SeedMuon tracks used for TrackingEff study
    const auto order = track.checkType( LHCb::Event::Enum::Track::Type::SeedMuon ) ? LHCb::Tag::State::AssumeUnordered
                                                                                   : LHCb::Tag::State::AssumeSorted;
    track.addToStates( newstates, order );
  } else {
    for ( LHCb::State* state : newstates ) delete state;
  }

  return sc;
}
