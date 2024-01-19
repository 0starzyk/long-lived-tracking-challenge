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
#include "Event/FitNode.h"
#include "Event/KalmanFitResult.h"
#include "Event/Track.h"
#include "Event/TrackUnitsConverters.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbMath/MatrixManip.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackInterfaces/ITrackInterpolator.h"
#include "TrackKernel/TrackFunctors.h"

using namespace Gaudi;
using namespace Gaudi::Math;
using namespace LHCb;

/** @class TrackInterpolator TrackInterpolator.h
 *
 *  This tool finds the two nearest nodes and interpolates between the nodes
 *  to get the best estimate of an intermediate state at the given z-position.
 *  It extrapolates the two filtered states to the intermediate z-position and
 *  calculated the weighted mean.
 *  The current implemtation also applies the Kalman filter step because only
 *  the result from the prediction step is stored in the node (not the result
 *  of the filtered step).
 *
 *  @author Jeroen van Tilburg
 *  @date   2006-10-06
 */

class TrackInterpolator : public extends<GaudiTool, ITrackInterpolator> {
public:
  /// Standard constructor
  using extends::extends;

  /// Interpolate between the two nearest nodes to get a state
  StatusCode interpolate( const LHCb::Track& track, double z, LHCb::State& state,
                          IGeometryInfo const& geometry ) const override;

private:
  /// extrapolator
  ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackMasterExtrapolator"};

  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_inverse_failure{this, "Failure inverting matrix in smoother",
                                                                          0};
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_propagate_upstream_failure{
      this, "Failure propagating upstream state", 0};
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_propagate_downstream_failure{
      this, "Failure propagating downstream state", 0};
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_out_of_range{
      this, "Failure extrapolating outside measurement range", 0};
};

//-----------------------------------------------------------------------------
// Implementation file for class : TrackInterpolator
//
// 2006-10-06 : Jeroen van Tilburg
//-----------------------------------------------------------------------------

DECLARE_COMPONENT( TrackInterpolator )

//=============================================================================
// Interpolate between the nearest nodes
//=============================================================================
StatusCode TrackInterpolator::interpolate( const Track& track, double z, State& state,
                                           IGeometryInfo const& geometry ) const {
  // Check that track was actally fitted. Otherwise quietly call
  // extrapolator.
  const LHCb::TrackFitResult* fr = fitResult( track );
  if ( !fr || fr->nodes().empty() ) return m_extrapolator->propagate( track, z, state, geometry );

  const auto& nodes = fr->nodes();

  // If we are between the first and last node with a measurement, we
  // interpolate. If not, we extrapolate from the closest 'inside'
  // node. (That's more stable than interpolation.) In the end this
  // needs to work both for upstream and downstream fits. I am not
  // sure that it works for either now.

  // first find the pair of iterators such that z is between 'prevnode' and 'nextnode'
  auto nextnode = nodes.begin();
  if ( nodes.front()->z() < nodes.back()->z() ) {
    nextnode = std::find_if( nextnode, nodes.end(), [z]( const auto& n ) { return n->z() >= z; } );
  } else {
    nextnode = std::find_if( nextnode, nodes.end(), [z]( const auto& n ) { return n->z() <= z; } );
  }

  // determine where we are wrt to nodes with (active) measurements
  // is there measurement in nodes < nextnode?
  bool foundprecedingmeasurement = std::any_of(
      nodes.begin(), nextnode, []( const auto& n ) { return n->type() == LHCb::FitNode::Type::HitOnTrack; } );
  // is there a measurement in nodes >= nextnode?
  bool foundprocedingmeasurement = std::any_of(
      nextnode, nodes.end(), []( const auto& n ) { return n->type() == LHCb::FitNode::Type::HitOnTrack; } );

  // we must find either of the two (there must be measurement nodes!)
  if ( !foundprecedingmeasurement && !foundprocedingmeasurement )
    return Error( "Did not find any measurement nodes on track!" );

  // This is not necessarily a valid iterator, but that should be
  // caught by the logic later on.
  auto prevnode = std::prev( nextnode );

  // interpolate only if we have measurements on both sides
  if ( !foundprecedingmeasurement || !foundprocedingmeasurement ) {
    const LHCb::FitNode* extrapolationnode = foundprocedingmeasurement ? *nextnode : *prevnode;
    state                                  = extrapolationnode->state();
    return m_extrapolator->propagate( state, z, geometry ).orElse( [&] {
      if ( msgLevel( MSG::DEBUG ) )
        debug() << "Failure with normal extrapolator: z_target = " << z << " track type = " << track.type() << std::endl
                << "state = " << extrapolationnode->state() << endmsg;
      ++m_out_of_range;
    } );
  }

  if ( ( z - ( *nextnode )->z() ) * ( z - ( *prevnode )->z() ) > 0 ) {
    error() << "logic failure in locating nodes: " << z << ", " << ( *prevnode )->z() << "," << ( *nextnode )->z()
            << endmsg;
    return StatusCode::FAILURE;
  }

  // bail out if we have actually reached our destination
  for ( const auto* node : {*nextnode, *prevnode} ) {
    if ( std::abs( node->z() - z ) < TrackParameters::propagationTolerance ) {
      state = node->state();
      return StatusCode::SUCCESS;
    }
  }

  // Get the filtered states
  State stateDown, stateUp;

  LHCb::FitNode* fnextnode = *nextnode;
  LHCb::FitNode* fprevnode = *prevnode;
  if ( fnextnode && fprevnode ) {
    // TrackMasterFit
    stateDown = fprevnode->filteredStateForward();
    stateUp   = fnextnode->filteredStateBackward();
  }

  // extrapolate the upstream and downstream states
  auto sc = m_extrapolator->propagate( stateDown, z, geometry );
  if ( sc.isFailure() ) {
    if ( msgLevel( MSG::DEBUG ) )
      debug() << "Error propagating downstream state to z = " << z << std::endl << "state = " << stateDown << endmsg;
    ++m_propagate_downstream_failure;
    return sc;
  }

  sc = m_extrapolator->propagate( stateUp, z, geometry );
  if ( sc.isFailure() ) {
    if ( msgLevel( MSG::DEBUG ) )
      debug() << "Error propagating upstream state to z = " << z << " tracktype = " << track.type() << std::endl
              << "state = " << stateUp << endmsg;
    ++m_propagate_upstream_failure;
    return sc;
  }

  // Get the predicted downstream state and invert the covariance matrix
  const TrackVector& stateDownX    = stateDown.stateVector();
  TrackSymMatrix     invStateDownC = stateDown.covariance();
  if ( !invStateDownC.InvertChol() ) {
    ++m_inverse_failure;
    return StatusCode::FAILURE;
  }

  // Get the predicted upstream state and invert the covariance matrix
  const TrackVector& stateUpX    = stateUp.stateVector();
  TrackSymMatrix     invStateUpC = stateUp.covariance();
  if ( !invStateUpC.InvertChol() ) {
    ++m_inverse_failure;
    return StatusCode::FAILURE;
  }

  // Add the inverted matrices
  TrackSymMatrix& stateC = state.covariance();
  stateC                 = invStateDownC + invStateUpC;
  if ( !stateC.InvertChol() ) {
    ++m_inverse_failure;
    return StatusCode::FAILURE;
  }

  // Get the state by calculating the weighted mean
  TrackVector& stateX = state.stateVector();
  stateX              = stateC * ( ( invStateDownC * stateDownX ) + ( invStateUpC * stateUpX ) );
  state.setZ( z );

  if ( msgLevel( MSG::DEBUG ) )
    debug() << "filteredstate A: " << stateUpX << std::endl
            << "filteredstate B: " << stateDownX << std::endl
            << "smoothed state A: " << ( *prevnode )->state() << "smoothed state B: " << ( *nextnode )->state()
            << "interpolated state: " << state << endmsg;

  return StatusCode::SUCCESS;
}
