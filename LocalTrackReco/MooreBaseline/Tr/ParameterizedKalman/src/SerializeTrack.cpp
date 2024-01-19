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
#include "SerializeTrack.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include <array>

// How many bytes a serialized state takes.
// 6 doubles for x,y,z,tx,t,z,qop and 5+4+..+1 for the symmetric covariance matrix.
constexpr int BytesPerState = sizeof( double ) * 6 + sizeof( double ) * 15;

using SerializedState = std::array<uint8_t, BytesPerState>;

namespace {
  void push( uint8_t*& ptr, double d ) {
    std::memcpy( ptr, &d, sizeof( double ) );
    ptr += sizeof( double );
  }
  double pop( uint8_t const*& ptr ) {
    double d = 0;
    std::memcpy( &d, ptr, sizeof( double ) );
    ptr += sizeof( double );
    return d;
  }
} // namespace

// Serialize one single state.
SerializedState SerializeState( const LHCb::State& state ) {
  SerializedState serializedState;
  auto            ptr = serializedState.begin();

  // Store state vector.
  push( ptr, state.x() );
  push( ptr, state.y() );
  push( ptr, state.z() );
  push( ptr, state.tx() );
  push( ptr, state.ty() );
  push( ptr, state.qOverP() );

  // Store covariance matrix.
  const auto& covMat = state.covariance();
  for ( int i = 0; i < 5; ++i ) {
    for ( int j = i; j < 5; ++j ) { push( ptr, covMat( i, j ) ); }
  }
  assert( ptr = serializedState.end() );
  return serializedState;
}

// Deserialize one single state.
LHCb::State DeserializeState( const SerializedState& serializedState ) {

  auto ptr = serializedState.begin();

  LHCb::State state;
  // Load state vector.
  double x      = pop( ptr );
  double y      = pop( ptr );
  double z      = pop( ptr );
  double tx     = pop( ptr );
  double ty     = pop( ptr );
  double qOverP = pop( ptr );
  state.setState( x, y, z, tx, ty, qOverP );

  // Load covariance matrix.
  Gaudi::SymMatrix5x5 covMat;
  for ( int i = 0; i < 5; ++i ) {
    for ( int j = i; j < 5; ++j ) { covMat( i, j ) = pop( ptr ); }
  }
  state.setCovariance( covMat );

  return state;
}

std::vector<uint8_t> SerializeStates( const std::vector<LHCb::State>& states ) {
  std::vector<uint8_t> binary;
  binary.reserve( states.size() * BytesPerState );

  for ( auto& state : states ) {
    SerializedState singleState = SerializeState( state );
    binary.insert( binary.end(), singleState.begin(), singleState.end() );
  }

  return binary;
}

std::vector<LHCb::State> DeserializeStates( const std::vector<uint8_t>& binary ) {
  std::vector<LHCb::State> states;
  size_t                   numStates = binary.size() / BytesPerState;
  assert( binary.size() % BytesPerState == 0 );
  states.reserve( numStates );

  auto binaryIt = binary.begin();
  for ( size_t i = 0; i < numStates; ++i ) {
    SerializedState singleState;
    for ( auto& v : singleState ) {
      v = *binaryIt;
      ++binaryIt;
    }
    states.push_back( DeserializeState( singleState ) );
  }

  return states;
}

StateVectorDifference CompareSingleStates( const LHCb::State& lhs, const LHCb::State& rhs ) {
  StateVectorDifference diff;

  // Compute distance of the two points.
  Gaudi::XYZPoint p1 = {lhs.x(), lhs.y(), lhs.z()};
  Gaudi::XYZPoint p2 = {rhs.x(), rhs.y(), rhs.z()};

  diff.avgPos = diff.maxPos = sqrt( ( p1 - p2 ).mag2() );

  // Compute L2 norm of the difference of the covariance matrices.
  double sumCov     = 0;
  auto   covMatDiff = lhs.covariance() - rhs.covariance();

  for ( int i = 0; i < 5; ++i ) {
    for ( int j = i; j < 5; ++j ) {
      // Divide diagonal by 2.
      sumCov += covMatDiff( i, j ) * covMatDiff( i, j ) / ( 1 + i == j );
    }
  }

  diff.avgCov = diff.maxCov = sqrt( sumCov * 2.0 / 25.0 );

  return diff;
}

StateVectorDifference CompareStates( const std::vector<LHCb::State>& lhs, const std::vector<LHCb::State>& rhs ) {
  StateVectorDifference diff;

  // Comparison makes no sense if sizes are not the same.
  // Possible to sort by Z and try to match states, but for now it doesn't matter.
  if ( lhs.size() != rhs.size() ) {
    diff.numStates = std::abs( intptr_t( lhs.size() ) - intptr_t( rhs.size() ) );
    return diff;
  }

  diff.numStates = 0;

  // Compare states one-by-one.
  for ( size_t i = 0; i < lhs.size(); ++i ) {
    StateVectorDifference singleDiff = CompareSingleStates( lhs[i], rhs[i] );
    diff.maxPos                      = std::max( diff.maxPos, singleDiff.maxPos );
    diff.maxCov                      = std::max( diff.maxCov, singleDiff.maxCov );
    diff.avgPos += singleDiff.avgPos;
    diff.avgCov += singleDiff.avgCov;
  }

  diff.avgPos /= lhs.size();
  diff.avgCov /= lhs.size();

  return diff;
}
