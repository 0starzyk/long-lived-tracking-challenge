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
#ifndef TRACKKERNEL_ZTRAJECTORY_H
#define TRACKKERNEL_ZTRAJECTORY_H 1

// Include files
#include "Event/State.h"
#include "Event/StateVector.h"
#include "Kernel/Trajectory.h"
#include <vector>

/** @class ZTrajectory ZTrajectory.h TrackKernel/ZTrajectory.h
 *  Interface for trajectories parameterized along Z. Gives access to statevectors.
 *
 *  @author Wouter HULSBERGEN
 *  @date   2007-11-29
 */

namespace LHCb {
  // TODO : Would be nice to be able to template the floating point type
  //        for State and StateVector
  template <typename FTYPE>
  struct StateVectorT {
    using Vector5 = typename ROOT::Math::SVector<FTYPE, 5>;

    Vector5 parameters;
    FTYPE   z;
  };

  template <typename FTYPE = double, unsigned int NBELMT = 1>
  class ZTrajectory : public Trajectory<FTYPE> {
  public:
    using Trajectory<FTYPE>::Trajectory;
    using Vector      = typename Trajectory<FTYPE>::Vector;
    using Point       = typename Trajectory<FTYPE>::Point;
    using Range       = typename Trajectory<FTYPE>::Range;
    using StateVector = typename std::conditional<std::is_floating_point<FTYPE>::value, LHCb::StateVector,
                                                  LHCb::StateVectorT<FTYPE>>::type;
    using AdaptState  = typename std::conditional<NBELMT == 1, LHCb::State, std::array<State, NBELMT>>::type;

    template <size_t... Is>
    constexpr AdaptState toStates( StateVector sv, std::index_sequence<Is...> ) const {
      if constexpr ( !std::is_floating_point<FTYPE>::value ) {
        alignas( 64 ) float data[6 * NBELMT];
        for ( unsigned int i = 0; i < 5; i++ ) sv.parameters[i].store( data + i * NBELMT );
        sv.z.store( data + 5 * NBELMT );
        return {State( LHCb::StateVector(
            {data[Is], data[NBELMT + Is], data[2 * NBELMT + Is], data[3 * NBELMT + Is], data[4 * NBELMT + Is]},
            data[5 * NBELMT + Is] ) )...};
      } else {
        return {State( sv )};
      }
    }

    /// Default constructor
    ZTrajectory() : Trajectory<FTYPE>( 0., 0. ) {}
    /// Constructor taking the values of mu that defined the valid range of the trajectory
    ZTrajectory( FTYPE begin, FTYPE end ) : Trajectory<FTYPE>( begin, end ) {}
    /// Constructor taking a range
    ZTrajectory( const Range& range ) : Trajectory<FTYPE>( range ) {}
    /// return stateVector at position mu
    virtual StateVector stateVector( FTYPE mu ) const = 0;
    /// return a state at position mu
    virtual AdaptState state( FTYPE mu ) const {
      return toStates( stateVector( mu ), std::make_index_sequence<NBELMT>() );
    }
    /// return the set of reference statevectors for this parameterization (if any)
    virtual std::vector<StateVector> refStateVectors() const { return {}; }
  };

} // namespace LHCb

#endif // EVENT_ZTRAJECTORY_H
