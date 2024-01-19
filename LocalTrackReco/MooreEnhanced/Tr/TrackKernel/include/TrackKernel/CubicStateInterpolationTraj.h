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
#ifndef TRACKKERNEL_CubicStateInterpolationTraj_H
#define TRACKKERNEL_CubicStateInterpolationTraj_H 1

#include "Event/State.h"
#include "TrackKernel/CubicStateVectorInterpolationTraj.h"
#include <memory>

// Include files

/** @class CubicStateInterpolationTraj CubicStateInterpolationTraj.h TrackKernel/CubicStateInterpolationTraj.h
 *
 * Cubic interpolating between two states. Used for TrackTraj.
 *
 *  @author Wouter Hulsbergen
 *  @date   2007-10-08
 */

namespace LHCb {
  class CubicStateInterpolationTraj : public CubicStateVectorInterpolationTraj {
  public:
    /// Constructor from two states
    CubicStateInterpolationTraj( const LHCb::State& begin, const LHCb::State& end );

    /// Constructor from two states and a bfield
    CubicStateInterpolationTraj( const LHCb::State& state, const Gaudi::XYZVector& bfield );

    /// Default constructor (does nothing)
    CubicStateInterpolationTraj() = default;

    /// Initialize with two states
    void init( const LHCb::State& begin, const LHCb::State& end );

    /// Initialize with single state and b field
    void init( const LHCb::State& state, const Gaudi::XYZVector& bfield );

    /// Covariance of state vector at location z
    Gaudi::TrackSymMatrix covariance( double z ) const;

    /// Calculate state at location z
    State state( double z ) const override final;

    /// Clone this trajectory
    std::unique_ptr<Trajectory> clone() const override {
      return std::make_unique<CubicStateInterpolationTraj>( *this );
    }

  private:
    Gaudi::TrackSymMatrix m_covbegin; ///< Covariance at first state
    Gaudi::TrackSymMatrix m_covend;   ///< Covariance at second state
  };

  /*************************************************************************************************/
  // inline functions
  /*************************************************************************************************/

  inline CubicStateInterpolationTraj::CubicStateInterpolationTraj( const LHCb::State& begin, const LHCb::State& end )
      : CubicStateVectorInterpolationTraj( begin, end )
      , m_covbegin( begin.covariance() )
      , m_covend( end.covariance() ) {}

  inline CubicStateInterpolationTraj::CubicStateInterpolationTraj( const LHCb::State&      state,
                                                                   const Gaudi::XYZVector& bfield )
      : CubicStateVectorInterpolationTraj( state, bfield ), m_covbegin( state.covariance() ), m_covend( m_covbegin ) {}

  inline void CubicStateInterpolationTraj::init( const LHCb::State& begin, const LHCb::State& end ) {
    CubicStateVectorInterpolationTraj::init( begin, end );
    m_covbegin = begin.covariance();
    m_covend   = end.covariance();
  }

  inline void CubicStateInterpolationTraj::init( const LHCb::State& state, const Gaudi::XYZVector& bfield ) {
    CubicStateVectorInterpolationTraj::init( state, bfield );
    m_covbegin = m_covend = state.covariance();
  }
} // namespace LHCb

#endif // EVENT_CUBICSTATEINTERPOLATION_H
