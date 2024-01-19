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
#ifndef LINEARSTATEZTRAJ_H
#define LINEARSTATEZTRAJ_H

#include "Event/ZTrajectory.h"

#include "GaudiKernel/SystemOfUnits.h"

namespace LHCb {
  // ==========================================================================
  /** @class LinearStateZTraj
   *
   *  ZTrajectory class for a simple line parametrized by a State. The
   *  difference with StateZTraj is that this one cannot account for
   *  bending in a field. However, it does provide a correct covariance
   *  matrix in the call to 'state(z)', which unfortunately doesn't.
   *
   *  @author Wouter Hulsbergen
   *  @date 2018-05-25
   */
  template <typename FTYPE = double, unsigned int NBELMT = 1>
  class LinearStateZTraj : public ZTrajectory<FTYPE, NBELMT> {
  private:
    const LHCb::State m_state; /// state on which trajectory is based
  public:
    /// Typedefs
    using Point  = typename LHCb::ZTrajectory<FTYPE, NBELMT>::Point;
    using Vector = typename LHCb::ZTrajectory<FTYPE, NBELMT>::Vector;
    /// Constructor from a state
    LinearStateZTraj( const LHCb::State& state ) : m_state{state} {}
    /// Clone
    std::unique_ptr<LHCb::Trajectory<FTYPE>> clone() const override final {
      return std::unique_ptr<LHCb::Trajectory<FTYPE>>{new LinearStateZTraj{*this}};
    }
    /// Return the cached state
    const LHCb::State& state() const { return m_state; }
    /// Position at location z
    Point position( FTYPE z ) const override final {
      return Point{m_state.x() + ( z - m_state.z() ) * m_state.tx(), m_state.y() + ( z - m_state.z() ) * m_state.ty(),
                   1.0};
    }
    /// First derivative of position to z
    Vector direction( FTYPE /*z*/ ) const override final { return Vector{m_state.tx(), m_state.ty(), 1.0}; }
    /// Second derivative of position to z
    Vector curvature( FTYPE /*z*/ ) const override final { return Vector{0, 0, 0}; }
    /// Distance in z until the deviation from the linear
    /// approximation differs from this trajectory by a given tolerance.
    FTYPE distTo1stError( FTYPE /*z*/, FTYPE /*tolerance*/, int /*pathDirection*/ ) const override final {
      return 10 * Gaudi::Units::km;
    }
    /// Distance in z until the deviation from the quadratic
    /// approximation differs from this trajectory by a given tolerance.
    FTYPE distTo2ndError( FTYPE /*z*/, FTYPE /*tolerance*/, int /*pathDirection*/ ) const override final {
      return 10 * Gaudi::Units::km;
    }
    /// Create a parabolic approximation to the trajectory
    void expansion( FTYPE z, Point& p, Vector& dp, Vector& ddp ) const override final {
      p   = position( z );
      dp  = direction( z );
      ddp = curvature( z );
    }
    /// Arclength of total trajectory
    using ZTrajectory<FTYPE, NBELMT>::arclength;
    /// Arclength between 2 z -locations
    FTYPE arclength( FTYPE z1, FTYPE z2 ) const override final { return ( z2 - z1 ) * direction( z1 ).R(); }
    /// Estimate for expansion parameter 'z' closest to point
    FTYPE muEstimate( const Point& p ) const override final { return p.z(); }
    /// return a state at position z
    virtual LHCb::State state( FTYPE z ) const override final {
      LHCb::State rc = m_state;
      rc.linearTransportTo( z );
      return rc;
    }
    /// return a state vector at position z
    virtual LHCb::StateVector stateVector( FTYPE z ) const override final {
      const double dz = z - m_state.z();
      return LHCb::StateVector{Gaudi::TrackVector{m_state.x() + dz * m_state.tx(), m_state.y() + dz * m_state.ty(),
                                                  m_state.tx(), m_state.ty(), m_state.qOverP()},
                               z};
    }
    /// return the set of reference statevectors for this parameterization (if any)
    virtual std::vector<StateVector> refStateVectors() const override final {
      return std::vector<StateVector>{StateVector{m_state.stateVector(), m_state.z()}};
    }
  };
} // namespace LHCb

#endif
