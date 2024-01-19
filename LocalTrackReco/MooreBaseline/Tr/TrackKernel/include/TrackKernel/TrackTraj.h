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
#pragma once

#include "Event/Track.h"
#include "Kernel/STLExtensions.h"
#include "Kernel/Trajectory.h"
#include "Magnet/DeMagnet.h"
#include "TrackKernel/CubicStateInterpolationTraj.h"
#include <algorithm>
#include <memory>
#include <vector>

namespace LHCb {

  class FitNode;

  /** @class TrackTraj TrackTraj.h TrackKernel/TrackTraj.h
   *
   * Trajectory representation of a LHCb::Track
   *
   * @author Wouter Hulsbergen
   * @date   15/10/2007
   */

  class TrackTraj : public ZTrajectory<double> {
  public:
    using ZTrajectory<double>::StateVector;

    /// Constructor from a track and (optionally) pointer to mag field service
    TrackTraj( const LHCb::Track& track, const DeMagnet* magfieldsvc = nullptr );

    /// Constructor from an unsorted set of states and (optionally) pointer to mag field service
    TrackTraj( span<const LHCb::State* const> states, const DeMagnet* magfieldsvc = nullptr );

    /// Constructor from a (sorted) Track::StateContainer and (optionally) pointer to mag field service
    TrackTraj( span<const LHCb::State* const> states, LHCb::Tag::State::AssumeSorted_tag,
               const DeMagnet*                magfieldsvc = nullptr );

    /// Constructor from a (sorted) std::vector<State> and (optionally) pointer to mag field service
    TrackTraj( span<const LHCb::State> states, LHCb::Tag::State::AssumeSorted_tag,
               const DeMagnet*         magfieldsvc = nullptr );

    /// Constructor from a (sorted) Track::NodesContainer and (optionally) pointer to mag field service
    TrackTraj( span<const FitNode* const> nodes, const DeMagnet* magfieldsvc = nullptr );

#if defined( __GNUC__ ) && ( __GNUC__ < 10 )
    // cppgsl 3.x gcc 9.x workaround. To be removed when gcc 9 or older no longer supported
    TrackTraj( const std::vector<LHCb::State*>& states, const DeMagnet* magfieldsvc = nullptr )
        : TrackTraj(
              span<const LHCb::State* const>{states.data(),
                                             static_cast<span<const LHCb::State* const>::size_type>( states.size() )},
              magfieldsvc ) {}
    TrackTraj( const std::vector<LHCb::State*>& states, LHCb::Tag::State::AssumeSorted_tag tag,
               const DeMagnet* magfieldsvc = nullptr )
        : TrackTraj(
              span<const LHCb::State* const>{states.data(),
                                             static_cast<span<const LHCb::State* const>::size_type>( states.size() )},
              tag, magfieldsvc ) {}
#endif

    /// Clone method
    std::unique_ptr<Trajectory<double>> clone() const override { return std::make_unique<TrackTraj>( *this ); }

    /// Point on the trajectory at arclength from the starting point
    Point position( double z ) const override final;

    /// First derivative of the trajectory at arclength from the starting point
    Vector direction( double z ) const override final;

    /// Second derivative of the trajectory at arclength from the starting point
    Vector curvature( double z ) const override final;

    /// Point on the trajectory at arclength from the starting point
    Vector momentum( double z ) const;

    /// State at given z
    State state( double z ) const override final;

    /// State at given z
    StateVector stateVector( double z ) const override final;

    /// Expand this track in z
    void expansion( double z, Point& p, Vector& dp, Vector& ddp ) const override final;

    /// distance where the deviation of the trajectory from the expansion
    /// reaches the given tolerance.
    double distTo1stError( double mu, double tolerance, int pathDirection = +1 ) const override final;

    /// distance where the deviation of the trajectory from the expansion
    /// reaches the given tolerance.
    double distTo2ndError( double mu, double tolerance, int pathDirection = +1 ) const override final;

    using ZTrajectory<double>::arclength;
    /// arclength between 2 coordinates on the track
    double arclength( double z1, double z2 ) const override final;

    /// Derivative of arclength to mu
    double dArclengthDMu( double z ) const;

    /// Estimate for mu which minimizes point poca
    double muEstimate( const Gaudi::XYZPoint& p ) const override final;

    /// return the set of reference states
    const auto& refStates() const { return m_states; }

    /// return the set of reference statevectors for this parameterization (if any)
    std::vector<StateVector> refStateVectors() const override final;

  protected:
    /// return the set of reference states
    auto& refStates() { return m_states; }

    /// invalidate the cache
    void invalidateCache() { m_cachedindex = InvalidCacheIndex; }

  private:
    void updatecache( double z ) const;
    void init( const DeMagnet* magfieldsvc );

  private:
    enum : unsigned int { InvalidCacheIndex = (unsigned int)( -1 ) };
    std::vector<const LHCb::State*>     m_states;              ///< Container of states
    Gaudi::XYZVector                    m_bfield;              ///< Bfield at upstream end of track
    mutable size_t                      m_cachedindex;         ///< Index for cached z-range
    mutable CubicStateInterpolationTraj m_cachedinterpolation; ///< Cached interpolation for z-range
  };

  /*************************************************************************************************/
  // inline functions
  /*************************************************************************************************/

  inline Trajectory<double>::Point TrackTraj::position( double z ) const {
    updatecache( z );
    return m_cachedinterpolation.position( z );
  }

  inline Trajectory<double>::Vector TrackTraj::momentum( double z ) const {
    updatecache( z );
    return m_cachedinterpolation.momentum( z );
  }

  inline Trajectory<double>::Vector TrackTraj::direction( double z ) const {
    updatecache( z );
    return m_cachedinterpolation.direction( z );
  }

  inline Trajectory<double>::Vector TrackTraj::curvature( double z ) const {
    updatecache( z );
    return m_cachedinterpolation.curvature( z );
  }

  inline LHCb::State TrackTraj::state( double z ) const {
    updatecache( z );
    return m_cachedinterpolation.CubicStateInterpolationTraj::state( z );
  }

  inline LHCb::StateVector TrackTraj::stateVector( double z ) const {
    updatecache( z );
    return m_cachedinterpolation.CubicStateInterpolationTraj::stateVector( z );
  }

  inline void TrackTraj::expansion( double z, Trajectory<double>::Point& p, Trajectory<double>::Vector& dp,
                                    Trajectory<double>::Vector& ddp ) const {
    updatecache( z );
    return m_cachedinterpolation.expansion( z, p, dp, ddp );
  }

  inline double TrackTraj::muEstimate( const Gaudi::XYZPoint& p ) const {
    updatecache( p.z() );
    return m_cachedinterpolation.muEstimate( p );
  }

  inline double TrackTraj::dArclengthDMu( double z ) const {
    updatecache( z );
    return m_cachedinterpolation.dArclengthDMu( z );
  }

} // namespace LHCb
