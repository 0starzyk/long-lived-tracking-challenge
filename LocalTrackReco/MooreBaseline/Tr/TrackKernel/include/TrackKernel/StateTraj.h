/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#ifndef TRACKFITEVENT_STATETRAJ_H
#define TRACKFITEVENT_STATETRAJ_H 1

// Include files
// -------------

// STL
#include <memory>

// from LHCbKernel
#include "Kernel/DifTraj.h"

// from TrackEvent
#include "Event/State.h"
#include "Event/StateVector.h"
#include "Event/TrackParameters.h"

namespace LHCb {

  /** @class StateTraj StateTraj.h TrackKernel/StateTraj.h
   *
   * Trajectory created from a State.
   *
   * @author Edwin Bos, Jeroen van Tilburg, Eduardo Rodrigues
   * @date   01/12/2005
   */

  class StateTraj : public DifTraj<5> {

  public:
    /// Enum providing number of colums in derivative matrix
    enum { kSize = 5 };

    /// get me another one of these!
    std::unique_ptr<Trajectory<double>> clone() const override;

    /// Default Destructor
    virtual ~StateTraj() {}

    /// Constructor from a State and the magnetic field at the State position
    StateTraj( const LHCb::State& state, const Gaudi::XYZVector& bField );

    /// Constructor from a StateVector and the magnetic field at the State position
    StateTraj( const LHCb::StateVector& state, const Gaudi::XYZVector& bField );

    /// Constructor from a StateVector and the magnetic field at State position
    StateTraj( const Gaudi::TrackVector& stateVector, double z, const Gaudi::XYZVector& bField );

    /// Point on trajectory where parabolic approximation is made
    virtual Gaudi::XYZPoint position( double arclength ) const override;

    /// First derivative of the trajectory at the approximation point
    virtual Gaudi::XYZVector direction( double arclength ) const override;

    /// Second derivative of the trajectory at the approximation point,
    /// used as the constant value of the curvature of the parabolic approximation
    virtual Gaudi::XYZVector curvature( double arclength ) const override;

    /// Create a parabolic approximation to the trajectory
    virtual void expansion( double arclength, Gaudi::XYZPoint& p, Gaudi::XYZVector& dp,
                            Gaudi::XYZVector& ddp ) const override;

    /// Retrieve the parameters of this traj...
    virtual Parameters parameters() const override;

    /// Update the parameters of this traj...
    virtual StateTraj& operator+=( const Parameters& delta ) override;

    /// Retrieve the derivative of the parabolic approximation to the
    /// trajectory with respect to the state parameters
    virtual Derivative derivative( double arclength ) const override;

    /// give arclength where this trajectory is closest to the
    /// specified point
    virtual double muEstimate( const Gaudi::XYZPoint& point ) const override;

    /// Number of arclengths until deviation of the trajectory from the expansion
    /// reaches the given tolerance (does not account for the curvature).
    virtual double distTo1stError( double arclength, double tolerance, int pathDirection = +1 ) const override;

    /// Number of arclengths until deviation of the trajectory from the
    /// expansion reaches the given tolerance (accounts for the curvature).
    virtual double distTo2ndError( double arclen, double tolerance, int pathDirection = +1 ) const override;

    using LHCb::DifTraj<5>::arclength;
    /// Distance, along the Trajectory, between position(mu1) and
    /// position(mu2). Trivial because StateTraj is parameterized in
    /// arclength.
    virtual double arclength( double mu1, double mu2 ) const override { return mu2 - mu1; }

  private:
    Gaudi::XYZPoint  m_pos;    ///< the position of the State
    Gaudi::XYZVector m_dir;    ///< the unit direction of the State
    double           m_qOverP; ///< the charge-over-momentum Q/P of the State
    Gaudi::XYZVector m_bField; ///< the magnetic field vector at the State position
    Gaudi::XYZVector m_curv;   ///< constant value of parabola's curvature

  }; // class StateTraj

} // namespace LHCb

#endif /// TRACKFITEVENT_STATETRAJ_H
