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
#ifndef TRACKKERNEL_STATEZTRAJ_H
#define TRACKKERNEL_STATEZTRAJ_H 1

// Include files
// -------------

// STL
#include <cmath>
#include <memory>
#include <type_traits>

// from LHCbKernel
#include "Event/ZTrajectory.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Kernel/DifTraj.h"

namespace LHCb {
  template <typename FTYPE, unsigned int NBELMT>
  inline auto upderiv( FTYPE deriv, FTYPE tolerance ) {
    if constexpr ( !std::is_floating_point<FTYPE>::value ) {
      return select( deriv == 0., FTYPE{10 * Gaudi::Units::km}, sqrt( abs( 2 * tolerance / deriv ) ) );
    } else {
      return deriv != 0. ? sqrt( fabs( 2 * tolerance / deriv ) ) : 10 * Gaudi::Units::km;
    }
  }

  /** @class StateZTraj StateZTraj.h TrackKernel/StateZTraj.h
   *
   * Trajectory created from a State, parameterized in z. This still needs some
   * work. It is a DifTraj, but because I cannot use MI, it actually does not
   * derive from DifTraj. Since DifTraj is a tenmplated base anyway, it wouldn't
   * help either, so nobody will care.
   *
   * @author Wouter Hulsbergen (after StateTraj by Edwin Bos, Jeroen van Tilburg, Eduardo Rodrigues)
   * @date   15/10/2007
   */
  template <typename FTYPE = double, unsigned int NBELMT = 1>
  class StateZTraj : public ZTrajectory<FTYPE, NBELMT> {
  public:
    /// Enum providing number of colums in derivative matrix
    enum { kSize = 5 };

    using ZTrajectory<FTYPE, NBELMT>::ZTrajectory;
    using Vector      = typename ZTrajectory<FTYPE, NBELMT>::Vector;
    using Vector5     = typename ROOT::Math::SVector<FTYPE, 5>;
    using Point       = typename ZTrajectory<FTYPE, NBELMT>::Point;
    using Range       = typename ZTrajectory<FTYPE, NBELMT>::Range;
    using StateVector = typename ZTrajectory<FTYPE, NBELMT>::StateVector;

    // typedefs
    typedef ROOT::Math::SMatrix<FTYPE, 3, kSize> Derivative;
    typedef ROOT::Math::SVector<FTYPE, kSize>    Parameters;

    /// get me another one of these!
    std::unique_ptr<Trajectory<FTYPE>> clone() const override { return std::make_unique<StateZTraj>( *this ); }

    /// Default Destructor
    ~StateZTraj() {}

    /// Constructor from the magnetic field at the States position
    StateZTraj( FTYPE x, FTYPE y, FTYPE tx, FTYPE ty, FTYPE qop, FTYPE z, const Vector& bField = {0., 0., 0.} );

    /// Constructor from a State or StateVector and the magnetic field at the States position
    template <class StateT>
    StateZTraj( const StateT& state, const Vector& bField );

    /// Point on trajectory where parabolic approximation is made
    Point position( FTYPE z ) const override { return Point( x( z ), y( z ), z ); }

    /// First derivative of the trajectory at the approximation point
    Vector direction( FTYPE z ) const override { return Vector( tx( z ), ty( z ), 1 ); }

    /// Second derivative of the trajectory at the approximation point,
    /// used as the constant value of the curvature of the parabolic approximation
    Vector curvature( FTYPE z ) const override { return Vector( omegax( z ), omegay( z ), 0. ); }

    /// Create a parabolic approximation to the trajectory
    void expansion( FTYPE z, Point& p, Vector& dp, Vector& ddp ) const override;

    /// Retrieve the parameters of this traj...
    Parameters parameters() const;

    /// Update the parameters of this traj...
    StateZTraj<FTYPE, NBELMT>& operator+=( const Parameters& delta );

    /// Retrieve the derivative of the parabolic approximation to the
    /// trajectory with respect to the state parameters
    template <typename U = FTYPE>
    Derivative derivative( FTYPE z,
                           typename std::enable_if<not std::is_floating_point<U>::value>::type* = nullptr ) const {
      Derivative deriv{ROOT::Math::SMatrixNoInit()};
      for ( unsigned int i = 0; i < 3 * kSize; ++i ) deriv.Array()[i] = 0.;
      FTYPE dz      = z - m_z;
      deriv( 0, 0 ) = deriv( 1, 1 ) = 1.;
      deriv( 0, 2 ) = deriv( 1, 3 ) = dz;

      // to speed this up, we only calculate the rest if the trajectory is indeed curved.
      auto isnonlinear = abs( m_cx[2] ) > 1e-10 || abs( m_cy[2] ) > 1e-10;

      if ( !none( isnonlinear ) ) {
        FTYPE tx     = m_cx[1];
        FTYPE ty     = m_cy[1];
        FTYPE omegax = m_cx[2];
        FTYPE omegay = m_cy[2];
        FTYPE n      = sqrt( 1 + tx * tx + ty * ty );
        FTYPE dndtx  = tx / n;
        FTYPE dndty  = ty / n;

        FTYPE      half_dz_sqr = 0.5 * dz * dz;
        const auto omegn       = omegax / n;

        deriv( 0, 2 ) += select( isnonlinear, half_dz_sqr * omegn * dndtx, 0. );
        deriv( 0, 3 ) +=
            select( isnonlinear, half_dz_sqr * ( omegn * dndty + n * m_qOverP * Gaudi::Units::c_light * m_Bz ), 0. );
        deriv( 0, 4 ) += select( isnonlinear, half_dz_sqr * omegax / m_qOverP, 0. );
        deriv( 1, 2 ) +=
            select( isnonlinear, half_dz_sqr * ( omegn * dndtx - n * m_qOverP * Gaudi::Units::c_light * m_Bz ), 0. );
        deriv( 1, 3 ) += select( isnonlinear, half_dz_sqr * omegn * dndty, 0. );
        deriv( 1, 4 ) += select( isnonlinear, half_dz_sqr * omegay / m_qOverP, 0. );
      }

      return deriv;
    }

    /// Retrieve the derivative of the parabolic approximation to the
    /// trajectory with respect to the state parameters
    template <typename U = FTYPE>
    Derivative derivative( FTYPE z, typename std::enable_if<std::is_floating_point<U>::value>::type* = nullptr ) const {
      Derivative deriv;
      FTYPE      dz = z - m_z;
      deriv( 0, 0 ) = deriv( 1, 1 ) = 1;
      deriv( 0, 2 ) = deriv( 1, 3 ) = dz;

      // to speed this up, we only calculate the rest if the trajectory is indeed curved.
      auto isnonlinear = fabs( m_cx[2] ) > 1e-10 || fabs( m_cy[2] ) > 1e-10;
      if ( isnonlinear ) {
        FTYPE tx     = m_cx[1];
        FTYPE ty     = m_cy[1];
        FTYPE omegax = m_cx[2];
        FTYPE omegay = m_cy[2];
        FTYPE n      = std::sqrt( 1 + tx * tx + ty * ty );
        FTYPE dndtx  = tx / n;
        FTYPE dndty  = ty / n;

        FTYPE half_dz_sqr = 0.5 * dz * dz;
        FTYPE omegn       = omegax / n;

        deriv( 0, 2 ) += half_dz_sqr * omegn * dndtx;
        deriv( 0, 3 ) += half_dz_sqr * ( omegn * dndty + n * m_qOverP * Gaudi::Units::c_light * m_Bz );
        deriv( 0, 4 ) += half_dz_sqr * omegax / m_qOverP;

        deriv( 1, 2 ) += half_dz_sqr * ( omegn * dndtx - n * m_qOverP * Gaudi::Units::c_light * m_Bz );
        deriv( 1, 3 ) += half_dz_sqr * omegn * dndty;
        deriv( 1, 4 ) += half_dz_sqr * omegay / m_qOverP;
      }

      return deriv;
    }

    /// give arclength where this trajectory is closest to the
    /// specified point. (linear only. can be improved)
    FTYPE muEstimate( const Point& point ) const override;

    /// Number of arclengths until deviation of the trajectory from the expansion
    /// reaches the given tolerance (does not account for the curvature).
    FTYPE distTo1stError( FTYPE arclength, FTYPE tolerance, int pathDirection = +1 ) const override;

    /// Number of arclengths until deviation of the trajectory from the
    /// expansion reaches the given tolerance (accounts for the curvature).
    FTYPE distTo2ndError( FTYPE arclen, FTYPE tolerance, int pathDirection = +1 ) const override;

    using ZTrajectory<FTYPE, NBELMT>::arclength;
    /// Distance, along the Trajectory, between position(mu1) and
    /// position(mu2). Trivial because StateZTraj is parameterized in
    /// arclength.
    FTYPE arclength( FTYPE mu1, FTYPE mu2 ) const override { return mu2 - mu1; }

    /// return stateVector at position mu
    StateVector stateVector( FTYPE z ) const override {
      if constexpr ( !std::is_floating_point<FTYPE>::value ) {
        const Vector5 parameters( x( z ), y( z ), tx( z ), ty( z ), m_qOverP );
        return StateVector{parameters, z};
      } else {
        const Gaudi::Vector5 parameters( x( z ), y( z ), tx( z ), ty( z ), m_qOverP );
        return StateVector{parameters, z};
      }
    }

    FTYPE x( FTYPE z ) const { return polyeval( z - m_z, m_cx ); }
    FTYPE tx( FTYPE z ) const { return poly1stderiveval( z - m_z, m_cx ); }
    FTYPE omegax( FTYPE z ) const { return poly2ndderiveval( z - m_z, m_cx ); }
    FTYPE y( FTYPE z ) const { return polyeval( z - m_z, m_cy ); }
    FTYPE ty( FTYPE z ) const { return poly1stderiveval( z - m_z, m_cy ); }
    FTYPE omegay( FTYPE z ) const { return poly2ndderiveval( z - m_z, m_cy ); }
    FTYPE polyeval( FTYPE dz, const std::array<FTYPE, 3>& c ) const { return c[0] + dz * ( c[1] + dz * c[2] ); }
    FTYPE poly1stderiveval( FTYPE dz, const std::array<FTYPE, 3>& c ) const { return c[1] + 2 * dz * c[2]; }
    FTYPE poly2ndderiveval( FTYPE /*dz*/, const std::array<FTYPE, 3>& c ) const { return 2 * c[2]; }

    FTYPE                m_z;      ///< z-position of this state
    std::array<FTYPE, 3> m_cx;     ///< Coefficients for parabola x(z)
    std::array<FTYPE, 3> m_cy;     ///< Coefficients for parabola y(z)
    FTYPE                m_qOverP; ///< the charge-over-momentum Q/P of the State
    FTYPE                m_Bz;     ///< z-component of B field (needed to calculate derivative)
  };                               // class StateZTraj

  /*************************************************************************************************/
  // inline functions
  /*************************************************************************************************/

  template <typename FTYPE, unsigned int NBELMT>
  template <class StateT>
  StateZTraj<FTYPE, NBELMT>::StateZTraj( const StateT& state, const Vector& bfield )
      : ZTrajectory<FTYPE, NBELMT>(), m_z( state.z() ), m_qOverP( state.qOverP() ), m_Bz( bfield.z() ) {
    FTYPE n = sqrt( 1 + state.tx() * state.tx() + state.ty() * state.ty() );

    m_cx[0] = state.x();
    m_cx[1] = state.tx();
    m_cx[2] = n * Gaudi::Units::c_light * m_qOverP * ( -bfield.y() + state.ty() * bfield.z() );

    m_cy[0] = state.y();
    m_cy[1] = state.ty();
    m_cy[2] = n * Gaudi::Units::c_light * m_qOverP * ( bfield.x() - state.tx() * bfield.z() );
  }

  template <typename FTYPE, unsigned int NBELMT>
  StateZTraj<FTYPE, NBELMT>::StateZTraj( FTYPE x, FTYPE y, FTYPE tx, FTYPE ty, FTYPE qop, FTYPE z,
                                         const Vector& bfield )
      : ZTrajectory<FTYPE, NBELMT>(), m_z( z ), m_qOverP( qop ), m_Bz( bfield.z() ) {
    FTYPE n = sqrt( 1 + tx * tx + ty * ty );

    m_cx[0] = x;
    m_cx[1] = tx;
    m_cx[2] = n * Gaudi::Units::c_light * m_qOverP * ( -bfield.y() + ty * bfield.z() );

    m_cy[0] = y;
    m_cy[1] = ty;
    m_cy[2] = n * Gaudi::Units::c_light * m_qOverP * ( bfield.x() - tx * bfield.z() );
  }

  template <typename FTYPE, unsigned int NBELMT>
  inline FTYPE StateZTraj<FTYPE, NBELMT>::distTo1stError( FTYPE z, FTYPE tolerance, int /*pathDirection*/ ) const {
    // look only at x (because it curves most)
    FTYPE deriv = poly2ndderiveval( z - m_z, m_cx );
#ifdef __INTEL_COMPILER             // Disable ICC remark
#  pragma warning( disable : 1572 ) // Floating-point equality and inequality comparisons are unreliable
#  pragma warning( push )
#endif
    deriv = upderiv<FTYPE, NBELMT>( deriv, tolerance );
#ifdef __INTEL_COMPILER // End disable ICC remark
#  pragma warning( pop )
#endif
    return deriv;
  }

  template <typename FTYPE, unsigned int NBELMT>
  inline FTYPE StateZTraj<FTYPE, NBELMT>::distTo2ndError( FTYPE /*z*/, FTYPE /*tolerance*/,
                                                          int /*pathDirection*/ ) const {
    return 10 * Gaudi::Units::km;
  }

  template <typename FTYPE, unsigned int NBELMT>
  inline void StateZTraj<FTYPE, NBELMT>::expansion( FTYPE z, typename StateZTraj<FTYPE, NBELMT>::Point& p,
                                                    typename StateZTraj<FTYPE, NBELMT>::Vector& dp,
                                                    typename StateZTraj<FTYPE, NBELMT>::Vector& ddp ) const {
    p.SetXYZ( x( z ), y( z ), z );
    dp.SetXYZ( tx( z ), ty( z ), 1 );
    ddp.SetXYZ( omegax( z ), omegay( z ), 0. );
  }

  template <typename FTYPE, unsigned int NBELMT>
  inline FTYPE StateZTraj<FTYPE, NBELMT>::muEstimate( const StateZTraj<FTYPE, NBELMT>::Point& p ) const {
    StateZTraj<FTYPE, NBELMT>::Vector dir  = direction( p.z() );
    StateZTraj<FTYPE, NBELMT>::Vector dx   = p - position( p.z() );
    FTYPE                             dir2 = dir.mag2();
    FTYPE                             det  = dir2 - curvature( p.z() ).Dot( dx );

    if constexpr ( !std::is_floating_point<FTYPE>::value ) {
      det = select( det > (FTYPE)0., dir2, det );
    } else {
      det = det > (FTYPE)0. ? dir2 : det;
    }
    return p.z() + dx.Dot( dir ) / det;
  }

  template <typename FTYPE, unsigned int NBELMT>
  inline typename StateZTraj<FTYPE, NBELMT>::Parameters StateZTraj<FTYPE, NBELMT>::parameters() const {
    return {m_cx[0], m_cy[0], m_cx[1], m_cy[1], m_qOverP};
  }

  template <typename FTYPE, unsigned int NBELMT>
  inline StateZTraj<FTYPE, NBELMT>& StateZTraj<FTYPE, NBELMT>::operator+=( const Parameters& /*delta*/ ) {
    // to implement this we need the full b-field.
    std::cerr << __FUNCTION__ << " not yet implemented." << std::endl;
    return *this;
  }

} // namespace LHCb

#endif /// TRACKFITEVENT_STATETRAJ_H
