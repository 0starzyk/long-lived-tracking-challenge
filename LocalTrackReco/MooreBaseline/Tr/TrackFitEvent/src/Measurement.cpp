/*****************************************************************************\
* (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "Event/Measurement.h"
#include "FTDet/DeFTMat.h"
#include "Kernel/ROOTExtensions.h"
#include "MuonDet/DeMuonChamber.h"
#include "UTDet/DeUTSector.h"
#include "VPDet/DeVPSensor.h"

Gaudi::XYZPoint LHCb::Measurement::toLocal( const Gaudi::XYZPoint& globalPoint ) const {
  return visit( [&globalPoint]( const UT& ut ) { return ut.sec().toLocal( globalPoint ); },
                [&globalPoint]( const FT& ft ) { return ft.mat().toLocal( globalPoint ); },
                [&globalPoint]( const Muon& muon ) { return muon.chamber().toLocal( globalPoint ); },
                [&globalPoint]( const VP& vp ) { return vp.sensor().toLocal( globalPoint ); },
                [&globalPoint]( const VP2D& vp ) { return vp.sensor().toLocal( globalPoint ); } );
}

Gaudi::XYZVector LHCb::Measurement::toLocal( const Gaudi::XYZVector& globalPoint ) const {
  return visit( [&globalPoint]( const UT& ut ) { return ut.sec().toLocal( globalPoint ); },
                [&globalPoint]( const FT& ft ) { return ft.mat().toLocal( globalPoint ); },
                [&globalPoint]( const Muon& muon ) { return muon.chamber().toLocal( globalPoint ); },
                [&globalPoint]( const VP& vp ) { return vp.sensor().toLocal( globalPoint ); },
                [&globalPoint]( const VP2D& vp ) { return vp.sensor().toLocal( globalPoint ); } );
}

Gaudi::XYZPoint LHCb::Measurement::toGlobal( const Gaudi::XYZPoint& globalPoint ) const {
  return visit( [&globalPoint]( const UT& ut ) { return ut.sec().toGlobal( globalPoint ); },
                [&globalPoint]( const FT& ft ) { return ft.mat().toGlobal( globalPoint ); },
                [&globalPoint]( const Muon& muon ) { return muon.chamber().toGlobal( globalPoint ); },
                [&globalPoint]( const VP& vp ) { return vp.sensor().toGlobal( globalPoint ); },
                [&globalPoint]( const VP2D& vp ) { return vp.sensor().toGlobal( globalPoint ); } );
}

Gaudi::XYZVector LHCb::Measurement::toGlobal( const Gaudi::XYZVector& globalPoint ) const {
  return visit( [&globalPoint]( const UT& ut ) { return ut.sec().toGlobal( globalPoint ); },
                [&globalPoint]( const FT& ft ) { return ft.mat().toGlobal( globalPoint ); },
                [&globalPoint]( const Muon& muon ) { return muon.chamber().toGlobal( globalPoint ); },
                [&globalPoint]( const VP& vp ) { return vp.sensor().toGlobal( globalPoint ); },
                [&globalPoint]( const VP2D& vp ) { return vp.sensor().toGlobal( globalPoint ); } );
}

std::string LHCb::Measurement::name() const {
  return visit( []( const UT& ut ) -> std::string { return ut.sec().name(); },
                []( const FT& ft ) -> std::string { return ft.mat().name(); },
                []( const Muon& muon ) -> std::string { return muon.chamber().name(); },
                []( const VP& vp ) -> std::string { return vp.sensor().name(); },
                []( const VP2D& vp ) -> std::string { return vp.sensor().name(); } );
}

bool LHCb::Measurement::isSameDetectorElement( const Measurement& other ) const {
  return std::visit(
      []( const auto& lhs, const auto& rhs ) {
        if constexpr ( std::is_same_v<decltype( lhs ), decltype( rhs )> ) {
          return lhs.detElem == rhs.detElem;
        } else
          return false;
      },
      m_sub, other.m_sub );
}

namespace {
  // trivial helpers to make code clearer...
  typedef Gaudi::Matrix1x3 DualVector;

  inline DualVector dual( const Gaudi::XYZVector& v ) {
    DualVector d;
    v.GetCoordinates( d.Array() );
    return d;
  }
} // namespace

LHCb::MinimizeResult LHCb::minimize( const LHCb::Measurement& m, const LHCb::StateZTraj<double>& refTraj,
                                     double zState ) {

  return m.visit(
      [&]( LHCb::Measurement::VP2D const& m ) -> MinimizeResult {
        // 2D minimize for VP
        auto                             unitPoca = refTraj.position( m.trajectory.Z() ) - m.trajectory;
        constexpr std::array<double, 10> d        = {1, 0, 0, 0, 0, 0, -1, 0, 0, 0};
        auto                             H        = Gaudi::Matrix2x5{d.begin(), d.end()};
        return LHCb::Minimize2DResult{zState, 0., unitPoca.Dot( unitPoca ), std::move( H ), std::move( unitPoca )};
      },
      [&]( auto const& m ) -> MinimizeResult {
        // Determine the two points on the measurement and reference trajectory
        // between which we minimize the distance between the trajectories
        // ---> Yields Distance of Closest Approach (DOCA)
        //
        // This function is mainly used from the track fit.
        // It provides extra output and handles potential local curvature of refTraj
        // If you only need the doca or the location of the two points for two straight trajectories
        // have a look what's defined in LHCb/Kernel/LHCbMath/GeomFun.h
        //
        // The curvature of the reference trajectory is evaluted at zState and treated as
        // constant in the proximity of zState. This assumption only holds,
        // if the final point on the reference trajectory is in fact close to zState.
        // The final point is usually so close that the curvature term doesn't matter much
        // which is why most (all?) configurations call this function with a refTraj that
        // doesn't account for the local B-field vector
        //
        // The distance of closest approach (doca) vector is defined as (p_r - p_m)
        // where p_r and p_m are points on the reference and measurement trajectory.
        // The below answers: if I start at a point p on each trajectory,
        // how far along their respective direction d do I have to go,
        // to reach the points p_r or p_m respectively.
        // The result to this called mu[0] and mu[1] below
        // This then defines the two points as:
        // p_r = p_r_zstate - mu[0] * d_r (+ potentially terms for curvature)
        // p_m = p_m_0 - mu[1] * d_m

        // this is the distance between our starting points p_r_zstate and p_m_0
        auto const dist = refTraj.position( zState ) - m.trajectory.position( 0 );
        auto const d_m  = m.trajectory.direction( 0 );
        auto const d_r  = refTraj.direction( zState );
        // extra term corresponding to curvature if the reference traj includes treatment of B-field
        // normally we don't do that, thus c_r is just a vector of 0s
        auto const c_r = refTraj.curvature( zState );
        // now we solve the linear system for mu
        auto const mat    = std::array{d_m.mag2(), -d_r.Dot( d_m ), d_r.mag2() - dist.Dot( c_r )};
        auto       mu     = std::array{-dist.Dot( d_m ), dist.Dot( d_r )};
        auto const decomp = ROOT::Math::CholeskyDecomp<double, 2>{mat.data()};
        if ( !decomp.Solve( mu ) ) throw std::runtime_error( "singular matrix" );

        // we move to the correct point on the reference trajectory
        auto const p_r = zState - mu[1];
        // we started at 0 so m_r = 0 - mu[0]
        auto const m_r = -mu[0];

        // distance vector between the two points of closest approach
        auto const doca_vec = refTraj.position( p_r ) - m.trajectory.position( m_r );

        // vector onto which we project doca_vec to obtain a signed real value for the doca.
        auto const unitPoca = d_m.Cross( refTraj.direction( p_r ) ).Unit();

        // compute the projection matrix from parameter space onto the (signed!) unit
        auto H = eval( dual( unitPoca ) * refTraj.derivative( p_r ) );

        return LHCb::Minimize1DResult{p_r, m_r, unitPoca.Dot( doca_vec ), std::move( H ), std::move( unitPoca )};
      } );
}
