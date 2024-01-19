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
#include "TrackProjector.h"
#include "Event/FitNode.h"
#include "Event/Measurement.h"
#include "Event/State.h"
#include "Event/StateVector.h"
#include "Event/StateZTraj.h"
#include "Kernel/ITrajPoca.h"

using namespace ROOT::Math;
using ROOT::Math::SMatrix;
using ROOT::Math::SVector;

DECLARE_COMPONENT( TrackProjector )

//-----------------------------------------------------------------------------
/// Initialize
//-----------------------------------------------------------------------------
StatusCode TrackProjector::initialize() {
  return GaudiTool::initialize().andThen( [&] { m_poca = tool<ITrajPoca>( "TrajPoca" ); } );
}

// trivial helpers to make code clearer...
namespace {
  typedef Gaudi::Matrix1x3 DualVector;

  DualVector dual( const Gaudi::XYZVector& v ) {
    DualVector d;
    v.GetCoordinates( d.Array() );
    return d;
  }

} // namespace

//-----------------------------------------------------------------------------
// internal project method, doing the actual work
//-----------------------------------------------------------------------------
TrackProjector::ProjectResult TrackProjector::internal_project( const LHCb::StateVector& statevector,
                                                                const LHCb::Measurement& meas ) const {
  // Project onto the reference. First create the StateTraj with or without BField information.
  Gaudi::XYZVector bfield( 0, 0, 0 );
  if ( m_useBField ) {
    auto& ctx = Gaudi::Hive::currentContext();
    bfield    = m_magneticfield.get( getConditionContext( ctx ) ).fieldVector( statevector.position() );
  }
  const LHCb::StateZTraj<double> refTraj( statevector, bfield );

  // minimize in measurement returns a variant
  LHCb::MinimizeResult minimizeResult = minimize( meas, refTraj, statevector.z() );

  return meas.visit(
      minimizeResult,
      [&]( const auto& m, LHCb::Minimize1DResult& mr ) -> ProjectResult {
        auto errMeasure = ROOT::Math::SVector<double, 1>( m.errMeasure );
        auto residual   = ROOT::Math::SVector<double, 1>( -mr.doca );
        return Project1DResult{mr.sMeas,   mr.doca,           residual,
                               errMeasure, std::move( mr.H ), std::move( mr.unitPocaVector )};
      },
      [&]( const LHCb::Measurement::VP2D& vp2d, LHCb::Minimize2DResult& mr ) -> ProjectResult {
        auto dist     = refTraj.position( mr.zState ) - vp2d.trajectory;
        auto residual = ROOT::Math::SVector<double, 2>( -dist.X(), dist.Y() ); // TODO minus sign? wrong def elsewhere?
        return Project2DResult{mr.sMeas,        mr.doca,           residual,
                               vp2d.errMeasure, std::move( mr.H ), std::move( mr.unitPocaVector )};
      },
      [&]( const LHCb::Measurement::VP2D&, LHCb::Minimize1DResult& ) -> ProjectResult {
        throw std::logic_error( "invalid combination of 1D and 2D variants" );
      },
      [&]( ... ) -> ProjectResult { throw std::logic_error( "invalid combination of 1D and 2D variants" ); } );
}

//-----------------------------------------------------------------------------
/// Project a state onto a measurement
//-----------------------------------------------------------------------------
std::tuple<StatusCode, TrackProjector::ProjectResult> TrackProjector::project( const LHCb::State&       state,
                                                                               const LHCb::Measurement& meas ) const {
  try {
    // Project onto the reference (prevent the virtual function call)
    auto projectResult = internal_project( LHCb::StateVector( state.stateVector(), state.z() ), meas );
    return {StatusCode{StatusCode::SUCCESS}, projectResult};
    // Calculate the error on the residual
    // m_errResidual = sqrt( m_errMeasure*m_errMeasure + Similarity( m_H, state.covariance() )(0,0) );
  } catch ( StatusCode sc ) { return {sc, ITrackProjector::Project1DResult{}}; }
}

//-----------------------------------------------------------------------------
/// Project the state vector in this fitnode and update projection matrix and reference residual
//-----------------------------------------------------------------------------
StatusCode TrackProjector::projectReference( LHCb::FitNode& node ) const {
  StatusCode sc = StatusCode::FAILURE;
  if ( node.hasMeasurement() ) {
    auto projectResult = internal_project( node.refVector(), node.measurement() );

    node.visit2(
        projectResult,
        [&]( LHCb::FitNode::DimInfos<LHCb::Enum::nDim::Type::one>& n, ITrackProjector::Project1DResult result ) {
          n.updateProjection( node, result.H, result.residual, result.errMeasure );
          node.setPocaVector( std::move( result.unitPocaVector ) );
          node.setDoca( result.doca );
          sc = StatusCode::SUCCESS;
        },
        [&]( LHCb::FitNode::DimInfos<LHCb::Enum::nDim::Type::two>& n, ITrackProjector::Project2DResult result ) {
          n.updateProjection( node, result.H, result.residual, result.errMeasure );
          node.setPocaVector( std::move( result.unitPocaVector ) );
          node.setDoca( result.doca );
          sc = StatusCode::SUCCESS;
        },
        [&]( ... ) { throw std::logic_error( "invalid combination of 1D and 2D variants" ); } );
  }
  return sc;
}

//-----------------------------------------------------------------------------
/// Derivatives wrt.the measurement's alignment...
//-----------------------------------------------------------------------------
TrackProjector::Derivatives TrackProjector::alignmentDerivatives( const LHCb::StateVector& statevector,
                                                                  const LHCb::Measurement& meas,
                                                                  const Gaudi::XYZPoint&   pivot ) const {
  auto       resultvar = internal_project( statevector, meas );
  DualVector unit = std::visit( [&]( const auto r ) -> DualVector { return dual( r.unitPocaVector ); }, resultvar );

  // Calculate the derivative of the poca on measTraj to alignment parameters.
  // Only non-zero elements:
  // TODO: implement alignment derivatives for 2D measurements.
  LHCb::LineTraj<double> trajectory =
      meas.visit( [&]( const auto& m ) { return m.trajectory; },
                  [&]( const LHCb::Measurement::VP2D& ) {
                    throw std::runtime_error( "Alignment derivatives not implemented for 2D measurements." );
                    LHCb::LineTraj<double> dummy;
                    return dummy;
                  } );

  auto sMeas = std::visit(
      [&]( const auto r ) -> auto { return r.sMeas; }, resultvar );
  Gaudi::XYZVector                  arm = trajectory.position( sMeas ) - pivot;
  ROOT::Math::SMatrix<double, 3, 6> dPosdAlpha;
  // Derivative to translation
  dPosdAlpha( 0, 0 ) = dPosdAlpha( 1, 1 ) = dPosdAlpha( 2, 2 ) = 1;
  // Derivative to rotation around x-axis
  dPosdAlpha( 1, 3 ) = -arm.z();
  dPosdAlpha( 2, 3 ) = arm.y();
  // Derivative to rotation around y-axis
  dPosdAlpha( 0, 4 ) = arm.z();
  dPosdAlpha( 2, 4 ) = -arm.x();
  // Derivative to rotation around z-axis
  dPosdAlpha( 0, 5 ) = -arm.y();
  dPosdAlpha( 1, 5 ) = arm.x();

  return unit * dPosdAlpha;
  // compute the projection matrix from parameter space onto the (signed!) unit
  // return unit*AlignTraj( meas.trajectory(), pivot ).derivative( m_sMeas );
}
