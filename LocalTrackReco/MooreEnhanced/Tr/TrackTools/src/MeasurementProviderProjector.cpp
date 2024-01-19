/*****************************************************************************\
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
/** @class MeasurementProviderProjector MeasurementProviderProjector.cpp
 *
 * Implementation of templated MeasurementProviderProjector tool
 * see interface header for description
 *
 *  @author A. Usachov
 *  @date   15/10/2019
 */

#include "MeasurementProviderProjector.h"

using namespace ROOT::Math;
using ROOT::Math::SMatrix;
using ROOT::Math::SVector;

MeasurementProviderProjector::ProjectResult
MeasurementProviderProjector::internal_project( const LHCb::StateVector& statevector,
                                                const LHCb::Measurement& meas ) const {
  // Project onto the reference. First create the StateTraj with or without BField information.
  Gaudi::XYZVector bfield = m_useBField ? m_magnet.get().fieldVector( statevector.position() ) : Gaudi::XYZVector{};
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
std::tuple<StatusCode, MeasurementProviderProjector::ProjectResult>
MeasurementProviderProjector::project( const LHCb::State& state, const LHCb::Measurement& meas ) const {
  try {
    // Project onto the reference (prevent the virtual function call)
    auto projectResult = internal_project( LHCb::StateVector( state.stateVector(), state.z() ), meas );
    return {StatusCode{StatusCode::SUCCESS}, projectResult};
    // Calculate the error on the residual
    // m_errResidual = sqrt( m_errMeasure*m_errMeasure + Similarity( m_H, state.covariance() )(0,0) );
  } catch ( StatusCode sc ) { return {sc, MeasurementProviderProjector::Project1DResult{}}; }
}

//-----------------------------------------------------------------------------
/// Project the state vector in this fitnode and update projection matrix and reference residual
//-----------------------------------------------------------------------------
StatusCode MeasurementProviderProjector::projectReference( LHCb::FitNode& node ) const {
  StatusCode sc = StatusCode::FAILURE;
  if ( node.hasMeasurement() ) {
    auto projectResult = internal_project( node.refVector(), node.measurement() );

    node.visit2(
        projectResult,
        [&]( LHCb::FitNode::DimInfos<LHCb::Enum::nDim::Type::one>& n, Project1DResult result ) {
          n.updateProjection( node, result.H, result.residual, result.errMeasure );
          node.setPocaVector( std::move( result.unitPocaVector ) );
          node.setDoca( result.doca );
          sc = StatusCode::SUCCESS;
        },
        [&]( LHCb::FitNode::DimInfos<LHCb::Enum::nDim::Type::two>& n, Project2DResult result ) {
          n.updateProjection( node, result.H, result.residual, result.errMeasure );
          node.setPocaVector( std::move( result.unitPocaVector ) );
          node.setDoca( result.doca );
          sc = StatusCode::SUCCESS;
        },
        [&]( ... ) { throw std::logic_error( "invalid combination of 1D and 2D variants" ); } );
  }
  return sc;
}
