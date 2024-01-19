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
#include "DetDesc/Material.h"
#include "Event/State.h"
#include "Event/TrackParameters.h"
#include "GaudiAlg/GaudiTool.h"
#include "TrackInterfaces/IStateCorrectionTool.h"
#include "vdt/log.h"

using namespace Gaudi::Units;

//-----------------------------------------------------------------------------
// Implementation file for class : StateThinMSCorrectionTool
//
// 2006-08-21 : Eduardo Rodrigues
//-----------------------------------------------------------------------------

/** @class StateThinMSCorrectionTool StateThinMSCorrectionTool.h
 *
 *  This state correction tool applies a multiple scattering correction
 *  in the approximation of a thin scatter
 *
 *  @author Eduardo Rodrigues
 *  @date   2006-08-21
 *  (Original code taken from the master extrapolator)
 */
class StateThinMSCorrectionTool final : public extends<GaudiTool, IStateCorrectionTool> {
public:
  /// Standard constructor
  using extends::extends;

  /// Correct a State for multiple scattering in the approximation of a thin scatter
  void correctState( LHCb::State& state, const MaterialPtr material, std::any& cache, double wallThickness,
                     bool upstream, double ) const override;

private:
  // Job options
  Gaudi::Property<double> m_msff2{this, "MSFudgeFactor2", 1.0}; ///< fudge factor for multiple scattering errors
};

// Declaration of the Tool Factory
DECLARE_COMPONENT( StateThinMSCorrectionTool )

//=============================================================================
// Correct a State for multiple scattering
// in the approximation of a thin scatter
//=============================================================================
void StateThinMSCorrectionTool::correctState( LHCb::State& state, const MaterialPtr material, std::any& /*cache*/,
                                              double wallThickness, bool, double ) const {
#ifdef USE_DD4HEP
  const double t = wallThickness / material.radiationLength();
#else
  const double t = wallThickness / material->radiationLength();
#endif
  const double norm2      = 1.0 + std::pow( state.tx(), 2 ) + std::pow( state.ty(), 2 );
  double       scatLength = 0.;
  if ( t > TrackParameters::lowTolerance ) {
    const double radThick = sqrt( norm2 ) * t;
    scatLength = radThick * std::pow( TrackParameters::moliereFactor * ( 1. + 0.038 * vdt::fast_log( radThick ) ), 2 );
  }

  // protect against zero momentum
  const double p = std::max( state.p(), 1.0 * MeV );

  const double norm2cnoise = norm2 * m_msff2 * scatLength / std::pow( p, 2 );

  // update covariance matrix C = C + Q
  // multiple scattering covariance matrix Q only has three elements...
  Gaudi::TrackSymMatrix& tC = state.covariance();
  tC( 2, 2 ) += norm2cnoise * ( 1. + std::pow( state.tx(), 2 ) );
  tC( 3, 3 ) += norm2cnoise * ( 1. + std::pow( state.ty(), 2 ) );
  tC( 3, 2 ) += norm2cnoise * state.tx() * state.ty();
}

//=============================================================================
