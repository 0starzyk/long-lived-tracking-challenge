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
#include "GaudiAlg/GaudiTool.h"
#include "TrackInterfaces/IStateCorrectionTool.h"
#include <cmath>

//-----------------------------------------------------------------------------
// Implementation file for class : StateElectronEnergyCorrectionTool
//
// 2006-08-18 : Eduardo Rodrigues
//-----------------------------------------------------------------------------

/** @class StateElectronEnergyCorrectionTool
 *
 *  This state correction tool applies a dE/dx energy loss correction
 *  to electrons
 *
 *  @author Eduardo Rodrigues
 *  @date   2006-08-18
 *  (Original code taken from the master extrapolator)
 */
class StateElectronEnergyCorrectionTool : public extends<GaudiTool, IStateCorrectionTool> {
public:
  /// Standard constructor
  using extends::extends;
  /// Correct a State for electron dE/dx energy losses
  void correctState( LHCb::State& state, const MaterialPtr material, std::any& cache, double wallThickness,
                     bool upstream, double ) const override;

private:
  // Job options
  Gaudi::Property<double> m_maxRadLength{this, "MaximumRadLength", 10.}; ///< maximum radiation length (in mm)
};

// Declaration of the Tool Factory
DECLARE_COMPONENT( StateElectronEnergyCorrectionTool )

//=============================================================================
// Correct a State for electron dE/dx energy losses
//=============================================================================
void StateElectronEnergyCorrectionTool::correctState( LHCb::State& state, const MaterialPtr material,
                                                      std::any& /*cache*/, double wallThickness, bool upstream,
                                                      double ) const {
  // hard energy loss for electrons
  double t =
#ifdef USE_DD4HEP
      wallThickness / material.radiationLength() * sqrt( 1. + std::pow( state.tx(), 2 ) + std::pow( state.ty(), 2 ) );
#else
      wallThickness / material->radiationLength() * sqrt( 1. + std::pow( state.tx(), 2 ) + std::pow( state.ty(), 2 ) );
#endif
  if ( !upstream ) t *= -1.;

  // protect against t too big
  if ( fabs( t ) > m_maxRadLength ) t = std::copysign( 1.0, t ) * m_maxRadLength;

  // apply correction
  Gaudi::TrackVector&    tX = state.stateVector();
  Gaudi::TrackSymMatrix& tC = state.covariance();

  tC( 4, 4 ) += std::pow( tX[4], 2 ) * ( exp( -t * log( 3.0 ) / log( 2.0 ) ) - exp( -2.0 * t ) );
  tX[4] *= exp( -t );
}

//=============================================================================
