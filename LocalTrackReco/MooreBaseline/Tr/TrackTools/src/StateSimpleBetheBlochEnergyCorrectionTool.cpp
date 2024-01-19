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
#include "GaudiKernel/SystemOfUnits.h"
#include "LHCbMath/LHCbMath.h"
#include "TrackInterfaces/IStateCorrectionTool.h"
#include <cmath>

using namespace Gaudi::Units;

//-----------------------------------------------------------------------------
// Implementation file for class : StateSimpleBetheBlochEnergyCorrectionTool
//
// 2006-08-18 : Eduardo Rodrigues
//-----------------------------------------------------------------------------

/** @class StateSimpleBetheBlochEnergyCorrectionTool
 *
 *  This state correction tool applies a dE/dx energy loss correction
 *  with a simplified version of the Bethe-Bloch equation.
 *
 *  @author Eduardo Rodrigues
 *  @date   2006-08-18
 *  (Original code taken from the master extrapolator)
 */
class StateSimpleBetheBlochEnergyCorrectionTool : public extends<GaudiTool, IStateCorrectionTool> {
public:
  /// Standard constructor
  using extends::extends;
  /// Correct a State for dE/dx energy losses with a simplified Bethe-Bloch equiaton
  void correctState( LHCb::State& state, const MaterialPtr material, std::any& cache, double wallThickness,
                     bool upstream, double ) const override;

private:
  // Job options
  Gaudi::Property<double> m_energyLossCorr{this, "EnergyLossFactor",
                                           354.1 * Gaudi::Units::MeV* Gaudi::Units::mm2 /
                                               Gaudi::Units::mole}; ///< tunable energy loss correction
  Gaudi::Property<double> m_maxEnergyLoss{this, "MaximumEnergyLoss",
                                          100. * Gaudi::Units::MeV}; ///< maximum energy loss in dE/dx correction
  Gaudi::Property<double> m_minMomentumAfterEnergyCorr{
      this, "MinMomentumAfterEnergyCorr", 10. * Gaudi::Units::MeV}; ///< minimum momentum after dE/dx correction
  Gaudi::Property<double> m_sqrtEError{this, "EnergyLossError",
                                       1.7 * sqrt( Gaudi::Units::MeV )}; // proportional to sqrt(E); ///<
                                                                         // sigma(dE)/sqrt(dE)
  Gaudi::Property<bool> m_useEnergyLossError{this, "UseEnergyLossError",
                                             false}; ///< flag to turn on using error on energy loss
};
// Declaration of the Tool Factory
DECLARE_COMPONENT( StateSimpleBetheBlochEnergyCorrectionTool )

//=============================================================================
// Correct a State for dE/dx energy losses
//=============================================================================
void StateSimpleBetheBlochEnergyCorrectionTool::correctState( LHCb::State& state, const MaterialPtr material,
                                                              std::any& /*cache*/, double wallThickness, bool upstream,
                                                              double ) const {
  double bbLoss = wallThickness * sqrt( 1. + std::pow( state.tx(), 2 ) + std::pow( state.ty(), 2 ) ) *
#ifdef USE_DD4HEP
                  m_energyLossCorr * material.Z() * material.density() / material.A();
#else
                  m_energyLossCorr * material->Z() * material->density() / material->A();
#endif
  bbLoss = std::min( m_maxEnergyLoss.value(), bbLoss );
  if ( !upstream ) bbLoss *= -1.;

  // apply correction - note for now only correct the state vector
  Gaudi::TrackVector& tX = state.stateVector();

  //  double minMomentumForEnergyCorrection = 10*Gaudi::Units::MeV;

  double qOverP        = 0.0;
  tX[4] < 0.0 ? qOverP = std::min( tX[4], -LHCb::Math::lowTolerance )
              : qOverP = std::max( tX[4], LHCb::Math::lowTolerance );

  double newP         = 0.0;
  qOverP < 0.0 ? newP = std::min( 1.0 / qOverP - bbLoss, -m_minMomentumAfterEnergyCorr.value() )
               : newP = std::max( 1.0 / qOverP + bbLoss, m_minMomentumAfterEnergyCorr.value() );

  tX[4] = 1.0 / newP;

  // correction on cov
  if ( m_useEnergyLossError && m_sqrtEError > 0 ) {
    // error on dE is proportional to the sqrt of dE:
    // double sigmadE = m_sqrtEError * std::sqrt(std::abs(bbLoss))
    double                 err2 = m_sqrtEError * m_sqrtEError * std::abs( bbLoss ) * tX[4] * tX[4];
    Gaudi::TrackSymMatrix& cov  = state.covariance();
    cov( 4, 4 ) += err2;
  }
}

//=============================================================================
