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
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Gaudi::Units;

namespace {
  template <typename T>
  inline decltype( auto ) pow_2( T x ) noexcept {
    return x * x;
  }
} // namespace
//-----------------------------------------------------------------------------
// Implementation file for class : StateThickMSCorrectionTool
//
// 2006-08-21 : Eduardo Rodrigues
//-----------------------------------------------------------------------------

/** @class StateThickMSCorrectionTool StateThickMSCorrectionTool.h
 *
 *  This state correction tool applies a multiple scattering correction
 *  for a thick scatter.
 *  As scattering formula is not linear, add possibility to omit log term
 *  using 'RossiAndGreisen', which should lead to better pulls.
 *
 * Parameters:
 *  - MSFudgeFactor2: Scale factor (squared) multiplied with the (13.6MeV)^2 in MS formula
 *  - UseRossiAndGreisen: Don't use the log-term in MS formula (this is recommended from 2015 on, as the log-term breaks
 * the linearity of many scatterers after each other)
 *  - RossiAndGreisenFact2: Constant (squared) in MS formula to replace the (13.6MeV)^2. Only used if larger than 0.
 * Default is '-1e42' and then (13.6MeV)^2 is used.
 *
 *  @author Eduardo Rodrigues
 *  @date   2006-08-21
 *  (Original code taken from the master extrapolator)
 *
 *  @author Miriam Hess, Michel De Cian
 *  @date 2015-07-10
 *
 */
class StateThickMSCorrectionTool : public extends<GaudiTool, IStateCorrectionTool> {
public:
  /// Standard constructor
  using extends::extends;

  StatusCode initialize() override;

  /// Correct a State for multiple scattering in the case of a thick scatter
  void correctState( LHCb::State& state, const MaterialPtr material, std::any& cache, double wallThickness = 0,
                     bool upstream = true, double pmass = 0 ) const override;

private:
  double m_msff2MoliereFact2 =
      std::numeric_limits<double>::quiet_NaN(); ///< auxiliary variable, precomputed; break things good and hard if
                                                ///< initialize isn't called

  // Job options
  Gaudi::Property<double> m_msff2{this, "MSFudgeFactor2", 1.0}; ///< fudge factor for multiple scattering errors
  Gaudi::Property<double> m_msff2RossiAndGreisenFact2{this, "RossiAndGreisenFact2",
                                                      -1e42}; ///< factor to use when log term is omitted
  Gaudi::Property<bool>   m_useRossiAndGreisen{this, "UseRossiAndGreisen", false}; ///< if true, not log term is used
};

// Declaration of the Tool Factory
DECLARE_COMPONENT( StateThickMSCorrectionTool )

//=============================================================================
// Initialize variables
//=============================================================================
StatusCode StateThickMSCorrectionTool::initialize() {
  StatusCode sc = GaudiTool::initialize();
  if ( sc.isFailure() ) return sc;
  // move computation out of the hot code path
  m_msff2MoliereFact2 = m_msff2 * pow_2( TrackParameters::moliereFactor );

  if ( m_msff2RossiAndGreisenFact2 < 0 ) m_msff2RossiAndGreisenFact2 = m_msff2MoliereFact2;
  if ( msgLevel( MSG::DEBUG ) ) {
    if ( m_useRossiAndGreisen ) {
      debug() << "Using Rossi & Greisen Factor squared of " << m_msff2RossiAndGreisenFact2 << endmsg;
    } else {
      debug() << "Using log term with Moliere factor squared of " << m_msff2MoliereFact2 << endmsg;
    }
  }

  return sc;
}
//=============================================================================
// Correct a State for multiple scattering in the case of a thick scatter
//=============================================================================
void StateThickMSCorrectionTool::correctState( LHCb::State& state, const MaterialPtr material, std::any& /*cache*/,
                                               double wallThickness, bool upstream, double ) const {
#ifdef USE_DD4HEP
  const auto t = wallThickness / material.radiationLength();
#else
  const auto t = wallThickness / material->radiationLength();
#endif
  // if t is below tolerance, all corrections end up zero anyway, so we can
  // stop early
  if ( t <= TrackParameters::lowTolerance ) return;

  const auto& stv   = state.stateVector(); // x, y, tx, ty, q/p
  const auto  norm2 = 1 + pow_2( stv[2] ) + pow_2( stv[3] );
  // protect against zero momentum
  static constexpr decltype( stv[4] ) iMeV           = 1. / MeV;
  const auto                          norm2cnoisetmp = norm2 * pow_2( std::min( std::abs( stv[4] ), iMeV ) );

  const auto radThick = t * std::sqrt( norm2 );
  // in a normal tracking run, for around 95% of cases, radThick is in the
  // interval [1e-7, 1e-1] - that's the region where the log turns nasty
  // because it decreases so rapidly, and it's not easily possible to
  // approximate it in a quick and dirty way
  //
  // be FMA friendly
  auto norm2cnoise = norm2cnoisetmp * radThick;

  // -- first one omits log-term, which 'solves' the problem of having a non-linear behaviour when going through many
  // layers
  if ( m_useRossiAndGreisen ) {
    norm2cnoise *= m_msff2RossiAndGreisenFact2;
  } else {
    norm2cnoise *= m_msff2MoliereFact2 * pow_2( .038 * vdt::fast_log( radThick ) + 1 );
  }

  const auto                                 covTxTx          = norm2cnoise * ( 1 + pow_2( stv[2] ) );
  const auto                                 covTyTy          = norm2cnoise * ( 1 + pow_2( stv[3] ) );
  const auto                                 covTxTy          = norm2cnoise * stv[2] * stv[3];
  const auto                                 wallThicknessD   = wallThickness * ( upstream ? -.5 : .5 );
  static constexpr decltype( wallThickness ) thirds           = 1. / 3.;
  const auto                                 wallThickness2_3 = pow_2( wallThickness ) * thirds;

  // update covariance matrix C = C + Q
  auto& cov = state.covariance();
  cov( 0, 0 ) += covTxTx * wallThickness2_3,     // don't care about the order or
      cov( 1, 0 ) += covTxTy * wallThickness2_3, // update so use commata here,
      cov( 1, 1 ) += covTyTy * wallThickness2_3, // and let compiler decide
      cov( 2, 0 ) += covTxTx * wallThicknessD, cov( 2, 1 ) += covTxTy * wallThicknessD, cov( 2, 2 ) += covTxTx,
      cov( 3, 0 ) += covTxTy * wallThicknessD, cov( 3, 1 ) += covTyTy * wallThicknessD, cov( 3, 2 ) += covTxTy,
      cov( 3, 3 ) += covTyTy;
}
