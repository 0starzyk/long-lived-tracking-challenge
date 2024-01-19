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
#include <Event/State.h>
#include <Event/TrackParameters.h>
#include <GaudiAlg/GaudiTool.h>
#include <GaudiKernel/SystemOfUnits.h>
#include <Kernel/ILHCbMagnetSvc.h>
#include <Kernel/STLExtensions.h>
#include <TrackInterfaces/ITrackMomentumEstimate.h>
#include <array>
#include <cmath>
#include <numeric>

//-----------------------------------------------------------------------------
// Implementation file for class : FastMomentumEstimate
//
// 2007-10-30 : S. Hansmann-Menzemer
//-----------------------------------------------------------------------------

/** @class FastMomentumEstimate FastMomentumEstimate.h
 *
 *  Calculate momentum of a seed track or a match track using the states and a parametrisation
 *  The parameters for VeloTCubic can be obtained using the "MatchFitParams" algorithm.
 *
 *  Parameters:
 *
 *  - ParamsTCubic: Coefficients for momentum calculation using a state in the T stations when using a cubic fit.
 *  - ParamsTParabola: Coefficients for momentum calculation using a state in the T stations when using a parbolic fit.
 *  - ParamsVeloTCubic: Coefficients for momentum calculation using a state in the Velo and the T stations when using a
 * cubic fit.
 *  - ParamsVeloTParabola: Coefficients for momentum calculation using a state in the Velo and the T stations when using
 * a parbolic fit.
 *  - TResolution: Resolution for q/p, given as: sigma(q/p) = TResolution * (q/p)
 *  - VeloPlusTResolution: Resolution for q/p, given as: sigma(q/p) = VeloPlusTResolution * (q/p)
 *
 *  Note: When giving the momentum parameters / coefficients in the options file, all 4 vectors have to be given.
 *
 *  @author Stephanie Hansmann-Menzemer
 *  @date   2007-10-30
 */
class FastMomentumEstimate : public extends<GaudiTool, ITrackMomentumEstimate> {
public:
  /// Standard constructor
  using extends::extends;

  StatusCode initialize() override;

  /// Estimate the momentum P of a State in T at ZAtMidT
  StatusCode calculate( const DeMagnet& magnet, const LHCb::State* tState, double& qOverP, double& sigmaQOverP,
                        bool tCubicFit ) const override;

  /// Estimate the momentum P of a velo State and a State in T (for matching)
  StatusCode calculate( const DeMagnet& magnet, const LHCb::State* veloState, const LHCb::State* tState, double& qOverP,
                        double& sigmaQOverP, bool tCubicFit ) const override;

  StatusCode calculate( const DeMagnet& magnet, const simd::float_v txT, const simd::float_v txV,
                        const simd::float_v tyV, simd::float_v& qOverP, simd::float_v& sigmaQOverP,
                        bool tCubicFit ) const override;

  StatusCode calculate( const DeMagnet& magnet, const simd::float_v txT, const simd::float_v tyT,
                        const simd::float_v xT, const simd::float_v zT, simd::float_v& qOverP,
                        simd::float_v& sigmaQOverP, bool tCubicFit ) const override;

private:
  template <typename floatType>
  StatusCode calculateImpl( const DeMagnet& magnet, const floatType txT, const floatType txV, const floatType tyV,
                            floatType& qOverP, floatType& sigmaQOverP, bool tCubicFit ) const;

  template <typename floatType>
  StatusCode calculateImpl( const DeMagnet& magnet, const floatType txT, const floatType tyT, const floatType xT,
                            const floatType zT, floatType& qOverP, floatType& sigmaQOverP, bool tCubicFit ) const;

  /// Momentum parameters (coefficients) for T state when using cubic fit
  Gaudi::Property<std::vector<double>> m_paramsTCubic{this, "ParamsTCubic", {0.0}};
  /// Momentum parameters (coefficients) for T state when using parabolic fit
  Gaudi::Property<std::vector<double>> m_paramsTParab{this, "ParamsTParabola", {0.0}};
  /// Momentum parameters (coefficients) for Velo and T state when using cubic fit
  Gaudi::Property<std::vector<double>> m_paramsVeloTCubic{this, "ParamsVeloTCubic", {0.0}};
  /// Momentum parameters (coefficients) for Velo and T state when using parabolic fit
  Gaudi::Property<std::vector<double>> m_paramsVeloTParab{this, "ParamsVeloTParabola", {0.0}};
  Gaudi::Property<double>              m_tResolution{this, "TResolution", 0.025};
  Gaudi::Property<double>              m_veloPlusTResolution{this, "VeloPlusTResolution", 0.015};

  bool m_useDefaultParametrization = false;

  LHCb::span<const double> getTParams( bool realMap, bool cubic ) const;
  LHCb::span<const double> getVeloParams( bool realMap, bool cubic ) const;
};

DECLARE_COMPONENT( FastMomentumEstimate )

namespace {
  static constexpr std::array<double, 4> DEFAULT_REAL_T_PARAB       = {-6.30991, -4.83533, -12.9192, 4.23025e-08};
  static constexpr std::array<double, 6> DEFAULT_REAL_VELO_PARAB    = {1.20812,  0.636694, -0.251334,
                                                                    0.414017, 2.87247,  -20.0982};
  static constexpr std::array<double, 4> DEFAULT_REAL_T_CUBIC       = {-6.34025, -4.85287, -12.4491, 4.25461e-08};
  static constexpr std::array<double, 6> DEFAULT_REAL_VELO_CUBIC    = {1.21352,  0.626691, -0.202483,
                                                                    0.426262, 2.47057,  -13.2917};
  static constexpr std::array<double, 4> DEFAULT_NO_REAL_T_PARAB    = {-6.3453, -4.77725, -14.9039, 3.13647e-08};
  static constexpr std::array<double, 6> DEFAULT_NO_REAL_VELO_PARAB = {1.21909,  0.627841, -0.235216,
                                                                       0.433811, 2.92798,  -21.3909};
  static constexpr std::array<double, 4> DEFAULT_NO_REAL_T_CUBIC    = {-6.31652, -4.46153, -16.694, 2.55588e-08};
  static constexpr std::array<double, 6> DEFAULT_NO_REAL_VELO_CUBIC = {1.21485,  0.64199, -0.27158,
                                                                       0.440325, 2.9191,  -20.4831};
} // namespace

//=============================================================================
// Initialization
//=============================================================================
StatusCode FastMomentumEstimate::initialize() {

  return extends::initialize().andThen( [&] {
    //@FIXME: why, if one of the properties is not properly set, are all
    //        four of them preset to these defaults???
    m_useDefaultParametrization = 4 != m_paramsTParab.size() || 4 != m_paramsTCubic.size() ||
                                  6 != m_paramsVeloTParab.size() || 6 != m_paramsVeloTCubic.size();

    return StatusCode::SUCCESS;
  } );
}

LHCb::span<const double> FastMomentumEstimate::getTParams( bool realMap, bool cubic ) const {
  if ( m_useDefaultParametrization ) {
    return realMap ? ( cubic ? DEFAULT_REAL_T_CUBIC : DEFAULT_REAL_T_PARAB )
                   : ( cubic ? DEFAULT_NO_REAL_T_CUBIC : DEFAULT_NO_REAL_T_PARAB );
  } else {
    return ( cubic ? m_paramsTCubic : m_paramsTParab ).value();
  }
}

LHCb::span<const double> FastMomentumEstimate::getVeloParams( bool realMap, bool cubic ) const {
  if ( m_useDefaultParametrization ) {
    return realMap ? ( cubic ? DEFAULT_REAL_VELO_CUBIC : DEFAULT_REAL_VELO_PARAB )
                   : ( cubic ? DEFAULT_NO_REAL_VELO_CUBIC : DEFAULT_NO_REAL_VELO_PARAB );
  } else {
    return ( cubic ? m_paramsVeloTCubic : m_paramsVeloTParab ).value();
  }
}

//=============================================================================
// Calculate momentum, given T state only
//=============================================================================
template <typename floatType>
StatusCode FastMomentumEstimate::calculateImpl( const DeMagnet& magnet, const floatType txT, const floatType tyT,
                                                const floatType xT, const floatType zT, floatType& qOverP,
                                                floatType& sigmaQOverP, bool tCubicFit ) const {
  const auto x0 = xT - txT * zT;

  const auto params = getTParams( magnet.useRealMap(), tCubicFit );

  const auto p = params[0] + params[1] * txT * txT + params[2] * tyT * tyT + params[3] * x0 * x0;

  const auto      scaleFactor = magnet.signedRelativeCurrent();
  const floatType denom       = p * scaleFactor * 1e6 * ( -1 );

  if ( std::abs( scaleFactor ) < 1e-6 ) {
    qOverP      = 0.01 / Gaudi::Units::GeV;
    sigmaQOverP = 1.0 / Gaudi::Units::MeV;
  } else {
    qOverP      = x0 / denom;
    sigmaQOverP = m_tResolution.value() * abs( qOverP );
  }

  return StatusCode::SUCCESS;
}

StatusCode FastMomentumEstimate::calculate( const DeMagnet& magnet, const LHCb::State* tState, double& qOverP,
                                            double& sigmaQOverP, bool tCubicFit ) const {
  const double tx = tState->tx();
  const double ty = tState->ty();
  const double x  = tState->x();
  const double z  = tState->z();
  return calculateImpl( magnet, tx, ty, x, z, qOverP, sigmaQOverP, tCubicFit );
}

StatusCode FastMomentumEstimate::calculate( const DeMagnet& magnet, const simd::float_v txT, const simd::float_v tyT,
                                            const simd::float_v xT, const simd::float_v zT, simd::float_v& qOverP,
                                            simd::float_v& sigmaQOverP, bool tCubicFit ) const {
  return calculateImpl( magnet, txT, tyT, xT, zT, qOverP, sigmaQOverP, tCubicFit );
}

//=============================================================================
// Calculate momentum, given T and Velo state (for matching)
//=============================================================================
template <typename floatType>
StatusCode FastMomentumEstimate::calculateImpl( const DeMagnet& magnet, const floatType txT, const floatType txV,
                                                const floatType tyV, floatType& qOverP, floatType& sigmaQOverP,
                                                bool tCubicFit ) const {

  const auto params = getVeloParams( magnet.useRealMap(), tCubicFit );

  const auto coef = params[0] + params[1] * txT * txT + params[2] * txT * txT * txT * txT + params[3] * txT * txV +
                    params[4] * tyV * tyV + params[5] * tyV * tyV * tyV * tyV;

  const auto proj        = sqrt( ( 1. + txV * txV + tyV * tyV ) / ( 1. + txV * txV ) );
  const auto scaleFactor = magnet.signedRelativeCurrent();

  if ( std::abs( scaleFactor ) < 1e-6 ) {
    qOverP      = 1.0 / Gaudi::Units::GeV;
    sigmaQOverP = 1.0 / Gaudi::Units::MeV;
  } else {
    qOverP      = ( txV - txT ) / ( coef * Gaudi::Units::GeV * proj * scaleFactor );
    sigmaQOverP = m_veloPlusTResolution.value() * abs( qOverP );
  }
  return StatusCode::SUCCESS;
}

StatusCode FastMomentumEstimate::calculate( const DeMagnet& magnet, const LHCb::State* veloState,
                                            const LHCb::State* tState, double& qOverP, double& sigmaQOverP,
                                            bool tCubicFit ) const {
  const auto txT = tState->tx();
  const auto txV = veloState->tx();
  const auto tyV = veloState->ty();

  return calculateImpl( magnet, txT, txV, tyV, qOverP, sigmaQOverP, tCubicFit );
}

StatusCode FastMomentumEstimate::calculate( const DeMagnet& magnet, const simd::float_v txT, const simd::float_v txV,
                                            const simd::float_v tyV, simd::float_v& qOverP, simd::float_v& sigmaQOverP,
                                            bool tCubicFit ) const {
  return calculateImpl( magnet, txT, txV, tyV, qOverP, sigmaQOverP, tCubicFit );
}
