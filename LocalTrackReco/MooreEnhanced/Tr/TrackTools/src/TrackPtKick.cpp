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
#include "Event/State.h"
#include "Event/TrackParameters.h"
#include "Event/TrackTypes.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Kernel/IBIntegrator.h"
#include "TrackInterfaces/ITrackMomentumEstimate.h"
#include <GaudiKernel/SystemOfUnits.h>
#include <cmath>

//-----------------------------------------------------------------------------
// Implementation file for class : TrackPtKick
//
// 2000-08-16 : M. Needham
// 2005-05-13 : J. Nardulli (adaptations to new track event model)
// 2006-07-24 : M Needham - tune for DC 06
//-----------------------------------------------------------------------------

/** @class TrackPtKick TrackPtKick.h TrackTools/TrackPtKick.h
 *
 *  @author M. Needham
 *  @date   2000-08-16
 */
class TrackPtKick : public extends<GaudiTool, ITrackMomentumEstimate> {
public:
  using extends::extends;

  StatusCode initialize() override;

  // Estimate the momentum P of a State in T at ZAtMidT
  StatusCode calculate( const DeMagnet& magnet, const LHCb::State* tState, double& qOverP, double& sigmaQOverP,
                        bool tCubicFit ) const override;

  // Estimate the momentum P of a velo State and a State in T at ZAtMidT
  StatusCode calculate( const DeMagnet& magnet, const LHCb::State* veloState, const LHCb::State* tState, double& qOverP,
                        double& sigmaQOverP,
                        bool    tCubicFit ) const override; // Estimate the momentum P of a State

  StatusCode calculate( const DeMagnet& magnet, const simd::float_v txT, const simd::float_v txV,
                        const simd::float_v tyV, simd::float_v& qOverP, simd::float_v& sigmaQOverP,
                        bool tCubicFit ) const override;

  StatusCode calculate( const DeMagnet& magnet, const simd::float_v txT, const simd::float_v tyT,
                        const simd::float_v xT, const simd::float_v zT, simd::float_v& qOverP,
                        simd::float_v& sigmaQOverP, bool tCubicFit ) const override;

private:
  PublicToolHandle<IBIntegrator> m_bIntegrator{this, "BFieldIntegrator", "BIntegrator"};

  /// Define the parameters of the Z dependance
  Gaudi::Property<std::vector<double>> m_ParabolicCorrection{this, "ParabolicCorrection", {1.04, 0.14}};
  Gaudi::Property<std::vector<double>> m_resParams{this, "resParams", {0.015, 0.29}};
  Gaudi::Property<double>              m_Constant{this, "ConstantCorrection", 0. * Gaudi::Units::MeV};
};

DECLARE_COMPONENT( TrackPtKick )

//=============================================================================
// Initialization
//=============================================================================
StatusCode TrackPtKick::initialize() {
  return extends::initialize().andThen( [&]() {
    info() << " Pt kick parameters(" << m_ParabolicCorrection.size() << ") ==" << m_ParabolicCorrection[0] << " + "
           << m_ParabolicCorrection[1] << " tx^2 " << endmsg;
  } );
}

//=============================================================================
// Estimate the momentum P of a State
//=============================================================================
StatusCode TrackPtKick::calculate( const DeMagnet&    magnet, const LHCb::State* /* veloState */,
                                   const LHCb::State* tState, double& qOverP, double& sigmaQOverP,
                                   bool tCubicFit ) const {
  return calculate( magnet, tState, qOverP, sigmaQOverP, tCubicFit );
}

//=============================================================================
// Estimate the momentum P of a State
//=============================================================================
StatusCode TrackPtKick::calculate( const DeMagnet& magnet, const LHCb::State* tState, double& qOverP,
                                   double& sigmaQOverP, bool /* tCubicFit */ ) const {
  // calculate intial estimate of track momentum assuming it came from
  // the primary vertex

  // scan in cm steps
  const Gaudi::XYZPoint begin{0., 0., 0.};
  Gaudi::XYZVector      bdl;
  double                zCenter;

  m_bIntegrator->calculateBdlAndCenter( magnet.fieldGrid(), begin, tState->position(), tState->tx(), tState->ty(),
                                        zCenter, bdl );

  double q = 0.;
  double p = 1e6 * Gaudi::Units::MeV;

  if ( fabs( bdl.x() ) > TrackParameters::hiTolerance ) {
    // can estimate momentum and charge

    // Rotate to the  0-0-z axis and do the ptkick
    const double tX      = tState->tx();
    const double xCenter = tState->x() + tX * ( zCenter - tState->z() );

    const double zeta_trk = -tX / sqrt( 1.0 + tX * tX );
    const double tx_vtx   = xCenter / zCenter;
    const double zeta_vtx = -tx_vtx / sqrt( 1.0 + tx_vtx * tx_vtx );

    // curvature
    const double curv = ( zeta_trk - zeta_vtx );

    // charge
    int sign = 1;
    if ( curv < TrackParameters::hiTolerance ) { sign *= -1; }
    if ( bdl.x() < TrackParameters::hiTolerance ) { sign *= -1; }
    q = -1. * ( magnet.isDown() ? -1 : 1 ) * sign;

    // momentum
    p = Gaudi::Units::eplus * Gaudi::Units::c_light * fabs( bdl.x() ) *
        sqrt( ( 1.0 + tX * tX + std::pow( tState->ty(), 2 ) ) / ( 1.0 + std::pow( tX, 2 ) ) ) / fabs( curv );

    //   Addition Correction factor for the angle of the track!
    if ( m_ParabolicCorrection.size() == 2u ) {
      // p*= (a + b*tx*tx )
      p += m_Constant;
      p *= ( m_ParabolicCorrection[0] + ( m_ParabolicCorrection[1] * tX * tX ) );
    }

  } else {
    // can't estimate momentum or charge
    error() << "B integral is 0!" << endmsg;
    return StatusCode::FAILURE;
  }

  qOverP      = q / p;
  sigmaQOverP = std::hypot( m_resParams[0], m_resParams[1] / p ) / p;

  return StatusCode::SUCCESS;
}

//=============================================================================
// These is here just to make interface happy
//=============================================================================
StatusCode TrackPtKick::calculate( const DeMagnet& /*magnet*/, const simd::float_v /* txT */,
                                   const simd::float_v /* txV */, const simd::float_v /* tyV */,
                                   simd::float_v& /* qOverP */, simd::float_v& /* sigmaQOverP */,
                                   bool /* tCubicFit */ ) const {
  return StatusCode::SUCCESS;
}
StatusCode TrackPtKick::calculate( const DeMagnet& /*magnet*/, const simd::float_v /* txT */,
                                   const simd::float_v /* tyT */, const simd::float_v /* xT */,
                                   const simd::float_v /* zT */, simd::float_v& /* qOverP */,
                                   simd::float_v& /* sigmaQOverP */, bool /* tCubicFit */ ) const {
  return StatusCode::SUCCESS;
}
