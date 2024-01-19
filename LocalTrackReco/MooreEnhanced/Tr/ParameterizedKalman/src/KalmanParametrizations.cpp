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
#include <fstream>
#include <iostream>

#include <cmath>
#include <string.h>
#include <vector>

#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/System.h"

#include "KalmanParametrizations.h"

//##################################################################################################
//
// Implementation file for class : KalmanParametrizations
//
// 2017-10-26: Simon Stemmle
//
//##################################################################################################

namespace {
  // Set a 5x5 diagonal matrix for later use
  std::array<double, 25> F_diag = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
} // namespace

////////////////////////////////////////////////////////////////////////////////////////////////////
// Set the parameters (if needed) for the given magnet polarity
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::SetParameters( const std::string& ParamFileLocation, const Polarity polarity,
                                            bool useOneParameterSet ) {
  if ( ( m_Polarity == polarity ) && paramsLoaded ) return;
  ////////////////////
  // Load parameters
  ////////////////////
  std::string pol = ( polarity == Polarity::Up ? "Up" : "Down" );

  // The down parameter set is the default one
  if ( useOneParameterSet ) pol = "Down";

  std::string parameterPath = ParamFileLocation + "/Mag" + pol;

  // read the parameters for parametrizations
  read_params( parameterPath + "/params_predictV.txt", Par_predictV );
  read_params( parameterPath + "/params_predictVUT.txt", Par_predictVUT );
  read_params( parameterPath + "/params_predictUT.txt", Par_predictUT );
  read_params( parameterPath + "/params_predictUTFUT.txt", Par_predictUTFUT );
  read_params( parameterPath + "/params_predictUTTF.txt", Par_predictUTTF );
  read_params( parameterPath + "/params_predictTFT.txt", Par_predictTFT );
  read_params( parameterPath + "/params_predictT.txt", Par_predictT );

  read_params( parameterPath + "/params_TLayer.txt", Par_TLayer );
  read_params( parameterPath + "/params_UTLayer.txt", Par_UTLayer );

  // Get the up parameters from the down parameters
  if ( useOneParameterSet && polarity == Polarity::Up ) {
    SwitchParamsForPolarity( Par_predictV, flip_Par_predictV );
    SwitchParamsForPolarity( Par_predictVUT, flip_Par_predictVUT );
    SwitchParamsForPolarity( Par_predictUT, flip_Par_predictUT );
    SwitchParamsForPolarity( Par_predictUTFUT, flip_Par_predictUTFUT );
    SwitchParamsForPolarity( Par_predictUTTF, flip_Par_predictUTTF );
    SwitchParamsForPolarity( Par_predictTFT, flip_Par_predictTFT );
    SwitchParamsForPolarity( Par_predictT, flip_Par_predictT );
  }

  read_params_UTT( ParamFileLocation + "/MagDown/v5r0_7957.tab" );
  if ( polarity == Polarity::Up ) m_qop_flip = true;

  m_Polarity   = polarity;
  paramsLoaded = true;
}

// This contains the new version of Pierres extrapolation
// It uses quadratic extrapolation between the x,y coefficients
#include "Extrapolations/propag_UT_SciFi_quad.icpp"

////////////////////////////////////////////////////////////////////////////////////////////////////
// This switches all parameters that linearly/cubicly/... depend on q/p
////////////////////////////////////////////////////////////////////////////////////////////////////
template <std::size_t SIZE, std::size_t BN, std::size_t BS>
void KalmanParametrizations::SwitchParamsForPolarity( KalmanParameters<BN, BS>&            params,
                                                      const std::array<unsigned int, SIZE> list ) {
  for ( unsigned int i = 0; i < params.batchN; i++ ) {
    for ( auto j : list ) params( i, j ) *= -1;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate inside the VELO
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateInV( double zFrom, double zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                                             Gaudi::SymMatrix5x5& Q ) const {
  // cache the old state
  Gaudi::Vector5 x_old = x;
  // step size in z
  double dz = zTo - zFrom;
  // which set of parameters should be used
  const auto par = Par_predictV[dz > 0 ? 0 : 1];

  // do not update if there is nothing to update
  if ( dz == 0 ) return;

  // parametrizations for state extrapolation
  // tx
  x[2] = x_old[2] + x_old[4] * par[4] * 1e-5 * dz * ( ( dz > 0 ? zFrom : zTo ) + par[5] * 1e3 );
  // x
  x[0] = x_old[0] + ( x[2] + x_old[2] ) * 0.5 * dz;
  // ty
  x[3] = x_old[3];
  // y
  x[1] = x_old[1] + x[3] * dz;
  // qop
  x[4] = x_old[4];

  // determine the Jacobian

  F.SetElements( F_diag.begin(), F_diag.end() );
  F( 0, 2 ) = dz;
  F( 1, 3 ) = dz;

  // tx
  F( 2, 4 ) = par[4] * 1e-5 * dz * ( ( dz > 0 ? zFrom : zTo ) + par[5] * 1e3 );

  // x
  F( 0, 4 ) = 0.5 * dz * F( 2, 4 );

  // Set noise matrix

  double sigt = par[1] * 1e-5 + par[2] * std::abs( x_old[4] );
  // sigma x/y
  double sigx = par[6] * sigt * std::abs( dz );
  // Correlation between x/y and tx/ty
  double corr = par[7];

  Q( 0, 0 ) = sigx * sigx;
  Q( 1, 1 ) = sigx * sigx;
  Q( 2, 2 ) = sigt * sigt;
  Q( 3, 3 ) = sigt * sigt;

  Q( 0, 2 ) = corr * sigx * sigt;
  Q( 1, 3 ) = corr * sigx * sigt;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate VELO <-> UT
////////////////////////////////////////////////////////////////////////////////////////////////////
bool KalmanParametrizations::ExtrapolateVUT( double zFrom, double zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                                             Gaudi::SymMatrix5x5& Q ) const {
  // cache the old state
  Gaudi::Vector5 x_old = x;
  // step size in z
  double dz = zTo - zFrom;
  // which set of parameters should be used
  const auto par = Par_predictVUT[dz > 0 ? 0 : 1];
  // extrapolate the current state and define noise
  if ( dz > 0 ) {
    // ty
    x[3] = x_old[3] + par[0] * std::copysign( 1.0, x[1] ) * x_old[4] * x_old[2];

    double tyErr = par[3] * std::fabs( x_old[4] );

    // y
    x[1] = x_old[1] + ( par[5] * x_old[3] + ( 1 - par[5] ) * x[3] ) * dz;

    double yErr = par[6] * std::abs( dz * x_old[4] );

    // tx
    double coeff = par[8] * 1e1 + par[9] * 1e-2 * zFrom + par[10] * 1e2 * x_old[3] * x_old[3];

    double a = x_old[2] / std::sqrt( 1.0 + x_old[2] * x_old[2] + x_old[3] * x_old[3] ) - x_old[4] * coeff;

    // Check that the track is not deflected
    if ( std::fabs( a ) >= 1 ) return false;

    x[2] = a * sqrt( 1.0 / ( 1.0 - a * a ) * ( 1.0 + x[3] * x[3] ) );

    double txErr = par[15] * std::fabs( x_old[4] );

    // x
    double zmag =
        par[16] * 1e3 + par[17] * zFrom + par[18] * 1e-5 * zFrom * zFrom + par[19] * 1e3 * x_old[3] * x_old[3];

    x[0] = x_old[0] + ( zmag - zFrom ) * x_old[2] + ( zTo - zmag ) * x[2];

    double xErr = par[20] * std::abs( dz * x_old[4] );

    // calculate jacobian
    // ty
    F( 3, 0 ) = 0;
    F( 3, 1 ) = 0;
    F( 3, 2 ) = par[0] * x_old[4];
    F( 3, 3 ) = 1;
    F( 3, 4 ) = par[0] * x_old[2];
    // y
    double DyDty = ( 1 - par[5] ) * dz;
    F( 1, 0 )    = 0.0;
    F( 1, 1 )    = 1.0;
    F( 1, 2 )    = DyDty * F( 3, 2 );
    F( 1, 3 )    = dz;
    F( 1, 4 )    = DyDty * F( 3, 4 );

    // tx
    double sqrtTmp = std::sqrt( ( 1 - a * a ) * ( 1 + x[3] * x[3] ) );
    double DtxDty  = a * x[3] * 1.0 / sqrtTmp;
    double DtxDa   = sqrtTmp / ( ( a * a - 1 ) * ( a * a - 1 ) );
    F( 2, 0 )      = 0;
    F( 2, 1 )      = 0;

    sqrtTmp   = std::sqrt( 1 + x_old[2] * x_old[2] + x_old[3] * x_old[3] );
    F( 2, 2 ) = DtxDa * ( 1 + x_old[3] * x_old[3] ) / ( sqrtTmp * ( 1 + x_old[2] * x_old[2] + x_old[3] * x_old[3] ) ) +
                DtxDty * F( 3, 2 );

    F( 2, 3 ) = DtxDa * ( -x_old[2] * x_old[3] / ( sqrtTmp * ( 1 + x_old[2] * x_old[2] + x_old[3] * x_old[3] ) ) -
                          x_old[4] * 2 * par[10] * 1e2 * x_old[3] ) +
                DtxDty * F( 3, 3 );

    F( 2, 4 ) = DtxDa * ( -coeff ) + DtxDty * F( 3, 4 );

    // x
    F( 0, 0 ) = 1;
    F( 0, 1 ) = 0;
    F( 0, 2 ) = ( zmag - zFrom ) + ( zTo - zmag ) * F( 2, 2 );

    F( 0, 3 ) = ( zTo - zmag ) * F( 2, 3 ) + ( x_old[2] - x[2] ) * 2 * par[19] * 1e3 * x_old[3];

    F( 0, 4 ) = ( zTo - zmag ) * F( 2, 4 );

    // qop
    F( 4, 0 ) = 0;
    F( 4, 1 ) = 0;
    F( 4, 2 ) = 0;
    F( 4, 3 ) = 0;
    F( 4, 4 ) = 1;

    // add noise
    Q( 0, 0 ) = xErr * xErr;
    Q( 0, 2 ) = par[4] * xErr * txErr;
    Q( 1, 1 ) = yErr * yErr;
    Q( 1, 3 ) = par[21] * yErr * tyErr;
    Q( 2, 2 ) = txErr * txErr;
    Q( 3, 3 ) = tyErr * tyErr;
  } else {
    // ty
    x[3] = x_old[3] + par[0] * std::copysign( 1.0, x[1] ) * x_old[4] * x_old[2];

    double tyErr = par[3] * std::fabs( x_old[4] );

    // y
    x[1] = x_old[1] + ( par[5] * x_old[3] + ( 1 - par[5] ) * x[3] ) * dz;

    double yErr = par[6] * std::abs( dz * x_old[4] );

    // tx
    double coeff = par[8] * 1e1 + par[9] * 1e-2 * zTo + par[10] * 1e2 * x_old[3] * x_old[3];

    double a = x_old[2] / std::sqrt( 1.0 + x_old[2] * x_old[2] + x_old[3] * x_old[3] ) - x_old[4] * coeff;

    // Check that the track is not deflected
    if ( std::fabs( a ) >= 1 ) return false;

    x[2]         = a * sqrt( 1.0 / ( 1.0 - a * a ) * ( 1.0 + x[3] * x[3] ) );
    double txErr = par[15] * std::fabs( x_old[4] );

    // x
    double zmag = par[16] * 1e3 + par[17] * zTo + par[18] * 1e-5 * zTo * zTo + par[19] * 1e3 * x_old[3] * x_old[3];

    x[0] = x_old[0] + ( zmag - zFrom ) * x_old[2] + ( zTo - zmag ) * x[2];

    double xErr = par[20] * std::abs( dz * x_old[4] );

    // calculate jacobian
    // ty
    F( 3, 0 ) = 0;
    F( 3, 1 ) = 0;
    F( 3, 2 ) = par[0] * x_old[4];
    F( 3, 3 ) = 1;
    F( 3, 4 ) = par[0] * x_old[2];
    // y
    double DyDty = ( 1 - par[5] ) * dz;
    F( 1, 0 )    = 0.0;
    F( 1, 1 )    = 1.0;
    F( 1, 2 )    = DyDty * F( 3, 2 );
    F( 1, 3 )    = dz;
    F( 1, 4 )    = DyDty * F( 3, 4 );

    // tx
    double sqrtTmp = std::sqrt( ( 1 - a * a ) * ( 1 + x[3] * x[3] ) );
    double DtxDty  = a * x[3] * 1.0 / sqrtTmp;
    double DtxDa   = sqrtTmp / ( ( a * a - 1 ) * ( a * a - 1 ) );
    F( 2, 0 )      = 0;
    F( 2, 1 )      = 0;

    sqrtTmp   = std::sqrt( 1 + x_old[2] * x_old[2] + x_old[3] * x_old[3] );
    F( 2, 2 ) = DtxDa * ( 1 + x_old[3] * x_old[3] ) / ( sqrtTmp * ( 1 + x_old[2] * x_old[2] + x_old[3] * x_old[3] ) ) +
                DtxDty * F( 3, 2 );

    F( 2, 3 ) = DtxDa * ( -x_old[2] * x_old[3] / ( sqrtTmp * ( 1 + x_old[2] * x_old[2] + x_old[3] * x_old[3] ) ) -
                          x_old[4] * 2 * par[10] * 1e2 * x_old[3] ) +
                DtxDty * F( 3, 3 );

    F( 2, 4 ) = DtxDa * ( -coeff ) + DtxDty * F( 3, 4 );

    // x
    F( 0, 0 ) = 1;
    F( 0, 1 ) = 0;
    F( 0, 2 ) = ( zmag - zFrom ) + ( zTo - zmag ) * F( 2, 2 );

    F( 0, 3 ) = ( zTo - zmag ) * F( 2, 3 ) + ( x_old[2] - x[2] ) * 2 * par[19] * 1e3 * x_old[3];

    F( 0, 4 ) = ( zTo - zmag ) * F( 2, 4 );

    // qop
    F( 4, 0 ) = 0;
    F( 4, 1 ) = 0;
    F( 4, 2 ) = 0;
    F( 4, 3 ) = 0;
    F( 4, 4 ) = 1;

    // add noise
    Q( 0, 0 ) = xErr * xErr;
    Q( 0, 2 ) = par[4] * xErr * txErr;
    Q( 1, 1 ) = yErr * yErr;
    Q( 1, 3 ) = par[21] * yErr * tyErr;
    Q( 2, 2 ) = txErr * txErr;
    Q( 3, 3 ) = tyErr * tyErr;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate VELO <-> UT (traj)
////////////////////////////////////////////////////////////////////////////////////////////////////
bool KalmanParametrizations::ExtrapolateVUT( double zFrom, ROOT::Math::XYZPoint point, ROOT::Math::XYZVector dir,
                                             double& zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                                             Gaudi::SymMatrix5x5& Q ) const {
  // next z position:
  // use the straigt line extrapolation in y

  double z0 = point.Z();
  double y0 = point.Y();

  double dy = dir.Y();
  double dz = dir.Z();
  zTo       = ( zFrom * x[3] * dz - z0 * dy - x[1] * dz + y0 * dz ) / ( x[3] * dz - dy );

  return ExtrapolateVUT( zFrom, zTo, x, F, Q );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Get noise for VELO <- UT
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::GetNoiseVUTBackw( double zFrom, double zTo, const Gaudi::Vector5& x,
                                               Gaudi::SymMatrix5x5& Q ) const {
  // step size in z
  double dz = zTo - zFrom;
  // which set of parameters should be used
  const auto par = Par_predictVUT[1];

  // ty
  double tyErr = par[3] * std::fabs( x[4] );

  // y
  double yErr = par[6] * std::abs( dz * x[4] );

  // tx
  double txErr = par[15] * std::fabs( x[4] );

  // x
  double xErr = par[20] * std::abs( dz * x[4] );

  // add noise
  Q( 0, 0 ) = xErr * xErr;
  Q( 0, 2 ) = par[4] * xErr * txErr;
  Q( 1, 1 ) = yErr * yErr;
  Q( 1, 3 ) = par[21] * yErr * tyErr;
  Q( 2, 2 ) = txErr * txErr;
  Q( 3, 3 ) = tyErr * tyErr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Predict UT <-> UT (traj)
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateInUT( double zFrom, int nLayer, ROOT::Math::XYZPoint point,
                                              ROOT::Math::XYZVector dir, double& zTo, Gaudi::Vector5& x,
                                              Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const {
  // next z position:
  // use the straigt line extrapolation in y
  double z0 = point.Z();
  double y0 = point.Y();

  double dy = dir.Y();
  double dz = dir.Z();
  zTo       = ( zFrom * x[3] * dz - z0 * dy - x[1] * dz + y0 * dz ) / ( x[3] * dz - dy );

  ExtrapolateInUT( zFrom, nLayer, zTo, x, F, Q );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Predict UT <-> UT
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateInUT( double zFrom, int nLayer, double zTo, Gaudi::Vector5& x,
                                              Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const {
  // In case no specific z-position is set the default z position is used
  if ( zFrom == zTo ) zTo = Par_UTLayer[0][nLayer];
  // cache the old state
  Gaudi::Vector5 x_old = x;
  // step size in z
  double dz = zTo - zFrom;
  // which set of parameters should be used
  const auto par = Par_predictUT[( dz > 0 ? nLayer - 1 : ( 5 - nLayer ) )];

  // extrapolate state vector
  // tx
  x[2] += dz * ( par[5] * 1.e-1 * x[4] + par[6] * 1.e3 * x[4] * x[4] * x[4] + par[7] * 1e-7 * x[1] * x[1] * x[4] );
  // x
  x[0] += dz * ( par[0] * x_old[2] + ( 1 - par[0] ) * x[2] );
  // ty
  x[3] += par[10] * x[4] * x[2] * std::copysign( 1.0, x[1] );
  // y
  x[1] += dz * ( par[3] * x_old[3] + ( 1 - par[3] ) * x[3] );

  F( 2, 0 ) = 0;
  F( 2, 1 ) = 2 * dz * par[7] * 1e-7 * x_old[1] * x[4];
  F( 2, 2 ) = 1;
  F( 2, 3 ) = 0;
  F( 2, 4 ) = dz * ( par[5] * 1.e-1 + 3 * par[6] * 1.e3 * x[4] * x[4] + par[7] * 1e-7 * x_old[1] * x_old[1] );

  F( 0, 0 ) = 1;
  F( 0, 1 ) = dz * ( 1 - par[0] ) * F( 2, 1 );
  F( 0, 2 ) = dz;
  F( 0, 3 ) = 0;
  F( 0, 4 ) = dz * ( 1 - par[0] ) * F( 2, 4 );

  F( 3, 0 ) = 0;
  F( 3, 1 ) = 0;
  F( 3, 2 ) = par[10] * x[4] * std::copysign( 1.0, x[1] );
  F( 3, 3 ) = 1;
  F( 3, 4 ) = par[10] * x[2] * std::copysign( 1.0, x[1] );

  F( 1, 0 ) = 0;
  F( 1, 1 ) = 1;
  F( 1, 2 ) = dz * ( 1 - par[3] ) * F( 3, 2 );
  F( 1, 3 ) = dz;
  F( 1, 4 ) = dz * ( 1 - par[3] ) * F( 3, 4 );

  F( 4, 0 ) = 0;
  F( 4, 1 ) = 0;
  F( 4, 2 ) = 0;
  F( 4, 3 ) = 0;
  F( 4, 4 ) = 1;

  // Define noise
  double xErr  = par[2] * std::fabs( dz * x_old[4] );
  double yErr  = par[4] * std::fabs( dz * x_old[4] );
  double txErr = par[12] * std::fabs( x_old[4] );
  double tyErr = par[15] * std::fabs( x_old[4] );

  // Add noise
  Q( 0, 0 ) = xErr * xErr;
  Q( 0, 2 ) = par[14] * xErr * txErr;
  Q( 1, 1 ) = yErr * yErr;
  Q( 1, 3 ) = par[17] * yErr * tyErr;
  Q( 2, 2 ) = txErr * txErr;
  Q( 3, 3 ) = tyErr * tyErr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate UT (fixed z) -> T (fixed z)
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateUTT( Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const {
  // cache old state
  Gaudi::Vector5 x_old = x;

  const auto par = Par_predictUTTF[0];

  // extrapolating from last UT layer (z=2642.5) to fixed z in T (z=7855)

  // determine the momentum at this state from the momentum saved in the state vector
  //(representing always the PV qop)
  double qopHere =
      x[4] + x[4] * 1e-4 * par[18] + x[4] * std::abs( x[4] ) * par[19]; // TODO make this a tuneable parameter

  // do the actual extrapolation
  double der_tx[4], der_ty[4], der_qop[4]; //, der_x[4], der_y[4];
  extrapUTT( UTTExtrBeginZ(), UTTExtrEndZ(), QUADRATICINTERPOLATION, x[0], x[1], x[2], x[3], qopHere, der_tx, der_ty,
             der_qop );

  // apply additional correction
  x[0] +=
      par[9] * x_old[4] * 1e2 + par[10] * x_old[4] * x_old[4] * 1e5 + par[11] * x_old[4] * x_old[4] * x_old[4] * 1e10;
  x[1] += par[3] * x_old[4] * 1e2;
  x[2] += par[6] * x_old[4] + par[7] * x_old[4] * x_old[4] * 1e5 + par[8] * x_old[4] * x_old[4] * x_old[4] * 1e8;
  x[3] += par[0] * x_old[4];

  // Set jacobian matrix
  // TODO study impact of der_x, der_y
  // ty
  F( 3, 0 ) = 0; // der_x[3];
  F( 3, 1 ) = 0; // der_y[3];
  F( 3, 2 ) = der_tx[3];
  F( 3, 3 ) = der_ty[3];
  F( 3, 4 ) = der_qop[3] * ( 1 + 2 * std::abs( x[4] ) * par[18] ) + par[0] + 2 * par[1] * x_old[4] * 1e5 +
              3 * par[2] * x_old[4] * x_old[4] * 1e8;
  // y
  F( 1, 0 ) = 0; // der_x[1];
  F( 1, 1 ) = 1; // der_y[1];
  F( 1, 2 ) = der_tx[1];
  F( 1, 3 ) = der_ty[1];
  F( 1, 4 ) = der_qop[1] * ( 1 + 2 * std::abs( x[4] ) * par[18] ) + par[3] * 1e2 + 2 * par[4] * x_old[4] * 1e5 +
              3 * par[5] * x_old[4] * x_old[4] * 1e8;

  // tx
  F( 2, 0 ) = 0; // der_x[2];
  F( 2, 1 ) = 0; // der_y[2];
  F( 2, 2 ) = der_tx[2];
  F( 2, 3 ) = der_ty[2];
  F( 2, 4 ) = der_qop[2] * ( 1 + 2 * std::abs( x[4] ) * par[18] ) + par[6] + 2 * par[7] * x_old[4] * 1e5 +
              3 * par[8] * x_old[4] * x_old[4] * 1e8;

  // x
  F( 0, 0 ) = 1; // der_x[0];
  F( 0, 1 ) = 0; // der_y[0];
  F( 0, 2 ) = der_tx[0];
  F( 0, 3 ) = der_ty[0];
  F( 0, 4 ) = der_qop[0] * ( 1 + 2 * std::abs( x[4] ) * par[18] ) + par[9] * 1e2 + 2 * par[10] * x_old[4] * 1e5 +
              3 * par[11] * x_old[4] * x_old[4] * 1e10;

  // qop
  F( 4, 0 ) = 0;
  F( 4, 1 ) = 0;
  F( 4, 2 ) = 0;
  F( 4, 3 ) = 0;
  F( 4, 4 ) = 1;

  // Define noise
  double xErr  = par[13] * 1e2 * std::abs( x_old[4] );
  double yErr  = par[16] * 1e2 * std::abs( x_old[4] );
  double txErr = par[12] * std::abs( x_old[4] );
  double tyErr = par[15] * std::abs( x_old[4] );

  // Add noise
  Q( 0, 0 ) = xErr * xErr;
  Q( 0, 2 ) = par[14] * xErr * txErr;
  Q( 1, 1 ) = yErr * yErr;
  Q( 1, 3 ) = par[17] * yErr * tyErr;
  Q( 2, 2 ) = txErr * txErr;
  Q( 3, 3 ) = tyErr * tyErr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Get noise for UT (fixed z) <- T (fixed z)
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::GetNoiseUTTBackw( const Gaudi::Vector5& x, Gaudi::SymMatrix5x5& Q ) const {
  const auto par = Par_predictUTTF[1];

  // Define noise
  double xErr  = par[13] * 1e2 * std::abs( x[4] );
  double yErr  = par[16] * 1e2 * std::abs( x[4] );
  double txErr = par[12] * std::abs( x[4] );
  double tyErr = par[15] * std::abs( x[4] );

  // Add noise
  Q( 0, 0 ) = xErr * xErr;
  Q( 0, 2 ) = par[14] * xErr * txErr;
  Q( 1, 1 ) = yErr * yErr;
  Q( 1, 3 ) = par[17] * yErr * tyErr;
  Q( 2, 2 ) = txErr * txErr;
  Q( 3, 3 ) = tyErr * tyErr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate T <-> T (traj)
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateInT( double zFrom, int nLayer, ROOT::Math::XYZPoint point,
                                             ROOT::Math::XYZVector dir, double& zTo, Gaudi::Vector5& x,
                                             Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const {
  // determine next z position:
  // use the straigt line extrapolation in y
  // if there is a hit in the first layer: use it's trajectory

  double z0 = point.Z();
  double y0 = point.Y();
  double dy = dir.Y();
  double dz = dir.Z();

  double numerator   = dz * ( zFrom * x[3] - x[1] + y0 ) - z0 * dy;
  double denominator = 1.0f / ( x[3] * dz - dy );

  zTo = numerator * denominator;
  // TODO use this derivatives: Tested: it does not help. Remove it at some point!
  double DzDy  = -dz * denominator;
  double DzDty = dz * ( zFrom * denominator - numerator * ( denominator * denominator ) );

  ExtrapolateInT( zFrom, nLayer, zTo, DzDy, DzDty, x, F, Q );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate T <-> T (No hit)
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateInT( double zFrom, int nLayer, double& zTo, Gaudi::Vector5& x,
                                             Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const {
  // determine next z position:
  // use the straigt line extrapolation in y
  // and calculate the intersection with the detector layer
  double z0   = Par_TLayer[0][nLayer];
  double y0   = 0;
  double dydz = Par_TLayer[1][nLayer];
  zTo         = ( zFrom * x[3] - z0 * dydz - x[1] + y0 ) / ( x[3] - dydz );
  // TODO use this derivatives: Tested: it does not help. Remove it at some point!
  double DzDy = -1.0 / ( x[3] - dydz );
  double DzDty =
      zFrom / ( x[3] - dydz ) - ( zFrom * x[3] - z0 * dydz - x[1] + y0 ) / ( ( x[3] - dydz ) * ( x[3] - dydz ) );

  ExtrapolateInT( zFrom, nLayer, zTo, DzDy, DzDty, x, F, Q );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate T <-> T
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateInT( double zFrom, int nLayer, double zTo, double DzDy, double DzDty,
                                             Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const {
  // cache the old state
  Gaudi::Vector5 x_old = x;
  // step size in z
  double dz = zTo - zFrom;
  // which set of parameters should be used
  // Reminder: backward T station label is different for the iPar definition
  // 44 42 40 38     36 34 32 30    28 26 24
  //|  |  |  |      |  |  |  |     |  |  |  |
  //|  |  |  |      |  |  |  |     |  |  |  |
  // 45 43 41 39     37 35 33 31    29 27 25
  int iPar = ( dz > 0 ? 2 * nLayer - 2 : ( 42 - 2 * nLayer ) );
  if ( x[1] < 0 ) iPar += 1;
  const auto par = Par_predictT[iPar];

  // predict state
  // tx
  x[2] += dz * ( par[5] * 1.e-1 * x[4] + par[6] * 1.e3 * x[4] * x[4] * x[4] + par[7] * 1e-7 * x[1] * x[1] * x[4] );
  // x
  x[0] += dz * ( par[0] * x_old[2] + ( 1 - par[0] ) * x[2] );
  // ty
  x[3] += par[10] * x[4] * x[4] * x[1];
  // y
  x[1] += dz * ( par[3] * x_old[3] + ( 1 - par[3] ) * x[3] );

  // calculate jacobian

  double dtxddz = par[5] * 1.e-1 * x[4] + par[6] * 1.e3 * x[4] * x[4] * x[4] + par[7] * 1e-7 * x[1] * x[1] * x[4];

  F( 2, 0 ) = 0;
  F( 2, 1 ) = 2 * dz * par[7] * 1e-7 * x_old[1] * x[4] + dtxddz * DzDy;
  F( 2, 2 ) = 1;
  F( 2, 3 ) = dtxddz * DzDty;
  F( 2, 4 ) = dz * ( par[5] * 1.e-1 + 3 * par[6] * 1.e3 * x[4] * x[4] + par[7] * 1e-7 * x_old[1] * x_old[1] );

  double dxddz = par[0] * x_old[2] + ( 1 - par[0] ) * x[2];
  F( 0, 0 )    = 1;
  F( 0, 1 )    = dz * ( 1 - par[0] ) * F( 2, 1 ) + dxddz * DzDy;
  F( 0, 2 )    = dz;
  F( 0, 3 )    = dz * ( 1 - par[0] ) * F( 2, 3 ) + dxddz * DzDty;
  F( 0, 4 )    = dz * ( 1 - par[0] ) * F( 2, 4 );

  F( 3, 0 ) = 0;
  F( 3, 1 ) = 0;
  F( 3, 2 ) = 0;
  F( 3, 3 ) = 1;
  F( 3, 4 ) = 2 * par[10] * x[4];

  F( 1, 0 ) = 0;
  F( 1, 1 ) = 1;
  F( 1, 2 ) = 0;
  F( 1, 3 ) = dz;
  F( 1, 4 ) = dz * ( 1 - par[3] ) * F( 3, 4 );

  F( 4, 0 ) = 0;
  F( 4, 1 ) = 0;
  F( 4, 2 ) = 0;
  F( 4, 3 ) = 0;
  F( 4, 4 ) = 1;

  // Define noise
  double xErr  = par[2] * std::fabs( dz * x_old[4] );
  double yErr  = par[4] * std::fabs( dz * x_old[4] );
  double txErr = par[12] * std::fabs( x_old[4] );
  double tyErr = par[15] * std::fabs( x_old[4] );

  Q( 0, 0 ) = xErr * xErr;
  Q( 0, 2 ) = par[14] * xErr * txErr;
  Q( 1, 1 ) = yErr * yErr;
  Q( 1, 3 ) = par[17] * yErr * tyErr;
  Q( 2, 2 ) = txErr * txErr;
  Q( 3, 3 ) = tyErr * tyErr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate UT(fixed z) <-> last UT layer (traj)
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateUTFUT( double& zFrom, ROOT::Math::XYZPoint point, ROOT::Math::XYZVector dir,
                                               Gaudi::Vector5& x, Gaudi::Matrix5x5& F ) const {
  // determine next z position:
  // use the straigt line extrapolation in y
  // if there is a hit in the first layer: use it's trajector

  double z0 = point.Z();
  double y0 = point.Y();

  double zTo;

  double dy = dir.Y();
  double dz = dir.Z();
  zTo       = ( zFrom * x[3] * dz - z0 * dy - x[1] * dz + y0 * dz ) / ( x[3] * dz - dy );

  ExtrapolateUTFUT( zFrom, zTo, x, F );
  zFrom = zTo;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate UT to start point of UTTF extrapolation
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateUTFUTDef( double& zFrom, Gaudi::Vector5& x, Gaudi::Matrix5x5& F ) const {
  // Use the start position of the UTTF extrapolation as default z value
  ExtrapolateUTFUT( zFrom, Par_UTLayer[0][3], x, F );
  zFrom = Par_UTLayer[0][3];
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate around last UT layer
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateUTFUT( double zFrom, double zTo, Gaudi::Vector5& x,
                                               Gaudi::Matrix5x5& F ) const {
  // cache the old state
  Gaudi::Vector5 x_old = x;
  // step size in z
  double dz = zTo - zFrom;
  // which parameters should be used?
  const auto par = Par_predictUTFUT[0];

  // do the extrapolation of the state vector
  // tx
  x[2] = x_old[2] + par[0] * x_old[4] * dz;
  // x
  x[0] = x_old[0] + ( x[2] + x_old[2] ) * 0.5 * dz;
  // y
  x[1] = x_old[1] + x_old[3] * dz;

  // Jacobian
  F.SetElements( F_diag.begin(), F_diag.end() );

  // tx
  F( 2, 4 ) = par[0] * dz;
  // x
  F( 0, 2 ) = dz;
  F( 0, 4 ) = 0.5 * dz * F( 2, 4 );
  // y
  F( 1, 3 ) = dz;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate T(fixed z) <-> first T layer (traj)
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateTFT( double zFrom, ROOT::Math::XYZPoint point, ROOT::Math::XYZVector dir,
                                             double& zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                                             Gaudi::SymMatrix5x5& Q ) const {
  // determine next z position:
  // use the straigt line extrapolation in y
  // if there is a hit in the first layer: use it's trajector

  double z0 = point.Z();
  double y0 = point.Y();
  double dy = dir.Y();
  double dz = dir.Z();

  zTo = ( zFrom * x[3] * dz - z0 * dy - x[1] * dz + y0 * dz ) / ( x[3] * dz - dy );
  ExtrapolateTFT( zFrom, zTo, x, F, Q );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate T(fixed z) <-> first T layer (no hit)
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateTFTDef( double zFrom, double& zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                                                Gaudi::SymMatrix5x5& Q ) const {
  // determine next z position:
  // use the straigt line extrapolation in y
  double z0   = UTTExtrEndZ();
  double y0   = 0;
  double dydz = Par_TLayer[1][0];
  zTo         = ( zFrom * x[3] - z0 * dydz - x[1] + y0 ) / ( x[3] - dydz );

  ExtrapolateTFT( zFrom, zTo, x, F, Q );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Extrapolate T(fixed z) <-> first T layer
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::ExtrapolateTFT( double zFrom, double zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                                             Gaudi::SymMatrix5x5& Q ) const {
  // cache the old state
  Gaudi::Vector5 x_old = x;
  // step size in z
  double dz = zTo - zFrom;
  // which parameters should be used?
  const auto par = Par_predictTFT[dz > 0 ? 0 : 1];

  // do the extrapolation of the state vector
  // tx
  x[2] = x_old[2] + par[5] * x_old[4] * dz + 1e4 * par[6] * x_old[4] * dz * x_old[4] * dz * x_old[4] * dz;
  // x
  x[0] = x_old[0] + ( ( 1 - par[8] ) * x[2] + par[8] * x_old[2] ) * dz;
  // ty
  x[3] = x_old[3] + par[0] * ( x_old[4] * dz ) * ( x_old[4] * dz );
  // y
  x[1] = x_old[1] + ( x[3] + x_old[3] ) * 0.5 * dz;
  // qop
  x[4] = x_old[4];

  // Jacobian
  F.SetElements( F_diag.begin(), F_diag.end() );
  F( 0, 2 ) = dz;
  F( 1, 3 ) = dz;

  // tx
  F( 2, 4 ) = par[5] * dz + 3 * 1e4 * par[6] * dz * dz * dz * x_old[4] * x_old[4];
  // x
  F( 0, 4 ) = ( 1 - par[8] ) * dz * F( 2, 4 );
  // ty
  F( 3, 4 ) = 2 * par[0] * x_old[4] * dz * dz;
  // y
  F( 1, 4 ) = 0.5 * dz * F( 3, 4 );

  // Set noise: none
  // Should be already initialized to 0
  Q( 0, 0 ) = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//  read parameters from file
////////////////////////////////////////////////////////////////////////////////////////////////////
template <std::size_t BN, std::size_t BS>
void KalmanParametrizations::read_params( std::string_view file, KalmanParameters<BN, BS>& params ) {
  // read new parameters
  std::string   line;
  std::ifstream myfile( std::string{file} );

  bool foundSet = false;
  if ( myfile.is_open() ) {
    int iSet = 0;
    while ( getline( myfile, line ) ) {
      // determine which parameterset the respective line of paramters belongs to
      foundSet = false;
      for ( unsigned int s = 0; s < params.batchN; s++ ) {
        std::stringstream ss;
        ss << "_" << s << "_";
        std::string str = ss.str();
        if ( line.find( str ) != std::string::npos ) {
          iSet     = s;
          foundSet = true;
        }
      }
      if ( !foundSet ) continue;
      // set the values
      std::istringstream iss( line );
      std::string        sub;
      iss >> sub;
      unsigned int p = 0;
      while ( iss >> sub && p < params.batchS ) {
        params( iSet, p ) = std::atof( sub.c_str() );
        p++;
      }
    }
    myfile.close();
  } else
    throw GaudiException( std::string{"Failed to set the parameters from file "}.append( file ),
                          "KalmanParametrizations", StatusCode::FAILURE );

  // std::cout << file << std::endl;
  // for(int i=0; i<params.batchN; i++){
  //  for(int j=0; j<params.batchS; j++){
  //    std::cout << params[i][j] << " ";
  //  }
  //  std::cout << std::endl;
  //}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//  read parameters from file - here for the extrapolation UT -> T
////////////////////////////////////////////////////////////////////////////////////////////////////
void KalmanParametrizations::read_params_UTT( std::string_view file ) {
  std::string   line;
  std::ifstream myfile( std::string{file} );
  if ( !myfile.is_open() )
    throw GaudiException( std::string{"Failed to set the parameters from file "}.append( file ),
                          "KalmanParametrizations", StatusCode::FAILURE );

  myfile >> ZINI >> ZFIN >> PMIN >> BENDX >> BENDX_X2 >> BENDX_Y2 >> BENDY_XY >> Txmax >> Tymax >> XFmax >> Dtxy;
  myfile >> Nbinx >> Nbiny >> XGridOption >> YGridOption >> DEGX1 >> DEGX2 >> DEGY1 >> DEGY2;

  for ( int ix = 0; ix < Nbinx; ix++ )
    for ( int iy = 0; iy < Nbiny; iy++ ) C[ix][iy].Read( myfile, DEGX1, DEGX2, DEGY1, DEGY2 );

  Xmax = ZINI * Txmax;
  Ymax = ZINI * Tymax;
}

double KalmanParametrizations::UTTExtrEndZ() const { return ZFIN; }

double KalmanParametrizations::UTTExtrBeginZ() const { return ZINI; }

double KalmanParametrizations::VUTExtrEndZ() const { return Par_UTLayer[0][0]; }
