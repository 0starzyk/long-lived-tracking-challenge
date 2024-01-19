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
// Include
#include "KalmanParametrizationsCoef.h"
#include <stdio.h>

//##################################################################################################
//
// Implementation file for class : KalmanParametrizationsCoef
//
// 2018-07-06: Pierre Billoir, Simon Stemmle
//
//##################################################################################################

int KalmanParametrizationsCoef::Read( std::istream& inp, int degx1, int degx2, int degy1, int degy2 ) {
  if ( degx1 > 10 || degx2 > 10 || degy1 > 10 || degy2 > 10 )
    throw GaudiException(
        "You have to increase the size of the internal arrays of the KalmanParametrizationsCoef class",
        "KalmanParametrizationsCoef", StatusCode::FAILURE );

  auto GetPar = []( std::istream& ss ) -> float {
    std::string s;
    ss >> s;
    // the following is only necessary at the moment because some of the parameters still have nonsense values
    // This will be fixed with a future version of the parameters
    // For the moment I think it is not worth releasing another DataPackage ParamFiles version. I am sorry for that..
    try {
      return std::stof( s );
    } catch ( const std::out_of_range& e ) { return 0; }
    return std::stof( s );
  };

  Degx1 = degx1;
  Degx2 = degx2;
  Degy1 = degy1;
  Degy2 = degy2;
  for ( int i = 0; i < degx2; i++ ) { x00[i] = GetPar( inp ); }
  for ( int i = 0; i < degx2; i++ ) { tx00[i] = 1e-3 * GetPar( inp ); }
  for ( int i = 0; i < degx1; i++ ) { x10[i] = GetPar( inp ); }
  for ( int i = 0; i < degx1; i++ ) { x01[i] = GetPar( inp ); }
  for ( int i = 0; i < degx1; i++ ) { tx10[i] = 1e-3 * GetPar( inp ); }
  for ( int i = 0; i < degx1; i++ ) { tx01[i] = 1e-3 * GetPar( inp ); }

  for ( int i = 0; i < degy2; i++ ) { y00[i] = GetPar( inp ); }
  for ( int i = 0; i < degy2; i++ ) { ty00[i] = 1e-3 * GetPar( inp ); }
  for ( int i = 0; i < degy1; i++ ) { y10[i] = GetPar( inp ); }
  for ( int i = 0; i < degy1; i++ ) { y01[i] = GetPar( inp ); }
  for ( int i = 0; i < degy1; i++ ) { ty10[i] = 1e-3 * GetPar( inp ); }
  for ( int i = 0; i < degy1; i++ ) { ty01[i] = 1e-3 * GetPar( inp ); }
  return 1;
}

KalmanParametrizationsCoef operator+( KalmanParametrizationsCoef a, KalmanParametrizationsCoef b ) {
  KalmanParametrizationsCoef c = a;
  for ( int i = 0; i < c.Degx2; i++ ) {
    c.x00[i] += b.x00[i];
    c.tx00[i] += b.tx00[i];
  }
  for ( int i = 0; i < c.Degx1; i++ ) {
    c.x10[i] += b.x10[i];
    c.x01[i] += b.x01[i];
    c.tx10[i] += b.tx10[i];
    c.tx01[i] += b.tx01[i];
  }
  for ( int i = 0; i < c.Degy2; i++ ) {
    c.y00[i] += b.y00[i];
    c.ty00[i] += b.ty00[i];
  }
  for ( int i = 0; i < c.Degy1; i++ ) {
    c.y10[i] += b.y10[i];
    c.y01[i] += b.y01[i];
    c.ty10[i] += b.ty10[i];
    c.ty01[i] += b.ty01[i];
  }
  return c;
}

KalmanParametrizationsCoef operator-( KalmanParametrizationsCoef a, KalmanParametrizationsCoef b ) {
  KalmanParametrizationsCoef c = a;
  for ( int i = 0; i < c.Degx2; i++ ) {
    c.x00[i] -= b.x00[i];
    c.tx00[i] -= b.tx00[i];
  }
  for ( int i = 0; i < c.Degx1; i++ ) {
    c.x10[i] -= b.x10[i];
    c.x01[i] -= b.x01[i];
    c.tx10[i] -= b.tx10[i];
    c.tx01[i] -= b.tx01[i];
  }
  for ( int i = 0; i < c.Degy2; i++ ) {
    c.y00[i] -= b.y00[i];
    c.ty00[i] -= b.ty00[i];
  }
  for ( int i = 0; i < c.Degy1; i++ ) {
    c.y10[i] -= b.y10[i];
    c.y01[i] -= b.y01[i];
    c.ty10[i] -= b.ty10[i];
    c.ty01[i] -= b.ty01[i];
  }
  return c;
}

KalmanParametrizationsCoef operator*( KalmanParametrizationsCoef a, double p ) {
  KalmanParametrizationsCoef c = a;
  for ( int i = 0; i < c.Degx2; i++ ) {
    c.x00[i] *= p;
    c.tx00[i] *= p;
  }
  for ( int i = 0; i < c.Degx1; i++ ) {
    c.x10[i] *= p;
    c.x01[i] *= p;
    c.tx10[i] *= p;
    c.tx01[i] *= p;
  }
  for ( int i = 0; i < c.Degy2; i++ ) {
    c.y00[i] *= p;
    c.ty00[i] *= p;
  }
  for ( int i = 0; i < c.Degy1; i++ ) {
    c.y10[i] *= p;
    c.y01[i] *= p;
    c.ty10[i] *= p;
    c.ty01[i] *= p;
  }
  return c;
}
