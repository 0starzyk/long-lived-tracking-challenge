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
#ifndef KALMANPARAMETRIZATIONSCOEF_H
#define KALMANPARAMETRIZATIONSCOEF_H 1

// Include files
#include <array>
#include <fstream>
#include <stdio.h>

#include "GaudiKernel/GaudiException.h"
#include "GaudiKernel/System.h"

/** @class KalmanParametrizationsCoef KalmanParametrizationsCoef.h
 *  Class to handle coefficients of the extrapolation of a state through the magnet
 *
 *
 *  @author Pierre Billoir, Simon Stemmle
 *  @date   2018-07-06
 */

class KalmanParametrizationsCoef {
public:
  int                   Degx1, Degx2, Degy1, Degy2;
  std::array<float, 10> x00, x10, x01, tx00, tx10, tx01, y00, y10, y01, ty00, ty10, ty01;
  // kacoef() {};
  //~Coef() {};
  int Read( std::istream& inFile, int degx1, int degx2, int degy1, int degy2 );
};

KalmanParametrizationsCoef operator+( KalmanParametrizationsCoef a, KalmanParametrizationsCoef b );

KalmanParametrizationsCoef operator-( KalmanParametrizationsCoef a, KalmanParametrizationsCoef b );

KalmanParametrizationsCoef operator*( KalmanParametrizationsCoef a, double p );

#endif
