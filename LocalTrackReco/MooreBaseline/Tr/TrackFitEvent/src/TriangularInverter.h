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
#ifndef TRIANGULARINVERTER_H
#define TRIANGULARINVERTER_H

#include "Event/TrackTypes.h"
#include <cmath>

/** Some tools to deal with triangular matrices
 * or almost triangular matrices (ie. one of the
 * element just under diagonal is non null )
 **/
namespace TriangularInversion {

  //  Fast inversion of upper triangular matrix
  void invertUpperTriangular( const Gaudi::TrackMatrix& F, Gaudi::TrackMatrix& Finv );

  //  Fast inversion of almost upper triangular matrix
  //  almost means one that the element under diagonal at line "anomaly"
  //  is non-null.
  void invertAlmostUpperTriangular( const Gaudi::TrackMatrix& F, Gaudi::TrackMatrix& Finv, const int& anomaly );

  //  Compute similarity matrix of with upper triangular matrix
  void similarityUpperTriangular( const Gaudi::TrackMatrix& F, const Gaudi::TrackMatrix& orig,
                                  Gaudi::TrackSymMatrix& target );

  //  Compute similarity matrix of with almost upper triangular matrix
  //  almost means one that the element under diagonal at line "anomaly"
  //  is non-null.
  void similarityAlmostUpperTriangular( const Gaudi::TrackMatrix& F, const Gaudi::TrackMatrix& orig,
                                        Gaudi::TrackSymMatrix& target, const int& anomaly );
} // namespace TriangularInversion

#endif // TRIANGULARINVERTER_H
