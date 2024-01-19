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
#ifndef TRACKINTERFACES_ITRACKCHI2CALCULATOR_H
#define TRACKINTERFACES_ITRACKCHI2CALCULATOR_H 1

// Include files
// -------------
#include "GaudiKernel/IAlgTool.h"

// From TrackEvent
#include "Event/TrackTypes.h"

/** @class ITrackChi2Calculator ITrackChi2Calculator.h TrackInterfaces/ITrackChi2Calculator.h
 *
 *  Interface class for the TrackChi2Calculator tool.
 *  This is used to calculate the chi2-distance between to track states.
 *
 *  @author Jeroen van Tilburg
 *  @date   2003-09-18
 */

struct ITrackChi2Calculator : extend_interfaces<IAlgTool> {
  DeclareInterfaceID( ITrackChi2Calculator, 2, 0 );

  /** Calculate the chi2 distance between two track vectors.
   *  The track vectors must be given in (x,y,tx,ty,q/p).
   *  @return StatusCode:   Failure if matrix inversion failed
   *  @param  trackVector1: input 1st track HepVector
   *  @param  trackCov1:    input covariance matrix corresponding to 1st vector
   *  @param  trackVector2: input 2nd track HepVector
   *  @param  trackCov2:    input covariance matrix corresponding to 2nd vector
   *  @param  chi2:         output chi2 distance between the two vectors
   */
  virtual StatusCode calculateChi2( const Gaudi::TrackVector& trackVector1, const Gaudi::TrackSymMatrix& trackCov1,
                                    const Gaudi::TrackVector& trackVector2, const Gaudi::TrackSymMatrix& trackCov2,
                                    double& chi2 ) const = 0;
};

#endif // TRACKINTERFACES_ITRACKCHI2CALCULATOR_H
