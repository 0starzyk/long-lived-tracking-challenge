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
#ifndef TRACKINTERFACES_ITRACKKALMANFILTER_H
#define TRACKINTERFACES_ITRACKKALMANFILTER_H

// Include files
// -------------

// from Gaudi
#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

// Forward declarations

/** @class ITrackKalmanFilter ITrackKalmanFilter.h TrackInterfaces/ITrackKalmanFilter.h
 *
 *  Interface for a track fitting tool.
 *
 *  @author Jose A. Hernando, Eduardo Rodrigues
 *  @date   2005-05-25
 *
 *  @author Rutger van der Eijk  07-04-1999
 *  @author Mattiew Needham
 */
struct ITrackKalmanFilter : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackKalmanFilter, 2, 0 );

  //! fit a track
  virtual StatusCode fit( LHCb::Track& track ) const = 0;
};
#endif // TRACKINTERFACES_ITRACKFITTER_H
