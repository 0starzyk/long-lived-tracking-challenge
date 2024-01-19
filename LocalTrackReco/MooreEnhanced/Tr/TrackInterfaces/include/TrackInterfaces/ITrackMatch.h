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
#ifndef TRACKINTERFACES_ITRACKMATCH_H
#define TRACKINTERFACES_ITRACKMATCH_H 1

// Include files
// from STL
#include <string>

// from Gaudi
#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"
#include <DetDesc/IGeometryInfo.h>

// forward declarations

/** @class ITrackMatch ITrackMatch.h TrackInterfaces/ITrackMatch.h
 *
 *
 *  @author Jose A. Hernando
 *  @date   2007-06-16
 */
struct ITrackMatch : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackMatch, 3, 0 );

  virtual StatusCode match( const LHCb::Track& track0, const LHCb::Track& track1, LHCb::Track& matchTrack,
                            IGeometryInfo const& geometry, double& quality, double& quality2 ) = 0;
};
#endif // TRACKINTERFACES_IMATCHTVELOTRACKS_H
