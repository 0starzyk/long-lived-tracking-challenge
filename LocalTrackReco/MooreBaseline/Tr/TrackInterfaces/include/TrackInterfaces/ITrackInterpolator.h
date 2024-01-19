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
#ifndef TRACKINTERFACES_ITRACKINTERPOLATOR_H
#define TRACKINTERFACES_ITRACKINTERPOLATOR_H 1

// Include files
// -------------
// from Gaudi
#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"
#include <DetDesc/IGeometryInfo.h>

// Forward declarations
namespace LHCb {
  class State;
}

/** @class ITrackInterpolator ITrackInterpolator.h TrackExtrapolators/ITrackInterpolator.h
 *
 *  Interface for track interpolators
 *
 *  @author Eduardo Rodrigues
 *  @date   2006-07-27
 */
struct ITrackInterpolator : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackInterpolator, 3, 0 );

  /// Interpolate a Track at a given z-position (the track may be re-fitted if needed!)
  virtual StatusCode interpolate( const LHCb::Track& track, double z, LHCb::State& state,
                                  IGeometryInfo const& geometry ) const = 0;
};
#endif // TRACKINTERFACES_ITRACKINTERPOLATOR_H
