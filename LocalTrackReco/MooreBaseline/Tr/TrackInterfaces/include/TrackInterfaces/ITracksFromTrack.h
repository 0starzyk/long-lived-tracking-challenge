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
#ifndef TRACKINTERFACES_ITRACKSFROMTRACK_H
#define TRACKINTERFACES_ITRACKSFROMTRACK_H 1

// Include files
// from STL
#include <vector>

// from Gaudi
#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

// Forward declarations

/** @class ITracksFromTrack ITracksFromTrack.h TrackInterfaces/ITracksFromTrack.h
 *  Interface to the forward pattern tool
 *
 *  @author David Hutchcroft
 *  @date   2007-05-24
 */
class ITracksFromTrack : public extend_interfaces<IAlgTool> {
public:
  DeclareInterfaceID( ITracksFromTrack, 2, 0 );

  /// Take an existing track and make new tracks from it (usually with hits from more detectors)
  virtual StatusCode tracksFromTrack( const LHCb::Track& seed, std::vector<LHCb::Track*>& tracks ) const = 0;
};
#endif // TRACKINTERFACES_ITRACKSFROMTRACK_H
