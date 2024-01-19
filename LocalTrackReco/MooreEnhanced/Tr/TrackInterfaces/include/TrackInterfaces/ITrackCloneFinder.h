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
#ifndef TRACKINTERFACES_ITRACKCLONEFINDER_H
#define TRACKINTERFACES_ITRACKCLONEFINDER_H 1

// Include files
// -------------

// from Gaudi
#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

// Forward declarations

/** @class ITrackCloneFinder ITrackCloneFinder.h TrackInterfaces/ITrackCloneFinder.h
 *
 *  Interface for the clone finder tool among two input tracks
 *
 *  @author Eduardo Rodrigues
 *  @date   2005-12-08
 *  Modified for speed reason
 *  @author Adrian Perieanu
 *  @date   2008-05-05
 */
struct ITrackCloneFinder : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackCloneFinder, 2, 0 );

  /** Compare two input Tracks and find whether one is a clone
   *  of the other based on some "overlap criteria".
   *  Note: the method ignores whether the Tracks themselves have been
   *        previously flagged as clones! It merely does a comparison.
   *  @param  track1 input 1st track
   *  @param  track2 input 2nd track
   *  @param  setFlag input parameter indicates whether the clone track
   *          is to be set as such (default = false)
   */
  virtual bool areClones( const LHCb::Track& track1, const LHCb::Track& track2 ) const = 0;

  /** Compare two input Tracks and find whether one is a clone
   *  of the other based on some "overlap criteria".
   *  The corresponding flag may be set accordingly (NOT DONE BY DEFAULT)
   *  depending on the value of the "setFlag" argument.
   *  Note: the method ignores whether the Tracks themselves have been
   *        previously flagged as clones! It merely does a comparison.
   *  @param  track1 input 1st track
   *  @param  track2 input 2nd track
   *  @param  setFlag input parameter indicates whether the clone track
   *          is to be set as such (default = false)
   */
  virtual bool flagClones( LHCb::Track& track1, LHCb::Track& track2 ) const = 0;
};
#endif // TRACKINTERFACES_ITRACKCLONEFINDER_H
