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
#ifndef PATALG_ITRACKSTATEINIT_H
#define PATALG_ITRACKSTATEINIT_H 1

// from Gaudi
#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"
#include "Kernel/TrackDefaultParticles.h"

/** @class ITrackStateInit ITrackStateInit.h
 *
 * An interface to the TrackStateInitTool
 *
 * @author Pavel Krokovny <krokovny@physi.uni-heidelberg.de>
 * @date   2009-03-02
 */

struct ITrackStateInit : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackStateInit, 2, 0 );

  /**
   * remove all states on the track, apply a simple fit only based on LHCbIds
   * and reinitalize a list of states. clearStates = false leaves the tracks
   * untouched for the moment , more features will be added here later.
   */
  virtual StatusCode fit( LHCb::Track& track, bool clearStates = true ) const = 0;
  /**
   * check if track states are initalized at a given list of z reference positions.
   * if list of state not complete, add additional states based on extrapolation of
   * states given on the track
   */
  virtual StatusCode initializeRefStates( LHCb::Track&        track,
                                          const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion() ) const = 0;
};
#endif
