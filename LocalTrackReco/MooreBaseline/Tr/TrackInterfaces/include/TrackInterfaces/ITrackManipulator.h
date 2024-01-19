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
#ifndef _ITrackManipulator_H
#define _ITrackManipulator_H

#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

/** @class ITrackManipulator
 *
 *  interface for setting the track reference parameters from the measurement info
 *
 *  @author M.Needham
 *  @date   16/06/2006
 */

struct ITrackManipulator : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackManipulator, 1, 0 );

  /** Add the reference information */
  virtual StatusCode execute( LHCb::Track& aTrack ) const = 0;
};

#endif
