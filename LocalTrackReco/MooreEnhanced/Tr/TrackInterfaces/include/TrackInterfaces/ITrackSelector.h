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
#ifndef _ITrackSelector_H
#define _ITrackSelector_H

#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

/** @class ITrackSelector
 *
 *  interface for selecting tracks....
 *
 *  @author M.Needham
 *  @date   31/05/2004
 */

struct ITrackSelector : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackSelector, 1, 0 );

  /// the method
  virtual bool accept( const LHCb::Track& aTrack ) const = 0;
};

#endif
