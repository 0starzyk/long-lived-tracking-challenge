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
#ifndef _ITrackFunctor_h
#define _ITrackFunctor_h

#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

/** @class ITrackFunctor
 *
 *  interface to compute some float given a Track
 */

struct ITrackFunctor : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackFunctor, 1, 0 );

  virtual float operator()( const LHCb::Track& aTrack ) const = 0;
};

#endif
