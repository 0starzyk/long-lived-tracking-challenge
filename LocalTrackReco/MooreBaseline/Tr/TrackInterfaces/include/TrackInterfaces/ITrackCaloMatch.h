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
#ifndef _ITrackCaloMatch_H
#define _ITrackCaloMatch_H

#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

/** @class ITrackCaloMatch
 *
 *  interface for getting energy deposited in calos associated to track
 *  returned value is the appropriately weighted sum of ecal, hcal and preshower
 *  zero indicates no energy found
 *
 *  @author M.Needham
 *  @date   31/05/2005
 */

struct ITrackCaloMatch : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackCaloMatch, 1, 0 );

  /// the method
  virtual double energy( const LHCb::Track& aTrack ) const = 0;
};

#endif
