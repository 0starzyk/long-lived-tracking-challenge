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
#pragma once

#include "Event/PrLongTracks.h"
#include "Event/State.h"
#include "GaudiKernel/IAlgTool.h"
#include "PrKernel/PrMutUTHits.h"
#include "boost/container/static_vector.hpp"

/** @class IPrAddUTHitsTool IPrAddUTHitsTool.h TrackInterfaces/IPrAddUTHitsTool.h
 *
 *  @author:  Michel De Cian
 *  @date:    13-11-2013
 */

class IPrAddUTHitsTool : public extend_interfaces<IAlgTool> {
  using Tracks = LHCb::Pr::Long::Tracks;

public:
  DeclareInterfaceID( IPrAddUTHitsTool, 2, 0 );

  /// Add UT clusters to Long tracks
  virtual void addUTHits( Tracks& tracks ) const = 0;
  virtual void
  getUTHits( const LHCb::State&                                                              state,
             boost::container::static_vector<LHCb::LHCbID, LHCb::Pr::TracksInfo::MaxUTHits>& utHits ) const = 0;
};
