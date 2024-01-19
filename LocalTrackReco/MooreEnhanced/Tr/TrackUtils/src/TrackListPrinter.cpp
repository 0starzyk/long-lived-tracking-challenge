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
#include "Event/Track.h"
#include "LHCbAlgs/Consumer.h"

/** @class TrackListPrinter TrackListPrinter.h
 *
 *  @author M.Needham
 *  @date   30/05/2006
 */

struct TrackListPrinter final : LHCb::Algorithm::Consumer<void( const LHCb::Track::Range& )> {
  TrackListPrinter( const std::string& name, ISvcLocator* pSvc )
      : Consumer{name, pSvc, {"InputLocation", LHCb::TrackLocation::Default}} {}

  void operator()( const LHCb::Track::Range& in ) const override {
    for ( const auto& t : in ) info() << *t << endmsg;
  }
};

DECLARE_COMPONENT( TrackListPrinter )
