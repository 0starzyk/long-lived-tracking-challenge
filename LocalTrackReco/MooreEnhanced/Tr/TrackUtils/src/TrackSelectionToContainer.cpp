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
/** @class TrackSelectionToContainer TrackSelectionToContainer.h
 *
 *  Convert LHCb::Tracks to LHCb::Track::Selection
 */

#include "Event/Track.h"
#include "LHCbAlgs/Transformer.h"
#include <string>

struct TrackSelectionToContainer final : LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Track::Selection& )> {
  TrackSelectionToContainer( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator, {"InputLocation", {}}, {"OutputLocation", {}} ) {}

  LHCb::Tracks operator()( const LHCb::Track::Selection& tracks ) const override {
    LHCb::Tracks out;
    for ( const auto& track : tracks ) { out.insert( new LHCb::Track( *track ) ); }
    return out;
  }
};

DECLARE_COMPONENT( TrackSelectionToContainer )
