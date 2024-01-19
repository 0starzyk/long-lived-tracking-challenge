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
/** @class TrackListFilter TrackListFilter.h
 *
 *  Algorithm to filter events in which a track list is not empty
 *
 *  @author W. Hulsbergen
 *  @date   2008
 */
#include "Event/Track.h"
#include "LHCbAlgs/FilterPredicate.h"

namespace LHCb {
  struct TrackListFilter final : Algorithm::FilterPredicate<bool( const LHCb::Track::Range& )> {
    TrackListFilter( const std::string& name, ISvcLocator* pSvcLocator )
        : FilterPredicate{name, pSvcLocator, {"inputLocation", LHCb::TrackLocation::Default}} {}
    bool operator()( const LHCb::Track::Range& r ) const override { return !r.empty(); }
  };

  DECLARE_COMPONENT_WITH_ID( TrackListFilter, "TrackListFilter" )
} // namespace LHCb
