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
#include "Event/Track.h"
#include <map>
#include <string>

namespace TrackMonitorMaps {

  struct HistoRange final {

    HistoRange( std::string id, double xMin, double xMax ) : fid( std::move( id ) ), fxMin( xMin ), fxMax( xMax ) {}
    std::string fid;
    double      fxMin;
    double      fxMax;
  };

  using InfoHistMap = std::map<LHCb::Track::AdditionalInfo, HistoRange>;
  const InfoHistMap& infoHistDescription();

} // namespace TrackMonitorMaps
