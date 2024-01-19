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
#include "Map.h"

namespace {
  using namespace TrackMonitorMaps;
  static const InfoHistMap s_map2 = {{LHCb::Track::AdditionalInfo::FitTChi2, HistoRange( "8", 0., 100. )},
                                     {LHCb::Track::AdditionalInfo::FitTNDoF, HistoRange( "9", 0., 50. )},
                                     {LHCb::Track::AdditionalInfo::FitVeloChi2, HistoRange( "10", 0., 100. )},
                                     {LHCb::Track::AdditionalInfo::FitVeloNDoF, HistoRange( "11", 0., 50. )},
                                     {LHCb::Track::AdditionalInfo::FitMatchChi2, HistoRange( "12", 0., 100. )}};
} // namespace

const TrackMonitorMaps::InfoHistMap& TrackMonitorMaps::infoHistDescription() { return s_map2; }
