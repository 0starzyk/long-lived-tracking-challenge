/*****************************************************************************\
* (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "Map.h"

const TrackMaps::TypeMap& TrackMaps::typeDescription() {
  static const TrackMaps::TypeMap s_map = {
      {"Velo", LHCb::Track::Types::Velo},         {"Long", LHCb::Track::Types::Long},
      {"Upstream", LHCb::Track::Types::Upstream}, {"Downstream", LHCb::Track::Types::Downstream},
      {"Ttrack", LHCb::Track::Types::Ttrack},     {"Muon", LHCb::Track::Types::Muon}};
  return s_map;
}

// const TrackMaps::InfoMap& TrackMaps::infoDescription()
//{
//   static TrackMaps::InfoMap s_map ;
//   if ( s_map.empty() ) {
//     s_map = boost::assign::map_list_of("Likelihood",LHCb::Track::Likelihood)
//                                       ("PatQuality",LHCb::Track::AdditionalInfo::PatQuality)
//                                       ("MatchChi2",LHCb::Track::AdditionalInfo::MatchChi2);
//   };
//   return s_map ;
//}

const TrackMaps::RecMap& TrackMaps::recDescription() {
  static const TrackMaps::RecMap s_map = {{"ChargedLong", IMCReconstructible::ChargedLong},
                                          {"ChargedDownstream", IMCReconstructible::ChargedDownstream},
                                          {"ChargedUpstream", IMCReconstructible::ChargedUpstream},
                                          {"ChargedTtrack", IMCReconstructible::ChargedTtrack},
                                          {"ChargedVelo", IMCReconstructible::ChargedVelo}};
  return s_map;
}

const TrackMaps::InfoHistMap& TrackMaps::infoHistDescription() {
  static const InfoHistMap s_map = {
      {LHCb::Track::AdditionalInfo::FitTChi2, HistoRange( "FitTChi2", 0., 100. )},
      {LHCb::Track::AdditionalInfo::FitTNDoF, HistoRange( "FitTNDof", 0., 50. )},
      {LHCb::Track::AdditionalInfo::FitVeloChi2, HistoRange( "FitVeloChi2", 0., 100. )},
      {LHCb::Track::AdditionalInfo::FitVeloNDoF, HistoRange( "FitVeloNDoF", 0., 50. )},
      {LHCb::Track::AdditionalInfo::FitMatchChi2, HistoRange( "FitMatchChi2", 0., 100. )}};
  return s_map;
}
