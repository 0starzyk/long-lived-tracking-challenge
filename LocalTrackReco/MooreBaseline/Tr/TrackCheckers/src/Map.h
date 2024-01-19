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
#ifndef _TrackMaps_H
#define _TrackMaps_H

#include <map>
#include <string>

#include "Event/Track.h"
#include "MCInterfaces/IMCReconstructible.h"

namespace TrackMaps {

  typedef std::map<std::string, LHCb::Track::Types> TypeMap;
  const TypeMap&                                    typeDescription();

  typedef std::map<std::string, LHCb::Track::AdditionalInfo> InfoMap;
  // const InfoMap& infoDescription() ;

  typedef std::map<std::string, IMCReconstructible::RecCategory> RecMap;
  const RecMap&                                                  recDescription();

  class HistoRange {
  public:
    HistoRange( std::string id, double xMin, double xMax ) : fid( id ), fxMin( xMin ), fxMax( xMax ) { ; }
    std::string fid;
    double      fxMin;
    double      fxMax;
  };

  typedef std::map<LHCb::Track::AdditionalInfo, HistoRange> InfoHistMap;
  const InfoHistMap&                                        infoHistDescription();

} // namespace TrackMaps

#endif
