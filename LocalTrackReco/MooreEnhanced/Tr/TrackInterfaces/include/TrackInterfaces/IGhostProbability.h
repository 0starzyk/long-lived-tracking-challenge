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
#include "GaudiKernel/IAlgTool.h"
#include <string_view>

/** @class IGhostProbability
 *
 *  interface for the ghost probability calculation
 *
 *  @author P.Seyfert
 *  @date   27/01/2015
 */

struct IGhostProbability : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IGhostProbability, 2, 0 );

  /** Add the ghost probability */
  virtual StatusCode execute( LHCb::Track& ) const = 0;

  /** reveal the variable names for a track type */
  virtual std::vector<std::string_view> variableNames( LHCb::Track::Types type ) const = 0;

  /** reveal the variable values for a track */
  virtual std::vector<float> netInputs( LHCb::Track& aTrack ) const = 0;
};
