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
#include <cstdint>
#include <vector>

///
std::vector<uint8_t> SerializeStates( const std::vector<LHCb::State>& states );

///
std::vector<LHCb::State> DeserializeStates( const std::vector<uint8_t>& binary );

///
struct StateVectorDifference {
  intptr_t numStates = 0;
  double   avgPos    = 0;
  double   maxPos    = 0;
  double   avgCov    = 0;
  double   maxCov    = 0;
};

///
StateVectorDifference CompareStates( const std::vector<LHCb::State>& lhs, const std::vector<LHCb::State>& rhs );
