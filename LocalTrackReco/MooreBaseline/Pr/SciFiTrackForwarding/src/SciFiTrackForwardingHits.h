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
#include "FTDAQ/FTInfo.h"
#include "Kernel/EventLocalAllocator.h"

#include <array>
#include <utility>
#include <vector>

// Custom span that allow for negative indices
struct SciFiTrackForwardingHits {
  using allocator_type = LHCb::Allocators::EventLocal<float>;
  using hits_t         = std::vector<float, allocator_type>;
  using IDs_t          = std::vector<unsigned, std::allocator_traits<allocator_type>::rebind_alloc<unsigned>>;
  std::array<std::pair<int, int>, LHCb::Detector::FT::NFTZones> zonerange{};
  hits_t                                                        hits;
  IDs_t                                                         IDs;
  SciFiTrackForwardingHits( allocator_type alloc = {} ) : hits{alloc}, IDs{alloc} {}
};
