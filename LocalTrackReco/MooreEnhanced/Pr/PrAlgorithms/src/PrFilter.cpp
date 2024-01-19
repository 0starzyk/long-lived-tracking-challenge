/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "PrAlgorithms/PrFilter.h"
#include "Event/Particle_v2.h"
#include "Event/PrFittedForwardTracks.h"
#include "Event/PrVeloTracks.h"
#include "Event/PrimaryVertices.h"
#include "PrKernel/PrSelection.h"

namespace Pr {

  // LHCb::Pr::Velo::Tracks -> LHCb::Pr::Velo::Tracks
  DECLARE_COMPONENT_WITH_ID( Filter<LHCb::Pr::Velo::Tracks>, "PrFilter__PrVeloTracks" )

  // LHCb::Pr::Fitted::Forward::Tracks -> LHCb::Pr::Fitted::Forward::Tracks
  DECLARE_COMPONENT_WITH_ID( Filter<LHCb::Pr::Fitted::Forward::Tracks>, "PrFilter__PrFittedForwardTracks" )

  // LHCb::Pr::Fitted::Forward::TracksWithPVs
  DECLARE_COMPONENT_WITH_ID( Filter<LHCb::Pr::Fitted::Forward::TracksWithPVs>,
                             "PrFilter__PrFittedForwardTracksWithPVs" )

  // LHCb::Pr::Fitted::Forward::TracksWithMuonID
  DECLARE_COMPONENT_WITH_ID( Filter<LHCb::Pr::Fitted::Forward::TracksWithMuonID>,
                             "PrFilter__PrFittedForwardTracksWithMuonID" )

  DECLARE_COMPONENT_WITH_ID( Filter<LHCb::Event::Composites>, "PrFilter__Composites" )

  DECLARE_COMPONENT_WITH_ID( Filter<LHCb::Event::PV::PrimaryVertexContainer>, "PrFilter__PV" )

} // namespace Pr
