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
#ifndef PVUTILS_H
#define PVUTILS_H 1

#include "Event/RecVertex.h"
#include "Event/State.h"
#include "Event/Track.h"

#include <vector>

class PVTrack final {
public:
  PVTrack() = default;

  const LHCb::Track* refTrack = nullptr;
  // Current state of the track at the current point
  LHCb::State stateG;
  // Normalized vector of slope
  Gaudi::XYZVector unitVect;
  // Flag if the track has been used in a previous vertex
  bool isUsed = false;

  // Result for impact parameter
  Gaudi::XYZVector   vd0;        // Impact parameter vector
  double             d0sq   = 0; // Impact parameter squared
  double             err2d0 = 0; // IP error squared
  double             chi2   = 0; // chi2 = d02 / d0err**2
  double             weight = 0; // Weight assigned to track
  LHCb::Track::Types type   = LHCb::Track::Types::Velo;

  double zClose() const {
    return stateG.z() - unitVect.z() * ( unitVect.x() * stateG.x() + unitVect.y() * stateG.y() ) /
                            ( 1.0 - std::pow( unitVect.z(), 2 ) );
  }
};

typedef std::vector<PVTrack>  PVTracks;
typedef std::vector<PVTrack*> PVTrackPtrs;

struct PVVertex final {
  PVTrackPtrs     pvTracks;
  LHCb::RecVertex primVtx{0};
};

typedef std::vector<PVVertex> PVVertices;

#endif // PVUTILS_H
