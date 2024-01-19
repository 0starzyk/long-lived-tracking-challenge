/*****************************************************************************\
* (c) Copyright 2020 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once

#include "DetDesc/ITransportSvc.h"

namespace LHCb {
  /**
   * See https://gitlab.cern.ch/lhcb/Rec/-/issues/153 about the numerical precision
   **/
  struct TrackMaterialIntersection final {
    double      z1 = 0.;
    double      z2 = 0.;
    double      tx = 0.;
    double      ty = 0.;
    MaterialPtr material{};
  };
} // namespace LHCb
