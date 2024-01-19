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

// from Kernel/LHCbDefinitions
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "LHCbMath/Line.h"

// local
#include "Line.h"

namespace Tf::Tsa {

  /* @typedef for a line in 3-D
   *
   *  @author M.Needham
   *  @date   31/05/2004
   */
  using Line3D = Gaudi::Math::Line<Gaudi::XYZPoint, Gaudi::XYZVector>;

  /// Create a Line3D from a point line and z reference point
  inline Line3D createLine3D( const Tsa::Line& xLine, const Tsa::Line& yLine, const double zRef ) {
    return Line3D{Gaudi::XYZPoint{xLine.value( zRef ), yLine.value( zRef ), zRef},
                  Gaudi::XYZVector{xLine.m(), yLine.m(), 1.}.Unit()};
  }

} // namespace Tf::Tsa
