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
#include "GaudiKernel/Point3DTypes.h"
#include "Kernel/STLExtensions.h"

#include <vector>

#include "DetDesc/IGeometryInfo.h"
namespace LHCb {
  class RecVertex;
}

struct IPVFitter : extend_interfaces<IAlgTool> {
  DeclareInterfaceID( IPVFitter, 4, 0 );
  virtual StatusCode fitVertex( const Gaudi::XYZPoint& seedPoint, LHCb::span<const LHCb::Track* const> tracks,
                                LHCb::RecVertex& vtx, std::vector<const LHCb::Track*>& tracks2remove,
                                IGeometryInfo const& geometry ) const = 0;
};
