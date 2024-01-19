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
#include <string>

#include "DetDesc/IGeometryInfo.h"

/** @class IPVSeeding IPVSeeding.h newtool/IPVSeeding.h
 *
 *
 *  @author Mariusz Witek
 *  @date   2008-05-19
 */
struct IPVSeeding : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IPVSeeding, 4, 0 );

  virtual std::vector<Gaudi::XYZPoint> getSeeds( LHCb::span<const LHCb::Track* const> inputTracks,
                                                 const Gaudi::XYZPoint&               beamspot,
                                                 IGeometryInfo const&                 geometry ) const = 0;
};
