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
#ifndef KERNEL_IVPCLUSTERPOSITION_H
#define KERNEL_IVPCLUSTERPOSITION_H 1

#include "VPDet/DeVP.h"

#include "Kernel/PixelPositionInfo.h"

#include "GaudiKernel/IAlgTool.h"

/**
 *  @author Victor Coco
 *  @date   2010-02-02
 */

namespace LHCb {
  class VPLightCluster;
}

struct IVPClusterPosition : extend_interfaces<IAlgTool> {
  DeclareInterfaceID( IVPClusterPosition, 3, 0 );

  /** Calculate position of a given VPCluster
   * @return struct containing coordinates and errors
   * The returned error estimate depends both on the pixel size and
   * the projected angle of a track.
   */
  virtual LHCb::VPPositionInfo position( const DeVP& det, const LHCb::VPLightCluster& cluster ) const = 0;
};
#endif // KERNEL_IVPCLUSTERPOSITION_H
