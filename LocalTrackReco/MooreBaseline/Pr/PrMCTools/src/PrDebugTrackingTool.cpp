/*****************************************************************************\
* (c) Copyright 2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "PrKernel/PrDebugTrackingToolBase.h"

#include <limits>
#include <vector>

namespace LHCb::Pr {
  struct PrDebugTrackingTool : public PrDebugTrackingToolBase {

    // inherit standard constructors
    using PrDebugTrackingToolBase::PrDebugTrackingToolBase;
    // dummy implementation because not needed at the moment
    int check( int = -1, int = -1, const std::vector<int>& = {} ) const override {
      return std::numeric_limits<int>::lowest();
    };
  };

  // Declaration of the Tool Factory
  DECLARE_COMPONENT_WITH_ID( PrDebugTrackingTool, "PrDebugTrackingTool" )

} // namespace LHCb::Pr
