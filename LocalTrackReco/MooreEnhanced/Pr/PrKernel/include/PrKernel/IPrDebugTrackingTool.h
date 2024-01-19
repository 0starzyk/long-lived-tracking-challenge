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
#pragma once

#include "GaudiKernel/IAlgTool.h"

#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "Kernel/STLExtensions.h"

struct IPrDebugTrackingTool : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IPrDebugTrackingTool, 1, 0 );

  using VariableDef =
      std::pair<std::string_view,
                std::variant<int, float, double, std::vector<int>, std::vector<float>, std::vector<double>>>;

  virtual int check( int = -1, int = -1, const std::vector<int>& = {} ) const = 0;

  virtual void storeData( LHCb::span<const VariableDef>, std::string_view = "Tuple" ) const = 0;
};
