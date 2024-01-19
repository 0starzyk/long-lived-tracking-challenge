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

#include <array>
#include <tuple>
#include <vector>

// Gaudi
#include "LHCbAlgs/Transformer.h"

// LHCb
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Detector/VP/VPChannelID.h"
#include "Event/RawEvent.h"
#include "Event/VPFullCluster.h"
#include "Kernel/STLExtensions.h"
#include "VPDet/DeVP.h"
#include "VPKernel/PixelModule.h"
#include "VPKernel/VeloPixelInfo.h"

// FIXME still need to handle nx and ny for the clusters

// Namespace for locations in TDS
namespace LHCb {
  namespace VPFullClusterLocation {
    inline const std::string Offsets = "Raw/VP/FullClustersOffsets";
  }
} // namespace LHCb

namespace LHCb::Pr::Velo {

  class VPClusFull : public LHCb::Algorithm::MultiTransformer<
                         std::tuple<std::vector<VPFullCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>>(
                             const EventContext&, const RawEvent&, const DeVP& ),
                         LHCb::DetDesc::usesConditions<DeVP>> {

  public:
    /// Standard constructor
    VPClusFull( const std::string& name, ISvcLocator* pSvcLocator );

    /// Algorithm execution
    std::tuple<std::vector<VPFullCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>>
    operator()( const EventContext&, const RawEvent&, const DeVP& ) const override;

  private:
    /// Maximum allowed cluster size (no effect when running on lite clusters).
    unsigned int m_maxClusterSize = 999999999;

    std::bitset<VP::NModules>                  m_modulesToSkipMask;
    Gaudi::Property<std::vector<unsigned int>> m_modulesToSkip{this,
                                                               "ModulesToSkip",
                                                               {},
                                                               [=]( auto& ) {
                                                                 m_modulesToSkipMask.reset();
                                                                 for ( auto i : m_modulesToSkip ) {
                                                                   m_modulesToSkipMask.set( i );
                                                                 }
                                                               },
                                                               Gaudi::Details::Property::ImmediatelyInvokeHandler{true},
                                                               "List of modules that should be skipped in decoding"};
  };
} // namespace LHCb::Pr::Velo
