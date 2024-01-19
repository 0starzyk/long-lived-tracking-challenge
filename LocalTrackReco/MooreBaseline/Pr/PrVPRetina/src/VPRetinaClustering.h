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
#include "Event/VPLightCluster.h"
#include "VPDet/DeVP.h"

// Rec
#include "PrKernel/PrPixelModule.h"
#include "PrKernel/VeloPixelInfo.h"

// local
#include "VPRetinaMatrix.h"

// FIXME still need to handle nx and ny for the clusters

// Namespace for locations in TDS
namespace LHCb {
  namespace VPClusterLocation {
    inline const std::string Offsets = "Raw/VP/LightClustersOffsets";
  }
} // namespace LHCb

class VPRetinaClustering
    : public LHCb::Algorithm::MultiTransformer<
          std::tuple<std::vector<LHCb::VPLightCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>>(
              const LHCb::RawEvent&, const DeVP& ),
          LHCb::DetDesc::usesConditions<DeVP>> {

public:
  /// Standard constructor
  VPRetinaClustering( const std::string& name, ISvcLocator* pSvcLocator );

  /// Algorithm execution
  std::tuple<std::vector<LHCb::VPLightCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>>
  operator()( const LHCb::RawEvent&, const DeVP& ) const override;

private:
  /// make RetinaClusters from bank
  std::vector<uint32_t> makeRetinaClusters( const uint32_t* bank ) const;

  std::bitset<VP::NModules>                  m_modulesToSkipMask;
  Gaudi::Property<std::vector<unsigned int>> m_modulesToSkip{this,
                                                             "ModulesToSkip",
                                                             {},
                                                             [=]( auto& ) {
                                                               // if( msgLevel(MSG::DEBUG)){
                                                               //   info() << "Modules to skip size: " <<
                                                               //   m_modulesToSkip.size() << endmsg;
                                                               // }
                                                               m_modulesToSkipMask.reset();
                                                               for ( auto i : m_modulesToSkip ) {
                                                                 // if( msgLevel(MSG::DEBUG)){
                                                                 //   info() << "Skipping module " << moduleNumber <<
                                                                 //   endmsg;
                                                                 // }
                                                                 m_modulesToSkipMask.set( i );
                                                               }
                                                             },
                                                             Gaudi::Details::Property::ImmediatelyInvokeHandler{true},
                                                             "List of modules that should be skipped in decoding"};
};
