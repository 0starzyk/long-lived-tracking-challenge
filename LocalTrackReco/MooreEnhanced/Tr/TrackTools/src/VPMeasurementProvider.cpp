/*****************************************************************************\
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "MeasurementProviderProjector.h"
#include <DetDesc/GenericConditionAccessorHolder.h>
#include <Event/TrackParameters.h>
#include <Event/VPCluster.h>
#include <Event/VPLightCluster.h>
#include <GaudiAlg/GaudiTool.h>
#include <GaudiKernel/DataObjectHandle.h>
#include <GaudiKernel/ToolHandle.h>
#include <Kernel/LineTraj.h>
#include <TrackInterfaces/IVPClusterPosition.h>
#include <VPDet/DeVP.h>
#include <mutex>
#include <type_traits>

class VPMeasurementProvider final : public MeasurementProviderProjector {
public:
  using MeasurementProviderProjector::MeasurementProviderProjector;

  void addToMeasurements( LHCb::span<LHCb::LHCbID> ids, std::vector<LHCb::Measurement>& measurements,
                          LHCb::ZTrajectory<double> const& reftraj ) const override;

  StatusCode load( LHCb::Track& ) const override {
    error() << "sorry, MeasurementProviderBase::load not implemented" << endmsg;
    return StatusCode::FAILURE;
  }

  void reset() override {}

private:
  using Cluster              = LHCb::VPLightCluster;
  using ClusterContainerType = LHCb::VPLightClusters;

  void addToMeasurements( DeVP const& det, LHCb::span<LHCb::LHCbID> ids,
                          std::vector<LHCb::Measurement>& measurements ) const;

  LHCb::Measurement measurement( DeVP const& det, Cluster const& cluster,
                                 LHCb::Measurement::VP::Projection proj ) const;

  DataObjectReadHandle<ClusterContainerType> m_clustersDH{this, "ClusterLocation", LHCb::VPClusterLocation::Light};

  Gaudi::Property<bool>          m_useReference{this, "UseReference", true, "[[deprecated]] ignored"};
  Gaudi::Property<bool>          m_use2D{this, "Use2D", false}; // use 2D measurements for VP
  ToolHandle<IVPClusterPosition> m_positiontool = {"VPClusterPosition"};

  ConditionAccessor<DeVP> m_vp{this, "DEVP", LHCb::Det::VP::det_path};
};

LHCb::Measurement VPMeasurementProvider::measurement( DeVP const& det, Cluster const& clus,
                                                      LHCb::Measurement::VP::Projection proj ) const {
#ifdef USE_DD4HEP
  const auto& sensor = det.sensor( clus.channelID().sensor() );
#else
  const auto* sensor = &det.sensor( clus.channelID().sensor() );
#endif

  LHCb::VPPositionInfo info = m_positiontool->position( det, clus );
  Gaudi::XYZPoint      position( info.x, info.y, clus.z() );

  if ( proj == LHCb::Measurement::VP::Projection::XY ) {
    // 2D measurement
    std::array<double, 2> r = {info.dx, info.dy};
    return LHCb::Measurement{clus.channelID(), clus.z(), position, ROOT::Math::SVector<double, 2>{r.begin(), r.end()},
                             sensor};
  } else {
    if ( proj == LHCb::Measurement::VP::Projection::Y ) {
      return LHCb::Measurement{clus.channelID(), clus.z(),
                               LHCb::LineTraj<double>{position, LHCb::Trajectory<double>::Vector{1, 0, 0},
                                                      LHCb::Trajectory<double>::Range{-info.dx, info.dx},
                                                      LHCb::Trajectory<double>::DirNormalized{true}},
                               info.dy, sensor};
    } else {
      return LHCb::Measurement{clus.channelID(), clus.z(),
                               LHCb::LineTraj<double>{position, LHCb::Trajectory<double>::Vector{0, 1, 0},
                                                      LHCb::Trajectory<double>::Range{-info.dy, info.dy},
                                                      LHCb::Trajectory<double>::DirNormalized{true}},
                               info.dx, sensor};
    }
  }
}
void VPMeasurementProvider::addToMeasurements( LHCb::span<LHCb::LHCbID>        ids,
                                               std::vector<LHCb::Measurement>& measurements,
                                               LHCb::ZTrajectory<double> const& ) const {
  addToMeasurements( m_vp.get(), ids, measurements );
}

void VPMeasurementProvider::addToMeasurements( DeVP const& det, LHCb::span<LHCb::LHCbID> ids,
                                               std::vector<LHCb::Measurement>& measurements ) const {
  measurements.reserve( measurements.size() + ( m_use2D ? 1 : 2 ) * ids.size() );
  auto to_clus = [clusters = m_clustersDH.get()]( LHCb::LHCbID id ) {
    const auto clus = id.isVP() ? std::find_if( clusters->begin(), clusters->end(),
                                                [cid = id.vpID().channelID()]( LHCb::VPLightCluster const& clus ) {
                                                  return clus.channelID() == cid;
                                                } )
                                : clusters->end();
    return clus != clusters->end() ? &*clus : nullptr;
  };
  std::for_each( ids.begin(), ids.end(), [&]( LHCb::LHCbID id ) {
    auto clus = to_clus( id );
    assert( clus != nullptr );
    if ( m_use2D ) {
      measurements.push_back( measurement( det, *clus, LHCb::Measurement::VP::Projection::XY ) );
    } else {
      measurements.push_back( measurement( det, *clus, LHCb::Measurement::VP::Projection::X ) );
      measurements.push_back( measurement( det, *clus, LHCb::Measurement::VP::Projection::Y ) );
    }
  } );
}

DECLARE_COMPONENT( VPMeasurementProvider )
