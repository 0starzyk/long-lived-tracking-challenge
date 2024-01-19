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
#include <Event/TrackParameters.h>
#include <Event/UTCluster.h>
#include <Event/UTLiteCluster.h>
#include <GaudiAlg/GaudiTool.h>
#include <GaudiKernel/DataObjectHandle.h>
#include <GaudiKernel/ToolHandle.h>
#include <Kernel/LineTraj.h>
#include <PrKernel/UTHitHandler.h>
#include <TrackInterfaces/IUTClusterPosition.h>
#include <UTDet/DeUTDetector.h>
#include <UTDet/DeUTSector.h>
#include <UTDet/DeUTSensor.h>
#include <mutex>
#include <type_traits>

class UTMeasurementProvider final : public MeasurementProviderProjector {
public:
  using MeasurementProviderProjector::MeasurementProviderProjector;

  void addToMeasurements( LHCb::span<LHCb::LHCbID> ids, std::vector<LHCb::Measurement>& measurements,
                          LHCb::ZTrajectory<double> const& reftraj ) const override;

  StatusCode load( LHCb::Track& ) const override {
    info() << "sorry, MeasurementProviderBase::load not implemented" << endmsg;
    return StatusCode::FAILURE;
  }

  void reset() override {}

private:
  using Cluster              = UT::Hit;
  using ClusterContainerType = UT::HitHandler;

  void addToMeasurements( DeUTDetector const& det, LHCb::span<LHCb::LHCbID> ids,
                          std::vector<LHCb::Measurement>& measurements ) const;

  LHCb::Measurement measurement( DeUTDetector const& det, Cluster const& cluster ) const;

  DataObjectReadHandle<ClusterContainerType> m_clustersDH{this, "ClusterLocation", UTInfo::HitLocation};

  Gaudi::Property<bool>          m_useReference{this, "UseReference", true, "[[deprecated]] ignored"};
  Gaudi::Property<bool>          m_use2D{this, "Use2D", false, "[[deprecated]] ignored"};
  ToolHandle<IUTClusterPosition> m_positiontool = {"UTOnlinePosition/UTLiteClusterPosition"};

  ConditionAccessor<DeUTDetector> m_ut{this, "DEUT", DeUTDetLocation::location()};
};

LHCb::Measurement UTMeasurementProvider::measurement( DeUTDetector const&                   det,
                                                      UTMeasurementProvider::Cluster const& clus ) const {
  auto utSector = det.findSector( clus.chanID() );
  return {clus.lhcbID(), utSector->globalCentre().z(), utSector->trajectory( clus.chanID(), clus.fracStrip() ),
          m_positiontool->error( clus.pseudoSize() ) * utSector->pitch(),
#ifdef USE_DD4HEP
          *utSector};
#else
          utSector};
#endif
}

void UTMeasurementProvider::addToMeasurements( LHCb::span<LHCb::LHCbID>        ids,
                                               std::vector<LHCb::Measurement>& measurements,
                                               LHCb::ZTrajectory<double> const& ) const {
  addToMeasurements( m_ut.get(), ids, measurements );
}

void UTMeasurementProvider::addToMeasurements( DeUTDetector const& det, LHCb::span<LHCb::LHCbID> ids,
                                               std::vector<LHCb::Measurement>& measurements ) const {
  measurements.reserve( measurements.size() + ids.size() );
  auto to_clus = [clusters = m_clustersDH.get()]( LHCb::LHCbID id ) { return clusters->hit( id.utID() ); };
  std::transform( ids.begin(), ids.end(), std::back_inserter( measurements ), [&]( LHCb::LHCbID id ) {
    auto clus = to_clus( id );
    assert( clus != nullptr );
    return this->measurement( det, *clus );
  } );
}

DECLARE_COMPONENT( UTMeasurementProvider )
