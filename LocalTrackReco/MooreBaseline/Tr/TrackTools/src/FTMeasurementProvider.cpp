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
#include <Event/FTLiteCluster.h>
#include <FTDet/DeFTDetector.h>

#include <GaudiAlg/GaudiTool.h>
#include <GaudiKernel/DataObjectHandle.h>

/** @class FTMeasurementProvider FTMeasurementProvider.cpp
 *
 * Implementation of FTMeasurementProvider
 * see interface header for description
 *
 *  @author Wouter Hulsbergen
 *  @date   30/12/2005
 */
class FTMeasurementProvider final : public MeasurementProviderProjector {
public:
  using MeasurementProviderProjector::MeasurementProviderProjector;

  void addToMeasurements( LHCb::span<LHCb::LHCbID> ids, std::vector<LHCb::Measurement>& measurements,
                          const LHCb::ZTrajectory<double>& reftraj ) const override;

  StatusCode load( LHCb::Track& ) const override {
    return Error( "sorry, MeasurementProviderBase::load not implemented" );
  }

private:
  ConditionAccessor<DeFT>                                   m_det{this, DeFTDetectorLocation::Default};
  DataObjectReadHandle<LHCb::FTLiteCluster::FTLiteClusters> m_clustersDh{this, "ClusterLocation",
                                                                         LHCb::FTLiteClusterLocation::Default};
};

//=============================================================================
// Declare to tool factory
//=============================================================================

DECLARE_COMPONENT( FTMeasurementProvider )

//-----------------------------------------------------------------------------
/// Create measurements for list of LHCbIDs
//-----------------------------------------------------------------------------

void FTMeasurementProvider::addToMeasurements( LHCb::span<LHCb::LHCbID>        ids,
                                               std::vector<LHCb::Measurement>& measurements,
                                               const LHCb::ZTrajectory<double>& ) const {
  const auto& det = m_det.get();

  measurements.reserve( measurements.size() + ids.size() );
  assert( std::all_of( ids.begin(), ids.end(), []( LHCb::LHCbID id ) { return id.isFT(); } ) );
  std::transform(
      ids.begin(), ids.end(), std::back_inserter( measurements ),
      [&, clusters = m_clustersDh.get()]( const LHCb::LHCbID& id ) {
        /// The clusters are not sorted anymore, so we can use a find_if
        /// to find the element corresponding to the channel ID
        const auto  ftMat = det.findMat( id.ftID() );
        const auto& c     = clusters->range( id.ftID().globalQuarterIdx() );

        auto itH =
            id.isFT()
                ? std::find_if( c.begin(), c.end(),
                                [&id]( const LHCb::FTLiteCluster clus ) { return clus.channelID() == id.ftID(); } )
                : c.end();
        if ( itH == c.end() ) {
          throw GaudiException( "Can not find FTLiteCluster for given lhcbID", __func__, StatusCode::FAILURE );
        }
#ifdef USE_DD4HEP
        return LHCb::Measurement{itH->channelID(), ftMat->globalZ(),
                                 ftMat->trajectory( itH->channelID(), itH->fraction() ),
                                 0.04 + 0.01 * itH->pseudoSize(), // FIXME need a better error parametrization
                                 *ftMat};
#else
        return LHCb::Measurement{itH->channelID(), ftMat->globalZ(),
                                 ftMat->trajectory( itH->channelID(), itH->fraction() ),
                                 0.04 + 0.01 * itH->pseudoSize(), // FIXME need a better error parametrization
                                 ftMat};
#endif
      } );
}
