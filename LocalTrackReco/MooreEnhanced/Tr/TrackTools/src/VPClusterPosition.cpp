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

#include "TrackInterfaces/IVPClusterPosition.h"

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/VPLightCluster.h"
#include "VPDet/DeVP.h"
#include "VPDet/DeVPSensor.h"

#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/Vector3DTypes.h"

/**
 *  @author Victor Coco
 *  @date   2010-02-02
 */
class VPClusterPosition : public extends<GaudiTool, IVPClusterPosition> {
public:
  using extends::extends;
  LHCb::VPPositionInfo position( const DeVP& det, const LHCb::VPLightCluster& cluster ) const override;

private:
  /// Average error for single-pixel clusters
  Gaudi::Property<double> m_errorSinglePixel{this, "ErrorSinglePixel", 12.5 * Gaudi::Units::micrometer};
  /// Average error for two-pixel clusters
  Gaudi::Property<double> m_errorTwoPixel{this, "ErrorTwoPixel", 8. * Gaudi::Units::micrometer};
  /// Average error (all cluster sizes)
  Gaudi::Property<double> m_errorAverage{this, "ErrorAverage", 12. * Gaudi::Units::micrometer};
};

DECLARE_COMPONENT( VPClusterPosition )

LHCb::VPPositionInfo VPClusterPosition::position( const DeVP& det, const LHCb::VPLightCluster& cluster ) const {
  // TODO: include track information, parameterise error as function of angle.
  LHCb::VPPositionInfo pos;
  // Get the position directly from the cluster.
  pos.x = cluster.x();
  pos.y = cluster.y();
  // Get the sensor.
  const LHCb::Detector::VPChannelID channel = cluster.channelID();
  const DeVPSensor&                 sensor  = det.sensor( channel.sensor() );
  // Initialise the local error estimate.
  double dx = m_errorAverage;
  double dy = m_errorAverage;
  // Get the inter-pixel fraction.
  const auto fx = cluster.xfraction();
  const auto fy = cluster.yfraction();
  // Get the cluster size.
  const unsigned int nPixels = 1;
  if ( nPixels == 1 ) {
    // Single-pixel clusters
    dx = m_errorSinglePixel;
    dy = m_errorSinglePixel;
    if ( sensor.isLong( channel ) ) dx *= 2.;
  } else if ( nPixels == 2 ) {
    // Two-pixel clusters
    if ( fx < 0.1 ) {
      dx = m_errorSinglePixel;
      if ( sensor.isLong( channel ) ) dx *= 2.;
    } else {
      dx = m_errorTwoPixel;
      dx *= sensor.xPitch()[channel.scol()] / sensor.xPitch()[0];
    }
    dy = fy < 0.1 ? m_errorSinglePixel : m_errorTwoPixel;
  } else if ( nPixels >= 9 ) {
    dx = 20. * Gaudi::Units::micrometer;
    dy = 20. * Gaudi::Units::micrometer;
  }
  // Transform the error estimate to the global frame.
  const auto sensorNumber = sensor.sensorNumber();
  pos.dx = sqrt( dx * dx * det.cos2OfSensor( sensorNumber ) + dy * dy * det.sin2OfSensor( sensorNumber ) );
  pos.dy = sqrt( dx * dx * det.sin2OfSensor( sensorNumber ) + dy * dy * det.cos2OfSensor( sensorNumber ) );
  return pos;
}
