/*****************************************************************************\
 * * (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
 * *                                                                             *
 * * This software is distributed under the terms of the GNU General Public      *
 * * Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
 * *                                                                             *
 * * In applying this licence, CERN does not waive the privileges and immunities *
 * * granted to it by virtue of its status as an Intergovernmental Organization  *
 * * or submit itself to any jurisdiction.                                       *
 * \*****************************************************************************/

//-----------------------------------------------------------------------------

/** @file VPExpectation.cpp
 *
 *  Implementation file for reconstruction tool : VPHitExpectation
 *
 *  @author Bhagyashree Pagare bhagyashree.pagare@cern.ch
 *  @date   22/03/2022
 */
//-----------------------------------------------------------------------------

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/State.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "Line3D.h"
#include "VPDet/DeVP.h"

#include "Event/Track.h"
#include "Event/TrackParameters.h"
#include "GaudiAlg/GaudiTool.h"
#include "Kernel/HitPattern.h"
#include "Kernel/LHCbID.h"
#include "TrackInterfaces/IDetailedHitExpectation.h"
#include "TrackInterfaces/IHitExpectation.h"

#include "GaudiAlg/FunctionalTool.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IBinder.h"

#include <bitset>
#include <string>

namespace LHCb {

  class VPHitExpectation : public Gaudi::Functional::ToolBinder<
                               Gaudi::Interface::Bind::Box<IDetailedHitExpectation>( DeVP const& ),
                               // Gaudi::Interface::Bind::Box<IHitExpectation, IDetailedHitExpectation>( DeVP const& ),
                               DetDesc::usesBaseAndConditions<DetDesc::ConditionAccessorHolder<GaudiTool>, DeVP>> {

    class BoundInstance final
        //: public Gaudi::Interface::Bind::Stub<IHitExpeactation, IDetailedHitExpectation> {
        : public Gaudi::Interface::Bind::Stub<IDetailedHitExpectation> {
      const VPHitExpectation* m_parent;
      const DeVP&             m_detector;

    public:
      BoundInstance( VPHitExpectation const* parent, DeVP const& det ) : m_parent{parent}, m_detector{det} {}
      // FIXME BoundInstnace and VPHitExpectation should also implement the IHitExpectation interface
      // However the Binder mecanism is not supporting it right now and this interface is
      // currently not used anywhere. Thus it has been commented for now
      /*
      unsigned int nExpected( const Track& aTrack, IGeometryInfo const& geometry ) const override {
        return m_parent->nExpected( aTrack, geometry, m_detector );
      }
      Info expectation( const Track& aTrack, IGeometryInfo const& geometry ) const override {
        return m_parent->expectation( aTrack, geometry, m_detector );
      }
      void collect( const Track& aTrack, std::vector<LHCbID>& ids, IGeometryInfo const& geometry ) const override {
        return m_parent->collect( aTrack, ids, geometry, m_detector );
      }
      */
      virtual DetailedInfo detailedExpectation( const Track& aTrack, IGeometryInfo const& geometry ) const override {
        return m_parent->detailedExpectation( aTrack, geometry, m_detector );
      }
    };

  public:
    VPHitExpectation( std::string type, std::string name, const IInterface* parent )
        : ToolBinder{std::move( type ),
                     std::move( name ),
                     parent,
                     {KeyValue{"VPDetectorLocation", DeVPLocation::Default}},
                     construct<BoundInstance>( this )} {}

    /**
     * Returns number of hits expected
     *
     * @param aTrack Reference to the Track to test
     *
     * @return number of hits expected
     */
    unsigned int nExpected( const Track& aTrack, IGeometryInfo const& geometry, const DeVP& det ) const;

    /** Returns number of hits expected
     *
     *  @param aTrack Reference to the Track to test
     *
     *  @return Info nexpected
     */

    IHitExpectation::Info expectation( const Track& aTrack, IGeometryInfo const& geometry, const DeVP& det ) const;

    /** Returns detailed info about the missed hits
     *
     *  @param aTrack Reference to the Track to test
     *
     *  @return Info
     */
    IDetailedHitExpectation::DetailedInfo detailedExpectation( const Track& aTrack, IGeometryInfo const& geometry,
                                                               const DeVP& det ) const;

    /** Collects the hits that are missed during the reconstruction
     *
     * @param aTrack:Reference to the Track to test and ids:reference to the track ids
     */
    void collect( const Track& aTrack, std::vector<LHCbID>& ids, IGeometryInfo const& geometry, const DeVP& det ) const;

  private:
    bool   isInside( const Track& aTrack, const DeVPSensor& sensor, const double zStart, const double zStop ) const;
    double zMin( const Track& aTrack, const DeVP& det ) const;
    double zMax( const Track& aTrack, const DeVP& det ) const;
    Detector::VPChannelID findIntersectingChannelID( const Track& aTrack, const DeVPSensor& sensor ) const;
  };

  DECLARE_COMPONENT_WITH_ID( VPHitExpectation, "VPHitExpectation" )

} // namespace LHCb

unsigned int LHCb::VPHitExpectation::nExpected( const Track& aTrack, IGeometryInfo const& geometry,
                                                const DeVP& det ) const {
  return expectation( aTrack, geometry, det ).nExpected;
}

IDetailedHitExpectation::DetailedInfo LHCb::VPHitExpectation::detailedExpectation( const Track&         aTrack,
                                                                                   IGeometryInfo const& geometry,
                                                                                   const DeVP&          det ) const {
  const auto geom = geometry;
  // work out the first and last z on the track
  const auto zStart = zMin( aTrack, det ) - 1e-3; // 1e-3 is chosen to ensure that even if the sensor is tilted we
                                                  // are counting the z cordinate of it
  const auto                            zStop = zMax( aTrack, det ) + 1e-3;
  IDetailedHitExpectation::DetailedInfo info;
  const auto&                           trkIDs = aTrack.lhcbIDs();
  std::vector<LHCbID>                   collectedIDs;
  collectedIDs.reserve( trkIDs.size() );
  info.reserve( 208 );
  std::copy_if( begin( trkIDs ), end( trkIDs ), back_inserter( collectedIDs ),
                []( auto i ) -> bool { return i.isVP(); } );

  det.runOnAllSensors( [&]( const DeVPSensor& sensor ) {
    if ( isInside( aTrack, sensor, zStart, zStop ) ) {
      const Detector::VPChannelID     channelID    = findIntersectingChannelID( aTrack, sensor );
      Detector::VPChannelID::SensorID sensorNumber = sensor.sensorNumber();
      const unsigned                  chipNumber   = unsigned( channelID.chip() );
      if ( chipNumber > 2 ) { Warning( format( "Chip number was %i", chipNumber ) ).ignore(); }
      const Detector::VPChannelID chanid( sensorNumber, channelID.chip(),
                                          static_cast<Detector::VPChannelID::ColumnID>( 0 ),
                                          static_cast<Detector::VPChannelID::RowID>( 0 ) );
      const LHCbID                lhcbid( chanid );
      const unsigned int          moduleNumber = sensor.module();

      info.emplace_back( lhcbid, moduleNumber, false );
      if ( std::find_if( begin( collectedIDs ), end( collectedIDs ), [sensorNumber, chipNumber]( auto i ) -> bool {
             return ( sensorNumber == i.vpID().sensor() ) && ( chipNumber == unsigned( i.vpID().chip() ) );
           } ) != end( collectedIDs ) ) {
        info.back().found = true;
      }
    }
  } );

  return info;
}

IHitExpectation::Info LHCb::VPHitExpectation::expectation( const Track& aTrack, IGeometryInfo const& geometry,
                                                           const DeVP& det ) const {
  auto                  detailed = detailedExpectation( aTrack, geometry, det );
  IHitExpectation::Info nHits;
  nHits.nExpected = detailed.size();
  nHits.nFound    = std::count_if( detailed.begin(), detailed.end(), []( const auto& ele ) { return ele.found; } );
  return nHits;
}

bool LHCb::VPHitExpectation::isInside( const LHCb::Track& aTrack, const DeVPSensor& sensor, const double zStart,
                                       const double zStop ) const {
  const double z = sensor.z();
  if ( z <= zStart || z >= zStop ) return false;
  const State&    state = aTrack.closestState( z );
  const auto      xLine = Tf::Tsa::Line( state.tx(), state.x(), state.z() );
  const auto      yLine = Tf::Tsa::Line( state.ty(), state.y(), state.z() );
  Gaudi::XYZPoint trackPos( xLine.value( z ), yLine.value( z ), z );
  return sensor.isInsideSensor( trackPos.x(), trackPos.y() ); // check if the point is in Active area of the sensor
}

double LHCb::VPHitExpectation::zMin( const Track& aTrack, const DeVP& det ) const {
  // get the hit at least z
  auto z = std::numeric_limits<double>::max();
  for ( auto id : aTrack.lhcbIDs() ) {
    if ( id.isVP() ) {
      const DeVPSensor& sensor = det.sensor( id.vpID() );
      z                        = std::min( z, sensor.z() );
    }
  } // loop ids
  return z;
}

void LHCb::VPHitExpectation::collect( const LHCb::Track& aTrack, std::vector<LHCb::LHCbID>& ids,
                                      IGeometryInfo const& geometry, const DeVP& det ) const {
  const auto          zStart = zMin( aTrack, det ) - 1e-3;
  const auto          zStop  = zMax( aTrack, det ) + 1e-3;
  const auto          geom   = geometry;
  std::vector<LHCbID> collectedIDs;
  const auto&         trkIDs = ids;
  std::copy_if( begin( trkIDs ), end( trkIDs ), back_inserter( collectedIDs ),
                []( auto i ) -> bool { return i.isVP(); } );
  det.runOnAllSensors( [&]( const DeVPSensor& sensor ) {
    if ( isInside( aTrack, sensor, zStart, zStop ) ) {
      Detector::VPChannelID::SensorID sensorNumber = sensor.sensorNumber();
      if ( std::find_if( begin( collectedIDs ), end( collectedIDs ), [sensorNumber]( auto i ) -> bool {
             return ( sensorNumber == i.vpID().sensor() );
           } ) == end( collectedIDs ) ) {
        const double z = sensor.z();
        if ( z >= zStart && z <= zStop ) {
          const Detector::VPChannelID channelID = findIntersectingChannelID( aTrack, sensor );
          collectedIDs.emplace_back( channelID );
        }
      }
    }
  } );
}

LHCb::Detector::VPChannelID LHCb::VPHitExpectation::findIntersectingChannelID( const Track&      aTrack,
                                                                               const DeVPSensor& sensor ) const {
  const double          z     = sensor.z();
  const State&          state = aTrack.closestState( z );
  const Tf::Tsa::Line   xLine( state.tx(), state.x(), state.z() );
  const Tf::Tsa::Line   yLine( state.ty(), state.y(), state.z() );
  const Gaudi::XYZPoint point( xLine.value( z ), yLine.value( z ), z );
  const bool            local = false;
  Detector::VPChannelID channelID;
  sensor.pointToChannel( point, local, channelID );
  return channelID;
}

double LHCb::VPHitExpectation::zMax( const Track& aTrack, const DeVP& det ) const {
  auto z = std::numeric_limits<double>::lowest();
  for ( auto id : aTrack.lhcbIDs() ) {
    if ( id.isVP() ) {
      const DeVPSensor& sensor = det.sensor( id.vpID() );
      z                        = std::max( z, sensor.z() );
    }
  } // loop ids
  return z;
}
