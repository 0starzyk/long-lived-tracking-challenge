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

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/State.h"
#include "Event/StateParameters.h"
#include "Event/StateVector.h"
#include "Event/Track.h"
#include "Kernel/IUTChannelIDSelector.h"
#include "LHCbMath/GeomFun.h"
#include "Line3D.h"
#include "TrackInterfaces/IHitExpectation.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackKernel/TrackFunctors.h"
#include "UTDet/DeUTDetector.h"
#include "UTDet/DeUTSensor.h"

#include "GaudiAlg/FunctionalTool.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IBinder.h"
#include "GaudiKernel/Plane3DTypes.h"
#include "GaudiKernel/SystemOfUnits.h"

#include <algorithm>
#include <string>

namespace {
  struct Cache {
    double zUTa;
    double zUTb;
  };
} // namespace

namespace LHCb {

  /**
   * Implementation of UTHitExpectation tool
   * see interface header for description
   *
   *  @author M.Needham
   *  @date   22/5/2007
   */
  class UTHitExpectation
      : public Gaudi::Functional::ToolBinder<
            Gaudi::Interface::Bind::Box<IHitExpectation>( DeUTDetector const& ),
            DetDesc::usesBaseAndConditions<DetDesc::ConditionAccessorHolder<GaudiTool>, DeUTDetector>> {

    class BoundInstance final : public Gaudi::Interface::Bind::Stub<IHitExpectation> {
      const UTHitExpectation* m_parent;
      const DeUTDetector&     m_detector;

    public:
      BoundInstance( UTHitExpectation const* parent, DeUTDetector const& det ) : m_parent{parent}, m_detector{det} {}
      unsigned int nExpected( const Track& aTrack, IGeometryInfo const& geometry ) const override {
        return m_parent->nExpected( aTrack, geometry, m_detector );
      }
      Info expectation( const Track& aTrack, IGeometryInfo const& geometry ) const override {
        return m_parent->expectation( aTrack, geometry, m_detector );
      }
      void collect( const Track& aTrack, std::vector<LHCbID>& ids, IGeometryInfo const& geometry ) const override {
        return m_parent->collect( aTrack, ids, geometry, m_detector );
      }
    };

  public:
    UTHitExpectation( std::string type, std::string name, const IInterface* parent )
        : ToolBinder{std::move( type ),
                     std::move( name ),
                     parent,
                     {KeyValue{"UTLocation", DeUTDetLocation::location()}},
                     construct<BoundInstance>( this )} {}

    /** intialize */
    StatusCode initialize() override;

    /** Returns number of hits expected, from zFirst to inf
     *
     *  @param aTrack Reference to the Track to test
     *
     *  @return number of hits expected
     */
    unsigned int nExpected( const Track& aTrack, IGeometryInfo const& geometry, DeUTDetector const& det ) const;

    /** Returns number of hits expected
     *
     *  @param aTrack Reference to the Track to test
     *
     *  @return Info info including likelihood
     */
    IHitExpectation::Info expectation( const Track& aTrack, IGeometryInfo const& geometry,
                                       DeUTDetector const& det ) const;

    /** Collect all the expected hits
     *
     * @param aTrack Reference to the Track to test
     * @param hits collected lhcbIDs
     *
     **/
    void collect( const Track& aTrack, std::vector<LHCbID>& ids, IGeometryInfo const& geometry,
                  DeUTDetector const& det ) const;

  private:
    void collectHits( std::vector<Detector::UT::ChannelID>& chans, StateVector stateVec, unsigned int station,
                      IGeometryInfo const& geometry, DeUTDetector const& det ) const;

    bool insideSensor( const DeUTSensor& sensor, const Tf::Tsa::Line3D& line ) const;

    Gaudi::XYZPoint intersection( const Tf::Tsa::Line3D& line, const Gaudi::Plane3D& aPlane ) const;

    ToolHandle<Gaudi::Interface::Bind::IBinder<IUTChannelIDSelector>> m_selector{this, "SelectorType",
                                                                                 "UTSelectChannelIDByElement"};
    PublicToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackParabolicExtrapolator"};

    ConditionAccessor<Cache> m_cache{this, name() + "_Cache"};
    Gaudi::Property<bool>    m_allStrips{this, "allStrips", false};
  };

  DECLARE_COMPONENT_WITH_ID( UTHitExpectation, "UTHitExpectation" )

} // namespace LHCb

bool LHCb::UTHitExpectation::insideSensor( const DeUTSensor& sensor, const Tf::Tsa::Line3D& line ) const {

  bool            isIn = false;
  Gaudi::XYZPoint point;
  double          mu;
  if ( Gaudi::Math::intersection( line, sensor.plane(), point, mu ) == true ) { isIn = sensor.globalInActive( point ); }
  return isIn;
}

StatusCode LHCb::UTHitExpectation::initialize() {
  return ToolBinder::initialize().andThen( [&]() -> StatusCode {
    addConditionDerivation( {DeUTDetLocation::location()}, m_cache.key(), [&]( DeUTDetector const& utDet ) -> Cache {
      if ( utDet.nStation() != 2u ) {
        throw GaudiException( "2 stations needed in UT", "UTHitExpectation::initialize", StatusCode::FAILURE );
      }
      return {utDet.station( 0 ).globalCentre().z(), utDet.station( 1 ).globalCentre().z()};
    } );
    return StatusCode::SUCCESS;
  } );
}

IHitExpectation::Info LHCb::UTHitExpectation::expectation( const LHCb::Track& aTrack, IGeometryInfo const& geometry,
                                                           DeUTDetector const& det ) const {

  IHitExpectation::Info info;
  info.likelihood = 0.0;
  info.nFound     = 0;
  info.nExpected  = nExpected( aTrack, geometry, det );
  return info;
}

unsigned int LHCb::UTHitExpectation::nExpected( const LHCb::Track& aTrack, IGeometryInfo const& geometry,
                                                DeUTDetector const& det ) const {

  // make a line at UTa and UTb
  auto&        cache    = m_cache.get();
  const State& UTaState = closestState( aTrack, cache.zUTa );
  StateVector  stateVectorUTa( UTaState.position(), UTaState.slopes() );

  const State& UTbState = closestState( aTrack, cache.zUTb );
  StateVector  stateVectorUTb( UTbState.position(), UTbState.slopes() );

  // determine which modules should be hit
  std::vector<Detector::UT::ChannelID> expectedHitsA;
  expectedHitsA.reserve( 4 );
  std::vector<Detector::UT::ChannelID> expectedHitsB;
  expectedHitsB.reserve( 4 );
  collectHits( expectedHitsA, stateVectorUTa, 1, geometry, det );
  collectHits( expectedHitsB, stateVectorUTb, 2, geometry, det );

  return expectedHitsA.size() + expectedHitsB.size();
}

void LHCb::UTHitExpectation::collect( const LHCb::Track& aTrack, std::vector<LHCb::LHCbID>& ids,
                                      IGeometryInfo const& geometry, DeUTDetector const& det ) const {
  // make a line at UTa and UTb
  auto&        cache    = m_cache.get();
  const State& UTaState = closestState( aTrack, cache.zUTa );
  StateVector  stateVectorUTa( UTaState.position(), UTaState.slopes() );

  const State& UTbState = closestState( aTrack, cache.zUTb );
  StateVector  stateVectorUTb( UTbState.position(), UTbState.slopes() );

  // determine which modules should be hit
  std::vector<Detector::UT::ChannelID> expectedHits;
  expectedHits.reserve( 8 );
  collectHits( expectedHits, stateVectorUTa, 1, geometry, det );
  collectHits( expectedHits, stateVectorUTb, 2, geometry, det );

  // convert to LHCb ids
  ids.reserve( expectedHits.size() );
  std::transform( expectedHits.begin(), expectedHits.end(), std::back_inserter( ids ),
                  []( Detector::UT::ChannelID chan ) { return LHCbID{chan}; } );
}

void LHCb::UTHitExpectation::collectHits( std::vector<LHCb::Detector::UT::ChannelID>& chans, LHCb::StateVector stateVec,
                                          unsigned int station, IGeometryInfo const& geometry,
                                          DeUTDetector const& det ) const {

  // loop over the sectors
  det.applyToAllSectors( [&]( DeUTSector const& sector ) {
    // propagate to z of sector
    m_extrapolator->propagate( stateVec, sector.globalCentre().z(), geometry )
        .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    const Detector::UT::ChannelID elemID = sector.elementID();
    if ( elemID.station() == station ) {
      // loop over sensors
      sector.applyToAllSensors( [&]( DeUTSensor const& sensor ) {
        Tf::Tsa::Line3D             aLine3D( stateVec.position(), stateVec.slopes() );
        auto                        box          = m_selector.bind( Gaudi::Hive::currentContext() );
        const IUTChannelIDSelector& selectorTool = box;
        if ( selectorTool.select( sector.elementID(), det ) && insideSensor( sensor, aLine3D ) == true ) {

          // get the list of strips that could be in the cluster
          Gaudi::XYZPoint globalEntry = intersection( aLine3D, sensor.entryPlane() );
          Gaudi::XYZPoint globalExit  = intersection( aLine3D, sensor.exitPlane() );
          Gaudi::XYZPoint localEntry  = sensor.toLocal( globalEntry );
          Gaudi::XYZPoint localExit   = sensor.toLocal( globalExit );

          unsigned int firstStrip = sensor.localUToStrip( localEntry.x() );
          unsigned int lastStrip  = sensor.localUToStrip( localExit.x() );

          // might have to swap...
          if ( firstStrip > lastStrip ) std::swap( firstStrip, lastStrip );

          // allow for capacitive coupling....
          if ( sensor.isStrip( firstStrip - 1 ) == true ) --firstStrip;
          if ( sensor.isStrip( lastStrip + 1 ) == true ) ++lastStrip;

          bool         found       = false;
          unsigned int middleStrip = ( firstStrip + lastStrip ) / 2;
          for ( unsigned int iStrip = firstStrip; iStrip != lastStrip; ++iStrip ) {
            const Detector::UT::ChannelID chan = sector.stripToChan( iStrip );
            if ( sector.isOKStrip( chan ) == true ) { // check it is alive
              if ( m_allStrips == true ) {
                chans.push_back( chan ); // take them all
              } else {
                found = true; // take just the sector
              }
            } // ok strip
          }   // loop strips

          if ( !m_allStrips && found == true ) {
            Detector::UT::ChannelID midChan( elemID.type(), elemID.station(), elemID.layer(), elemID.detRegion(),
                                             elemID.sector(), middleStrip );
            chans.push_back( midChan );
          }
        }
      } );
    }
  } );
}

Gaudi::XYZPoint LHCb::UTHitExpectation::intersection( const Tf::Tsa::Line3D& line,
                                                      const Gaudi::Plane3D&  aPlane ) const {
  // make a plane
  Gaudi::XYZPoint inter;
  double          mu = 0;
  Gaudi::Math::intersection( line, aPlane, inter, mu );
  return inter;
}
