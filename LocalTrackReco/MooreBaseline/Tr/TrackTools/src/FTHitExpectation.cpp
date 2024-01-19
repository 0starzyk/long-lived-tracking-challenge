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
/**
 * Implementation of FTHitExpectation tool
 * see interface header for description
 *
 *  @author D.Milanes
 *  @date   20/13/2013
 */

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/State.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "FTDet/DeFTDetector.h"
#include "LHCbMath/GeomFun.h"
#include "TrackInterfaces/IHitExpectation.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackKernel/TrackFunctors.h"

#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IMagneticFieldSvc.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/SystemOfUnits.h"

#include "Line.h"
#include "Line3D.h"
#include "Parabola.h"

#include "GaudiAlg/FunctionalTool.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IBinder.h"

#include <algorithm>
#include <utility>

namespace LHCb {

  class FTHitExpectation : public Gaudi::Functional::ToolBinder<
                               Gaudi::Interface::Bind::Box<IHitExpectation>( DeFT const& ),
                               DetDesc::usesBaseAndConditions<DetDesc::ConditionAccessorHolder<GaudiTool>, DeFT>> {

    class BoundInstance final : public Gaudi::Interface::Bind::Stub<IHitExpectation> {
      const FTHitExpectation* m_parent;
      const DeFT&             m_detector;

    public:
      BoundInstance( FTHitExpectation const* parent, DeFT const& det ) : m_parent{parent}, m_detector{det} {}
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
    FTHitExpectation( std::string type, std::string name, const IInterface* parent )
        : ToolBinder{std::move( type ),
                     std::move( name ),
                     parent,
                     {KeyValue{"FTLocation", DeFTDetectorLocation::Default}},
                     construct<BoundInstance>( this )} {}

    /** Returns number of hits expected, from zFirst to inf
     *
     *  @param aTrack Reference to the Track to test
     *
     *  @return number of hits expected
     */
    unsigned int nExpected( const Track& aTrack, IGeometryInfo const& geometry, DeFT const& det ) const;

    /** Returns number of hits expected
     *
     *  @param aTrack Reference to the Track to test
     *
     *  @return Info info including likelihood
     */
    IHitExpectation::Info expectation( const Track& aTrack, IGeometryInfo const& geometry, DeFT const& det ) const;

    /** Collect all the expected hits
     *
     * @param aTrack Reference to the Track to test
     * @param hits collected lhcbIDs
     *
     **/
    void collect( const Track& aTrack, std::vector<LHCbID>& ids, IGeometryInfo const& geometry, DeFT const& det ) const;

  private:
    bool findIntersectingMat( const DeFTLayer& layer, const Track& aTrack, const DeFTMat*& mat,
                              Gaudi::XYZPoint& intersectionPoint, IGeometryInfo const& geometry ) const;

    Gaudi::XYZPoint intersection( const Tf::Tsa::Line3D& line, const Gaudi::Plane3D& aPlane ) const;

    double            curvature( const State& aState ) const;
    Tf::Tsa::Parabola xParabola( const Track& aTrack, const double z, IGeometryInfo const& geometry ) const;
    Tf::Tsa::Line     yLine( const Track& aTrack, const double z, IGeometryInfo const& geometry ) const;

    ServiceHandle<IMagneticFieldSvc> m_pIMF{this, "MagneticFieldSvc", "MagneticFieldSvc"};
    ToolHandle<ITrackExtrapolator>   m_extrapolator{this, "extrapolatorName", "TrackParabolicExtrapolator"};
  };

  DECLARE_COMPONENT_WITH_ID( FTHitExpectation, "FTHitExpectation" )

} // namespace LHCb

IHitExpectation::Info LHCb::FTHitExpectation::expectation( const LHCb::Track& aTrack, IGeometryInfo const& geometry,
                                                           DeFT const& det ) const {

  IHitExpectation::Info info;
  info.likelihood = 0.0;
  info.nFound     = 0;
  info.nExpected  = 0;

  const auto&         ids = aTrack.lhcbIDs();
  std::vector<LHCbID> ftHits;
  ftHits.reserve( ids.size() );
  std::copy_if( ids.begin(), ids.end(), std::back_inserter( ftHits ), []( const LHCbID& id ) { return id.isFT(); } );

  // Loop over all layers
  auto func = [this, &aTrack, &geometry, &info, &ftHits]( const DeFTLayer& layer ) {
    const DeFTMat*  mat{nullptr};
    Gaudi::XYZPoint intersectionPoint;
    bool            expectHit = findIntersectingMat( layer, aTrack, mat, intersectionPoint, geometry );

    if ( expectHit ) {
      ++( info.nExpected );

      // Check of this layer is present in the hits
      auto itIter = std::find_if( ftHits.begin(), ftHits.end(), [&mat]( const LHCbID& id ) {
        return id.ftID().globalLayerID() == mat->elementID().globalLayerID();
      } );
      if ( itIter != ftHits.end() ) ++info.nFound;
    }
  };
  det.applyToAllLayers( func );

  return info;
}

bool LHCb::FTHitExpectation::findIntersectingMat( const DeFTLayer& layer, const LHCb::Track& aTrack,
                                                  const DeFTMat*& mat, Gaudi::XYZPoint& intersectionPoint,
                                                  IGeometryInfo const& geometry ) const {

  // make a Line3D from the track
  double            layerZ  = layer.globalZ();
  Tf::Tsa::Line     line    = yLine( aTrack, layerZ, geometry );
  Tf::Tsa::Parabola aParab  = xParabola( aTrack, layerZ, geometry );
  Tf::Tsa::Line     tanLine = aParab.tangent( layerZ );
  Tf::Tsa::Line3D   aLine3D = Tf::Tsa::createLine3D( tanLine, line, layerZ );

  // find intersection point of track and plane of layer
  const auto layerPlane = layer.plane();
  intersectionPoint     = intersection( aLine3D, layerPlane );

  // find the module that is hit
  const auto module = layer.findModule( intersectionPoint );
  if ( !module ) return false;

  // find intersection point of track and plane of module
  // (to account for misalignments between module and layer)
  tanLine                = aParab.tangent( intersectionPoint.z() );
  aLine3D                = Tf::Tsa::createLine3D( tanLine, line, intersectionPoint.z() );
  const auto modulePlane = module->plane();
  intersectionPoint      = intersection( aLine3D, modulePlane );

  // Find the corresponding mat
  const auto foundMat = module->findMat( intersectionPoint );
  // check if intersection point is inside the fibres of the module
  if ( !foundMat ) return false;
  //  mat = static_cast<DeFTMat*>(foundMat);//FIXME
  mat = &( *foundMat ); // FIXME
  return foundMat->isInside( intersectionPoint );
}

void LHCb::FTHitExpectation::collect( const LHCb::Track& aTrack, std::vector<LHCb::LHCbID>& ids,
                                      IGeometryInfo const& geometry, DeFT const& det ) const {
  // Loop over all layers
  // FIXME
  auto func = [this, &aTrack, &geometry, &ids]( const DeFTLayer& layer ) {
    const DeFTMat*  mat{nullptr};
    Gaudi::XYZPoint intersectionPoint;
    bool            expectHit = findIntersectingMat( layer, aTrack, mat, intersectionPoint, geometry );

    if ( expectHit ) {
      // Find the channel that is closest
      Gaudi::XYZPoint localP    = mat->toLocal( intersectionPoint );
      const auto [ftChan, frac] = mat->calculateChannelAndFrac( localP.x() );

      // Add the channels
      // JvT: Without the fraction the bare FTChannelID is pretty useless...
      if ( std::abs( frac ) <= 0.5f ) ids.push_back( LHCbID( ftChan ) );
    }
  };
  det.applyToAllLayers( func );
}

unsigned int LHCb::FTHitExpectation::nExpected( const LHCb::Track& aTrack, IGeometryInfo const& geometry,
                                                DeFT const& det ) const {
  return expectation( aTrack, geometry, det ).nExpected;
}

Gaudi::XYZPoint LHCb::FTHitExpectation::intersection( const Tf::Tsa::Line3D& line,
                                                      const Gaudi::Plane3D&  aPlane ) const {
  Gaudi::XYZPoint inter;
  double          mu = 0;
  Gaudi::Math::intersection( line, aPlane, inter, mu );
  return inter;
}

Tf::Tsa::Parabola LHCb::FTHitExpectation::xParabola( const LHCb::Track& aTrack, const double z,
                                                     IGeometryInfo const& geometry ) const {
  // find the closest state
  const State& aState = closestState( aTrack, z );
  StateVector  stateVector( aState.position(), aState.slopes() );
  m_extrapolator->propagate( stateVector, StateParameters::ZMidT, geometry )
      .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

  // make a parabola from this...
  const double zS = aState.z();
  const double a  = curvature( aState );
  const double b  = aState.tx() - 2 * zS * a;
  const double c  = aState.x() - zS * ( b + a * zS );
  return Tf::Tsa::Parabola( a, b, c );
}

double LHCb::FTHitExpectation::curvature( const LHCb::State& aState ) const {
  Gaudi::XYZPoint  P = aState.position();
  Gaudi::XYZVector B;
  m_pIMF->fieldVector( P, B ).ignore();
  const double tx   = aState.tx();
  const double ty   = aState.ty();
  auto         tmp  = 1.0 + std::pow( tx, 2 );
  const double nTx  = sqrt( tmp );
  const double norm = sqrt( tmp + std::pow( ty, 2 ) );

  return -0.5 * norm * nTx * Gaudi::Units::c_light * B.y() * aState.qOverP();
}

Tf::Tsa::Line LHCb::FTHitExpectation::yLine( const LHCb::Track& aTrack, const double z,
                                             IGeometryInfo const& geometry ) const {
  // find the closest state
  State       aState = closestState( aTrack, z );
  StateVector stateVector( aState.position(), aState.slopes() );
  m_extrapolator->propagate( stateVector, StateParameters::ZMidT, geometry )
      .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

  const double m = aState.ty();
  const double c = aState.y() - m * aState.z();
  return Tf::Tsa::Line( m, c );
}
