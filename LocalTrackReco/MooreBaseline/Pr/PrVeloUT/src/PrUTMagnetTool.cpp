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
#include "PrUTMagnetTool.h"
#include <Event/State.h>
#include <GaudiKernel/Point3DTypes.h>
#include <GaudiKernel/SystemOfUnits.h>
#include <GaudiKernel/ThreadLocalContext.h>
#include <GaudiKernel/Vector3DTypes.h>
#include <Kernel/LUTForFunction.h>
#include <Kernel/STLExtensions.h>

//-----------------------------------------------------------------------------
// Implementation file for class : PrUTMagnetTool
//
// 2006-09-25 : Mariusz Witek
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_COMPONENT_WITH_ID( LHCb::Pr::UTMagnetTool, "PrUTMagnetTool" )

namespace LHCb::Pr {
  namespace {
    struct bdl_t {
      float BdlTrack, zHalfBdlTrack;
    };

    bdl_t f_bdl( const DeMagnet& magnet, float slopeY, float zOrigin, float zStart, float zStop ) {

      if ( zStart > zStop ) return {0, 0};

      float Bdl      = 0.0;
      float zHalfBdl = 0.0;

      Gaudi::XYZPoint  aPoint( 0., 0., 0. );
      Gaudi::XYZVector bField;

      constexpr int np = 500;
      float         dz = ( zStop - zStart ) / np;
      float         dy = dz * slopeY;
      float         z  = zStart + dz / 2.;
      float         y  = slopeY * ( zStart - zOrigin );

      // vectors to calculate z of half Bdl
      std::vector<float> bdlTmp, zTmp;
      bdlTmp.reserve( np + 1 );
      zTmp.reserve( np + 1 );
      while ( z < zStop ) {

        aPoint.SetY( y );
        aPoint.SetZ( z );

        bField = magnet.fieldVector( aPoint );
        Bdl += dy * bField.z() - dz * bField.y();
        if ( z > 100. * Gaudi::Units::cm ) {
          bdlTmp.push_back( Bdl );
          zTmp.push_back( z );
        }
        z += dz;
        y += dy;
      }

      float bdlhalf = std::abs( Bdl ) / 2.;

      for ( unsigned int i = 5; i < bdlTmp.size() - 5; i++ ) {
        if ( std::abs( bdlTmp[i] ) > bdlhalf ) {
          float zrat = ( Bdl / 2. - bdlTmp[i - 1] ) / ( bdlTmp[i] - bdlTmp[i - 1] );
          zHalfBdl   = zTmp[i - 1] + dz * zrat;
          break;
        }
      }

      return {Bdl, zHalfBdl};
    }

  } // namespace

  const UTMagnetTool::Cache& UTMagnetTool::cache() const {
    return m_cache.get( getConditionContext( Gaudi::Hive::currentContext() ) );
  }

  StatusCode UTMagnetTool::initialize() {
    return ConditionAccessorHolder::initialize().andThen( [&]() {
      addConditionDerivation( {LHCb::standard_geometry_top, LHCb::Det::Magnet::det_path, DeUTDetLocation::location()},
                              m_cache.key(),
                              [&]( const LHCb::Detector::DeLHCb& lhcb, const DeMagnet& magnet,
                                   const DeUTDetector& utdet ) { return makeCache( lhcb, magnet, utdet ); } );
    } );
  }

  //=========================================================================
  // Callback function for cache creation
  //=========================================================================
  UTMagnetTool::Cache UTMagnetTool::makeCache( const LHCb::Detector::DeLHCb& lhcb, const DeMagnet& magnet,
                                               const DeUTDetector& utdet ) const {
    UTMagnetTool::Cache c;

    prepareBdlTables( magnet, utdet, c );

    prepareDeflectionTables( *lhcb.geometry(), magnet, c );

    // check whether B=0
    auto bdl  = f_bdl( magnet, 0., 0., 400. * Gaudi::Units::mm, c.zCenterUT );
    c.noField = ( std::abs( bdl.BdlTrack ) < 10e-4 );

    if ( c.noField ) {
      info() << " No B field detected." << endmsg;
      // override computed values with hardcoded constants
      c.dist2mom  = s_averageDist2mom_NoB;
      c.zMidField = s_zMidField_NoB;
    } else {
      c.zMidField = c.lutZHalfBdl->getInterpolatedValue( {0.05, 0.0, 0.0} );
    }

    return c;
  }

  //=========================================================================
  // prepareBdlTables
  //=========================================================================
  void UTMagnetTool::prepareBdlTables( const DeMagnet& magnet, const DeUTDetector& utdet,
                                       UTMagnetTool::Cache& c ) const {

    if ( msgLevel( MSG::DEBUG ) ) debug() << "Start generation of VeloUT Bdl LUTs" << endmsg;
    // prepare table with Bdl integrations
    // Bdl integral depends on 3 track parameters
    //  slopeY     - y slope of the track
    //  zOrigin    - z of track intersection with z axis (in YZ projection)
    //  zVeloEnd   - z of the track at which slopeY is given
    //                      slopeY    zOrigin    zVeloEnd
    // m_zCenterUT is a normalization plane which should be close to middle of UT.
    // It is used to normalize dx deflection at different UT layers.
    // No need to update with small UT movement up to +- 5 cm.

    c.zCenterUT     = 2484.6;
    float zCenterUT = 0.;

    assert( utdet.nLayer() == 4 );
    unsigned il = 0;
    utdet.applyToAllLayers( [&il, &zCenterUT, &c]( DeUTLayer const& layer ) {
      float zlay      = layer.sector( 0 ).globalCentre().Z();
      c.zLayers[il++] = zlay;
      zCenterUT += zlay;
    } );
    zCenterUT /= c.zLayers.size();

    if ( std::abs( c.zCenterUT - zCenterUT ) > 50. ) {
      warning() << "Calculated center of UT station far away from nominal value: " << zCenterUT << " wrt nominal "
                << c.zCenterUT << endmsg;
      warning() << " Calculated value taken: " << zCenterUT << endmsg;
      c.zCenterUT = zCenterUT;
    }
    // warning: layers a-priori not in order of increasing z!
    std::sort( c.zLayers.begin(), c.zLayers.end() );

    c.lutBdl->fillTable( [&]( LHCb::span<const float, 3> var ) {
      return f_bdl( magnet, var[0], var[1], var[2], c.zCenterUT ).BdlTrack;
    } );

    c.lutZHalfBdl->fillTable( [&]( LHCb::span<const float, 3> var ) {
      return f_bdl( magnet, var[0], var[1], var[2], c.zCenterUT ).zHalfBdlTrack;
    } );

    if ( msgLevel( MSG::DEBUG ) ) debug() << "Generation of VeloUT Bdl LUTs finished" << endmsg;
  }

  //=========================================================================
  // prepareDeflectionTables
  //=========================================================================

  void UTMagnetTool::prepareDeflectionTables( IGeometryInfo const& geometry, const DeMagnet& magnet,
                                              UTMagnetTool::Cache& c ) const {

    if ( msgLevel( MSG::DEBUG ) ) debug() << "Start generation of VeloUT deflection LUTs" << endmsg;

    // tmp state
    LHCb::State tmpState;
    float       qpBeg = 1. / ( 10. * Gaudi::Units::GeV );
    tmpState.setState( 0., 0., 0., 0., 0., qpBeg );
    // set dummy covariance matrix
    Gaudi::TrackSymMatrix cov = Gaudi::TrackSymMatrix();
    cov( 0, 0 )               = 0.1;
    cov( 1, 1 )               = 0.1;
    cov( 2, 2 )               = 0.1;
    cov( 3, 3 )               = 0.1;
    cov( 4, 4 )               = 0.1;
    tmpState.setCovariance( cov );

    // determine normalization factors for deflections in different UT layers
    // wrt center of UT
    auto dxLay = [&]( LHCb::span<const float, 2> lutVar ) {
      auto idLay   = int( lutVar[0] + 0.000001 );
      auto dydzBeg = lutVar[1];
      auto zLay    = c.zLayers[idLay];

      if ( msgLevel( MSG::DEBUG ) ) {
        debug() << lutVar[0] << " " << lutVar[1] << " "
                << " idlay " << idLay << " " << c.zLayers.size() << " zlay" << zLay << " " << endmsg;
      }

      tmpState.setState( 0., 0., 0., 0., dydzBeg, qpBeg );

      LHCb::State stalin_mid = tmpState;
      LHCb::State stapar_mid = tmpState;
      // extrapolate state to middle of UT
      // extrapolate state to middle of UT
      StatusCode sc1 =
          m_linear->propagate( stalin_mid, c.zCenterUT, geometry, LHCb::Tr::PID::Pion(), magnet.fieldGrid() );
      StatusCode sc2 =
          m_parabolic->propagate( stapar_mid, c.zCenterUT, geometry, LHCb::Tr::PID::Pion(), magnet.fieldGrid() );

      if ( sc1.isFailure() || sc2.isFailure() ) { Warning( "Extrapolation failed ", StatusCode::SUCCESS, 0 ).ignore(); }

      float       ratio      = 0.;
      LHCb::State stalin_lay = tmpState;
      LHCb::State stapar_lay = tmpState;

      StatusCode sc3 = m_linear->propagate( stalin_lay, zLay, geometry, LHCb::Tr::PID::Pion(), magnet.fieldGrid() );
      StatusCode sc4 = m_parabolic->propagate( stapar_lay, zLay, geometry, LHCb::Tr::PID::Pion(), magnet.fieldGrid() );

      if ( sc3.isFailure() || sc4.isFailure() ) {
        Warning( "Extrapolation failed ", StatusCode::SUCCESS, 0 ).ignore();
      } else {
        auto dx_mid = stapar_mid.x() - stalin_mid.x();
        auto dx_lay = stapar_lay.x() - stalin_lay.x();
        if ( std::abs( dx_mid ) > 1.e-8 ) ratio = dx_mid / dx_lay;
      }
      return ratio;
    };

    c.lutDxLay->fillTable( dxLay );

    // distance to momentum table (depends on y slope only)
    auto dxToMom = [&]( LHCb::span<const float, 1> var ) {
      tmpState.setState( 0., 0., 0., 0., var[0], qpBeg );
      LHCb::State stalin_mid = tmpState;
      LHCb::State stapar_mid = tmpState;
      // extrapolate state to middle of UT
      StatusCode sc1 =
          m_linear->propagate( stalin_mid, c.zCenterUT, geometry, LHCb::Tr::PID::Pion(), magnet.fieldGrid() );
      StatusCode sc2 =
          m_parabolic->propagate( stapar_mid, c.zCenterUT, geometry, LHCb::Tr::PID::Pion(), magnet.fieldGrid() );

      float dx2mom = 0.;
      if ( sc1.isFailure() || sc2.isFailure() ) {
        Warning( "Extrapolation failed ", StatusCode::SUCCESS, 0 ).ignore();
      } else {
        float dx = stapar_mid.x() - stalin_mid.x();
        if ( std::abs( dx ) > 1e-8 ) dx2mom = qpBeg / dx;
      }
      return dx2mom;
    };

    c.lutDxToMom->fillTable( dxToMom );

    // determine distToMomentum parameter
    c.dist2mom = c.lutDxToMom->getValue( {0.05} );

    if ( msgLevel( MSG::DEBUG ) ) debug() << "Generation of VeloUT deflection LUTs finished" << endmsg;
  }
} // namespace LHCb::Pr
//****************************************************************************
