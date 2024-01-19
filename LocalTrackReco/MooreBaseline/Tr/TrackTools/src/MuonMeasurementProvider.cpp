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
#include <Kernel/LineTraj.h>
#include <Kernel/STLExtensions.h>
#include <LHCbMath/LHCbMath.h>
#include <MuonDet/DeMuonDetector.h>
#include <MuonDet/MuonNamespace.h>
#include <memory>

namespace {
  // to calculate the coordinate of the cluster
  struct pos_t {
    double pos;
    double dpos;
  };
  pos_t recomputePos( const std::vector<std::array<double, 2>>& data ) {
    int    np  = 0;
    double sum = 0., sum2 = 0.;
    for ( auto i = data.begin(); i != data.end(); ++i ) {
      // check that this position is not already the same of a previous pad
      if ( std::none_of( data.begin(), i, [ref = ( *i )[0], step = 0.5 * ( *i )[1]]( const std::array<double, 2>& x ) {
             return std::abs( x[0] - ref ) < step;
           } ) ) {
        np++;
        sum += ( *i )[0];
        sum2 += std::pow( ( *i )[0], 2 );
      }
    }
    double av = sum / np;
    return {av, np > 1 ? sqrt( ( sum2 - sum * av ) / ( np - 1 ) ) : data.back()[1]};
  }

  std::array<LHCb::Measurement, 2> makeMeasurements( LHCb::LHCbID id, const DeMuonDetector& det, Gaudi::XYZPoint p,
                                                     double dx, double dy ) {
    LHCb::Detector::Muon::TileID muid = id.muonID();
    // they have promised to fix the const
#ifdef USE_DD4HEP
    const DeMuonChamber chamber = det.getChamberFromTile( muid );
#else
    const DeMuonChamber* chamber =
        &det.getChamber( muid.station(), muid.region(), det.Tile2FirstChamberNumber( muid ) );
#endif
    return {LHCb::Measurement( id, p.z(),
                               LHCb::LineTraj<double>{p, Gaudi::XYZVector{0, 1, 0}, std::pair{-dy, dy},
                                                      LHCb::Trajectory<double>::DirNormalized{true}},
                               2. * dx * LHCb::Math::inv_sqrt_12, chamber ),
            LHCb::Measurement( id, p.z(),
                               LHCb::LineTraj<double>{p, Gaudi::XYZVector{1, 0, 0}, std::pair{-dx, dx},
                                                      LHCb::Trajectory<double>::DirNormalized{true}},
                               2. * dy * LHCb::Math::inv_sqrt_12, chamber )};
  }
} // namespace

/** @class MuonMeasurementProvider MuonMeasurementProvider.cpp
 *
 *  @author Wouter Hulsbergen + Stefania Vecchi
 *  @date   30/12/2005            18/12/2009
 */
class MuonMeasurementProvider final : public MeasurementProviderProjector {
public:
  using MeasurementProviderProjector::MeasurementProviderProjector;

  void addToMeasurements( LHCb::span<LHCb::LHCbID> ids, std::vector<LHCb::Measurement>& measurements,
                          const LHCb::ZTrajectory<double>& reftraj ) const override;

  StatusCode load( LHCb::Track& ) const override {
    info() << "sorry, MeasurementProviderBase::load not implemented" << endmsg;
    return StatusCode::FAILURE;
  }

private:
  /// measurement for single hits
  std::array<LHCb::Measurement, 2> measurement( const DeMuonDetector& det, const LHCb::LHCbID& id ) const;
  /// measurement for cluster of hits
  std::array<LHCb::Measurement, 2> measurement( const DeMuonDetector& det, LHCb::span<LHCb::LHCbID> ids ) const;

  // pointer to detector
  ConditionAccessor<DeMuonDetector> m_det{this, DeMuonLocation::Default};
  Gaudi::Property<bool>             m_clusterize{this, "clusterize", false};
};

//=============================================================================
// Declare to tool factory
//=============================================================================

DECLARE_COMPONENT( MuonMeasurementProvider )

//-----------------------------------------------------------------------------
/// Create a measurement
//-----------------------------------------------------------------------------

std::array<LHCb::Measurement, 2> MuonMeasurementProvider::measurement( const DeMuonDetector& det,
                                                                       const LHCb::LHCbID&   id ) const {
  if ( !id.isMuon() ) {
    error() << "Not a Muon measurement" << endmsg;
    throw GaudiException( "LHCbID provided is not a muon ID", __func__, StatusCode::FAILURE );
  }

  auto pos = det.position( id.muonID() );
  if ( !pos ) {
    Warning( "Failed to get x,y,z of tile " ).ignore();
    if ( msgLevel( MSG::DEBUG ) ) debug() << "Failed to get x,y,z of tile " << id.muonID() << endmsg;
    throw GaudiException( "Failed to get tile positions", __func__, StatusCode::FAILURE );
  }

  if ( msgLevel( MSG::DEBUG ) )
    debug() << " Created muon measurement! " << id.muonID() << "  " << pos->position() << " dx " << pos->dX() << " dy "
            << pos->dY() << " dz " << pos->dZ() << endmsg;

  return makeMeasurements( id, det, pos->position(), pos->dX(), pos->dY() );
}
//-----------------------------------------------------------------------------
/// Create a measurement
//-----------------------------------------------------------------------------

std::array<LHCb::Measurement, 2> MuonMeasurementProvider::measurement( const DeMuonDetector&    det,
                                                                       LHCb::span<LHCb::LHCbID> ids ) const {
  using std::pow;
  std::vector<std::array<double, 2>> padx;
  std::vector<std::array<double, 2>> pady;
  std::vector<std::array<double, 2>> padz;
  double                             hit_minx = 100000;
  double                             hit_maxx = -100000;
  double                             hit_miny = 100000;
  double                             hit_maxy = -100000;
  double                             hit_minz = 100000;
  double                             hit_maxz = -100000;

  padx.reserve( ids.size() );
  pady.reserve( ids.size() );
  padz.reserve( ids.size() );
  for ( auto id : ids ) {
    if ( !id.isMuon() ) {
      error() << "Not a Muon measurement" << endmsg;
      continue;
    }
    auto pos = det.position( id.muonID() );
    if ( !pos ) {
      Warning( "Failed to get x,y,z of tile " ).ignore();
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Failed to get x,y,z of tile " << id.muonID() << endmsg;
      continue;
    }
    double x  = pos->position().X();
    double y  = pos->position().Y();
    double z  = pos->position().Z();
    double dx = pos->dX();
    double dy = pos->dY();
    double dz = pos->dZ();
    padx.push_back( {x, dx} );
    pady.push_back( {y, dy} );
    padz.push_back( {z, dz} );
    if ( ( x - dx ) < hit_minx ) hit_minx = x - dx;
    if ( ( x + dx ) > hit_maxx ) hit_maxx = x + dx;
    if ( ( y - dy ) < hit_miny ) hit_miny = y - dy;
    if ( ( y + dy ) > hit_maxy ) hit_maxy = y + dy;
    if ( ( z - dz ) < hit_minz ) hit_minz = z - dz;
    if ( ( z + dz ) > hit_maxz ) hit_maxz = z + dz;
  }

  auto [Cx, Cdx] = recomputePos( padx );
  auto [Cy, Cdy] = recomputePos( pady );
  auto [Cz, Cdz] = recomputePos( padz );

  LHCb::LHCbID Cid;
  double       min_dist2 = 100000000;
  for ( auto id : ids ) {
    if ( id.isMuon() ) {
      LHCb::Detector::Muon::TileID muid = id.muonID();

      auto pos = det.position( muid );
      if ( pos ) {
        double dist2 = std::pow( pos->position().X() - Cx, 2 ) + pow( pos->position().Y() - Cy, 2 ) +
                       pow( pos->position().Z() - Cz, 2 );
        if ( min_dist2 > dist2 ) {
          min_dist2 = dist2;
          Cid       = id; // choose the id closest to the cluster center
        }
      }
    }
  }
  if ( !Cid.muonID().isValid() ) {
    error() << " IMPOSSIBLE to Create muon measurement from a cluster of " << padx.size() << " hits ! " << Cid.muonID()
            << endmsg;
    error() << " x " << Cx << " +/- " << Cdx << " in the range [" << hit_minx << "," << hit_maxx << "] " << endmsg;
    error() << " y " << Cy << " +/- " << Cdy << " in the range [" << hit_miny << "," << hit_maxy << "] " << endmsg;
    error() << " z " << Cz << " +/- " << Cdz << " in the range [" << hit_minz << "," << hit_maxz << "] " << endmsg;
  }

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << " Created muon measurement from a cluster of " << padx.size() << " hits ! " << Cid.muonID() << endmsg;
    debug() << " x " << Cx << " +/- " << Cdx << " in the range [" << hit_minx << "," << hit_maxx << "] " << endmsg;
    debug() << " y " << Cy << " +/- " << Cdy << " in the range [" << hit_miny << "," << hit_maxy << "] " << endmsg;
    debug() << " z " << Cz << " +/- " << Cdz << " in the range [" << hit_minz << "," << hit_maxz << "] " << endmsg;
  }

  return makeMeasurements( Cid, det, Gaudi::XYZPoint( Cx, Cy, Cz ), Cdx, Cdy );
}

//-----------------------------------------------------------------------------
/// Create measurements for list of LHCbIDs
//-----------------------------------------------------------------------------

void MuonMeasurementProvider::addToMeasurements( LHCb::span<LHCb::LHCbID>        ids,
                                                 std::vector<LHCb::Measurement>& measurements,
                                                 const LHCb::ZTrajectory<double>& ) const {
  const auto& det = m_det.get();
  if ( !m_clusterize ) {

    if ( msgLevel( MSG::DEBUG ) ) debug() << " MuonMeasurementProvider makes measurements for each LHCbID" << endmsg;

    std::for_each( ids.begin(), ids.end(), [&]( const LHCb::LHCbID& id ) {
      auto m = measurement( det, id );
      measurements.push_back( m[0] );
      measurements.push_back( m[1] );
    } );

  } else {

    if ( msgLevel( MSG::DEBUG ) )
      debug() << " MuonMeasurementProvider makes measurements for each CLUSTER of LHCbIDs" << endmsg;

    auto first = ids.begin();
    auto last  = std::stable_partition( first, ids.end(), []( const LHCb::LHCbID& id ) { return id.isMuon(); } );

    for ( unsigned iMS = 0; iMS < (unsigned)det.stations(); ++iMS ) {
      auto pivot = std::stable_partition(
          first, last, [&]( const LHCb::LHCbID& id ) { return id.isMuon() && iMS == id.muonID().station(); } );
      if ( msgLevel( MSG::DEBUG ) )
        debug() << " Station M" << iMS << " lhcbIDS " << std::distance( first, pivot ) << endmsg;
      if ( std::distance( first, pivot ) != 0 ) {
        auto m = measurement( det, make_span( first, pivot ) );
        measurements.push_back( m[0] );
        measurements.push_back( m[1] );
      }
      first = pivot;
    }
  }
}
