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

#include "Event/Track.h"
#include "Event/TrackTypes.h"
#include "Kernel/TrackDefaultParticles.h"

#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/Plane3DTypes.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "GaudiKernel/Vector3DTypes.h"

#include "DetDesc/IDetectorElement.h"
#include "DetDesc/IGeometryInfo.h"

namespace LHCb {
  class State;
  class StateVector;
  namespace Event {
    namespace v3 {
      struct States;
    }
  } // namespace Event
  namespace Magnet {
    class MagneticFieldGrid;
  }
} // namespace LHCb

/**
 *  Interface for track extrapolator tools
 *
 *  @author Edwin Bos (added method)
 *  @date 06/07/2005
 *  @author Eduardo Rodrigues (changes and new features for new track event model)
 *  @date   14/12/2004
 *  @author Marco Cattaneo
 *  @date   09/01/2002
 *
 *  Based on TrackExtrapolator ABS by Rutger van der Eijk, 07-04-1999
 */

struct ITrackExtrapolator : extend_interfaces<IAlgTool> {
  /// Return the interface ID
  DeclareInterfaceID( ITrackExtrapolator, 5, 0 );

  /// Whether it uses or not grid interpolation for the magnetic field
  virtual bool usesGridInterpolation() const { return true; }

  /// Propagate a state vector from zOld to zNew
  StatusCode propagate( Gaudi::TrackVector& stateVec, double zOld, double zNew, IGeometryInfo const& geometry,
                        const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                        const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    return propagate( stateVec, zOld, zNew, nullptr, geometry, pid, grid );
  }

  /// Propagate a state vector from zOld to zNew
  /// Transport matrix is calulated when transMat pointer is not NULL
  virtual StatusCode propagate( Gaudi::TrackVector& stateVec, double zOld, double zNew, Gaudi::TrackMatrix* transMat,
                                IGeometryInfo const& geometry, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  virtual LHCb::Event::v3::States propagate( const LHCb::Event::v3::States& states, double zNew,
                                             IGeometryInfo const&                   geometry,
                                             const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                                             const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  /// Propagate a statevector
  StatusCode propagate( LHCb::StateVector& state, double z, IGeometryInfo const& geometry,
                        Gaudi::TrackMatrix* transportmatrix = nullptr, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                        const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    return propagate( state.parameters(), state.z(), z, transportmatrix, geometry, pid, grid ).andThen( [&] {
      state.setZ( z );
    } );
  }

  /// Propagate a track to a given z-position
  virtual StatusCode propagate( const LHCb::Track& track, double z, LHCb::State& state, IGeometryInfo const& geometry,
                                const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  /// Propagate a track to a given z-position
  virtual StatusCode propagate( const LHCb::Track& track, double z, LHCb::StateVector& statevector,
                                IGeometryInfo const& geometry, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  /// Propagate a state to a given z-position
  virtual StatusCode propagate( LHCb::State& state, double z, IGeometryInfo const& geometry,
                                const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  /// Propagate a state to a given z-position
  /// Transport matrix is calulated when transMat pointer is not NULL
  virtual StatusCode propagate( LHCb::State& state, double z, Gaudi::TrackMatrix* transMat,
                                IGeometryInfo const& geometry, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  /// Propagate a track to the closest point to the specified point
  virtual StatusCode propagate( const LHCb::Track& track, const Gaudi::XYZPoint& point, LHCb::State& state,
                                IGeometryInfo const& geometry, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  /// Propagate a state to the closest point to the specified point
  virtual StatusCode propagate( LHCb::State& state, const Gaudi::XYZPoint& point, IGeometryInfo const& geometry,
                                const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  /// Propagate a track to within tolerance of a plane (default = 10 microns)
  virtual StatusCode propagate( const LHCb::Track& track, const Gaudi::Plane3D& plane, LHCb::State& state,
                                IGeometryInfo const& geometry, double tolerance = 0.01,
                                const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  /// Propagate a state to within tolerance of a plane (default = 10 microns)
  virtual StatusCode propagate( LHCb::State& state, const Gaudi::Plane3D& plane, IGeometryInfo const& geometry,
                                double tolerance = 0.01, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const = 0;

  //--------- ACCESS METHODS ---------------------------------------

  /// Retrieve the position and momentum vectors and the corresponding
  /// 6D covariance matrix (pos:1->3,mom:4-6) of a track at a given z-position
  StatusCode positionAndMomentum( const LHCb::Track& track, double z, Gaudi::XYZPoint& pos, Gaudi::XYZVector& mom,
                                  Gaudi::SymMatrix6x6& cov6D, IGeometryInfo const& geometry,
                                  const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                                  const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] {
      tmpState.positionAndMomentum( pos, mom, cov6D );
    } );
  }

  /// Retrieve the position and momentum vectors of a track at a given z-position
  StatusCode positionAndMomentum( const LHCb::Track& track, double z, Gaudi::XYZPoint& pos, Gaudi::XYZVector& mom,
                                  IGeometryInfo const& geometry, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                                  const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] {
      pos = tmpState.position();
      mom = tmpState.momentum();
    } );
  }

  /// Retrieve the 3D-position vector and error matrix of a track at a given z-position
  StatusCode position( const LHCb::Track& track, double z, Gaudi::XYZPoint& pos, Gaudi::SymMatrix3x3& errPos,
                       IGeometryInfo const& geometry, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                       const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] {
      pos    = tmpState.position();
      errPos = tmpState.errPosition();
    } );
  }

  /// Retrieve the 3D-position vector of a track at a given z-position
  StatusCode position( const LHCb::Track& track, double z, Gaudi::XYZPoint& pos, IGeometryInfo const& geometry,
                       const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                       const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] { pos = tmpState.position(); } );
  }

  /// Retrieve the slopes (dx/dz,dy/dz,1) and error matrix of a track at a given z-position
  StatusCode slopes( const LHCb::Track& track, double z, Gaudi::XYZVector& slopes, IGeometryInfo const& geometry,
                     Gaudi::SymMatrix3x3& errSlopes, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                     const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] {
      slopes    = tmpState.slopes();
      errSlopes = tmpState.errSlopes();
    } );
  }

  /// Retrieve the slopes (dx/dz,dy/dz,1) of a track at a given z-position
  StatusCode slopes( const LHCb::Track& track, double z, Gaudi::XYZVector& slopes, IGeometryInfo const& geometry,
                     const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                     const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] { slopes = tmpState.slopes(); } );
  }

  /// Retrieve the momentum of a track at a given z-position
  StatusCode p( const LHCb::Track& track, double z, double& p, IGeometryInfo const& geometry,
                const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] { p = tmpState.p(); } );
  }

  /// Retrieve the transverse momentum of a track at a given z-position
  StatusCode pt( const LHCb::Track& track, double z, double& pt, IGeometryInfo const& geometry,
                 const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                 const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] { pt = tmpState.pt(); } );
  }

  /// Retrieve the momentum vector and error matrix of a track at a given z-position
  StatusCode momentum( const LHCb::Track& track, double z, Gaudi::XYZVector& mom, Gaudi::SymMatrix3x3& errMom,
                       IGeometryInfo const& geometry, const LHCb::Tr::PID pid = LHCb::Tr::PID::Pion(),
                       const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] {
      mom    = tmpState.momentum();
      errMom = tmpState.errMomentum();
    } );
  }

  /// Retrieve the momentum vector of a track at a given z-position
  StatusCode momentum( const LHCb::Track& track, double z, Gaudi::XYZVector& mom, IGeometryInfo const& geometry,
                       const LHCb::Tr::PID                    pid  = LHCb::Tr::PID::Pion(),
                       const LHCb::Magnet::MagneticFieldGrid* grid = nullptr ) const {
    LHCb::State tmpState;
    return propagate( track, z, tmpState, geometry, pid, grid ).andThen( [&] { mom = tmpState.momentum(); } );
  }
};
