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

#include "DetDesc/ILVolume.h"
#include "Event/TrackMaterialIntersection.h"
#include "Event/TrackTypes.h"
#include "Event/ZTrajectory.h"
#include "Kernel/TrackDefaultParticles.h"

#include "GaudiKernel/IAlgTool.h"

#include <any>

#include "DetDesc/IGeometryInfo.h"
namespace LHCb {
  class State;
  class StateVector;
} // namespace LHCb

/** @class IMaterialLocatorLocator
 *
 *  Interface for tools that locate materil intersections on a trajectory
 *
 *  @author Wouter Hulsbergen
 *  @date   2006-05-16
 */

struct IMaterialLocator : extend_interfaces<IAlgTool> {
  /// interface ID
  DeclareInterfaceID( IMaterialLocator, 3, 0 );

  /// embedded class representing intersection
  using Intersection = LHCb::TrackMaterialIntersection;

  /// container of intersections
  typedef std::vector<Intersection> Intersections;

  /// Create an instance of the accelerator cache
  virtual std::any createCache() const = 0;

  /// Intersect a trajectory with volumes in the geometry
  virtual Intersections intersect( const LHCb::ZTrajectory<double>& traj, std::any& accelCache,
                                   IGeometryInfo const& geometry ) const = 0;

  /// Apply material corrections using material in intersepts
  virtual void applyMaterialCorrections( LHCb::State& stateAtTarget, const Intersections& intersepts, double zorigin,
                                         const LHCb::Tr::PID pid, bool applyScatteringCorrection = true,
                                         bool applyELossCorrection = true ) const = 0;
};
