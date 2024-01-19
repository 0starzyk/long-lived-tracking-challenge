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
#include "Event/TrackParameters.h"

#include "GaudiKernel/IAlgTool.h"

#include "DetDesc/IGeometryInfo.h"
namespace LHCb {
  class State;
  class TrackTraj;
} // namespace LHCb

/**
 *  Interface for track stateprovider. TrackStateProvider provides the
 *  state of a track at a certain z position using several methods. To
 *  make this reasonably fast for use in DaVinci it caches every state
 *  it has ever computed for a given track.
 *
 *  @author Wouter Hulsbergen
 *  @date 16/08/2010
 */

struct ITrackStateProvider : extend_interfaces<IAlgTool> {
  DeclareInterfaceID( ITrackStateProvider, 5, 0 );

  /// Compute the state of the track at position z.  The third
  /// argument is the tolerance: if an existing state is found within
  /// a z-distance 'tolerance', that state is returned.
  virtual StatusCode state( LHCb::State& state, const LHCb::Track& track, double z, IGeometryInfo const& geometry,
                            double ztolerance = TrackParameters::propagationTolerance ) const = 0;

  /// Compute the state of the track at a certain z position by using
  /// the trajectory approximation.
  virtual StatusCode stateFromTrajectory( LHCb::State& state, const LHCb::Track& track, double z,
                                          IGeometryInfo const& geometry ) const = 0;

  /// Return the trajectory approximation of this track. Ownership
  /// stays with the tool. A zero pointer indicates an error.
  virtual const LHCb::TrackTraj* trajectory( const LHCb::Track& track, IGeometryInfo const& geometry ) const = 0;

  /// Clear the cache
  virtual void clearCache() const = 0;

  /// Clear the cache for a particular track
  virtual void clearCache( const LHCb::Track& track ) const = 0;
};
