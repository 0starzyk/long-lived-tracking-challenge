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
#include "Kernel/STLExtensions.h"
#include "Kernel/TrackDefaultParticles.h"

#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/SystemOfUnits.h"

#include <vector>

#include "DetDesc/IGeometryInfo.h"

/**
 *  Interface for a track fitting tool.
 *
 *  @author Jose A. Hernando, Eduardo Rodrigues
 *  @date   2005-05-25
 *
 *  @author Rutger van der Eijk  07-04-1999
 *  @author Mattiew Needham
 */
struct ITrackFitter : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackFitter, 6, 0 );

  // Fit a track
  virtual StatusCode operator()( LHCb::Track& track, IGeometryInfo const& geometry,
                                 LHCb::Tr::PID const& pid = LHCb::Tr::PID::Pion() ) const = 0;

  // Fit a batch of tracks
  // Note: Returns SUCCESS if all tracks are Fitted
  virtual StatusCode operator()( LHCb::span<LHCb::Track> tracks, IGeometryInfo const& geometry,
                                 LHCb::Tr::PID const& pid = LHCb::Tr::PID::Pion() ) const = 0;

  /// reset internal state of the fitter
  virtual void reset() = 0;
};
