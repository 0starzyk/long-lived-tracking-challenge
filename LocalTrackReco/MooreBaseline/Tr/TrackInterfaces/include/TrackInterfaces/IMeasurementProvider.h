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
#include "Event/ZTrajectory.h"
#include "GaudiKernel/IAlgTool.h"
#include "Kernel/LHCbID.h"
#include "Kernel/STLExtensions.h"
#include <vector>

// Forward declarations
namespace LHCb {
  class Measurement;
  class StateVector;
} // namespace LHCb

/** @class IMeasurementProvider IMeasurementProvider.h TrackInterfaces/IMeasurementProvider.h
 *
 *  Interface for the measurement provider tool
 *
 *  @author Jose Hernando
 *  @author Eduardo Rodrigues
 *  @date   2005-06-28
 */
struct IMeasurementProvider : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IMeasurementProvider, 5, 0 );

  /** Load (=create) all the Measurements from the list of LHCbIDs
   *  on the input Track
   */
  virtual StatusCode load( LHCb::Track& track ) const = 0;

  /** create measurements for a set of LHCbIDs **/
  virtual void addToMeasurements( LHCb::span<LHCb::LHCbID> ids, std::vector<LHCb::Measurement>&,
                                  const LHCb::ZTrajectory<double>& ) const = 0;
};
