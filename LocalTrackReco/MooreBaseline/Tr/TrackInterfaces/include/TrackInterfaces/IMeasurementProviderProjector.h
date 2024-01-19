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
#ifndef TRACKINTERFACES_IMEASUREMENTPROVIDERPROJECTOR_H
#define TRACKINTERFACES_IMEASUREMENTPROVIDERPROJECTOR_H 1

// Include files
// -------------
// from Gaudi
#include "GaudiKernel/IAlgTool.h"
#include <vector>

#include "Kernel/STLExtensions.h"

#include "Event/Track.h"
#include "Event/ZTrajectory.h"
#include "Kernel/LHCbID.h"

// Geometry definitions
#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/Point3DTypes.h"

// From TrackEvent
#include "Event/TrackTypes.h"

// Forward declarations
namespace LHCb {
  class Measurement;
  // from projectors
  class Node;
  class FitNode;
} // namespace LHCb

/** @class IMeasurementProviderProjector IMeasurementProviderProjector.h TrackInterfaces/IMeasurementProviderProjector.h
 *
 *  Interface for the measurement provider and projector tool
 *
 *  @author Andrii Usachov
 *  @date   2019-10-15
 */
struct IMeasurementProviderProjector : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IMeasurementProviderProjector, 5, 0 );

  /** Load (=create) all the Measurements from the list of LHCbIDs
   *  on the input Track
   */
  virtual StatusCode load( LHCb::Track& track ) const = 0;

  /** create measurements for a set of LHCbIDs **/
  virtual void addToMeasurements( LHCb::span<LHCb::LHCbID> ids, std::vector<LHCb::Measurement>&,
                                  const LHCb::ZTrajectory<double>& ) const = 0;

  /// Project the state vector in this fitnode and update projection matrix and reference residual
  virtual StatusCode projectReference( LHCb::FitNode& node ) const = 0;

  /// reset internal state of the provider, if any
  virtual void reset() = 0;

  typedef Gaudi::Matrix1x6 Derivatives;
};
#endif // TRACKINTERFACES_IMEASUREMENTPROVIDERPROJECTOR_H
