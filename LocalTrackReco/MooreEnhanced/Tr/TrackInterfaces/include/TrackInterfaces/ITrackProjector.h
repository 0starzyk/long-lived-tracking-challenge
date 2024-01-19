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
#ifndef TRACKINTERFACES_ITRACKPROJECTOR_H
#define TRACKINTERFACES_ITRACKPROJECTOR_H 1

// Include files
// -------------
// from Gaudi
#include "GaudiKernel/IAlgTool.h"

// Geometry definitions
#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/Vector3DTypes.h"

// From TrackEvent
#include "Event/TrackTypes.h"

// Forward declarations
namespace LHCb {
  class State;
  class StateVector;
  class Measurement;
  class Node;
  class FitNode;
} // namespace LHCb

/** @class ITrackProjector ITrackProjector.h
 *
 *  Interface for tracking projector tools
 *
 *  @author Jose Hernando
 *  @author Eduardo Rodrigues
 *  @author Sebastien Ponce
 */
struct ITrackProjector : extend_interfaces<IAlgTool> {
  // Return the interface ID
  DeclareInterfaceID( ITrackProjector, 4, 0 );

  /// Project a state onto a measurement.
  struct Project1DResult final {
    double                         sMeas;
    double                         doca;
    ROOT::Math::SVector<double, 1> residual;
    ROOT::Math::SVector<double, 1> errMeasure;
    Gaudi::TrackProjectionMatrix1D H;
    Gaudi::XYZVector               unitPocaVector;
  };
  struct Project2DResult final {
    double                         sMeas;
    double                         doca;
    ROOT::Math::SVector<double, 2> residual;
    ROOT::Math::SVector<double, 2> errMeasure;
    Gaudi::TrackProjectionMatrix2D H;
    Gaudi::XYZVector               unitPocaVector;
  };
  using ProjectResult = std::variant<Project1DResult, Project2DResult>;

  virtual std::tuple<StatusCode, ProjectResult> project( const LHCb::State&       state,
                                                         const LHCb::Measurement& meas ) const = 0;
  /// Project the state vector in this fitnode and update projection matrix and reference residual
  virtual StatusCode projectReference( LHCb::FitNode& node ) const = 0;

  /// Retrieve the derivative of the residual wrt. the alignment parameters
  /// of the measurement. The details of the alignment transformation are
  /// defined in AlignTraj.
  typedef Gaudi::Matrix1x6 Derivatives;
  virtual Derivatives      alignmentDerivatives( const LHCb::StateVector& state, const LHCb::Measurement& meas,
                                                 const Gaudi::XYZPoint& pivot ) const = 0;

  /// Simple getter to know whether we have to use BField
  virtual bool useBField() const = 0;
};

//==============================================================================
//   end of class
//==============================================================================

#endif // TRACKINTERFACES_ITRACKPROJECTOR_H
