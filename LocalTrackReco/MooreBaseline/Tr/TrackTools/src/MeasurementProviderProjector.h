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

#pragma once
#include <DetDesc/GenericConditionAccessorHolder.h>
#include <Event/FitNode.h>
#include <Event/Measurement.h>
#include <Event/State.h>
#include <Event/StateVector.h>
#include <Event/StateZTraj.h>
#include <Event/TrackParameters.h>
#include <Event/ZTrajectory.h>
#include <GaudiAlg/GaudiTool.h>
#include <GaudiKernel/ToolHandle.h>
#include <Magnet/DeMagnet.h>
#include <TrackInterfaces/IMeasurementProviderProjector.h>
#include <type_traits>

/** @class MeasurementProviderProjector MeasurementProviderProjector.h
 *
 * Implementation of templated MeasurementProviderProjector tool
 * see interface header for description
 *
 *  @author A. Usachov
 *  @date   15/10/2019
 */
class MeasurementProviderProjector
    : public LHCb::DetDesc::ConditionAccessorHolder<extends<GaudiTool, IMeasurementProviderProjector>> {
public:
  using ConditionAccessorHolder::ConditionAccessorHolder;

  /// Project the state vector in this fitnode and update projection matrix and reference residual
  StatusCode projectReference( LHCb::FitNode& node ) const override;

  /// reset internal state of the provider, if any
  void reset() override{};

  virtual bool useBField() const { return m_useBField; }

protected:
  /**
   * Helper struct storing the result of calls to internal_project
   */
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

  virtual ProjectResult internal_project( const LHCb::StateVector& statevector, const LHCb::Measurement& meas ) const;

  Gaudi::Property<bool> m_useBField{this, "UseBField", false}; /// Create StateTraj with or without BField information.

public:
  /// Project a state onto a measurement.
  virtual std::tuple<StatusCode, ProjectResult> project( const LHCb::State&       state,
                                                         const LHCb::Measurement& meas ) const;

  ConditionAccessor<DeMagnet> m_magnet{this, "Magnet", LHCb::Det::Magnet::det_path};
};
