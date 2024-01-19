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
#ifndef TRACKPROJECTORS_TRACKPROJECTOR_H
#define TRACKPROJECTORS_TRACKPROJECTOR_H 1

// Include files

#include "Event/Measurement.h"

// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/Vector3DTypes.h"

// from TrackInterfaces
#include "Magnet/DeMagnet.h"
#include "TrackInterfaces/ITrackProjector.h"
#include <DetDesc/GenericConditionAccessorHolder.h>
struct ITrajPoca;

/** @class TrackProjector TrackProjector.h TrackProjectors/TrackProjector.h
 *
 *  TrackProjector is the base class implementing methods
 *  from the ITrackProjector interface.
 *
 *  @author Jose Hernando
 *  @author Eduardo Rodrigues
 *  @author Sebastien Ponce
 */

class TrackProjector : public LHCb::DetDesc::ConditionAccessorHolder<extends<GaudiTool, ITrackProjector>> {

public:
  using ConditionAccessorHolder::ConditionAccessorHolder;

  /// ::initialize
  StatusCode initialize() override;

  /// Project the state vector in this fitnode and update projection matrix and reference residual
  StatusCode projectReference( LHCb::FitNode& node ) const override;

  /// Retrieve the derivative of the residual wrt. the alignment parameters
  /// of the measurement. The details of the alignment transformation are
  /// defined in AlignTraj.
  Derivatives alignmentDerivatives( const LHCb::StateVector& state, const LHCb::Measurement& meas,
                                    const Gaudi::XYZPoint& pivot ) const override final;

  bool useBField() const override { return m_useBField; }

protected:
  /**
   * Helper struct storing the result of calls to internal_project
   */

  using ProjectResult = ITrackProjector::ProjectResult;

  virtual ProjectResult internal_project( const LHCb::StateVector& statevector, const LHCb::Measurement& meas ) const;

  Gaudi::Property<bool> m_useBField{this, "UseBField", false}; /// Create StateTraj with or without BField information.
  Gaudi::Property<double>     m_tolerance{this, "Tolerance",
                                      0.0005 * Gaudi::Units::mm}; ///< Required accuracy of the projection
  ConditionAccessor<DeMagnet> m_magneticfield{this, "Magnet", LHCb::Det::Magnet::det_path};
  ITrajPoca*                  m_poca = nullptr; ///< Pointer to the ITrajPoca interface

public:
  /// Project a state onto a measurement.
  virtual std::tuple<StatusCode, ProjectResult> project( const LHCb::State&       state,
                                                         const LHCb::Measurement& meas ) const override;
};
#endif // TRACKPROJECTORS_TRACKPROJECTOR_H
