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

// Include files
// -------------
#include "boost/container/static_vector.hpp"
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/ToolHandle.h"

// from LHCbKernel
#include "Kernel/LHCbID.h"

// from TrackInterfaces
#include "TrackInterfaces/IMeasurementProviderProjector.h"
#include "TrackInterfaces/IUTClusterPosition.h"

// track kernel
#include "TrackKernel/TrackFunctors.h"
#include "TrackKernel/TrackTraj.h"

// Event
#include "Event/TrackFitResult.h"

class MeasurementProvider : public extends<GaudiTool, IMeasurementProviderProjector> {
public:
  /** standard tool constructor */
  MeasurementProvider( const std::string& type, const std::string& name, const IInterface* parent );

  /** initialize tool */
  StatusCode initialize() override;

  /** finalize tool */
  StatusCode finalize() override;

  /** See interface class */
  StatusCode load( LHCb::Track& track ) const override;

  /** See interface class */
  void addToMeasurements( LHCb::span<LHCb::LHCbID> ids, std::vector<LHCb::Measurement>& measurements,
                          const LHCb::ZTrajectory<double>& reftraj ) const override;

  StatusCode projectReference( LHCb::FitNode& node ) const override;
  void       reset() override;

private:
  // Handles to actual measurement providers
  ToolHandle<IMeasurementProviderProjector> m_vpProvider   = {"VPMeasurementProvider", this};
  ToolHandle<IMeasurementProviderProjector> m_utProvider   = {"UTMeasurementProvider", this};
  ToolHandle<IMeasurementProviderProjector> m_ftProvider   = {"FTMeasurementProvider", this};
  ToolHandle<IMeasurementProviderProjector> m_muonProvider = {"MuonMeasurementProvider", this};

  Gaudi::Property<bool> m_ignoreVP{this, "IgnoreVP",
                                   false}; // VP does not exist in default detector   ///< Ignore VP hits
  Gaudi::Property<bool> m_ignoreUT{this, "IgnoreUT", false};    ///< Ignore UT hits
  Gaudi::Property<bool> m_ignoreFT{this, "IgnoreFT", false};    ///< Ignore FT hits
  Gaudi::Property<bool> m_ignoreMuon{this, "IgnoreMuon", true}; ///< Ignore Muon hits
  Gaudi::Property<bool> m_initializeReference{
      this, "InitializeReference", true}; ///< Initialize measurement reference vector with closest state on track

  boost::container::static_vector<std::pair<LHCb::Measurement::Type, IMeasurementProviderProjector*>, 9> m_providers;
};
