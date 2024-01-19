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
#ifndef PVOFFLINERECALCULATE_H
#define PVOFFLINERECALCULATE_H 1

// Include files
// from Gaudi
#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

static const InterfaceID IID_PVOfflineRecalculate( "PVOfflineRecalculate", 1, 0 );
/** @class PVOfflineRecalculate PVOfflineRecalculate.h
 *
 *
 *  @author Mariusz Witek
 *  @date   2010-10-05
 */
class PVOfflineRecalculate : public GaudiTool {
public:
  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_PVOfflineRecalculate; }

  /// Standard constructor
  PVOfflineRecalculate( const std::string& type, const std::string& name, const IInterface* parent );

  // Initialization
  StatusCode initialize() override;

  // Finalize
  StatusCode finalize() override;

  StatusCode RecalculateVertex( const LHCb::RecVertex* pvin, const std::vector<const LHCb::Track*>& tracks2remove,
                                LHCb::RecVertex& pvout, IGeometryInfo const& geometry ) const;

  void print_stats();

private:
  Gaudi::Property<bool> m_updatePVTracks{this, "UpdatePVTracks", false};
  // Extrapolators
  ToolHandle<ITrackExtrapolator> m_linExtrapolator{this, "LinearExtrapolator", "TrackLinearExtrapolator"};
  ToolHandle<ITrackExtrapolator> m_fullExtrapolator{this, "FullExtrapolator", "TrackMasterExtrapolator"};

  mutable std::array<int, 9> m_counter_count;

private:
  bool remove_track( const LHCb::RecVertex* vtx, Gaudi::SymMatrix3x3& hess, ROOT::Math::SVector<double, 3>& d0vec,
                     double& dchi2, const LHCb::Track* vtrack, double wg, IGeometryInfo const& geometry ) const;

  void printRat( const std::string& mes, int a, int b );
};
#endif // PVOFFLINERECALCULATE_H
