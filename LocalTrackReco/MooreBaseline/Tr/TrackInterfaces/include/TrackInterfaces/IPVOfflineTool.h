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

#include "Event/RecVertex.h"
#include "Event/Track.h"

#include "GaudiKernel/AlgTool.h"

#include <string>
#include <vector>

#include "DetDesc/IGeometryInfo.h"

struct IPVOfflineTool : extend_interfaces<IAlgTool> {
  // Retrieve interface ID
  DeclareInterfaceID( IPVOfflineTool, 4, 0 );

  /** Reconstructs single PV excluding tracks from a given vector.
   *
   *  The full reconstruction of single PV is performed for a given
   *  seed point with the tracks from the default location
   *  excluding tracks in tracks2exclude vector.
   *  It is not a refit. Procedure can fail due to tracks excluded.
   */
  virtual StatusCode reDoSinglePV( const LHCb::Tracks& inputTracks, const Gaudi::XYZPoint xyzseed,
                                   std::vector<const LHCb::Track*>& tracks2exclude, LHCb::RecVertex& outvtx,
                                   IGeometryInfo const& geometry ) const = 0;

  /** Reconstructs all PVs excluding tracks from a given vector and returns new PV.
   *
   *  The full reconstruction of all PVs is performed with the tracks
   *  from the default location excluding tracks in tracks2exclude vector.
   *  The returned outvtx matches invtx by comparing track contents.
   *  It is not a refit. Procedure can fail due to tracks excluded.
   */
  virtual StatusCode reDoMultiPV( const LHCb::Tracks& inputTracks, const LHCb::RecVertex& invtx,
                                  std::vector<const LHCb::Track*>& tracks2exclude, LHCb::RecVertex& outvtx,
                                  IGeometryInfo const& geometry ) const = 0;

  /// Reconstructs single PV with tracks specified in a vector.
  virtual StatusCode reconstructSinglePVFromTracks( const Gaudi::XYZPoint                  xyzseed,
                                                    const std::vector<const LHCb::Track*>& tracks2use,
                                                    LHCb::RecVertex& outvtx, IGeometryInfo const& geometry ) const = 0;

  /// Reconstructs all PVs with tracks specified in a vector. Return weights.
  virtual StatusCode reconstructMultiPVFromTracks( std::vector<const LHCb::Track*>& tracks2use,
                                                   std::vector<LHCb::RecVertex>&    outvtxVec,
                                                   IGeometryInfo const&             geometry ) const = 0;

  /// Reconstructs all PVs with tracks from default location.
  virtual StatusCode reconstructMultiPV( const std::vector<LHCb::Track>& inputTracks,
                                         std::vector<LHCb::RecVertex>&   outvtxVec,
                                         IGeometryInfo const&            geometry ) const = 0;

  /// Reconstructs single PV for a given seed with tracks from default location.
  virtual StatusCode reconstructSinglePV( const LHCb::Tracks& inputTracks, const Gaudi::XYZPoint xyzseed,
                                          LHCb::RecVertex& outvtx, IGeometryInfo const& geometry ) const = 0;

  /// Remove tracks from PV and recalculate PV parameters without a fit
  virtual StatusCode removeTracksAndRecalculatePV( const LHCb::RecVertex*                 pvin,
                                                   const std::vector<const LHCb::Track*>& tracks2remove,
                                                   LHCb::RecVertex& vtx, IGeometryInfo const& geometry ) const = 0;
};
