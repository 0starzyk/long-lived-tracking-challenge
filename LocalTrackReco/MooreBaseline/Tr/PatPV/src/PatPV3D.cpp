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

#include "TrackInterfaces/IPVOfflineTool.h"

#include "DetDesc/DetectorElement.h"
#include "Event/RecVertex.h"

#include "GaudiKernel/SystemOfUnits.h"
#include "LHCbAlgs/Transformer.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PatPV3D
//
// 2004-02-17 : Eduardo Rodrigues
//-----------------------------------------------------------------------------

/** @class PatPV3D PatPV3D.h
 *  Algorithm to find the primary vertices at the HLT.
 *
 *  @author Eduardo Rodrigues
 *  @author Sebastien Ponce
 */

//-----------------------------------------------------------------------------

class PatPV3D : public LHCb::Algorithm::MultiTransformerFilter<std::tuple<std::vector<LHCb::RecVertex>>(
                    std::vector<LHCb::Track> const&, DetectorElement const& )> {
public:
  /// Standard constructor
  PatPV3D( const std::string& name, ISvcLocator* pSvcLocator );

  /// Algorithm execution
  std::tuple<bool, std::vector<LHCb::RecVertex>> operator()( std::vector<LHCb::Track> const&,
                                                             DetectorElement const& lhcb ) const override;

private:
  Gaudi::Property<bool> m_refitpv{this, "RefitPV", false, "Flag to refit PVs when converting to type PrimaryVertex"};
  ToolHandle<IPVOfflineTool> m_pvsfit{"PVOfflineTool", this}; // PV fitting tool

  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_nbPVsCounter{this, "Nb PVs"};
  mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING>     m_failed{this, "reconstructMultiPV failed!"};
};

DECLARE_COMPONENT( PatPV3D )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PatPV3D::PatPV3D( const std::string& name, ISvcLocator* pSvcLocator )
    : MultiTransformerFilter( name, pSvcLocator,
                              {KeyValue{"InputTracks", LHCb::TrackLocation::Default},
                               KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}},
                              KeyValue( "OutputVerticesName", LHCb::RecVertexLocation::Velo3D ) ) {}

//=============================================================================
// Main execution
//=============================================================================
std::tuple<bool, std::vector<LHCb::RecVertex>> PatPV3D::operator()( std::vector<LHCb::Track> const& inputTracks,
                                                                    DetectorElement const&          lhcb ) const {

  if ( msgLevel( MSG::DEBUG ) ) { debug() << "==> Execute" << endmsg; }

  std::vector<LHCb::RecVertex> rvts;
  rvts.reserve( 20 );
  bool filter = false;
  m_pvsfit->reconstructMultiPV( inputTracks, rvts, *lhcb.geometry() )
      .andThen( [&] { filter = !rvts.empty(); } )
      .orElse( [&] { ++m_failed; } )
      .ignore();

  m_nbPVsCounter += rvts.size();
  return std::make_tuple( filter, std::move( rvts ) );
}
