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
#include "Event/FTLiteCluster.h"
#include "Event/MCParticle.h"
#include "GaudiAlg/GaudiTool.h"
#include "Kernel/LHCbID.h"
#include "Linker/LinkedTo.h"
#include "PrKernel/IPrDebugTool.h" // Interface

//-----------------------------------------------------------------------------
// Implementation file for class : PrTStationDebugTool
//
// 2012-03-22 : Olivier Callot
//-----------------------------------------------------------------------------

/** @class PrTStationDebugTool PrTStationDebugTool.h
 *
 *
 *  @author Olivier Callot
 *  @date   2012-03-22
 */
class PrTStationDebugTool : public extends<GaudiTool, IPrDebugTool> {
public:
  /// Standard constructor
  using extends::extends;

  bool matchKey( LHCb::LHCbID id, int key ) const override;

  void printKey( MsgStream& msg, LHCb::LHCbID id ) const override;
};
// Declaration of the Tool Factory
DECLARE_COMPONENT( PrTStationDebugTool )

//=========================================================================
//  Check if a given LHCbID is associated to the MCParticle of the specified key
//=========================================================================
bool PrTStationDebugTool::matchKey( LHCb::LHCbID id, int key ) const {
  if ( id.isFT() ) {
    auto links =
        SmartDataPtr<LHCb::LinksByKey>{evtSvc(), LHCb::LinksByKey::linkerName( LHCb::FTLiteClusterLocation::Default )};
    auto r = LinkedTo<LHCb::MCParticle>{links}.range( id.ftID() );
    return std::any_of( r.begin(), r.end(), [key]( const auto& i ) { return i.key() == key; } );
  }
  return false;
}
//=========================================================================
//  Print the list of MCParticle keys associated to the specified LHCbID
//=========================================================================
void PrTStationDebugTool::printKey( MsgStream& msg, LHCb::LHCbID id ) const {
  if ( id.isFT() ) {
    auto links =
        SmartDataPtr<LHCb::LinksByKey>{evtSvc(), LHCb::LinksByKey::linkerName( LHCb::FTLiteClusterLocation::Default )};
    auto r = LinkedTo<LHCb::MCParticle>{links}.range( id.ftID() );
    if ( !r.empty() ) msg << " MC:";
    for ( const auto& p : r ) msg << " " << p.key();
  }
}
//=============================================================================
