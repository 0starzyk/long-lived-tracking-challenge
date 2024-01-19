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
#include "Event/MCParticle.h"
#include "Event/VPCluster.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IRegistry.h"
#include "Kernel/LHCbID.h"
#include "Linker/LinkedTo.h"
#include "PrKernel/IPrDebugTool.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrPixelDebugTool
//
// 2007-10-25 : Olivier Callot
//-----------------------------------------------------------------------------

/** @class PrPixelDebugTool PrPixelDebugTool.h
 *  Debug Pixel processing using MC truth
 *
 *  @author Olivier Callot
 *  @date   2007-10-25
 */
class PrPixelDebugTool final : public extends<GaudiTool, IPrDebugTool> {
public:
  /// Standard constructor
  using extends::extends;

  bool matchKey( LHCb::LHCbID id, int key ) const override;
  void printKey( MsgStream& msg, LHCb::LHCbID id ) const override;

private:
  double xTrue( int key, double z );
  double yTrue( int key, double z );
};

DECLARE_COMPONENT( PrPixelDebugTool )

//=============================================================================
// Check if a given LHCbID is associated to the MCParticle of the specified key
//=============================================================================
bool PrPixelDebugTool::matchKey( LHCb::LHCbID id, int key ) const {
  auto links =
      SmartDataPtr<LHCb::LinksByKey>{evtSvc(), LHCb::LinksByKey::linkerName( LHCb::VPClusterLocation::Default )};
  LHCb::Detector::VPChannelID idV = id.vpID();
  auto                        r   = LinkedTo<LHCb::MCParticle>{links}.range( idV );
  return std::any_of( r.begin(), r.end(), [key]( const auto& p ) { return p.key() == key; } );
}

//=========================================================================
// Print the list of MCParticle keys associated to the specified LHCbID
//=========================================================================
void PrPixelDebugTool::printKey( MsgStream& msg, LHCb::LHCbID id ) const {
  auto links =
      SmartDataPtr<LHCb::LinksByKey>{evtSvc(), LHCb::LinksByKey::linkerName( LHCb::VPClusterLocation::Default )};
  for ( auto const& part : LinkedTo<LHCb::MCParticle>{links}.range( id.vpID() ) ) msg << " " << part.key();
}

//=========================================================================
// Calculate the x position of a particle at a given z
//=========================================================================
double PrPixelDebugTool::xTrue( int key, double z ) {
  LHCb::MCParticles const* parts = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );
  LHCb::MCParticle const*  part  = parts->object( key );
  if ( !part ) return -999.;
  const double    tx     = part->momentum().px() / part->momentum().pz();
  Gaudi::XYZPoint origin = part->originVertex()->position();
  return origin.x() + tx * ( z - origin.z() );
}

//=========================================================================
// Calculate the y position of a particle at a given z
//=========================================================================
double PrPixelDebugTool::yTrue( int key, double z ) {

  LHCb::MCParticles* parts = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );
  LHCb::MCParticle*  part  = parts->object( key );
  if ( NULL == part ) return -999.;
  const double    ty     = part->momentum().py() / part->momentum().pz();
  Gaudi::XYZPoint origin = part->originVertex()->position();
  return origin.y() + ty * ( z - origin.z() );
}
