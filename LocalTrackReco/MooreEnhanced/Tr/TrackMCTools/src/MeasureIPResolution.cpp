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
#include "Event/MCVertex.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ParticleProperty.h"
#include "Linker/LinkedTo.h"
#include <string>

//-----------------------------------------------------------------------------
// Implementation file for class : MeasureIPResolution
//
// 2010-10-05 : Olivier Callot
//-----------------------------------------------------------------------------

/** @class MeasureIPResolution MeasureIPResolution.h
 *  Measure the IP resolution of Velo tracks
 *
 *  @author Olivier Callot
 *  @date   2010-10-05
 */
class MeasureIPResolution : public GaudiAlgorithm {
public:
  /// Standard constructor
  MeasureIPResolution( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution
  StatusCode finalize() override;   ///< Algorithm finalization

protected:
  void printMCParticle( const LHCb::MCParticle* part );

private:
  std::string                         m_containerName;
  mutable LHCb::IParticlePropertySvc* m_ppSvc; ///< Pointer to particle property service

  int    m_nTracks;
  double m_averX;
  double m_averY;
  int    m_nbInCore;
  double m_sumRInCore;
  double m_sumR2InCore;
  double m_sumIPS;
  double m_sumIPS2;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( MeasureIPResolution )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
MeasureIPResolution::MeasureIPResolution( const std::string& name, ISvcLocator* pSvcLocator )
    : GaudiAlgorithm( name, pSvcLocator ) {
  declareProperty( "ContainerName", m_containerName = LHCb::TrackLocation::Velo );
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode MeasureIPResolution::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;              // error printed already by GaudiAlgorithm

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Initialize" << endmsg;

  m_ppSvc = svc<LHCb::IParticlePropertySvc>( "LHCb::ParticlePropertySvc", true );

  m_nTracks     = 0;
  m_averX       = 0.;
  m_averY       = 0.;
  m_nbInCore    = 0;
  m_sumRInCore  = 0.;
  m_sumR2InCore = 0.;
  m_sumIPS      = 0.;
  m_sumIPS2     = 0.;
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode MeasureIPResolution::execute() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  //== Get the MC primary vertices
  std::vector<Gaudi::XYZPoint> pvPosition;
  LHCb::MCVertices*            verts = get<LHCb::MCVertices>( LHCb::MCVertexLocation::Default );
  for ( const auto& vtx : *verts ) {
    if ( vtx->isPrimary() ) {
      pvPosition.push_back( vtx->position() );
      debug() << format( "Found PV at %6.3f %6.3f %8.3f", vtx->position().x(), vtx->position().y(),
                         vtx->position().z() )
              << endmsg;
    }
  }

  const LHCb::Tracks* tracks      = get<LHCb::Tracks>( m_containerName );
  auto                track_links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( m_containerName ) );

  for ( LHCb::Track const* trk : *tracks ) {
    if ( trk->isVeloBackward() ) continue;
    LHCb::MCParticle const* part = LinkedTo<LHCb::MCParticle>{track_links}.range( trk ).try_front();
    /*
    if ( 0 == part ) continue;
    if ( !part->originVertex()->isPrimary() ) {
      bool isClose = false;
      for ( std::vector<Gaudi::XYZPoint>::iterator itV = pvPosition.begin(); pvPosition.end() != itV; ++itV ) {
        if ( fabs( (*itV).x() - part->originVertex()->position().x() ) > 0.01 ) continue;
        if ( fabs( (*itV).y() - part->originVertex()->position().y() ) > 0.01 ) continue;
        if ( fabs( (*itV).z() - part->originVertex()->position().z() ) > 0.1 ) continue;
        isClose = true;
        break;
      }
      if ( !isClose) continue;
    }
    //if ( part->pt() < 300. ) continue;
    */

    LHCb::State const& state   = trk->firstState();
    double             minIP   = 1.e9;
    double             dxAtMin = 1.e9;
    double             dyAtMin = 1.e9;
    for ( const Gaudi::XYZPoint& vtx : pvPosition ) {
      double dx = state.x() + state.tx() * ( vtx.z() - state.z() ) - vtx.x();
      double dy = state.y() + state.ty() * ( vtx.z() - state.z() ) - vtx.y();
      double ip = sqrt( dx * dx + dy * dy );
      if ( ip < minIP ) {
        minIP   = ip;
        dxAtMin = dx;
        dyAtMin = dy;
      }
    }
    m_nTracks++;
    m_averX += dxAtMin;
    m_averY += dyAtMin;
    double signedIP = minIP;
    if ( dxAtMin * state.tx() + dyAtMin * state.ty() > 0. ) signedIP = -signedIP;
    double ipS = signedIP / sqrt( state.errX2() + state.errY2() );
    if ( minIP < 0.5 ) {
      m_nbInCore++;
      m_sumRInCore += signedIP;
      m_sumR2InCore += minIP * minIP;
      m_sumIPS += ipS;
      m_sumIPS2 += ipS * ipS;
    }
    if ( msgLevel( MSG::DEBUG ) ) {
      if ( 1. < ipS ) {
        info() << format( "Track %4d : signedIP %7.3f ipS %7.3f ", trk->key(), signedIP, ipS );
        if ( !part ) {
          info() << endmsg;
        } else {
          printMCParticle( part );
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode MeasureIPResolution::finalize() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Finalize" << endmsg;

  if ( m_nTracks > 0 ) {
    m_averX /= m_nTracks;
    m_averY /= m_nTracks;
    info() << format( "On %7d tracks, <x> %8.3f  <y> %8.3f", m_nTracks, m_averX, m_averY ) << endmsg;
    if ( m_nbInCore > 0 ) {
      double aver  = m_sumRInCore / m_nbInCore;
      double sigma = sqrt( m_sumR2InCore / m_nbInCore - aver * aver );
      double sIPS  = sqrt( m_sumIPS2 * m_nbInCore - m_sumIPS * m_sumIPS ) / m_nbInCore;
      info() << format( "   In core: %7d tracks, <> %8.3f sigma %8.3f  sigmaIPS %8.3f outside%7d", m_nbInCore, aver,
                        sigma, sIPS, m_nTracks - m_nbInCore )
             << endmsg;
    }
  }
  return GaudiAlgorithm::finalize(); // must be called after all other actions
}

//=========================================================================
//  Print the info on the particle
//=========================================================================
void MeasureIPResolution::printMCParticle( const LHCb::MCParticle* part ) {
  const LHCb::MCParticle* mother = part;
  const LHCb::MCVertex*   vert   = part->originVertex();
  double                  p      = double( int( part->p() ) / 1000. );
  info() << "MC: [" << p << " GeV]";
  while ( 0 != mother ) {
    const LHCb::ParticleProperty* pp = m_ppSvc->find( mother->particleID() );
    if ( 0 == pp ) {
      info() << mother->key() << "[" << mother->particleID().pid() << "]";
    } else {
      info() << mother->key() << "[" << pp->particle() << "]";
    }
    vert = mother->originVertex();
    if ( 0 == vert ) {
      mother = 0;
    } else {
      info() << format( " <-(z=%7.2f)", vert->position().z() );
      mother = vert->mother();
    }
  }
  info() << endmsg;
}
//=============================================================================
