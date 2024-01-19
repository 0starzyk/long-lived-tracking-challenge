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
#include "Event/MCTrackInfo.h"
#include "Event/ProcStatus.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ParticleProperty.h"
#include "Linker/LinkedFrom.h"
#include "Linker/LinkedTo.h"
#include <vector>

//-----------------------------------------------------------------------------
// Implementation file for class : PrDebugTrackingLosses
//
// 2009-03-26 : Olivier Callot
//-----------------------------------------------------------------------------

/** @class PrDebugTrackingLosses PrDebugTrackingLosses.h
 *  Debug which MCParticles are not reconstructed.
 *
 *  @author Olivier Callot
 *  @date   2009-03-26
 */
class PrDebugTrackingLosses : public GaudiAlgorithm {
public:
  /// Standard constructor
  PrDebugTrackingLosses( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution
  StatusCode finalize() override;   ///< Algorithm finalization

private:
  void printMCParticle( const LHCb::MCParticle* part );

  mutable LHCb::IParticlePropertySvc* m_ppSvc = nullptr; ///< Pointer to particle property service

  bool                            m_velo;
  bool                            m_forward;
  bool                            m_seed;
  bool                            m_clone;
  bool                            m_ghost;
  bool                            m_fromStrange;
  bool                            m_fromBeauty;
  double                          m_minMomentum;
  bool                            m_saveList;
  std::string                     m_veloName;
  std::string                     m_forwardName;
  std::string                     m_seedName;
  std::vector<std::array<int, 3>> m_badGuys;
  int                             m_eventNumber = 0;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( PrDebugTrackingLosses )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrDebugTrackingLosses::PrDebugTrackingLosses( const std::string& name, ISvcLocator* pSvcLocator )
    : GaudiAlgorithm( name, pSvcLocator ) {
  declareProperty( "Velo", m_velo = false );
  declareProperty( "Forward", m_forward = false );
  declareProperty( "Seed", m_seed = false );
  declareProperty( "Ghost", m_ghost = false );
  declareProperty( "Clone", m_clone = false );
  declareProperty( "FromStrange", m_fromStrange = false );
  declareProperty( "FromBeauty", m_fromBeauty = false );
  declareProperty( "MinMomentum", m_minMomentum = 5000. );
  declareProperty( "SaveList", m_saveList = false );
  declareProperty( "VeloName", m_veloName = LHCb::TrackLocation::Velo );
  declareProperty( "ForwardName", m_forwardName = LHCb::TrackLocation::Forward );
  declareProperty( "SeedName", m_seedName = LHCb::TrackLocation::Seed );
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrDebugTrackingLosses::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;              // error printed already by GaudiAlgorithm

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Initialize" << endmsg;

  m_ppSvc = svc<LHCb::IParticlePropertySvc>( "LHCb::ParticlePropertySvc", true );
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode PrDebugTrackingLosses::execute() {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  if ( !m_velo && !m_forward && !m_seed ) return StatusCode::SUCCESS;

  ++m_eventNumber;

  const LHCb::MCParticles* partCont = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );

  LHCb::ProcStatus* procStat = getOrCreate<LHCb::ProcStatus, LHCb::ProcStatus>( LHCb::ProcStatusLocation::Default );

  if ( procStat->aborted() ) {
    debug() << "** Processing aborted. Don't analyse losses! " << endmsg;
    return StatusCode::SUCCESS;
  }

  auto velo_links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( m_veloName ) );

  const auto trackInfo = MCTrackInfo{*get<LHCb::MCProperty>( LHCb::MCPropertyLocation::TrackInfo )};

  for ( LHCb::MCParticle* part : *partCont ) {
    if ( 0 == trackInfo.fullInfo( part ) ) continue;
    if ( m_seed ) {
      if ( !trackInfo.hasT( part ) ) continue;
    } else {
      if ( !trackInfo.hasVeloAndT( part ) ) continue;
    }
    if ( abs( part->particleID().pid() ) == 11 ) continue; // reject electron
    if ( m_fromStrange || m_fromBeauty ) {
      bool                    isStrange = false;
      bool                    isBeauty  = false;
      const LHCb::MCParticle* mother    = part;
      while ( 0 != mother->originVertex() ) {
        mother = mother->originVertex()->mother();
        if ( 0 == mother ) break;
        if ( mother->particleID().pid() == 310 ) isStrange = true;
        if ( mother->particleID().pid() == 3122 ) isStrange = true;
        if ( mother->particleID().pid() == -3122 ) isStrange = true;
        if ( mother->particleID().hasBottom() && ( mother->particleID().isMeson() || mother->particleID().isBaryon() ) )
          isBeauty = true;
      }
      if ( m_fromStrange && !isStrange ) continue;
      if ( m_fromBeauty && !isBeauty ) continue;
    }
    if ( m_minMomentum > fabs( part->p() ) ) continue;

    auto vr = LinkedFrom<LHCb::Track>{velo_links}.range( part );
    if ( ( m_velo ) && vr.empty() && !m_clone ) {
      info() << "Missed Velo for MCParticle " << part->key() << " ";
      printMCParticle( part );
    } else if ( ( m_velo ) && m_clone && ++vr.begin() != vr.end() ) {
      info() << "Velo clone for particle " << part->key() << " ";
      printMCParticle( part );
    }

    if ( m_forward && !vr.empty() ) {
      auto links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( m_forwardName ) );
      auto r     = LinkedFrom<LHCb::Track>{links}.range( part );
      if ( r.empty() ) { // the original code had firstP(part)!=NULL, which I think should have been == instead...
        if ( !m_clone && !m_ghost ) {
          info() << "Missed Forward (Velo " << vr.front().key() << ") for MCParticle " << part->key() << " ";
          printMCParticle( part );
          if ( m_saveList ) { m_badGuys.emplace_back( std::array{m_eventNumber, vr.front().key(), part->key()} ); }
        }
      } else if ( m_clone && ++r.begin() != r.end() ) {
        info() << "Forward clone (Velo " << vr.front().key() << ") for MCParticle " << part->key() << " ";
        printMCParticle( part );
      }
    }

    if ( m_seed ) {
      auto links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( m_seedName ) );
      auto r     = LinkedFrom<LHCb::Track>{links}.range( part );
      if ( r.empty() ) {
        if ( !m_clone && !m_ghost ) {
          info() << "Missed Seed for MCParticle " << part->key() << " ";
          printMCParticle( part );
          if ( m_saveList ) { m_badGuys.emplace_back( std::array{m_eventNumber, 0, part->key()} ); }
        }
      } else if ( m_clone && ++r.begin() != r.end() ) {
        info() << "Seed clone for MCParticle " << part->key() << " ";
        printMCParticle( part );
      }
    }
  }

  if ( m_ghost ) {
    auto velo_links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( m_veloName ) );
    auto id_links   = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( "Pr/LHCbID" ) );

    std::string location = m_forwardName;
    if ( m_velo ) location = m_veloName;
    if ( m_seed ) location = m_seedName;
    auto loc_links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( location ) );
    for ( const auto& itT : *get<LHCb::Tracks>( location ) ) {
      if ( LinkedTo<LHCb::MCParticle>{loc_links}.range( itT ).empty() ) {
        info() << "Ghost track, nb " << itT->key();
        int vKey = -1;
        for ( const auto& itA : itT->ancestors() ) {
          info() << " from Velo " << itA->key() << " ";
          vKey = itA->key();
        }
        if ( vKey >= 0 ) {
          LHCb::MCParticle const* part = LinkedTo<LHCb::MCParticle>{velo_links}.range( vKey ).try_front();
          if ( part ) {
            printMCParticle( part );
          } else {
            info() << endmsg;
          }
        } else {
          info() << endmsg;
        }

        std::map<LHCb::MCParticle const*, int> listKeys;
        for ( auto id : itT->lhcbIDs() ) {
          if ( id.isVP() ) {
            LHCb::Detector::VPChannelID idV = id.vpID();
            info() << format( "   Velo Sensor %3d chip %3d col %3d row %3d", idV.sensor(), idV.chip(), idV.col(),
                              idV.row() );
          } else if ( id.isFT() ) {
            LHCb::Detector::FTChannelID ftID = id.ftID();
            info() << format( "    FT St%2d La%2d Pm%2d Cel%4d    ", ftID.station(), ftID.layer(), ftID.sipm(),
                              ftID.channel() );
          }
          for ( const auto& part : LinkedTo<LHCb::MCParticle>{id_links}.range( id.lhcbID() ) ) {
            info() << " " << part.key();
            listKeys[&part] += 1;
          }
          info() << endmsg;
        }
        for ( const auto& i : listKeys ) printMCParticle( i.first );
      }
    }
  }

  return StatusCode::SUCCESS;
}

//=========================================================================
//
//=========================================================================
void PrDebugTrackingLosses::printMCParticle( const LHCb::MCParticle* part ) {
  const LHCb::MCParticle* mother = part;
  double                  p      = double( int( part->p() ) / 1000. );
  info() << "MC: [" << p << " GeV]";
  while ( 0 != mother ) {
    const LHCb::ParticleProperty* pp = m_ppSvc->find( mother->particleID() );
    if ( 0 == pp ) {
      info() << mother->key() << "[" << mother->particleID().pid() << "]";
    } else {
      info() << mother->key() << "[" << pp->particle() << "]";
    }
    const LHCb::MCVertex* vert = mother->originVertex();
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
//  Finalize
//=============================================================================
StatusCode PrDebugTrackingLosses::finalize() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Finalize" << endmsg;

  if ( m_saveList ) {
    for ( const auto& i : m_badGuys ) {
      info() << format( "Event %4d  Velo %4d  MCParticle %4d", i[0], i[1], i[2] ) << endmsg;
    }
  }

  return GaudiAlgorithm::finalize(); // must be called after all other actions
}

//=============================================================================
