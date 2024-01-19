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
//-----------------------------------------------------------------------------
// Implementation file for class : PatKShortChecker
//
// 2002-11-23 : Olivier Callot
//-----------------------------------------------------------------------------
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/MCVertex.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Linker/LinkerTable.h"
#include <string>

/** @class PatKShortChecker PatKShortChecker.h
 *  Monitor the KShort in an event
 *
 *  @author Olivier Callot
 *  @date   2002-11-23
 *  @adapt to A-Team framework 2007-08-20 SHM
 */

class PatKShortChecker : public GaudiAlgorithm {
public:
  /// Standard constructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution
  StatusCode finalize() override;   ///< Algorithm finalization

protected:
  std::string trackType( const LHCb::Track* tr );

  bool isKChild( const LHCb::MCParticle* part );

private:
  Gaudi::Property<std::string> m_inputLocation{this, "InputLocation", LHCb::TrackLocation::Default};

  std::vector<int> m_counter;

  std::vector<double> m_cntSeed;
  std::vector<double> m_cntDown;
};

DECLARE_COMPONENT( PatKShortChecker )

//=============================================================================
// Initialisation. Check parameters
//=============================================================================
StatusCode PatKShortChecker::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;              // error printed already by GaudiAlgorithm

  debug() << "==> Initialize" << endmsg;

  unsigned int kk;
  for ( kk = 0; 20 > kk; kk++ ) m_counter.push_back( 0 );

  for ( kk = 0; 10 > kk; kk++ ) {
    m_cntSeed.push_back( 0. );
    m_cntDown.push_back( 0. );
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode PatKShortChecker::execute() {

  debug() << "==> Execute" << endmsg;

  if ( !exist<LHCb::Tracks>( m_inputLocation ) ) return StatusCode::SUCCESS;

  LHCb::Tracks const* tracks = get<LHCb::Tracks>( m_inputLocation );
  LHCb::Tracks const* seeds  = get<LHCb::Tracks>( LHCb::TrackLocation::Seed );
  LHCb::Tracks const* downs  = get<LHCb::Tracks>( LHCb::TrackLocation::Downstream );

  LHCb::MCParticles const* parts   = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );
  LHCb::LinksByKey*        trlinks = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( m_inputLocation ) );
  LHCb::LinksByKey* seedlinks      = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::TrackLocation::Seed ) );
  LHCb::LinksByKey* downlinks =
      get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::TrackLocation::Downstream ) );

  // try and load the linker tables
  if ( trlinks ) trlinks->resolveLinks( evtSvc() );
  if ( seedlinks ) seedlinks->resolveLinks( evtSvc() );
  if ( downlinks ) downlinks->resolveLinks( evtSvc() );
  const LinkerTable<LHCb::Track, LHCb::MCParticle> trTable{trlinks};
  const LinkerTable<LHCb::Track, LHCb::MCParticle> seedTable{seedlinks};
  const LinkerTable<LHCb::Track, LHCb::MCParticle> downTable{downlinks};

  const auto trackInfo = MCTrackInfo{*get<LHCb::MCProperty>( LHCb::MCPropertyLocation::TrackInfo )};

  if ( seedTable ) {
    debug() << "Start counting seeds" << endmsg;

    for ( auto const& itT : *seeds ) {
      auto range = seedTable.relations( itT );
      if ( range.empty() ) {
        m_cntSeed[0] += 1.;
        verbose() << "No truth for track " << itT->key() << endmsg;
      } else {
        m_cntSeed[1] += 1.;
        verbose() << "Truth for track " << itT->key() << " : ";
        for ( auto it = range.begin(); range.end() != it; it++ ) {
          verbose() << it->to()->key() << " ";

          if ( trackInfo.hasUT( it->to() ) ) {
            m_cntSeed[2] += 1.;
            bool KChild = isKChild( it->to() );
            if ( KChild ) { m_cntSeed[4] += 1.; }
            if ( 5. * Gaudi::Units::GeV < it->to()->p() ) {
              m_cntSeed[3] += 1.;
              if ( KChild ) m_cntSeed[5] += 1.;
            }
          }
        }
        verbose() << endmsg;
      }
    }
  }

  if ( downTable ) {
    debug() << "Start counting Downstream " << endmsg;

    for ( auto const& itT : *downs ) {
      auto range = downTable.relations( itT );
      if ( range.empty() ) {
        m_cntDown[0] += 1.;
      } else {
        m_cntDown[1] += 1.;
        for ( auto const& it : range ) {
          if ( trackInfo.hasUT( it.to() ) ) {
            m_cntDown[2] += 1.;
            bool KChild = isKChild( it.to() );
            if ( KChild ) { m_cntDown[4] += 1.; }
            if ( 5. * Gaudi::Units::GeV < it.to()->p() ) {
              m_cntDown[3] += 1.;
              if ( KChild ) m_cntDown[5] += 1.;
            }
          }
        }
      }
    }
  }

  //== Count the topology of true KShorts

  debug() << "Start counting true KShorts " << endmsg;

  m_counter[0]++;

  for ( LHCb::MCParticle const* part : *parts ) {

    if ( part->particleID().pid() != 310 ) continue;
    const LHCb::MCVertex* vert = part->originVertex();
    if ( !vert ) continue;
    const LHCb::MCParticle* mother = vert->mother();
    if ( !mother ) continue;
    if ( !( mother->particleID().hasBottom() &&
            ( mother->particleID().isMeson() || mother->particleID().isBaryon() ) ) )
      continue;

    m_counter[1]++;

    std::vector<LHCb::MCParticle const*> children;

    unsigned int nWithVelo = 0;

    for ( auto itV : part->endVertices() ) {
      for ( auto& decay : itV->products() ) {
        if ( !trackInfo.hasT( decay ) ) continue;
        if ( !trackInfo.hasUT( decay ) ) continue;
        if ( trackInfo.hasVelo( decay ) ) nWithVelo++;
        children.push_back( decay.target() );
      }
    }

    if ( 2 == children.size() ) {
      m_counter[2]++;
      m_counter[12 + nWithVelo]++;

      debug() << "== KShort from B == N with velo : " << nWithVelo << endmsg;
      unsigned int nbReco = 0;
      unsigned int nbLong = 0;
      unsigned int nbDown = 0;
      unsigned int nbSeed = 0;
      for ( auto itPi = children.begin(); children.end() != itPi; itPi++ ) {

        bool found      = false;
        bool longTrack  = false;
        bool downstream = false;
        bool seed       = false;
        if ( trTable ) {
          for ( auto itT = tracks->begin(); tracks->end() != itT; itT++ ) {
            for ( auto const& it : trTable.relations( *itT ) ) {
              if ( *itPi == it.to() ) {
                found = true;
                if ( LHCb::Track::Types::Long == ( *itT )->type() ) longTrack = true;
                if ( LHCb::Track::Types::Downstream == ( *itT )->type() ) downstream = true;
                if ( LHCb::Track::Types::Ttrack == ( *itT )->type() ) seed = true;

                debug() << "  -- child found as " << ( *itT )->key() << " type " << trackType( ( *itT ) ) << endmsg;
              }
            }
          }
        }
        if ( found ) {
          nbReco++;
          if ( longTrack ) {
            nbLong++;
          } else if ( downstream ) {
            nbDown++;
          } else if ( seed ) {
            nbSeed++;
          }
        }
      }
      if ( children.size() == nbReco ) {
        m_counter[3]++;
        if ( 2 == nbReco ) {
          m_counter[4]++;
          if ( 2 == nbLong ) {
            m_counter[5]++;
          } else if ( 1 == nbLong && 1 == nbDown ) {
            m_counter[6]++;
          } else if ( 2 == nbDown ) {
            m_counter[7]++;
          } else if ( 1 == nbSeed && 1 == nbLong ) {
            m_counter[8]++;
          } else if ( 1 == nbSeed && 1 == nbDown ) {
            m_counter[9]++;
          } else if ( 2 == nbSeed ) {
            m_counter[10]++;
          } else {
            m_counter[11]++;
          }
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PatKShortChecker::finalize() {

  if ( 0 != m_counter[0] ) {
    info() << format( "Nb events          %6d", m_counter[0] ) << endmsg;
    double frac = m_counter[1] / (double)m_counter[0];
    info() << format( "Nb KShort from B   %6d %6.2f per event", m_counter[1], frac ) << endmsg;
    frac = m_counter[2] / (double)m_counter[0];
    info() << format( "Nb reconstructible %6d %6.2f per event", m_counter[2], frac ) << endmsg;
    if ( 0 != m_counter[2] ) {
      double den = 100. / m_counter[2];

      frac = m_counter[12] * den;
      info() << format( "   with 0 Velo    %6d %6.2f %%", m_counter[12], frac ) << endmsg;
      frac = m_counter[13] * den;
      info() << format( "   with 1 Velo    %6d %6.2f %%", m_counter[13], frac ) << endmsg;
      frac = m_counter[14] * den;
      info() << format( "   with 2 Velo    %6d %6.2f %%", m_counter[14], frac ) << endmsg;

      frac = m_counter[3] / (double)m_counter[2];
      info() << format( "Nb with tracks    %6d %6.2f per reconstructible", m_counter[3], frac ) << endmsg;
      frac = m_counter[4] / (double)m_counter[2];
      info() << format( "Nb with 2 tracks  %6d %6.2f per reconstructible", m_counter[4], frac ) << endmsg;
    }
    if ( 0 != m_counter[4] ) {
      double den = 100. / m_counter[4];

      info() << format( "  2 Long          %6d %6.2f %% of 2-tracks", m_counter[5], m_counter[5] * den ) << endmsg;
      info() << format( "  1 Long  1 Dwnst %6d %6.2f %% of 2-tracks", m_counter[6], m_counter[6] * den ) << endmsg;
      info() << format( "  2 downstream    %6d %6.2f %% of 2-tracks", m_counter[7], m_counter[7] * den ) << endmsg;
      info() << format( "  1 Long 1 Seed   %6d %6.2f %% of 2-tracks", m_counter[8], m_counter[8] * den ) << endmsg;
      info() << format( "  1 Dwnst 1 Seed  %6d %6.2f %% of 2-tracks", m_counter[9], m_counter[9] * den ) << endmsg;
      info() << format( "  2 Seed          %6d %6.2f %% of 2-tracks", m_counter[10], m_counter[10] * den ) << endmsg;
      info() << format( "  Other...        %6d %6.2f %% of 2-tracks", m_counter[11], m_counter[11] * den ) << endmsg;
    }
  }

  if ( 0 < m_cntSeed[1] ) {
    info() << "                     "
           << "Ghosts      Found          WithUT        >5GeV "
           << "   UT+KChild    UT+K+>5GeV" << endmsg;
    double frac = 100. * m_cntSeed[0] / ( m_cntSeed[0] + m_cntSeed[1] );
    info() << format( "Seed           %7.0f %5.1f%7.0f      ", m_cntSeed[0], frac, m_cntSeed[1] );

    double mult = 100. / m_cntSeed[1];
    info() << format( "%7.0f %5.1f", m_cntSeed[2], mult * m_cntSeed[2] );
    info() << format( "%7.0f %5.1f", m_cntSeed[3], mult * m_cntSeed[3] );
    info() << format( "%7.0f %5.1f", m_cntSeed[4], mult * m_cntSeed[4] );
    info() << format( "%7.0f %5.1f", m_cntSeed[5], mult * m_cntSeed[5] );
    info() << endmsg;

    if ( 0 < m_cntDown[1] ) {
      // frac = 100. * m_cntDown[0] / (m_cntDown[0] + m_cntDown[1]);
      // info() << format( "Downstream   %7.0f %5.1f %7.0f",
      //               m_cntDown[0], frac, m_cntDown[1] );

      //     double mult = 100. /  m_cntDown[1];
      // info() << format( "%7.0f %5.1f", m_cntDown[2], mult * m_cntDown[2]);
      // info() << format( "%7.0f %5.1f", m_cntDown[3], mult * m_cntDown[3]);
      // info() << format( "%7.0f %5.1f", m_cntDown[4], mult * m_cntDown[4]);
      // info() << format( "%7.0f %5.1f", m_cntDown[5], mult * m_cntDown[5]);
      // info() << endmsg;

      //=== Ratios, Downstream efficiency...

      info() << format( "Downstream/Seed%7.0f      ", m_cntDown[0] );
      for ( unsigned int kk = 1; 6 > kk; kk++ ) {
        if ( 0 < m_cntSeed[kk] ) {
          frac = 100. * m_cntDown[kk] / m_cntSeed[kk];
        } else {
          frac = 0.;
        }
        info() << format( "%7.0f %5.1f", m_cntDown[kk], frac );
      }
      info() << endmsg;
    }
  }

  return GaudiAlgorithm::finalize();
}

//=========================================================================
// returns a string with the type of the track
//=========================================================================
std::string PatKShortChecker::trackType( const LHCb::Track* tr ) {
  std::string result = "";
  if ( !tr->checkFlag( LHCb::Track::Flags::Clone ) ) { result = "unique "; }
  switch ( tr->type() ) {
  case LHCb::Track::Types::Velo:
    return result + "Velo ";
  case LHCb::Track::Types::Long:
    return result + "Long ";
  case LHCb::Track::Types::Ttrack:
    return result + "TTrack ";
  case LHCb::Track::Types::Upstream:
    return result + "Upstream ";
  case LHCb::Track::Types::Downstream:
    return result + "Downstream ";
  default:
    return result;
  }
}

//=========================================================================
//  Return if there is a K0S in the motherhood
//=========================================================================
bool PatKShortChecker::isKChild( const LHCb::MCParticle* part ) {
  const LHCb::MCVertex* vert = part->originVertex();
  if ( !vert ) return false;
  auto mother = vert->mother();
  return mother && mother->particleID().pid() == 310;
}
//=============================================================================
