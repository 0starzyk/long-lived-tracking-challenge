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
/** @file TrackSelector.cpp
 *
 *  Implementation file for reconstruction tool : TrackSelector
 *
 *  @author M.Needham Matt.Needham@cern.ch
 *  @author Chris Jones   Christopher.Rob.Jones@cern.ch
 *  @date   30/12/2005
 */
//-----------------------------------------------------------------------------

#include "GaudiKernel/SystemOfUnits.h"

// Tsa
#include "Event/PrKalmanFitResult.h"
#include "Event/TrackFitResult.h"
#include "Kernel/HitPattern.h"
#include "TrackKernel/TrackFunctors.h"
#include "TrackSelector.h"

using namespace LHCb;

DECLARE_COMPONENT( TrackSelector )

//-----------------------------------------------------------------------------

bool TrackSelector::accept( const Track& aTrack ) const {

  // Use a try block to catch exceptions from Track and/or State classes
  try {

    if ( msgLevel( MSG::VERBOSE ) ) {
      verbose() << "Trying Track " << aTrack.key() << " " << aTrack.type();
      if ( !aTrack.states().empty() ) verbose() << " P=" << aTrack.p() << " Pt=" << aTrack.pt();
      verbose() << endmsg;
    }

    // NDOF
    const int ndof = aTrack.nDoF();
    if ( ndof < m_minNDoF || ndof > m_maxNDoF ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> NDoF " << aTrack.nDoF() << " failed cut" << endmsg;
      return false;
    }

    // chi-squared
    double chi2 = aTrack.chi2PerDoF();
    if ( ( m_maxChi2Cut >= 0 && ( chi2 > m_maxChi2Cut || aTrack.nDoF() <= 0 ) ) ||
         ( m_minChi2Cut >= 0 && ( chi2 < m_minChi2Cut || aTrack.nDoF() <= 0 ) ) ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Chi^2 " << chi2 << " failed cut" << endmsg;
      return false;
    }

    const auto* fit   = dynamic_cast<const TrackFitResult*>( aTrack.fitResult() );
    const auto* prfit = dynamic_cast<const PrKalmanFitResult*>( aTrack.fitResult() );

    if ( ( m_maxChi2Velo > 0 ) || ( m_maxChi2Upstream > 0 ) || ( m_maxChi2Downstream > 0 ) || ( m_maxChi2Match > 0 ) ) {
      if ( !fit && !prfit )
        throw GaudiException( "TrackSelector", "Unknown or empty TrackFitResult", StatusCode::FAILURE );
    }

    if ( m_maxChi2Velo > 0 &&
         ( chi2 = fit ? fit->chi2Velo().chi2PerDoF() : prfit->chi2Velo().chi2PerDoF() ) > m_maxChi2Velo ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Velo Chi^2 " << chi2 << " failed cut" << endmsg;
      return false;
    }
    if ( m_maxChi2Upstream > 0 &&
         ( chi2 = fit ? fit->chi2Upstream().chi2PerDoF() : prfit->chi2Upstream().chi2PerDoF() ) > m_maxChi2Upstream ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Upstream Chi^2 " << chi2 << " failed cut" << endmsg;
      return false;
    }
    if ( m_maxChi2Downstream > 0 && ( chi2 = fit ? fit->chi2Downstream().chi2PerDoF()
                                                 : prfit->chi2Downstream().chi2PerDoF() ) > m_maxChi2Downstream ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Downstream Chi^2 " << chi2 << " failed cut" << endmsg;
      return false;
    }
    if ( m_maxChi2Match > 0 &&
         ( chi2 = fit ? fit->chi2Match().chi2PerDoF() : prfit->chi2Match().chi2PerDoF() ) > m_maxChi2Match ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Match Chi^2 " << chi2 << " failed cut" << endmsg;
      return false;
    }

    // cut p
    const double p = aTrack.p();
    if ( p < m_minPCut || ( m_maxPCut > 0 && p > m_maxPCut ) ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> P " << aTrack.p() << " failed cut" << endmsg;
      return false;
    }

    // cut on pt
    const double pt = aTrack.pt();
    if ( pt < m_minPtCut || ( m_maxPtCut > 0 && pt > m_maxPtCut ) ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Pt " << aTrack.pt() << " failed cut" << endmsg;
      return false;
    }

    // track types
    if ( !checkTrackType( aTrack ) ) return false;

    // number of hits
    if ( m_minNVeloHits > 0 ) {
      auto numVelo = std::count_if( std::begin( aTrack.lhcbIDs() ), std::end( aTrack.lhcbIDs() ),
                                    []( const auto& id ) { return id.isVP(); } );
      if ( numVelo < m_minNVeloHits ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #velo " << numVelo << " failed cut" << endmsg;
        return false;
      }
    }

    if ( m_minNTHits > 0 ) {
      auto numT = std::count_if( std::begin( aTrack.lhcbIDs() ), std::end( aTrack.lhcbIDs() ),
                                 []( const auto& id ) { return id.isFT(); } );
      if ( numT < m_minNTHits ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #T hits " << numT << " failed cut" << endmsg;
        return false;
      }
    }

    if ( m_maxNVeloHoles < 99 || m_maxNTHoles < 99 || m_minNVeloLayers > 0 || m_minNVeloALayers > 0 ||
         m_minNVeloCLayers > 0 || m_minNTLayers > 0 || m_minNUTLayers > 0 || m_minNVeloOverlap > 0 ) {
      LHCb::HitPattern hitpattern( aTrack.lhcbIDs() );
      const int        numVeloHoles = hitpattern.numVeloHoles();
      if ( numVeloHoles > m_maxNVeloHoles ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #VeloHoles " << numVeloHoles << " failed cut" << endmsg;
        return false;
      }
      const int numVeloLayers = hitpattern.numVeloA() + hitpattern.numVeloC();
      if ( numVeloLayers < m_minNVeloLayers ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #VeloLayers " << numVeloLayers << " failed cut" << endmsg;
        return false;
      }
      const int numVeloALayers = hitpattern.numVeloA();
      if ( numVeloALayers < m_minNVeloALayers ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #VeloALayers " << numVeloALayers << " failed cut" << endmsg;
        return false;
      }
      const int numVeloCLayers = hitpattern.numVeloC();
      if ( numVeloCLayers < m_minNVeloCLayers ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #VeloCLayers " << numVeloCLayers << " failed cut" << endmsg;
        return false;
      }
      const int numVeloOverlap = std::min( hitpattern.numVeloA(), hitpattern.numVeloC() );
      if ( numVeloOverlap < m_minNVeloOverlap ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #VeloOverlap " << numVeloOverlap << " failed cut" << endmsg;
        return false;
      }
      const int numTHoles = hitpattern.numFTHoles();
      if ( numTHoles > m_maxNTHoles ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #THoles " << numTHoles << " failed cut" << endmsg;
        return false;
      }
      const int numTLayers = hitpattern.numFT();
      if ( numTLayers < m_minNTLayers ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #TLayers " << numTLayers << " failed cut" << endmsg;
        return false;
      }
      const int numUTLayers = hitpattern.numUT();
      if ( numUTLayers < m_minNUTLayers ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #UTLayers " << numUTLayers << " failed cut" << endmsg;
        return false;
      }
    }

    // eta
    const double eta = aTrack.pseudoRapidity();
    if ( eta < m_minEtaCut || eta > m_maxEtaCut ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #eta " << eta << " failed cut" << endmsg;
      return false;
    }

    // phi
    const double phi = aTrack.phi();
    if ( phi < m_minPhiCut || phi > m_maxPhiCut ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> #phi " << phi << " failed cut" << endmsg;
      return false;
    }

    // Clones
    const double cloneDist = aTrack.info( LHCb::Track::AdditionalInfo::CloneDist, 9e99 );
    if ( !m_acceptClones &&
         ( aTrack.checkFlag( LHCb::Track::Flags::Clone ) || cloneDist < m_minCloneCut || cloneDist > m_maxCloneCut ) ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Track failed clone rejection" << endmsg;
      return false;
    }

    if ( aTrack.likelihood() < m_minLikCut || aTrack.likelihood() > m_maxLikCut ) {
      if ( msgLevel( MSG::VERBOSE ) )
        verbose() << " -> Track Likelihood " << aTrack.likelihood() << " failed cut" << endmsg;
      return false;
    }

    // GhostProbability
    if ( aTrack.ghostProbability() < m_minGhostProb || aTrack.ghostProbability() > m_maxGhostProb ) {
      if ( msgLevel( MSG::VERBOSE ) )
        verbose() << " -> Track GhostProbability " << aTrack.ghostProbability() << " failed cut" << endmsg;
      return false;
    }

    // if get here track is selected !
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Track selected" << endmsg;

  } catch ( const GaudiException& excpt ) {
    // print the exception message as a warning and reject the track
    std::ostringstream mess;
    mess << "GaudiException caught " << excpt.message() << " " << excpt.tag() << " -> Track rejected";
    Warning( mess.str() ).ignore();
    return false;
  } catch ( const std::exception& excpt ) {
    // print the exception message as a warning and reject the track
    std::ostringstream mess;
    mess << "std::exception caught " << excpt.what() << " -> Track rejected";
    Warning( mess.str() ).ignore();
    return false;
  }

  // return OK
  return true;
}
