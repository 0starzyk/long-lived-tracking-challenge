/*****************************************************************************\
* (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "TrackCheckerBase.h"
#include "Event/Track.h"
#include "Map.h"

//=============================================================================
// Initialization. Check parameters
//=============================================================================
StatusCode TrackCheckerBase::initialize() {

  static const std::string histoDir = "Track/";
  if ( "" == histoTopDir() ) setHistoTopDir( histoDir );

  // Mandatory initialization of GaudiAlgorithm
  return GaudiHistoAlg::initialize().andThen( [&] {
    const TrackMaps::RecMap& theMap = TrackMaps::recDescription();
    m_recCat                        = theMap.find( m_selectionCriteria )->second;
  } );
}

const LHCb::MCParticle* TrackCheckerBase::mcTruth( const LHCb::Track& track, const LHCb::MCParticles& mcParts,
                                                   const LHCb::LinksByKey& links ) const {
  const LHCb::MCParticle* mcparticle{nullptr};
  links.applyToLinks( track.key(), [&mcparticle, &mcParts, this]( unsigned int, unsigned int mcPartKey, float ) {
    if ( !mcparticle ) {
      mcparticle = static_cast<const LHCb::MCParticle*>( mcParts.containedObject( mcPartKey ) );
      if ( mcparticle && mcparticle->particleID().threeCharge() == 0 ) {
        this->Warning( "Linker table for track contains pointer to particle with zero charge", StatusCode::SUCCESS, 0 )
            .ignore();
        mcparticle = nullptr;
      }
    }
  } );
  return mcparticle;
}

bool TrackCheckerBase::bAncestorWithReconstructibleDaughters( const LHCb::MCParticle* mcPart ) const {
  // loop back and see if there is a B in the history
  bool                    fromB  = false;
  const LHCb::MCParticle* mother = mcPart->mother();
  while ( mother && !fromB ) {
    fromB = mother->particleID().hasBottom() && ( mother->particleID().isMeson() || mother->particleID().isBaryon() );
    if ( fromB && !allDaughtersReconstructible( mother ) ) return false;
    mother = mother->mother();
  } // loop
  return fromB;
}

bool TrackCheckerBase::bAncestor( const LHCb::MCParticle* mcPart ) const {
  // loop back and see if there is a B in the history
  bool                    fromB  = false;
  const LHCb::MCParticle* mother = mcPart->mother();
  while ( mother && !fromB ) {
    fromB  = mother->particleID().hasBottom() && ( mother->particleID().isMeson() || mother->particleID().isBaryon() );
    mother = mother->mother();
  } // loop
  return fromB;
}

bool TrackCheckerBase::ksLambdaAncestor( const LHCb::MCParticle* mcPart ) const {
  // loop back and see if there is a B in the history
  bool                    fromKsL = false;
  const LHCb::MCParticle* mother  = mcPart->mother();
  while ( mother && !fromKsL ) {
    if ( abs( mother->particleID().pid() ) == 310 || abs( mother->particleID().pid() ) == 3122 ) fromKsL = true;
    mother = mother->mother();
  } // loop
  return fromKsL;
}

bool TrackCheckerBase::allDaughtersReconstructible( const LHCb::MCParticle* mcPart ) const {
  const SmartRefVector<LHCb::MCVertex>& vtx = mcPart->endVertices();

  for ( SmartRefVector<LHCb::MCVertex>::const_iterator i = vtx.begin(); i != vtx.end(); ++i ) {
    const SmartRefVector<LHCb::MCParticle>& ch = ( *i )->products();
    for ( SmartRefVector<LHCb::MCParticle>::const_iterator j = ch.begin(); j != ch.end(); ++j ) {

      if ( ( abs( ( *j )->particleID().pid() ) == 321 || abs( ( *j )->particleID().pid() ) == 211 ||
             abs( ( *j )->particleID().pid() ) == 13 || abs( ( *j )->particleID().pid() ) == 11 ||
             abs( ( *j )->particleID().pid() ) == 2212 ) ) {
        if ( !selected( *j ) && ( *j )->mother()->particleID().pid() != 22 &&
             ( *j )->mother()->particleID().pid() != -99000000 && ( *j )->mother()->particleID().pid() != 130 &&
             ( *j )->mother()->particleID().pid() != 310 && ( *j )->mother()->particleID().pid() != 3122 ) {
          return false;
        }
      } else if ( !allDaughtersReconstructible( *j ) )
        return false;
    }
  }

  return true;
}
