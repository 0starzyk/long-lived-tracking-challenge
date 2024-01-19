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
// from Gaudi
#include "Event/FTCluster.h"
#include "Event/MCParticle.h"
#include "Event/MCVertex.h"
#include "Event/MuonCoord.h"
#include "Event/Track.h"
#include "Event/UTCluster.h"
#include "Event/VPCluster.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "Linker/LinkedTo.h"
#include <optional>

/** @class TrackAssociator TrackAssociator.h
 *
 *  This algorithm computes the link between a Track and a MCParticle.
 *  The requirement is a match of both the VP and the Seed part of the
 *  Track. If there are not enough coordinates, the match is assumed so that
 *  a VP only or a Seed only are matched properly.
 *  The required fraction of hits is a jobOption 'FractionOK', default 0.70.
 *
 *  Adapted to the new Track Event Model using Linkers
 *  @author Edwin Bos
 *  @date   2005-11-14
 *
 *  Based on the Tr/TrFitAssociators package by :
 *  @author Olivier Callot
 *  @date   2002-07-01
 */

class TrackAssociator : public GaudiAlgorithm {
public:
  // Standard constructor
  TrackAssociator( const std::string& name, ISvcLocator* pSvcLocator );

  // Algorithm initialization
  StatusCode initialize() override;

  // Algorithm execution
  StatusCode execute() override;

private:
  // jobOptions
  std::string m_tracksInContainer; //< Name of the input Tracks container
  std::string m_linkerOutTable;    //< Name of the output Linker table
  double      m_fractionOK = 0;    //< minimal good matching fraction
  bool        m_decideUsingMuons;  //< use muon hits in link decision
};

using namespace LHCb;

DECLARE_COMPONENT( TrackAssociator )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TrackAssociator::TrackAssociator( const std::string& name, ISvcLocator* pSvcLocator )
    : GaudiAlgorithm( name, pSvcLocator ) {
  declareProperty( "TracksInContainer", m_tracksInContainer = TrackLocation::Default );
  declareProperty( "LinkerOutTable", m_linkerOutTable = "" );
  declareProperty( "FractionOK", m_fractionOK = 0.70 );
  declareProperty( "DecideUsingMuons", m_decideUsingMuons = false );
}

//=============================================================================
// Initialisation. Check parameters
//=============================================================================
StatusCode TrackAssociator::initialize() {

  // Mandatory initialization of GaudiAlgorithm
  StatusCode sc = GaudiAlgorithm::initialize();
  if ( sc.isFailure() ) { return sc; }

  // Set the path for the linker table Track - MCParticle
  if ( m_linkerOutTable == "" ) m_linkerOutTable = m_tracksInContainer;

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode TrackAssociator::execute() {

  // Vector containing the MCParticles which
  // have a Measurement of any type associated to them
  std::vector<const LHCb::MCParticle*> parts;

  // Number of VP Measurements assigned to the MCParticle which
  // has the same vector index in parts
  std::vector<double> nVP;

  // Number of UT Measurements assigned to the MCParticle which
  // has the same vector index in parts
  std::vector<double> nUT1;

  // Number of IT+OT Measurements assigned to the MCParticle which
  // has the same vector index in parts
  std::vector<double> nSeed;

  // Number of Muon Measurements assigned to the MCParticle which
  // has the same vector index in parts
  std::vector<double> nMuon;

  double nTotVP   = 0; // Total number of VP hits
  double nTotUT1  = 0; // Total number of UT (UT) hits
  double nTotSeed = 0; // Total number of IT+OT hits
  double nTotMuon = 0; // Total number of Muon hits

  // Retrieve the Tracks
  LHCb::Track::Range tracks = getIfExists<LHCb::Track::Range>( m_tracksInContainer );

  // Retrieve the MCParticles
  MCParticles* mcParts = get<MCParticles>( MCParticleLocation::Default );

  // Create the Linker table from Track to MCParticle
  // Linker table is stored in "Link/" + m_linkerOutTable
  // Sorted by decreasing weight, so first retrieved has highest weight
  auto* myLinker = LHCb::LinksByKey::fetchFromOrCreateOnTES<Track, LHCb::MCParticle>(
      evtSvc(), m_linkerOutTable, LHCb::LinksByKey::Order::decreasingWeight );
  myLinker->reset(); //== reset it, in case we work from a DST where associations are already there.

  // Get the linker table VpCluster => MCParticle

  LinkedTo<MCParticle> vpLink{
      get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::VPClusterLocation::Default ) )};

  // Get the linker table TTCluster (or UTCluster) => MCParticle
  LinkedTo<LHCb::MCParticle> utLink{
      get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::UTClusterLocation::UTClusters ) )};

  // Get the linker table FTCluster => MCParticle
  LinkedTo<LHCb::MCParticle> ftLink{
      get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::FTClusterLocation::Default ) )};

  std::optional<LinkedTo<LHCb::MCParticle>> muonLinks;
  // Get the linker table MuonCoord => MCParticle
  if ( m_decideUsingMuons ) {
    muonLinks.emplace( get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::MuonCoordLocation::MuonCoords ) ) );
  }

  // Loop over the Tracks
  for ( const Track* tr : tracks ) {
    nTotVP   = 0.;
    nTotUT1  = 0.;
    nTotSeed = 0.;
    nTotMuon = 0.;
    parts.clear();
    nVP.clear();
    nUT1.clear();
    nSeed.clear();
    nMuon.clear();

    //=============================================================================
    // Adjust the counters for this MCParticle, create one if needed
    //=============================================================================
    auto countMCPart = [&]( const MCParticle& part, double incVP, double incUT1, double incSeed, double incMuon ) {
      // Loop over the MCParticles vector to see whether part is already in parts
      auto i = std::find( parts.begin(), parts.end(), &part );
      if ( i != parts.end() ) {
        auto jj = i - parts.begin();
        nVP[jj] += incVP;
        nUT1[jj] += incUT1;
        nSeed[jj] += incSeed;
        nMuon[jj] += incMuon;
      } else {
        // If part an element of parts if it was not already and update counters
        parts.push_back( &part );
        nVP.push_back( incVP );
        nUT1.push_back( incUT1 );
        nSeed.push_back( incSeed );
        nMuon.push_back( incMuon );
      }
    };

    // Loop over the LHCbIDs of the Track
    static int nMeas;
    nMeas = 0;
    for ( const auto& id : tr->lhcbIDs() ) {
      if ( id.isVP() ) {
        ++nMeas;
        // Count number of VP hits
        nTotVP += 1.;
        auto mcParticles = vpLink.range( id.vpID() );
        if ( msgLevel( MSG::DEBUG ) && mcParticles.empty() ) {
          debug() << "No MCParticle linked with VPCluster " << id.vpID() << endmsg;
        }
        for ( const auto& mcParticle : mcParticles ) {
          if ( mcParts != mcParticle.parent() ) {
            if ( msgLevel( MSG::DEBUG ) ) debug() << " (other BX " << mcParticle.key() << ")" << endmsg;
          } else {
            countMCPart( mcParticle, 1., 0., 0., 0. );
          }
        }
        continue;
      }
      if ( id.isUT() ) { // UT
        ++nMeas;
        nTotUT1 += 1.;
        auto mcParticles = utLink.range( id.utID() );
        if ( msgLevel( MSG::DEBUG ) && mcParticles.empty() ) {
          debug() << "No MCParticle linked with UTCluster " << id.utID() << endmsg;
        }
        for ( const auto& mcParticle : mcParticles ) {
          if ( mcParts != mcParticle.parent() ) {
            if ( msgLevel( MSG::DEBUG ) ) debug() << " (other BX " << mcParticle.key() << ")" << endmsg;
          } else {
            countMCPart( mcParticle, 0., 1., 0., 0. );
          }
        }
        continue;
      }
      if ( m_decideUsingMuons ) {
        if ( id.isMuon() ) {
          ++nMeas;
          // Count number of Muon hits
          nTotMuon += 1.;
          auto mcParticles = muonLinks->range( id.muonID() );
          if ( msgLevel( MSG::DEBUG ) && mcParticles.empty() ) {
            debug() << "No MCParticle linked with MuonCoord " << id.muonID() << endmsg;
          }
          for ( const auto& mcParticle : mcParticles ) {
            if ( mcParts != mcParticle.parent() ) {
              debug() << " (other BX " << mcParticle.key() << ")" << endmsg;
            } else {
              countMCPart( mcParticle, 0., 0., 0., 1. );
            }
          }
        }
      }
    }

    //== Handling of decays between VP and T: Associate the hits also to the mother

    for ( unsigned int j1 = 0; parts.size() > j1; ++j1 ) {
      const MCVertex* vOrigin = parts[j1]->originVertex();
      if ( 0 == vOrigin ) continue;
      const MCParticle* mother = vOrigin->mother();
      if ( 0 == mother ) continue;
      if ( vOrigin->type() != MCVertex::MCVertexType::DecayVertex ) continue;
      for ( unsigned int j2 = 0; parts.size() > j2; ++j2 ) {
        if ( mother == parts[j2] ) {
          if ( msgLevel( MSG::DEBUG ) )
            debug() << "  *** Particle " << parts[j1]->key() << "[" << parts[j1]->particleID().pid() << "] (" << nVP[j1]
                    << "," << nUT1[j1] << "," << nSeed[j1] << ")"
                    << " is daughter at z=" << vOrigin->position().z() << " of " << parts[j2]->key() << "["
                    << parts[j2]->particleID().pid() << "] (" << nVP[j2] << "," << nUT1[j2] << "," << nSeed[j2] << ")"
                    << ". Merge hits to tag both." << endmsg;

          //== Daughter hits are added to mother.
          nVP[j2] += nVP[j1];
          nUT1[j2] += nUT1[j1];
          nSeed[j2] += nSeed[j1];

          //== Mother hits overwrite Daughter hits
          nVP[j1]   = nVP[j2];
          nUT1[j1]  = nUT1[j2];
          nSeed[j1] = nSeed[j2];
        }
      }
    }

    // Add parent muon hits to daughter MCParticle
    if ( m_decideUsingMuons ) {
      for ( unsigned j1 = 0; parts.size() > j1; ++j1 ) {
        if ( nMuon[j1] < m_fractionOK * nTotMuon ) { continue; }
        const MCVertex* vOrigin = parts[j1]->originVertex();
        while ( 0 != vOrigin ) {
          const MCParticle* mother = vOrigin->mother();
          if ( 0 == mother ) break;
          for ( unsigned j2 = 0; parts.size() > j2; ++j2 ) {
            if ( nSeed[j2] < m_fractionOK * nTotSeed ) { continue; }
            if ( mother == parts[j2] ) {
              nMuon[j2] += nMuon[j1];
              nMuon[j1] = nMuon[j2];
            }
          }
          vOrigin = mother->originVertex();
        }
      }
    }

    // For each MCParticle with a Measurement associated,
    //  Relate the Track to the MCParticles matching the below criteria,
    //  using as weight the ratio of ( total # Measurements associated
    //                                 to this MCParticle )
    //                           and ( total # Measurements )
    //   if total # VP hits > 2,
    //   or if 0.70 <= ( # VP hits associated to the MCParticle /
    //                   total # VP hits )
    //  AND
    //   if total # IT+OT hits > 2,
    //   or if 0.70 <= ( # IT+OT hits associated to the MCParticle /
    //                   total # IT+OT hits )
    //  AND
    //   if # associated UT hits > ( total # UT hits - 2 ) ??
    //   or if total # VP hits > 2 AND total # IT+OT hits > 2
    // When taking Muon hits into account, also:
    //  AND
    //   if 0.70 <= ( # Muon hits associated to the MCParticle /
    //                total # Muon hits )
    double ratio;
    for ( unsigned int jj = 0; parts.size() > jj; ++jj ) {
      bool vpOK = true;
      if ( 2 < nTotVP ) {
        vpOK  = false;
        ratio = nVP[jj] / nTotVP;
        if ( m_fractionOK <= ratio ) { vpOK = true; }
      }
      bool seedOK = true;
      if ( 2 < nTotSeed ) {
        seedOK = false;
        ratio  = nSeed[jj] / nTotSeed;
        if ( m_fractionOK <= ratio ) { seedOK = true; }
      }
      bool UT1OK = nUT1[jj] > nTotUT1 - 2;
      if ( 2 < nTotVP && 2 < nTotSeed ) { UT1OK = true; }
      bool muonOK = true;
      if ( m_decideUsingMuons ) {
        ratio = nMuon[jj] / nTotMuon;
        if ( m_fractionOK > ratio ) { muonOK = false; }
      }

      if ( msgLevel( MSG::VERBOSE ) )
        verbose() << " MC " << parts[jj]->key() << " VP " << nVP[jj] << "/" << nTotVP << " UT1 " << nUT1[jj] << "/"
                  << nTotUT1 << " Seed " << nSeed[jj] << "/" << nTotSeed << endmsg;

      //=== Decision. Use UT1. Fill Linker
      if ( vpOK && seedOK && UT1OK && muonOK ) {
        double muons = 0.;
        if ( m_decideUsingMuons ) { muons = nMuon[jj]; }
        ratio = ( nVP[jj] + nSeed[jj] + nUT1[jj] + muons ) / nMeas;
        myLinker->link( tr, parts[jj], ratio );
        if ( msgLevel( MSG::DEBUG ) )
          debug() << "Track: " << tr->key() << " was associated with ratio: " << ratio << endmsg;
      }
    }
  } // End loop over Tracks

  return StatusCode::SUCCESS;
}
