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
// Implementation file for class : PrGhostDumper
//
// 2021-02 : Grzegorz Golaszewski
// Based on PrDownstreamDumper
//-----------------------------------------------------------------------------

// Include files
#include "Associators/Associators.h"
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/MCVertex.h"
#include "Event/ODIN.h"
#include "Event/Track.h"
#include "Event/VPLightCluster.h"
#include "LHCbAlgs/Consumer.h"
#include "Linker/LinkedTo.h"

#include "PrKernel/PrFTHitHandler.h"
#include "PrKernel/PrHit.h"
#include "PrKernel/UTHit.h"

// ROOT
#include "TROOT.h"
#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include <boost/filesystem.hpp>
#include <utility>

class PrGhostDumper
    : public LHCb::Algorithm::Consumer<void( const LHCb::Tracks&, const LHCb::ODIN&, const LHCb::LinksByKey& )> {
public:
  using Track = LHCb::Track;
  /// Standard constructor
  PrGhostDumper( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override;
  void operator()( const LHCb::Tracks& SeedTracks, const LHCb::ODIN& odin, const LHCb::LinksByKey& ) const override;

private:
  Gaudi::Property<std::string> m_seed_location{this, "TrackLocation", LHCb::TrackLocation::Seed,
                                               "Container in which are stored reconstructed T-seeds linked to MC "
                                               "particles"};
  Gaudi::Property<std::string> m_output_directory{this, "OutputDirectory", ".", "Directory to store output hTuples"};
};

// local
namespace {

  namespace fs = boost::filesystem;
} // namespace

// Declaration of the Algorithm Factory

DECLARE_COMPONENT( PrGhostDumper )

PrGhostDumper::PrGhostDumper( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator,
                {KeyValue{"SeedTrackLocation", LHCb::Event::TrackLocation::Default},
                 KeyValue{"ODINLocation", LHCb::ODINLocation::Default},
                 KeyValue{"SeedLinkerLocation", LHCb::LinksByKey::linkerName( LHCb::TrackLocation::Default )}} ) {}

StatusCode PrGhostDumper::initialize() {
  info() << "Initizaling" << endmsg;
  StatusCode sc = Consumer::initialize();
  if ( sc.isFailure() ) return sc;

  auto dir = fs::path{m_output_directory.value()};

  if ( !fs::exists( dir ) ) {
    boost::system::error_code ec;
    bool                      success = fs::create_directories( dir, ec );
    success &= !ec;
    if ( !success ) {
      error() << "Failed to create directory " << dir.string() << ": " << ec.message() << endmsg;
      return StatusCode::FAILURE;
    }
  }

  return sc;
}

//=============================================================================
// operator() Iterates over all FT and UT hits and over all MC Particles,
// which left hits in UT or FT. Stores their data into nTuple.
//=============================================================================
void PrGhostDumper::operator()( const LHCb::Tracks& SeedTracks, const LHCb::ODIN& odin,
                                const LHCb::LinksByKey& seed_links ) const {
  std::string filename = ( "GhostDumper_runNb_" + std::to_string( odin.runNumber() ) + "_evtNb_" +
                           std::to_string( odin.eventNumber() ) + ".root" );
  TFile       file( filename.c_str(), "UPDATE" );

  int eventID = 0;

  double track_charge         = 0;
  double track_chi2           = 0;
  double track_chi2PerDoF     = 0;
  double track_nLHCbIDs       = 0;
  double track_p              = 0;
  double track_phi            = 0;
  double track_position_phi   = 0;
  double track_position_r     = 0;
  double track_position_x     = 0;
  double track_position_y     = 0;
  double track_position_z     = 0;
  double track_pt             = 0;
  double track_tx             = 0;
  double track_ty             = 0;
  double track_pseudoRapidity = 0;

  auto tree = TTree( "Particle_Data", "Particle_Data" );

  tree.Branch( "eventID", &eventID );
  tree.Branch( "track_charge", &track_charge );
  tree.Branch( "track_chi2", &track_chi2 );
  tree.Branch( "track_chi2PerDoF", &track_chi2PerDoF );
  tree.Branch( "track_nLHCbIDs", &track_nLHCbIDs );
  tree.Branch( "track_p", &track_p );
  tree.Branch( "track_phi", &track_phi );
  tree.Branch( "track_position_phi", &track_position_phi );
  tree.Branch( "track_position_r", &track_position_r );
  tree.Branch( "track_position_x", &track_position_x );
  tree.Branch( "track_position_y", &track_position_y );
  tree.Branch( "track_position_z", &track_position_z );
  tree.Branch( "track_pt", &track_pt );
  tree.Branch( "track_tx", &track_tx );
  tree.Branch( "track_ty", &track_ty );
  tree.Branch( "track_pseudoRapidity", &track_pseudoRapidity );

  eventID = odin.eventNumber();

  for ( const auto* track : SeedTracks ) {
    const LHCb::MCParticle* mcSeedPart = LinkedTo<LHCb::MCParticle>{&seed_links}.range( track->key() ).try_front();
    if ( !mcSeedPart && track && track->position().z() > 7400 ) {
      track_charge = track->charge();
      track_chi2   = track->chi2();

      track_chi2PerDoF = track->chi2PerDoF();

      track_nLHCbIDs = track->nLHCbIDs();
      track_p        = track->p();
      track_phi      = track->phi();
      track_tx       = track->slopes().x();
      track_ty       = track->slopes().y();

      track_pseudoRapidity = track->pseudoRapidity();
      track_position_phi   = track->position().phi();
      track_position_r =
          std::sqrt( track->position().x() * track->position().x() + track->position().y() * track->position().y() );
      track_position_x = track->position().x();
      track_position_y = track->position().y();
      track_position_z = track->position().z();
      track_pt         = track->pt();

      tree.Fill();
    }
  }

  file.Write();
  file.Close();
}
