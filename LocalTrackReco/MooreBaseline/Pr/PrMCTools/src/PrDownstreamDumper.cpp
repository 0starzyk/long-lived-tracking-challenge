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
// Implementation file for class : PrDownstreamDumper
//
// 2021-01 : Grzegorz Golaszewski - adapting
// 2019-03 : Artur Kucia
// Based on PrTrackerDumper2 by Renato Quagliani and Giulia Tuci.
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
#include "Linker/LinkedFrom.h"

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"

#include "PrKernel/PrFTHitHandler.h"
#include "PrKernel/PrHit.h"
#include "PrKernel/UTHit.h"
#include "PrKernel/UTHitHandler.h"

#include "MCInterfaces/IIdealStateCreator.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

// ROOT
#include "TROOT.h"
#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include <boost/filesystem.hpp>
#include <utility>

/** @class PrDownsteramDumper PrDownstreamDumper.cpp
 */

class PrDownstreamDumper
    : public LHCb::Algorithm::Consumer<void( const LHCb::MCParticles&, const UT::HitHandler&, const LHCb::ODIN&,
                                             const LHCb::LinksByKey&, DetectorElement const& ),
                                       LHCb::DetDesc::usesConditions<DetectorElement>> {
public:
  /// Standard constructor
  PrDownstreamDumper( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override;
  void       operator()( const LHCb::MCParticles& MCParticles, const UT::HitHandler& utHits, const LHCb::ODIN& odin,
                   const LHCb::LinksByKey& links, DetectorElement const& ) const override;
  const ToolHandle<ITrackExtrapolator> extrapolator() const { return m_extrapolator; }
  const ToolHandle<IIdealStateCreator> stateCreator() const { return m_stateCreator; }

private:
  Gaudi::Property<std::string> m_seed_location{this, "TrackLocation", LHCb::Event::v2::TrackLocation::Seed,
                                               "Container in which are stored reconstructed T-seeds linked to MC "
                                               "particles"};
  Gaudi::Property<std::string> m_output_directory{this, "OutputDirectory", ".", "Directory to store output hTuples"};

  ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackMasterExtrapolatro"};

  ToolHandle<IIdealStateCreator> m_stateCreator{this, "IdealStateCreator", "IdealStateCreator"};

  typedef std::map<const LHCb::MCParticle*, std::vector<UT::Hit>> ParticleToUTHits;

  std::tuple<PrDownstreamDumper::ParticleToUTHits, std::vector<UT::Hit>>
  SplitUTHits( const UT::HitHandler&, const InputLinks<ContainedObject, LHCb::MCParticle>& ) const;
};

DECLARE_COMPONENT( PrDownstreamDumper )

namespace {
  namespace fs = boost::filesystem;
} // namespace

PrDownstreamDumper::PrDownstreamDumper( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator,
                {KeyValue{"MCParticlesLocation", LHCb::MCParticleLocation::Default},
                 KeyValue{"UTHitsLocation", UTInfo::HitLocation}, KeyValue{"ODINLocation", LHCb::ODINLocation::Default},
                 KeyValue{"LinkerLocation", Links::location( "Pr/LHCbID" )},
                 KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}} ) {}

// initialize creates a .root file with tree for saving hits data.
StatusCode PrDownstreamDumper::initialize() {
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

  gROOT->ProcessLine( "#include <vector>" );
  gROOT->ProcessLine( "#include <tuple>" );
  return sc;
}

//=============================================================================
// operator() Iterates over all FT and UT hits and over all MC Particles,
// which left hits in UT or FT. Stores their data into nTuple.
//=============================================================================
void PrDownstreamDumper::operator()( const LHCb::MCParticles& MCParticles, const UT::HitHandler& prUTHitHandler,
                                     const LHCb::ODIN& odin, const LHCb::LinksByKey& links,
                                     DetectorElement const& lhcb ) const {
  std::string filename = ( "DownstreamDumper_runNb_" + std::to_string( odin.runNumber() ) + "_evtNb_" +
                           std::to_string( odin.eventNumber() ) + ".root" );
  TFile       file( filename.c_str(), "UPDATE" );

  int    eventID                = 0;
  bool   particle_fullInfo      = false;
  bool   particle_hasSciFi      = false;
  bool   particle_hasUT         = false;
  bool   particle_hasVelo       = false;
  bool   particle_isDown        = false;
  bool   particle_isDown_noVelo = false;
  bool   particle_isLong        = false;
  bool   particle_isLong_andUT  = false;
  bool   particle_isElectron    = false;
  double particle_eta           = 0;
  double particle_ovtx_x        = 0;
  double particle_ovtx_y        = 0;
  double particle_ovtx_z        = 0;
  double particle_evtx_x        = 0;
  double particle_evtx_y        = 0;
  double particle_evtx_z        = 0;
  double particle_p             = 0;
  double particle_phi           = 0;
  double particle_pt            = 0;
  double particle_px            = 0;
  double particle_py            = 0;
  double particle_pz            = 0;
  int    particle_pid           = 0;
  int    particle_key           = 0;

  double track_charge           = 0;
  double track_chi2             = 0;
  double track_chi2PerDoF       = 0;
  double track_nLHCbIDs         = 0;
  double track_p                = 0;
  double track_phi              = 0;
  double track_position_phi     = 0;
  double track_position_r       = 0;
  double track_position_x       = 0;
  double track_position_y       = 0;
  double track_position_z       = 0;
  double track_pt               = 0;
  double track_tx               = 0;
  double track_ty               = 0;
  double track_ghostProbability = 0;
  double track_pseudoRapidity   = 0;

  std::vector<float> track_stateVector;
  std::vector<float> track_stateCovariance;
  std::vector<float> track_stateZ;

  std::vector<float> track_residue;

  std::vector<float> UT_cos;
  std::vector<float> UT_cosT;
  std::vector<float> UT_dxDy;
  std::vector<float> UT_sinT;
  std::vector<float> UT_tanT;

  std::vector<float> UT_weight;
  std::vector<float> UT_xAtYEq0;
  std::vector<float> UT_xAtYMid;
  std::vector<float> UT_xT;
  std::vector<float> UT_zAtYEq0;

  std::vector<int> UT_lhcbID;
  std::vector<int> UT_particle_key;
  std::vector<int> UT_planeCode;

  auto tree = TTree( "Particle_Data", "Particle_Data" );
  tree.Branch( "eventID", &eventID );

  tree.Branch( "particle_key", &particle_key );
  tree.Branch( "particle_eta", &particle_eta );
  tree.Branch( "particle_fullInfo", &particle_fullInfo );
  tree.Branch( "particle_hasScifi", &particle_hasSciFi );
  tree.Branch( "particle_hasUT", &particle_hasUT );
  tree.Branch( "particle_hasVelo", &particle_hasVelo );
  tree.Branch( "particle_isDown", &particle_isDown );
  tree.Branch( "particle_isDown_noVelo", &particle_isDown_noVelo );
  tree.Branch( "particle_isLong", &particle_isLong );
  tree.Branch( "particle_isLong_andUT", &particle_isLong_andUT );
  tree.Branch( "partiole_isElectron", &particle_isElectron );
  tree.Branch( "particle_ovtx_x", &particle_ovtx_x );
  tree.Branch( "particle_ovtx_y", &particle_ovtx_y );
  tree.Branch( "particle_ovtx_z", &particle_ovtx_z );
  tree.Branch( "particle_evtx_x", &particle_evtx_x );
  tree.Branch( "particle_evtx_y", &particle_evtx_y );
  tree.Branch( "particle_evtx_z", &particle_evtx_z );
  tree.Branch( "particle_p", &particle_p );
  tree.Branch( "particle_phi", &particle_phi );
  tree.Branch( "particle_pid", &particle_pid );
  tree.Branch( "particle_pt", &particle_pt );
  tree.Branch( "particle_px", &particle_px );
  tree.Branch( "particle_py", &particle_py );
  tree.Branch( "particle_pz", &particle_pz );

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
  tree.Branch( "track_ghostProbability", &track_ghostProbability );

  tree.Branch( "track_stateVector", &track_stateVector );
  tree.Branch( "track_stateCovariance", &track_stateCovariance );
  tree.Branch( "track_stateZ", &track_stateZ );
  tree.Branch( "track_residue", &track_residue );

  tree.Branch( "UT_cos", &UT_cos );
  tree.Branch( "UT_cosT", &UT_cosT );
  tree.Branch( "UT_dxDy", &UT_dxDy );
  tree.Branch( "UT_lhcbID", &UT_lhcbID );
  tree.Branch( "UT_planeCode", &UT_planeCode );
  tree.Branch( "UT_sinT", &UT_sinT );
  tree.Branch( "UT_tanT", &UT_tanT );
  tree.Branch( "UT_weight", &UT_weight );
  tree.Branch( "UT_xAtYEq0", &UT_xAtYEq0 );
  tree.Branch( "UT_xAtYMid", &UT_xAtYMid );
  tree.Branch( "UT_xT", &UT_xT );
  tree.Branch( "UT_zAtYEq0", &UT_zAtYEq0 );

  eventID = odin.eventNumber();

  InputLinks<ContainedObject, LHCb::MCParticle> HitMCParticleLinks( links );
  auto seedLinks = SmartDataPtr<LHCb::LinksByKey>{evtSvc(), LHCb::LinksByKey::linkerName( m_seed_location )};
  auto [mcparticles_to_ut_hits, non_assoc_ut_hits] = SplitUTHits( prUTHitHandler, links );

  MCTrackInfo track_info = make_MCTrackInfo( evtSvc(), msgSvc() );

  for ( const auto* particle : MCParticles ) {
    const LHCb::Track* seed_track = LinkedFrom<LHCb::Track>{seedLinks}.range( particle ).try_front();

    UT_cos.clear();
    UT_cosT.clear();
    UT_dxDy.clear();
    UT_sinT.clear();
    UT_tanT.clear();
    UT_weight.clear();
    UT_xAtYEq0.clear();
    UT_xAtYMid.clear();
    UT_xT.clear();
    UT_zAtYEq0.clear();
    UT_lhcbID.clear();
    UT_planeCode.clear();

    track_stateVector.clear();
    track_stateCovariance.clear();
    track_stateZ.clear();
    track_residue.clear();

    if ( seed_track ) {

      for ( auto const& [UTparticle, UThits] : mcparticles_to_ut_hits ) {
        if ( particle->key() == UTparticle->key() ) {

          for ( auto const& hit : UThits ) {
            UT_cos.push_back( hit.cos() );
            UT_cosT.push_back( hit.cosT() );
            UT_dxDy.push_back( hit.dxDy() );
            UT_sinT.push_back( hit.sinT() );
            UT_tanT.push_back( hit.tanT() );
            UT_weight.push_back( hit.weight() );
            UT_xAtYEq0.push_back( hit.xAtYEq0() );
            UT_xAtYMid.push_back( hit.xAtYMid() );
            UT_xT.push_back( hit.xT() );
            UT_zAtYEq0.push_back( hit.zAtYEq0() );
            UT_lhcbID.push_back( hit.chanID().channelID() );
            UT_planeCode.push_back( hit.planeCode() );

            LHCb::State state;
            LHCb::State trueState;
            StatusCode  extrapolator_sc =
                extrapolator()->propagate( *seed_track, hit.zAtYEq0(), state, *lhcb.geometry() );
            StatusCode statecreator_sc =
                stateCreator()->createState( particle, hit.zAtYEq0(), trueState, *lhcb.geometry() );

            if ( extrapolator_sc.isFailure() ) continue;

            track_stateZ.push_back( hit.zAtYEq0() );
            for ( int i = 0; i < 5; ++i ) {
              track_stateVector.push_back( state.stateVector().At( i ) );

              if ( statecreator_sc.isSuccess() )
                track_residue.push_back( state.stateVector().At( i ) - trueState.stateVector().At( i ) );

              for ( int j = 0; j < 5; ++j ) { track_stateCovariance.push_back( state.covariance().At( i, j ) ); }
            }
          }
        }
      }

      particle_fullInfo      = track_info.fullInfo( particle );
      particle_hasSciFi      = track_info.hasT( particle );
      particle_hasUT         = track_info.hasUT( particle );
      particle_hasVelo       = track_info.hasVelo( particle );
      particle_isDown        = particle_hasSciFi && particle_hasUT;
      particle_isDown_noVelo = particle_isDown && !( particle_hasVelo );
      particle_isLong        = particle_hasSciFi && particle_hasVelo;
      particle_isLong_andUT  = particle_isLong && particle_hasUT;
      particle_isElectron    = std::abs( particle->particleID().pid() ) == 11;

      particle_p   = particle->p();
      particle_px  = particle->momentum().Px();
      particle_py  = particle->momentum().Py();
      particle_pz  = particle->momentum().Pz();
      particle_pt  = particle->pt();
      particle_eta = particle->momentum().Eta();
      particle_phi = particle->momentum().phi();
      particle_pid = particle->particleID().pid();

      particle_key = particle->key();

      if ( nullptr != particle->originVertex() ) {
        particle_ovtx_x = particle->originVertex()->position().x();
        particle_ovtx_y = particle->originVertex()->position().y();
        particle_ovtx_z = particle->originVertex()->position().z();
      }
      // Fill end vertex
      if ( 1 == particle->endVertices().size() ) {
        particle_evtx_x = particle->endVertices().at( 0 )->position().x();
        particle_evtx_y = particle->endVertices().at( 0 )->position().y();
        particle_evtx_z = particle->endVertices().at( 0 )->position().z();
      } else if ( 1 < particle->endVertices().size() ) {
        auto minV =
            std::min_element( particle->endVertices().begin(), particle->endVertices().end(),
                              []( auto& v1, auto& v2 ) { return ( *v1 ).position().z() < ( *v2 ).position().z(); } );
        particle_evtx_x = ( *minV )->position().x();
        particle_evtx_y = ( *minV )->position().y();
        particle_evtx_z = ( *minV )->position().z();
      }

      track_charge     = seed_track->charge();
      track_chi2       = seed_track->chi2();
      track_chi2PerDoF = seed_track->chi2PerDoF();
      track_nLHCbIDs   = seed_track->nLHCbIDs();
      track_p          = seed_track->p();
      track_phi        = seed_track->phi();
      track_tx         = seed_track->slopes().x();
      track_ty         = seed_track->slopes().y();

      track_ghostProbability = seed_track->ghostProbability();
      track_pseudoRapidity   = seed_track->pseudoRapidity();
      track_position_phi     = seed_track->position().phi();
      track_position_r       = std::sqrt( seed_track->position().x() * seed_track->position().x() +
                                    seed_track->position().y() * seed_track->position().y() );
      track_position_x       = seed_track->position().x();
      track_position_y       = seed_track->position().y();
      track_position_z       = seed_track->position().z();
      track_pt               = seed_track->pt();

      LHCb::State st;
      StatusCode  statecreator_sc =
          stateCreator()->createState( particle, seed_track->firstState().z(), st, *lhcb.geometry() );

      track_stateZ.push_back( seed_track->firstState().z() );
      for ( int i = 0; i < 5; ++i ) {
        track_stateVector.push_back( seed_track->firstState().stateVector().At( i ) );
        if ( statecreator_sc.isSuccess() )
          track_residue.push_back( seed_track->firstState().stateVector().At( i ) - st.stateVector().At( i ) );
        for ( int j = 0; j < 5; ++j ) {
          track_stateCovariance.push_back( seed_track->firstState().covariance().At( i, j ) );
        }
      }
      tree.Fill();
    }
  }
  file.Write();
  file.Close();
}

// SplitUTHits splits reads an input UT hit handler and returns a map of MC
// particles matched with hits and a vector of hits that were not matched with
// any MC particle
std::tuple<PrDownstreamDumper::ParticleToUTHits, std::vector<UT::Hit>>
PrDownstreamDumper::SplitUTHits( const UT::HitHandler&                                prUTHitHandler,
                                 const InputLinks<ContainedObject, LHCb::MCParticle>& HitMCParticleLinks ) const {
  ParticleToUTHits     particle_to_ut_hits;
  std::vector<UT::Hit> non_assoc_ut_hits;
  for ( int iStation = 1; iStation < 3; ++iStation ) {
    for ( int iLayer = 1; iLayer < 3; ++iLayer ) {
      for ( int iRegion = 1; iRegion < 4; ++iRegion ) {
        for ( int iSector = 1; iSector < 99; ++iSector ) {
          for ( const auto& hit : prUTHitHandler.hits( iStation, iLayer, iRegion, iSector ) ) {
            LHCb::LHCbID lhcbid               = hit.lhcbID();
            auto         mcparticlesrelations = HitMCParticleLinks.from( lhcbid.lhcbID() );
            if ( mcparticlesrelations.empty() ) { non_assoc_ut_hits.push_back( hit ); }
            for ( const auto& mcp : mcparticlesrelations ) {
              auto MCP = mcp.to();
              particle_to_ut_hits[MCP].push_back( hit );
            }
          }
        }
      }
    }
  }
  return {particle_to_ut_hits, non_assoc_ut_hits};
}
