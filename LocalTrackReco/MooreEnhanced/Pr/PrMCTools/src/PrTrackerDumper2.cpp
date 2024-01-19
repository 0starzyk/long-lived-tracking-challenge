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
#include "Associators/Associators.h"
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/MCVertex.h"
#include "Event/ODIN.h"
#include "Event/Track.h"
#include "Event/VPLightCluster.h"
#include "LHCbAlgs/Consumer.h"
#include "Linker/LinkedFrom.h"
#include "Linker/LinkedTo.h"
#include "PrKernel/PrFTHitHandler.h"
#include "PrKernel/PrHit.h"
#include "PrKernel/UTHit.h"
#include "PrKernel/UTHitHandler.h"

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <mutex>
#include <utility>

//-----------------------------------------------------------------------------
// Implementation file for class : PrTrackerDumper2
//
// 2017-11-06 : Renato Quagliani
// 2018-11-21 : Giulia Tuci
// Based on PrTrackerDumper. Dumped also information of reconstructed tracks associated to a MC particle

// Branches produced regarding the MCParticle information:
/*
  fullInfo   : has full info MCParticle
  hasSciFi   : reconstructible in SciFi
  hasUT      : reconstructible in UT
  hasVelo    : reconstructible in Velo
  isDown     : reconstructible in UT and SciFi
  isDown_noVelo : reconstructible in UT and SciFi but not in Velo
  isLong        : reconstructible in Velo and SciFi
  isLong_andUT  : reconstructible in Velo and SciFi and UT
  p             : track momentum (cut for p>0 & at least having 1 hit in one of sub-detector)
 [among MCPartcile there are also the intermediate particles]
  pt            : track trsansverse momentum
  pid           : PID of the particle (can distinguish among muon/electrons/pion/kaons/protons etc...)
  eta           : track pseudorapidity
  ovtx_x        : track origin position X
  ovtx_y        : track origin position Y
  ovtx_z        : track origin position Z
  fromBeautyDecay : the track belongs to a decay chain with a b-quark hadron
  fromCharmDecay  : the track belongs to a decay chain with a c-quark hadron
  fromStrangeDecay : the track belongs to a decay chain with a s-quark hadron
  DecayOriginMother_pid : it store the PID of the head particle in the decay chain if found ,
  you can filter based on the simulated sample. If for instance you run over Bs->PhiPhi,
  you can filter the 4 kaons among all tracks requiring Bs PID for this variable
*/
/*  VELO related part
  "nVeloHits" : Number of VeloHits associated to the MCParticle
  Velo_x       : vector of x position for Velo  hits (size = nVeloHits)
  Velo_y       : vector of y position for Velo  hits (size = nVeloHits)
  Velo_z       : vector of z position for Velo  hits (size = nVeloHits)
  Velo_Module  : vector of ModuleID Velo  hits (size = nVeloHits)
  Velo_Sensor  : vector of SensorID Velo  hits (size = nVeloHits)
  Velo_Station : vector of StationID Velo  hits (size = nVeloHits)
  Velo_lhcbID  : vector of lhcbID Velo  hits (size = nVeloHits)
*/
/*  SciFi related part
  nFTHits   : Number of FTHits associated to the MCParticle
  FT_x      : vector of x(y=0) position for SciFi
  FT_z      : vector of z(y=0) position for SciFi hits
  FT_w      : vector of weight error   for SciFi hits
  FT_dxdy   : vector of slopes dxdy for SciFi hits
  FT_YMin   : vector of yMin  for SciFi hits
  FT_YMax   : vector of yMax  for SciFi hits
  FT_hitPlaneCode : vector of planeCode  for SciFi hits
  FT_hitzone      : vector of hitzone (up/down)  for SciFi hits
  FT_lhcbID       : vector of lhcbID  for SciFi hits
*/
/* UT related part
   nUTHits   : Number of UTHits associated to the MCParticle
  //---- see private members of UT:Hit in PrKernel package
  UT_cos
  UT_cosT
  UT_dxDy
  UT_lhcbID
  UT_planeCode
  UT_sinT
  UT_size
  UT_tanT
  UT_weight
  UT_xAtYEq0
  UT_xAtYMid
  UT_xMax
  UT_xMin
  UT_xT
  UT_yBegin
  UT_yEnd
  UT_yMax
  UT_yMid
  UT_yMin
  UT_zAtYEq0
*/
//-----------------------------------------------------------------------------
namespace {

  int computeNumberUTHits( const UT::HitHandler& prUTHitHandler ) {
    int nbHits = 0;
    for ( int iStation = 1; iStation < 3; ++iStation ) {
      for ( int iLayer = 1; iLayer < 3; ++iLayer ) {
        for ( int iRegion = 1; iRegion < 4; ++iRegion ) {
          for ( int iSector = 1; iSector < 99; ++iSector ) {
            nbHits += prUTHitHandler.hits( iStation, iLayer, iRegion, iSector ).size();
          }
        }
      }
    }
    return nbHits;
  }
} // namespace
/** @class PrTrackerDumper PrTrackerDumper2.h
 *  Algorithm to store tracking hits to a single root file
 *  (PrTrackerDumper stores events in separate files)
 *
 *  @author Renato Quagliani
 *  @date   2017-11-06
 */
/*

*/

class PrTrackerDumper2
    : public LHCb::Algorithm::Consumer<void(
          const LHCb::MCParticles&, const LHCb::VPLightClusters&, const PrFTHitHandler<PrHit>&, const UT::HitHandler&,
          const LHCb::ODIN&, const LHCb::LinksByKey&, const LHCb::MCProperty&, const LHCb::LinksByKey& )> {
public:
  PrTrackerDumper2( const std::string& name, ISvcLocator* pSvcLocator );
  StatusCode initialize() override;
  StatusCode finalize() override;
  void       operator()( const LHCb::MCParticles& MCParticles, const LHCb::VPLightClusters& VPClusters,
                   const PrFTHitHandler<PrHit>& ftHits, const UT::HitHandler& utHits, const LHCb::ODIN& odin,
                   const LHCb::LinksByKey& links, const LHCb::MCProperty&,
                   const LHCb::LinksByKey& tracklinks ) const override;

private:
  std::unique_ptr<TFile> file;
  std::unique_ptr<TTree> tree;

  mutable std::mutex m_mutex;
  mutable int        eventID = 0;

  mutable bool   fullInfo      = false;
  mutable bool   hasSciFi      = false;
  mutable bool   hasUT         = false;
  mutable bool   hasVelo       = false;
  mutable bool   isDown        = false;
  mutable bool   isDown_noVelo = false;
  mutable bool   isLong        = false;
  mutable bool   isLong_andUT  = false;
  mutable double p             = 0;
  mutable double px = 0, py = 0, pz = 0;
  mutable double pt  = 0;
  mutable double eta = 0, phi = 0;
  mutable double eta_track = 0, phi_track = 0, chi2_track = 0;
  // vertex origin of the particle
  mutable double                    ovtx_x                = 0;
  mutable double                    ovtx_y                = 0;
  mutable double                    ovtx_z                = 0;
  mutable int                       pid                   = 0;
  mutable bool                      fromBeautyDecay       = false;
  mutable bool                      fromCharmDecay        = false;
  mutable bool                      fromStrangeDecay      = false;
  mutable int                       DecayOriginMother_pid = 0;
  mutable int                       key                   = 0;
  mutable int                       nVeloHits_track       = 0;
  mutable std::vector<float>        Velo_x_track;
  mutable std::vector<float>        Velo_y_track;
  mutable std::vector<float>        Velo_z_track;
  mutable std::vector<int>          Velo_Module_track;
  mutable std::vector<int>          Velo_Sensor_track;
  mutable std::vector<int>          Velo_Station_track;
  mutable std::vector<unsigned int> Velo_lhcbID_track;
  mutable std::vector<unsigned int> Velo_index_track;

  mutable int                       nVeloHits = 0;
  mutable std::vector<float>        Velo_x;
  mutable std::vector<float>        Velo_y;
  mutable std::vector<float>        Velo_z;
  mutable std::vector<int>          Velo_Module;
  mutable std::vector<int>          Velo_Sensor;
  mutable std::vector<int>          Velo_Station;
  mutable std::vector<unsigned int> Velo_lhcbID;
  mutable std::vector<unsigned int> Velo_index;

  mutable std::vector<float>        FT_hitx;
  mutable std::vector<float>        FT_hitz;
  mutable std::vector<float>        FT_hitw;
  mutable std::vector<float>        FT_hitDXDY;
  mutable std::vector<float>        FT_hitYMin;
  mutable std::vector<float>        FT_hitYMax;
  mutable std::vector<int>          FT_hitPlaneCode;
  mutable std::vector<int>          FT_hitzone;
  mutable std::vector<unsigned int> FT_lhcbID;

  mutable int nFTHits         = 0;
  mutable int nbHits_in_UT    = 0;
  mutable int nbHits_in_SciFi = 0;

  mutable std::vector<float>        UT_cos;
  mutable std::vector<float>        UT_cosT;
  mutable std::vector<float>        UT_dxDy;
  mutable std::vector<unsigned int> UT_lhcbID;
  mutable std::vector<int>          UT_planeCode;
  mutable std::vector<float>        UT_sinT;
  mutable std::vector<int>          UT_size;
  mutable std::vector<float>        UT_tanT;
  mutable std::vector<float>        UT_weight;
  mutable std::vector<float>        UT_xAtYEq0;
  mutable std::vector<float>        UT_xAtYMid;
  mutable std::vector<float>        UT_xMax;
  mutable std::vector<float>        UT_xMin;
  mutable std::vector<float>        UT_xT;
  mutable std::vector<float>        UT_yBegin;
  mutable std::vector<float>        UT_yEnd;
  mutable std::vector<float>        UT_yMax;
  mutable std::vector<float>        UT_yMid;
  mutable std::vector<float>        UT_yMin;
  mutable std::vector<float>        UT_zAtYEq0;

  mutable int nUTHits = 0;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( PrTrackerDumper2 )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrTrackerDumper2::PrTrackerDumper2( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator,
                {KeyValue{"MCParticlesLocation", LHCb::MCParticleLocation::Default},
                 KeyValue{"VPLightClusterLocation", LHCb::VPClusterLocation::Light},
                 KeyValue{"FTHitsLocation", PrFTInfo::FTHitsLocation}, KeyValue{"UTHitsLocation", UTInfo::HitLocation},
                 KeyValue{"ODINLocation", LHCb::ODINLocation::Default},
                 KeyValue{"LinkerLocation", Links::location( "Pr/LHCbID" )},
                 KeyValue{"MCTrackInfo", LHCb::MCPropertyLocation::TrackInfo},
                 KeyValue{"TrackLinks", LHCb::LinksByKey::linkerName( "Rec/Track/Keyed/Velo" )}} ) {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrTrackerDumper2::initialize() {
  return Consumer::initialize().andThen( [&] {
    std::ostringstream oss;
    oss << "DumperFTUTHits_runNb_"
        //<< std::fixed << std::setfill('0') << std::setw(6) << std::to_string( odin.runNumber())
        //  << "_evtNb_" << std::setfill('0') << std::setw(6) << std::to_string(odin.eventNumber())
        << ".root";
    TString filename = oss.str();
    // TFile *
    file = std::make_unique<TFile>( filename.Data(), "RECREATE" );

    // TTree *
    tree = std::make_unique<TTree>( "Hits_detectors", "Hits_detectors" );

    tree->Branch( "eventID", &eventID );

    tree->Branch( "fullInfo", &fullInfo );
    tree->Branch( "hasScifi", &hasSciFi );
    tree->Branch( "hasUT", &hasUT );
    tree->Branch( "hasVelo", &hasVelo );
    tree->Branch( "isDown", &isDown );
    tree->Branch( "isDown_noVelo", &isDown_noVelo );
    tree->Branch( "isLong", &isLong );
    tree->Branch( "isLong_andUT", &isLong_andUT );
    tree->Branch( "p", &p );
    tree->Branch( "pt", &pt );
    tree->Branch( "px", &px );
    tree->Branch( "py", &py );
    tree->Branch( "pz", &pz );
    tree->Branch( "eta", &eta );
    tree->Branch( "phi", &phi );
    tree->Branch( "ovtx_x", &ovtx_x );
    tree->Branch( "ovtx_y", &ovtx_y );
    tree->Branch( "ovtx_z", &ovtx_z );
    tree->Branch( "pid", &pid );
    tree->Branch( "DecayOriginMother_pid", &DecayOriginMother_pid );
    tree->Branch( "key", &key );
    tree->Branch( "fromBeautyDecay", &fromBeautyDecay );
    tree->Branch( "fromCharmDecay", &fromCharmDecay );
    tree->Branch( "fromStrangeDecay", &fromStrangeDecay );
    tree->Branch( "eta_track", &eta_track );
    tree->Branch( "phi_track", &phi_track );
    tree->Branch( "chi2_track", &chi2_track );

    tree->Branch( "nVeloHits_track", &nVeloHits_track );
    tree->Branch( "Velo_x_track", &Velo_x_track );
    tree->Branch( "Velo_y_track", &Velo_y_track );
    tree->Branch( "Velo_z_track", &Velo_z_track );
    tree->Branch( "Velo_Module_track", &Velo_Module_track );
    tree->Branch( "Velo_Sensor_track", &Velo_Sensor_track );
    tree->Branch( "Velo_Station_track", &Velo_Station_track );
    tree->Branch( "Velo_lhcbID_track", &Velo_lhcbID_track );
    tree->Branch( "Velo_index_track", &Velo_index_track );

    tree->Branch( "nVeloHits", &nVeloHits );
    tree->Branch( "Velo_x", &Velo_x );
    tree->Branch( "Velo_y", &Velo_y );
    tree->Branch( "Velo_z", &Velo_z );
    tree->Branch( "Velo_Module", &Velo_Module );
    tree->Branch( "Velo_Sensor", &Velo_Sensor );
    tree->Branch( "Velo_Station", &Velo_Station );
    tree->Branch( "Velo_lhcbID", &Velo_lhcbID );
    tree->Branch( "Velo_index", &Velo_index );

    tree->Branch( "nbHits_in_UT", &nbHits_in_UT );
    tree->Branch( "nbHits_in_SciFi", &nbHits_in_SciFi );

    // SciFi

    tree->Branch( "nFTHits", &nFTHits );
    tree->Branch( "FT_x", &FT_hitx );
    tree->Branch( "FT_z", &FT_hitz );
    tree->Branch( "FT_w", &FT_hitw );
    tree->Branch( "FT_dxdy", &FT_hitDXDY );
    tree->Branch( "FT_YMin", &FT_hitYMin );
    tree->Branch( "FT_YMax", &FT_hitYMax );
    tree->Branch( "FT_hitPlaneCode", &FT_hitPlaneCode );
    tree->Branch( "FT_hitzone", &FT_hitzone );
    tree->Branch( "FT_lhcbID", &FT_lhcbID );

    // UT info

    tree->Branch( "nUTHits", &nUTHits );
    tree->Branch( "UT_cos", &UT_cos );
    tree->Branch( "UT_cosT", &UT_cosT );
    tree->Branch( "UT_dxDy", &UT_dxDy );
    tree->Branch( "UT_lhcbID", &UT_lhcbID );
    tree->Branch( "UT_planeCode", &UT_planeCode );
    tree->Branch( "UT_sinT", &UT_sinT );
    tree->Branch( "UT_size", &UT_size );
    tree->Branch( "UT_tanT", &UT_tanT );
    tree->Branch( "UT_weight", &UT_weight );
    tree->Branch( "UT_xAtYEq0", &UT_xAtYEq0 );
    tree->Branch( "UT_xAtYMid", &UT_xAtYMid );
    tree->Branch( "UT_xMax", &UT_xMax );
    tree->Branch( "UT_xMin", &UT_xMin );
    tree->Branch( "UT_xT", &UT_xT );
    tree->Branch( "UT_yBegin", &UT_yBegin );
    tree->Branch( "UT_yEnd", &UT_yEnd );
    tree->Branch( "UT_yMax", &UT_yMax );
    tree->Branch( "UT_yMid", &UT_yMid );
    tree->Branch( "UT_yMin", &UT_yMin );
    tree->Branch( "UT_zAtYEq0", &UT_zAtYEq0 );
    return StatusCode::SUCCESS;
  } );
}

//=============================================================================
// Finalization
//=============================================================================
StatusCode PrTrackerDumper2::finalize() {
  return Consumer::finalize().andThen( [&] {
    file->Write();
    file->Close();
  } );
}

//=============================================================================
// operator()
//=============================================================================
void PrTrackerDumper2::operator()( const LHCb::MCParticles& MCParticles, const LHCb::VPLightClusters& VPClusters,
                                   const PrFTHitHandler<PrHit>& prFTHitHandler, const UT::HitHandler& prUTHitHandler,
                                   const LHCb::ODIN& odin, const LHCb::LinksByKey& links,
                                   const LHCb::MCProperty& mcProperty, const LHCb::LinksByKey& seed_links ) const {

  verbose() << "Starting to dump..." << endmsg;
  auto lock = std::scoped_lock{m_mutex};

  // Look for associated MC  particle to the hit
  InputLinks<ContainedObject, LHCb::MCParticle> HitMCParticleLinks( links );

  verbose() << "SciFi" << endmsg;
  // SciFi
  std::map<const LHCb::MCParticle*, std::vector<PrHit>> FTHits_on_MCParticles;
  std::vector<PrHit>                                    non_Assoc_FTHits;

  for ( unsigned int zone = 0; LHCb::Detector::FT::nbZones() > zone; ++zone ) {
    for ( const auto& hit : prFTHitHandler.hits( zone ) ) {
      // get the LHCbID from the PrHit
      LHCb::LHCbID lhcbid = hit.id();
      // Get the linking to the MCParticle given the LHCbID
      auto mcparticlesrelations = HitMCParticleLinks.from( lhcbid.lhcbID() );

      if ( mcparticlesrelations.empty() ) { non_Assoc_FTHits.push_back( hit ); }
      for ( const auto& mcp : mcparticlesrelations ) {
        // MCP is MCParticle*
        auto MCP = mcp.to();
        //---> weightassociation = mcp.weight();
        FTHits_on_MCParticles[MCP].push_back( hit );
      }
    }
  }

  verbose() << "UT" << endmsg;
  // UT detector. loop over all hits in detector, extract for each MCParticle the vector<Hit> , then
  // See Pr/PrKernel/UTHit definitions to know the info to store
  std::map<const LHCb::MCParticle*, std::vector<UT::Hit>> UTHits_on_MCParticles;
  std::vector<UT::Hit>                                    non_Assoc_UTHits;
  for ( int iStation = 1; iStation < 3; ++iStation ) {
    for ( int iLayer = 1; iLayer < 3; ++iLayer ) {
      for ( int iRegion = 1; iRegion < 4; ++iRegion ) {
        for ( int iSector = 1; iSector < 99; ++iSector ) {
          for ( auto& hit : prUTHitHandler.hits( iStation, iLayer, iRegion, iSector ) ) {
            LHCb::LHCbID lhcbid = hit.lhcbID();

            auto mcparticlesrelations = HitMCParticleLinks.from( lhcbid.lhcbID() );
            if ( mcparticlesrelations.empty() ) {
              non_Assoc_UTHits.push_back( hit );
            } else {

              for ( const auto& mcp : mcparticlesrelations ) {
                auto MCP = mcp.to();
                //---> weightassociation = mcp.weight();
                UTHits_on_MCParticles[MCP].push_back( hit );
              }
            }
          }
        }
      }
    }
  }

  verbose() << "VP" << endmsg;
  // VP Detector
  std::map<const LHCb::MCParticle*, std::vector<std::pair<LHCb::VPLightCluster, unsigned int>>> VPHits_on_MCParticles;
  std::vector<std::pair<LHCb::VPLightCluster, unsigned int>>                                    non_Assoc_VPHits;
  std::cout << "Nb Velo Clusters in TES = " << VPClusters.size() << std::endl;

  // adding index within the module (bank) assuming that cluster are put on TES in order (true for
  // VPRetinaClusterCreator)
  unsigned int idx_cluster = 0, prev_module = 999;

  for ( const auto& vpclus : VPClusters ) {
    unsigned int the_module = vpclus.channelID().module();
    if ( the_module != prev_module ) {
      prev_module = the_module;
      idx_cluster = 0;
    }

    LHCb::LHCbID lhcbid               = LHCb::LHCbID( vpclus.channelID() );
    auto         mcparticlesrelations = HitMCParticleLinks.from( lhcbid.lhcbID() );
    if ( mcparticlesrelations.empty() ) {
      non_Assoc_VPHits.emplace_back( vpclus, idx_cluster );
    } else {
      for ( const auto& mcp : mcparticlesrelations ) {
        VPHits_on_MCParticles[mcp.to()].emplace_back( vpclus, idx_cluster );
      }
    }
    ++idx_cluster;
  }

  //---- We use trackInfo for a given MCParticle to know if the particle is reconstructible or not
  verbose() << "Track" << endmsg;

  eventID = odin.eventNumber();

  const auto trackInfo = MCTrackInfo{mcProperty};

  nbHits_in_UT    = computeNumberUTHits( prUTHitHandler );
  nbHits_in_SciFi = prFTHitHandler.hits().size();

  verbose() << "Loop on particles" << endmsg;

  // std::string track_location = "Rec/Track/Keyed/Velo";
  for ( const auto* mcparticle : MCParticles ) {
    //---- We can speed up things if we filter only tracks which are either reconstructible in Velo or UT or SciFi,
    //---- Here is very inefficient, we go through ALL MCParticles in the chain, even the non-final states one
    /*
      if( ! ( trackInfo.hasVelo( mcparticle) || trackInfo.hasT( mcparticle) || trackInfo.hasSciFiT( mcparticle) ) ){
        continue;
      };
    */

    const LHCb::Track* mcSeedPart = LinkedFrom<LHCb::Track>{&seed_links}.range( mcparticle ).try_front();

    // Velo
    nVeloHits = 0;
    Velo_x.clear();
    Velo_y.clear();
    Velo_z.clear();
    Velo_Module.clear();
    Velo_Sensor.clear();
    Velo_Station.clear();
    Velo_lhcbID.clear();
    Velo_index.clear();

    if ( VPHits_on_MCParticles.find( mcparticle ) != VPHits_on_MCParticles.end() ) {
      nVeloHits = VPHits_on_MCParticles[mcparticle].size();
      for ( auto& vphit_pair : VPHits_on_MCParticles[mcparticle] ) {
        auto vphit = vphit_pair.first;
        Velo_x.push_back( vphit.x() );
        Velo_y.push_back( vphit.y() );
        Velo_z.push_back( vphit.z() );
        Velo_Module.push_back( to_unsigned( vphit.channelID().sensor() ) / 4 );
        Velo_Sensor.push_back( to_unsigned( vphit.channelID().sensor() ) );
        Velo_Station.push_back( vphit.channelID().station() );
        Velo_lhcbID.push_back( LHCb::LHCbID( vphit.channelID() ).lhcbID() );
        Velo_index.push_back( vphit_pair.second );
      }
    }

    nFTHits = 0;
    // SciFi
    FT_hitz.clear();
    FT_hitx.clear();
    FT_hitw.clear();
    FT_hitPlaneCode.clear();
    FT_hitzone.clear();
    FT_hitDXDY.clear();
    FT_hitYMin.clear();
    FT_hitYMax.clear();
    FT_lhcbID.clear();

    if ( FTHits_on_MCParticles.find( mcparticle ) != FTHits_on_MCParticles.end() ) {
      nFTHits = (int)FTHits_on_MCParticles[mcparticle].size();
      for ( auto& fthit : FTHits_on_MCParticles[mcparticle] ) {
        FT_hitz.push_back( fthit.z() );
        FT_hitx.push_back( fthit.x() );
        FT_hitw.push_back( fthit.w() );
        FT_hitPlaneCode.push_back( fthit.planeCode() );
        FT_hitzone.push_back( fthit.zone() );
        FT_hitDXDY.push_back( fthit.dxDy() );
        FT_hitYMin.push_back( fthit.yMin() );
        FT_hitYMax.push_back( fthit.yMax() );
        FT_lhcbID.push_back( fthit.id().lhcbID() );
      }
    }

    nUTHits = 0;
    UT_cos.clear();
    UT_cosT.clear();
    UT_dxDy.clear();
    UT_lhcbID.clear();
    UT_planeCode.clear();
    UT_sinT.clear();
    UT_size.clear();
    UT_tanT.clear();
    UT_weight.clear();
    UT_xAtYEq0.clear();
    UT_xAtYMid.clear();
    UT_xMax.clear();
    UT_xMin.clear();
    UT_xT.clear();
    UT_yBegin.clear();
    UT_yEnd.clear();
    UT_yMax.clear();
    UT_yMid.clear();
    UT_yMin.clear();
    UT_zAtYEq0.clear();

    if ( UTHits_on_MCParticles.find( mcparticle ) != UTHits_on_MCParticles.end() ) {
      nUTHits = (int)UTHits_on_MCParticles[mcparticle].size();
      for ( const auto& uthit : UTHits_on_MCParticles[mcparticle] ) {
        UT_cos.push_back( uthit.cos() );
        UT_cosT.push_back( uthit.cosT() );
        UT_dxDy.push_back( uthit.dxDy() );
        UT_lhcbID.push_back( uthit.chanID().channelID() );
        UT_planeCode.push_back( uthit.planeCode() );
        UT_sinT.push_back( uthit.sinT() );
        UT_size.push_back( uthit.size() );
        UT_tanT.push_back( uthit.tanT() );
        UT_weight.push_back( uthit.weight() );
        UT_xAtYEq0.push_back( uthit.xAtYEq0() );
        UT_xAtYMid.push_back( uthit.xAtYMid() );
        UT_xMax.push_back( uthit.xMax() );
        UT_xMin.push_back( uthit.xMin() );
        UT_xT.push_back( uthit.xT() );
        UT_yBegin.push_back( uthit.yBegin() );
        UT_yEnd.push_back( uthit.yEnd() );
        UT_yMax.push_back( uthit.yMax() );
        UT_yMid.push_back( uthit.yMid() );
        UT_yMin.push_back( uthit.yMin() );
        UT_zAtYEq0.push_back( uthit.zAtYEq0() );
      }
    }

    verbose() << "Skipped hits" << endmsg;
    // skip the MC particles without any hits in the tracking system
    if ( nFTHits == 0 && nVeloHits == 0 && nUTHits == 0 ) continue;
    // probably if fullInfo ==0 skip does the same , to check
    fullInfo      = trackInfo.fullInfo( mcparticle );
    hasSciFi      = trackInfo.hasT( mcparticle );
    hasUT         = trackInfo.hasUT( mcparticle );
    hasVelo       = trackInfo.hasVelo( mcparticle );
    isDown        = hasSciFi && hasUT;
    isDown_noVelo = hasSciFi && hasUT && !( hasVelo );
    isLong        = hasSciFi && hasVelo;
    isLong_andUT  = hasSciFi && hasVelo && hasUT;
    p             = mcparticle->p();
    px            = mcparticle->momentum().Px();
    py            = mcparticle->momentum().Py();
    pz            = mcparticle->momentum().Pz();
    pt            = mcparticle->pt();
    eta           = mcparticle->momentum().Eta();
    phi           = mcparticle->momentum().phi();
    pid = mcparticle->particleID().pid(); // offline you want to match the PID eventually to the e+, e- or whatever
    fromBeautyDecay       = false;
    fromCharmDecay        = false;
    fromStrangeDecay      = false;
    DecayOriginMother_pid = -999999;
    ovtx_x                = std::numeric_limits<double>::min();
    ovtx_y                = std::numeric_limits<double>::min();
    ovtx_z                = std::numeric_limits<double>::min();
    key                   = mcparticle->key();

    // Added code to dump also informations of reconstructed track associated to MC particle

    if ( mcSeedPart ) {

      Velo_x_track.clear();
      Velo_y_track.clear();
      Velo_z_track.clear();
      Velo_Module_track.clear();
      Velo_Sensor_track.clear();
      Velo_Station_track.clear();
      Velo_lhcbID_track.clear();
      Velo_index_track.clear();

      eta_track       = mcSeedPart->pseudoRapidity();
      phi_track       = mcSeedPart->phi();
      chi2_track      = mcSeedPart->chi2();
      auto ids        = mcSeedPart->lhcbIDs();
      nVeloHits_track = ids.size();
      for ( auto& id : ids ) {
        if ( id.isVP() ) {
          auto vp_ID = id.vpID();

          bool foundID = false;
          for ( auto& vphit : VPClusters )
            if ( vphit.channelID() == vp_ID.channelID() ) {
              foundID = true;
              Velo_x_track.push_back( vphit.x() );
              Velo_y_track.push_back( vphit.y() );
              Velo_z_track.push_back( vphit.z() );
              Velo_Module_track.push_back( vp_ID.module() );
              Velo_Sensor_track.push_back( to_unsigned( vp_ID.sensor() ) );
              Velo_Station_track.push_back( vp_ID.station() );

              Velo_lhcbID_track.push_back( vp_ID.channelID() );
              Velo_index_track.push_back( 0 );
              break;
            }
          if ( !foundID ) error() << "Hit not found: " << vp_ID.channelID() << endmsg;
        }
      }
    } else {
      eta_track       = 0;
      phi_track       = 0;
      chi2_track      = 0;
      nVeloHits_track = 0;
    }

    // navigate decay back to mother origin
    if ( nullptr != mcparticle->originVertex() ) {
      // store the mcparticle origin vertex information , and navigate back to mother of the particle!

      ovtx_x                         = mcparticle->originVertex()->position().x();
      ovtx_y                         = mcparticle->originVertex()->position().y();
      ovtx_z                         = mcparticle->originVertex()->position().z();
      const LHCb::MCParticle* mother = mcparticle->originVertex()->mother();
      if ( nullptr != mother ) {
        if ( nullptr != mother->originVertex() ) {
          double rOrigin = mother->originVertex()->position().rho();
          if ( fabs( rOrigin ) < 5. ) { // radial origin position of the mother within 5 mm from beam pipe
            int pid = abs( mother->particleID().pid() );
            if ( 130 == pid ||  // K0L
                 310 == pid ||  // K0S
                 3122 == pid || // Lambda
                 3222 == pid || // Sigma+
                 3212 == pid || // Sigma0
                 3112 == pid || // Sigma-
                 3322 == pid || // Xsi0
                 3312 == pid || // Xsi-
                 3334 == pid    // Omega-
            ) {
              fromStrangeDecay = true;
            }
          }
        }
      }
      while ( mother ) {
        if ( mother->particleID().hasBottom() &&
             ( mother->particleID().isMeson() || mother->particleID().isBaryon() ) ) {
          DecayOriginMother_pid = mother->particleID().pid();
          fromBeautyDecay       = true;
        }
        if ( mother->particleID().hasCharm() &&
             ( mother->particleID().isMeson() || mother->particleID().isBaryon() ) ) {
          fromCharmDecay        = true;
          DecayOriginMother_pid = mother->particleID().pid();
        }
        mother = mother->originVertex()->mother();
      }
    }
    tree->Fill();
  }

  verbose() << "Last track (fake)" << endmsg;

  // We filled the tree with hits having a MCParticle linked to [ no filter done if the
  // particle is reconstructible or not in Velo/UT/SciFi]
  //---- Offline, to grab the hits on reconstructible tracks plot the ones having the
  // flag hasUT or hasT or hasVelo or combine the flags to your preference

  // Empty the vectors of info before fillong the remaining non-associated hits [offline
  // you want to check uniqueness of lhcbID info to have the actual hits to use for tracking, since
  // 1 hit can be associated to more MCParticles]
  FT_hitz.clear();
  FT_hitx.clear();
  FT_hitw.clear();
  FT_hitPlaneCode.clear();
  FT_hitzone.clear();
  FT_hitDXDY.clear();
  FT_hitYMin.clear();
  FT_hitYMax.clear();
  FT_lhcbID.clear();
  // store all remaining hits in a dummy tuple, non associated ones in FT for the event!
  nFTHits = non_Assoc_FTHits.size();
  for ( const auto& fthit : non_Assoc_FTHits ) {
    FT_hitz.push_back( fthit.z() );
    FT_hitx.push_back( fthit.x() );
    FT_hitw.push_back( fthit.w() );
    FT_hitPlaneCode.push_back( fthit.planeCode() );
    FT_hitzone.push_back( fthit.zone() );
    FT_hitDXDY.push_back( fthit.dxDy() );
    FT_hitYMin.push_back( fthit.yMin() );
    FT_hitYMax.push_back( fthit.yMax() );
    FT_lhcbID.push_back( fthit.id().lhcbID() );
  }

  nUTHits = non_Assoc_UTHits.size();
  // Empty the vectors of info before fillong the remaining non-associated hits [offline
  // you want to check uniqueness of lhcbID info to have the actual hits to use for tracking, since
  // 1 hit can be associated to more MCParticles]
  UT_cos.clear();
  UT_cosT.clear();
  UT_dxDy.clear();
  UT_lhcbID.clear();
  UT_planeCode.clear();
  UT_sinT.clear();
  UT_size.clear();
  UT_tanT.clear();
  UT_weight.clear();
  UT_xAtYEq0.clear();
  UT_xAtYMid.clear();
  UT_xMax.clear();
  UT_xMin.clear();
  UT_xT.clear();
  UT_yBegin.clear();
  UT_yEnd.clear();
  UT_yMax.clear();
  UT_yMid.clear();
  UT_yMin.clear();
  UT_zAtYEq0.clear();
  for ( const auto& uthit : non_Assoc_UTHits ) {
    UT_cos.push_back( uthit.cos() );
    UT_cosT.push_back( uthit.cosT() );
    UT_dxDy.push_back( uthit.dxDy() );
    UT_lhcbID.push_back( uthit.lhcbID().lhcbID() );
    UT_planeCode.push_back( uthit.planeCode() );
    UT_sinT.push_back( uthit.sinT() );
    UT_size.push_back( uthit.size() );
    UT_tanT.push_back( uthit.tanT() );
    UT_weight.push_back( uthit.weight() );
    UT_xAtYEq0.push_back( uthit.xAtYEq0() );
    UT_xAtYMid.push_back( uthit.xAtYMid() );
    UT_xMax.push_back( uthit.xMax() );
    UT_xMin.push_back( uthit.xMin() );
    UT_xT.push_back( uthit.xT() );
    UT_yBegin.push_back( uthit.yBegin() );
    UT_yEnd.push_back( uthit.yEnd() );
    UT_yMax.push_back( uthit.yMax() );
    UT_yMid.push_back( uthit.yMid() );
    UT_yMin.push_back( uthit.yMin() );
    UT_zAtYEq0.push_back( uthit.zAtYEq0() );
  }

  // velo part
  Velo_x.clear();
  Velo_y.clear();
  Velo_z.clear();
  Velo_Module.clear();
  Velo_Sensor.clear();
  Velo_Station.clear();
  Velo_lhcbID.clear();
  Velo_index.clear();
  nVeloHits = non_Assoc_VPHits.size();
  for ( const auto& vphit_pair : non_Assoc_VPHits ) {
    auto vphit = vphit_pair.first;
    Velo_x.push_back( vphit.x() );
    Velo_y.push_back( vphit.y() );
    Velo_z.push_back( vphit.z() );
    Velo_Module.push_back( to_unsigned( vphit.channelID().sensor() ) / 4 );
    Velo_Sensor.push_back( to_unsigned( vphit.channelID().sensor() ) );
    Velo_Station.push_back( vphit.channelID().station() );
    Velo_lhcbID.push_back( LHCb::LHCbID( vphit.channelID() ).lhcbID() );
    Velo_index.push_back( vphit_pair.second );
  }

  fullInfo              = false;
  hasSciFi              = false;
  hasUT                 = false;
  hasVelo               = false;
  isDown                = false;
  isDown_noVelo         = false;
  isLong                = false;
  isLong_andUT          = false;
  p                     = -9999999999999.;
  pt                    = -9999999999999.;
  eta                   = -9999999999999.;
  pid                   = -999999;
  fromBeautyDecay       = false;
  fromCharmDecay        = false;
  fromStrangeDecay      = false;
  DecayOriginMother_pid = -999999;
  ovtx_x                = -9999999999999.;
  ovtx_y                = -9999999999999.;
  ovtx_z                = -9999999999999.;
  key                   = -999999;

  tree->Fill();
}
