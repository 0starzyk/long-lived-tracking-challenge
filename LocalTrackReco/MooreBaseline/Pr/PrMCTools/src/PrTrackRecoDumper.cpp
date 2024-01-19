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
#include "Linker/LinkedTo.h"
#include "PrKernel/PrFTHitHandler.h"
#include "PrKernel/PrHit.h"
#include "PrKernel/UTHit.h"
#include "PrKernel/UTHitHandler.h"
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <utility>
//-----------------------------------------------------------------------------
// Implementation file for class : PrTrackRecoDumper
//
// 07-2018 Riccardo Cenci
// 21-11-2018 Giulia Tuci

// Based on PrTrackerDumper. Here, instead of MC information, we dump information of reconstructed tracks.
// If the variable isMatched==1, then the reconstructed track is correctly associated to a MC particle.

// Reconstructed tracks information:
/*
  p             : track momentum
  pt            : track trsansverse momentum
  eta           : track pseudorapidity
  ovtx_x        : track origin position X
  ovtx_y        : track origin position Y
  ovtx_z        : track origin position Z
  fromBeautyDecay : the track belongs to a decay chain with a b-quark hadron
  fromCharmDecay  : the track belongs to a decay chain with a c-quark hadron
  fromStrangeDecay : the track belongs to a decay chain with a s-quark hadron

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

/** @class PrTraackRecoDumper PrTrackRecoDumper.h
 *  Dumping of reconstructed tracks and truth of MC particle associated to these tracks
 *
 */
/*

 */

// typedef std::vector<LHCb::Track> Tracks;
typedef LHCb::Tracks Tracks;

class PrTrackRecoDumper
    : public LHCb::Algorithm::Consumer<void( const Tracks&, const LHCb::VPLightClusters&, const LHCb::ODIN&,
                                             const PrFTHitHandler<PrHit>&, const UT::HitHandler&,
                                             const LHCb::LinksByKey&, const LHCb::MCParticles& )> {
public:
  /// Standard constructor
  PrTrackRecoDumper( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override;
  StatusCode finalize() override;

  void operator()( const Tracks& recTracks, const LHCb::VPLightClusters& VPClusters, const LHCb::ODIN& odin,
                   const PrFTHitHandler<PrHit>& ftHits, const UT::HitHandler& utHits, const LHCb::LinksByKey& links,
                   const LHCb::MCParticles& mcParts ) const override;

private:
  mutable std::mutex m_mutex;

  TFile* file;
  TTree* tree;

  mutable int eventID;

  mutable double p;
  mutable double px, py, pz;
  mutable double pt;
  mutable double eta, phi;
  mutable bool   isMatched;
  // vertex origin of the particle
  mutable double ovtx_x;
  mutable double ovtx_y;
  mutable double ovtx_z;
  mutable int    pid;
  mutable bool   fromBeautyDecay;
  mutable bool   fromCharmDecay;
  mutable bool   fromStrangeDecay;
  mutable int    key;
  mutable double ghostProb;
  mutable double chi2;
  mutable int    ndof;

  mutable int                       nVeloHits;
  mutable std::vector<float>        Velo_x;
  mutable std::vector<float>        Velo_y;
  mutable std::vector<float>        Velo_z;
  mutable std::vector<int>          Velo_Module;
  mutable std::vector<int>          Velo_Sensor;
  mutable std::vector<int>          Velo_Station;
  mutable std::vector<unsigned int> Velo_lhcbID;

  mutable std::vector<float>        FT_hitx;
  mutable std::vector<float>        FT_hitz;
  mutable std::vector<float>        FT_hitw;
  mutable std::vector<float>        FT_hitDXDY;
  mutable std::vector<float>        FT_hitYMin;
  mutable std::vector<float>        FT_hitYMax;
  mutable std::vector<int>          FT_hitPlaneCode;
  mutable std::vector<int>          FT_hitzone;
  mutable std::vector<unsigned int> FT_lhcbID;
  mutable int                       nFTHits;

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
  mutable int                       nUTHits;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( PrTrackRecoDumper )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================

PrTrackRecoDumper::PrTrackRecoDumper( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator,
                {KeyValue{"TrackLocation", "Rec/Track/Velo"},
                 KeyValue{"VPLightClusterLocation", LHCb::VPClusterLocation::Light},
                 KeyValue{"ODINLocation", LHCb::ODINLocation::Default},
                 KeyValue{"FTHitsLocation", PrFTInfo::FTHitsLocation}, KeyValue{"UTHitsLocation", UTInfo::HitLocation},
                 KeyValue{"LinksLocation", ""}, KeyValue{"MCParticlesLocation", ""}} ) {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrTrackRecoDumper::initialize() {

  StatusCode sc = Consumer::initialize();
  if ( sc.isFailure() ) return sc;

  std::ostringstream oss;
  oss << "Dumper_recTracks"
      << ".root";
  TString filename = oss.str();

  // Reserve space for the vectors
  const unsigned int maxVPhits = 50;
  const unsigned int maxFThits = 1000;
  const unsigned int maxUThits = 1000;

  Velo_x.reserve( maxVPhits );
  Velo_y.reserve( maxVPhits );
  Velo_z.reserve( maxVPhits );
  Velo_Module.reserve( maxVPhits );
  Velo_Sensor.reserve( maxVPhits );
  Velo_Station.reserve( maxVPhits );
  Velo_lhcbID.reserve( maxVPhits );

  FT_hitx.reserve( maxFThits );
  FT_hitz.reserve( maxFThits );
  FT_hitw.reserve( maxFThits );
  FT_hitDXDY.reserve( maxFThits );
  FT_hitYMin.reserve( maxFThits );
  FT_hitYMax.reserve( maxFThits );
  FT_hitPlaneCode.reserve( maxFThits );
  FT_hitzone.reserve( maxFThits );
  FT_lhcbID.reserve( maxFThits );

  UT_cos.reserve( maxUThits );
  UT_cosT.reserve( maxUThits );
  UT_dxDy.reserve( maxUThits );
  UT_lhcbID.reserve( maxUThits );
  UT_planeCode.reserve( maxUThits );
  UT_sinT.reserve( maxUThits );
  UT_size.reserve( maxUThits );
  UT_tanT.reserve( maxUThits );
  UT_weight.reserve( maxUThits );
  UT_xAtYEq0.reserve( maxUThits );
  UT_xAtYMid.reserve( maxUThits );
  UT_xMax.reserve( maxUThits );
  UT_xMin.reserve( maxUThits );
  UT_xT.reserve( maxUThits );
  UT_yBegin.reserve( maxUThits );
  UT_yEnd.reserve( maxUThits );
  UT_yMax.reserve( maxUThits );
  UT_yMid.reserve( maxUThits );
  UT_yMin.reserve( maxUThits );
  UT_zAtYEq0.reserve( maxUThits );

  // TFile *
  file = new TFile( filename.Data(), "RECREATE" );

  // TTree *
  tree = new TTree( "Hits_detectors", "Hits_detectors" );

  tree->Branch( "eventID", &eventID );
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
  tree->Branch( "key", &key );
  tree->Branch( "ghostProb", &ghostProb );
  tree->Branch( "chi2", &chi2 );
  tree->Branch( "ndof", &ndof );
  tree->Branch( "fromBeautyDecay", &fromBeautyDecay );
  tree->Branch( "fromCharmDecay", &fromCharmDecay );
  tree->Branch( "fromStrangeDecay", &fromStrangeDecay );
  tree->Branch( "isMatched", &isMatched );

  tree->Branch( "nVeloHits", &nVeloHits );
  tree->Branch( "Velo_x", &Velo_x );
  tree->Branch( "Velo_y", &Velo_y );
  tree->Branch( "Velo_z", &Velo_z );
  tree->Branch( "Velo_Module", &Velo_Module );
  tree->Branch( "Velo_Sensor", &Velo_Sensor );
  tree->Branch( "Velo_Station", &Velo_Station );
  tree->Branch( "Velo_lhcbID", &Velo_lhcbID );

  nFTHits = int( 0 );

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

  nUTHits = int( 0 );

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

  return sc;
}

//=============================================================================
// Finalization
//=============================================================================
StatusCode PrTrackRecoDumper::finalize() {
  return Consumer::finalize().andThen( [&] {
    file->Write();
    file->Close();
  } );
}

//=============================================================================
// operator()
//=============================================================================
void PrTrackRecoDumper::operator()( const Tracks& recTracks, const LHCb::VPLightClusters& VPClusters,
                                    const LHCb::ODIN& odin, const PrFTHitHandler<PrHit>& prFTHitHandler,
                                    const UT::HitHandler& prUTHitHandler, const LHCb::LinksByKey& links,
                                    const LHCb::MCParticles& mcParts ) const {
  std::scoped_lock lock( m_mutex );
  verbose() << "Starting to dump..." << endmsg;

  verbose() << "Track" << endmsg;

  eventID = odin.eventNumber();

  verbose() << "Loop on tracks" << endmsg;

  for ( const auto& track : recTracks ) {
    // Here we retireve the link between MC particles and reconstructed tracks
    const LHCb::MCParticle* mcSeedPart{nullptr};
    double                  maxWeight( 0 );
    links.applyToLinks( track->key(),
                        [&maxWeight, &mcSeedPart, &mcParts]( unsigned int, unsigned int mcPartKey, float weight ) {
                          if ( weight > maxWeight ) {
                            maxWeight  = weight;
                            mcSeedPart = static_cast<const LHCb::MCParticle*>( mcParts.containedObject( mcPartKey ) );
                          }
                        } );
    isMatched = !( mcSeedPart == nullptr );

    // Information of reconstructed track
    p   = track->p();
    px  = track->momentum().x();
    py  = track->momentum().y();
    pz  = track->momentum().z();
    pt  = track->pt();
    eta = track->pseudoRapidity();
    phi = track->phi();
    pid = 0; // track->particleID().pid(); //offline you want to match the PID eventually to the e+, e- or whatever
    fromBeautyDecay  = false;
    fromCharmDecay   = false;
    fromStrangeDecay = false;
    ovtx_x           = track->firstState().position().x();
    ovtx_y           = track->firstState().position().y();
    ovtx_z           = track->firstState().position().z();
    key              = track->key();
    ghostProb        = track->ghostProbability();
    chi2             = track->chi2();
    ndof             = track->nDoF();

    // Velo
    nVeloHits = 0;
    Velo_x.clear();
    Velo_y.clear();
    Velo_z.clear();
    Velo_Module.clear();
    Velo_Sensor.clear();
    Velo_Station.clear();
    Velo_lhcbID.clear();

    // SciFi
    nFTHits = 0;
    FT_hitz.clear();
    FT_hitx.clear();
    FT_hitw.clear();
    FT_hitPlaneCode.clear();
    FT_hitzone.clear();
    FT_hitDXDY.clear();
    FT_hitYMin.clear();
    FT_hitYMax.clear();
    FT_lhcbID.clear();

    // UT
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

    auto ids = track->lhcbIDs();
    for ( auto& id : ids ) {
      if ( id.isVP() ) {
        auto vp_ID = id.vpID();
        nVeloHits++;
        bool foundID = false;
        for ( auto& vphit : VPClusters )
          if ( vphit.channelID() == vp_ID.channelID() ) {
            foundID = true;
            Velo_x.push_back( vphit.x() );
            Velo_y.push_back( vphit.y() );
            Velo_z.push_back( vphit.z() );
            Velo_Module.push_back( vp_ID.module() );
            Velo_Sensor.push_back( to_unsigned( vp_ID.sensor() ) );
            Velo_Station.push_back( vp_ID.station() );

            Velo_lhcbID.push_back( vp_ID.channelID() );
            break;
          }

        if ( !foundID ) error() << "Hit not found: " << vp_ID.channelID() << endmsg;
      }
      if ( id.isFT() ) {
        auto ft_ID   = id.ftID();
        nFTHits      = nFTHits + 1;
        bool foundID = false;
        for ( unsigned int zone = 0; LHCb::Detector::FT::nbZones() > zone; ++zone ) {
          for ( const auto& fthit : prFTHitHandler.hits( zone ) ) {
            if ( fthit.id().channelID() == ft_ID.channelID() ) {
              foundID = true;
              FT_hitz.push_back( fthit.z() );
              FT_hitx.push_back( fthit.x() );
              FT_hitw.push_back( fthit.w() );
              FT_hitPlaneCode.push_back( fthit.planeCode() );
              FT_hitzone.push_back( fthit.zone() );
              FT_hitDXDY.push_back( fthit.dxDy() );
              FT_hitYMin.push_back( fthit.yMin() );
              FT_hitYMax.push_back( fthit.yMax() );
              FT_lhcbID.push_back( fthit.id().channelID() );
              break;
            }
          }
        }
        if ( !foundID ) error() << "Hit not found: " << ft_ID.channelID() << endmsg;
      }

      if ( id.isUT() ) {
        auto ut_ID   = id.utID();
        bool foundID = false;
        nUTHits      = nUTHits + 1;
        for ( int iStation = 1; iStation < 3; ++iStation ) {
          for ( int iLayer = 1; iLayer < 3; ++iLayer ) {
            for ( int iRegion = 1; iRegion < 4; ++iRegion ) {
              for ( int iSector = 1; iSector < 99; ++iSector ) {
                for ( auto& uthit : prUTHitHandler.hits( iStation, iLayer, iRegion, iSector ) ) {
                  if ( uthit.chanID().channelID() == ut_ID.channelID() ) {
                    foundID = true;
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
                    break;
                  }
                }
              }
            }
          }
        }

        if ( !foundID ) error() << "Hit not found: " << ut_ID.channelID() << endmsg;
      }
    }

    tree->Fill();
  }
}
