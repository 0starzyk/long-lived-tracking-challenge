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
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/LinksByKey.h"
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/RawEvent.h"
#include "Event/VPDigit.h"
#include "Event/VPFullCluster.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "LHCbAlgs/Consumer.h"
#include "VPDet/DeVP.h"

#include <math.h>

/** @class VPClusterEfficiency VPClusterEfficiency.h
 *
 * Checks the VPCluster efficiency for simulation
 */

class VPClusterEfficiency
    : public LHCb::Algorithm::Consumer<void( const LHCb::RawEvent&, const std::vector<LHCb::VPFullCluster>&,
                                             const LHCb::MCHits&, const LHCb::MCParticles&, const LHCb::LinksByKey&,
                                             const LHCb::MCProperty&, const DeVP& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DeVP>> {
public:
  /// Standard constructor
  VPClusterEfficiency( const std::string& name, ISvcLocator* pSvcLocator );

  /// Consumer operator: takes VP clusters & MCHits to make efficiency plots
  void operator()( const LHCb::RawEvent&, const std::vector<LHCb::VPFullCluster>&, const LHCb::MCHits&,
                   const LHCb::MCParticles&, const LHCb::LinksByKey&, const LHCb::MCProperty&,
                   const DeVP& ) const override;

private:
  mutable Gaudi::Accumulators::SigmaCounter<>    m_num_clusters{this, "# clusters per event"};
  mutable Gaudi::Accumulators::SigmaCounter<>    m_num_pix_clu{this, "# pixels per cluster"};
  mutable Gaudi::Accumulators::SigmaCounter<>    m_num_pix_hit{this, "# pixels per MCHit"};
  mutable Gaudi::Accumulators::BinomialCounter<> m_efficiency{this, "Cluster Efficiency"};
  mutable Gaudi::Accumulators::SigmaCounter<>    m_residual_x{this, "Residuals x [mm]"};
  mutable Gaudi::Accumulators::SigmaCounter<>    m_residual_y{this, "Residuals y [mm]"};
  mutable Gaudi::Accumulators::SigmaCounter<>    m_purity{this, "Purity"};
};

DECLARE_COMPONENT( VPClusterEfficiency )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
VPClusterEfficiency::VPClusterEfficiency( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer{name,
               pSvcLocator,
               {KeyValue{"RawEventLocation", LHCb::RawEventLocation::Default},
                KeyValue{"VPClusterLocation", LHCb::VPFullClusterLocation::Default},
                KeyValue{"MCHitLocation", LHCb::MCHitLocation::VP},
                KeyValue{"MCParticleLocation", LHCb::MCParticleLocation::Default},
                KeyValue{"VPDigit2MCHitLinksLocation", LHCb::VPDigitLocation::Default + "2MCHits"},
                KeyValue{"MCProperty", LHCb::MCPropertyLocation::TrackInfo},
                KeyValue{"DeVP", LHCb::Det::VP::det_path}}} {}

//=============================================================================
// Main execution
//=============================================================================
void VPClusterEfficiency::operator()( const LHCb::RawEvent& rawEvent, const std::vector<LHCb::VPFullCluster>& clusters,
                                      const LHCb::MCHits&     mcHits, const LHCb::MCParticles&,
                                      const LHCb::LinksByKey& links, const LHCb::MCProperty& mcprop,
                                      const DeVP& det ) const {
  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Found " << clusters.size() << " VP clusters" << endmsg;
    debug() << "Found " << mcHits.size() << " mcHits" << endmsg;
  }

  // convert Superpixels into digits
  std::vector<LHCb::Detector::VPChannelID> digits;
  digits.reserve( clusters.size() * 2 ); // guess as to max pixels to clusters ratio
  const auto& tBanks = rawEvent.banks( LHCb::RawBank::VP );
  if ( tBanks.empty() ) { GaudiException( "Missing VP banks", "No VP superpixel banks", StatusCode::FAILURE ); }
  // offsets from SP address of component pixels
  /*
   * row,y
   *
   *  ^  37
   *  |  26
   *  |  15
   *  |  04
   *  +---> col,x
   */
  constexpr uint32_t CHIP_COLUMNS = 256; // hard code number of columns per chip
  // Loop over VP RawBanks
  for ( const auto& bank : tBanks ) {
    const auto beforeSize = digits.size();
    const auto sensor     = LHCb::Detector::VPChannelID::SensorID( bank->sourceID() ); // sensor = source ID
    auto       data       = bank->range<uint32_t>();

    assert( data.size() == data.front() + 1 );
    for ( const uint32_t sp_word :
          data.subspan( 1 ) ) {     // note: the subspan(1) will skip the first entry, which is `nsp`
      uint8_t sp = sp_word & 0xFFU; // mask out pattern [bit 0-7]
      if ( 0 == sp ) continue;      // protect against empty super pixels.

      // bit 8-13   Super Pixel Row (0-63), bit 14-22  Super Pixel Column (0-383)
      const uint32_t sp_addr = ( sp_word & 0x007FFF00U ) >> 8; // get bits and left shift align
      const uint32_t sp_row  = sp_addr & 0x3FU;                // pick out row
      const uint32_t sp_col  = ( sp_addr >> 6 );               // left shift to align column val
      for ( uint32_t off = 0; off < 8; ++off ) {
        const uint8_t mask = 0x1U << off;
        if ( sp & mask ) {                               // test each of the 8 bits in turn
          const uint32_t cx = sp_col * 2 + ( off >> 2 ); // 0->3 no offset, 4->7 offset = 1
          const auto     cy =
              LHCb::Detector::VPChannelID::RowID{sp_row * 4 + ( off & 0b11U )}; // 0->3 add off, 4->7 add off-4
          const auto chip = LHCb::Detector::VPChannelID::ChipID{cx / CHIP_COLUMNS};
          const auto ccol = LHCb::Detector::VPChannelID::ColumnID{cx % CHIP_COLUMNS};
          digits.emplace_back( sensor, chip, ccol, cy );
        }
      }
    }
    if ( msgLevel( MSG::VERBOSE ) ) {
      verbose() << "Found " << digits.size() - beforeSize << " pixels in sensor " << to_unsigned( sensor ) << endmsg;
    }
  }

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Found " << digits.size() << " VP Digits from SuperPixel banks" << endmsg;
  }

  // Table linking a LHCb::Detector::VPChannelID* to std::vector<LHCb::MCHit>
  std::map<const unsigned int, std::vector<LHCb::MCHit const*>> MCHitForchannelId;
  links.applyToAllLinks( [&MCHitForchannelId, &mcHits]( unsigned int channelId, unsigned int mcHitKey, float ) {
    MCHitForchannelId[channelId].emplace_back( mcHits[mcHitKey] );
  } );

  // Table linking a LHCb::MCHit* to std::vector<LHCb::Detector::VPChannelID>> --orig
  std::map<const LHCb::MCHit*, std::vector<unsigned int>> channelIdForMCHit;
  links.applyToAllLinks( [&channelIdForMCHit, &mcHits]( unsigned int channelId, unsigned int mcHitKey, float ) {
    channelIdForMCHit[mcHits[mcHitKey]].emplace_back( channelId );
  } );

  // split MCHits into modules
  std::array<std::vector<LHCb::MCHit const*>, 52> hitsInModules;
  for ( unsigned int i = 0; i < 52; ++i ) hitsInModules[i].reserve( mcHits.size() / 52 ); // make some space
  for ( auto& mcH : mcHits ) {
    unsigned int module = mcH->sensDetID() / 4; // 4 sensors per module, numbered consecutively
    hitsInModules[module].push_back( mcH );
  }

  if ( msgLevel( MSG::VERBOSE ) ) {
    for ( unsigned int i = 0; i < 52; ++i ) {
      verbose() << "Found " << hitsInModules[i].size() << " MCHits in module " << i << endmsg;
    }
  }

  for ( const auto& hitsInModule : hitsInModules ) {
    std::string strModNum;
    if ( hitsInModule.size() != 0 ) strModNum = std::to_string( hitsInModule[0]->sensDetID() / 4 );
    for ( const auto& mcH : hitsInModule ) {
      // all MCHits
      plot2D( mcH->midPoint().x(), mcH->midPoint().y(), "All MCHit x,y VP module " + strModNum, -60, 60, -60, 60, 240,
              240 );
      plot2D( mcH->midPoint().x(), mcH->midPoint().y(), "All MCHit x,y VP all modules ", -60, 60, -60, 60, 240, 240 );
      // distance particle traveled in sensor & energy
      auto dist = mcH->pathLength();
      plot1D( dist, "Distance travelled by MCHit in VP sensor (mm)", 0, 1, 100 );
      // energy
      auto energy = mcH->energy();
      plot1D( energy, "Energy depotisted by MCHit in VP sensor (MeV)", 0, 1, 100 );
      // were MCHits found
      double efficiency = ( channelIdForMCHit.find( mcH ) != channelIdForMCHit.end() );

      profile2D( mcH->midPoint().x(), mcH->midPoint().y(), efficiency,
                 "Efficiency MCHit -> Pixel x,y VP module " + strModNum, -60, 60, -60, 60, 240, 240 );
      profile2D( mcH->midPoint().x(), mcH->midPoint().y(), efficiency, "Efficiency MCHit -> Pixel x,y VP all modules",
                 -60, 60, -60, 60, 240, 240 );

      profile1D( dist, efficiency, "Efficiency (MCHit->Pixel) v distance travelled by MCHit in VP sensor (mm)", 0, 1,
                 100 );
      profile1D( energy, efficiency, "Efficiency (MCHit->Pixel) v energy depotisted by MCHit in VP sensor (MeV)", 0, 1,
                 100 );
    }
  }

  // what pixels are in the event -- copy them all now, then delete as they are found in clusters
  auto pixUsed = std::set<LHCb::Detector::VPChannelID>( digits.begin(), digits.end() );

  // plot info on clusters
  for ( auto& clus : clusters ) {
    unsigned int           mainPix  = 0;
    unsigned int           otherPix = 0;
    std::set<unsigned int> mcKeysClus;
    for ( auto& channelID : clus.pixels() ) {

      unsigned int numMC = ( MCHitForchannelId.find( channelID.channelID() )->second.size() );
      plot1D( (double)numMC, "Number of MCHits per VP pixel", -0.5, 9.5, 10 );
      // count MCParticles contributing to Cluster
      std::set<unsigned int> mcKeysPix;
      for ( auto& mcHit : ( *MCHitForchannelId.find( channelID.channelID() ) ).second ) {
        auto particle = mcHit->mcParticle();
        // Check if the hit originates from a delta ray.
        while ( particle && particle->originVertex() &&
                particle->originVertex()->type() == LHCb::MCVertex::MCVertexType::DeltaRay ) {
          particle = particle->mother();
        }
        mcKeysPix.insert( particle->key() );
      }
      unsigned int numMCP = mcKeysPix.size();
      plot1D( (double)numMCP, "Number of MCParticles per VP pixel", -0.5, 9.5, 10 );
      plot2D( (double)numMCP, (double)channelID.module(), "Number of MCParticles per VP Pixel v module", -0.5, 9.5,
              -0.5, 51.5, 10, 52 );

      // check if pixel is from MCHit
      if ( numMC != 0 ) {
        ++mainPix;
      } else {
        ++otherPix;
      }
      pixUsed.erase( channelID );                             // leave only "unused" pixels
      for ( auto& key : mcKeysPix ) mcKeysClus.insert( key ); // add to list of cluster MCParticles
    }
    double numMCParticles = (double)mcKeysClus.size();
    plot1D( numMCParticles, "Number of MCParticles per VP Cluster", -0.5, 9.5, 10 );
    plot2D( numMCParticles, (double)clus.channelID().module(), "Number of MCParticles per VP Cluster v module", -0.5,
            9.5, -0.5, 51.5, 10, 52 );
    double fracOther = ( (double)otherPix ) / ( (double)( otherPix + mainPix ) );
    plot1D( fracOther, "Fraction pixels spill or noise per cluster", 0., 1.04, 26 );
    plot2D( fracOther, clus.pixels().size(), "Fraction pixels spill or noise per cluster v cluster size", 0., 1.04, 0.5,
            50.5, 26, 50 );

    // plot positions of "good" clusters i.e. > 70% from MCHits
    unsigned int module    = clus.channelID().module();
    auto         strModNum = std::to_string( module );
    if ( fracOther < 0.3 ) {
      plot2D( clus.x(), clus.y(), "Good (>70% true) Clusters pos all modules", -60, 60, -60, 60, 240, 240 );
      plot2D( clus.x(), clus.y(), "Good (>70% true) Clusters pos module " + strModNum, -60, 60, -60, 60, 240, 240 );
    } else {
      plot2D( clus.x(), clus.y(), "Bad (<70% true) Clusters pos all modules", -60, 60, -60, 60, 240, 240 );
      plot2D( clus.x(), clus.y(), "Bad (<70% true) Clusters pos module " + strModNum, -60, 60, -60, 60, 240, 240 );
    }
  }

  // all pixels
  for ( auto& channelID : digits ) {
    const DeVPSensor& sensor      = det.sensor( channelID.sensor() );
    Gaudi::XYZPoint   pointGlobal = sensor.channelToGlobalPoint( channelID );
    plot2D( pointGlobal.x(), pointGlobal.y(), "Pixel xy all sensors", -60, 60, -60, 60, 240, 240 );
    auto strModNum = std::to_string( channelID.module() );
    plot2D( pointGlobal.x(), pointGlobal.y(), "Pixel xy sensor " + strModNum, -60, 60, -60, 60, 240, 240 );
  }

  // lost pixels
  double nPixNotUsed = (double)( pixUsed.size() );
  if ( msgLevel( MSG::DEBUG ) ) { debug() << "Found " << nPixNotUsed << " VP Digits unused in clusters" << endmsg; }

  plot1D( nPixNotUsed, "Number of pixels not in clusters per event", 0., 10000, 100 );
  if ( digits.size() > 0 ) {
    plot1D( nPixNotUsed / (double)digits.size(), "Fraction of pixels not in clusters per event", 0., 1., 100 );
  }
  for ( auto& channelID : pixUsed ) {
    const DeVPSensor& sensor      = det.sensor( channelID.sensor() );
    Gaudi::XYZPoint   pointGlobal = sensor.channelToGlobalPoint( channelID );
    plot2D( pointGlobal.x(), pointGlobal.y(), "Pixel not in cluster xy all sensors", -60, 60, -60, 60, 240, 240 );
    auto strModNum = std::to_string( channelID.module() );
    plot2D( pointGlobal.x(), pointGlobal.y(), "Pixel not in cluster xy sensor " + strModNum, -60, 60, -60, 60, 240,
            240 );
  }

  // cluster efficiency
  const auto trackInfo = MCTrackInfo{mcprop};

  std::vector<LHCb::MCHit*> hitsReconstructible;
  std::vector<LHCb::MCHit*> hitsMissed;
  std::vector<LHCb::MCHit*> hitsReconstructibleVeloReco;
  std::vector<LHCb::MCHit*> hitsMissedVeloReco;

  std::array<std::vector<LHCb::MCHit*>, 208> hitsMissedInSensors;
  for ( unsigned int i = 0; i < 208; ++i ) hitsMissedInSensors[i].reserve( mcHits.size() / 208 ); // make some space
  std::array<std::vector<LHCb::MCHit*>, 208> hitsReconstructibleInSensors;
  for ( unsigned int i = 0; i < 208; ++i )
    hitsReconstructibleInSensors[i].reserve( mcHits.size() / 208 ); // make some space

  std::array<std::vector<LHCb::MCHit*>, 208> hitsMissedInSensorsMu;
  for ( unsigned int i = 0; i < 208; ++i ) hitsMissedInSensorsMu[i].reserve( mcHits.size() / 208 ); // make some space
  std::array<std::vector<LHCb::MCHit*>, 208> hitsReconstructibleInSensorsMu;
  for ( unsigned int i = 0; i < 208; ++i )
    hitsReconstructibleInSensorsMu[i].reserve( mcHits.size() / 208 ); // make some space

  for ( auto& mcH : mcHits ) {
    if ( channelIdForMCHit.find( mcH ) != channelIdForMCHit.end() ) { // check that hit has at least a digit associated
      hitsReconstructible.push_back( mcH );
      hitsMissed.push_back( mcH );
      hitsMissedInSensors[mcH->sensDetID()].push_back( mcH );
      hitsReconstructibleInSensors[mcH->sensDetID()].push_back( mcH );
      if ( trackInfo.hasVelo( mcH->mcParticle() ) ) { // make residuals for VELO reconstructible tracks
        hitsReconstructibleVeloReco.push_back( mcH );
        hitsMissedVeloReco.push_back( mcH );
      }
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 13 &&
           sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ) < 7 ) {
        hitsReconstructibleInSensorsMu[mcH->sensDetID()].push_back( mcH );
        hitsMissedInSensorsMu[mcH->sensDetID()].push_back( mcH );
      }
    }
  }

  int unmatched_clusters = 0;
  for ( auto& clus : clusters ) {
    bool matched = false;
    for ( auto& channelID : clus.pixels() ) {
      for ( auto& mcHit : ( *MCHitForchannelId.find( channelID.channelID() ) ).second ) {
        matched = true;
        hitsMissed.erase( std::remove( hitsMissed.begin(), hitsMissed.end(), mcHit ), hitsMissed.end() );
        hitsMissedInSensors[mcHit->sensDetID()].erase( std::remove( hitsMissedInSensors[mcHit->sensDetID()].begin(),
                                                                    hitsMissedInSensors[mcHit->sensDetID()].end(),
                                                                    mcHit ),
                                                       hitsMissedInSensors[mcHit->sensDetID()].end() );
        if ( trackInfo.hasVelo( mcHit->mcParticle() ) ) {
          hitsMissedVeloReco.erase( std::remove( hitsMissedVeloReco.begin(), hitsMissedVeloReco.end(), mcHit ),
                                    hitsMissedVeloReco.end() );
        }
        hitsMissedInSensorsMu[mcHit->sensDetID()].erase( std::remove( hitsMissedInSensorsMu[mcHit->sensDetID()].begin(),
                                                                      hitsMissedInSensorsMu[mcHit->sensDetID()].end(),
                                                                      mcHit ),
                                                         hitsMissedInSensorsMu[mcHit->sensDetID()].end() );
      }
    }
    if ( !matched ) { unmatched_clusters += 1; }
  }

  // plots for missed MCHits
  for ( auto& mcH : hitsMissed ) {
    plot1D( mcH->mcParticle()->momentum().Eta(), "# MCHits not found - eta", -5.5, 5.5, 70 );
    plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits not found - phi", -M_PI, M_PI, 21 );
    if ( nullptr != mcH->mcParticle()->originVertex()->mother() ) {
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 321 &&
           abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == 333 ) {
        plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits not found - phi - phikk", -M_PI, M_PI, 21 );
      }
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 11 &&
           abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == 22 ) {
        plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits not found - phi - gammaee", -M_PI, M_PI, 21 );
      }
    }
    plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
            "# MCHits not found - r", 5.1, 50, 50 );
    plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
            "# MCHits not found - r<7mm", 5.1, 7, 14 );
    plot1D( mcH->mcParticle()->p() / 1000.0, "# MCHits not found - p", 0, 100, 120 );
    plot1D( mcH->mcParticle()->pt() / 1000.0, "# MCHits not found - pt", 0, 10, 120 );
    plot1D( mcH->midPoint().z(), "# MCHits not found - z", -300, 800, 90 );
    plot1D( mcH->sensDetID() / 4, "# MCHits not found - module", -0.5, 51.5, 52 );
  }

  // plots for reconstructed MCHits
  for ( auto& mcH : hitsReconstructible ) {
    plot1D( mcH->mcParticle()->momentum().Eta(), "# MCHits reconstructible - eta", -5.5, 5.5, 70 );
    plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits reconstructible - phi", -M_PI, M_PI, 21 );
    if ( nullptr != mcH->mcParticle()->originVertex()->mother() ) {
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 321 &&
           abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == 333 ) {
        plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits reconstructible - phi - phikk", -M_PI, M_PI, 21 );
      }
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 11 &&
           abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == 22 ) {
        plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits reconstructible - phi - gammaee", -M_PI, M_PI, 21 );
      }
    }
    plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
            "# MCHits reconstructible - r", 5.1, 50, 50 );
    plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
            "# MCHits reconstructible - r<7mm", 5.1, 7, 14 );
    plot1D( mcH->mcParticle()->p() / 1000.0, "# MCHits reconstructible - p", 0, 100, 120 );
    plot1D( mcH->mcParticle()->pt() / 1000.0, "# MCHits reconstructible - pt", 0, 10, 120 );
    plot1D( mcH->midPoint().z(), "# MCHits reconstructible - z", -300, 800, 90 );
    plot1D( mcH->sensDetID() / 4, "# MCHits reconstructible - module", -0.5, 51.5, 52 );
  }

  // plots for missed MCHits from Velo reconstructible track
  for ( auto& mcH : hitsMissedVeloReco ) {
    plot1D( mcH->mcParticle()->momentum().Eta(), "# MCHits not found Velo reco - eta", -5.5, 5.5, 70 );
    plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits not found Velo reco - phi", -M_PI, M_PI, 21 );
    if ( nullptr != mcH->mcParticle()->originVertex()->mother() ) {
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 321 &&
           abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == 333 ) {
        plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits not found Velo reco - phi - phikk", -M_PI, M_PI, 21 );
      }
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 11 &&
           abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == 22 ) {
        plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits not found Velo reco - phi - gammaee", -M_PI, M_PI, 21 );
      }
    }
    plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
            "# MCHits not found Velo reco - r", 5.1, 50, 50 );
    plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
            "# MCHits not found Velo reco - r<7mm", 5.1, 7, 14 );
    plot1D( mcH->mcParticle()->p() / 1000.0, "# MCHits not found Velo reco - p", 0, 100, 120 );
    plot1D( mcH->mcParticle()->pt() / 1000.0, "# MCHits not found Velo reco - pt", 0, 10, 120 );
    plot1D( mcH->midPoint().z(), "# MCHits not found Velo reco - z", -300, 800, 90 );
    plot1D( mcH->sensDetID() / 4, "# MCHits not found Velo reco - module", -0.5, 51.5, 52 );
  }

  // plots for reconstructed MCHits from Velo reconstructible track
  for ( auto& mcH : hitsReconstructibleVeloReco ) {
    plot1D( mcH->mcParticle()->momentum().Eta(), "# MCHits reconstructible Velo reco - eta", -5.5, 5.5, 70 );
    plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits reconstructible Velo reco - phi", -M_PI, M_PI, 21 );
    if ( nullptr != mcH->mcParticle()->originVertex()->mother() ) {
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 321 &&
           abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == 333 ) {
        plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits reconstructible Velo reco - phi - phikk", -M_PI, M_PI,
                21 );
      }
      if ( abs( mcH->mcParticle()->particleID().pid() ) == 11 &&
           abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == 22 ) {
        plot1D( mcH->mcParticle()->momentum().phi(), "# MCHits reconstructible Velo reco - phi - gammaee", -M_PI, M_PI,
                21 );
      }
    }
    plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
            "# MCHits reconstructible Velo reco - r", 5.1, 50, 50 );
    plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
            "# MCHits reconstructible Velo reco - r<7mm", 5.1, 7, 14 );
    plot1D( mcH->mcParticle()->p() / 1000.0, "# MCHits reconstructible Velo reco - p", 0, 100, 120 );
    plot1D( mcH->mcParticle()->pt() / 1000.0, "# MCHits reconstructible Velo reco - pt", 0, 10, 120 );
    plot1D( mcH->midPoint().z(), "# MCHits reconstructible Velo reco - z", -300, 800, 90 );
    plot1D( mcH->sensDetID() / 4, "# MCHits reconstructible Velo reco - module", -0.5, 51.5, 52 );
  }

  // plots for reconstructed MCHits from muons
  for ( auto& hitsInSensor : hitsReconstructibleInSensorsMu ) {
    plot1D( hitsInSensor.size(), "# muons r<7mm", -0.5, 10.5, 11 );
    if ( hitsInSensor.size() > 1 ) {
      auto   sensor = hitsInSensor[0]->sensDetID();
      double efficiency;
      if ( hitsInSensor.size() - hitsMissedInSensorsMu[sensor].size() > 1 ) {
        efficiency = 1.0;
      } else {
        efficiency = 0.0;
      }
      profile1D( hitsInSensor.size(), efficiency, "Efficiency vs # #mu (at least 2 #mu)", -0.5, 10.5, 11 );
      plot1D( sensor, "Sensors w at least 2 #mu", -0.5, 207.5, 208 );
    }
  }

  // cluster residual plots
  for ( auto& clus : clusters ) {
    m_num_pix_clu += clus.pixels().size();
    std::set<LHCb::MCHit const*> associated_hits;
    auto                         strModNum = std::to_string( clus.channelID().module() );
    for ( auto& channelID : clus.pixels() ) {
      for ( auto& mcHit : ( *MCHitForchannelId.find( channelID.channelID() ) ).second ) {
        if ( channelIdForMCHit.find( mcHit ) != channelIdForMCHit.end() ) { associated_hits.insert( mcHit ); }
      }
    }
    for ( auto& mcHit : associated_hits ) {
      double x_dist = clus.x() - mcHit->midPoint().x();
      double y_dist = clus.y() - mcHit->midPoint().y();

      if ( abs( x_dist ) < 0.1 ) { m_residual_x += x_dist; }
      if ( abs( y_dist ) < 0.1 ) { m_residual_y += y_dist; }

      plot1D( x_dist, "Residuals along x [mm]", -0.2, 0.2, 200 );
      plot1D( y_dist, "Residuals along y [mm]", -0.2, 0.2, 200 );
      plot1D( x_dist, "Residuals along x [mm] - module" + strModNum, -0.2, 0.2, 200 );
      plot1D( y_dist, "Residuals along y [mm] - module" + strModNum, -0.2, 0.2, 200 );

      if ( trackInfo.hasVelo( mcHit->mcParticle() ) ) { // make residuals for VELO reconstructible tracks
        plot1D( x_dist, "Residuals along x - VELO reco [mm]", -0.2, 0.2, 200 );
        plot1D( y_dist, "Residuals along y - VELO reco [mm]", -0.2, 0.2, 200 );
      }

      if ( trackInfo.hasT( mcHit->mcParticle() ) && trackInfo.hasVelo( mcHit->mcParticle() ) ) { // make residuals for
                                                                                                 // long tracks
        plot1D( x_dist, "Residuals along x - long reco [mm]", -0.2, 0.2, 200 );
        plot1D( y_dist, "Residuals along y - long reco [mm]", -0.2, 0.2, 200 );
      }

      if ( clus.channelID().col() == LHCb::Detector::VPChannelID::ColumnID{255} ||
           clus.channelID().col() == LHCb::Detector::VPChannelID::ColumnID{256} ||
           clus.channelID().col() == LHCb::Detector::VPChannelID::ColumnID{511} ||
           clus.channelID().col() == LHCb::Detector::VPChannelID::ColumnID{512} ) { // make residuals for long pixels
        plot1D( x_dist, "Residuals along x - long pixels [mm]", -0.2, 0.2, 200 );
        plot1D( y_dist, "Residuals along y - long pixels [mm]", -0.2, 0.2, 200 );
      }

      if ( clus.channelID().col() < LHCb::Detector::VPChannelID::ColumnID{3} ||
           clus.channelID().col() > LHCb::Detector::VPChannelID::ColumnID{764} ||
           clus.channelID().row() < LHCb::Detector::VPChannelID::RowID{3} ||
           clus.channelID().row() > LHCb::Detector::VPChannelID::RowID{252} ) { // make residuals for sensor edges
        plot1D( x_dist, "Residuals along x - sensor edges [mm]", -0.2, 0.2, 200 );
        plot1D( y_dist, "Residuals along y - sensor edges [mm]", -0.2, 0.2, 200 );
      }
    }
  }

  // matrix edge check
  for ( const auto& hitsMissedInSensor : hitsMissedInSensors ) {
    for ( const auto& mcH : hitsMissedInSensor ) {
      double x_dist     = 9999;
      double y_dist     = 9999;
      double dist       = 9999;
      double distmcHmcH = 9999;

      for ( auto& clus : clusters ) {
        if ( clus.channelID().sensor() == LHCb::Detector::VPChannelID::SensorID( mcH->sensDetID() ) ) { // Velo sensor
                                                                                                        // are non
                                                                                                        // negative
          if ( sqrt( ( mcH->midPoint().x() - clus.x() ) * ( mcH->midPoint().x() - clus.x() ) +
                     ( mcH->midPoint().y() - clus.y() ) * ( mcH->midPoint().y() - clus.y() ) ) < dist ) {
            dist   = sqrt( ( mcH->midPoint().x() - clus.x() ) * ( mcH->midPoint().x() - clus.x() ) +
                         ( mcH->midPoint().y() - clus.y() ) * ( mcH->midPoint().y() - clus.y() ) );
            x_dist = ( mcH->midPoint().x() - clus.x() );
            y_dist = ( mcH->midPoint().y() - clus.y() );
          }
        }
      }
      if ( x_dist != 9999 && y_dist != 9999 ) {
        plot2D( x_dist, y_dist, "Non reconstructed MCHit - closest cluster", -0.99, 0.99, -0.99, 0.99, 36, 36 );
        plot1D( dist, "Non reconstructed MCHit - closest cluster - 1D", 0, 0.66, 48 );
      }

      for ( auto& mcHReco : hitsReconstructibleInSensors[mcH->sensDetID()] ) {

        if ( ( mcHReco != mcH ) && ( mcHReco->mcParticle() != mcH->mcParticle() ) ) {
          if ( sqrt( ( mcH->midPoint().x() - mcHReco->midPoint().x() ) *
                         ( mcH->midPoint().x() - mcHReco->midPoint().x() ) +
                     ( mcH->midPoint().y() - mcHReco->midPoint().y() ) *
                         ( mcH->midPoint().y() - mcHReco->midPoint().y() ) ) < distmcHmcH ) {
            distmcHmcH = sqrt(
                ( mcH->midPoint().x() - mcHReco->midPoint().x() ) * ( mcH->midPoint().x() - mcHReco->midPoint().x() ) +
                ( mcH->midPoint().y() - mcHReco->midPoint().y() ) * ( mcH->midPoint().y() - mcHReco->midPoint().y() ) );
          }
        }
      }

      if ( distmcHmcH != 9999 && distmcHmcH > 0 ) {
        plot1D( distmcHmcH, "Non reconstructed MCHit - closest reconstructible MCHit", 0, 0.66, 48 );
      }
    }
  }

  for ( const auto& hitsReconstructibleInSensor : hitsReconstructibleInSensors ) {
    for ( const auto& mcH : hitsReconstructibleInSensor ) {
      double distmcHmcH = 9999;

      for ( auto& mcHReco : hitsReconstructibleInSensors[mcH->sensDetID()] ) {
        if ( ( mcHReco != mcH ) && ( mcHReco->mcParticle() != mcH->mcParticle() ) ) {
          if ( sqrt( ( mcH->midPoint().x() - mcHReco->midPoint().x() ) *
                         ( mcH->midPoint().x() - mcHReco->midPoint().x() ) +
                     ( mcH->midPoint().y() - mcHReco->midPoint().y() ) *
                         ( mcH->midPoint().y() - mcHReco->midPoint().y() ) ) < distmcHmcH ) {
            distmcHmcH = sqrt(
                ( mcH->midPoint().x() - mcHReco->midPoint().x() ) * ( mcH->midPoint().x() - mcHReco->midPoint().x() ) +
                ( mcH->midPoint().y() - mcHReco->midPoint().y() ) * ( mcH->midPoint().y() - mcHReco->midPoint().y() ) );
          }
        }
      }
      if ( distmcHmcH != 9999 && distmcHmcH > 0 ) {
        plot1D( distmcHmcH, "Reconstructible MCHit - closest reconstructible MCHit", 0, 0.66, 48 );
        plot1D( distmcHmcH, "Reconstructible MCHit - closest reconstructible MCHit - extended", 0, 52.8, 3840 );
      }
    }
  }

  // matrix edge check - phi->kk
  constexpr auto KPlus    = 321;
  constexpr auto phi_1020 = 333;
  for ( const auto& hitsMissedInSensor : hitsMissedInSensors ) {
    for ( const auto& mcHMiss : hitsMissedInSensor ) {
      double dist = 9999;
      if ( nullptr != mcHMiss->mcParticle()->originVertex()->mother() ) {
        if ( abs( mcHMiss->mcParticle()->particleID().pid() ) == KPlus &&
             abs( mcHMiss->mcParticle()->originVertex()->mother()->particleID().pid() ) == phi_1020 ) {
          for ( auto& mcH : hitsReconstructibleInSensors[mcHMiss->sensDetID()] ) {
            if ( mcH != mcHMiss && nullptr != mcH->mcParticle()->originVertex()->mother() ) {
              if ( abs( mcH->mcParticle()->particleID().pid() ) == KPlus &&
                   abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == phi_1020 ) {
                if ( sqrt( ( mcH->midPoint().x() - mcHMiss->midPoint().x() ) *
                               ( mcH->midPoint().x() - mcHMiss->midPoint().x() ) +
                           ( mcH->midPoint().y() - mcHMiss->midPoint().y() ) *
                               ( mcH->midPoint().y() - mcHMiss->midPoint().y() ) ) < dist ) {
                  dist = sqrt( ( mcH->midPoint().x() - mcHMiss->midPoint().x() ) *
                                   ( mcH->midPoint().x() - mcHMiss->midPoint().x() ) +
                               ( mcH->midPoint().y() - mcHMiss->midPoint().y() ) *
                                   ( mcH->midPoint().y() - mcHMiss->midPoint().y() ) );
                }
              }
            }
          }
          if ( dist != 9999 && dist > 0 ) {
            plot1D( dist, "Non reconstructed MCHit - closest reconstructible MCHit - phikk", 0, 0.66, 48 );
          }
        }
      }
    }
  }

  for ( const auto& hitsReconstructibleInSensor : hitsReconstructibleInSensors ) {
    for ( const auto& mcHReco : hitsReconstructibleInSensor ) {
      double distmcHmcH = 9999;
      if ( nullptr != mcHReco->mcParticle()->originVertex()->mother() ) {
        if ( abs( mcHReco->mcParticle()->particleID().pid() ) == KPlus &&
             abs( mcHReco->mcParticle()->originVertex()->mother()->particleID().pid() ) == phi_1020 ) {
          for ( auto& mcH : hitsReconstructibleInSensors[mcHReco->sensDetID()] ) {
            if ( mcHReco != mcH && nullptr != mcH->mcParticle()->originVertex()->mother() ) {
              if ( abs( mcH->mcParticle()->particleID().pid() ) == KPlus &&
                   abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == phi_1020 ) {
                if ( sqrt( ( mcH->midPoint().x() - mcHReco->midPoint().x() ) *
                               ( mcH->midPoint().x() - mcHReco->midPoint().x() ) +
                           ( mcH->midPoint().y() - mcHReco->midPoint().y() ) *
                               ( mcH->midPoint().y() - mcHReco->midPoint().y() ) ) < distmcHmcH ) {
                  distmcHmcH = sqrt( ( mcH->midPoint().x() - mcHReco->midPoint().x() ) *
                                         ( mcH->midPoint().x() - mcHReco->midPoint().x() ) +
                                     ( mcH->midPoint().y() - mcHReco->midPoint().y() ) *
                                         ( mcH->midPoint().y() - mcHReco->midPoint().y() ) );
                }
              }
            }
          }
          if ( distmcHmcH != 9999 && distmcHmcH > 0 ) {
            plot1D( distmcHmcH, "Reconstructible MCHit - closest reconstructible MCHit - phikk", 0, 0.66, 48 );
          }
        }
      }
    }
  }

  // matrix edge check - gamma->ee
  constexpr auto posit = 11;
  constexpr auto gamma = 22;
  for ( const auto& hitsMissedInSensor : hitsMissedInSensors ) {
    for ( const auto& mcHMiss : hitsMissedInSensor ) {
      double dist = 9999;
      if ( nullptr != mcHMiss->mcParticle()->originVertex()->mother() ) {
        if ( abs( mcHMiss->mcParticle()->particleID().pid() ) == posit &&
             abs( mcHMiss->mcParticle()->originVertex()->mother()->particleID().pid() ) == gamma ) {
          for ( auto& mcH : hitsReconstructibleInSensors[mcHMiss->sensDetID()] ) {
            if ( mcH != mcHMiss && nullptr != mcH->mcParticle()->originVertex()->mother() ) {
              if ( abs( mcH->mcParticle()->particleID().pid() ) == posit &&
                   abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == gamma ) {
                if ( sqrt( ( mcH->midPoint().x() - mcHMiss->midPoint().x() ) *
                               ( mcH->midPoint().x() - mcHMiss->midPoint().x() ) +
                           ( mcH->midPoint().y() - mcHMiss->midPoint().y() ) *
                               ( mcH->midPoint().y() - mcHMiss->midPoint().y() ) ) < dist ) {
                  dist = sqrt( ( mcH->midPoint().x() - mcHMiss->midPoint().x() ) *
                                   ( mcH->midPoint().x() - mcHMiss->midPoint().x() ) +
                               ( mcH->midPoint().y() - mcHMiss->midPoint().y() ) *
                                   ( mcH->midPoint().y() - mcHMiss->midPoint().y() ) );
                }
              }
            }
          }
          if ( dist != 9999 && dist > 0 ) {
            plot1D( dist, "Non reconstructed MCHit - closest reconstructible MCHit - gammaee", 0, 0.66, 48 );
          }
        }
      }
    }
  }

  for ( const auto& hitsReconstructibleInSensor : hitsReconstructibleInSensors ) {
    for ( const auto& mcHReco : hitsReconstructibleInSensor ) {
      double distmcHmcH = 9999;
      if ( nullptr != mcHReco->mcParticle()->originVertex()->mother() ) {
        if ( abs( mcHReco->mcParticle()->particleID().pid() ) == posit &&
             abs( mcHReco->mcParticle()->originVertex()->mother()->particleID().pid() ) == gamma ) {
          for ( auto& mcH : hitsReconstructibleInSensors[mcHReco->sensDetID()] ) {
            if ( mcHReco != mcH && nullptr != mcH->mcParticle()->originVertex()->mother() ) {
              if ( abs( mcH->mcParticle()->particleID().pid() ) == posit &&
                   abs( mcH->mcParticle()->originVertex()->mother()->particleID().pid() ) == gamma ) {
                if ( sqrt( ( mcH->midPoint().x() - mcHReco->midPoint().x() ) *
                               ( mcH->midPoint().x() - mcHReco->midPoint().x() ) +
                           ( mcH->midPoint().y() - mcHReco->midPoint().y() ) *
                               ( mcH->midPoint().y() - mcHReco->midPoint().y() ) ) < distmcHmcH ) {
                  distmcHmcH = sqrt( ( mcH->midPoint().x() - mcHReco->midPoint().x() ) *
                                         ( mcH->midPoint().x() - mcHReco->midPoint().x() ) +
                                     ( mcH->midPoint().y() - mcHReco->midPoint().y() ) *
                                         ( mcH->midPoint().y() - mcHReco->midPoint().y() ) );
                }
              }
            }
          }
          if ( distmcHmcH != 9999 && distmcHmcH > 0 ) {
            plot1D( distmcHmcH, "Reconstructible MCHit - closest reconstructible MCHit - gammaee", 0, 0.66, 48 );
          }
        }
      }
    }
  }

  // plot distribution number of clusters
  m_num_clusters += clusters.size();
  plot1D( clusters.size(), "# clusters", 0, 10000, 10000 );
  plot1D( unmatched_clusters, "# unmatched clusters", 0, 10000, 100 );
  for ( auto& clus : clusters ) {
    auto strModNum = std::to_string( clus.channelID().module() );
    plot1D( sqrt( clus.x() * clus.x() + clus.y() * clus.y() ), "# cluster distribution - r", 5.1, 50, 50 );
    plot1D( sqrt( clus.x() * clus.x() + clus.y() * clus.y() ), "# cluster distribution - r - module" + strModNum, 5.1,
            50, 50 );
    plot1D( clus.channelID().module(), "# cluster distribution - module", -0.5, 51.5, 52 );
    plot1D( clus.pixels().size(), "# pixel per cluster distribution", -0.5, 50.5, 51 );
    plot1D( clus.pixels().size(), "# pixel per cluster distribution - module" + strModNum, -0.5, 50.5, 51 );
  }

  // purity plots
  for ( auto& clus : clusters ) {
    std::set<LHCb::MCHit const*> associated_hits;
    std::vector<unsigned int>    ids_clu;
    for ( auto& channelID : clus.pixels() ) {
      ids_clu.push_back( channelID.channelID() );
      for ( auto& mcHit : ( *MCHitForchannelId.find( channelID.channelID() ) ).second ) {
        associated_hits.insert( mcHit );
      }
    }
    for ( auto& mcHit : associated_hits ) {
      std::vector<unsigned int> ids_hit;
      for ( auto& pix : ( *channelIdForMCHit.find( mcHit ) ).second ) { ids_hit.push_back( pix ); }
      std::sort( ids_clu.begin(), ids_clu.end() );
      std::sort( ids_hit.begin(), ids_hit.end() );
      std::vector<int> common;
      set_intersection( ids_clu.begin(), ids_clu.end(), ids_hit.begin(), ids_hit.end(), back_inserter( common ) );
      plot1D( ( (double)common.size() ) / ids_hit.size(), "Purity distribution", 0, 1.1, 110 );
      m_purity += ( (double)common.size() ) / ids_hit.size();
      plot1D( ids_hit.size(), "# pixel distribution - seen", 0, 100, 20 );
      profile1D( ids_hit.size(), ( (double)common.size() ) / ids_hit.size(), "Purity vs # pixels", 0, 100, 20 );
    }
    plot1D( ids_clu.size(), "Number of pixels in cluster", 0, 32, 33 );
  }

  // pixel stat plots
  for ( auto& mcHit : hitsReconstructible ) {
    std::set<int> associated_clusters;
    for ( auto& clus : clusters ) {
      if ( LHCb::Detector::VPChannelID::SensorID( mcHit->sensDetID() ) == clus.channelID().sensor() ) {
        for ( auto& channelID : clus.pixels() ) {
          for ( auto& mcHit_ass : ( *MCHitForchannelId.find( channelID.channelID() ) ).second ) {
            if ( mcHit_ass == mcHit ) { associated_clusters.insert( clus.channelID().channelID() ); }
          }
        }
      }
    }
    plot1D( associated_clusters.size(), "# associated cluster to hit", -0.5, 10.5, 11 );
    plot1D( channelIdForMCHit.find( mcHit )->second.size(), "# pixel distribution", 0, 100, 20 );
    m_efficiency += !associated_clusters.empty();
    m_num_pix_hit += ( channelIdForMCHit.find( mcHit )->second.size() );
  }

  // efficiency vs occupancy
  for ( const auto& hitsReconstructibleInSensor : hitsReconstructibleInSensors ) {
    for ( const auto& mcH : hitsReconstructibleInSensor ) {
      // reconstructible MCHits
      const auto sensor = mcH->sensDetID();
      plot2D( mcH->midPoint().x(), mcH->midPoint().y(),
              "Reconstructible MCHit x,y VP sensor " + std::to_string( sensor ), -60, 60, -60, 60, 240, 240 );
      plot2D( mcH->midPoint().x(), mcH->midPoint().y(), "Reconstructible MCHit x,y VP all sensors ", -60, 60, -60, 60,
              240, 240 );
      if ( trackInfo.hasVelo( mcH->mcParticle() ) ) {
        plot2D( mcH->midPoint().x(), mcH->midPoint().y(),
                "Reconstructible MCHit x,y VP sensor - VELO " + std::to_string( sensor ), -60, 60, -60, 60, 240, 240 );
        plot2D( mcH->midPoint().x(), mcH->midPoint().y(), "Reconstructible MCHit x,y VP all sensors - VELO ", -60, 60,
                -60, 60, 240, 240 );
      }
    }
  }

  for ( const auto& hitsMissedInSensor : hitsMissedInSensors ) {

    std::string strSenNum;
    if ( hitsMissedInSensor.size() != 0 ) strSenNum = std::to_string( hitsMissedInSensor[0]->sensDetID() );

    for ( auto& mcH : hitsMissedInSensor ) {
      // missed MCHits
      plot2D( mcH->midPoint().x(), mcH->midPoint().y(), "Missed MCHit x,y VP sensor " + strSenNum, -60, 60, -60, 60,
              240, 240 );
      plot2D( mcH->midPoint().x(), mcH->midPoint().y(), "Missed MCHit x,y VP all sensors ", -60, 60, -60, 60, 240,
              240 );
      if ( trackInfo.hasVelo( mcH->mcParticle() ) ) {
        plot2D( mcH->midPoint().x(), mcH->midPoint().y(), "Missed MCHit x,y VP sensor - VELO " + strSenNum, -60, 60,
                -60, 60, 240, 240 );
        plot2D( mcH->midPoint().x(), mcH->midPoint().y(), "Missed MCHit x,y VP all sensors - VELO ", -60, 60, -60, 60,
                240, 240 );
      }
    }
  }

  for ( const auto& hitsReconstructibleInSensor : hitsReconstructibleInSensors ) {

    std::string strModNum;
    if ( hitsReconstructibleInSensor.size() != 0 )
      strModNum = std::to_string( hitsReconstructibleInSensor[0]->sensDetID() / 4 );

    for ( const auto& mcH : hitsReconstructibleInSensor ) {
      // reconstructible MCHits
      plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
              "# MCHits reconstructible - r - module" + strModNum, 5.1, 50, 50 );
      plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
              "# MCHits reconstructible - r<7mm - module" + strModNum, 5.1, 7, 14 );
      if ( trackInfo.hasVelo( mcH->mcParticle() ) ) {
        plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
                "# MCHits reconstructible Velo reco - r - module" + strModNum, 5.1, 50, 50 );
        plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
                "# MCHits reconstructible Velo reco - r<7mm - module" + strModNum, 5.1, 7, 14 );
      }
    }
  }

  for ( const auto& hitsMissedInSensor : hitsMissedInSensors ) {

    std::string strModNum;
    if ( hitsMissedInSensor.size() != 0 ) strModNum = std::to_string( hitsMissedInSensor[0]->sensDetID() / 4 );

    for ( auto& mcH : hitsMissedInSensor ) {
      // missed MCHits
      plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
              "# MCHits not found - r - module" + strModNum, 5.1, 50, 50 );
      plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
              "# MCHits not found - r<7mm - module" + strModNum, 5.1, 7, 14 );
      if ( trackInfo.hasVelo( mcH->mcParticle() ) ) {
        plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
                "# MCHits not found Velo reco - r - module" + strModNum, 5.1, 50, 50 );
        plot1D( sqrt( mcH->midPoint().x() * mcH->midPoint().x() + mcH->midPoint().y() * mcH->midPoint().y() ),
                "# MCHits not found Velo reco - r<7mm - module" + strModNum, 5.1, 7, 14 );
      }
    }
  }
}
