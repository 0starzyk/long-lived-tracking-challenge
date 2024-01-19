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
// Include files

// from Gaudi
#include "LHCbAlgs/Transformer.h"

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/FTLiteCluster.h"
#include "Event/LinksByKey.h"
#include "Event/MCTrackInfo.h"
#include "Event/PrSeedTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "PrKernel/PrFTHitHandler.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrCheatedSciFiTracking
//
// 2015-03-23 : Michel De Cian
//-----------------------------------------------------------------------------

/** @class PrCheatedSciFiTracking PrCheatedSciFiTracking.h
 *
 *  Cheated track reconstruction in the SciFi.
 *  It creates tracks by getting all Clusters associated to a reconstructible MCParticle.
 *  Cuts can be set on the minimum number of total hits, hits in x layers and hits in stereo layers.
 *  Top and bottom modules are not mixed in the x layers.
 *  Beware: This just produces a "container of LHCbIDs", and no directional information of the track, i.e. no meaningful
 * state.
 *
 * - NumZones: Number of zones (normally 2 x number of layers if no y-segmentation)
 * - MinXHits: Minimum number of required hits in x layers to make a track.
 * - MinStereoHits: Minimum number of required hits in stereo layers to make a track.
 * - MinTotHits: Minimum number of total hits to make a track.
 *
 *
 *  @author Michel De Cian
 *  @date   2015-03-23
 */

using SeedTracks = LHCb::Pr::Seeding::Tracks;
using SeedTag    = LHCb::Pr::Seeding::Tag;

class PrCheatedSciFiTracking
    : public LHCb::Algorithm::Transformer<SeedTracks( const PrFTHitHandler<PrHit>&, const LHCb::MCParticles&,
                                                      const LHCb::MCProperty&, const LHCb::LinksByKey&,
                                                      DetectorElement const& ),
                                          LHCb::DetDesc::usesConditions<DetectorElement>> {
public:
  using base_t = LHCb::Algorithm::Transformer<SeedTracks( const PrFTHitHandler<PrHit>&, const LHCb::MCParticles&,
                                                          const LHCb::MCProperty&, const LHCb::LinksByKey&,
                                                          DetectorElement const& ),
                                              LHCb::DetDesc::usesConditions<DetectorElement>>;

  /// Standard constructor
  PrCheatedSciFiTracking( const std::string& name, ISvcLocator* pSvcLocator );

  /// make cheated tracks by getting the clusters matched to an MCParticle
  SeedTracks operator()( const PrFTHitHandler<PrHit>&, const LHCb::MCParticles&, const LHCb::MCProperty&,
                         const LHCb::LinksByKey&, DetectorElement const& ) const override;

private:
  Gaudi::Property<int> m_numZones      = {this, "NumZones", 24};
  Gaudi::Property<int> m_minXHits      = {this, "MinXHits", 5};
  Gaudi::Property<int> m_minStereoHits = {this, "MinStereoHits", 5};
  Gaudi::Property<int> m_minTotHits    = {this, "MinTotHits", 10};

  /// The track extrapolator
  ToolHandle<ITrackExtrapolator> m_extrapolator = {this, "ReferenceExtrapolator", "TrackMasterExtrapolator"};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( PrCheatedSciFiTracking )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrCheatedSciFiTracking::PrCheatedSciFiTracking( const std::string& name, ISvcLocator* pSvcLocator )
    : base_t( name, pSvcLocator,
              {KeyValue{"FTHitsLocation", PrFTInfo::FTHitsLocation},
               KeyValue{"MCParticleLocation", LHCb::MCParticleLocation::Default},
               KeyValue{"MCPropertyLocation", LHCb::MCPropertyLocation::TrackInfo},
               KeyValue{"LinkLocation", "Link/Raw/FT/LiteClusters"},
               KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}},
              KeyValue{"OutputName", LHCb::TrackLocation::Seed} ) {}

//=============================================================================
// make cheated tracks by getting the clusters matched to an MCParticle
//=============================================================================
SeedTracks PrCheatedSciFiTracking::operator()( const PrFTHitHandler<PrHit>& FTHitHandler,
                                               const LHCb::MCParticles& mcParts, const LHCb::MCProperty& mcProps,
                                               const LHCb::LinksByKey& links, DetectorElement const& lhcb ) const {
  SeedTracks result;

  MCTrackInfo trackInfo( mcProps );

  for ( const LHCb::MCParticle* mcPart : mcParts ) {

    const bool isSeed = trackInfo.hasT( mcPart );
    if ( !isSeed ) continue;

    auto outTrack = result.emplace_back<SIMDWrapper::InstructionSet::Scalar>();

    std::vector<int> firedXLayers( m_numZones.value(), 0 );
    std::vector<int> firedStereoLayers( m_numZones.value(), 0 );
    int              totHits = 0;

    // -- loop over all zones
    const double qOverP = ( mcPart->particleID().threeCharge() / 3 ) / mcPart->p();
    const double x      = mcPart->originVertex()->position().X();
    const double y      = mcPart->originVertex()->position().Y();
    const double z      = mcPart->originVertex()->position().Z();
    const double tx     = mcPart->momentum().X() / mcPart->momentum().Z();
    const double ty     = mcPart->momentum().Y() / mcPart->momentum().Z();
    LHCb::State  oState;
    oState.setLocation( LHCb::State::Location::ClosestToBeam );
    oState.setState( x, y, z, tx, ty, qOverP );
    for ( int iZone = 0; iZone < m_numZones.value(); ++iZone ) {
      // -- loop over all hits in a zone
      for ( const auto& hit : FTHitHandler.hits( iZone ) ) {
        bool found = false;

        links.applyToLinks( hit.id().ftID(), [&]( unsigned int, unsigned int index, float ) {
          const LHCb::MCParticle* linkedMCPart =
              static_cast<const LHCb::MCParticle*>( mcParts.containedObject( index ) );
          if ( mcPart == linkedMCPart ) found = true;
        } );

        if ( hit.isX() && found ) {
          if ( firedXLayers[iZone] == 0 ) firedXLayers[iZone]++;
        }
        if ( !( hit.isX() ) && found ) {
          if ( firedStereoLayers[iZone] == 0 ) firedStereoLayers[iZone]++;
        }
        if ( found ) {
          if ( totHits < static_cast<int>( LHCb::Pr::TracksInfo::MaxFTHits ) ) {
            outTrack.field<SeedTag::FTHits>().resize( totHits + 1 );
            outTrack.field<SeedTag::FTHits>()[totHits].template field<SeedTag::Index>().set( 0 ); // TODO
            outTrack.field<SeedTag::FTHits>()[totHits].template field<SeedTag::LHCbID>().set(
                LHCb::Event::lhcbid_v<SIMDWrapper::scalar::types>( hit.id().lhcbID() ) );
          }
          totHits++;
        }
      }
    }

    int iState = 0;
    for ( auto z : {StateParameters::ZBegT, StateParameters::ZMidT, StateParameters::ZEndT} ) {
      auto tState = outTrack.field<SeedTag::States>( iState );
      tState.setQOverP( static_cast<float>( qOverP ) );
      LHCb::State state = oState;
      StatusCode  sc    = m_extrapolator->propagate( state, z, *lhcb.geometry() );

      iState++;
      if ( !sc ) {
        tState.setPosition( 0.0f, 0.0f, static_cast<float>( z ) );
        tState.setDirection( 0.0f, 0.0f );
        debug() << "extrapolation failed. p: " << std::abs( 1 / qOverP ) << endmsg;
        continue;
      }

      tState.setPosition( static_cast<float>( state.x() ), static_cast<float>( state.y() ), static_cast<float>( z ) );
      tState.setDirection( static_cast<float>( state.tx() ), static_cast<float>( state.ty() ) );
    }

    outTrack.field<SeedTag::Chi2PerDoF>().set( 1.0f );

    int sumLowerX = 0;
    int sumUpperX = 0;
    int sumStereo = 0;
    for ( int i = 0; i < m_numZones.value(); i += 2 ) { sumLowerX += firedXLayers[i]; }
    for ( int i = 1; i < m_numZones.value(); i += 2 ) { sumUpperX += firedXLayers[i]; }
    for ( int i = 0; i < m_numZones.value(); i++ ) { sumStereo += firedStereoLayers[i]; }
    debug() << "sumLowerX: " << sumLowerX << " sumUpperX " << sumUpperX << " sumStereo " << sumStereo << " totHits "
            << totHits << endmsg;

    if ( ( sumLowerX < m_minXHits.value() && sumUpperX < m_minXHits.value() ) || sumStereo < m_minStereoHits.value() ||
         totHits < m_minTotHits.value() ) {
      continue;
    }
  }
  return result;
}
