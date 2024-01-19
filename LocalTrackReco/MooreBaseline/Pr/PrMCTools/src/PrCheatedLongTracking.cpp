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
#include "Event/FTLiteCluster.h"
#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/MCVertex.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Event/TrackTags.h"
#include "LHCbAlgs/Transformer.h"
#include "MCInterfaces/IIdealStateCreator.h"
#include "PrKernel/PrFTHitHandler.h"

using Track = LHCb::Event::v2::Track;

/** @class PrCheatedLongTracking
 *
 *  Ideal pattern recognition for long (v2) tracks. Uses the links created by the
 *  PrLHCbID2MCParticle algorithm to construct tracks consisting only out of LHCbIDs
 *  associated to the MCParticle.
 *
 *  Hits in the VP are added, and a possibility exists to add hits in the UT
 *  as well, given a certain threshold (e.g. at least 4 hits).
 *
 *  Following the cheated PR of the SciFi: cuts can be set on the minimum number
 *  of total hits, hits in x layers and hits in stereo layers.
 *  Top and bottom modules are not mixed in the x layers.
 *
 *  Beware: This just produces a "container of LHCbIDs", and no directional information
 *  of the track, i.e. no meaningful state covariance matrix.
 *  The track parameters at creation are filled in the closest to beam state.
 *
 *  The standard use-case for this is the tuning of the track fit,
 *  which should result in well estimated uncertainties with this ideal
 *  PR.
 *
 *  Note that an explicit dependency exists on MCVertices, while it's not
 *  being used directly. This is the result if the indirect use through
 *  MCParticles: we call MCParticle::originVertex().
 *
 *  @todo Ideally we would merge all CheatedPR algorithms.
 *
 *  @author Laurent Dufour
 */
class PrCheatedLongTracking : public LHCb::Algorithm::Transformer<std::vector<Track>(
                                  const PrFTHitHandler<PrHit>&, const LHCb::MCParticles&, const LHCb::MCVertices&,
                                  const LHCb::MCProperty&, const LHCb::LinksByKey& )> {
public:
  /// Standard constructor
  PrCheatedLongTracking( const std::string& name, ISvcLocator* pSvcLocator );

  /// make cheated tracks by getting the clusters matched to an MCParticle
  std::vector<Track> operator()( const PrFTHitHandler<PrHit>&, const LHCb::MCParticles&, const LHCb::MCVertices&,
                                 const LHCb::MCProperty&, const LHCb::LinksByKey& ) const override;

private:
  Gaudi::Property<int> m_numFTZones      = {this, "NumFTZones", 24};
  Gaudi::Property<int> m_minFTXHits      = {this, "MinFTXHits", 3};
  Gaudi::Property<int> m_minFTStereoHits = {this, "MinFTStereoHits", 3};
  Gaudi::Property<int> m_minFTHits       = {this, "MinFTHits", 6};

  // Minimum number of VP clusters to be found for the track to be created
  Gaudi::Property<size_t> m_minVPHits = {this, "MinVPHits", 3};

  // UT Hits are only added if the number of associated hits is great or equal than this
  Gaudi::Property<size_t> m_thresholdUTHits = {this, "ThresholdUTHits", 4};

  // UT Hits are only added if this is true
  Gaudi::Property<bool> m_addUTHits = {this, "AddUTHits", true};

  Gaudi::Property<bool>          m_add_ideal_states = {this, "AddIdealStates", false};
  ToolHandle<IIdealStateCreator> m_ideal_state_creator{this, "IdealStateCreator", "IdealStateCreator"};
};

DECLARE_COMPONENT( PrCheatedLongTracking )

PrCheatedLongTracking::PrCheatedLongTracking( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator,
                   {
                       KeyValue{"FTHitsLocation", PrFTInfo::FTHitsLocation},
                       KeyValue{"MCParticleLocation", LHCb::MCParticleLocation::Default},
                       KeyValue{"MCVerticesLocation", LHCb::MCVertexLocation::Default},
                       KeyValue{"MCPropertyLocation", LHCb::MCPropertyLocation::TrackInfo},
                       KeyValue{"LHCbIdLinkLocation", "Link/Pr/LHCbID"},
                   },
                   KeyValue{"OutputName", "Rec/Track/CheatedLong"} ) {}

std::vector<Track> PrCheatedLongTracking::operator()( const PrFTHitHandler<PrHit>& FTHitHandler,
                                                      const LHCb::MCParticles&     mcParts, const LHCb::MCVertices&,
                                                      const LHCb::MCProperty&      mcProps,
                                                      const LHCb::LinksByKey&      links ) const {
  std::vector<Track> result;
  MCTrackInfo        trackInfo( mcProps );

  for ( const LHCb::MCParticle* mcPart : mcParts ) {
    const bool isLong = trackInfo.hasVeloAndT( mcPart );
    if ( !isLong ) continue;

    Track newTrack;
    newTrack.setType( Track::Type::Long );
    newTrack.setHistory( Track::History::TrackIdealPR );
    newTrack.setPatRecStatus( Track::PatRecStatus::PatRecIDs );
    std::vector<LHCb::LHCbID> utIds;
    const bool                isAddUTHits = m_addUTHits.value();
    size_t                    nVPHits     = 0;

    // add the VP hits & UT hits
    links.applyToAllLinks( [&mcParts, &mcPart, &newTrack, &utIds, &isAddUTHits,
                            &nVPHits]( unsigned int srcKey, unsigned int mcPartKey, float ) {
      const LHCb::MCParticle* linkedMCPart =
          static_cast<const LHCb::MCParticle*>( mcParts.containedObject( mcPartKey ) );
      LHCb::LHCbID theId( srcKey );

      if ( mcPart == linkedMCPart ) {
        if ( theId.isVP() ) {
          newTrack.addToLhcbIDs( theId );
          ++nVPHits;
        }

        if ( theId.isUT() && isAddUTHits ) { utIds.push_back( theId ); }
      }
    } );

    if ( nVPHits < m_minVPHits.value() ) continue;

    /**
     * Addition of FT hits. Makes use of the FT detector information
     * (zones) to ensure bottom and top modules are not mixed.
     **/
    std::vector<int> firedXLayers( m_numFTZones.value(), 0 );
    std::vector<int> firedStereoLayers( m_numFTZones.value(), 0 );
    int              ftHits = 0;

    // -- loop over all zones
    for ( int iZone = 0; iZone < m_numFTZones.value(); ++iZone ) {
      // -- loop over all hits in a zone
      for ( const auto& hit : FTHitHandler.hits( iZone ) ) {
        bool found = false;

        links.applyToLinks( hit.id().lhcbID(), [&]( unsigned int, unsigned int index, float ) {
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
          ftHits++;
          newTrack.addToLhcbIDs( hit.id() );
        }
      }
    }

    int sumLowerX = 0;
    int sumUpperX = 0;
    int sumStereo = 0;
    for ( int i = 0; i < m_numFTZones.value(); i += 2 ) { sumLowerX += firedXLayers[i]; }
    for ( int i = 1; i < m_numFTZones.value(); i += 2 ) { sumUpperX += firedXLayers[i]; }
    for ( int i = 0; i < m_numFTZones.value(); i++ ) { sumStereo += firedStereoLayers[i]; }
    debug() << "sumLowerX: " << sumLowerX << " sumUpperX " << sumUpperX << " sumStereo " << sumStereo << " totHits "
            << ftHits << endmsg;

    if ( ( sumLowerX < m_minFTXHits.value() && sumUpperX < m_minFTXHits.value() ) ||
         sumStereo < m_minFTStereoHits.value() || ftHits < m_minFTHits.value() ) {
      continue; // Don't even make the track if the number of FT hits is too little.
    }

    /**
     * Add UT hits
     **/
    if ( utIds.size() >= m_thresholdUTHits.value() ) {
      for ( const auto& utId : utIds ) newTrack.addToLhcbIDs( utId );
    }

    const double qOverP = ( mcPart->particleID().threeCharge() / 3 ) / mcPart->p();
    const double x      = mcPart->originVertex()->position().X();
    const double y      = mcPart->originVertex()->position().Y();
    const double z      = mcPart->originVertex()->position().Z();
    const double tx     = mcPart->momentum().X() / mcPart->momentum().Z();
    const double ty     = mcPart->momentum().Y() / mcPart->momentum().Z();

    LHCb::State stateClosestToBeam;
    stateClosestToBeam.setLocation( LHCb::State::Location::ClosestToBeam );
    stateClosestToBeam.setState( x, y, z, tx, ty, qOverP );
    newTrack.addToStates( stateClosestToBeam );

    if ( m_add_ideal_states ) {
      std::vector<LHCb::State> newstates{};
      newstates.reserve( 30 );
      m_ideal_state_creator->getMCHitStates( *mcPart, newstates ).ignore();
      newTrack.addToStates( newstates, LHCb::Tag::Unordered_tag{} );
    }
    result.emplace_back( newTrack );
  }

  return result;
}
