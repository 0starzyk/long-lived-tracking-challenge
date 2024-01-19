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
#include "AIDA/IHistogram1D.h"
#include "Event/GhostTrackInfo.h"
#include "Event/MCParticle.h"
#include "Event/Track.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "LHCbAlgs/Consumer.h"
#include "Map.h"
#include "TrackCheckerBase.h"

/** @class TrackEffChecker TrackEffChecker.h
 *
 * Class for track monitoring
 *  @author M. Needham.
 *  @date   6-5-2007
 */

class TrackEffChecker : public LHCb::Algorithm::Consumer<void( const LHCb::Tracks&, const LHCb::MCParticles&,
                                                               const LHCb::LinksByKey&, const LHCb::LinksByKey& ),
                                                         Gaudi::Functional::Traits::BaseClass_t<TrackCheckerBase>> {

public:
  /** Standard construtor */
  TrackEffChecker( const std::string& name, ISvcLocator* pSvcLocator );

  /** Algorithm execute */
  void operator()( const LHCb::Tracks&, const LHCb::MCParticles&, const LHCb::LinksByKey&,
                   const LHCb::LinksByKey& ) const override;

  /** Algorithm finalize */
  StatusCode finalize() override;

private:
  Gaudi::Property<bool> m_requireLong{this, "RequireLongTrack", false};

  void ghostInfo( const LHCb::Tracks&, const LHCb::MCParticles&, const LHCb::LinksByKey&, unsigned int ) const;

  void effInfo( const LHCb::Tracks&, const LHCb::MCParticles&, const LHCb::LinksByKey&,
                const std::vector<std::vector<int>>& ) const;

  void plots( const std::string& type, const LHCb::Track* track ) const;

  void plots( const std::string& type, const LHCb::MCParticle* part ) const;

  double weightedMeasurementSum( const LHCb::Track* aTrack ) const;

  mutable Gaudi::Accumulators::StatCounter<> m_nTrackCounter{this, "nTrack"};
  mutable Gaudi::Accumulators::StatCounter<> m_nGhostCounter{this, "nGhost"};

  mutable Gaudi::Accumulators::StatCounter<> m_nClone{this, "nClone"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneG5{this, "nCloneG5"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneB{this, "nCloneB"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneG5B{this, "nCloneG5B"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneBAll{this, "nCloneBAll"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneG5BAll{this, "nCloneG5BAll"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneKsL{this, "nCloneKsL"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneG5KsL{this, "nCloneG5KsL"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneIP{this, "nCloneIP"};
  mutable Gaudi::Accumulators::StatCounter<> m_nCloneG5IP{this, "nCloneG5IP"};

  mutable Gaudi::Accumulators::StatCounter<> m_nToFind{this, "nToFind"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFound{this, "nFound"};
  mutable Gaudi::Accumulators::StatCounter<> m_nToFindB{this, "nToFindB"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundB{this, "nFoundB"};
  mutable Gaudi::Accumulators::StatCounter<> m_nToFindBAll{this, "nToFindBAll"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundBAll{this, "nFoundBAll"};
  mutable Gaudi::Accumulators::StatCounter<> m_nToFindKsL{this, "nToFindKsL"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundKsL{this, "nFoundKsL"};
  mutable Gaudi::Accumulators::StatCounter<> m_nToFindIP{this, "nToFindIP"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundIP{this, "nFoundIP"};

  mutable Gaudi::Accumulators::StatCounter<> m_nToFindG5{this, "nToFindG5"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundG5{this, "nFoundG5"};
  mutable Gaudi::Accumulators::StatCounter<> m_nToFindG5B{this, "nToFindG5B"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundG5B{this, "nFoundG5B"};
  mutable Gaudi::Accumulators::StatCounter<> m_nToFindG5BAll{this, "nToFindG5BAll"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundG5BAll{this, "nFoundG5BAll"};
  mutable Gaudi::Accumulators::StatCounter<> m_nToFindG5KsL{this, "nToFindG5KsL"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundG5KsL{this, "nFoundG5KsL"};
  mutable Gaudi::Accumulators::StatCounter<> m_nToFindG5IP{this, "nToFindG5IP"};
  mutable Gaudi::Accumulators::StatCounter<> m_nFoundG5IP{this, "nFoundG5IP"};

  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiency{this, "efficiency"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyG5{this, "efficiencyG5"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyB{this, "efficiencyB"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyG5B{this, "efficiencyG5B"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyBAll{this, "efficiencyBAll"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyG5BAll{this, "efficiencyG5BAll"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyKsL{this, "efficiencyKsL"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyG5KsL{this, "efficiencyG5KsL"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyIP{this, "efficiencyIP"};
  mutable Gaudi::Accumulators::SigmaCounter<> m_efficiencyG5IP{this, "efficiencyG5IP"};
};

DECLARE_COMPONENT( TrackEffChecker )

TrackEffChecker::TrackEffChecker( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator,
                {KeyValue{"TracksInContainer", LHCb::TrackLocation::Default},
                 KeyValue{"MCParticleInContainer", LHCb::MCParticleLocation::Default},
                 KeyValue{"LinkerInTable", "Link/" + LHCb::TrackLocation::Default},
                 KeyValue{"AllLinksLocation", "Link/Pat/LHCbID"}} ) {}

void TrackEffChecker::operator()( const LHCb::Tracks& tracks, const LHCb::MCParticles& mcParts,
                                  const LHCb::LinksByKey& links, const LHCb::LinksByKey& allIds ) const {

  std::vector<std::vector<int>> linkedIds;
  allIds.applyToAllLinks( [&linkedIds]( unsigned int srcKey, unsigned int mcPartKey, float ) {
    while ( linkedIds.size() <= mcPartKey ) {
      std::vector<int> dum;
      linkedIds.push_back( dum );
    }
    linkedIds[mcPartKey].push_back( srcKey );
  } );

  unsigned int nTracksInThisEvent =
      std::count_if( tracks.begin(), tracks.end(), [requireLong = m_requireLong.value()]( const auto* track ) {
        return track->type() == LHCb::Track::Types::Long || !requireLong;
      } );

  m_nTrackCounter += nTracksInThisEvent;

  // we want to count ghosts etc
  ghostInfo( tracks, mcParts, links, nTracksInThisEvent );

  // then we want to look at efficiencies
  effInfo( tracks, mcParts, links, linkedIds );
}

void TrackEffChecker::ghostInfo( const LHCb::Tracks& tracks, const LHCb::MCParticles& mcParts,
                                 const LHCb::LinksByKey& links, unsigned int nTracksInThisEvent ) const {

  unsigned int nGhost = 0;
  std::string  type   = "";
  for ( const auto* track : tracks ) {

    if ( track->type() != LHCb::Track::Types::Long && m_requireLong.value() ) continue;

    const LHCb::MCParticle* particle = mcTruth( *track, mcParts, links );

    splitByAlgorithm() == true ? type = Gaudi::Utils::toString( track->history() ) : type = "ALL";
    if ( !particle ) {
      ++nGhost;
      if ( fullDetail() == true ) {
        plots( type + "/ghost", track );
        LHCb::GhostTrackInfo gInfo;
        StatusCode           sc = ghostClassification()->info( *track, gInfo );
        if ( sc.isSuccess() ) {
          const LHCb::GhostTrackInfo::Classification& gtype = gInfo.classification();
          plot( gtype, "ghost classification", -0.5, 30.5, 31 );
          ++counter( Gaudi::Utils::toString( gtype ) );
        }
      }
    } else {
      if ( fullDetail() == true ) plots( type + "/real", track );
    }
  }

  // counter for ghost rate
  m_nGhostCounter += nGhost;

  // plot the event ghost rate
  if ( nTracksInThisEvent != 0 ) plot( nGhost / double( nTracksInThisEvent ), "ghost rate", -0.01, 1.01, 51 );
  if ( fullDetail() == true ) {
    // ghost rate versus # interactions
    long nVert = visPrimVertTool()->countVertices();
    if ( nTracksInThisEvent != 0 )
      plot2D( nVert, 100 * nGhost / double( nTracksInThisEvent ), "ghost rate vs visible", -0.5, 20.5, -0.5, 100.5, 21,
              101 );
  }
}

void TrackEffChecker::effInfo( const LHCb::Tracks& tracks, const LHCb::MCParticles& mcParts,
                               const LHCb::LinksByKey& links, const std::vector<std::vector<int>>& linkedIds ) const {

  double efficiency       = 0;
  double efficiencyG5     = 0;
  double efficiencyB      = 0;
  double efficiencyG5B    = 0;
  double efficiencyBAll   = 0;
  double efficiencyG5BAll = 0;
  double efficiencyKsL    = 0;
  double efficiencyG5KsL  = 0;
  double efficiencyIP     = 0;
  double efficiencyG5IP   = 0;

  unsigned int nToFind   = 0u;
  unsigned int nFound    = 0u;
  unsigned int nToFindG5 = 0u;
  unsigned int nFoundG5  = 0u;

  unsigned int nToFindB   = 0u;
  unsigned int nFoundB    = 0u;
  unsigned int nToFindG5B = 0u;
  unsigned int nFoundG5B  = 0u;

  unsigned int nToFindBAll   = 0u;
  unsigned int nFoundBAll    = 0u;
  unsigned int nToFindG5BAll = 0u;
  unsigned int nFoundG5BAll  = 0u;

  unsigned int nToFindKsL   = 0u;
  unsigned int nFoundKsL    = 0u;
  unsigned int nToFindG5KsL = 0u;
  unsigned int nFoundG5KsL  = 0u;

  unsigned int nToFindIP   = 0u;
  unsigned int nFoundIP    = 0u;
  unsigned int nToFindG5IP = 0u;
  unsigned int nFoundG5IP  = 0u;

  // Build a map of reconstructed tracks for each MCParticle
  std::map<int, TrackCheckerBase::LinkInfo> reconstructedTracks;
  links.applyToAllLinks(
      [&reconstructedTracks, &tracks]( unsigned int trackKey, unsigned int mcPartKey, float weight ) {
        const LHCb::Track* track = static_cast<const LHCb::Track*>( tracks.containedObject( trackKey ) );
        auto [it, inserted]      = reconstructedTracks.try_emplace( mcPartKey, track, 0, weight );
        if ( !inserted ) ++it->second.clone;
      } );

  for ( const auto* mcPart : mcParts ) {
    bool reconstructible     = false;
    bool reconstructibleB    = false;
    bool reconstructibleBAll = false;
    bool reconstructibleKsL  = false;
    bool largeIP             = false;

    reconstructible     = selected( mcPart );
    reconstructibleBAll = reconstructible && bAncestorWithReconstructibleDaughters( mcPart );
    reconstructibleB    = reconstructible && bAncestor( mcPart );
    reconstructibleKsL  = reconstructible && ksLambdaAncestor( mcPart );

    if ( reconstructible == true ) {

      double x  = mcPart->originVertex()->position4vector().x();
      double y  = mcPart->originVertex()->position4vector().y();
      double tx = mcPart->momentum().x() / mcPart->momentum().z();
      double ty = mcPart->momentum().y() / mcPart->momentum().z();

      double IP2 = ( x * x + y * y ) - ( tx * x + ty * y ) * ( tx * x + ty * y ) / ( tx * tx + ty * ty );
      largeIP    = reconstructibleB && IP2 > 1. && IP2 < 4.; // 1 to 2 mm

      double nTrue = 0;

      TrackCheckerBase::LinkInfo info = reconstructedTracks[mcPart->key()];

      if ( info.track ) {

        std::vector<LHCb::LHCbID> ids;
        if ( linkedIds.size() > (unsigned int)mcPart->key() ) {
          for ( const auto id : linkedIds[mcPart->key()] ) {
            LHCb::LHCbID temp;
            temp.setDetectorType( LHCb::LHCbID::channelIDtype( id >> 28 ) ); // create LHCbId from int. Clumsy !
            temp.setID( id );
            ids.push_back( temp );
          }
        }

        for ( const LHCb::LHCbID& id : ids ) {
          if ( id.isVP() ) {
            if ( info.track->type() == LHCb::Track::Types::Ttrack ||
                 info.track->type() == LHCb::Track::Types::Downstream )
              continue;
            nTrue += 1.;
          } else if ( id.isVP() ) {
            if ( info.track->type() == LHCb::Track::Types::Ttrack ||
                 info.track->type() == LHCb::Track::Types::Downstream )
              continue;
            nTrue += 1.;
          } else if ( id.isUT() ) {
            if ( info.track->type() != LHCb::Track::Types::Ttrack && info.track->type() != LHCb::Track::Types::Velo )
              nTrue += 1.;
          }
        }
      }

      plots( "mcSelected", mcPart );
      ++nToFind;
      if ( mcPart->p() > 5 * Gaudi::Units::GeV ) ++nToFindG5;

      if ( reconstructibleB ) {
        ++nToFindB;
        if ( mcPart->p() > 5 * Gaudi::Units::GeV ) ++nToFindG5B;
      }

      if ( reconstructibleBAll ) {
        ++nToFindBAll;
        if ( mcPart->p() > 5 * Gaudi::Units::GeV ) ++nToFindG5BAll;
      }

      if ( reconstructibleKsL ) {
        ++nToFindKsL;
        if ( mcPart->p() > 5 * Gaudi::Units::GeV ) ++nToFindG5KsL;
      }

      if ( largeIP ) {
        ++nToFindIP;
        if ( mcPart->p() > 5 * Gaudi::Units::GeV ) ++nToFindG5IP;
      }

      if ( info.track != 0 && ( info.track->type() == LHCb::Track::Types::Long || !m_requireLong.value() ) ) {
        ++nFound;
        efficiency += info.purity * info.track->lhcbIDs().size() / nTrue;

        if ( mcPart->p() > 5 * Gaudi::Units::GeV ) {
          ++nFoundG5;
          efficiencyG5 += info.purity * info.track->lhcbIDs().size() / nTrue;
          ;
        }

        plots( "selected", info.track );
        plots( "mcSelectedAndRec", mcPart );

        plot( info.purity, "hitpurity", -0.01, 1.01, 51 );
        if ( mcPart->p() > 5 * Gaudi::Units::GeV ) plot( info.purity, "hitpurityG5", -0.01, 1.01, 51 );

        m_nClone += info.clone;
        if ( mcPart->p() > 5 * Gaudi::Units::GeV ) m_nCloneG5 += info.clone;

        if ( reconstructibleB ) {

          if ( info.track != 0 && ( info.track->type() == LHCb::Track::Types::Long || !m_requireLong.value() ) ) {
            ++nFoundB;
            efficiencyB += info.purity * info.track->lhcbIDs().size() / nTrue;
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) {
              ++nFoundG5B;
              efficiencyG5B += info.purity * info.track->lhcbIDs().size() / nTrue;
            }

            plot( info.purity, "hitpurityB", -0.01, 1.01, 51 );
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) plot( info.purity, "hitpurityG5B", -0.01, 1.01, 51 );

            m_nCloneB += info.clone;
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) m_nCloneG5B += info.clone;
          }
        }

        if ( reconstructibleBAll ) {

          if ( info.track != 0 && ( info.track->type() == LHCb::Track::Types::Long || !m_requireLong.value() ) ) {
            ++nFoundBAll;
            efficiencyBAll += info.purity * info.track->lhcbIDs().size() / nTrue;
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) {
              ++nFoundG5BAll;
              efficiencyG5BAll += info.purity * info.track->lhcbIDs().size() / nTrue;
            }

            plot( info.purity, "hitpurityBAll", -0.01, 1.01, 51 );
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) plot( info.purity, "hitpurityG5BAll", -0.01, 1.01, 51 );

            m_nCloneBAll += info.clone;
            if ( mcPart->p() > 5 * Gaudi::Units::GeV )
              m_nCloneG5BAll += info.clone * info.track->lhcbIDs().size() / nTrue;
          }
        }

        if ( reconstructibleKsL ) {

          if ( info.track != 0 && ( info.track->type() == LHCb::Track::Types::Long || !m_requireLong.value() ) ) {
            ++nFoundKsL;
            efficiencyKsL += info.purity * info.track->lhcbIDs().size() / nTrue;

            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) {
              ++nFoundG5KsL;
              efficiencyG5KsL += info.purity * info.track->lhcbIDs().size() / nTrue;
            }

            plot( info.purity, "hitpurityKsL", -0.01, 1.01, 51 );
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) plot( info.purity, "hitpurityG5KsL", -0.01, 1.01, 51 );

            m_nCloneKsL += info.clone;
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) m_nCloneG5KsL += info.clone;
          }
        }

        if ( largeIP ) {

          if ( info.track != 0 && ( info.track->type() == LHCb::Track::Types::Long || !m_requireLong.value() ) ) {
            ++nFoundIP;
            efficiencyIP += info.purity * info.track->lhcbIDs().size() / nTrue;

            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) {
              ++nFoundG5IP;
              efficiencyG5IP += info.purity * info.track->lhcbIDs().size() / nTrue;
            }

            plot( info.purity, "hitpurityIP", -0.01, 1.01, 51 );
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) plot( info.purity, "hitpurityG5IP", -0.01, 1.01, 51 );

            m_nCloneIP += info.clone;
            if ( mcPart->p() > 5 * Gaudi::Units::GeV ) m_nCloneG5IP += info.clone;
          }
        }
      }
    }
  } // loop particles

  // update counters
  m_nToFind += nToFind;
  m_nFound += nFound;
  m_nToFindB += nToFindB;
  m_nFoundB += nFoundB;
  m_nToFindBAll += nToFindBAll;
  m_nFoundBAll += nFoundBAll;
  m_nToFindKsL += nToFindKsL;
  m_nFoundKsL += nFoundKsL;
  m_nToFindIP += nToFindIP;
  m_nFoundIP += nFoundIP;

  m_nToFindG5 += nToFindG5;
  m_nFoundG5 += nFoundG5;
  m_nToFindG5B += nToFindG5B;
  m_nFoundG5B += nFoundG5B;
  m_nToFindG5BAll += nToFindG5BAll;
  m_nFoundG5BAll += nFoundG5BAll;
  m_nToFindG5KsL += nToFindG5KsL;
  m_nFoundG5KsL += nFoundG5KsL;
  m_nToFindG5IP += nToFindG5IP;
  m_nFoundG5IP += nFoundG5IP;

  m_efficiency += efficiency;
  m_efficiencyG5 += efficiencyG5;
  m_efficiencyB += efficiencyB;
  m_efficiencyG5B += efficiencyG5B;
  m_efficiencyBAll += efficiencyBAll;
  m_efficiencyG5BAll += efficiencyG5BAll;
  m_efficiencyKsL += efficiencyKsL;
  m_efficiencyG5KsL += efficiencyG5KsL;
  m_efficiencyIP += efficiencyIP;
  m_efficiencyG5IP += efficiencyG5IP;

  // event efficiency
  if ( nToFind != 0u ) plot( nFound / double( nToFind ), "eff", -0.01, 1.01, 51 );
  if ( nToFindG5 != 0u ) plot( nFoundG5 / double( nToFindG5 ), "effG5", -0.01, 1.01, 51 );
  if ( nToFindB != 0u ) plot( nFoundB / double( nToFindB ), "effB", -0.01, 1.01, 51 );
  if ( nToFindG5B != 0u ) plot( nFoundG5B / double( nToFindG5B ), "effG5B", -0.01, 1.01, 51 );
  if ( nToFindBAll != 0u ) plot( nFoundBAll / double( nToFindBAll ), "effBAll", -0.01, 1.01, 51 );
  if ( nToFindG5BAll != 0u ) plot( nFoundG5BAll / double( nToFindG5BAll ), "effG5BAll", -0.01, 1.01, 51 );
  if ( nToFindKsL != 0u ) plot( nFoundKsL / double( nToFindKsL ), "effKsL", -0.01, 1.01, 51 );
  if ( nToFindG5KsL != 0u ) plot( nFoundG5KsL / double( nToFindG5KsL ), "effG5KsL", -0.01, 1.01, 51 );

  if ( fullDetail() == true ) {
    // ghost rate versus # interactions
    long nVert = visPrimVertTool()->countVertices();
    if ( nToFind != 0u )
      plot2D( nVert, 100 * nFound / double( nToFind ), "eff vs visible", -0.5, 20.5, -0.5, 100.5, 21, 101 );
  }
}

void TrackEffChecker::plots( const std::string& type, const LHCb::MCParticle* part ) const {

  plot( part->p() / Gaudi::Units::GeV, type + "/p", "p", 0., 50., 25 );
  plot( part->pt() / Gaudi::Units::GeV, type + "/pt", "pt", 0., 10., 20 );
  plot( part->pseudoRapidity(), type + "/eta", "eta", 1., 6., 25 );
}

void TrackEffChecker::plots( const std::string& type, const LHCb::Track* track ) const {

  //  const double nMeas = weightedMeasurementSum(track);

  if ( track->type() != LHCb::Track::Types::Velo && track->history() != LHCb::Track::History::PrPixel ) {
    plot( track->pt() / Gaudi::Units::GeV, type + "/pt", "pt", 0., 10., 100 );
    plot( track->p() / Gaudi::Units::GeV, type + "/p", "p", 0., 50., 25 );
  }
  plot( track->chi2PerDoF(), type + "/chi2", "chi2", 0., 500., 100 );
  plot( track->probChi2(), type + "/probchi2", "probChi2", 0., 1., 100 );
  plot( track->pseudoRapidity(), type + "/eta", "eta", 1., 6., 25 );

  // information from the extra info list
  const LHCb::Track::ExtraInfo&          info     = track->extraInfo();
  LHCb::Track::ExtraInfo::const_iterator iterInfo = info.begin();
  for ( ; iterInfo != info.end(); ++iterInfo ) {
    LHCb::Track::AdditionalInfo            infoName = LHCb::Track::AdditionalInfo( iterInfo->first );
    std::string                            title    = Gaudi::Utils::toString( infoName );
    const TrackMaps::InfoHistMap&          histMap  = TrackMaps::infoHistDescription();
    TrackMaps::InfoHistMap::const_iterator iterM    = histMap.find( infoName );
    if ( iterM != histMap.end() ) {
      const TrackMaps::HistoRange range = histMap.find( infoName )->second;
      plot( iterInfo->second, type + "/info/" + range.fid, title, range.fxMin, range.fxMax, 100 );
    }
  } // iterInfo
}

double TrackEffChecker::weightedMeasurementSum( const LHCb::Track* aTrack ) const {

  double                           wSum = 0;
  const std::vector<LHCb::LHCbID>& ids  = aTrack->lhcbIDs();
  for ( std::vector<LHCb::LHCbID>::const_iterator iter = ids.begin(); iter != ids.end(); ++iter ) {
    if ( iter->isVP() == true ) {
      wSum += 1;
    } else if ( iter->isFT() == true ) {
      wSum += 1;
    } else {
      // nothing
    }
  }
  return wSum;
}

StatusCode TrackEffChecker::finalize() {

  // ghost rate
  std::string         histName = "ghost rate";
  AIDA::IHistogram1D* hist     = histo1D( histName );
  double              eGhost   = 0;
  if ( hist != 0 ) eGhost = hist->mean();
  const double tGhost =
      m_nTrackCounter.sum() == 0 ? 0.0 : double( m_nGhostCounter.sum() ) / double( m_nTrackCounter.sum() );

  histName   = "hitpurity";
  hist       = histo1D( histName );
  double pur = 0;
  if ( hist != 0 ) pur = hist->mean();

  histName     = "hitpurityG5";
  hist         = histo1D( histName );
  double purG5 = 0;
  if ( hist != 0 ) purG5 = hist->mean();

  histName    = "hitpurityB";
  hist        = histo1D( histName );
  double purB = 0;
  if ( hist != 0 ) purB = hist->mean();

  histName      = "hitpurityG5B";
  hist          = histo1D( histName );
  double purG5B = 0;
  if ( hist != 0 ) purG5B = hist->mean();

  histName       = "hitpurityBAll";
  hist           = histo1D( histName );
  double purBAll = 0;
  if ( hist != 0 ) purBAll = hist->mean();

  histName         = "hitpurityG5BAll";
  hist             = histo1D( histName );
  double purG5BAll = 0;
  if ( hist != 0 ) purG5BAll = hist->mean();

  histName      = "hitpurityKsL";
  hist          = histo1D( histName );
  double purKsL = 0;
  if ( hist != 0 ) purKsL = hist->mean();

  histName        = "hitpurityG5KsL";
  hist            = histo1D( histName );
  double purG5KsL = 0;
  if ( hist != 0 ) purG5KsL = hist->mean();

  histName     = "hitpurityIP";
  hist         = histo1D( histName );
  double purIP = 0;
  if ( hist != 0 ) purIP = hist->mean();

  histName       = "hitpurityG5IP";
  hist           = histo1D( histName );
  double purG5IP = 0;
  if ( hist != 0 ) purG5IP = hist->mean();

  const double tEff = m_nToFind.nEntries() == 0 ? 0.0 : double( m_nFound.nEntries() ) / double( m_nToFind.nEntries() );

  const double tEffG5 =
      m_nToFindG5.nEntries() == 0 ? 0.0 : double( m_nFoundG5.nEntries() ) / double( m_nToFindG5.nEntries() );

  const double hEff =
      m_nFound.nEntries() == 0 ? 0.0 : double( m_efficiency.nEntries() ) / double( m_nFound.nEntries() );

  const double hEffG5 =
      m_nFoundG5.nEntries() == 0 ? 0.0 : double( m_efficiencyG5.nEntries() ) / double( m_nFoundG5.nEntries() );

  const double tCloneG5 =
      m_nToFindG5.nEntries() == 0 ? 0.0 : double( m_nCloneG5.nEntries() ) / double( m_nToFindG5.nEntries() );

  const double tClone =
      m_nToFind.nEntries() == 0 ? 0.0 : double( m_nClone.nEntries() ) / double( m_nToFind.nEntries() );

  const double tEffB =
      m_nToFindB.nEntries() == 0 ? 0.0 : double( m_nFoundB.nEntries() ) / double( m_nToFindB.nEntries() );

  const double tEffG5B =
      m_nToFindG5B.nEntries() == 0 ? 0.0 : double( m_nFoundG5B.nEntries() ) / double( m_nToFindG5B.nEntries() );

  const double hEffB =
      m_nFoundB.nEntries() == 0 ? 0.0 : double( m_efficiencyB.nEntries() ) / double( m_nFoundB.nEntries() );

  const double hEffG5B =
      m_nFoundG5B.nEntries() == 0 ? 0.0 : double( m_efficiencyG5B.nEntries() ) / double( m_nFoundG5B.nEntries() );

  const double tCloneG5B =
      m_nToFindG5B.nEntries() == 0 ? 0.0 : double( m_nCloneG5B.nEntries() ) / double( m_nToFindG5B.nEntries() );

  const double tCloneB =
      m_nToFindB.nEntries() == 0 ? 0.0 : double( m_nCloneB.nEntries() ) / double( m_nToFindB.nEntries() );

  const double tEffBAll =
      m_nToFindBAll.nEntries() == 0 ? 0.0 : double( m_nFoundBAll.nEntries() ) / double( m_nToFindBAll.nEntries() );

  const double tEffG5BAll = m_nToFindG5BAll.nEntries() == 0
                                ? 0.0
                                : double( m_nFoundG5BAll.nEntries() ) / double( m_nToFindG5BAll.nEntries() );

  const double hEffBAll =
      m_nFoundBAll.nEntries() == 0 ? 0.0 : double( m_efficiencyBAll.nEntries() ) / double( m_nFoundBAll.nEntries() );

  const double hEffG5BAll = m_nFoundG5BAll.nEntries() == 0
                                ? 0.0
                                : double( m_efficiencyG5BAll.nEntries() ) / double( m_nFoundG5BAll.nEntries() );

  const double tCloneG5BAll = m_nToFindG5BAll.nEntries() == 0
                                  ? 0.0
                                  : double( m_nCloneG5BAll.nEntries() ) / double( m_nToFindG5BAll.nEntries() );

  const double tCloneBAll =
      m_nToFindBAll.nEntries() == 0 ? 0.0 : double( m_nCloneBAll.nEntries() ) / double( m_nToFindBAll.nEntries() );

  const double tEffKsL =
      m_nToFindKsL.nEntries() == 0 ? 0.0 : double( m_nFoundKsL.nEntries() ) / double( m_nToFindKsL.nEntries() );

  const double tEffG5KsL =
      m_nToFindG5KsL.nEntries() == 0 ? 0.0 : double( m_nFoundG5KsL.nEntries() ) / double( m_nToFindG5KsL.nEntries() );

  const double hEffKsL =
      m_nFoundKsL.nEntries() == 0 ? 0.0 : double( m_efficiencyKsL.nEntries() ) / double( m_nFoundKsL.nEntries() );

  const double hEffG5KsL =
      m_nFoundG5KsL.nEntries() == 0 ? 0.0 : double( m_efficiencyG5KsL.nEntries() ) / double( m_nFoundG5KsL.nEntries() );

  const double tCloneG5KsL =
      m_nToFindG5KsL.nEntries() == 0 ? 0.0 : double( m_nCloneG5KsL.nEntries() ) / double( m_nToFindG5KsL.nEntries() );

  const double tCloneKsL =
      m_nToFindKsL.nEntries() == 0 ? 0.0 : double( m_nCloneKsL.nEntries() ) / double( m_nToFindKsL.nEntries() );

  const double tEffIP =
      m_nToFindIP.nEntries() == 0 ? 0.0 : double( m_nFoundIP.nEntries() ) / double( m_nToFindIP.nEntries() );

  const double tEffG5IP =
      m_nToFindG5IP.nEntries() == 0 ? 0.0 : double( m_nFoundG5IP.nEntries() ) / double( m_nToFindG5IP.nEntries() );

  const double hEffIP =
      m_nFoundIP.nEntries() == 0 ? 0.0 : double( m_efficiencyIP.nEntries() ) / double( m_nFoundIP.nEntries() );

  const double hEffG5IP =
      m_nFoundG5IP.nEntries() == 0 ? 0.0 : double( m_efficiencyG5IP.nEntries() ) / double( m_nFoundG5IP.nEntries() );

  const double tCloneG5IP =
      m_nToFindG5IP.nEntries() == 0 ? 0.0 : double( m_nCloneG5IP.nEntries() ) / double( m_nToFindG5IP.nEntries() );

  const double tCloneIP =
      m_nToFindIP.nEntries() == 0 ? 0.0 : double( m_nCloneIP.nEntries() ) / double( m_nToFindIP.nEntries() );

  info() << "             **************            "
         << format( "       %8.0f tracks including %8.0f ghosts [%5.1f %%] Event average %5.1f  ****",
                    double( m_nTrackCounter.sum() ), double( m_nGhostCounter.sum() ), tGhost * 100, eGhost * 100 )
         << endmsg;

  info() << "     all  long                         "
         << format( " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                    double( m_nFound.nEntries() ), double( m_nToFind.nEntries() ), tEff * 100, tClone * 100, pur * 100,
                    hEff * 100 )
         << endmsg;

  info() << "     long, p > 5 GeV                   "
         << format( " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiecy",
                    double( m_nFoundG5.nEntries() ), double( m_nToFindG5.nEntries() ), tEffG5 * 100, tCloneG5 * 100,
                    purG5 * 100, hEffG5 * 100 )
         << endmsg;

  if ( m_nToFindB.nEntries() > 0 ) {
    info() << "     all  long B daughers              "
           << format(
                  " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                  double( m_nFoundB.nEntries() ), double( m_nToFindB.nEntries() ), tEffB * 100, tCloneB * 100,
                  purB * 100, hEffB * 100 )
           << endmsg;

    info() << "   long B daughters, p > 5 GeV         "
           << format(
                  " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                  double( m_nFoundG5B.nEntries() ), double( m_nToFindG5B.nEntries() ), tEffG5B * 100, tCloneG5B * 100,
                  purG5B * 100, hEffG5B * 100 )
           << endmsg;
  }

  if ( m_nToFindBAll.nEntries() > 0 ) {
    info() << "     all  long good B daughers         "
           << format(
                  " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                  double( m_nFoundBAll.nEntries() ), double( m_nToFindBAll.nEntries() ), tEffBAll * 100,
                  tCloneBAll * 100, purBAll * 100, hEffBAll * 100 )
           << endmsg;

    info() << "   long good B daughters , p > 5 GeV   "
           << format(
                  " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                  double( m_nFoundG5BAll.nEntries() ), double( m_nToFindG5BAll.nEntries() ), tEffG5BAll * 100,
                  tCloneG5BAll * 100, purG5BAll * 100, hEffG5BAll * 100 )
           << endmsg;
  }

  if ( m_nToFindKsL.nEntries() > 0 ) {
    info() << "     all  long  Ks/Lambda              "
           << format(
                  " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                  double( m_nFoundKsL.nEntries() ), double( m_nToFindKsL.nEntries() ), tEffKsL * 100, tCloneKsL * 100,
                  purKsL * 100, hEffKsL * 100 )
           << endmsg;

    info() << "   long Ks/Lambda, p > 5 GeV           "
           << format(
                  " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                  double( m_nFoundG5KsL.nEntries() ), double( m_nToFindG5KsL.nEntries() ), tEffG5KsL * 100,
                  tCloneG5KsL * 100, purG5KsL * 100, hEffG5KsL * 100 )
           << endmsg;
  }

  if ( m_nToFindIP.nEntries() > 0 ) {
    info() << "     all  long large IP                "
           << format(
                  " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                  double( m_nFoundIP.nEntries() ), double( m_nToFindIP.nEntries() ), tEffIP * 100, tCloneIP * 100,
                  purIP * 100, hEffIP * 100 )
           << endmsg;

    info() << "   long largeIP, p > 5 GeV             "
           << format(
                  " :     %8.0f from %8.0f [%5.1f %%]; %5.1f %% clones; %5.1f %% hit purity; %5.1f %% hit efficiency",
                  double( m_nFoundG5IP.nEntries() ), double( m_nToFindG5IP.nEntries() ), tEffG5IP * 100,
                  tCloneG5IP * 100, purG5IP * 100, hEffG5IP * 100 )
           << endmsg;
  }

  return TrackCheckerBase::finalize();
}
