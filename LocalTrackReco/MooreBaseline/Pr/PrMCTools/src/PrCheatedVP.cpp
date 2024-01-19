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
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Event/VPCluster.h"
#include "LHCbAlgs/Transformer.h"
#include "Linker/LinkedFrom.h"
#include "PrFitParams/IPrFitTool.h"

using namespace LHCb;

/** @class PrCheatedVP PrCheatedVP.h
 *  Cheated pattern recognition for the upgraded VELO
 *
 *  @author Olivier Callot
 *  @date   2012-07-26
 */

namespace {
  LHCb::Tracks getTracks( const Gaudi::Algorithm& alg, const MCTrackInfo& trackInfo, const IPrFitTool& fitTool,
                          const LHCb::MCParticles& particles, LHCb::LinksByKey const& my_links,
                          const std::function<void( LHCb::Track&, std::vector<Gaudi::XYZPoint>&,
                                                    const LHCb::MCParticle*, LinkedFrom<LHCb::VPCluster> )>
                              getPoints ) {
    LHCb::Tracks tracks;

    LinkedFrom<LHCb::VPCluster> link{&my_links};

    constexpr double zVelo = 0.;

    for ( const LHCb::MCParticle* const particle : particles ) {
      // Skip particles without track info.
      if ( 0 == trackInfo.fullInfo( particle ) ) continue;
      // Skip particles not linked to a VELO track.
      if ( !trackInfo.hasVelo( particle ) ) continue;
      // Skip electrons.
      if ( abs( particle->particleID().pid() ) == 11 ) continue;

      auto                         track = std::make_unique<LHCb::Track>();
      std::vector<Gaudi::XYZPoint> points;
      getPoints( *track, points, particle, link );

      // Make a straight-line fit of the track.
      const auto xResult = fitTool.fitLine( points, IPrFitTool::XY::X, zVelo );
      const auto yResult = fitTool.fitLine( points, IPrFitTool::XY::Y, zVelo );
      if ( !xResult || !yResult ) {
        alg.err() << "Fit matrix is singular" << endmsg;
        continue;
      }

      const auto& [x0, x1] = *xResult;
      const auto& [y0, y1] = *yResult;

      LHCb::State state;
      state.setLocation( LHCb::State::Location::ClosestToBeam );
      state.setState( x0, y0, zVelo, x1, y1, 0. );
      track->addToStates( state );
      if ( 0 > particle->momentum().z() ) {
        track->setType( LHCb::Track::Types::VeloBackward );
        // Cut out backwards tracks.
        // continue;
      } else {
        track->setType( LHCb::Track::Types::Velo );
      }
      tracks.insert( track.release() );
    }

    return tracks;
  }
} // namespace

template <bool useMCHits>
class PrCheatedVPBase
    : public LHCb::Algorithm::Transformer<std::conditional_t<
          useMCHits,
          LHCb::Tracks( const LHCb::MCParticles&, const LHCb::MCHits&, const LHCb::MCProperty&,
                        const LHCb::LinksByKey& ),
          LHCb::Tracks( const LHCb::MCParticles&, const LHCb::MCProperty&, const LHCb::LinksByKey& )>> {
public:
  /// Using Transfomer's constructor
  using PrCheatedVPBase::Transformer::Transformer;

protected:
  ToolHandle<const IPrFitTool> m_fitTool{"PrFitTool", this};
};

struct PrCheatedVP final : PrCheatedVPBase<false> {

  PrCheatedVP( const std::string& name, ISvcLocator* pSvcLocator )
      : PrCheatedVPBase( name, pSvcLocator,
                         {KeyValue{"MCParticles", LHCb::MCParticleLocation::Default},
                          KeyValue{"MCTrackInfo", LHCb::MCPropertyLocation::TrackInfo},
                          KeyValue{"MCLinks", LHCb::LinksByKey::linkerName( LHCb::VPClusterLocation::Default )}},
                         KeyValue{"Tracks", LHCb::TrackLocation::Velo} ) {}

  LHCb::Tracks operator()( const LHCb::MCParticles& particles, const LHCb::MCProperty& trackInfo,
                           const LHCb::LinksByKey& my_links ) const override {
    return getTracks( *this, MCTrackInfo{trackInfo}, *m_fitTool, particles, my_links,
                      []( LHCb::Track& track, std::vector<Gaudi::XYZPoint>& points,
                          const LHCb::MCParticle* const particle, LinkedFrom<LHCb::VPCluster> link ) {
                        for ( const LHCb::VPCluster& cluster : link.range( particle ) ) {
                          track.addToLhcbIDs( LHCb::LHCbID( cluster.channelID() ) );
                          points.emplace_back( cluster.x(), cluster.y(), cluster.z() );
                        }
                      } );
  }
};
DECLARE_COMPONENT( PrCheatedVP )

struct PrCheatedVPMCHits final : PrCheatedVPBase<true> {
  PrCheatedVPMCHits( const std::string& name, ISvcLocator* pSvcLocator )
      : PrCheatedVPBase( name, pSvcLocator,
                         {KeyValue{"MCParticles", LHCb::MCParticleLocation::Default},
                          KeyValue{"MCHits", LHCb::MCHitLocation::VP},
                          KeyValue{"MCTrackInfo", LHCb::MCPropertyLocation::TrackInfo},
                          KeyValue{"MCLinks", LHCb::LinksByKey::linkerName( LHCb::VPClusterLocation::Default )}},
                         KeyValue{"Tracks", LHCb::TrackLocation::Velo} ) {}

  LHCb::Tracks operator()( const LHCb::MCParticles& particles, const LHCb::MCHits& hits,
                           const LHCb::MCProperty& trackInfo, const LHCb::LinksByKey& my_links ) const override {
    return getTracks( *this, MCTrackInfo{trackInfo}, *m_fitTool, particles, my_links,
                      [&hits]( LHCb::Track& track, std::vector<Gaudi::XYZPoint>& points,
                               const LHCb::MCParticle* const particle, LinkedFrom<LHCb::VPCluster> link ) {
                        for ( const auto id : link.keyRange( particle ) ) {
                          track.addToLhcbIDs( LHCb::LHCbID( LHCb::Detector::VPChannelID( id ) ) );
                        }
                        for ( const LHCb::MCHit* const hit : hits ) {
                          if ( hit->mcParticle() == particle ) { points.emplace_back( hit->midPoint() ); }
                        }
                      } );
  }
};
DECLARE_COMPONENT( PrCheatedVPMCHits )
