/*****************************************************************************\
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
// Include files

#include "DetDesc/GenericConditionAccessorHolder.h"

#include "GaudiKernel/DataObjectHandle.h"
#include "GaudiKernel/IUpdateManagerSvc.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "Kernel/IBIntegrator.h"
#include "LHCbAlgs/Transformer.h"

#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Event/Track_v1.h"

#include "Event/PrHits.h"
#include "LHCbMath/Utils.h"
#include "Magnet/DeMagnet.h"
#include "MuonDet/DeMuonDetector.h"
#include "MuonDet/MuonNamespace.h"
#include "vdt/vdtMath.h"

#include "StandaloneMuonTrack.h"

#include <algorithm>
#include <array>
#include <cmath>

//-----------------------------------------------------------------------------
// Implementation file for class : StandaloneMuonRec
//
// 2004-10-06 : Alessia Satta
//
//  Removed from Hlt/HltMuon and ported to Tr/TrackTools
//
// 2011-03-03 : Paul Seyfert
/*
// Input: MuonHitsLocation   - Output: V1 Muon Tracks, with momentum estimate
// muonSearch: M5 hits as seed, to search for hits in M4, M3, M2 with LineExtrapolation
// secondLoop: search from M5 to M3, M2 if M4 is not found, or search from M4 to M3, M2
// findCoincidence: find the hits within search window
// recMomentum: estimate the momentum of muon track
// deleteClone: remove clones
*/
//-----------------------------------------------------------------------------

namespace {
  // -- initialize the pad size. Hardwired to speed up.
  constexpr std::array<float, 12> m_Xmax = {{                          //   R1  R2   R3   R4
                                             100., 200., 300., 400.,   // M2
                                             100., 200., 300., 400.,   // M3
                                             400., 400., 400., 400.}}; // M4

  constexpr std::array<float, 12> m_Ymax = {{                         //  R1   R2   R3   R4
                                             60., 120., 180., 240.,   // M2
                                             60., 120., 240., 480.,   // M3
                                             60., 120., 240., 480.}}; // M4

  // -- Set tolerances for hit search in region
  constexpr std::array<float, 4> m_tolForRegion{{2.0, 4.0, 8.0, 10.0}};

  class Cache {
  public:
    std::array<float, 4> stationZ{};
    Gaudi::XYZVector     bdl;
    double               zCenter{};
    Cache(){};
    Cache( DeMuonDetector const& det ) {
      for ( int s = 0; s != det.stations(); ++s ) { stationZ[s] = det.getStationZ( s ); }
    }
  };
} // namespace

class StandaloneMuonRec : public LHCb::Algorithm::Transformer<LHCb::Tracks( const MuonHitContainer&, const Cache& ),
                                                              LHCb::DetDesc::usesConditions<Cache>> {

public:
  using base_class_t = LHCb::Algorithm::Transformer<LHCb::Tracks( const MuonHitContainer&, const Cache& ),
                                                    LHCb::DetDesc::usesConditions<Cache>>;
  using base_class_t::addConditionDerivation;

  /// Standard constructor
  StandaloneMuonRec( const std::string& name, ISvcLocator* pSvcLocator )
      : base_class_t( name, pSvcLocator,
                      {KeyValue{"MuonHitsLocation", MuonHitContainerLocation::Default},
                       KeyValue{"ConditionsCache", "StandaloneMuonAlg-" + name + "-ConditionsCache"}},
                      KeyValue{"OutputMuonTracks", "Rec/Track/Muon"} ) {}

  StatusCode initialize() override {
    return base_class_t::initialize().andThen( [&] {
      this->addConditionDerivation( {DeMuonLocation::Default, LHCb::Det::Magnet::det_path},
                                    this->inputLocation<Cache>(),
                                    [&]( const DeMuonDetector& det, const DeMagnet& magnet ) {
                                      Cache                 cache{det};
                                      const Gaudi::XYZPoint begin( 0., 0., 0. );
                                      const Gaudi::XYZPoint end( 0., 0., cache.stationZ[M2] );
                                      m_bIntegrator->calculateBdlAndCenter( magnet.fieldGrid(), begin, end, 0.0001, 0.,
                                                                            cache.zCenter, cache.bdl );
                                      debug()
                                          << "Integrated B field is " << cache.bdl.x() << " Tm"
                                          << "  centered at Z=" << cache.zCenter / Gaudi::Units::m << " m" << endmsg;
                                      return cache;
                                    } );
    } );
  }

  LHCb::Tracks operator()( const MuonHitContainer& hitContainer, const Cache& cache ) const override;

  ToolHandle<IBIntegrator> m_bIntegrator{this, "BIntegrator", "BIntegrator"};

private:
  enum { M2 = 0, M3, M4, M5 };
  Gaudi::Property<bool>               m_cloneKiller{this, "CloneKiller", true};
  Gaudi::Property<bool>               m_chi2Cut{this, "Chi2Cut", false};
  Gaudi::Property<float>              m_maxchi2Cut{this, "MaxChi2Cut", 1.0};
  Gaudi::Property<bool>               m_strongCloneKiller{this, "StrongCloneKiller", false};
  Gaudi::Property<bool>               m_secondLoop{this, "SecondLoop", false};
  Gaudi::Property<std::vector<float>> m_ParabolicCorrection{this, "ParabolicCorrection", {1.04, 0.14}};
  Gaudi::Property<std::vector<float>> m_resParams{this, "m_resParams", {0.015, 0.29}};
  Gaudi::Property<float>              m_Constant{this, "ConstantCorrection", 0., "In MeV"};

  std::vector<StandaloneMuonTrack> muonSearch( const MuonHitContainer& hitContainer, const Cache& cache ) const;
  bool findCoincidence( const float x, const float y, const unsigned int station, const unsigned int regionBefore,
                        const CommonMuonHitRange& hits, CommonMuonHit& hitcand ) const;
  void findmuonTrack( const MuonHitContainer& hitContainer, const Cache& cache,
                      std::array<CommonMuonHit, 4>& bestCandidates, const int seed,
                      std::vector<StandaloneMuonTrack>& muonTracks ) const;
  void recMomentum( StandaloneMuonTrack& muonTrack, const Cache& cache, LHCb::Track& track ) const;
  void detectClone( std::vector<StandaloneMuonTrack>& muonTracks, const Cache& cache ) const;

  // counters
  mutable Gaudi::Accumulators::Counter<>                    m_countEvents{this, "nEvents"};
  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_countMuCandidates{this, "nMuonTrackCandidates"};

  mutable Gaudi::Accumulators::MsgCounter<MSG::INFO> m_failed_linearfit{this, "Linear Fit Failed!"};
  mutable Gaudi::Accumulators::MsgCounter<MSG::INFO> m_error_zeroBint{this, "B integral is 0!!"};
};

DECLARE_COMPONENT( StandaloneMuonRec )

//=============================================================================
// Main execution
//=============================================================================
LHCb::Tracks StandaloneMuonRec::operator()( const MuonHitContainer& hitContainer, const Cache& cache ) const {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  LHCb::Tracks outputTracks;
  outputTracks.reserve( 1024 );

  ++m_countEvents;

  auto muonTracks = muonSearch( hitContainer, cache );
  if ( m_cloneKiller ) { detectClone( muonTracks, cache ); }

  for ( auto& muTrack : muonTracks ) {
    if ( muTrack.isClone() ) continue;

    auto sc = muTrack.linearFit();
    if ( !sc ) {
      ++m_failed_linearfit;
      continue;
    }

    auto trgTr = std::make_unique<LHCb::Track>();
    recMomentum( muTrack, cache, *trgTr );
    if ( m_chi2Cut && trgTr->chi2PerDoF() > m_maxchi2Cut.value() ) continue;
    outputTracks.insert( trgTr.release() );
  }
  m_countMuCandidates += outputTracks.size();

  if ( msgLevel( MSG::DEBUG ) ) debug() << " stored candidates " << outputTracks.size() << endmsg;

  return outputTracks;
}

bool StandaloneMuonRec::findCoincidence( const float x, const float y, const unsigned int station,
                                         const unsigned int regionBefore, const CommonMuonHitRange& hits,
                                         CommonMuonHit& hitcand ) const {

  const auto tol       = m_tolForRegion[regionBefore];
  float      deltaYmin = 9999.;
  float      deltaXmin = 9999.;
  bool       findCand  = false;
  for ( auto& hit : hits ) {
    float deltaX = fabs( x - hit.x() );
    float deltaY = fabs( y - hit.y() );
    //-- Check if the hit is within the FOI (copy from MuonCombRec); in any case a xFOI >1000 is not considered
    if ( deltaX < m_Xmax[station * 4 + regionBefore] && deltaY < m_Ymax[station * 4 + regionBefore] &&
         ( deltaY < deltaYmin - tol ||
           ( deltaY < deltaYmin + tol && ( deltaX < deltaXmin - tol || fabs( deltaXmin - deltaX ) < 0.1 ) ) ) ) {
      deltaXmin = deltaX;
      deltaYmin = deltaY;
      hitcand   = hit;
      findCand  = true;
    }
  }
  return findCand;
}
void StandaloneMuonRec::findmuonTrack( const MuonHitContainer& hitContainer, const Cache& cache,
                                       std::array<CommonMuonHit, 4>& bestCandidates, const int seed,
                                       std::vector<StandaloneMuonTrack>& muonTracks ) const {
  float        xseed      = bestCandidates[seed].x();
  float        yseed      = bestCandidates[seed].y();
  unsigned int hit_region = bestCandidates[seed].region();
  float        x          = xseed * cache.stationZ[seed - 1] / cache.stationZ[seed];
  float        y          = yseed * cache.stationZ[seed - 1] / cache.stationZ[seed];

  bool findAll = false;
  for ( int ista = seed - 1; ista > -1; ista-- ) {
    auto findCandidate = findCoincidence( x, y, ista, hit_region, hitContainer.hits( ista ), bestCandidates[ista] );
    if ( !findCandidate && ( !m_secondLoop || ista < M4 ) ) break;
    if ( findCandidate && ista == M2 ) {
      findAll = true;
      break;
    }

    // if no hit found in M4 but m_secondLoop enable
    if ( m_secondLoop && ista == M4 && !findCandidate ) {
      bestCandidates[ista] = bestCandidates[seed];
      x                    = bestCandidates[ista + 1].x() * cache.stationZ[ista - 1] / cache.stationZ[ista + 1];
      y                    = bestCandidates[ista + 1].y() * cache.stationZ[ista - 1] / cache.stationZ[ista + 1];
    } else {
      x = -1.0 * ( bestCandidates[ista + 1].x() - bestCandidates[ista].x() ) + bestCandidates[ista].x();
      y = bestCandidates[ista].y() * cache.stationZ[ista - 1] / cache.stationZ[ista];
      if ( fabs( cache.bdl.x() ) < 0.1 && ista == M3 ) {
        x = bestCandidates[ista].x() * cache.stationZ[ista - 1] / cache.stationZ[ista];
      }
    }
    hit_region = bestCandidates[ista].region();
  }
  if ( findAll ) {
    // create the muon track
    StandaloneMuonTrack muon;
    muon.setPoint( 0, bestCandidates[M2] );
    muon.setPoint( 1, bestCandidates[M3] );
    if ( m_secondLoop && bestCandidates[M4].station() == bestCandidates[M5].station() ) {
      muon.setPoint( 2, bestCandidates[M4] );
      muon.setnHits( 3 );
    } else {
      muon.setPoint( 2, bestCandidates[M4] );
      muon.setPoint( 3, bestCandidates[M5] );
      muon.setnHits( 4 );
    }
    muonTracks.push_back( muon );
  }
}
std::vector<StandaloneMuonTrack> StandaloneMuonRec::muonSearch( const MuonHitContainer& hitContainer,
                                                                const Cache&            cache ) const {
  std::vector<StandaloneMuonTrack> muonTracks;
  muonTracks.reserve( 48 );

  const auto&                  hitsM5 = hitContainer.hits( M5 );
  std::array<CommonMuonHit, 4> bestCandidates;
  for ( auto& hit : hitsM5 ) {
    bestCandidates[3] = hit;
    findmuonTrack( hitContainer, cache, bestCandidates, M5, muonTracks );
  }
  ///---second loop from M4
  if ( m_secondLoop ) {
    const auto& hitsM4 = hitContainer.hits( M4 );
    for ( auto& hit : hitsM4 ) {
      // To remove these hits of M4 used in the first round of search
      auto used = std::find_if( muonTracks.begin(), muonTracks.end(),
                                [&hit]( auto mutrack ) { return hit.tile() == mutrack.point( M4 ).tile(); } );
      if ( used != muonTracks.end() ) continue;
      bestCandidates[3] = hit;
      bestCandidates[2] = hit;
      findmuonTrack( hitContainer, cache, bestCandidates, M4, muonTracks );
    }
  }
  return muonTracks;
}

// estimate the momentum of muonTrack
void StandaloneMuonRec::recMomentum( StandaloneMuonTrack& track, const Cache& cache, LHCb::Track& lbtrack ) const {

  const float bdlX          = cache.bdl.x();
  const float FieldPolarity = ( bdlX > 0.0 ? 1 : -1 );

  // create a state at the Z of M2
  const auto      Zfirst = cache.stationZ[M2];
  Gaudi::XYZPoint trackPos( track.bx() + track.sx() * Zfirst, track.by() + track.sy() * Zfirst, Zfirst );
  LHCb::State     state( LHCb::StateVector( trackPos, Gaudi::XYZVector( track.sx(), track.sy(), 1.0 ), 1. / 10000. ) );

  const auto      Zend = cache.stationZ[M5];
  Gaudi::XYZPoint endtrackPos( track.bx() + track.sx() * Zfirst, track.by() + track.sy() * Zend, Zend );
  LHCb::State     endstate(
      LHCb::StateVector( endtrackPos, Gaudi::XYZVector( track.sx(), track.sy(), 1.0 ), 1. / 10000. ) );

  // copied from the MuonTrackMomRec
  // double q = 0.;
  // double p = 1e6 * Gaudi::Units::MeV;

  // can't estimate momentum or charge
  if ( fabs( bdlX ) < TrackParameters::hiTolerance ) { ++m_error_zeroBint; }

  // Rotate to the 0-0-z zixs and do the ptkick
  const auto tX      = state.tx();
  const auto xCenter = state.x() + tX * ( cache.zCenter - state.z() );

  const auto zeta_trk = -tX / sqrt( 1.0 + tX * tX );
  const auto tx_vtx   = xCenter / cache.zCenter;
  const auto zeta_vtx = -tx_vtx / sqrt( 1.0 + tx_vtx * tx_vtx );

  // curvature
  const double curv = ( zeta_trk - zeta_vtx );

  // charge
  int sign = 1;
  if ( curv < TrackParameters::hiTolerance ) { sign *= -1; }
  if ( bdlX < TrackParameters::hiTolerance ) { sign *= -1; }
  const auto q = -1. * FieldPolarity * sign;
  // momentum
  const auto p = Gaudi::Units::eplus * Gaudi::Units::c_light * fabs( bdlX ) *
                 sqrt( ( 1.0 + tX * tX + std::pow( state.ty(), 2 ) ) / ( 1.0 + tX * tX ) ) / fabs( curv );

  /// from Run 2 tunning, commented out for the moment
  /*
  if ( m_ParabolicCorrection.size() == 2u ) {
    // p*= (a + b*tx*tx )
    p += m_Constant;
    p *= ( m_ParabolicCorrection[0] + ( m_ParabolicCorrection[1] * tX * tX ) );
  }
  */

  const double qOverP       = q / p;
  const double err2         = std::pow( m_resParams[0], 2 ) + std::pow( m_resParams[1] / p, 2 );
  double       sigmaQOverP2 = err2 / std::pow( p, 2 );

  state.setQOverP( qOverP );
  endstate.setQOverP( qOverP );

  Gaudi::TrackSymMatrix seedCov;
  seedCov( 0, 0 ) = track.errbx() * track.errbx();
  seedCov( 2, 2 ) = track.errsx() * track.errsx();
  seedCov( 1, 1 ) = track.errby() * track.errby();
  seedCov( 3, 3 ) = track.errsy() * track.errsy();
  seedCov( 4, 4 ) = sigmaQOverP2;
  state.setCovariance( seedCov );
  endstate.setCovariance( seedCov );

  state.setLocation( LHCb::State::Location::Muon );
  endstate.setLocation( LHCb::State::Location::LastMeasurement );

  debug() << "Muon state = " << state << endmsg;

  lbtrack.clearStates();
  lbtrack.addToStates( state );
  lbtrack.addToStates( endstate );
  lbtrack.setChi2PerDoF( track.chi2x() + track.chi2y() );
  lbtrack.setNDoF( track.nHits() - 2 );

  for ( int i = 0; i < track.nHits(); i++ ) {
    const auto Tile = track.point( i ).tile();
    lbtrack.addToLhcbIDs( ( LHCb::LHCbID )( Tile ) );
    debug() << " Muon Hit " << i << " tile " << Tile << " tiles in station " << track.point( i ).station() << endmsg;
  }

  lbtrack.setPatRecStatus( LHCb::Track::PatRecStatus::PatRecIDs );
  lbtrack.setType( LHCb::Track::Types::Muon );
}

void StandaloneMuonRec::detectClone( std::vector<StandaloneMuonTrack>& muonTracks, const Cache& cache ) const {

  for ( auto itMuonTrackFirst = muonTracks.begin(); itMuonTrackFirst < muonTracks.end(); itMuonTrackFirst++ ) {
    for ( auto itMuonTrackSecond = itMuonTrackFirst + 1; itMuonTrackSecond < muonTracks.end(); itMuonTrackSecond++ ) {
      bool sameM2 = false;
      bool sameM3 = false;
      bool sameM4 = false;
      if ( itMuonTrackFirst->point( 0 ).tile() == itMuonTrackSecond->point( 0 ).tile() ) sameM2 = true;
      if ( itMuonTrackFirst->point( 1 ).tile() == itMuonTrackSecond->point( 1 ).tile() ) sameM3 = true;
      if ( itMuonTrackFirst->point( 2 ).tile() == itMuonTrackSecond->point( 2 ).tile() ) sameM4 = true;
      if ( ( sameM2 && sameM3 && sameM4 ) ||
           ( ( m_strongCloneKiller || itMuonTrackFirst->nHits() == 3 ) && sameM2 && sameM3 ) ) {
        const auto x_extra5 =
            -( itMuonTrackFirst->point( 1 ).x() - itMuonTrackFirst->point( 2 ).x() ) + itMuonTrackFirst->point( 2 ).x();
        unsigned int ihit = 2;
        if ( itMuonTrackFirst->nHits() == 4 && itMuonTrackFirst->nHits() == 4 ) { ihit = 3; }
        const auto y_extra5 = itMuonTrackFirst->point( ihit - 1 ).y() * cache.stationZ[ihit] / cache.stationZ[ihit - 1];

        const auto distuno =
            ( itMuonTrackFirst->point( ihit ).x() - x_extra5 ) * ( itMuonTrackFirst->point( ihit ).x() - x_extra5 ) +
            ( itMuonTrackFirst->point( ihit ).y() - y_extra5 ) * ( itMuonTrackFirst->point( ihit ).y() - y_extra5 );
        const auto distdue =
            ( itMuonTrackSecond->point( ihit ).x() - x_extra5 ) * ( itMuonTrackSecond->point( ihit ).x() - x_extra5 ) +
            ( itMuonTrackSecond->point( ihit ).y() - y_extra5 ) * ( itMuonTrackSecond->point( ihit ).y() - y_extra5 );

        if ( distuno > distdue ) {
          itMuonTrackFirst->setClone();
        } else {
          itMuonTrackSecond->setClone();
        }
      }
    }
  }
}
