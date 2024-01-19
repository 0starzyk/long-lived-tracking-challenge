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
#include "Event/RecVertex.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "IPVFitter.h"
#include "IPVSeeding.h"
#include "Kernel/STLExtensions.h"
#include "PVOfflineRecalculate.h"
#include "TrackInterfaces/IPVOfflineTool.h"
#include "VPDet/DeVP.h"

#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/ISequencerTimerTool.h"
#include "GaudiKernel/IUpdateManagerSvc.h"
#include "GaudiKernel/SystemOfUnits.h"

#include <functional>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

class PVOfflineTool : public extends<LHCb::DetDesc::ConditionAccessorHolder<GaudiTool>, IPVOfflineTool> {
public:
  // Standard constructor
  using extends::extends;
  // Destructor
  StatusCode initialize() override;
  // PV fitting

  StatusCode reDoSinglePV( const LHCb::Tracks& inputTracks, const Gaudi::XYZPoint xyzseed,
                           std::vector<const LHCb::Track*>& tracks2exclude, LHCb::RecVertex& outvtx,
                           IGeometryInfo const& geometry ) const override;

  StatusCode reDoMultiPV( const LHCb::Tracks& inputTracks, const LHCb::RecVertex& invtx,
                          std::vector<const LHCb::Track*>& tracks2exclude, LHCb::RecVertex& outvtx,
                          IGeometryInfo const& geometry ) const override;

  StatusCode reconstructSinglePVFromTracks( const Gaudi::XYZPoint                  xyzseed,
                                            const std::vector<const LHCb::Track*>& tracks2use, LHCb::RecVertex& outvtx,
                                            IGeometryInfo const& geometry ) const override;

  StatusCode reconstructMultiPVFromTracks( std::vector<const LHCb::Track*>& tracks2use,
                                           std::vector<LHCb::RecVertex>&    outvtxVec,
                                           IGeometryInfo const&             geometry ) const override;

  StatusCode reconstructMultiPV( const std::vector<LHCb::Track>& inputTracks, std::vector<LHCb::RecVertex>& outvtxVec,
                                 IGeometryInfo const& geometry ) const override;

  StatusCode reconstructSinglePV( const LHCb::Tracks& inputTracks, const Gaudi::XYZPoint xyzseed,
                                  LHCb::RecVertex& outvtx, IGeometryInfo const& geometry ) const override;

  StatusCode removeTracksAndRecalculatePV( const LHCb::RecVertex*                 pvin,
                                           const std::vector<const LHCb::Track*>& tracks2remove, LHCb::RecVertex& vtx,
                                           IGeometryInfo const& geometry ) const override;

private:
  Gaudi::Property<bool> m_requireVelo{this, "RequireVelo", true, "Option to use tracks with VELO segment only"};
  Gaudi::Property<bool> m_saveSeedsAsPV{this, "SaveSeedsAsPV", false, "Save seeds as PVs (for monitoring"};

  // Tools
  ToolHandle<IPVFitter>            m_pvfit{this, "PVFitterName", "AdaptivePV3DFitter"}; // PV fitting tool
  ToolHandle<IPVSeeding>           m_pvSeedTool{this, "PVSeedingName", "PVSeed3DTool"}; // Seeding tool
  ToolHandle<PVOfflineRecalculate> m_pvRecalc{this, "PVOfflineRecalculate", "PVOfflineRecalculate"};

  Gaudi::Property<double> m_pvsChi2Separation{this, "PVsChi2Separation", 25.};
  Gaudi::Property<double> m_pvsChi2SeparationLowMult{this, "PVsChi2SeparationLowMult", 91.};

  Gaudi::Property<bool>         m_useBeamSpotRCut{this, "UseBeamSpotRCut", false};
  Gaudi::Property<double>       m_beamSpotRCut{this, "BeamSpotRCut", 0.2};
  Gaudi::Property<double>       m_beamSpotRCutHMC{this, "BeamSpotRHighMultiplicityCut", 0.4};
  Gaudi::Property<unsigned int> m_beamSpotRMT{this, "BeamSpotRMultiplicityTreshold", 10};

  ConditionAccessor<DeVP> m_vp{this, "DeVP", LHCb::Det::VP::det_path};

  // Member functions
  LHCb::RecVertex* matchVtxByTracks( const LHCb::RecVertex& invtx, std::vector<LHCb::RecVertex>& outvtxvec ) const;

  std::vector<const LHCb::Track*> readTracks( const LHCb::Tracks& inputTracks ) const;

  void removeTracksByLHCbIDs( std::vector<const LHCb::Track*>& tracks,
                              LHCb::span<const LHCb::Track*>   tracks2remove ) const;

  void removeTracksUsedByVertex( std::vector<const LHCb::Track*>& tracks, LHCb::RecVertex& rvtx ) const;

  // timing
  Gaudi::Property<bool>           m_doTiming{this, "TimingMeasurement", false};
  ToolHandle<ISequencerTimerTool> m_timerTool{"SequencerTimerTool/Timer", this};
  std::array<int, 3>              m_timer;

  enum class timers_t { Total = 0, Seeding, Fitting };
  class TimerGuard {
    ISequencerTimerTool* m_tool;
    int                  m_timer;

  public:
    TimerGuard( const ISequencerTimerTool* t, int i ) : m_tool( const_cast<ISequencerTimerTool*>( t ) ), m_timer( i ) {
      m_tool->start( m_timer );
    }
    ~TimerGuard() { m_tool->stop( m_timer ); }
  };
  std::optional<TimerGuard> make_timeguard( timers_t type ) const {
    if ( !m_doTiming ) return {};
    return TimerGuard{m_timerTool.get(), m_timer[static_cast<int>( type )]};
  }
  // trivial accessor to minimum allowed chi2
  double minAllowedChi2( const LHCb::RecVertex& rvtx ) const {
    return ( rvtx.tracks().size() < 7 ? std::max( m_pvsChi2Separation, m_pvsChi2SeparationLowMult )
                                      : m_pvsChi2Separation );
  }
};

DECLARE_COMPONENT( PVOfflineTool )

namespace {

  template <typename Iterator, typename Projection, typename Threshold>
  std::pair<Iterator, Threshold> max_above_threshold( Iterator begin, Iterator end, Projection proj, Threshold mx ) {
    auto r = std::make_pair( end, std::move( mx ) );
    for ( ; begin != end; ++begin ) {
      auto val = std::invoke( proj, *begin );
      if ( val > r.second ) {
        r.second = std::move( val );
        r.first  = begin;
      }
    }
    return r;
  }

  class velo_overlap_with {
    std::vector<LHCb::LHCbID> m_ref;

  public:
    velo_overlap_with( LHCb::span<const LHCb::LHCbID> ids ) {
      m_ref.reserve( ids.size() );
      std::copy_if( ids.begin(), ids.end(), std::back_inserter( m_ref ),
                    []( const LHCb::LHCbID& id ) { return id.isVP(); } );
    }
    velo_overlap_with( const LHCb::Track& trk ) : velo_overlap_with( trk.lhcbIDs() ) {}

    std::pair<int, int> operator()( LHCb::span<const LHCb::LHCbID> ids ) const {
      auto first1 = m_ref.begin();
      auto last1  = m_ref.end();
      auto first2 = ids.begin();
      auto last2  = ids.end();
      auto n      = std::make_pair( 0, 0 );
      while ( first1 != last1 && first2 != last2 ) {
        if ( *first1 < *first2 ) {
          ++first1;
        } else if ( *first2 < *first1 ) {
          if ( first2->isVP() ) ++n.second;
          ++first2;
        } else {
          ++n.first;
          ++n.second;
          ++first1;
          ++first2;
        }
      }
      return n;
    }
    std::pair<int, int> operator()( const LHCb::Track& trk ) const { return ( *this )( trk.lhcbIDs() ); }
  };

  double zCloseBeam( const LHCb::Track& track ) {
    Gaudi::XYZVector   unitVect = track.firstState().slopes().Unit();
    const LHCb::State& stateG   = track.firstState();
    return stateG.z() - unitVect.z() * ( unitVect.x() * stateG.x() + unitVect.y() * stateG.y() ) /
                            ( 1.0 - std::pow( unitVect.z(), 2 ) );
  }
  //=============================================================================
  // removeTracks
  //=============================================================================
  void removeTracks( std::vector<const LHCb::Track*>& tracks, LHCb::span<const LHCb::Track*> tracks2remove ) {
    auto firstToErase = std::remove_if( begin( tracks ), end( tracks ), [tracks2remove]( auto trk ) {
      return std::find( tracks2remove.begin(), tracks2remove.end(), trk ) != tracks2remove.end();
    } );
    tracks.erase( firstToErase, std::end( tracks ) );
  }

  void removeTracks( std::vector<const LHCb::Track*>& tracks, const SmartRefVector<LHCb::Track>& tracks2remove ) {
    auto firstToErase = std::remove_if( begin( tracks ), end( tracks ), [&tracks2remove]( auto trk ) {
      return std::find_if( tracks2remove.begin(), tracks2remove.end(),
                           [&trk]( auto item ) { return item.target() == trk; } ) != tracks2remove.end();
    } );
    tracks.erase( firstToErase, std::end( tracks ) );
  }

  //=============================================================================
  // Store dummy vertices
  //=============================================================================
  std::vector<LHCb::RecVertex> storeDummyVertices( const std::vector<Gaudi::XYZPoint>& seeds,
                                                   LHCb::span<const LHCb::Track*>      rtracks ) {
    std::vector<LHCb::RecVertex> out;
    out.reserve( seeds.size() );
    for ( const auto& seed : seeds ) {
      out.emplace_back();
      auto& tVertex = out.back();
      tVertex.setPosition( seed );
      Gaudi::SymMatrix3x3 errMat;
      for ( int i = 0; i < 3; i++ ) {
        for ( int j = 0; j < 3; j++ ) { errMat( i, j ) = 1.0; }
      }
      tVertex.setCovMatrix( errMat );
      tVertex.setNDoF( 1 );
      tVertex.setChi2( 99999.0 );
      // Fill close tracks
      auto is_close = [z = seed.Z()]( const auto* trk ) {
        return std::abs( zCloseBeam( *trk ) - z ) < 3.0 * Gaudi::Units::mm;
      };
      for ( const auto& trk : rtracks ) {
        if ( is_close( trk ) ) tVertex.addToTracks( trk );
      }
    }
    return out;
  }
} // namespace

//=========================================================================
// Initialize
//=========================================================================
StatusCode PVOfflineTool::initialize() {
  return extends::initialize().andThen( [&] {
    if ( m_doTiming ) {
      m_timerTool.retrieve().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      m_timer[static_cast<int>( timers_t::Total )] = m_timerTool->addTimer( "PatPV total" );
      m_timerTool->increaseIndent();
      m_timer[static_cast<int>( timers_t::Seeding )] = m_timerTool->addTimer( "PatPV seeding" );
      m_timer[static_cast<int>( timers_t::Fitting )] = m_timerTool->addTimer( "PatPV fitting" );
      m_timerTool->decreaseIndent();
    } else {
      m_timerTool.disable();
    }
    return StatusCode::SUCCESS;
  } );
}

//=============================================================================
// reconstruct PV for a given seed
//=============================================================================
StatusCode PVOfflineTool::reconstructSinglePV( const LHCb::Tracks& inputTracks, const Gaudi::XYZPoint xyzseed,
                                               LHCb::RecVertex& outvtx, IGeometryInfo const& geometry ) const {
  std::vector<const LHCb::Track*> tracks2remove;
  return m_pvfit->fitVertex( xyzseed, readTracks( inputTracks ), outvtx, tracks2remove, geometry );
}
//=============================================================================
// reconstruct PV for a given seed using tracks from the list
//=============================================================================
StatusCode PVOfflineTool::reconstructSinglePVFromTracks( const Gaudi::XYZPoint                  xyzseed,
                                                         const std::vector<const LHCb::Track*>& tracks2use,
                                                         LHCb::RecVertex&                       outvtx,
                                                         IGeometryInfo const&                   geometry ) const {
  std::vector<const LHCb::Track*> tracks2remove;
  return m_pvfit->fitVertex( xyzseed, tracks2use, outvtx, tracks2remove, geometry );
}

//=============================================================================
// reconstruct PV for a given seed with a list of tracks to be excluded
//=============================================================================
StatusCode PVOfflineTool::reDoSinglePV( const LHCb::Tracks& inputTracks, const Gaudi::XYZPoint xyzseed,
                                        std::vector<const LHCb::Track*>& tracks2exclude, LHCb::RecVertex& outvtx,
                                        IGeometryInfo const& geometry ) const {
  auto rtracks = readTracks( inputTracks );
  if ( !tracks2exclude.empty() ) removeTracksByLHCbIDs( rtracks, tracks2exclude );
  std::vector<const LHCb::Track*> tracks2remove;
  return m_pvfit->fitVertex( xyzseed, rtracks, outvtx, tracks2remove, geometry );
}

//=============================================================================
// multi vtx search and fit. Return new vtx after track removal
//=============================================================================
StatusCode PVOfflineTool::reDoMultiPV( const LHCb::Tracks& inputTracks, const LHCb::RecVertex& invtx,
                                       std::vector<const LHCb::Track*>& tracks2exclude, LHCb::RecVertex& outvtx,
                                       IGeometryInfo const& geometry ) const {
  auto rtracks = readTracks( inputTracks );

  if ( !tracks2exclude.empty() ) removeTracksByLHCbIDs( rtracks, tracks2exclude );

  std::vector<LHCb::RecVertex> outvtxVec;
  StatusCode                   scvfit = reconstructMultiPVFromTracks( rtracks, outvtxVec, geometry );

  // check which vtx corresponds to input vtx
  if ( scvfit.isFailure() ) return scvfit;
  auto pvtx = matchVtxByTracks( invtx, outvtxVec );
  if ( !pvtx ) return StatusCode::FAILURE;
  outvtx = *pvtx;
  return StatusCode::SUCCESS;
}

//=============================================================================
// multi vtx search and fit with tracks from default location
//=============================================================================
StatusCode PVOfflineTool::reconstructMultiPV( const std::vector<LHCb::Track>& inputTracks,
                                              std::vector<LHCb::RecVertex>&   outvtxvec,
                                              IGeometryInfo const&            geometry ) const {
  std::vector<const LHCb::Track*> rtracks;

  std::for_each( inputTracks.begin(), inputTracks.end(), [&]( const LHCb::Track& trk ) {
    if ( !m_requireVelo || trk.hasVelo() ) { rtracks.push_back( &trk ); }
  } );

  if ( msgLevel( MSG::DEBUG ) ) debug() << "readTracks: " << rtracks.size() << endmsg;

  return reconstructMultiPVFromTracks( rtracks, outvtxvec, geometry );
}

namespace {
  using RecVertexSpan = LHCb::span<const LHCb::RecVertex>;
  bool isChi2Separated( const LHCb::RecVertex& rvtx, RecVertexSpan outvtxvec, double minAllowedChi2 ) {
    return std::none_of(
        outvtxvec.begin(), outvtxvec.end(),
        [rz = rvtx.position().z(), sigma2z = rvtx.covMatrix()( 2, 2 ), minAllowedChi2]( const LHCb::RecVertex& v ) {
          return std::pow( rz - v.position().z(), 2 ) / ( sigma2z + v.covMatrix()( 2, 2 ) ) < minAllowedChi2;
        } );
  }
} // namespace

//=============================================================================
// multi vtx search and fit with tracks specified
//=============================================================================
StatusCode PVOfflineTool::reconstructMultiPVFromTracks( std::vector<const LHCb::Track*>& rtracks,
                                                        std::vector<LHCb::RecVertex>&    outvtxvec,
                                                        IGeometryInfo const&             geometry ) const {
  auto            totaltime_guard = make_timeguard( timers_t::Total );
  DeVP const&     deVP            = m_vp.get();
  Gaudi::XYZPoint beamSpot        = deVP.beamSpot();
  outvtxvec.clear();

  if ( m_saveSeedsAsPV ) {
    auto seeds = m_pvSeedTool->getSeeds( rtracks, beamSpot, geometry );
    outvtxvec  = storeDummyVertices( seeds, rtracks );
    return StatusCode::SUCCESS;
  }

  bool goOn = true;
  while ( goOn ) {
    goOn       = false;
    auto seeds = [&]() {
      auto seedingtime_guard = make_timeguard( timers_t::Seeding );
      // seeding
      return m_pvSeedTool->getSeeds( rtracks, beamSpot, geometry );
    }();

    if ( msgLevel( MSG::DEBUG ) ) { debug() << " seeds  " << seeds.size() << " pos: " << seeds << endmsg; }
    for ( const auto& seed : seeds ) {
      if ( msgLevel( MSG::VERBOSE ) ) verbose() << "ready to fit" << endmsg;
      LHCb::RecVertex&                recvtx = outvtxvec.emplace_back(); // create in place in case we will keep it
      std::vector<const LHCb::Track*> tracks2remove;
      {
        auto fittime_guard = make_timeguard( timers_t::Fitting );
        // fitting
        StatusCode scvfit = m_pvfit->fitVertex( seed, rtracks, recvtx, tracks2remove, geometry );
        if ( !scvfit.isSuccess() ) {
          outvtxvec.pop_back();
          continue;
        }
      }
      assert( !outvtxvec.empty() ); // code below implicitly assumes not empty !!
      if ( !isChi2Separated( outvtxvec.back(),
                             {outvtxvec.data(), static_cast<RecVertexSpan::size_type>( outvtxvec.size() - 1 )},
                             minAllowedChi2( outvtxvec.back() ) ) ) {
        outvtxvec.pop_back();
        continue;
      }
      if ( m_useBeamSpotRCut.value() && deVP.veloClosed() ) {
        const auto& pos = recvtx.position();
        auto        r2  = std::pow( pos.x() - beamSpot.x(), 2 ) + std::pow( pos.y() - beamSpot.y(), 2 );
        auto        r =
            ( recvtx.tracks().size() < m_beamSpotRMT.value() ? m_beamSpotRCut.value() : m_beamSpotRCutHMC.value() );
        if ( r2 > r * r ) {
          outvtxvec.pop_back();
          continue;
        }
      }
      goOn = true;
      removeTracks( rtracks, tracks2remove );
    } // iterate on seeds
  }   // iterate on vtx

  return StatusCode::SUCCESS;
}

//=============================================================================
// Read tracks
//=============================================================================
std::vector<const LHCb::Track*> PVOfflineTool::readTracks( const LHCb::Tracks& inputTracks ) const {
  std::vector<const LHCb::Track*> rtracks;
  rtracks.reserve( inputTracks.size() );

  std::copy_if( inputTracks.begin(), inputTracks.end(), std::back_inserter( rtracks ),
                [=]( const LHCb::Track* trk ) { return !m_requireVelo || trk->hasVelo(); } );

  if ( msgLevel( MSG::DEBUG ) ) debug() << "readTracks: " << rtracks.size() << endmsg;

  return rtracks;
}

//=============================================================================
// removeTracks
//=============================================================================
void PVOfflineTool::removeTracksByLHCbIDs( std::vector<const LHCb::Track*>& tracks,
                                           LHCb::span<const LHCb::Track*>   tracks2remove ) const {
  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "removeTracksByLHCbIDs. Input number of tracks: " << tracks.size() << endmsg;
  }

  for ( const auto* ptr1 : tracks2remove ) {
    auto compute_overlap = velo_overlap_with( *ptr1 );
    auto i               = std::find_if( tracks.begin(), tracks.end(), [&]( const LHCb::Track* trk ) {
      auto n = compute_overlap( *trk );
      return 1. * n.first > 0.99 * n.second;
    } );
    if ( i != tracks.end() ) tracks.erase( i );
  } // over tracks2remove

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "removeTracksByLHCbIDs. Output number of tracks: " << tracks.size() << endmsg;
  }
}

//=============================================================================
// removeTracksUsedByVertex
//=============================================================================
void PVOfflineTool::removeTracksUsedByVertex( std::vector<const LHCb::Track*>& tracks, LHCb::RecVertex& rvtx ) const {
  removeTracks( tracks, rvtx.tracks() );
}

//=============================================================================
// Match vtx in vector of vtx by matching tracks
//=============================================================================
LHCb::RecVertex* PVOfflineTool::matchVtxByTracks( const LHCb::RecVertex&        invtx,
                                                  std::vector<LHCb::RecVertex>& outvtxvec ) const {
  auto is_in = []( const auto& haystack ) {
    return [&]( const LHCb::Track* needle ) {
      return std::any_of( haystack.begin(), haystack.end(), [&]( const LHCb::Track* i ) { return i == needle; } );
    };
  };
  const auto& tracksIn = invtx.tracks();
  auto        best     = max_above_threshold(
      outvtxvec.begin(), outvtxvec.end(),
      [&]( const auto& ovtx ) { return std::count_if( tracksIn.begin(), tracksIn.end(), is_in( ovtx.tracks() ) ); },
      std::lround( std::floor( 0.3 * tracksIn.size() ) ) );
  if ( best.first == outvtxvec.end() ) return nullptr;
  // if tracksIn.empty(), we cannot get here... so dividing by tracksIn.size() will always be possible
  if ( msgLevel( MSG::DEBUG ) )
    debug() << " vtx succesfully matched at tracks rate: " << double( best.second ) / tracksIn.size() << endmsg;
  return &*best.first;
}

StatusCode PVOfflineTool::removeTracksAndRecalculatePV( const LHCb::RecVertex*                 pvin,
                                                        const std::vector<const LHCb::Track*>& tracks2remove,
                                                        LHCb::RecVertex& vtx, IGeometryInfo const& geometry ) const {
  return m_pvRecalc->RecalculateVertex( pvin, tracks2remove, vtx, geometry );
}
