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

#include "DetDesc/Condition.h"
#include "Event/RecVertex.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "LHCbAlgs/Transformer.h"
#include "VPDet/DeVP.h"

#include "GaudiAlg/ISequencerTimerTool.h"
#include "GaudiKernel/SystemOfUnits.h"

#include <assert.h>
#include <string>
#include <vector>

namespace {

  // auxiliary class for searching of clusters of tracks
  struct vtxCluster final {
    double z          = 0; // z of the cluster
    double sigsq      = 0; // sigma**2 of the cluster
    double sigsqmin   = 0; // minimum sigma**2 of the tracks forming cluster
    int    ntracks    = 0; // number of tracks in the cluster
    int    not_merged = 0; // flag for iterative merging
  };

  class PVTrack final {
  public:
    using Track           = LHCb::Event::v1::Track;
    const Track* refTrack = nullptr;
    // Current state of the track at the current point
    LHCb::State stateG;
    // Normalized vector of slope
    Gaudi::XYZVector unitVect;
    // Flag if the track has been used in a previous vertex
    bool isUsed = false;

    // Result for impact parameter
    Gaudi::XYZVector vd0;        // Impact parameter vector
    double           d0sq   = 0; // Impact parameter squared
    double           err2d0 = 0; // IP error squared
    double           chi2   = 0; // chi2 = d02 / d0err**2
    double           weight = 0; // Weight assigned to track
    Track::Types     type   = Track::Types::Velo;

    double zClose() const {
      return ( stateG.z() - ( unitVect.z() * ( unitVect.x() * stateG.x() + unitVect.y() * stateG.y() ) /
                              ( 1.0 - std::pow( unitVect.z(), 2 ) ) ) );
    }
  };

  using PVTracks    = std::vector<PVTrack>;
  using PVTrackPtrs = std::vector<PVTrack*>;

  struct PVVertex final {
    PVTrackPtrs     pvTracks;
    LHCb::RecVertex primVtx;
  };

} // namespace

/**
 *  Algorithm to find the primary vertices at the HLT.
 *
 *  This version is temporary and should be dropped and replaced
 *  with TrackBeamLineVertexFinder when it is ready. Its only purpose
 *  is to provide a PatPV using Tracks v2 and RecVertex v2 in the
 *  mean time.
 *  Note that this class integrates directly the code of the tools
 *  previously used by PatPV3D
 */
using TrackContainer = LHCb::Event::v1::Track::Range;
class PatPV3DFuture : public LHCb::Algorithm::Transformer<LHCb::RecVertices( const TrackContainer&, DeVP const& ),
                                                          LHCb::DetDesc::usesConditions<DeVP>> {
public:
  using RecVertex   = LHCb::RecVertex;
  using RecVertices = LHCb::RecVertices;
  using Track       = LHCb::Event::v1::Track;

  PatPV3DFuture( const std::string& name, ISvcLocator* pSvcLocator );

  RecVertices operator()( const TrackContainer&, DeVP const& ) const override;

private:
  ////////  From PVOfflineTool
  PatPV3DFuture::RecVertices reconstructMultiPV( const TrackContainer&  inputTracks,
                                                 Gaudi::XYZPoint const& beamSpot ) const;

  PatPV3DFuture::RecVertices reconstructMultiPVFromTracks( std::vector<const PatPV3DFuture::Track*>& tracks2use,
                                                           Gaudi::XYZPoint const&                    beamSpot ) const;

  void removeTracks( std::vector<const PatPV3DFuture::Track*>&       tracks,
                     const std::vector<const PatPV3DFuture::Track*>& tracks2remove ) const;

  ////////  From AdaptivePV3DFitter
  StatusCode fitVertex( const Gaudi::XYZPoint& seedPoint, const std::vector<const Track*>& tracks,
                        PatPV3DFuture::RecVertices& outputVtxVec, std::vector<const Track*>& tracks2remove ) const;

  // Get Tukey's weight
  double getTukeyWeight( double trchi2, int iter ) const;

  ////////  From PVSeedTool
  std::vector<Gaudi::XYZPoint> getSeeds( const std::vector<const PatPV3DFuture::Track*>& inputTracks,
                                         const Gaudi::XYZPoint&                          beamspot ) const;

  std::vector<double> findClusters( std::vector<vtxCluster>& vclus ) const;
  double              errorForPVSeedFinding( double tx, double ty ) const;

private:
  ////////  From PVOfflineTool
  Gaudi::Property<bool>   m_refitpv{this, "RefitPV", false, "Flag to refit PVs when converting to type PrimaryVertex"};
  Gaudi::Property<bool>   m_requireVelo{this, "RequireVelo", true, "Option to use tracks with VELO segment only"};
  Gaudi::Property<double> m_pvsChi2Separation{this, "PVsChi2Separation", 25.};
  Gaudi::Property<double> m_pvsChi2SeparationLowMult{this, "PVsChi2SeparationLowMult", 91.};
  Gaudi::Property<bool>   m_useBeamSpotRCut{this, "UseBeamSpotRCut", false};
  Gaudi::Property<double> m_beamSpotRCut{this, "BeamSpotRCut", 0.2};
  Gaudi::Property<double> m_beamSpotRCutHMC{this, "BeamSpotRHighMultiplicityCut", 0.4};
  Gaudi::Property<unsigned int> m_beamSpotRMT{this, "BeamSpotRMultiplicityTreshold", 10};
  Gaudi::Property<double>       m_resolverBound{this, "ResolverBound", 5 * Gaudi::Units::mm};
  bool                          m_veloClosed = false;

  ////////  From AdaptivePV3DFitter
  Gaudi::Property<size_t> m_minTr{this, "MinTracks", 4, "Minimum number of tracks to make a vertex"};
  Gaudi::Property<int>    m_Iterations{this, "Iterations", 20, "Number of iterations for minimisation"};
  Gaudi::Property<int>    m_minIter{this, "minIter", 5, "Min number of iterations"};
  Gaudi::Property<double> m_maxDeltaZ{this, "maxDeltaZ", 0.0005 * Gaudi::Units::mm,
                                      "Fit convergence condition"}; // 0.001 * Gaudi::Units::mm)
  Gaudi::Property<double> m_minTrackWeight{this, "minTrackWeight", 0.00000001,
                                           "Minimum Tukey's weight to accept a track"}; // 0.00001)
  Gaudi::Property<double> m_TrackErrorScaleFactor{this, "TrackErrorScaleFactor", 1.0};
  Gaudi::Property<double> m_maxChi2{this, "maxChi2", 400.0, "max chi2 of track-to-vtx to be considered for fit"};
  Gaudi::Property<double> m_trackMaxChi2{this, "trackMaxChi2", 12.};
  double                  m_trackChi; // sqrt of trackMaxChi2
  Gaudi::Property<double> m_trackMaxChi2Remove{
      this,
      "trackMaxChi2Remove",
      25.,
      [this]( auto& ) { this->m_trackChi = std::sqrt( this->m_trackMaxChi2 ); },
      Gaudi::Details::Property::ImmediatelyInvokeHandler{true}, // insure m_trackChi is in sync with m_trackMaxChi2
      "Max chi2 tracks to be removed from next PV search"};
  Gaudi::Property<double> m_maxDeltaZCache{this, "maxDeltaZCache", 1 * Gaudi::Units::mm,
                                           "Update derivatives if distance of reference is more than this"};

  ////////  From PVSeedTool
  // steering parameters for merging procedure
  Gaudi::Property<double> m_maxChi2Merge{this, "maxChi2Merge", 25.};
  Gaudi::Property<double> m_factorToIncreaseErrors{this, "factorToIncreaseErrors", 15.};
  // steering parameters for final cluster selection
  Gaudi::Property<int>    m_minClusterMult{this, "minClusterMult", 3};
  Gaudi::Property<double> m_dzCloseTracksInCluster{this, "dzCloseTracksInCluster", 5. * Gaudi::Units::mm};
  Gaudi::Property<int>    m_minCloseTracksInCluster{this, "minCloseTracksInCluster", 3};
  Gaudi::Property<int>    m_highMult{this, "highMult", 10};
  Gaudi::Property<double> m_ratioSig2HighMult{this, "ratioSig2HighMult", 1.0};
  Gaudi::Property<double> m_ratioSig2LowMult{this, "ratioSig2LowMult", 0.9};
  double                  m_scatCons = 0; // calculated from m_x0MS
  Gaudi::Property<double> m_x0MS{this, "x0MS", 0.01,
                                 [=]( auto& ) {
                                   double X0        = this->m_x0MS;
                                   this->m_scatCons = ( 13.6 * sqrt( X0 ) * ( 1. + 0.038 * log( X0 ) ) );
                                 },
                                 Gaudi::Details::Property::ImmediatelyInvokeHandler{true}}; // X0 (tunable) of MS to add
                                                                                            // for extrapolation of
                                                                                            // track parameters to PV

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
  double minAllowedChi2( const RecVertex& rvtx ) const {
    return ( rvtx.tracks().size() < 7 ? std::max( m_pvsChi2Separation, m_pvsChi2SeparationLowMult )
                                      : m_pvsChi2Separation );
  }

  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_nbPVsCounter{this, "Nb PVs"};
};

DECLARE_COMPONENT( PatPV3DFuture )

PatPV3DFuture::PatPV3DFuture( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator,
                   {KeyValue{"InputTracks", LHCb::TrackLocation::Default}, KeyValue{"DeVP", LHCb::Det::VP::det_path}},
                   KeyValue( "OutputVerticesName", LHCb::RecVertexLocation::Velo3D ) ) {}

namespace {

  using RecVertexSpan = LHCb::span<const PatPV3DFuture::RecVertex>;
  bool isChi2Separated( const PatPV3DFuture::RecVertex& rvtx, const LHCb::RecVertices& outvtxvec,
                        double minAllowedChi2 ) {
    return std::none_of( outvtxvec.begin(), outvtxvec.end() - 1,
                         [rz = rvtx.position().z(), sigma2z = rvtx.covMatrix()( 2, 2 ),
                          minAllowedChi2]( const PatPV3DFuture::RecVertex* v ) {
                           return std::pow( rz - v->position().z(), 2 ) / ( sigma2z + v->covMatrix()( 2, 2 ) ) <
                                  minAllowedChi2;
                         } );
  }

  bool vtxcomp( vtxCluster* first, vtxCluster* second ) { return first->z < second->z; }
  bool multcomp( vtxCluster* first, vtxCluster* second ) { return first->ntracks > second->ntracks; }

  // auxiliary class for merging procedure of tracks/clusters
  struct pair_to_merge final {
    vtxCluster* first    = nullptr; // pointer to first cluster to be merged
    vtxCluster* second   = nullptr; // pointer to second cluster to be merged
    double      chi2dist = 10.e+10; // a chi2dist = zdistance**2/(sigma1**2+sigma2**2)
    pair_to_merge( vtxCluster* f, vtxCluster* s, double chi2 ) : first( f ), second( s ), chi2dist( chi2 ) {}
  };

  bool paircomp( const pair_to_merge& first, const pair_to_merge& second ) { return first.chi2dist < second.chi2dist; }

  constexpr static const int s_p2mstatic = 5000;

  double zCloseBeamOfflineTool( const PatPV3DFuture::Track* track, const Gaudi::XYZPoint& beamspot,
                                const double maxr2 ) {
    Gaudi::XYZPoint  tpoint   = track->position();
    Gaudi::XYZVector tdir     = track->slopes();
    double           wx       = ( 1. + tdir.x() * tdir.x() ) / track->firstState().errX2();
    double           wy       = ( 1. + tdir.y() * tdir.y() ) / track->firstState().errY2();
    double           x0       = tpoint.x() - tpoint.z() * tdir.x() - beamspot.X();
    double           y0       = tpoint.y() - tpoint.z() * tdir.y() - beamspot.Y();
    double           den      = wx * tdir.x() * tdir.x() + wy * tdir.y() * tdir.y();
    double           zAtBeam  = -( wx * x0 * tdir.x() + wy * y0 * tdir.y() ) / den;
    double           xb       = tpoint.x() + tdir.x() * ( zAtBeam - tpoint.z() ) - beamspot.X();
    double           yb       = tpoint.y() + tdir.y() * ( zAtBeam - tpoint.z() ) - beamspot.Y();
    double           r2AtBeam = xb * xb + yb * yb;
    return r2AtBeam < maxr2 ? zAtBeam : 10e8;
  }

} // namespace

PatPV3DFuture::RecVertices PatPV3DFuture::operator()( TrackContainer const& inputTracks, DeVP const& vpdet ) const {
  auto rvts = reconstructMultiPV( inputTracks, vpdet.beamSpot() );
  m_nbPVsCounter += rvts.size();
  return rvts;
}

PatPV3DFuture::RecVertices PatPV3DFuture::reconstructMultiPV( const TrackContainer&  inputTracks,
                                                              Gaudi::XYZPoint const& beamSpot ) const {
  std::vector<const PatPV3DFuture::Track*> rtracks;
  std::for_each( inputTracks.begin(), inputTracks.end(), [&]( const PatPV3DFuture::Track* trk ) {
    if ( !m_requireVelo || trk->hasVelo() ) { rtracks.push_back( trk ); }
  } );

  return reconstructMultiPVFromTracks( rtracks, beamSpot );
}

void PatPV3DFuture::removeTracks( std::vector<const PatPV3DFuture::Track*>&       tracks,
                                  const std::vector<const PatPV3DFuture::Track*>& tracks2remove ) const {
  auto firstToErase = std::remove_if( begin( tracks ), end( tracks ), [&tracks2remove]( auto trk ) {
    return std::find( begin( tracks2remove ), end( tracks2remove ), trk ) != tracks2remove.end();
  } );
  tracks.erase( firstToErase, std::end( tracks ) );
}

namespace {
  class AdaptivePVTrack final {
  public:
    using Track = LHCb::Event::v1::Track;
    AdaptivePVTrack( const Track& track, const Gaudi::XYZPoint& vtx );
    void                       updateCache( const Gaudi::XYZPoint& vtx );
    double                     weight() const { return m_weight; }
    void                       setWeight( double w ) { m_weight = w; }
    const Gaudi::SymMatrix3x3& halfD2Chi2DX2() const { return m_halfD2Chi2DX2; }
    const Gaudi::Vector3&      halfDChi2DX() const { return m_halfDChi2DX; }
    double                     chi2() const { return m_chi2; }
    inline double              chi2( const Gaudi::XYZPoint& vtx ) const;
    const Track*               track() const { return m_track; }

  private:
    void                              invertCov();
    double                            m_weight;
    const Track*                      m_track;
    LHCb::State                       m_state;
    Gaudi::SymMatrix2x2               m_invcov;
    Gaudi::SymMatrix3x3               m_halfD2Chi2DX2;
    Gaudi::Vector3                    m_halfDChi2DX;
    double                            m_chi2;
    ROOT::Math::SMatrix<double, 3, 2> m_H;
  };

  AdaptivePVTrack::AdaptivePVTrack( const LHCb::Event::v1::Track& track, const Gaudi::XYZPoint& vtx )
      : m_track( &track ) {
    // get the state
    m_state = track.firstState();
    if ( m_state.location() != LHCb::State::Location::ClosestToBeam ) {
      const LHCb::State* closestToBeam = track.stateAt( LHCb::State::Location::ClosestToBeam );
      if ( closestToBeam ) m_state = *closestToBeam;
    }
    // do here things we could evaluate at z_seed. may add cov matrix here, which'd save a lot of time.
    m_H( 0, 0 ) = m_H( 1, 1 ) = 1;
    m_H( 2, 0 )               = -m_state.tx();
    m_H( 2, 1 )               = -m_state.ty();
    // update the cache
    updateCache( vtx );
  }

  void AdaptivePVTrack::invertCov() {
    auto invdet      = 1 / ( m_invcov( 0, 0 ) * m_invcov( 1, 1 ) - m_invcov( 0, 1 ) * m_invcov( 1, 0 ) );
    auto new00       = m_invcov( 1, 1 ) * invdet;
    m_invcov( 1, 1 ) = m_invcov( 0, 0 ) * invdet;
    m_invcov( 0, 0 ) = new00;
    m_invcov( 0, 1 ) = m_invcov( 1, 0 ) = -m_invcov( 0, 1 ) * invdet;
  }

  void AdaptivePVTrack::updateCache( const Gaudi::XYZPoint& vtx ) {
    // transport to vtx z
    m_state.linearTransportTo( vtx.z() );

    // invert cov matrix
    m_invcov = m_state.covariance().Sub<Gaudi::SymMatrix2x2>( 0, 0 );
    invertCov();

    // The following can all be written out, omitting the zeros, once
    // we know that it works.
    Gaudi::Vector2                    res{vtx.x() - m_state.x(), vtx.y() - m_state.y()};
    ROOT::Math::SMatrix<double, 3, 2> HW = m_H * m_invcov;
    ROOT::Math::AssignSym::Evaluate( m_halfD2Chi2DX2, HW * ROOT::Math::Transpose( m_H ) );
    // m_halfD2Chi2DX2 = ROOT::Math::Similarity(H, invcov ) ;
    m_halfDChi2DX = HW * res;
    m_chi2        = ROOT::Math::Similarity( res, m_invcov );
  }

  inline double AdaptivePVTrack::chi2( const Gaudi::XYZPoint& vtx ) const {
    double         dz = vtx.z() - m_state.z();
    Gaudi::Vector2 res{vtx.x() - ( m_state.x() + dz * m_state.tx() ), vtx.y() - ( m_state.y() + dz * m_state.ty() )};
    return ROOT::Math::Similarity( res, m_invcov );
  }
} // namespace

StatusCode PatPV3DFuture::fitVertex( const Gaudi::XYZPoint& seedPoint, const std::vector<const Track*>& rTracks,
                                     PatPV3DFuture::RecVertices& outputVtxVec,
                                     std::vector<const Track*>&  tracks2remove ) const {
  tracks2remove.clear();
  // position at which derivatives are evaluated
  Gaudi::XYZPoint refpos = seedPoint;
  // prepare tracks
  std::vector<AdaptivePVTrack> pvTracks;
  pvTracks.reserve( rTracks.size() );
  for ( const auto& track : rTracks )
    if ( track->hasVelo() ) {
      pvTracks.emplace_back( *track, refpos );
      if ( pvTracks.back().chi2() >= m_maxChi2 ) pvTracks.pop_back();
    }
  if ( pvTracks.size() < m_minTr ) { return StatusCode::FAILURE; }

  // current vertex position
  Gaudi::XYZPoint vtxpos = refpos;
  // vertex covariance matrix
  Gaudi::SymMatrix3x3 vtxcov;
  bool                converged = false;
  double              maxdz     = m_maxDeltaZ;
  int                 nbIter    = 0;
  double              chi2      = 0;
  size_t              ntrin     = 0;

  while ( ( nbIter < m_minIter ) || ( !converged && nbIter < m_Iterations ) ) {
    ++nbIter;
    Gaudi::SymMatrix3x3 halfD2Chi2DX2;
    Gaudi::Vector3      halfDChi2DX;
    // update cache if too far from reference position. this is the slow part.
    if ( std::abs( refpos.z() - vtxpos.z() ) > m_maxDeltaZCache ) {
      refpos = vtxpos;
      for ( auto& trk : pvTracks ) trk.updateCache( refpos );
    }

    // add contribution from all tracks
    chi2  = 0;
    ntrin = 0;
    for ( auto& trk : pvTracks ) {
      // compute weight
      double trkchi2 = trk.chi2( vtxpos );
      double weight  = getTukeyWeight( trkchi2, nbIter );
      trk.setWeight( weight );
      // add the track
      if ( weight > m_minTrackWeight ) {
        ++ntrin;
        halfD2Chi2DX2 += weight * trk.halfD2Chi2DX2();
        halfDChi2DX += weight * trk.halfDChi2DX();
        chi2 += weight * trk.chi2();
      }
    }

    // check nr of tracks that entered the fit
    if ( ntrin < m_minTr ) { return StatusCode::FAILURE; }

    // compute the new vertex covariance
    vtxcov = halfD2Chi2DX2;
    if ( !vtxcov.InvertChol() ) { return StatusCode::FAILURE; }
    // compute the delta w.r.t. the reference
    Gaudi::Vector3 delta = -1.0 * vtxcov * halfDChi2DX;

    // note: this is only correct if chi2 was chi2 of reference!
    chi2 += ROOT::Math::Dot( delta, halfDChi2DX );

    // deltaz needed for convergence
    const double deltaz = refpos.z() + delta( 2 ) - vtxpos.z();

    // update the position
    vtxpos.SetX( refpos.x() + delta( 0 ) );
    vtxpos.SetY( refpos.y() + delta( 1 ) );
    vtxpos.SetZ( refpos.z() + delta( 2 ) );

    // loose convergence criteria if close to end of iterations
    if ( 1. * nbIter > 0.8 * m_Iterations ) maxdz = 10. * m_maxDeltaZ;
    converged = std::abs( deltaz ) < maxdz;

  } // end iteration loop
  if ( !converged ) return StatusCode::FAILURE;

  LHCb::RecVertex* vtx = new LHCb::RecVertex();
  vtx->setPosition( vtxpos );
  vtx->setCovMatrix( vtxcov );
  vtx->setChi2AndDoF( chi2, 2 * ntrin - 3 );

  for ( const auto& trk : pvTracks ) {
    if ( trk.weight() > m_minTrackWeight ) vtx->addToTracks( trk.track(), trk.weight() );
    // remove track for next PV search
    if ( trk.chi2( vtxpos ) < m_trackMaxChi2Remove ) tracks2remove.push_back( trk.track() );
  }
  vtx->setTechnique( RecVertex::RecVertexType::Primary );
  outputVtxVec.insert( vtx );
  return StatusCode::SUCCESS;
}

double PatPV3DFuture::getTukeyWeight( double trchi2, int iter ) const {
  if ( iter < 1 ) return 1.;
  double ctrv = m_trackChi * std::max( m_minIter - iter, 1 );
  auto   sq   = std::pow( ctrv * m_TrackErrorScaleFactor, 2 );
  return sq > trchi2 ? std::pow( 1. - ( trchi2 / sq ), 2 ) : 0.;
}

double PatPV3DFuture::errorForPVSeedFinding( double tx, double ty ) const {
  // the seeding results depend weakly on this eror parametrization
  double invPMean2    = 1. / ( 3000. * Gaudi::Units::MeV * 3000. * Gaudi::Units::MeV );
  double tanTheta2    = tx * tx + ty * ty;
  double invSinTheta2 = 1. / tanTheta2 + 1.;
  // assume that first hit in VD at 8 mm
  double distr        = 8. * Gaudi::Units::mm;
  double dist2        = distr * distr * invSinTheta2;
  double sigma_ms2    = m_scatCons * m_scatCons * dist2 * invPMean2;
  double fslope2      = 0.0005 * 0.0005;
  double sigma_slope2 = fslope2 * dist2;
  return ( sigma_ms2 + sigma_slope2 ) * invSinTheta2;
}

std::vector<double> PatPV3DFuture::findClusters( std::vector<vtxCluster>& vclus ) const {
  std::vector<double> zclusters;
  if ( vclus.empty() ) return zclusters;
  std::vector<vtxCluster*> pvclus;
  pvclus.reserve( vclus.size() );
  for ( auto& itvtx : vclus ) {
    itvtx.sigsq *= m_factorToIncreaseErrors * m_factorToIncreaseErrors; // blow up errors
    itvtx.sigsqmin = itvtx.sigsq;
    pvclus.push_back( &itvtx );
  }
  std::sort( pvclus.begin(), pvclus.end(), vtxcomp );
  bool more_merging = true;
  while ( more_merging ) {
    // find pair of clusters for merging
    // refresh flag for this iteration
    for ( auto ivc = pvclus.begin(); ivc != pvclus.end(); ivc++ ) { ( *ivc )->not_merged = 1; }
    // for a given cluster look only up to a few consequtive ones to merge
    // "a few" might become a property?
    auto                       n_consequtive = std::min( 5, static_cast<int>( pvclus.size() ) );
    auto                       ivc2up        = std::next( pvclus.begin(), n_consequtive );
    std::vector<pair_to_merge> vecp2m;
    vecp2m.reserve( std::min( static_cast<int>( pvclus.size() ) * n_consequtive, s_p2mstatic ) );
    for ( auto ivc1 = pvclus.begin(); ivc1 != pvclus.end() - 1; ivc1++ ) {
      if ( ivc2up != pvclus.end() ) ++ivc2up;
      for ( auto ivc2 = std::next( ivc1 ); ivc2 != ivc2up; ivc2++ ) {
        double zdist    = ( *ivc1 )->z - ( *ivc2 )->z;
        double chi2dist = zdist * zdist / ( ( *ivc1 )->sigsq + ( *ivc2 )->sigsq );
        if ( chi2dist < m_maxChi2Merge && vecp2m.size() < s_p2mstatic ) {
          vecp2m.emplace_back( *ivc1, *ivc2, chi2dist );
        }
      }
    }
    // done if no more pairs to merge
    if ( vecp2m.empty() ) {
      more_merging = false;
    } else {
      // sort if number of pairs reasonable. Sorting increases efficency.
      if ( vecp2m.size() < 100 ) std::sort( vecp2m.begin(), vecp2m.end(), paircomp );
      // merge pairs
      for ( auto itp2m = vecp2m.begin(); itp2m != vecp2m.end(); itp2m++ ) {
        vtxCluster* pvtx1 = itp2m->first;
        vtxCluster* pvtx2 = itp2m->second;
        if ( pvtx1->not_merged == 1 && pvtx2->not_merged == 1 ) {
          double z1       = pvtx1->z;
          double z2       = pvtx2->z;
          double s1       = pvtx1->sigsq;
          double s2       = pvtx2->sigsq;
          double s1min    = pvtx1->sigsqmin;
          double s2min    = pvtx2->sigsqmin;
          double sigsqmin = s1min;
          if ( s2min < s1min ) sigsqmin = s2min;
          double w_inv    = ( s1 * s2 / ( s1 + s2 ) );
          double zmerge   = w_inv * ( z1 / s1 + z2 / s2 );
          pvtx1->z        = zmerge;
          pvtx1->sigsq    = w_inv;
          pvtx1->sigsqmin = sigsqmin;
          pvtx1->ntracks += pvtx2->ntracks;
          pvtx2->ntracks    = 0; // mark second cluster as used
          pvtx1->not_merged = 0;
          pvtx2->not_merged = 0;
        }
      }
      // remove clusters which where used
      pvclus.erase(
          std::remove_if( pvclus.begin(), pvclus.end(), []( const vtxCluster* cl ) { return cl->ntracks < 1; } ),
          pvclus.end() );
    }
  }
  // End of clustering.
  // Sort according to multiplicity
  std::sort( pvclus.begin(), pvclus.end(), multcomp );
  // Select good clusters.
  for ( auto ivc = pvclus.begin(); ivc != pvclus.end(); ivc++ ) {
    int n_tracks_close = 0;
    for ( auto itvtx = vclus.begin(); itvtx != vclus.end(); itvtx++ ) {
      if ( fabs( itvtx->z - ( *ivc )->z ) < m_dzCloseTracksInCluster ) n_tracks_close++;
    }
    double dist_to_closest = 1000000.;
    if ( pvclus.size() > 1 ) {
      for ( auto ivc1 = pvclus.begin(); ivc1 != pvclus.end(); ivc1++ ) {
        if ( ivc != ivc1 && ( fabs( ( *ivc1 )->z - ( *ivc )->z ) < dist_to_closest ) ) {
          dist_to_closest = fabs( ( *ivc1 )->z - ( *ivc )->z );
        }
      }
    }
    // ratio to remove clusters made of one low error track and many large error ones
    double rat     = ( *ivc )->sigsq / ( *ivc )->sigsqmin;
    bool   igood   = false;
    int    ntracks = ( *ivc )->ntracks;
    if ( ntracks >= m_minClusterMult ) {
      if ( dist_to_closest > 10. && rat < 0.95 ) igood = true;
      if ( ntracks >= m_highMult && rat < m_ratioSig2HighMult ) igood = true;
      if ( ntracks < m_highMult && rat < m_ratioSig2LowMult ) igood = true;
    }
    // veto
    if ( n_tracks_close < m_minCloseTracksInCluster ) igood = false;
    if ( igood ) zclusters.push_back( ( *ivc )->z );
  }
  //  print_clusters(pvclus);
  return zclusters;
}

PatPV3DFuture::RecVertices
PatPV3DFuture::reconstructMultiPVFromTracks( std::vector<const PatPV3DFuture::Track*>& rtracks,
                                             Gaudi::XYZPoint const&                    beamSpot ) const {
  PatPV3DFuture::RecVertices outvtxvec;
  outvtxvec.reserve( 20 );
  auto totaltime_guard = make_timeguard( timers_t::Total );
  bool goOn            = true;

  while ( goOn ) {
    goOn       = false;
    auto seeds = [&]() {
      auto seedingtime_guard = make_timeguard( timers_t::Seeding );
      // seeding
      return getSeeds( rtracks, beamSpot );
    }();
    for ( const auto& seed : seeds ) {
      std::vector<const Track*> tracks2remove;
      {
        auto fittime_guard = make_timeguard( timers_t::Fitting );
        // fitting
        StatusCode scvfit = fitVertex( seed, rtracks, outvtxvec, tracks2remove );
        if ( !scvfit.isSuccess() ) { continue; }
      }

      assert( !outvtxvec.empty() ); // code below implicitly assumes not empty !!
      if ( !isChi2Separated( **( outvtxvec.end() - 1 ), outvtxvec, minAllowedChi2( **( outvtxvec.end() - 1 ) ) ) ) {
        outvtxvec.erase( *( outvtxvec.end() - 1 ) );
        continue;
      }
      if ( m_useBeamSpotRCut.value() && m_veloClosed ) {
        const RecVertex& recvtx = **( outvtxvec.end() - 1 );
        const auto&      pos    = recvtx.position();
        auto             r2     = std::pow( pos.x() - beamSpot.x(), 2 ) + std::pow( pos.y() - beamSpot.y(), 2 );
        auto             r =
            ( recvtx.tracks().size() < m_beamSpotRMT.value() ? m_beamSpotRCut.value() : m_beamSpotRCutHMC.value() );
        if ( r2 > r * r ) {
          outvtxvec.erase( *( outvtxvec.end() - 1 ) );
          continue;
        }
      }
      goOn = true;
      removeTracks( rtracks, tracks2remove );
    } // iterate on seeds
  }   // iterate on vtx
  return outvtxvec;
}

std::vector<Gaudi::XYZPoint> PatPV3DFuture::getSeeds( const std::vector<const PatPV3DFuture::Track*>& inputTracks,
                                                      const Gaudi::XYZPoint&                          beamspot ) const {
  std::vector<Gaudi::XYZPoint> seeds;
  if ( inputTracks.size() < 3 ) return seeds;
  std::vector<vtxCluster> vclusters;
  for ( const auto& trk : inputTracks ) {
    double sigsq;
    double zclu;

    const double maxr2 = m_useBeamSpotRCut ? 0.5 * 0.5 : 1e6;
    zclu               = zCloseBeamOfflineTool( trk, beamspot, maxr2 );
    sigsq              = errorForPVSeedFinding( trk->firstState().tx(), trk->firstState().ty() );

    if ( fabs( zclu ) > 2000. ) continue;
    vtxCluster clu;
    clu.z        = zclu;
    clu.sigsq    = sigsq;
    clu.sigsqmin = clu.sigsq;
    clu.ntracks  = 1;
    vclusters.push_back( clu );
  }
  auto zseeds = findClusters( vclusters );
  seeds.reserve( zseeds.size() );
  std::transform( zseeds.begin(), zseeds.end(), std::back_inserter( seeds ), [&]( double z ) {
    return Gaudi::XYZPoint{beamspot.X(), beamspot.Y(), z};
  } );
  return seeds;
}
