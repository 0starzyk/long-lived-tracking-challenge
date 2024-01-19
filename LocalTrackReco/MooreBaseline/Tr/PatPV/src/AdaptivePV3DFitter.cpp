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
#include "Event/RecVertex.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "IPVFitter.h"

class AdaptivePV3DFitter : public extends<GaudiTool, IPVFitter> {

public:
  // Standard constructor
  using extends::extends;
  // Fitting
  StatusCode fitVertex( const Gaudi::XYZPoint& seedPoint, LHCb::span<const LHCb::Track* const> tracks,
                        LHCb::RecVertex& vtx, std::vector<const LHCb::Track*>& tracks2remove,
                        IGeometryInfo const& geometry ) const override;

private:
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

  // Get Tukey's weight
  double getTukeyWeight( double trchi2, int iter ) const;
};

DECLARE_COMPONENT( AdaptivePV3DFitter )

namespace {
  class AdaptivePVTrack final {
  public:
    AdaptivePVTrack( const LHCb::Track& track, const Gaudi::XYZPoint& vtx );
    void                       updateCache( const Gaudi::XYZPoint& vtx );
    double                     weight() const { return m_weight; }
    void                       setWeight( double w ) { m_weight = w; }
    const Gaudi::SymMatrix3x3& halfD2Chi2DX2() const { return m_halfD2Chi2DX2; }
    const Gaudi::Vector3&      halfDChi2DX() const { return m_halfDChi2DX; }
    double                     chi2() const { return m_chi2; }
    inline double              chi2( const Gaudi::XYZPoint& vtx ) const;
    const LHCb::Track*         track() const { return m_track; }

  private:
    void                              invertCov();
    double                            m_weight;
    const LHCb::Track*                m_track;
    LHCb::State                       m_state;
    Gaudi::SymMatrix2x2               m_invcov;
    Gaudi::SymMatrix3x3               m_halfD2Chi2DX2;
    Gaudi::Vector3                    m_halfDChi2DX;
    double                            m_chi2;
    ROOT::Math::SMatrix<double, 3, 2> m_H;
  };

  AdaptivePVTrack::AdaptivePVTrack( const LHCb::Track& track, const Gaudi::XYZPoint& vtx ) : m_track( &track ) {
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

//=============================================================================
// Least square adaptive fitting method
//=============================================================================
StatusCode AdaptivePV3DFitter::fitVertex( const Gaudi::XYZPoint&               seedPoint,
                                          LHCb::span<const LHCb::Track* const> rTracks, LHCb::RecVertex& vtx,
                                          std::vector<const LHCb::Track*>& tracks2remove, IGeometryInfo const& ) const {
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

  if ( pvTracks.size() < m_minTr ) {
    if ( msgLevel( MSG::DEBUG ) ) debug() << "Too few tracks to fit PV" << endmsg;
    return StatusCode::FAILURE;
  }

  // current vertex position
  Gaudi::XYZPoint vtxpos = refpos;
  // vertex covariance matrix
  Gaudi::SymMatrix3x3 vtxcov;
  bool                converged = false;
  double              maxdz     = m_maxDeltaZ;
  int                 nbIter    = 0;
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
    double chi2( 0 );
    size_t ntrin( 0 );
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
    if ( ntrin < m_minTr ) {
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Too few tracks after PV fit" << endmsg;
      return StatusCode::FAILURE;
    }

    // compute the new vertex covariance
    vtxcov = halfD2Chi2DX2;
    if ( !vtxcov.InvertChol() ) {
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Error inverting hessian matrix" << endmsg;
      return StatusCode::FAILURE;
    }
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
    vtx.setChi2AndDoF( chi2, 2 * ntrin - 3 );

    // loose convergence criteria if close to end of iterations
    if ( 1. * nbIter > 0.8 * m_Iterations ) maxdz = 10. * m_maxDeltaZ;
    converged = std::abs( deltaz ) < maxdz;

  } // end iteration loop
  if ( !converged ) return StatusCode::FAILURE;

  // set position and covariance
  vtx.setPosition( vtxpos );
  vtx.setCovMatrix( vtxcov );
  // Set tracks. Compute final chi2.
  vtx.clearTracks();
  for ( const auto& trk : pvTracks ) {
    if ( trk.weight() > m_minTrackWeight ) vtx.addToTracks( trk.track(), trk.weight() );
    // remove track for next PV search
    if ( trk.chi2( vtxpos ) < m_trackMaxChi2Remove ) tracks2remove.push_back( trk.track() );
  }
  vtx.setTechnique( LHCb::RecVertex::RecVertexType::Primary );
  return StatusCode::SUCCESS;
}

//=============================================================================
// Get Tukey's weight
//=============================================================================
double AdaptivePV3DFitter::getTukeyWeight( double trchi2, int iter ) const {
  if ( iter < 1 ) return 1.;
  double ctrv = m_trackChi * std::max( m_minIter - iter, 1 );
  auto   sq   = std::pow( ctrv * m_TrackErrorScaleFactor, 2 );
  return sq > trchi2 ? std::pow( 1. - ( trchi2 / sq ), 2 ) : 0.;
}
