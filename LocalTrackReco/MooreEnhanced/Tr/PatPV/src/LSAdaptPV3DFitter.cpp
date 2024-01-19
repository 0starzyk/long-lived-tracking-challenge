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
#include "PVUtils.h"

class LSAdaptPV3DFitter : public extends<GaudiTool, IPVFitter> {

public:
  // Standard constructor
  using extends::extends;
  // Fitting
  StatusCode fitVertex( const Gaudi::XYZPoint& seedPoint, LHCb::span<const LHCb::Track* const> tracks,
                        LHCb::RecVertex& vtx, std::vector<const LHCb::Track*>& tracks2remove,
                        IGeometryInfo const& geometry ) const override;

private:
  Gaudi::Property<double> m_maxDeltaZ{this, "maxDeltaZ", 0.0005 * Gaudi::Units::mm, "Fit convergence condition"};
  Gaudi::Property<double> m_minTrackWeight{this, "minTrackWeight", 0.00000001,
                                           "Minimum Tukey's weight to accept a track"}; // 0.00001
  Gaudi::Property<double> m_TrackErrorScaleFactor{this, "TrackErrorScaleFactor", 1.0};
  Gaudi::Property<double> m_trackMaxChi2Remove{this, "trackMaxChi2Remove", 25.,
                                               "Max chi2 tracks to be removed from next PV search"};
  Gaudi::Property<double> m_maxChi2{this, "maxChi2", 400.0,
                                    "maximal chi2 to accept track"}; // Chi2 of completely wrong tracks;       //
  Gaudi::Property<double> m_zVtxShift{this, "zVtxShift", 0.0};
  Gaudi::Property<int>    m_minTr{this, "MinTracks", 5, "Minimum number of tracks to make a vertex"};
  Gaudi::Property<int>    m_Iterations{this, "Iterations", 20, "Number of iterations for minimisation"};
  Gaudi::Property<int>    m_minIter{this, "minIter", 5,
                                 "Minimum number of iterations"}; // iterate at least m_minIter times
  Gaudi::Property<bool>   m_AddMultipleScattering{this, "AddMultipleScattering", true,
                                                "add multiple scattering calculation to not fitted tracks"};

  Gaudi::Property<bool>   m_CalculateMultipleScattering{this, "CalculateMultipleScattering", true,
                                                      "calculate multiple scattering"};
  Gaudi::Property<bool>   m_UseFittedTracks{this, "UseFittedTracks", false, "use tracks fitted by kalman fit"};
  double                  m_scatCons = 0.; // calculated from m_x0MS and 3 GeV
  double                  m_scatConsNoMom; // calculated from m_x0MS to be divided by momemntum [MeV]
  Gaudi::Property<double> m_x0MS{this, "x0MS", 0.02,
                                 [=]( auto& ) {
                                   double X0             = this->m_x0MS;
                                   this->m_scatConsNoMom = ( 13.6 * std::sqrt( X0 ) * ( 1. + 0.038 * log( X0 ) ) );
                                   this->m_scatCons      = this->m_scatConsNoMom / ( 3.0 * Gaudi::Units::GeV );
                                 },
                                 Gaudi::Details::Property::ImmediatelyInvokeHandler{true}}; // X0 (tunable) of MS to add
                                                                                            // for extrapolation of
                                                                                            // track parameters to PV
  double                  m_trackChi = 3.;                                                  // sqrt of m_trackMaxChi2
  Gaudi::Property<double> m_trackMaxChi2{this,
                                         "trackMaxChi2",
                                         9.,
                                         [=]( auto& ) { this->m_trackChi = std::sqrt( m_trackMaxChi2 ); },
                                         Gaudi::Details::Property::ImmediatelyInvokeHandler{true},
                                         "maximum chi2 track to accept track in PV fit"};

  // Add track for PV
  void addTrackForPV( const LHCb::Track* str, PVTracks& pvTracks, const Gaudi::XYZPoint& seed ) const;

  double err2d0( const LHCb::Track* track, const Gaudi::XYZPoint& seed ) const;
  // Get Tukey's weight
  double getTukeyWeight( double trchi2, int iter ) const {
    if ( iter < 1 ) return 1.;
    auto ctrv = m_trackChi * std::max( m_minIter - iter, 1 );
    auto cT2  = trchi2 / std::pow( ctrv * m_TrackErrorScaleFactor, 2 );
    return cT2 < 1 ? std::pow( 1 - cT2, 2 ) : 0.;
  }
};
namespace {
  constexpr const double s_myZero = 1E-12; // myzero=1E-12 small number

  Gaudi::XYZVector impactParameterVector( const Gaudi::XYZPoint& vertex, const Gaudi::XYZPoint& point,
                                          const Gaudi::XYZVector& udir ) {
    auto distance = point - vertex;
    // return udir.Cross(distance.Cross(udir));
    // Ax(BxC) = (A.C)B - (B.C)A
    return distance - udir * ( distance.Dot( udir ) );
  }
} // namespace
DECLARE_COMPONENT( LSAdaptPV3DFitter )

//=============================================================================
// Least square adaptive fitting method
//=============================================================================
StatusCode LSAdaptPV3DFitter::fitVertex( const Gaudi::XYZPoint& seedPoint, LHCb::span<const LHCb::Track* const> rTracks,
                                         LHCb::RecVertex& vtx, std::vector<const LHCb::Track*>& tracks2remove,
                                         IGeometryInfo const& ) const {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "================Test==================" << endmsg;

  tracks2remove.clear();

  Gaudi::XYZPoint xyzvtx = seedPoint;
  // prepare tracks
  PVTracks pvTracks;
  pvTracks.reserve( rTracks.size() );
  for ( const auto& track : rTracks ) {
    if ( !( track->hasVelo() ) ) continue;
    addTrackForPV( track, pvTracks, seedPoint );
  }

  if ( (int)pvTracks.size() < m_minTr ) {
    if ( msgLevel( MSG::DEBUG ) ) { verbose() << "Too few tracks to fit PV" << endmsg; }
    return StatusCode::FAILURE;
  }

  Gaudi::SymMatrix3x3            hess;
  ROOT::Math::SVector<double, 3> d0vec;

  bool converged = false;
  int  nbIter    = 0;
  while ( ( nbIter < m_minIter ) || ( !converged && nbIter < m_Iterations ) ) {
    ++nbIter;

    std::fill( d0vec.begin(), d0vec.end(), 0. );
    std::fill( hess.begin(), hess.end(), 0. );

    int ntrin = 0;
    for ( auto& pvTrack : pvTracks ) {
      std::array<double, 3> unitVect;
      pvTrack.unitVect.GetCoordinates( unitVect.data() );

      pvTrack.vd0  = impactParameterVector( xyzvtx, pvTrack.stateG.position(), pvTrack.unitVect );
      pvTrack.d0sq = pvTrack.vd0.Mag2();

      if ( pvTrack.err2d0 < s_myZero ) {
        if ( msgLevel( MSG::DEBUG ) ) debug() << "fitVertex: pvTrack error is too small for computation" << endmsg;
        pvTrack.chi2 = ( pvTrack.d0sq ) / s_myZero;
      } else
        pvTrack.chi2 = ( pvTrack.d0sq ) / pvTrack.err2d0;

      if ( msgLevel( MSG::DEBUG ) )
        debug() << format( "itr %d new track position: chi2 d0 w zc %7.2f %7.2f %7.2f %7.2f", nbIter, pvTrack.chi2,
                           std::sqrt( pvTrack.d0sq ), pvTrack.err2d0, pvTrack.stateG.z() )
                << endmsg;

      pvTrack.weight = getTukeyWeight( pvTrack.chi2, nbIter );

      if ( pvTrack.weight > m_minTrackWeight ) {
        ntrin++;

        double invs = 1.0;
        if ( pvTrack.err2d0 < s_myZero ) {
          invs = ( 2.0 / s_myZero ) * pvTrack.weight;
        } else
          invs = ( 2.0 / pvTrack.err2d0 ) * pvTrack.weight;

        std::array<double, 3> ipVec;
        pvTrack.vd0.GetCoordinates( ipVec.data() );
        for ( int k = 0; k < 3; ++k ) {
          for ( int l = 0; l < 3; ++l ) {
            double delta_kl = ( k == l ) ? 1.0 : 0.0;
            double val      = ( delta_kl - unitVect[k] * unitVect[l] ) * invs;
            d0vec[k] += ipVec[l] * val;
            if ( l <= k ) hess( k, l ) += val;
          }
        }
      }
    }

    // check nr of tracks that entered the fit
    if ( ntrin < 2 ) return StatusCode::FAILURE;

    if ( !hess.InvertChol() ) {
      if ( msgLevel( MSG::DEBUG ) ) {
        debug() << "Error inverting hessian matrix" << endmsg;
        // typically 2 tracks (clones) left with identical slopes
        // dump all tracks used in hessian calculation
        int dntrin = 0;
        for ( const auto& fitTr : pvTracks ) {
          if ( fitTr.weight > m_minTrackWeight ) {
            ++dntrin;
            info() << format( "chi2 d0 w %12.6e %12.6e %12.6e xyz %12.6e %12.6e %12.6e tx ty %12.6e %12.6e", fitTr.chi2,
                              std::sqrt( fitTr.d0sq ), fitTr.weight, fitTr.refTrack->firstState().x(),
                              fitTr.refTrack->firstState().y(), fitTr.refTrack->firstState().z(),
                              fitTr.refTrack->firstState().tx(), fitTr.refTrack->firstState().ty() )
                   << endmsg;
          }
        }
        debug() << "nr of track in hessian failure: " << dntrin++ << endmsg;
      }
      return StatusCode::FAILURE;
    }
    ROOT::Math::SVector<double, 3> delta = hess * d0vec;

    if ( msgLevel( MSG::DEBUG ) ) { debug() << " delta " << delta[0] << " " << delta[1] << " " << delta[2] << endmsg; }

    xyzvtx += Gaudi::XYZVector{delta[0], delta[1], delta[2]};

    // loose convergence criteria if close to end of iterations
    auto maxdz = ( 1. * nbIter > 0.8 * m_Iterations ? 10. * m_maxDeltaZ : m_maxDeltaZ.value() );
    if ( fabs( delta[2] ) < maxdz ) converged = true;

  } // end iteration loop
  if ( !converged ) return StatusCode::FAILURE;
  int    outTracks = 0;
  int    nDoF      = -3;
  double chi2      = 0.0;
  for ( const auto& fitTr : pvTracks ) {
    if ( fitTr.weight > m_minTrackWeight ) {
      ++outTracks;
      chi2 += fitTr.chi2 * fitTr.weight;
      nDoF += 2;
    }
  }

  if ( outTracks < m_minTr ) {
    if ( msgLevel( MSG::DEBUG ) ) { debug() << "Too few tracks after PV fit" << endmsg; }
    return StatusCode::FAILURE;
  }

  // accepted PV
  vtx.setChi2( chi2 );
  vtx.setNDoF( nDoF );

  xyzvtx.SetZ( xyzvtx.z() + m_zVtxShift );
  vtx.setPosition( xyzvtx );
  hess *= 2;
  vtx.setCovMatrix( hess );
  // Set tracks
  vtx.clearTracks();
  for ( const auto& fitTr : pvTracks ) {
    if ( fitTr.weight > m_minTrackWeight ) { vtx.addToTracks( fitTr.refTrack, fitTr.weight ); }
  }
  vtx.setTechnique( LHCb::RecVertex::RecVertexType::Primary );

  // fill tracks to remove from next PV search
  for ( const auto& fitTr : pvTracks ) {
    if ( fitTr.chi2 < m_trackMaxChi2Remove ) { tracks2remove.push_back( fitTr.refTrack ); }
  }

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Vertex" << endmsg;
    debug() << "===================" << endmsg;
    debug() << format( "chi2/ndof %7.2f  %7.2f / %5d xyz %7.3f %7.3f %7.3f  err xyz: %7.4f %7.4f %7.4f",
                       vtx.chi2() / vtx.nDoF(), vtx.chi2(), vtx.nDoF(), vtx.position().X(), vtx.position().Y(),
                       vtx.position().Z(), std::sqrt( vtx.covMatrix()( 0, 0 ) ), std::sqrt( vtx.covMatrix()( 1, 1 ) ),
                       std::sqrt( vtx.covMatrix()( 2, 2 ) ) )
            << endmsg << endmsg;
    debug() << "Tracks in this vertex" << endmsg;
    debug() << "---------------------" << endmsg;
    for ( const auto& fitTr : pvTracks ) {
      int invx = 0;
      if ( fitTr.weight > m_minTrackWeight ) invx = 1;
      debug() << format( "chi2 d0 w zc %7.2f %7.2f %7.2f %7.2f %7.2f %3d", fitTr.chi2, std::sqrt( fitTr.d0sq ),
                         fitTr.weight, fitTr.refTrack->firstState().z(), fitTr.refTrack->pt(), invx )
              << endmsg;
    }
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
// Add track for PV
//=============================================================================
void LSAdaptPV3DFitter::addTrackForPV( const LHCb::Track* pvtr, PVTracks& pvTracks,
                                       const Gaudi::XYZPoint& seed ) const {
  // Add new PVTrack
  PVTrack pvtrack;
  pvtrack.isUsed   = false;
  pvtrack.stateG   = pvtr->firstState();
  pvtrack.unitVect = pvtrack.stateG.slopes().Unit();
  pvtrack.vd0      = impactParameterVector( seed, pvtrack.stateG.position(), pvtrack.unitVect );
  pvtrack.d0sq     = pvtrack.vd0.Mag2();

  if ( m_UseFittedTracks ) {
    Gaudi::XYZVector               vd0Unit = pvtrack.vd0 / std::sqrt( pvtrack.d0sq );
    ROOT::Math::SVector<double, 2> xyVec( vd0Unit.x(), vd0Unit.y() );
    pvtrack.err2d0 = Similarity( pvtrack.stateG.covariance().Sub<Gaudi::SymMatrix2x2>( 0, 0 ), xyVec );
  } else {
    pvtrack.err2d0 = err2d0( pvtr, seed );
  }
  pvtrack.chi2 = 0;
  if ( pvtrack.err2d0 > s_myZero ) {
    pvtrack.chi2 = ( pvtrack.d0sq ) / pvtrack.err2d0;
  } else {
    if ( msgLevel( MSG::DEBUG ) ) debug() << "addTrackForPV: pvTrack error is too small for computation" << endmsg;
    pvtrack.chi2 = std::sqrt( pvtrack.d0sq ) / s_myZero;
  }
  // Keep reference to the original track
  pvtrack.refTrack = pvtr;
  if ( pvtrack.chi2 < m_maxChi2 ) pvTracks.push_back( pvtrack );
}

//=============================================================================
//  err2d0
//=============================================================================
double LSAdaptPV3DFitter::err2d0( const LHCb::Track* track, const Gaudi::XYZPoint& seed ) const {

  // fast parametrization of track parameters
  double z  = track->firstState().z();
  double tx = track->firstState().tx();
  double ty = track->firstState().ty();

  double ex2 = track->firstState().errX2();
  double ey2 = track->firstState().errY2();

  double tr2   = tx * tx + ty * ty;
  double fcos2 = 1. / ( 1. + tr2 );
  double fsin2 = tr2 / ( 1. + tr2 );
  double err2  = ( ex2 + ey2 ) * fcos2;

  bool backward = track->isVeloBackward();

  double z_seed = seed.z();
  double dz_pv  = z - z_seed;
  if ( backward ) dz_pv *= -1.;
  double et2 = 0.;
  if ( dz_pv > 0. ) {
    double etx2 = track->firstState().errTx2();
    double ety2 = track->firstState().errTy2();
    et2         = dz_pv * dz_pv * ( etx2 + ety2 );
    err2 += et2;
  }

  // add multiple scattering if track errors don't have it (HLT case)
  double corr2   = 0.;
  double l_scat2 = 0.;

  if ( m_AddMultipleScattering ) {
    if ( m_CalculateMultipleScattering ) {
      // mutliple scattering from RF foil at r = 10 mm

      double r_ms = 10.;
      l_scat2     = r_ms * r_ms / fsin2;
      corr2       = m_scatCons * m_scatCons * l_scat2;

      if ( msgLevel( MSG::DEBUG ) ) {
        debug() << " Track printout " << track->type() << " " << track->firstState().qOverP() << " " << endmsg;
      }
      err2 += corr2;
    } else {
      err2 += 0.05 * 0.05;
    }

    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << " Track printout " << track->type() << " " << track->firstState().qOverP() << " " << corr2 << " "
              << err2 << " " << endmsg;
    }
  }
  return err2;
}
