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
#pragma once

#include <vector>

// LHCb
#include "Event/PrHits.h"
#include "Event/PrVeloTracks.h"
#include <iomanip>

/**
 * Velo only Kalman fit helpers
 */

namespace VeloKalmanParam {
  constexpr float err = 0.0125f;
  constexpr float wx  = err * err;
  constexpr float wy  = wx;

  constexpr float scatterSensorParameters[4] = {0.54772f, 1.478845f, 0.626634f, -0.78f};
  constexpr float scatterFoilParameters[2]   = {1.67f, 20.f};
} // namespace VeloKalmanParam

template <typename T>
class FittedState {
public:
  T x, y, z;
  T tx, ty;
  T covXX, covXTx, covTxTx;
  T covYY, covYTy, covTyTy;

  FittedState() {}

  FittedState( Vec3<T> pos, Vec3<T> dir, T covXX, T covXTx, T covTxTx, T covYY, T covYTy, T covTyTy )
      : x( pos.x )
      , y( pos.y )
      , z( pos.z )
      , tx( dir.x )
      , ty( dir.y )
      , covXX( covXX )
      , covXTx( covXTx )
      , covTxTx( covTxTx )
      , covYY( covYY )
      , covYTy( covYTy )
      , covTyTy( covTyTy ) {}

  inline Vec3<T> pos() const { return Vec3<T>( x, y, z ); }
  inline Vec3<T> dir() const { return Vec3<T>( tx, ty, 1.f ); }
  inline Vec3<T> covX() const { return Vec3<T>( covXX, covXTx, covTxTx ); }
  inline Vec3<T> covY() const { return Vec3<T>( covYY, covYTy, covTyTy ); }

  inline T zBeam() const {
    const T x0    = x - z * tx;
    const T y0    = y - z * ty;
    T       denom = tx * tx + ty * ty;
    return select( denom < 0.001f * 0.001f, z, -( x0 * tx + y0 * ty ) / denom );
  }

  inline void transportTo( const T& toZ ) {
    const T dz  = toZ - z;
    const T dz2 = dz * dz;

    x = x + dz * tx;
    y = y + dz * ty;
    z = toZ;

    covXX  = covXX + dz2 * covTxTx + 2.f * dz * covXTx;
    covXTx = covXTx + dz * covTxTx;
    covYY  = covYY + dz2 * covTyTy + 2.f * dz * covYTy;
    covYTy = covYTy + dz * covTyTy;
  }
};

template <typename M, typename F>
inline __attribute__( ( always_inline ) ) void filter( const M mask, const F z, F& x, F& tx, F& covXX, F& covXTx,
                                                       F& covTxTx, const F zhit, const F xhit, const F winv ) {
  // compute prediction
  const F dz    = zhit - z;
  const F predx = x + dz * tx;

  const F dz_t_covTxTx = dz * covTxTx;
  const F predcovXTx   = covXTx + dz_t_covTxTx;
  const F dz_t_covXTx  = dz * covXTx;

  const F predcovXX = covXX + 2.f * dz_t_covXTx + dz * dz_t_covTxTx;

  // compute the gain matrix
  const F R   = 1.0f / ( winv + predcovXX );
  const F Kx  = predcovXX * R;
  const F KTx = predcovXTx * R;

  // update the state vector
  const F r = xhit - predx;
  x         = select( mask, predx + Kx * r, x );
  tx        = select( mask, tx + KTx * r, tx );

  // update the covariance matrix
  covTxTx = select( mask, R * ( covTxTx * ( winv + covXX ) - covXTx * covXTx ), covTxTx );
  covXTx  = select( mask, winv * KTx, covXTx );
  covXX   = select( mask, winv * Kx, covXX );
}

template <typename F, typename I, typename M>
inline __attribute__( ( always_inline ) ) FittedState<F>
fitBackward( const M track_mask, I& nHits, const LHCb::Pr::VP::Hits& hits, Vec3<F>& dir,
             std::array<I, LHCb::Pr::Velo::Tracks::MaxVPHits>& vp_index ) {
  int        maxHits    = nHits.hmax( track_mask );
  I          idxHit0    = vp_index[0];
  const auto hitproxy   = hits.simd();
  auto       hit_gather = hitproxy.gather( idxHit0, track_mask );
  Vec3<F>    pos        = {select( track_mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().x(), 0.f ),
                 select( track_mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().y(), 0.f ),
                 select( track_mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().z(), 0.f )};

  FittedState<F> s =
      FittedState<F>( pos, dir, F( VeloKalmanParam::wx ), 0.f, 0.01f, F( VeloKalmanParam::wy ), 0.f, 0.01f );

  // Parameters for kalmanfit scattering. calibrated on MC, shamelessly hardcoded:
  const F noise2PerLayer = 1e-8f + 7e-6f * ( s.tx * s.tx + s.ty * s.ty );

  for ( int i = 1; i < maxHits; i++ ) {
    auto    mask       = track_mask && ( I( i ) < nHits );
    I       idxHit     = vp_index[i];
    auto    hit_gather = hitproxy.gather( idxHit, mask );
    Vec3<F> hit        = {select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().x(), 0.f ),
                   select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().y(), 0.f ),
                   select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().z(), 0.f )};

    s.covTxTx = select( mask, s.covTxTx + noise2PerLayer, s.covTxTx );
    s.covTyTy = select( mask, s.covTyTy + noise2PerLayer, s.covTyTy );

    filter( mask, s.z, s.x, s.tx, s.covXX, s.covXTx, s.covTxTx, hit.z, hit.x, F( VeloKalmanParam::wx ) );
    filter( mask, s.z, s.y, s.ty, s.covYY, s.covYTy, s.covTyTy, hit.z, hit.y, F( VeloKalmanParam::wy ) );
    s.z = select( mask, hit.z, s.z );
  }

  s.covTxTx = s.covTxTx + noise2PerLayer;
  s.covTyTy = s.covTyTy + noise2PerLayer;

  return s;
}

template <typename F, typename I, typename M>
inline __attribute__( ( always_inline ) ) FittedState<F>
fitForward( const M track_mask, I& nHits, const LHCb::Pr::VP::Hits& hits, Vec3<F>& dir,
            std::array<I, LHCb::Pr::Velo::Tracks::MaxVPHits>& vp_index ) {
  int        maxHits    = nHits.hmax( track_mask );
  auto       mask       = track_mask && I( maxHits - 1 ) < nHits;
  I          idxHit0    = vp_index[maxHits - 1];
  const auto hitproxy   = hits.simd();
  auto       hit_gather = hitproxy.gather( idxHit0, mask );
  Vec3<F>    pos        = {select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().x(), 0.f ),
                 select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().y(), 0.f ),
                 select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().z(), 0.f )};

  FittedState<F> s =
      FittedState<F>( pos, dir, F( VeloKalmanParam::wx ), 0.f, 0.01f, F( VeloKalmanParam::wy ), 0.f, 0.01f );

  // Parameters for kalmanfit scattering. calibrated on MC, shamelessly hardcoded:
  const F noise2PerLayer = 1e-8f + 7e-6f * ( s.tx * s.tx + s.ty * s.ty );

  for ( int i = maxHits - 2; i >= 0; i-- ) {
    auto    mask       = track_mask && ( I( i ) < nHits );
    I       idxHit     = vp_index[i];
    auto    hit_gather = hitproxy.gather( idxHit, mask );
    Vec3<F> hit        = {select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().x(), 0.f ),
                   select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().y(), 0.f ),
                   select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().z(), 0.f )};

    s.covTxTx = select( mask, s.covTxTx + noise2PerLayer, s.covTxTx );
    s.covTyTy = select( mask, s.covTyTy + noise2PerLayer, s.covTyTy );

    filter( mask, s.z, s.x, s.tx, s.covXX, s.covXTx, s.covTxTx, hit.z, hit.x, F( VeloKalmanParam::wx ) );
    filter( mask, s.z, s.y, s.ty, s.covYY, s.covYTy, s.covTyTy, hit.z, hit.y, F( VeloKalmanParam::wy ) );
    s.z = select( mask, hit.z, s.z );
  }

  s.covTxTx = s.covTxTx + noise2PerLayer;
  s.covTyTy = s.covTyTy + noise2PerLayer;

  return s;
}

template <typename M, typename F>
inline F filterWithMomentum( const M mask, const F z, F& x, F& tx, F& covXX, F& covXTx, F& covTxTx, const F zhit,
                             const F xhit, const F winv, const F qop ) {
  // compute prediction
  const F dz    = zhit - z;
  const F predx = x + dz * tx;

  const F dz_t_covTxTx = dz * covTxTx;
  const F dz_t_covXTx  = dz * covXTx;

  // Add noise
  const F par1 = VeloKalmanParam::scatterSensorParameters[0];
  const F par2 = VeloKalmanParam::scatterSensorParameters[1];
  const F par6 = VeloKalmanParam::scatterSensorParameters[2];
  const F par7 = VeloKalmanParam::scatterSensorParameters[3];

  const F sigTx = par1 * 1e-5f + par2 * abs( qop );
  const F sigX  = par6 * sigTx * abs( dz );
  const F corr  = par7;

  const F eXX   = sigX * sigX;
  const F eXTx  = corr * sigX * sigTx;
  const F eTxTx = sigTx * sigTx;

  const F predcovXX  = covXX + 2.f * dz_t_covXTx + dz * dz_t_covTxTx + eXX;
  const F predcovXTx = covXTx + dz_t_covTxTx + eXTx;
  // compute the gain matrix
  const F R   = 1.0f / ( winv + predcovXX );
  const F Kx  = predcovXX * R;
  const F KTx = predcovXTx * R;

  // update the state vector
  const F r = xhit - predx;
  x         = select( mask, predx + Kx * r, x );
  tx        = select( mask, tx + KTx * r, tx );

  // update the covariance matrix
  /*
    Linearisation of the expression to avoid absorbtion:

    covTxTx = predcovTxTx - KTx * predcovXTx
    covTxTx = predcovTxTx - predcovXTx^2 / ( winv + predcovXX )
    covTxTx = eTxTx + (covTxTx * ( winv + predcovXX ) - predcovXTx^2) / ( winv + predcovXX )
    covTxTx = eTxTx + (covTxTx * ( winv + predcovXX ) - predcovXTx^2) / ( winv + predcovXX )
    ((((((
    predcovXTx^2 = (covXTx + dz*covTxTx + eXTx)^2
                = covXTx^2 + (dz*covTxTx)^2 + eXTx^2 + 2*covXTx*dz*covTxTx + 2*covXTx*eXTx + 2*dz*covTxTx*eXTx
    covTxTx * ( winv + predcovXX ) = covTxTx * ( winv + covXX + 2*dz*covXTx + dz^2*covTxTx + eXX )
                                   = covTxTx * ( winv + covXX) + 2*dz*covXTx*covTxTx + (dz*covTxTx)^2 + eXX*covTxTx
    ))))))
    covTxTx = eTxTx + (covTxTx * ( winv + covXX) - covXTx^2 + eXX*covTxTx - eXTx*(eXTx + 2*(covXTx + dz*covTxTx))) / (
    winv + predcovXX )
   */
  covTxTx = select( mask,
                    eTxTx + R * ( covTxTx * ( winv + covXX ) - covXTx * covXTx + eXX * covTxTx -
                                  eXTx * ( eXTx + 2.f * ( covXTx + dz_t_covTxTx ) ) ),
                    covTxTx );
  covXTx  = select( mask, winv * KTx, covXTx );
  covXX   = select( mask, winv * Kx, covXX );
  // return the chi2
  return r * r * R;
}

template <typename F, typename I, typename M>
inline __attribute__( ( always_inline ) ) std::tuple<FittedState<F>, F, I>
fitBackwardWithMomentum( const M track_mask, const LHCb::Pr::Velo::Tracks& tracks, const I idxVP, const F qop,
                         const LHCb::Pr::VP::Hits& hits, const int state_id ) {
  auto const velotracks = tracks.simd();
  auto const track      = velotracks.gather( idxVP, track_mask );
  I          nHits      = track.nHits();
  int        maxHits    = nHits.hmax( track_mask );
  I          idxHit0    = track.vp_index( 0 );
  Vec3<F>    dir        = track.StateDir( state_id );
  auto       hitproxy   = hits.simd();
  auto       hit_gather = hitproxy.gather( idxHit0, track_mask );
  Vec3<F>    pos        = {select( track_mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().x(), 0.f ),
                 select( track_mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().y(), 0.f ),
                 select( track_mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().z(), 0.f )};

  FittedState<F> s =
      FittedState<F>( pos, dir, F( VeloKalmanParam::wx ), 0.f, 0.01f, F( VeloKalmanParam::wy ), 0.f, 0.01f );

  F chi2 = 0.f;

  for ( int i = 1; i < maxHits; i++ ) {
    auto    mask       = track_mask && ( I( i ) < nHits );
    I       idxHit     = track.vp_index( i );
    auto    hit_gather = hitproxy.gather( idxHit, mask );
    Vec3<F> hit        = {select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().x(), 0.f ),
                   select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().y(), 0.f ),
                   select( mask, hit_gather.template get<LHCb::Pr::VP::VPHitsTag::pos>().z(), 0.f )};

    chi2 = select( mask,
                   chi2 + filterWithMomentum( mask, s.z, s.x, s.tx, s.covXX, s.covXTx, s.covTxTx, hit.z, hit.x,
                                              F( VeloKalmanParam::wx ), qop ),
                   chi2 );
    chi2 = select( mask,
                   chi2 + filterWithMomentum( mask, s.z, s.y, s.ty, s.covYY, s.covYTy, s.covTyTy, hit.z, hit.y,
                                              F( VeloKalmanParam::wy ), qop ),
                   chi2 );
    s.z  = select( mask, hit.z, s.z );
  }

  // Convert state at first measurement to state at closest to beam
  const F t2 = s.dir().rho();

  const F scat2RFFoil =
      VeloKalmanParam::scatterFoilParameters[0] * ( 1.0 + VeloKalmanParam::scatterFoilParameters[1] * t2 ) * qop * qop;
  s.covTxTx = s.covTxTx + scat2RFFoil;
  s.covTyTy = s.covTyTy + scat2RFFoil;

  s.transportTo( s.zBeam() );

  return {s, chi2, 2 * nHits - 4};
}
