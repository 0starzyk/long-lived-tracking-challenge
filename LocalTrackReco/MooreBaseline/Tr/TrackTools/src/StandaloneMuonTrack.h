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
#ifndef STANDALONEMUONTRACK_H
#define STANDALONEMUONTRACK_H 1

// Include files
#include "MuonDAQ/CommonMuonHit.h"
// from Gaudi
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include <string>

/** @class StandaloneMuonTrack StandaloneMuonTrack.h
 *
 *  @author Alessia Satta
 *  @date   2004-10-08
 *
 * Removed from Hlt/HltMuon and ported to Tr/TrackTools
 *
 *  @author Paul Seyfert
 *  @date   2011-03-03
 */
namespace {
  void LinearFitXZ( bool XZ, std::array<CommonMuonHit, 4>& m_points, int nHits, float& slope, float& trunc,
                    float& chi2ndof, float& err_slope, float& err_trunc, float& cov ) {
    float sz2, sz, s0, sxz, sx, sx2;
    sz2 = sz = s0 = sxz = sx = sx2 = 0;
    for ( int i = 0; i < nHits; i++ ) {
      const auto hit  = m_points[i];
      const auto z    = hit.z();
      auto       x    = hit.x();
      auto       xerr = 2. * hit.dx();
      if ( !XZ ) {
        x    = hit.y();
        xerr = 2 * hit.dy();
      }
      sz2 += z * z / xerr / xerr;
      sz += z / xerr / xerr;
      s0 += 1.0 / xerr / xerr;
      sxz += z * x / xerr / xerr;
      sx += x / xerr / xerr;
      sx2 += x * x / xerr / xerr;
    }
    float xdet = sz2 * s0 - sz * sz;
    if ( xdet != 0 ) {
      slope = ( sxz * s0 - sx * sz ) / xdet;
      trunc = ( sx * sz2 - sxz * sz ) / xdet;

      err_trunc = sqrt( sz2 / xdet );
      err_slope = sqrt( s0 / xdet );
      cov       = -sz / xdet;

      chi2ndof = ( sx2 + slope * slope * sz2 + trunc * trunc * s0 - 2. * slope * sxz - 2. * trunc * sx +
                   2 * slope * trunc * sz ) /
                 ( nHits - 2 );
    }
  }
} // namespace

class StandaloneMuonTrack final {
public:
  /// Standard constructor
  StandaloneMuonTrack() { m_clone = 0; };

  virtual ~StandaloneMuonTrack(){}; ///< Destructor
  void setPoint( unsigned int station, CommonMuonHit point ) { m_points[station] = point; };

  CommonMuonHit point( unsigned int station ) const { return m_points[station]; };

  void setClone() { m_clone = 1; };
  bool isClone() const { return m_clone > 0; }

  double slopeX( int stationFirst, int stationSecond, double zFirst, double zSecond ) const {
    return ( m_points[stationFirst].x() - m_points[stationSecond].x() ) / ( zFirst - zSecond );
  };
  double slopeY( int stationFirst, int stationSecond, double zFirst, double zSecond ) const {
    return ( m_points[stationFirst].y() - m_points[stationSecond].y() ) / ( zFirst - zSecond );
  };

  // Fit with a min chi^2 in the 2 projections xz and yz
  bool linearFit() {
    double dof = nHits() - 2.;
    if ( dof < 0 ) return false;
    LinearFitXZ( true, m_points, nHits(), m_sx, m_bx, m_chi2x, m_errsx, m_errbx, m_covbsx );
    LinearFitXZ( false, m_points, nHits(), m_sy, m_by, m_chi2y, m_errsy, m_errby, m_covbsy );
    if ( m_chi2x > -1.f && m_chi2y > -1.f )
      return true;
    else
      return false;
  };

  inline double chi2x() const { return m_chi2x; } /// chi2/dof XZ
  inline double chi2y() const { return m_chi2y; } /// chi2/dof YZ
  // slope XZ
  inline double sx() const { return m_sx; }
  // intercept XZ
  inline double bx() const { return m_bx; }
  // slope YZ
  inline double sy() const { return m_sy; }
  // intercept YZ
  inline double by() const { return m_by; }

  // errors on the above parameters
  inline double errsx() const { return m_errsx; }
  inline double errbx() const { return m_errbx; }
  inline double covbsx() const { return m_covbsx; }
  inline double errsy() const { return m_errsy; }
  inline double errby() const { return m_errby; }
  inline double covbsy() const { return m_covbsy; }

  inline void setnHits( int n ) { m_nHits = n; }
  inline int  nHits() const { return m_nHits; }

private:
  std::array<CommonMuonHit, 4> m_points;
  unsigned int                 m_clone;
  unsigned int                 m_nHits;

  float m_chi2x  = -1.f;
  float m_chi2y  = -1.f;
  float m_sx     = -1.f;
  float m_bx     = -1.f;
  float m_sy     = -1.f;
  float m_by     = -1.f;
  float m_errsx  = -1.f;
  float m_errbx  = -1.f;
  float m_covbsx = -1.f;
  float m_errsy  = -1.f;
  float m_errby  = -1.f;
  float m_covbsy = -1.f;
};
#endif // STANDALONEMUONTRACK_H
