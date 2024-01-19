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

// Include files
#include "GaudiKernel/Point3DTypes.h"

#include "Event/PrHits.h"
#include "Event/State.h"
#include "Event/StateParameters.h"

#include "boost/container/small_vector.hpp"

namespace Downstream {
  struct Hit {
    const LHCb::Pr::UT::Hits* hits;
    int                       hit;
    float                     x, z;
    float                     projection;

    using F = SIMDWrapper::scalar::types::float_v;
    using I = SIMDWrapper::scalar::types::int_v;

    Hit( const LHCb::Pr::UT::Hits* hits, const int hit, float x, float z, float proj )
        : hits( hits ), hit( hit ), x( x ), z( z ), projection( proj ) {}

    [[nodiscard]] auto lhcbID() const {
      const auto myHits = hits->scalar();
      const auto mH     = myHits[hit];
      const auto chanID = mH.get<LHCb::Pr::UT::UTHitsTag::channelID>().cast();
      return bit_cast<int>( LHCb::LHCbID( LHCb::Detector::UT::ChannelID( chanID ) ).lhcbID() );
    }
    [[nodiscard]] int planeCode() const {
      const auto myHits  = hits->scalar();
      const auto mH      = myHits[hit];
      auto       lhcbid  = mH.get<LHCb::Pr::UT::UTHitsTag::channelID>().cast();
      auto       station = ( lhcbid & static_cast<int>( UTInfo::MasksBits::StationMask ) ) >>
                     static_cast<int>( UTInfo::MasksBits::StationBits );
      auto layer = ( lhcbid & static_cast<int>( UTInfo::MasksBits::LayerMask ) ) >>
                   static_cast<int>( UTInfo::MasksBits::LayerBits );
      return 2 * ( station - 1 ) + ( layer - 1 );
    }
    [[nodiscard]] auto weight() const {
      const auto myHits = hits->scalar();
      const auto mH     = myHits[hit];
      return mH.get<LHCb::Pr::UT::UTHitsTag::weight>().cast();
    }

    [[nodiscard]] auto sin() const {
      const auto myHits = hits->scalar();
      const auto mH     = myHits[hit];
      return -mH.get<LHCb::Pr::UT::UTHitsTag::dxDy>().cast() * mH.get<LHCb::Pr::UT::UTHitsTag::cos>().cast();
    }

    [[nodiscard]] auto zAtYEq0() const {
      const auto myHits = hits->scalar();
      const auto mH     = myHits[hit];
      return mH.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>().cast();
    }
    [[nodiscard]] bool isYCompatible( const float y, const float tol ) const {
      const auto myHits = hits->scalar();
      const auto mH     = myHits[hit];
      const auto yMin =
          std::min( mH.get<LHCb::Pr::UT::UTHitsTag::yBegin>().cast(), mH.get<LHCb::Pr::UT::UTHitsTag::yEnd>().cast() );
      const auto yMax =
          std::max( mH.get<LHCb::Pr::UT::UTHitsTag::yBegin>().cast(), mH.get<LHCb::Pr::UT::UTHitsTag::yEnd>().cast() );
      return yMin - tol <= y && y <= yMax + tol;
    }
    [[nodiscard]] auto xAt( const float y ) const {
      const auto myHits = hits->scalar();
      const auto mH     = myHits[hit];
      return mH.get<LHCb::Pr::UT::UTHitsTag::xAtYEq0>().cast() + y * mH.get<LHCb::Pr::UT::UTHitsTag::dxDy>().cast();
    }
  };

  using Hits = std::vector<Hit, LHCb::Allocators::EventLocal<Hit>>;

  inline constexpr auto IncreaseByProj = []( const Hit& lhs, const Hit& rhs ) {
    if ( lhs.projection < rhs.projection ) return true;
    if ( rhs.projection < lhs.projection ) return false;
    return lhs.lhcbID() < rhs.lhcbID();
  };
} // namespace Downstream

/** @class PrDownTrack PrDownTrack.h
 *  Track helper for Downstream track search
 *  Adapted from Pat/PatKShort package
 *  Further adapted for use with PrLongLivedTracking
 *
 *  @author Olivier Callot
 *  @date   2007-10-18
 *
 *  @author Adam Davis
 *  @date   2016-04-10
 *
 *  @author Christoph Hasse (new framework)
 *  @date   2017-03-01
 */

class PrDownTrack final {
public:
  using Hits = boost::container::small_vector<Downstream::Hit, 12, LHCb::Allocators::EventLocal<Downstream::Hit>>;
  // using Hits = boost::container::static_vector<Downstream::Hit, 20>;
  // Until we can put a bound on the number of hits, use a small_vector

  PrDownTrack( Gaudi::TrackVector stateVector, double stateZ, double zUT, LHCb::span<const double, 7> magnetParams,
               LHCb::span<const double> yParams, LHCb::span<const double, 3> momPar, double magnetScale )
      : m_stateVector( stateVector ), m_stateZ( stateZ ), m_zUT( zUT ) {
    const auto tx2  = stateTx() * stateTx();
    const auto ty2  = stateTy() * stateTy();
    m_momentumParam = ( momPar[0] + momPar[1] * tx2 + momPar[2] * ty2 ) * magnetScale;

    // -- See PrFitKsParams to see how these coefficients are derived.
    double zMagnet = magnetParams[0] + magnetParams[1] * ty2 + magnetParams[2] * tx2 +
                     magnetParams[3] * std::abs( stateQoP() ) + /// this is where the old one stopped.
                     magnetParams[4] * std::abs( stateX() ) + magnetParams[5] * std::abs( stateY() ) +
                     magnetParams[6] * std::abs( stateTy() );

    const double dz      = zMagnet - stateZ;
    double       xMagnet = stateX() + dz * stateTx();
    m_slopeX             = xMagnet / zMagnet;
    const double dSlope  = std::abs( m_slopeX - stateTx() );
    const double dSlope2 = dSlope * dSlope;

    double by = stateY() / ( stateZ + ( yParams[0] * fabs( stateTy() ) * zMagnet + yParams[1] ) * dSlope2 );
    m_slopeY  = by * ( 1. + yParams[0] * fabs( by ) * dSlope2 );

    const double yMagnet = stateY() + dz * by - yParams[1] * by * dSlope2;

    // -- These resolutions are semi-empirical and are obtained by fitting residuals
    // -- with MCHits and reconstructed tracks
    // -- See Tracking &Alignment meeting, 19.2.2015, for the idea
    double errXMag = dSlope2 * 15.0 + dSlope * 15.0 + 3.0;
    double errYMag = dSlope2 * 80.0 + dSlope * 10.0 + 4.0;

    // -- Assume better resolution for SciFi than for OT
    // -- obviously this should be properly tuned...
    errXMag /= 2.0;
    errYMag /= 1.5;

    // errXMag = 0.5  + 5.3*dSlope + 6.7*dSlope2;
    // errYMag = 0.37 + 0.7*dSlope - 4.0*dSlope2 + 11*dSlope2*dSlope;

    m_weightXMag = 1.0 / ( errXMag * errXMag );
    m_weightYMag = 1.0 / ( errYMag * errYMag );

    m_magnet = Gaudi::XYZPoint( xMagnet, yMagnet, zMagnet );

    //=== Save for reference
    m_displX = 0.;
    m_displY = 0.;

    //=== Initialize all other data members
    m_chi2 = 0.;
  }

  /// getters
  double      stateX() const { return m_stateVector[0]; }
  double      stateY() const { return m_stateVector[1]; }
  double      stateZ() const { return m_stateZ; }
  double      stateTx() const { return m_stateVector[2]; }
  double      stateTy() const { return m_stateVector[3]; }
  double      stateQoP() const { return m_stateVector[4]; }
  Hits&       hits() { return m_hits; }
  const Hits& hits() const { return m_hits; }
  double      xMagnet() const { return m_magnet.x(); }
  double      yMagnet() const { return m_magnet.y(); }
  double      zMagnet() const { return m_magnet.z(); }
  double      slopeX() const { return m_slopeX; }
  double      slopeY() const { return m_slopeY; }
  double      weightXMag() const { return m_weightXMag; }
  double      weightYMag() const { return m_weightYMag; }
  double      chi2() const { return m_chi2; }

  /// setters
  void setSlopeX( double slopeX ) noexcept { m_slopeX = slopeX; }
  void setChi2( double chi2 ) noexcept { m_chi2 = chi2; }

  // functions
  double xAtZ( double z ) const noexcept {
    const double curvature = 1.6e-5 * ( stateTx() - m_slopeX );
    return xMagnet() + ( z - zMagnet() ) * m_slopeX + curvature * ( z - m_zUT ) * ( z - m_zUT );
  }

  double yAtZ( double z ) const noexcept { return yMagnet() + m_displY + ( z - zMagnet() ) * slopeY(); }

  void updateX( double dx, double dsl ) noexcept {
    m_displX += dx;
    m_magnet = Gaudi::XYZPoint( m_magnet.x() + dx, m_magnet.y(), m_magnet.z() );
    m_slopeX += dsl;
  }
  void updateY( double dy ) noexcept { m_displY += dy; }

  double dxMagnet() const noexcept { return -m_displX; }

  double initialChi2() const noexcept {
    return m_displX * m_displX * m_weightXMag + m_displY * m_displY * m_weightYMag;
  }

  double momentum() const noexcept { return m_momentumParam / ( stateTx() - m_slopeX ); }

  double pt() const noexcept {
    const double tx2      = slopeX() * slopeX();
    const double ty2      = slopeY() * slopeY();
    const double sinTrack = sqrt( 1. - 1. / ( 1. + tx2 + ty2 ) );
    return sinTrack * std::abs( momentum() );
  }

  double distance( const Downstream::Hit& hit ) const noexcept { return hit.x - xAtZ( hit.z ); }

  void sortFinalHits() noexcept {
    std::sort( m_hits.begin(), m_hits.end(), []( const Downstream::Hit& lhs, const Downstream::Hit& rhs ) {
      return std::make_tuple( lhs.z, lhs.lhcbID() ) < std::make_tuple( rhs.z, rhs.lhcbID() );
    } );
  }

private:
  Gaudi::TrackVector m_stateVector;
  double             m_stateZ;
  Gaudi::XYZPoint    m_magnet;

  double m_momentumParam;
  double m_zUT;
  double m_slopeX;
  double m_slopeY;
  double m_displX;
  double m_displY;
  double m_weightXMag;
  double m_weightYMag;
  double m_chi2;

  Hits m_hits; /// working list of hits on this track
};

// -- A typedef for a collection of downstream tracks... From PatDownTrack
typedef std::vector<PrDownTrack> PrDownTracks;
