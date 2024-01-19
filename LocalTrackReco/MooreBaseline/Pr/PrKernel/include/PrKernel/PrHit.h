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
#ifndef PRKERNEL_PRHIT_H
#define PRKERNEL_PRHIT_H 1

// Include files
#include "Kernel/EventLocalAllocator.h"
#include "Kernel/LHCbID.h"

/** @class PrHit PrHit.h PrKernel/PrHit.h
 *  Hits to be used in the pattern in the T/TT stations
 *
 *  @author Olivier Callot
 *  @date   2012-03-13
 *  @author Thomas Nikodem
 *  @date   2016-04-11
 */
class PrHit final {
public:
  /// Standard constructor
  PrHit( const LHCb::LHCbID id, const float x0, const float z0, const float dxDy, const float dzDy, const float yMin,
         const float yMax, const float werr, const float w,
         const uint info )
      : m_x0( x0 )
      , ///< x coordinate at y = 0
      m_info( info )
      , m_yMin( yMin )
      , m_yMax( yMax )
      , m_werr( werr )
      , m_w( w )
      , m_z0( z0 )
      , ///< z coordinate at y = 0
      m_dxDy( dxDy )
      , ///< Slope x vs y, typically 0 for x layers
      m_dzDy( dzDy )
      , ///< Slope z vs y, as detectors are vertical while Z axis is not horizontal
      m_id( id ){};

  /// Standard constructor
  PrHit( const LHCb::LHCbID id, const float x0, const float z0, const float dxDy, const float dzDy, const float yMin,
         const float yMax, const float errX, const int zone,
         const uint planeCode )
      : m_x0( x0 )
      , ///< x coordinate at y = 0
      m_info( ( zone & 1u ) | ( ( planeCode & 31u ) << 1 ) | ( ( std::abs( dxDy ) < 0.001f ) << 7 ) )
      , m_yMin( yMin )
      , m_yMax( yMax )
      , m_werr( 1.f / errX )
      , m_w( m_werr * m_werr )
      , m_z0( z0 )
      , ///< z coordinate at y = 0
      m_dxDy( dxDy )
      , ///< Slope x vs y, typically 0 for x layers
      m_dzDy( dzDy )
      , ///< Slope z vs y, as detectors are vertical while Z axis is not horizontal
      m_id( id ){};

  void setHit( const LHCb::LHCbID id, const float x0, const float z0, const float dxDy, const float dzDy,
               const float yMin, const float yMax, const float errX, const int zone, const int planeCode ) {
    *this = PrHit( id, x0, z0, dxDy, dzDy, yMin, yMax, errX, zone, planeCode );
  }

  LHCb::LHCbID id() const { return m_id; }
  float        x() const { return m_x0; }
  float        x( float y ) const { return m_x0 + y * m_dxDy; }
  float        z() const { return m_z0; }
  float        z( float y ) const { return m_z0 + y * m_dzDy; }
  float        werr() const { return m_werr; }
  float        w() const { return m_w; }
  float        yMin() const { return m_yMin; }
  float        yMax() const { return m_yMax; }
  float        yOnTrack( float y0, float dyDz ) const { return ( y0 + dyDz * m_z0 ) / ( 1.f - dyDz * m_dzDy ); }
  float        dxDy() const { return m_dxDy; } // used for deltaY
  float        dzDy() const { return m_dzDy; }

  float distance( const float x_track, const float y_track ) const { return x( y_track ) - x_track; }
  float distanceXHit( const float x_track ) const { return m_x0 - x_track; }

  int  planeCode() const { return ( m_info & 63u ) >> 1; }
  int  zone() const { return ( m_info & 1u ); } // only needed in printHit...
  bool isX() const { return ( m_info & 128u ) >> 7; }

  struct LowerByX0 {
    bool operator()( const PrHit& lhs, const PrHit& rhs ) const { return lhs.m_x0 < rhs.m_x0; }
    bool operator()( const PrHit* lhs, const PrHit* rhs ) const { return lhs->x() < rhs->x(); }
  };
  struct LowerByZ {
    bool operator()( const PrHit* lhs, const PrHit* rhs ) const { return lhs->z() < rhs->z(); }
    bool operator()( const PrHit& lhs, const PrHit& rhs ) const { return lhs.z() < rhs.z(); }
  };

private:
  // sort according to access
  float m_x0;   /// x coordinate at y = 0
  uint  m_info; /// several infos are stored in this variable, for definition look at code above
  float m_yMin; /// minimum y coordinate along this segment
  float m_yMax; /// maximum y coordinate along this segment
  float m_werr; /// Add for line fitter vectorised
  float m_w;    /// 1/error^2 of hit

  // Detector segment
  float        m_z0;   ///  z coordinate at y = 0
  float        m_dxDy; ///  Slope x vs y, typically 0 for x layers
  float        m_dzDy; ///  Slope z vs y, as detectors are vertical while Z axis is not horizontal
  LHCb::LHCbID m_id;
};

using PrHits = std::vector<const PrHit*>;

// struct to capsul modifiable information
struct ModPrHit final {
  ModPrHit() = default;
  ModPrHit( float c, size_t i ) : coord{c}, fullDex{i} {}
  float  coord   = 0.f;
  size_t fullDex = 0;
  bool   isValid() const {
    return coord != -std::numeric_limits<float>::max();
  }; //---LoH: crucial for the logic of the HybridSeeding that it is a large negative
  void setInvalid() {
    coord = -std::numeric_limits<float>::max();
  }; //---LoH: crucial for the logic of the HybridSeeding that it is a large negative
  friend bool operator==( const ModPrHit& lhs, const ModPrHit& rhs ) {
    return lhs.fullDex == rhs.fullDex;
  } // FIXME: this is _not_ equality, but (at best) equivalence... one can have ( a==b &&  a<b ) evaluate to true --
    // which is confusing!
  friend bool operator<( const ModPrHit& lhs, const ModPrHit& rhs ) { return lhs.coord < rhs.coord; }
};
using ModPrHits            = std::vector<ModPrHit, LHCb::Allocators::EventLocal<ModPrHit>>;
using ModPrHitIter         = ModPrHits::iterator;
using ModPrHitConstIter    = ModPrHits::const_iterator;
using ModPrHitConstRevIter = ModPrHits::const_reverse_iterator;

#endif // PRKERNEL_PRHIT_H
