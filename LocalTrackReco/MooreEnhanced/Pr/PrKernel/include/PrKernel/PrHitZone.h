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
#include "Kernel/DetectorSegment.h"
#include "PrKernel/PrHit.h"

/** @class PrHitZone PrHitZone.h PrKernel/PrHitZone.h
 *  Store the information of a zone in the T stations
 *  A zone is a part of a layer with boundaries, so that all hits are measuring
 *  one coordinate in this zone.
 *  For FT, this is the top or bottom parts of each layer.
 *  The zone has a number, a plane code, slopes (dxDy, dzDy) and boundaries
 *
 *  @author Olivier Callot
 *  @date   2012-03-13
 */
class PrHitZone final {
public:
  enum class Side { Lower, Upper };

  PrHitZone() = default;
  // All geometry informations casted here!
  void setZone( unsigned int number, DetectorSegment& seg, float xMin, float xMax, float yMin, float yMax ) {
    m_number    = number;
    m_planeCode = number / 2;
    m_z         = seg.z( 0.f );
    m_dxDy      = ( seg.x( 100.f ) - seg.x( 0.f ) ) / 100.f;
    m_dzDy      = ( seg.z( 100.f ) - seg.z( 0.f ) ) / 100.f;
    m_isX       = std::abs( m_dxDy ) < 0.010f;
    m_xMin      = xMin;
    m_xMax      = xMax;
    m_yMin      = yMin;
    m_yMax      = yMax;
  }
  unsigned int number() const { return m_number; }
  unsigned int planeCode() const { return m_planeCode; }
  float        z( float y = 0.f ) const { return m_z + m_dzDy * y; }
  float        dxDy() const { return m_dxDy; }
  float        dzDy() const { return m_dzDy; }
  bool         isX() const { return m_isX; }

  /// Is the point in the surrounding box?
  bool isInside( float x, float y ) const { return y < m_yMax && y > m_yMin && x > m_xMin && x < m_xMax; }

  bool isInsideY( float y ) const { return y < m_yMax && y > m_yMin; }

  float dxOnAFibre() const { return ( m_yMax - m_yMin ) * m_dxDy; }

  // for debugging
  friend std::ostream& operator<<( std::ostream& os, const PrHitZone& zone ) {
    os << "PrHitZone " << zone.m_number << ": z=" << zone.m_z << ", dxDy=" << zone.m_dxDy << ", dzDy=" << zone.m_dzDy
       << ", isX=" << zone.m_isX << std::endl;
    os << "xEdges=[" << zone.m_xMin << ", " << zone.m_xMax << "] "
       << "yEdges=[" << zone.m_yMin << ", " << zone.m_yMax << "]" << std::endl;
    return os;
  }

private:
  unsigned int m_number{};
  unsigned int m_planeCode{};
  float        m_z{};
  float        m_dxDy{};
  float        m_dzDy{};
  bool         m_isX{};
  float        m_xMin{};
  float        m_xMax{};
  float        m_yMin{};
  float        m_yMax{};
};
