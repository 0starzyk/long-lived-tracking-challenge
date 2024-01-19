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
//============================================================================
/** @file Track1DTabFunc.cpp
 *
 *  Implementation file for class : Track::TabulatedFunction1D
 *
 *  @author Paul Seyfert   Paul.Seyfert@cern.ch
 *  @date   2016-12-20
 */
//============================================================================

// STL
#include <algorithm>

// GaudiKernel
#include "GaudiKernel/GaudiException.h"

// local
#include "TrackKernel/Track1DTabFunc.h"

using namespace Track;

//============================================================================

// Constructor
TabulatedFunction1D::TabulatedFunction1D( std::initializer_list<float> x ) : m_xedges( x ) {
  if ( m_xedges.size() < 2 ) {
    throw GaudiException( "TabulatedFunction1D() : must be initialized with more than one bin edge",
                          "*Track::TabulatedFunction1D*", StatusCode::FAILURE );
  }
  m_width = 1.f / ( m_xedges.size() - 1 );
  if ( !std::is_sorted( m_xedges.begin(), m_xedges.end() ) ) {
    throw GaudiException( "TabulatedFunction1D() : must be initialized with sorted braced initializer list",
                          "*Track::TabulatedFunction1D*", StatusCode::FAILURE );
  }
}

//============================================================================

// evaluation function
float TabulatedFunction1D::value( float x ) const {
  // out of range check
  if ( x <= m_xedges.front() ) return 0.f;
  if ( x >= m_xedges.back() ) return 1.f;

  // iterator to the first element that is not smaller than x
  // may be end() - if x is larger than the last element
  // may be begin() - if x is smaller than the first element
  auto up = std::lower_bound( m_xedges.begin(), m_xedges.end(), x );
  // iterator to the last element that is smaller than x
  // (may be out of range)
  auto low = up - 1;

  // y-value for the lower edge of the x-bin we're in
  float edge = ( low - m_xedges.begin() );

  // by what fraction did we enter the x bin
  float add = ( x - ( *low ) ) / ( ( *up ) - ( *low ) );

  return m_width * ( edge + add );
}

//============================================================================
