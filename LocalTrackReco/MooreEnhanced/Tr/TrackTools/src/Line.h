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

namespace Tf {
  namespace Tsa {

    /** @class Line Line.h Line.h
     *
     *  y = mx + c !
     *
     *  @author M.Needham
     *  @date   31/05/2004
     */

    class Line final {
      double m_m;
      double m_c;

    public:
      Line( const double y1, const double y2, const double x1, const double x2 )
          : m_m{( y2 - y1 ) / ( x2 - x1 )}, m_c{y1 - ( m_m * x1 )} {}
      Line( const double m, const double c ) : m_m( m ), m_c( c ) {}
      Line( const double m, const double y, const double x ) : m_m( m ), m_c( y - ( m * x ) ) {}
      double value( const double x ) const { return ( m_c + ( m_m * x ) ); }
      double m() const { return m_m; }
      double c() const { return m_c; }
    };

  } // namespace Tsa
} // namespace Tf
