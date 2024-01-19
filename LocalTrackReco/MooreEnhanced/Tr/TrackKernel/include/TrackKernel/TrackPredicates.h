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
#include "Event/Track.h"

namespace TrackPredicates {
  struct Type {
    LHCb::Track::Types m_type;
    Type( LHCb::Track::Types atype ) : m_type( atype ) {}
    bool operator()( const LHCb::Track* track ) const { return track->type() == m_type; }
  };

  struct Flag {
    LHCb::Track::Flags m_flag;
    bool               m_positive;
    Flag( LHCb::Track::Flags flag, bool positive = true ) : m_flag( flag ), m_positive( positive ) {}
    bool operator()( const LHCb::Track* track ) const {
      return m_positive ? track->checkFlag( m_flag ) : !track->checkFlag( m_flag );
    }
  };

  struct VeloSide {
    int m_sign;
    VeloSide( int asign ) : m_sign( asign ) {}
    bool operator()( const LHCb::Track* track ) const {
      return track->firstState().tx() * m_sign * ( track->isVeloBackward() ? -1 : 1 ) > 0;
    }
  };

  struct MaxChisqPerDoF {
    double m_maxchisq;
    MaxChisqPerDoF( double maxChisqPerDof ) : m_maxchisq( maxChisqPerDof ) {}
    bool operator()( const LHCb::Track* track ) const { return track->chi2PerDoF() < m_maxchisq; }
  };
} // namespace TrackPredicates
