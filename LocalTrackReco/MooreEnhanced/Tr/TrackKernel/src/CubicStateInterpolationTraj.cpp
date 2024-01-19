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

#include "TrackKernel/CubicStateInterpolationTraj.h"

#include "LHCbMath/Similarity.h"

namespace LHCb {

  LHCb::State CubicStateInterpolationTraj::state( double z ) const {
    return LHCb::State( {x( z ), y( z ), tx( z ), ty( z ), qop( z )}, covariance( z ), z,
                        LHCb::State::Location::LocationUnknown );
  }

  namespace {
    void transport( Gaudi::TrackSymMatrix& cov, const double dz ) {
      const double dz2 = dz * dz;
      cov( 0, 0 ) += dz2 * cov( 2, 2 ) + 2 * dz * cov( 2, 0 );
      cov( 1, 0 ) += dz2 * cov( 3, 2 ) + dz * ( cov( 3, 0 ) + cov( 2, 1 ) );
      cov( 2, 0 ) += dz * cov( 2, 2 );
      cov( 3, 0 ) += dz * cov( 3, 2 );
      cov( 4, 0 ) += dz * cov( 4, 2 );
      cov( 1, 1 ) += dz2 * cov( 3, 3 ) + 2 * dz * cov( 3, 1 );
      cov( 2, 1 ) += dz * cov( 3, 2 );
      cov( 3, 1 ) += dz * cov( 3, 3 );
      cov( 4, 1 ) += dz * cov( 4, 3 );
    }
  } // namespace

  Gaudi::TrackSymMatrix CubicStateInterpolationTraj::covariance( double z ) const {
    Gaudi::TrackSymMatrix cov;
    if ( z <= zbegin() ) {
      cov = m_covbegin;
      transport( cov, z - zbegin() );
    } else if ( z >= zend() ) {
      cov = m_covend;
      transport( cov, z - zend() );
    } else {
      // what is the right weight? FIXME!
      const double zfrac  = ( z - zbegin() ) / ( zend() - zbegin() );
      const double weight = ( 1 - zfrac ); // linear
      // double weight = 0.5 - 4*std::pow(zfrac-0.5,3) ; // cubic
      Gaudi::TrackSymMatrix covA = m_covbegin;
      transport( covA, z - zbegin() );
      Gaudi::TrackSymMatrix covB = m_covend;
      transport( covB, z - zend() );
      cov = weight * covA + ( 1 - weight ) * covB;
    }
    return cov;
  }

} // namespace LHCb
