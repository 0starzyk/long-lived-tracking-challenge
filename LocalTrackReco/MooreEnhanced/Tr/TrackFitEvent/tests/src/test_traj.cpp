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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_traj

#include "Event/StateZTraj.h"
#include <boost/test/unit_test.hpp>

#include "LHCbMath/SIMDWrapper.h"

#include <array>
#include <cstring>

using simd   = SIMDWrapper::best::types;
using scalar = SIMDWrapper::scalar::types;

namespace {

  struct SoAData {
    constexpr static unsigned int nb_entries = 5;
    constexpr static unsigned int nb_vec     = ( nb_entries + simd::size - 1 ) / simd::size;
    constexpr static unsigned int line_size  = nb_vec * simd::size;
    constexpr static unsigned int max_size   = line_size * 6;
    SOA_ACCESSOR( x, &data[0 * line_size] )
    SOA_ACCESSOR( y, &data[1 * line_size] )
    SOA_ACCESSOR( tx, &data[2 * line_size] )
    SOA_ACCESSOR( ty, &data[3 * line_size] )
    SOA_ACCESSOR( qop, &data[4 * line_size] )
    SOA_ACCESSOR( z, &data[5 * line_size] )
    alignas( 64 ) float data[max_size]{0};
    SoAData( std::array<const float, nb_entries> xs, std::array<const float, nb_entries> ys,
             std::array<const float, nb_entries> txs, std::array<const float, nb_entries> tys,
             std::array<const float, nb_entries> qops, std::array<const float, nb_entries> zs ) {
      memcpy( &data[0 * line_size], xs.data(), nb_entries * sizeof( float ) );
      memcpy( &data[1 * line_size], ys.data(), nb_entries * sizeof( float ) );
      memcpy( &data[2 * line_size], txs.data(), nb_entries * sizeof( float ) );
      memcpy( &data[3 * line_size], tys.data(), nb_entries * sizeof( float ) );
      memcpy( &data[4 * line_size], qops.data(), nb_entries * sizeof( float ) );
      memcpy( &data[5 * line_size], zs.data(), nb_entries * sizeof( float ) );
    }
  };

  SoAData data{
      {-36.1555, -1.01366, 787.603, -5.59514, -5.59514},                 // x
      {-273.555, -7.13817, -1087.02, -39.4269, -39.4269},                // y
      {0.000175838, 0.000175838, 0.000179528, 0.000175838, 0.000175838}, // tx
      {-0.0128336, -0.0166972, 0.213602, -0.0165059, -0.0165059},        // ty
      {-0.116527, -0.116762, -0.111832, -0.117043, -0.117043},           // qop
      {2320.28, 36.919, 9403, 313.081, 313.081}                          // z
  };

  std::array<LHCb::StateZTraj<float>, SoAData::nb_entries> aost;
  // initialization of AoS data
  struct InitAosZTraj {
    InitAosZTraj() {
      for ( unsigned int i = 0; i < SoAData::nb_entries; i++ ) {
        aost[i] = {data.x<scalar::float_v>( i ).cast(),   data.y<scalar::float_v>( i ).cast(),
                   data.tx<scalar::float_v>( i ).cast(),  data.ty<scalar::float_v>( i ).cast(),
                   data.qop<scalar::float_v>( i ).cast(), data.z<scalar::float_v>( i ).cast()};
      }
    }
  };

  std::array<LHCb::StateZTraj<simd::float_v, simd::size>, SoAData::nb_vec> soat;
  // initialization of SoA data
  struct InitSoaZTraj {
    InitSoaZTraj() {
      for ( unsigned int i = 0, j = 0; i < SoAData::nb_entries; i += simd::size, j++ ) {
        soat[j] = {data.x<simd::float_v>( i ),  data.y<simd::float_v>( i ),   data.tx<simd::float_v>( i ),
                   data.ty<simd::float_v>( i ), data.qop<simd::float_v>( i ), data.z<simd::float_v>( i )};
      }
    }
  };

} // namespace

BOOST_TEST_GLOBAL_FIXTURE( InitAosZTraj );
BOOST_TEST_GLOBAL_FIXTURE( InitSoaZTraj );

BOOST_AUTO_TEST_CASE( CHECK_OMEGA_Y ) {
  float pz = 2103.28564;
  for ( unsigned int i = 0; i < SoAData::nb_vec; ++i ) {
    auto                res = soat[i].omegay( pz );
    alignas( 64 ) float soa_oy[simd::size];
    res.store( soa_oy );
    for ( unsigned int j = 0; j < simd::size; j++ ) {
      auto k = i * simd::size + j;
      if ( k == SoAData::nb_entries ) break;
      auto aos_oy = aost[k].omegay( pz );
      BOOST_TEST_MESSAGE( " traj " << k );
      BOOST_CHECK_MESSAGE( aos_oy == soa_oy[j], aos_oy << " == " << soa_oy[j] );
    }
  }
}

BOOST_AUTO_TEST_CASE( CHECK_OMEGA_X ) {
  float pz = 5413.28564;
  for ( unsigned int i = 0; i < SoAData::nb_vec; ++i ) {
    auto                res = soat[i].omegax( pz );
    alignas( 64 ) float soa_ox[simd::size];
    res.store( soa_ox );
    for ( unsigned int j = 0; j < simd::size; j++ ) {
      auto k = i * simd::size + j;
      if ( k == SoAData::nb_entries ) break;
      auto aos_ox = aost[k].omegax( pz );
      BOOST_TEST_MESSAGE( " traj " << k );
      BOOST_CHECK_MESSAGE( aos_ox == soa_ox[j], aos_ox << " == " << soa_ox[j] );
    }
  }
}

BOOST_AUTO_TEST_CASE( CHECK_ATT ) {
  for ( unsigned int i = 0; i < SoAData::nb_vec; ++i ) {
    alignas( 64 ) float soadata[8 * simd::size]; // holds z, qOverP, cx[0], cx[1], cx[2], cy[0], cy[1], cy[2] in SoA
    auto                soa_z = soadata;
    soat[i].m_z.store( soa_z );
    auto soa_qOverP = soadata + simd::size;
    soat[i].m_qOverP.store( soa_qOverP );
    auto soa_cx = soadata + 2 * simd::size;
    auto soa_cy = soadata + 5 * simd::size;
    for ( unsigned int l = 0; l < 3; ++l ) {
      soat[i].m_cx[l].store( soa_cx + l * simd::size );
      soat[i].m_cy[l].store( soa_cy + l * simd::size );
    }
    for ( unsigned int j = 0; j < simd::size; ++j ) {
      auto k = i * simd::size + j;
      if ( k == SoAData::nb_entries ) break;
      BOOST_TEST_MESSAGE( " traj " << i );
      BOOST_CHECK_MESSAGE( aost[k].m_z == soa_z[j], aost[k].m_z << "==" << soa_z[j] );
      BOOST_CHECK_MESSAGE( aost[k].m_qOverP == soa_qOverP[j], aost[k].m_qOverP << "==" << soa_qOverP[j] );
      for ( unsigned int l = 0; l < 3; ++l ) {
        BOOST_CHECK_MESSAGE( aost[k].m_cx[l] == soa_cx[l * simd::size + j],
                             aost[k].m_cx[l] << "==" << soa_cx[l * simd::size + j] );
        BOOST_CHECK_MESSAGE( aost[k].m_cy[l] == soa_cy[l * simd::size + j],
                             aost[k].m_cy[l] << "==" << soa_cy[l * simd::size + j] );
      }
    }
  }
}

BOOST_AUTO_TEST_CASE( CHECK_POSITION ) {
  float pz = 4453.28;
  for ( unsigned int i = 0; i < SoAData::nb_vec; ++i ) {
    auto                soa_pos = soat[i].position( pz );
    alignas( 64 ) float soadata[3 * simd::size]; // holds X, Y, Z in SoA
    auto                soa_X = soadata;
    soa_pos.X().store( soa_X );
    auto soa_Y = soadata + simd::size;
    soa_pos.Y().store( soa_Y );
    auto soa_Z = soadata + 2 * simd::size;
    soa_pos.Z().store( soa_Z );
    for ( unsigned int j = 0; j < simd::size; ++j ) {
      auto k = i * simd::size + j;
      if ( k == SoAData::nb_entries ) break;
      auto aos_pos = aost[k].position( pz );
      BOOST_TEST_MESSAGE( " traj " << k );
      BOOST_CHECK_MESSAGE( aos_pos.X() == soa_X[j], "X: " << aos_pos.X() << " == " << soa_X[j] );
      BOOST_CHECK_MESSAGE( aos_pos.Y() == soa_Y[j], "Y: " << aos_pos.Y() << " == " << soa_Y[j] );
      BOOST_CHECK_MESSAGE( aos_pos.Z() == soa_Z[j], "Z: " << aos_pos.Z() << " == " << soa_Z[j] );
    }
  }
}

BOOST_AUTO_TEST_CASE( CHECK_MU_ESTIMATE ) {
  // Point pp { 2.9, 5.18, 9.2 };
  for ( unsigned int i = 0; i < SoAData::nb_vec; ++i ) {
    auto                res = soat[i].muEstimate( {2.9, 5.18, 9.2} );
    alignas( 64 ) float soa_mue[simd::size];
    res.store( soa_mue );
    for ( unsigned int j = 0; j < simd::size; ++j ) {
      auto k = i * simd::size + j;
      if ( k == SoAData::nb_entries ) break;
      auto aos_mue = aost[k].muEstimate( {2.9, 5.18, 9.2} );
      BOOST_TEST_MESSAGE( " traj " << k );
      BOOST_CHECK_MESSAGE( aos_mue == soa_mue[j], aos_mue << " == " << soa_mue[j] );
    }
  }
}

BOOST_AUTO_TEST_CASE( CHECK_DERIVATIVE ) {
  float pz = 3564.28;
  for ( unsigned int i = 0; i < SoAData::nb_vec; ++i ) {
    const auto          soa_deriv = soat[i].derivative( pz );
    alignas( 64 ) float soadata[6 * simd::size]; // holds deriv(a,b) for a=0,1 and b=2,3,4
    for ( unsigned int a = 0; a <= 1; ++a ) {
      for ( unsigned int b = 2; b <= 4; ++b ) { soa_deriv( a, b ).store( soadata + ( a * 3 + b - 2 ) * simd::size ); }
    }
    for ( unsigned int j = 0; j < simd::size; ++j ) {
      auto k = i * simd::size + j;
      if ( k == SoAData::nb_entries ) break;
      auto aos_deriv = aost[k].derivative( pz );
      BOOST_TEST_MESSAGE( " traj " << k );
      for ( unsigned int l = 2; l <= 4; ++l ) {
        BOOST_CHECK_MESSAGE( aos_deriv( 0, l ) == soadata[( l - 2 ) * simd::size + j],
                             "0:" << l << " " << aos_deriv( 0, l ) << " == " << soadata[( l - 2 ) * simd::size + j] );
        BOOST_CHECK_MESSAGE( aos_deriv( 1, l ) == soadata[( 1 + l ) * simd::size + j],
                             "1:" << l << " " << aos_deriv( 1, l ) << " == " << soadata[( 1 + l ) * simd::size + j] );
      }
    }
  }
}
