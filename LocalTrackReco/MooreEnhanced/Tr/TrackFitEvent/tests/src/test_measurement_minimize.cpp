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
#define BOOST_TEST_MODULE test_measurement_minimize

#include <array>
#include <boost/test/unit_test.hpp>

#include "Event/Measurement.h"
#include "Event/State.h"
#include "Event/StateZTraj.h"
#include "Kernel/LHCbID.h"
#include "Kernel/LineTraj.h"
#include "VPDet/DeVPSensor.h"

namespace {

  // input are begin and end point of measurement trajectory
  LHCb::LineTraj<double> const meas_trajectory( {-0x1.3a22b4p+11, -0x1.814f77ffffff5p+3, 0x1.23a876p+13},
                                                {-0x1.1fba36p+11, -0x1.2f59a2p+11, 0x1.2362e4p+13} );

  // need this for the constructor but it doesn't need to actually point to something
#ifdef USE_DD4HEP
  DeVPSensor det{};
#else
  DeVPSensor* det = nullptr;
#endif
  // input: ID, z, trajectory, err, SubInfo
  LHCb::Measurement meas{LHCb::LHCbID(), 0x1.23a898p+13, meas_trajectory, 0, det};

  //  input is: x,y,tx,ty,qop,z, Bfield-vec
  LHCb::StateZTraj<double> const reference_no_b_field( -0x1.35675a91ee8dbp+11, -0x1.c68b550c7a396p+8,
                                                       -0x1.113aba118e781p-1, -0x1.989f13e3421ddp-5, -0x1.5b1b54p-12,
                                                       0x1.23a898p+13, {0, 0, 0} );
  //  input is: x,y,tx,ty,qop,z, Bfield-vec
  LHCb::StateZTraj<double> const reference_with_b_field( -0x1.35675a91ee8dbp+11, -0x1.c68b550c7a396p+8,
                                                         -0x1.113aba118e781p-1, -0x1.989f13e3421ddp-5, -0x1.5b1b54p-12,
                                                         0x1.23a898p+13,
                                                         {0x1.65d6b4p-17, -0x1.a162dp-15, -0x1.828402p-16} );
} // namespace

BOOST_AUTO_TEST_CASE( CHECK_MINIMIZE_NO_BFIELD, *boost::unit_test::tolerance( double{1e-12} ) ) {
  BOOST_TEST_MESSAGE( " Testing minimize with reference trajectory without B-field " );

  double constexpr zState_true{0x1.239bbc7ea9dbep+13};
  double constexpr sMeas_true{-0x1.7ff114977e4aep+9};
  double constexpr doca_true{-0x1.79e1f9add5f0cp-8};
  std::array<double, 3> const unitPoca_true{-0x1.c1a1059dce72ap-1, -0x1.33bbb15201d7dp-4, -0x1.e3ba40fe8ecd2p-2};
  std::array<double, 5> constexpr H_true{-0x1.c1a1059dce72ap-1, -0x1.33bbb15201d7dp-4, 0x1.695142ca30c3dp+0,
                                         0x1.ee952b35424e4p-4, -0x0p+0};

  // std::visit([&](const auto& min) {
  std::visit( Gaudi::overload(
                  [&]( const LHCb::Minimize1DResult& min ) {
                    auto const [zState, sMeas, doca, H, unitPocaVector] = min;

                    // convert Gaudi::xyzvector vector because we can't index it...
                    std::array<double, 3> unitpoca{unitPocaVector.x(), unitPocaVector.y(), unitPocaVector.z()};

                    BOOST_TEST( zState == zState_true, zState << " == " << zState_true );
                    BOOST_TEST( sMeas == sMeas_true, sMeas << " == " << sMeas_true );
                    BOOST_TEST( doca == doca_true, doca << " == " << doca_true );
                    for ( size_t i{0}; i < H_true.size(); ++i ) {
                      BOOST_TEST_CONTEXT( "H index " << i )
                      BOOST_TEST( H( 0, i ) == H_true[i] );
                    }

                    for ( size_t i{0}; i < unitPoca_true.size(); ++i ) {
                      BOOST_TEST_CONTEXT( "UnitPocaVector index " << i )
                      BOOST_TEST( unitpoca[i] == unitPoca_true[i] );
                    }
                  },
                  [&]( const LHCb::Minimize2DResult& ) {
                    throw std::logic_error( "Test is 1D but 2D result returned by minimize." );
                  } ),
              minimize( meas, reference_no_b_field, 0x1.23a898p+13 ) );
}

BOOST_AUTO_TEST_CASE( CHECK_MINIMIZE_WITH_BFIELD, *boost::unit_test::tolerance( double{1e-12} ) ) {
  BOOST_TEST_MESSAGE( " Testing minimize with reference trajectory with B-field " );

  double constexpr zState_true{0x1.239bb9c9bd50ap+13};
  double constexpr sMeas_true{-0x1.7ff11492f8902p+9};
  double constexpr doca_true{-0x1.790870d158c79p-8};
  std::array<double, 3> const unitPoca_true{-0x1.c1a1e63a91dbbp-1, -0x1.33bc5aa823285p-4, -0x1.e3b6f721652d1p-2};
  std::array<double, 5> constexpr H_true{-0x1.c1a1e63a91dbbp-1, -0x1.33bc5aa823285p-4, 0x1.699ddaebbe0a1p+0,
                                         0x1.eefaf977e9a83p-4, -0x1.4167176a60aacp-6};

  std::visit( Gaudi::overload(
                  [&]( const LHCb::Minimize1DResult& min ) {
                    auto const [zState, sMeas, doca, H, unitPocaVector] = min;

                    // convert Gaudi::xyzvector vector because we can't index it...
                    std::array<double, 3> unitpoca{unitPocaVector.x(), unitPocaVector.y(), unitPocaVector.z()};

                    BOOST_TEST( zState == zState_true, zState << " == " << zState_true );
                    BOOST_TEST( sMeas == sMeas_true, sMeas << " == " << sMeas_true );
                    BOOST_TEST( doca == doca_true, doca << " == " << doca_true );

                    for ( size_t i{0}; i < H_true.size(); ++i ) {
                      BOOST_TEST_CONTEXT( "H index " << i )
                      BOOST_TEST( H( 0, i ) == H_true[i] );
                    }

                    for ( size_t i{0}; i < unitPoca_true.size(); ++i ) {
                      BOOST_TEST_CONTEXT( "UnitPocaVector index " << i )
                      BOOST_TEST( unitpoca[i] == unitPoca_true[i] );
                    }
                  },
                  [&]( const LHCb::Minimize2DResult& ) {
                    throw std::logic_error( "Test is 1D but 2D result returned by minimize." );
                  } ),
              minimize( meas, reference_with_b_field, 0x1.23a898p+13 ) );
}
