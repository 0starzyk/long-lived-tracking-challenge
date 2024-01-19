/*****************************************************************************\
* (c) Copyright 2019-20 CERN for the benefit of the LHCb Collaboration        *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
// Make sure these are turned on for these tests, even if they are disabled
// elsewhere in the build
#define ZIPPING_SEMANTIC_CHECKS

#include "Event/GeneratePrFittedForwardTracks.h"
#include "Event/PrFittedForwardTracks.h"
#include "Event/RecVertex_v2.h"
#include "Event/StateParameters.h"
#include "Event/UniqueIDGenerator.h"
#include "GaudiKernel/SerializeSTL.h"
#include "SelKernel/VertexRelation.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE utestZipInfrastructure
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <random>

using Tracks   = LHCb::Pr::Fitted::Forward::Tracks;
using Vertices = LHCb::Event::v2::RecVertices;

static const LHCb::UniqueIDGenerator unique_id_gen;

auto make_covmatrix( double diag = 1e-2 ) {
  Gaudi::SymMatrix3x3 cov;
  cov( 0, 0 ) = cov( 1, 1 ) = cov( 2, 2 ) = diag;
  return cov;
}

Vertices make_vertices() {
  Vertices vertices;
  auto constexpr nvertices = 3;
  std::mt19937                    gen( 24 ); // Random engine with fixed seed
  std::normal_distribution<float> xy_dist{0.f, 0.05f}, z_dist{0.f, 1.f};
  for ( auto i = 0; i < nvertices; ++i ) {
    auto x = xy_dist( gen ), y = xy_dist( gen ), z = z_dist( gen );
    vertices.emplace_back( Gaudi::XYZPoint{x, y, z}, make_covmatrix(), LHCb::Event::v2::Track::Chi2PerDoF{1.f, 42} );
  }
  return vertices;
}

BOOST_AUTO_TEST_CASE( test_dummy_tracks ) {
  auto ntracks         = 20;
  auto tracks          = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto iterable_tracks = LHCb::Event::make_zip( tracks );
  for ( auto const track_chunk : iterable_tracks ) { popcount( track_chunk.loop_mask() ); }
}

auto make_relations( Tracks const& tracks, Vertices const& vertices = make_vertices() ) {
  return Sel::calculateBestVertices( tracks, vertices );
}

BOOST_AUTO_TEST_CASE( test_making_vertex_relations ) {
  auto ntracks = 20;
  make_relations( LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen ) );
}

BOOST_AUTO_TEST_CASE( test_zipping_relations ) {
  auto ntracks   = 20;
  auto tracks    = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto relations = make_relations( tracks );
  auto zipped    = LHCb::Event::make_zip( tracks, relations );
  for ( auto const& zipped_chunk : zipped ) {
    static_cast<void>( zipped_chunk.bestPV() );
    static_cast<void>( zipped_chunk.closestToBeamState() );
  }
}

// Check that we can go from the original containers to a zip and back again
BOOST_AUTO_TEST_CASE( test_decomposing_zip ) {
  // Make the owning containers
  auto ntracks   = 20;
  auto tracks    = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto relations = make_relations( tracks );
  // Zip them together
  auto zip = LHCb::Event::make_zip( tracks, relations );
  // Get references to the owning containers
  auto const& tracks_ref    = zip.get<decltype( tracks )>();
  auto const& relations_ref = zip.get<decltype( relations )>();
  // Check that those references do indeed point to the same containers
  BOOST_CHECK( &tracks == &tracks_ref );
  BOOST_CHECK( &relations == &relations_ref );
}
BOOST_AUTO_TEST_CASE( test_using_zipped_relations ) {
  auto ntracks   = 20;
  auto tracks    = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto vertices  = make_vertices();
  auto relations = make_relations( tracks, vertices );
  auto iterable  = LHCb::Event::make_zip( tracks );
  auto zipped    = LHCb::Event::make_zip( tracks, relations );
  // For the scalar case when masks are bools
  for ( auto const& chunk : zipped ) {
    // This should use the zipped relations
    auto bestPVs = Sel::getBestPVRel( chunk, vertices );
    // This should calculate the same thing without using the zipped result
    auto calcrel = Sel::calculateBestVertex( chunk, vertices );
    // Check they're either the same or the elements are out of range
    BOOST_CHECK( all( ( bestPVs == calcrel ) || !chunk.loop_mask() ) );
  }
  // For the first chunk, explicitly check that we get the same results whether
  // or not the containers are zipped
  BOOST_CHECK( all( Sel::getBestPVRel( zipped[0], vertices ) == Sel::getBestPVRel( iterable[0], vertices ) ) );
}

BOOST_AUTO_TEST_CASE( filter_tracks_using_zip ) {
  auto ntracks   = 20;
  auto tracks    = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto relations = make_relations( tracks );
  auto zipped    = LHCb::Event::make_zip( tracks, relations );
  using dType    = decltype( zipped )::default_simd_t;
  // Produce some filtered tracks using information from both parts of the zip
  Tracks output{nullptr, unique_id_gen, LHCb::Event::Enum::Track::History::Unknown};
  for ( auto const& chunk : zipped ) {
    auto loop_mask = chunk.loop_mask();
    auto filt_mask = ( chunk.pt() > 400.f ) && ( chunk.bestPV().ipchi2() > 1.f );
    output.copy_back<dType>( tracks, chunk.offset(), loop_mask && filt_mask );
  }
}

BOOST_AUTO_TEST_CASE( new_struct_from_zip ) {
  constexpr bool print{false};
  auto           ntracks   = 20;
  auto           tracks    = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto           relations = make_relations( tracks );

  // Make an iterable (non-owning) zip
  auto zipped = LHCb::Event::make_zip( tracks, relations );

  // Check we can iterate over it, and manually check the retention of the cut
  // that we're about to apply
  std::size_t passing_cut{0};
  for ( auto const& chunk : zipped ) {
    auto pt     = chunk.pt();
    auto ipchi2 = chunk.bestPV().ipchi2();
    auto index  = chunk.bestPV().index();
    passing_cut += popcount( chunk.loop_mask() && ( pt > 400.f ) && ( ipchi2 > 1.f ) );
    if ( print ) { std::cout << "pt " << pt << " ipchi2 " << ipchi2 << " index " << index << std::endl; }
  }

  // Try and make a new structure containing the fields from both 'tracks' and
  // 'relations', applying some selection
  auto new_data =
      zipped.filter( []( auto const& chunk ) { return ( chunk.pt() > 400.f ) && ( chunk.bestPV().ipchi2() > 1.f ); } );

  // Make a new non-owning iterable view into this new structure
  auto new_iterable = LHCb::Event::make_zip( new_data );

  // Check we retained the right number
  BOOST_CHECK( passing_cut == new_iterable.size() );

  // The handling of the merged data type returned by `filter` should be such
  // that the iterable version of it is the same type as the original zip.
  static_assert( std::is_same_v<decltype( new_iterable ), decltype( zipped )> );

  // Check we can form a loop over this one too [this is maybe redundant...]
  for ( auto const& chunk : new_iterable ) {
    auto pt     = chunk.pt();
    auto ipchi2 = chunk.bestPV().ipchi2();
    auto index  = chunk.bestPV().index();
    auto mask   = chunk.loop_mask();
    if ( print ) {
      std::cout << "pt " << pt << " ipchi2 " << ipchi2 << " index " << index << " mask " << mask << std::endl;
    }
  }

  // Try printing with a scalar loop too
  for ( auto const& chunk : new_iterable.with<SIMDWrapper::InstructionSet::Scalar>() ) {
    auto pt     = chunk.pt();
    auto ipchi2 = chunk.bestPV().ipchi2();
    auto index  = chunk.bestPV().index();
    if ( print ) { std::cout << "pt " << pt << " ipchi2 " << ipchi2 << " index " << index << std::endl; }
  }
}

BOOST_AUTO_TEST_CASE( test_semantic_check ) {
  // Make two incompatible zip containers
  auto ntracks      = 20;
  auto tracks       = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto other_tracks = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto relations    = make_relations( other_tracks );

  // These two containers have different zip identifiers, so trying to zip them
  // together should throw.
  BOOST_CHECK_THROW( LHCb::Event::make_zip( tracks, relations ), GaudiException );
}

BOOST_AUTO_TEST_CASE( test_size_check ) {
  // Make two compatible containers
  auto ntracks   = 20;
  auto tracks    = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto relations = make_relations( tracks );

  // Change the size of one of them, so zipping them becomes invalid
  tracks = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks - 1, unique_id_gen );

  // Check that we actually get an error
  BOOST_CHECK_THROW( LHCb::Event::make_zip( tracks, relations ), GaudiException );
}

BOOST_AUTO_TEST_CASE( test_growing_zip ) {
  // Make two compatible containers
  auto ntracks   = 20;
  auto tracks    = LHCb::Pr::Fitted::Forward::generate_tracks( ntracks, unique_id_gen );
  auto relations = make_relations( tracks );

  // Make some iterable tracks
  auto iterable_tracks = LHCb::Event::make_zip( tracks );

  // Check that calling make_zip again is a no-op
  auto iterable_tracks_2 = LHCb::Event::make_zip( iterable_tracks );
  static_assert( std::is_same_v<decltype( iterable_tracks ), decltype( iterable_tracks_2 )> );
  BOOST_CHECK( iterable_tracks == iterable_tracks_2 );

  // Extend this zip with relations
  auto tracks_with_rels = LHCb::Event::make_zip( iterable_tracks, relations );

  // Make the same zip directly
  auto tracks_with_rels_2 = LHCb::Event::make_zip( tracks, relations );

  // And with the opposite argument order
  auto tracks_with_rels_3 = LHCb::Event::make_zip( relations, tracks );

  // Check that we get the same type via all three routes
  static_assert( std::is_same_v<decltype( tracks_with_rels ), decltype( tracks_with_rels_2 )> );
  static_assert( std::is_same_v<decltype( tracks_with_rels ), decltype( tracks_with_rels_3 )> );

  // Do some simple, scalar checking
  auto scalar_1 = LHCb::Event::make_zip<SIMDWrapper::Scalar>( tracks_with_rels );
  auto scalar_2 = LHCb::Event::make_zip<SIMDWrapper::Scalar>( tracks_with_rels_2 );

  // Check they ended up with the same underlying pointers
  BOOST_CHECK( scalar_1 == scalar_2 );

  // Do some basic checking that the two zips yield the right values
  for ( auto i = 0; i < (int)tracks.size(); ++i ) {
    BOOST_CHECK( all( scalar_1[i].pt() == scalar_2[i].pt() ) );
    BOOST_CHECK( all( scalar_1[i].bestPV().index() == scalar_2[i].bestPV().index() ) );
  }
}
