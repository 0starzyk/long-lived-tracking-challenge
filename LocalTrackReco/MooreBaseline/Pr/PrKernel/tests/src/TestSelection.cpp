/*****************************************************************************\
* (c) Copyright 2018 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE utestSelection
#include "Event/Track_v2.h"
#include "PrKernel/PrSelection.h"
#include <array>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>

using Pr::Selection;

namespace {
  static constexpr std::array test_array          = {1, 3, 9, 14, 18};
  static constexpr std::array test_array_clone    = test_array;
  static constexpr std::array test_even_values    = {test_array[3], test_array[4]};
  static constexpr std::array test_even_addresses = {&test_array[3], &test_array[4]};

  struct Stubborn {
    Stubborn( Stubborn&& )      = delete;
    Stubborn( Stubborn const& ) = delete;
    Stubborn& operator=( Stubborn&& ) & = delete;
    Stubborn& operator=( Stubborn const& ) & = delete;
  };

  static constexpr std::array<Stubborn, 2> stubborn_array = {Stubborn{}, Stubborn{}};

  bool is_even( int num ) { return num % 2 == 0; }
  bool is_odd( int num ) { return num % 2 == 1; }
} // namespace

BOOST_AUTO_TEST_CASE( test_selection_set_operations ) {
  Selection all_sel{test_array}, even_sel{test_array, is_even}, odd_sel{test_array, is_odd};
  auto      union_sel = set_union( even_sel, odd_sel );
  BOOST_CHECK( union_sel == all_sel );               // same indices into the same container
  BOOST_CHECK( equal_values( union_sel, all_sel ) ); // => this also holds

  Selection clone_sel{test_array_clone};
  // selections point to different underlying storage => should throw
  BOOST_CHECK_THROW( includes( all_sel, clone_sel ), GaudiException );
  BOOST_CHECK_THROW( set_union( all_sel, clone_sel ), GaudiException );
  BOOST_CHECK_THROW( set_difference( all_sel, clone_sel ), GaudiException );
  BOOST_CHECK_THROW( set_intersection( all_sel, clone_sel ), GaudiException );
  BOOST_CHECK_THROW( set_symmetric_difference( all_sel, clone_sel ), GaudiException );

  // evens are a subset of the full list
  BOOST_CHECK( includes( all_sel, even_sel ) );

  auto all_that_are_not_even = set_difference( all_sel, even_sel );
  BOOST_CHECK( all_that_are_not_even == odd_sel );

  // overlap of even and odd
  auto both_even_and_odd = set_intersection( odd_sel, even_sel );
  BOOST_CHECK( both_even_and_odd.empty() );

  Selection        at_least_8{test_array, []( auto x ) { return x >= 8; }};
  auto             odd_and_at_least_8     = set_intersection( odd_sel, at_least_8 );
  std::vector<int> odd_and_at_least_8_ref = {9};
  BOOST_CHECK( equal_values( odd_and_at_least_8, odd_and_at_least_8_ref ) );

  // even or odd but not both
  auto even_or_odd_not_both = set_symmetric_difference( odd_sel, even_sel );
  BOOST_CHECK( even_or_odd_not_both == all_sel );

  auto             odd_or_at_least_8_but_not_both     = set_symmetric_difference( odd_sel, at_least_8 );
  std::vector<int> odd_or_at_least_8_but_not_both_ref = {1, 3, 14, 18};
  BOOST_CHECK( equal_values( odd_or_at_least_8_but_not_both, odd_or_at_least_8_but_not_both_ref ) );
}

BOOST_AUTO_TEST_CASE( test_selection_from_vector ) {
  std::vector<int> input = {0, 1, 2};
  Selection        selection{input}; // select everything
  BOOST_CHECK( equal_values( selection, input ) );
  BOOST_CHECK( equal_values( input, selection ) ); // check the arg ordering doesn't matter
}

BOOST_AUTO_TEST_CASE( test_selection_from_array ) {
  Selection selection{test_array};
  BOOST_CHECK( equal_values( selection, test_array ) );
}

BOOST_AUTO_TEST_CASE( test_selection_from_span ) {
  LHCb::span<int const> input = test_array;
  Selection             selection{input};
  BOOST_CHECK( equal_values( selection, input ) );
  BOOST_CHECK( equal_values( selection, test_array ) );
}

BOOST_AUTO_TEST_CASE( test_selection_of_immovable_objects ) {
  Selection<Stubborn> sel{stubborn_array};
  // for ( auto x : sel ) { } is not allowed
}

template <typename T>
void dump( T const& cont ) {
  for ( auto const& x : cont ) { std::cout << x << " "; }
  std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE( test_even_integer_selection ) {
  std::vector<int> integers = {0, 1, 2, 3, 4, 5}, even_integers_ref = {0, 2, 4};

  // Select the even integers
  Selection even_integers{integers, is_even};
  BOOST_CHECK( equal_values( even_integers, even_integers_ref ) );

  // Test iteration too
  std::vector<int> range_based_even, index_based_even;
  for ( auto const& value : even_integers ) range_based_even.push_back( value );
  for ( auto i = 0u; i < even_integers.size(); ++i ) index_based_even.push_back( even_integers[i] );

  BOOST_CHECK( range_based_even == index_based_even );

  // Add another selection on top of the selection
  Selection mod4_integers{even_integers, []( auto i ) { return i % 4 == 0; }};
  BOOST_CHECK( mod4_integers.size() == 2 );
  std::array const mod4_integers_ref = {0, 4};
  BOOST_CHECK( equal_values( mod4_integers, mod4_integers_ref ) );

  // Add an incompatible selection
  Selection even_and_odd{mod4_integers, is_odd};
  BOOST_CHECK( even_and_odd.empty() );
  BOOST_CHECK( even_and_odd.size() == 0 );
}

BOOST_AUTO_TEST_CASE( test_addresses_selection ) {
  Selection even_test_ints{test_array, is_even};
  BOOST_CHECK( test_even_addresses.size() == even_test_ints.size() );
  BOOST_CHECK( equal_values( even_test_ints, test_even_values ) );

  // Check that the selection refers to the same objects as we started with
  for ( auto i = 0u; i < even_test_ints.size(); ++i ) BOOST_CHECK( &even_test_ints[i] == test_even_addresses[i] );
}

BOOST_AUTO_TEST_CASE( test_selection_equality ) {
  // test array is {odd, odd, odd, even, even}
  Selection even_test_ints1{test_array, is_even}, even_test_ints2{test_array, is_even};

  // same predicate applied to the same input => should equal equal
  BOOST_CHECK( even_test_ints1 == even_test_ints2 );

  Selection even_test_ints_clone{test_array_clone, is_even};
  // same predicate applied to a different input sequence with the same values => should not compare equal
  BOOST_CHECK( even_test_ints1 != even_test_ints_clone );

  LHCb::span<int const> truncated_test_array{test_array.begin(), std::prev( test_array.end() )};
  Selection             odd_test_ints{test_array, is_odd}, truncated_odd_test_ints{truncated_test_array, is_odd};
  // same predicate applied to different numbers of elements starting in the same place => ??? :-)
  BOOST_CHECK( odd_test_ints == truncated_odd_test_ints );

  // the values should match though, given what we chose test_array to be
  BOOST_CHECK( equal_values( odd_test_ints, truncated_odd_test_ints ) );
}

using Track          = LHCb::Event::v2::Track;
using TrackSelection = Selection<Track>;
BOOST_AUTO_TEST_CASE( test_selection_constructors ) {
  // Shouldn't be able to default construct, a Selection should always be bound to a container
  static_assert( !std::is_default_constructible_v<TrackSelection> );

  static_assert( std::is_move_constructible_v<TrackSelection> );
  static_assert( std::is_nothrow_move_constructible_v<TrackSelection> );

  static_assert( std::is_copy_constructible_v<TrackSelection> );
  // static_assert(std::is_nothrow_copy_constructible_v<Track>); // not possible, indices are heap-allocated
}

BOOST_AUTO_TEST_CASE( test_selection_assignments ) {
  static_assert( std::is_move_assignable_v<TrackSelection> );
  static_assert( std::is_nothrow_move_assignable_v<TrackSelection> );

  static_assert( std::is_copy_assignable_v<TrackSelection> );
  // static_assert(std::is_nothrow_copy_assignable_v<Track>); // not possible, indices are heap-allocated

  static_assert( std::is_assignable_v<TrackSelection&, TrackSelection> );
  static_assert( std::is_nothrow_assignable_v<TrackSelection&, TrackSelection> );

  static_assert( !std::is_assignable_v<LHCb::span<Track const>&, TrackSelection> );
}
