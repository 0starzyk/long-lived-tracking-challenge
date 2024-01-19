/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
/** @class TrackContainersMerger TrackContainersMerger.h
 *
 *  Merge different track containers into new one.
 *
 *  @author Andrii Usachov
 *  @date   29/05/2020
 */

#include "Event/Track.h"
#include "Event/Track_v3.h"
#include "Event/UniqueIDGenerator.h"
#include "GaudiKernel/GaudiException.h"
#include "LHCbAlgs/MergingTransformer.h"
#include "LHCbMath/SIMDWrapper.h"
#include "fmt/format.h"
#include <stdexcept>
#include <string>

template <typename TrackListType, typename TrackLists = Gaudi::Functional::vector_of_const_<TrackListType>>
struct TrackContainersMergerT final : LHCb::Algorithm::MergingTransformer<TrackListType( const TrackLists& )> {
  TrackContainersMergerT( const std::string& name, ISvcLocator* pSvcLocator )
      : LHCb::Algorithm::MergingTransformer<TrackListType( const TrackLists& )>(
            name, pSvcLocator, {"InputLocations", {}}, {"OutputLocation", {}} ) {}

  TrackListType operator()( const TrackLists& lists ) const override {

    TrackListType out;
    for ( const auto& list : lists ) {
      for ( const auto& track : list ) {
        // make sure the track is not yet there!
        if ( std::find( out.begin(), out.end(), track ) == out.end() ) { out.insert( new LHCb::Track( *track ) ); }
      }
    }
    return out;
  }
};

namespace LHCb::Event::v3 {

  using UniqueIDGeneratorLists = Gaudi::Functional::vector_of_const_<LHCb::UniqueIDGenerator>;
  template <typename TrackListType, typename TrackLists = Gaudi::Functional::vector_of_const_<TrackListType>>
  struct TrackContainersMerger final
      : LHCb::Algorithm::MergingTransformer<TrackListType( const TrackLists&, const UniqueIDGeneratorLists& )> {

    using base_class =
        LHCb::Algorithm::MergingTransformer<TrackListType( const TrackLists&, const UniqueIDGeneratorLists& )>;
    using KeyValues = typename base_class::KeyValues;
    using KeyValue  = typename base_class::KeyValue;

    TrackContainersMerger( const std::string& name, ISvcLocator* pSvcLocator )
        : LHCb::Algorithm::MergingTransformer<TrackListType( const TrackLists&, const UniqueIDGeneratorLists& )>(
              name, pSvcLocator,
              {KeyValues{"InputLocations", {}},
               KeyValues{"InputUniqueIDGenerators", {LHCb::UniqueIDGeneratorLocation::Default}}},
              KeyValue{"OutputLocation", ""} ) {}

    TrackListType operator()( const TrackLists& lists, const UniqueIDGeneratorLists& unique_id_list ) const override {

      // -- The MergingTransformer only allows lists as input, at least make sure there is not more than
      // -- one element in the unique_id_list
      if ( unique_id_list.size() != 1 )
        throw GaudiException( "Only one Unique ID generator can be used", "TrackContainersMerger",
                              StatusCode::FAILURE );

      if ( lists.size() == 0 ) { return TrackListType{LHCb::Event::v3::TrackType::Unknown, *unique_id_list.begin()}; }

      // -- Get the track type of the first container in the list.
      const auto type = ( *lists.begin() ).type();

      TrackListType out{type, *unique_id_list.begin(), lists[0].zipIdentifier()};

      for ( const auto& list : lists ) {
        auto listType = list.type();

        // -- Only one track type is allowed in the final container
        if ( listType != type ) {
          std::string err = fmt::format( "{} does not correspond to required track type {}",
                                         static_cast<int>( listType ), static_cast<int>( type ) );

          throw GaudiException( err, "TrackContainersMerger", StatusCode::FAILURE );
        }

        using dType = SIMDWrapper::best::types;
        for ( auto const& track : list.simd() ) {
          auto       loop_mask = track.loop_mask();
          auto const t         = track.offset();
          out.template copy_back<dType>( list, t, loop_mask );
        }
      }
      return out;
    }
  };
} // namespace LHCb::Event::v3

DECLARE_COMPONENT_WITH_ID( TrackContainersMergerT<LHCb::Tracks>, "TrackContainersMerger" )
DECLARE_COMPONENT_WITH_ID( LHCb::Event::v3::TrackContainersMerger<LHCb::Event::v3::Tracks>, "TrackContainersMergerV3" )
