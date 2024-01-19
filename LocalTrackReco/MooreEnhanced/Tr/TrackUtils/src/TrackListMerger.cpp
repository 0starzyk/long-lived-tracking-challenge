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
/** @class TrackListMerger TrackListMerger.h
 *
 *  Merge different track lists.
 *
 *  @author Wouter Hulsbergen
 *  @date   05/01/2010
 */

#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "LHCbAlgs/MergingTransformer.h"
#include <string>

namespace LHCb {

  template <typename InputListType, typename OutputListType = InputListType>
  struct TrackListMergerT final
      : Algorithm::MergingTransformer<OutputListType( const Gaudi::Functional::vector_of_const_<InputListType>& )> {
    TrackListMergerT( const std::string& name, ISvcLocator* pSvcLocator )
        : Algorithm::MergingTransformer<OutputListType( const Gaudi::Functional::vector_of_const_<InputListType>& )>{
              name, pSvcLocator, {"InputLocations", {}}, {"OutputLocation", {}}} {}

    OutputListType operator()( const Gaudi::Functional::vector_of_const_<InputListType>& lists ) const override {

      OutputListType out;

      for ( const auto& list : lists ) {

        if ( list.empty() ) continue;

        for ( const auto& track : list ) {
          // make sure the track is not yet there!
          if ( std::find( out.begin(), out.end(), track ) == out.end() ) { out.insert( track ); }
        }
      }
      return out;
    }
  };

  /// TODO (wh,10/2022): we can replace TrackListMerger by TrackSelectionMerge, but that's for next time.
  using TrackSelectionMerger = TrackListMergerT<Track::Range, Track::Selection>;
  DECLARE_COMPONENT_WITH_ID( TrackSelectionMerger, "TrackSelectionMerger" )
  DECLARE_COMPONENT_WITH_ID( TrackListMergerT<Tracks>, "TrackListMerger" )
  using RecVertexListMerger = TrackListMergerT<RecVertex::Range, RecVertex::Selection>;
  DECLARE_COMPONENT_WITH_ID( RecVertexListMerger, "RecVertexListMerger" )
} // namespace LHCb
