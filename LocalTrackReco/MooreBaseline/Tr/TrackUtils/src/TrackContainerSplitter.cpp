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

/** @class TrackContainerSplitter TrackContainerSplitter.h
 *
 *  Split container to two (passed and rejected) wrt to given selection criteria
 *
 *  @author Andrii Usachov
 *  @date   30/04/2020
 */

#include "Event/Track.h"
#include "Functors/Function.h"
#include "Functors/with_functors.h"
#include "LHCbAlgs/Transformer.h"

namespace {
  /// The output data
  using OutConts = std::tuple<LHCb::Tracks, LHCb::Tracks>;

  struct TrackPredicate {
    using Signature                    = bool( const LHCb::Track& );
    static constexpr auto PropertyName = "Code";
  };
} // namespace

class TrackContainerSplitter final
    : public with_functors<LHCb::Algorithm::MultiTransformer<OutConts( const LHCb::Tracks& )>, TrackPredicate> {
public:
  TrackContainerSplitter( const std::string& name, ISvcLocator* pSvcLocator )
      : with_functors( name, pSvcLocator, {KeyValue{"TracksInContainer", ""}},
                       {KeyValue{"PassedContainer", ""}, KeyValue{"RejectedContainer", ""}} ) {}

  OutConts operator()( const LHCb::Tracks& input_tracks ) const override {
    OutConts out;
    auto& [passed, rejected] = out;
    m_inputs += input_tracks.size();
    auto const& pred = getFunctor<TrackPredicate>();
    for ( auto& trk : input_tracks ) {
      auto& container = ( pred( *trk ) ? passed : rejected );
      container.insert( new LHCb::Track( *trk ) );
    }
    m_passed += passed.size();

    return out;
  }

private:
  mutable Gaudi::Accumulators::StatCounter<> m_inputs{this, "#inputs"};
  mutable Gaudi::Accumulators::StatCounter<> m_passed{this, "#passed"};
};

DECLARE_COMPONENT( TrackContainerSplitter )
