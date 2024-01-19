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

/** @class TrackListRefiner TrackListRefiner.h
 *
 *  Make a subselection of a track list
 *
 *  @author Wouter Hulsbergen
 *  @date   05/01/2010
 */

#include "Event/Track.h"
#include "Functors/with_functors.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Transformer.h"

namespace {
  struct TrackPredicate {
    using Signature                    = bool( const LHCb::Track& );
    static constexpr auto PropertyName = "Code";
  };
} // namespace

class TrackListRefiner
    : public with_functors<LHCb::Algorithm::Transformer<LHCb::Track::Selection( const LHCb::Track::Range& )>,
                           TrackPredicate> {
public:
  TrackListRefiner( const std::string& name, ISvcLocator* pSvcLocator )
      : with_functors( name, pSvcLocator, KeyValue{"inputLocation", {}}, KeyValue{"outputLocation", {}} ) {}

  LHCb::Track::Selection operator()( const LHCb::Track::Range& tracksin ) const override {
    m_seeds += tracksin.size();
    // TODO: can we use std::transform -- i.e. is there an 'inserter' for LHCb::Track::Selection?
    LHCb::Track::Selection tracksout;
    auto const&            pred = getFunctor<TrackPredicate>();
    for ( const auto& trk : tracksin ) {
      if ( pred( *trk ) ) tracksout.insert( trk );
    }
    m_passed += tracksout.size();
    return tracksout;
  }

private:
  mutable Gaudi::Accumulators::StatCounter<> m_seeds{this, "#seeds"};
  mutable Gaudi::Accumulators::StatCounter<> m_passed{this, "#passed"};
};

DECLARE_COMPONENT( TrackListRefiner )
