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
#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "GaudiAlg/FilterPredicate.h"
#include "PrKernel/PrSelection.h"
#include <vector>

template <typename Container>
class PrCheckEmptyFilter : public Gaudi::Functional::FilterPredicate<bool( const Container& )> {

public:
  // Standard Constructor
  PrCheckEmptyFilter( const std::string& name, ISvcLocator* pSvcLocator )
      : Gaudi::Functional::FilterPredicate<bool( const Container& )>( name, pSvcLocator, {"inputLocation", ""} ) {}

  bool operator()( const Container& inputs ) const override {
    if ( this->msgLevel( MSG::DEBUG ) )
      this->debug() << this->inputLocation() << " is " << ( inputs.empty() ? "empty" : "not empty" ) << endmsg;
    return !inputs.empty();
  }
};

namespace {
  using Selection__Track_v1 = Pr::Selection<LHCb::Event::v1::Track>;
} // namespace

DECLARE_COMPONENT_WITH_ID( PrCheckEmptyFilter<LHCb::Track::Selection>, "PrCheckEmptyTracks" )
DECLARE_COMPONENT_WITH_ID( PrCheckEmptyFilter<Selection__Track_v1>, "PrCheckEmpty__Track_v1" )
DECLARE_COMPONENT_WITH_ID( PrCheckEmptyFilter<LHCb::RecVertex::Selection>, "VertexListFilter" )
