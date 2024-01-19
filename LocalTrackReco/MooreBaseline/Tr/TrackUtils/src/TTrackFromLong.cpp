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
/** @class TTrackFromLong TTrackFromLong.h
 *
 *  Fake a T seed from a long track
 *
 *  @author M.Needham
 *  @date   30/05/2006
 */

#include "Event/Track.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "LHCbAlgs/ScalarTransformer.h"
#include <optional>
#include <string>
using namespace Gaudi::Units;

namespace {
  auto isT = []( const LHCb::LHCbID& id ) { return id.isFT(); };
}

struct TTrackFromLong : LHCb::Algorithm::ScalarTransformer<TTrackFromLong, LHCb::Tracks( const LHCb::Tracks& )> {
  TTrackFromLong( const std::string& name, ISvcLocator* pSvcLocator )
      : ScalarTransformer( name, pSvcLocator, {KeyValue( "inputLocation", LHCb::TrackLocation::Forward )},
                           KeyValue( "outputLocation", LHCb::TrackLocation::Seed ) ) {}

  using ScalarTransformer::  operator();
  std::optional<LHCb::Track> operator()( const LHCb::Track& trk ) const // Note: this is NOT a virtual function!!
  {
    auto count = std::count_if( begin( trk.lhcbIDs() ), end( trk.lhcbIDs() ), isT );
    if ( count < 5 ) return {};

    LHCb::Track seed;
    seed.setType( LHCb::Track::Types::Ttrack );
    seed.setHistory( LHCb::Track::History::PrSeeding );
    const auto& lastState = trk.closestState( 9000. * mm );
    LHCb::State tState;
    tState.setLocation( LHCb::State::Location::AtT );
    tState.setState( lastState.stateVector() );
    tState.setZ( lastState.z() );
    tState.setCovariance( lastState.covariance() );
    seed.addToStates( tState );
    seed.setPatRecStatus( LHCb::Track::PatRecStatus::PatRecIDs );
    for ( const auto& id : trk.lhcbIDs() ) {
      if ( isT( id ) ) seed.addToLhcbIDs( id );
    }
    return seed;
  };
};

DECLARE_COMPONENT( TTrackFromLong )
