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
#pragma once
//-----------------------------------------------------------------------------
/** @class TrackSelectorBase TrackSelectorBase.h
 *
 *  Common base class for track selectors
 *
 *  @author C. Jones  Christopher.Rob.Jones@cern.ch
 *
 *  @date   30/06/2009
 */
//-----------------------------------------------------------------------------

#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "TrackInterfaces/ITrackSelector.h"
#include <bitset>

using TrackType = LHCb::Event::Enum::Track::Type;

namespace details {
  constexpr int max_value_of_type() {
    int  i     = -1;
    auto types = LHCb::Event::Enum::Track::members_of<TrackType>();
    for ( auto type : types ) i = ( static_cast<int>( type ) > i ? static_cast<int>( type ) : i );
    return i;
  }

  class TrackTypeChecker {
    std::bitset<max_value_of_type() + 1> m_types{}; ///< Mapping linking track types
  public:
    template <typename Container>
    TrackTypeChecker( Container const& ts ) {
      for ( TrackType t : ts ) m_types.set( static_cast<int>( t ) );
    }
    TrackTypeChecker( std::initializer_list<TrackType> const& ts ) {
      for ( TrackType t : ts ) m_types.set( static_cast<int>( t ) );
    }

    bool operator()( TrackType t ) const { return m_types.test( static_cast<int>( t ) ); }

    friend StatusCode parse( TrackTypeChecker& ch, std::string const& s ) {
      std::vector<TrackType> types;
      using Gaudi::Parsers::parse;
      return parse( types, s ).andThen( [&] { ch = TrackTypeChecker{types}; } );
    }
    friend std::string toString( TrackTypeChecker const& tr ) {
      std::vector<std::string> strs;
      for ( auto type : LHCb::Event::Enum::Track::members_of<TrackType>() ) {
        if ( tr( type ) ) strs.push_back( toString( type ) );
      }
      using Gaudi::Utils::toString;
      return toString( strs );
    }
    friend std::ostream& toStream( TrackTypeChecker const& tr, std::ostream& os ) { return os << toString( tr ); }
  };
} // namespace details

class TrackSelectorBase : public extends<GaudiTool, ITrackSelector> {
public:
  /// constructer
  using extends::extends;

protected:
  // Check track type
  bool checkTrackType( const LHCb::Track& aTrack ) const {
    if ( m_selTypes.value()( aTrack.type() ) ) { return m_onlyBackwardTracks ? aTrack.isVeloBackward() : true; }
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << " -> Track type " << aTrack.type() << " is rejected" << endmsg;
    return false;
  }

private:
  /// Track types to accept
  Gaudi::Property<details::TrackTypeChecker> m_selTypes{this, "TrackTypes",
                                                        LHCb::Event::Enum::Track::members_of<TrackType>()};
  // select backward tracks
  Gaudi::Property<bool> m_onlyBackwardTracks{this, "OnlyBackwardTracks", false};
};
