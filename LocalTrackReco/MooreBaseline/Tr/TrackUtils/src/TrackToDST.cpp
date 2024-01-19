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
/** @class TrackToDST TrackToDST.h
 *
 *  @author M.Needham
 *  @date   30/05/2006
 */

#include "Event/State.h"
#include "Event/Track.h"
#include "GaudiKernel/SerializeSTL.h"
#include <map>
#include <string>
#include <vector>

namespace Gaudi::Parsers {
  StatusCode parse( std::vector<LHCb::State::Location>& result, const std::string& input );
}
namespace Gaudi::Utils {
  // the default implemenation of `toStream` doesn't quote the individual elements...
  std::ostream& toStream( const std::vector<LHCb::State::Location>& value, std::ostream& os ) {
    return GaudiUtils::details::ostream_joiner( os << "[ ", value, ", ",
                                                []( std::ostream& os, LHCb::State::Location loc ) -> std::ostream& {
                                                  return os << std::quoted( LHCb::State::LocationToString( loc ) );
                                                } )
           << " ]";
  }
} // namespace Gaudi::Utils

#include "GaudiAlg/GaudiAlgorithm.h"

class TrackToDST : public GaudiAlgorithm {

public:
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode execute() override;

private:
  // CRJ : Orignal list
  //   m_veloStrings = list_of("ClosestToBeam");
  //   m_longStrings = list_of("ClosestToBeam")("FirstMeasurement")("BegRich1")("BegRich2")("V0Vertex");
  //   m_tStrings = list_of("FirstMeasurement")( "BegRich2");
  //   m_downstreamStrings = list_of("BegRich1")("FirstMeasurement")("BegRich2")("V0Vertex");
  //   m_upstreamStrings = list_of("ClosestToBeam")("FirstMeasurement")("BegRich1");
  //   m_muonStrings     = list_of("ClosestToBeam")("BegRich1")("BegRich2")("Muon");
  using SLocations = std::vector<LHCb::State::Location>;
  // (Slightly) reduced list
  Gaudi::Property<SLocations> m_veloStates{this, "veloStates", {LHCb::State::Location::ClosestToBeam}};
  Gaudi::Property<SLocations> m_longStates{this,
                                           "longStates",
                                           {LHCb::State::Location::ClosestToBeam,
                                            LHCb::State::Location::FirstMeasurement, LHCb::State::Location::BegRich2,
                                            LHCb::State::Location::V0Vertex}};
  Gaudi::Property<SLocations> m_upstreamStates{
      this, "upstreamStates", {LHCb::State::Location::ClosestToBeam, LHCb::State::Location::FirstMeasurement}};
  Gaudi::Property<SLocations> m_downstreamStates{
      this,
      "downstreamStates",
      {LHCb::State::Location::FirstMeasurement, LHCb::State::Location::BegRich2, LHCb::State::Location::V0Vertex}};
  Gaudi::Property<SLocations> m_tStates{
      this, "TTrackStates", {LHCb::State::Location::FirstMeasurement, LHCb::State::Location::BegRich2}};
  Gaudi::Property<SLocations> m_muonStates{this, "muonStates", {}};

  void cleanStates( LHCb::Track& aTrack, LHCb::span<const LHCb::State::Location> loc ) const;

  Gaudi::Property<std::string> m_inputLocation{this, "TracksInContainer", LHCb::TrackLocation::Default};
  Gaudi::Property<bool>        m_storeAllStates{this, "StoreAllStates", false};
};

StatusCode Gaudi::Parsers::parse( std::vector<LHCb::State::Location>& result, const std::string& input ) {
  std::vector<std::string> temp;
  return parse( temp, input ).andThen( [&] {
    result.reserve( temp.size() );
    std::transform( temp.begin(), temp.end(), std::back_inserter( result ), LHCb::State::LocationToType );
  } );
}

DECLARE_COMPONENT( TrackToDST )

StatusCode TrackToDST::execute() {
  // loop
  for ( auto& trk : *get<LHCb::Tracks>( m_inputLocation ) ) {
    // remove the necessary States on the Track
    if ( !m_storeAllStates.value() ) {
      // done in an ugly way for now - will be easier with the new
      // jobOptions parser
      const LHCb::Track::Types type = trk->type();
      switch ( type ) {
      case LHCb::Track::Types::Velo:
        cleanStates( *trk, m_veloStates.value() );
        break;
      case LHCb::Track::Types::Long:
        cleanStates( *trk, m_longStates.value() );
        break;
      case LHCb::Track::Types::Upstream:
        cleanStates( *trk, m_upstreamStates.value() );
        break;
      case LHCb::Track::Types::Downstream:
        cleanStates( *trk, m_downstreamStates.value() );
        break;
      case LHCb::Track::Types::Ttrack:
        cleanStates( *trk, m_tStates.value() );
        break;
      case LHCb::Track::Types::Muon:
        cleanStates( *trk, m_muonStates.value() );
        break;
      default:
        Warning( format( "Unknown track type %i", type ), StatusCode::SUCCESS, 1 ).ignore();
        break;
      } // switch
    }   // if
    // set the appropriate flag!
    trk->setPatRecStatus( LHCb::Track::PatRecStatus::PatRecIDs );
  }

  return StatusCode::SUCCESS;
}

void TrackToDST::cleanStates( LHCb::Track& aTrack, LHCb::span<const LHCb::State::Location> loc ) const {
  if ( msgLevel( MSG::VERBOSE ) ) {
    verbose() << "Analysing Track key=" << aTrack.key() << " type=" << aTrack.type() << " : " << aTrack.states().size()
              << " States at : z =";
    for ( const auto& s : aTrack.states() ) {
      if ( s ) verbose() << " (" << s->z() << " " << s->location() << ")";
    }
    verbose() << endmsg;
  }

  std::vector<LHCb::State*> tempCont;
  for ( const auto& l : loc ) {
    const LHCb::State* state = aTrack.stateAt( l );
    if ( state ) {
      tempCont.push_back( new LHCb::State{*state} );
    } else if ( l != LHCb::State::Location::V0Vertex ) {
      Warning( "Failed to find state - more info in DEBUG", StatusCode::SUCCESS, 1 ).ignore();
      if ( msgLevel( MSG::DEBUG ) ) {
        debug() << "Missing state at " << l << " on track " << aTrack.key() << " of type " << aTrack.type() << endmsg;
      }
    }
  } // loca

  aTrack.clearStates();
  aTrack.addToStates( tempCont );
}
