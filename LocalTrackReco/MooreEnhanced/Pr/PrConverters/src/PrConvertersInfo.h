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
#include "Event/PrDownstreamTracks.h"
#include "Event/PrFittedForwardTracks.h"
#include "Event/PrLongTracks.h"
#include "Event/PrUpstreamTracks.h"
#include "Event/PrVeloTracks.h"
#include "Event/Track_v1.h"

namespace LHCb::Pr::ConversionInfo {

  struct Downstream {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::Downstream;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrDownstream;
    using Ancestor1                                            = LHCb::Pr::Downstream::Tag::trackSeed;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 1> StateLocations = {
        std::make_pair( LHCb::State::Location::AtTT, 1 )};
    static constexpr bool                       AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 1> AncestorLocations     = {
        "SeedTracksLocation"}; // std::string does not work prior to C++20, as it cannot be constexpr'd
  };

  struct Upstream {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::Upstream;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrVeloUT;
    using Ancestor1                                            = LHCb::Pr::Upstream::Tag::trackVP;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 0> StateLocations = {{}}; // UT does not add
                                                                                                 // states in this case
    static constexpr bool                       AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 1> AncestorLocations     = {"VeloTracksLocation"};
  };

  struct Match {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::Long;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrMatch;
    using Ancestor1                                            = LHCb::Pr::Long::Tag::trackVP;
    using Ancestor2                                            = LHCb::Pr::Long::Tag::trackSeed;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 0> StateLocations        = {{}};
    static constexpr bool                                                 AddStatesFromAncestor = true;
    static constexpr std::array<const char*, 2> AncestorLocations = {"VeloTracksLocation", "SeedTracksLocation"};
  };

  struct Forward {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::Long;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrForward;
    using Ancestor1                                            = LHCb::Pr::Long::Tag::trackVP;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 1> StateLocations = {
        std::make_pair( LHCb::State::Location::AtT, 1 )};
    static constexpr bool                       AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 1> AncestorLocations     = {"VeloTracksLocation"};
  };

  struct ForwardFromVeloUT {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::Long;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrForward;
    using Ancestor1                                            = LHCb::Pr::Long::Tag::trackUT;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 1> StateLocations = {
        std::make_pair( LHCb::State::Location::AtT, 1 )};
    static constexpr bool                       AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 1> AncestorLocations     = {"UpstreamTracksLocation"};
  };

  struct FittedForward {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::Long;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrForward;
    using Ancestor1                                            = LHCb::Pr::Fitted::Forward::Tag::trackSeed;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 0> StateLocations        = {{}};
    static constexpr bool                                                 AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 1>                           AncestorLocations = {"ForwardTracksLocation"};
  };

  struct VeloForward {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::Velo;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrPixel;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 2> StateLocations = {
        std::make_pair( LHCb::State::Location::ClosestToBeam, 0 ), std::make_pair( LHCb::State::Location::EndVelo, 1 )};
    static constexpr bool                       AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 0> AncestorLocations     = {};
  };

  struct VeloBackward {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::VeloBackward;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrPixel;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 1> StateLocations = {
        std::make_pair( LHCb::State::Location::ClosestToBeam, 0 )};
    static constexpr bool                       AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 0> AncestorLocations     = {{}};
  };

  struct Velo {
    static LHCb::Event::v1::Track::Types Type( bool backward = false ) {
      if ( backward )
        return LHCb::Event::v1::Track::Types::VeloBackward;
      else
        return LHCb::Event::v1::Track::Types::Velo;
    }
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrPixel;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 2> StateLocations = {
        std::make_pair( LHCb::State::Location::ClosestToBeam, 0 ), std::make_pair( LHCb::State::Location::EndVelo, 1 )};
    static constexpr bool                       AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 0> AncestorLocations     = {};
  };

  struct Seeding {
    static constexpr LHCb::Event::v1::Track::Types   Type      = LHCb::Event::v1::Track::Types::Ttrack;
    static constexpr LHCb::Event::v1::Track::History PrHistory = LHCb::Event::v1::Track::History::PrSeeding;
    static constexpr std::array<std::pair<LHCb::State::Location, int>, 3> StateLocations = {
        std::make_pair( LHCb::State::Location::FirstMeasurement, 0 ), std::make_pair( LHCb::State::Location::AtT, 1 ),
        std::make_pair( LHCb::State::Location::LastMeasurement, 2 )};
    static constexpr bool                       AddStatesFromAncestor = false;
    static constexpr std::array<const char*, 0> AncestorLocations     = {{}};
  };

} // namespace LHCb::Pr::ConversionInfo
