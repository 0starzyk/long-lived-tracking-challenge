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

/** @class TrackUsedLHCbID TrackUsedLHCbID.h
 *
 * Implementation of TrackUsedLHCbID
 * check if an LHCbID is used
 *
 * @author M.Needham
 * @date   2/08/2006
 *
 * @author M. Schiller
 * @date 2015-02-21
 *  - use BloomFilters to achieve O(1) lookup instead of O(log(nHits))
 */

#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "Kernel/IUsedLHCbID.h"
#include "Kernel/LHCbID.h"
#include "LHCbMath/BloomFilter.h"
#include "TrackInterfaces/ITrackSelector.h"
#include <algorithm>
#include <cstdint>
#include <exception>
#include <string>
#include <vector>

namespace {
  /// calculate number of occupied channels with safety factor
  constexpr uint64_t occCh( uint64_t nch, uint64_t occnumer, uint64_t occdenom, uint64_t safety = 4 ) {
    return ( safety * nch * occnumer ) / occdenom;
  }
  /// calculate BloomFilter capacity for detector with given number of channels and occupancy
  constexpr uint64_t cap( uint64_t nch, uint64_t occnumer, uint64_t occdenom ) {
    return nch < occCh( nch, occnumer, occdenom ) ? nch : occCh( nch, occnumer, occdenom );
  }
} // namespace

class TrackUsedLHCbID : public extends<GaudiTool, IUsedLHCbID, IIncidentListener> {
public:
  /** constructor */
  using extends::extends;

  /** intialize */
  StatusCode initialize() override;

  /** Test if the LHCbID is used
   * @param id to be test tested
   *  @return true if used
   */
  bool used( const LHCb::LHCbID id ) const override;

  /** Implement the handle method for the Incident service.
   *  This is used to nform the tool of software incidents.
   *
   *  @param incident The incident identifier
   */
  void handle( const Incident& incident ) override;

private:
  void initEvent() const;

  typedef std::vector<std::string>     TrackContainers;
  typedef std::vector<ITrackSelector*> Selectors;
  /** Define containers and corresponding selectors in same order.
   *  E.g. inputContainers = "Rec/Track/Forward" and selectorNames = "ForwardSelector".
   */
  Gaudi::Property<TrackContainers> m_inputs{this, "inputContainers", {""}};
  // for track selection
  Gaudi::Property<std::vector<std::string>> m_names{this, "selectorNames", {""}};

  Selectors m_selectors;

  // set maximum number of hits expected in the BloomFilter:
  // if there's more hits in the detector than that, the number of collisions
  // (one hit mistaken for another) will rise above the threshold (1e-4)
  //
  // strategy: max(safety * maxocc * nChannels, nChannels) with safety = 4
  static constexpr uint64_t s_MaxVPHits    = cap( 4100000u, 125u, 100000u );
  static constexpr uint64_t s_MaxUTHits    = cap( 540000u, 18u, 1000u );
  static constexpr uint64_t s_MaxFTHits    = cap( 300000u, 2u, 100u );
  static constexpr uint64_t s_MaxOtherHits = 1024u; // up to 1024 which don't fit above
  static constexpr uint64_t s_denom        = 1u << 20u;
  static constexpr uint64_t s_numer        = 1u * s_denom / 10000u;

  /// flag bits
  enum { Initialized = 1, VP = 16, UT = 32, FT = 64, Other = 128 };
  mutable unsigned m_flags = 0; ///< flags
  // since current and Upgrade detectors are never used in the same job, we
  // can eke out some memory by putting corresponding detectors in a union
  // each, so they "share" the memory by overlapping in the same physical
  // memory location
  mutable BloomFilter<LHCb::LHCbID, s_MaxVPHits, s_numer, s_denom> m_vtx; ///< vertex detector hits
  mutable BloomFilter<LHCb::LHCbID, s_MaxUTHits, s_numer, s_denom> m_bmg; ///< hits in tracking detectors before the
                                                                          ///< magnet
  mutable BloomFilter<LHCb::LHCbID, s_MaxFTHits, s_numer, s_denom> m_amg; ///< hits in tracking detectors after the
                                                                          ///< magnet
  mutable BloomFilter<LHCb::LHCbID, s_MaxOtherHits, s_numer, s_denom> m_otherHits;
};

DECLARE_COMPONENT( TrackUsedLHCbID )

StatusCode TrackUsedLHCbID::initialize() {
  StatusCode sc = GaudiTool::initialize();
  if ( sc.isFailure() ) return sc;

  if ( m_names.value().size() != m_inputs.size() ) {
    if ( m_names.value().size() > m_inputs.size() ) {
      Warning( "More selector names than input locations, discarding excess selectors", StatusCode::FAILURE, 1 )
          .ignore();
    } else {
      Warning( "More input locations than selectors, always accepting tracks where selectors are missing.",
               StatusCode::FAILURE, 1 )
          .ignore();
    }
    m_names.value().resize( m_inputs.size() );
  }

  // make the selector tools
  std::transform(
      std::begin( m_names.value() ), std::end( m_names.value() ), std::back_inserter( m_selectors ),
      [this]( const std::string& name ) { return name.empty() ? nullptr : tool<ITrackSelector>( name, this ); } );

  // make sure we start from empty BloomFilters
  m_vtx.clear(), m_bmg.clear(), m_amg.clear(), m_otherHits.clear();
  m_flags = 0;

  incSvc()->addListener( this, IncidentType::BeginEvent );

  if ( msgLevel( MSG::DEBUG ) ) {
    // printout to announce size of differnt BloomFilters
    debug() << "BloomFilters for hits initialised, sizes: [ VP " << sizeof( m_vtx ) << " UT " << sizeof( m_bmg )
            << " FT " << sizeof( m_amg ) << " Other " << sizeof( m_otherHits ) << "], class size is " << sizeof( *this )
            << endmsg;
  }

  return StatusCode::SUCCESS;
}

bool TrackUsedLHCbID::used( const LHCb::LHCbID id ) const {
  // get the input - seeds
  if ( !( m_flags & Initialized ) ) initEvent();
  switch ( id.detectorType() ) {
  case LHCb::LHCbID::channelIDtype::VP:
    return ( ( ~m_flags ) & VP ) ? false : m_vtx.find( id );
  case LHCb::LHCbID::channelIDtype::UT:
    return ( ( ~m_flags ) & UT ) ? false : m_bmg.find( id );
  case LHCb::LHCbID::channelIDtype::FT:
    return ( ( ~m_flags ) & FT ) ? false : m_amg.find( id );
  default:
    return ( ( ~m_flags ) & Other ) ? false : m_otherHits.find( id );
  };
  return false;
}

void TrackUsedLHCbID::handle( const Incident& incident ) {
  if ( IncidentType::BeginEvent == incident.type() ) {
    // reset Initialized bit to trigger reading new tracks
    m_flags &= ~Initialized;
  }
}

void TrackUsedLHCbID::initEvent() const {
  // only clear the BloomFilters which actually need clearing
  if ( m_flags & VP ) m_vtx.clear();
  if ( m_flags & UT ) m_bmg.clear();
  if ( m_flags & FT ) m_amg.clear();
  if ( m_flags & Other ) m_otherHits.clear();
  m_flags = 0;

  // loop over tracks locations
  auto iterSelector = m_selectors.begin();
  for ( auto iterS = m_inputs.begin(); iterS != m_inputs.end(); ++iterS, ++iterSelector ) {
    // get selection tool
    ITrackSelector* tSelector = *iterSelector;
    // get the containers and extract the ids from the track
    auto tCont = getIfExists<LHCb::Track::Range>( *iterS );
    if ( tCont.empty() ) {
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Track container '" << *iterS << "' does not exist" << endmsg;
      continue;
    }
    for ( auto iterTrack : tCont ) { // loop over tracks in container
      if ( tSelector && !( tSelector->accept( *iterTrack ) ) ) continue;
      // put hits on track into BloomFilters
      for ( const LHCb::LHCbID id : iterTrack->lhcbIDs() ) {
        switch ( id.detectorType() ) {
        case LHCb::LHCbID::channelIDtype::VP:
          m_flags |= VP;
          m_vtx.insert( id );
          break;
        case LHCb::LHCbID::channelIDtype::UT:
          m_flags |= UT;
          m_bmg.insert( id );
          break;
        case LHCb::LHCbID::channelIDtype::FT:
          m_flags |= FT;
          m_amg.insert( id );
          break;
        default:
          m_flags |= Other;
          m_otherHits.insert( id );
          break;
        };
      }
    } // iterTrack
  }   // iterS

  // tracks all read, set Initialized bit
  m_flags |= Initialized;
}
