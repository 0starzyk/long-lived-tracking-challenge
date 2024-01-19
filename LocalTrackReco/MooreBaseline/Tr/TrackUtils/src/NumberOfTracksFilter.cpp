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
#include "Event/Track.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include <numeric>

//-----------------------------------------------------------------------------
// Implementation file for class : NumberOfTracksFilter
//
// 2008-03-05 : Patrick Koppenburg
//-----------------------------------------------------------------------------

/** @class NumberOfTracksFilter NumberOfTracksFilter.h
 *
 *  Looks into TES locations foir tracks and returns filterPassed
 *  true or false depending on the numebr of tracks.
 *  Default is at least 0 tracks (i.e. do nothing)
 *
 *  @code
 *  NumberOfTracksFilter.TrackLocations = { "Rec/Track/Best" };
 *  NumberOfTracksFilter.MinTracks = 2 ;
 *  NumberOfTracksFilter.MaxTracks = 1000 ;
 *  @endcode
 *
 *  This will require there are at least 2 and at most 1000 tracks in
 *  Rec/Track/Best to constinue the sequence. MaxTracks = -1 (default)
 *  does nothing.
 *
 *  @author Patrick Koppenburg
 *  @date   2008-03-05
 */
class NumberOfTracksFilter : public GaudiAlgorithm {
public:
  /// Standard constructor
  NumberOfTracksFilter( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution

private:
  Gaudi::Property<std::vector<std::string>> m_tracksPath{this, "TrackLocations", {}}; ///< locations
  Gaudi::Property<int>                      m_minTracks{this, "MinTracks", 0};        ///< min number of tracks
  Gaudi::Property<int>                      m_maxTracks{this, "MaxTracks", -1};       ///< max number of tracks
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( NumberOfTracksFilter )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
NumberOfTracksFilter::NumberOfTracksFilter( const std::string& name, ISvcLocator* pSvcLocator )
    : GaudiAlgorithm( name, pSvcLocator ) {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode NumberOfTracksFilter::initialize() {
  return GaudiAlgorithm::initialize().andThen( [&] {
    if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Initialize" << endmsg;
    if ( m_tracksPath.empty() ) {
      if ( context() != "HLT" ) m_tracksPath.value().push_back( LHCb::TrackLocation::Default );
    }
    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << "Tracks will be taken from ";
      for ( const auto& t : m_tracksPath ) debug() << t << " ";
      debug() << endmsg;
      debug() << "Will require at least " << m_minTracks << " tracks" << endmsg;
      if ( m_maxTracks > -1 ) debug() << "Will require at most " << m_maxTracks << " tracks" << endmsg;
    }
  } );
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode NumberOfTracksFilter::execute() {

  if ( msgLevel( MSG::VERBOSE ) ) verbose() << "==> Execute" << endmsg;

  int nT = std::accumulate( m_tracksPath.begin(), m_tracksPath.end(), 0, [&]( int n, const std::string& p ) {
    LHCb::Track::Container* inTracks = getIfExists<LHCb::Track::Container>( p );
    if ( !inTracks ) {
      Warning( "No tracks at " + p, StatusCode::SUCCESS ).ignore();
    } else {
      if ( msgLevel( MSG::VERBOSE ) )
        verbose() << "Container " << p << " contains " << inTracks->size() << " Tracks" << endmsg;
      n += inTracks->size();
    }
    return n;
  } );

  if ( msgLevel( MSG::DEBUG ) ) debug() << "Found " << nT << " tracks" << endmsg;

  if ( nT < m_minTracks )
    setFilterPassed( false ); // bad
  else if ( ( m_maxTracks > -1 ) & ( nT > m_maxTracks ) )
    setFilterPassed( false ); // bad
  else
    setFilterPassed( true ); // good

  return StatusCode::SUCCESS;
}
