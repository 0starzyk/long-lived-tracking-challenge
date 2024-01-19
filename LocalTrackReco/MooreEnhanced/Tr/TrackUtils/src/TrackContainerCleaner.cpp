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
/** @class TrackContainerCleaner TrackContainerCleaner.h
 *
 *  Clean out tracks with some criteria from the container
 *
 *  @author M.Needham
 *  @date   30/05/2006
 */

#include "Event/Track.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "TrackInterfaces/ITrackSelector.h"
#include <string>
#include <vector>

using namespace LHCb;

class TrackContainerCleaner : public GaudiAlgorithm {

public:
  // Constructors and destructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode execute() override;

private:
  Gaudi::Property<std::string>           m_inputLocation{this, "inputLocation", TrackLocation::Default};
  ToolHandle<ITrackSelector>             m_selector{this, "Selector", "TrackSelector/Selector"};
  mutable Gaudi::Accumulators::Counter<> m_in{this, "# input tracks"};
  mutable Gaudi::Accumulators::Counter<> m_removed{this, "# removed tracks"};
};

DECLARE_COMPONENT( TrackContainerCleaner )

StatusCode TrackContainerCleaner::execute() {

  Tracks* trackCont = getIfExists<Tracks>( m_inputLocation );
  if ( !trackCont ) {
    return Warning( "Input container " + m_inputLocation + " does not exist", StatusCode::SUCCESS, 20 );
  }

  // loop and select bad tracks
  std::vector<Track*> tVec;
  std::copy_if( trackCont->begin(), trackCont->end(), std::back_inserter( tVec ),
                [&]( const Track* trk ) { return !m_selector->accept( *trk ); } );

  // remove from the container and delete the bad tracks
  std::size_t nRemoved = 0;
  for ( auto iter = tVec.rbegin(); iter != tVec.rend(); ++iter ) {
    trackCont->erase( *iter );
    ++nRemoved;
  }

  m_in += trackCont->size();
  m_removed += nRemoved;

  return StatusCode::SUCCESS;
}
