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
#include "TrackInterfaces/IGhostProbability.h"

//-----------------------------------------------------------------------------
// Implementation file for class : TrackAddNNGhostId
//
// 2009-10-06 : Johannes Albrecht
//-----------------------------------------------------------------------------
using namespace LHCb;

/** @class TrackAddNNGhostId TrackAddNNGhostId.h
 *
 *
 *  @author Johannes Albrecht
 *  @date   2009-10-06
 */

class TrackAddNNGhostId : public GaudiAlgorithm {
public:
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode execute() override; ///< Algorithm execution

private:
  ToolHandle<IGhostProbability> m_ghostTool{this, "GhostIdTool", "UpgradeGhostId"};
  Gaudi::Property<std::string>  m_inputLocation{this, "inputLocation", TrackLocation::Default};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TrackAddNNGhostId )

//=============================================================================
// Main execution
//=============================================================================
StatusCode TrackAddNNGhostId::execute() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  Tracks* inCont = getIfExists<Tracks>( m_inputLocation );
  if ( !inCont ) {
    if ( msgLevel( MSG::DEBUG ) ) debug() << "no tracks at " << m_inputLocation << endmsg;
    return StatusCode::SUCCESS;
  }

  for ( auto& t : *inCont ) m_ghostTool->execute( *t ).ignore();

  return StatusCode::SUCCESS;
}

//=============================================================================
