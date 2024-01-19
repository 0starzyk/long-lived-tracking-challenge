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
#include "GaudiKernel/SharedObjectsContainer.h"
#include <string>
#include <vector>

//-----------------------------------------------------------------------------
// Implementation file for class : CreateFastTrackCollection
//
// 2009-02-25 : Manuel Tobias Schiller <schiller@physi.uni-heidelberg.de>
//-----------------------------------------------------------------------------

/** @class CreateFastTrackCollection CreateFastTrackCollection.h
 * given a list of input track containers, this algorithm creates a fast
 * GaudiSharedObjectsContainer containing pointers to the tracks in the
 * input containers given
 *
 * @author Manuel Tobias Schiller <schiller@physi.uni-heidelberg.de>
 * @date   2009-02-25
 */
class CreateFastTrackCollection : public GaudiAlgorithm {
public:
  /// Standard Constructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm event execution

private:
  Gaudi::Property<std::vector<std::string>> m_inputLocations{this, "InputLocations", {}}; ///< input locations
  Gaudi::Property<std::string>              m_outputLocation{this, "OutputLocation", {}}; ///< output location
  Gaudi::Property<bool> m_slowContainer{this, "SlowContainer", false}; ///< optionally deep-copy tracks into keyed cont.
};

DECLARE_COMPONENT( CreateFastTrackCollection )

//=============================================================================
// Initialization
//=============================================================================
StatusCode CreateFastTrackCollection::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;              // error printed already by GaudiAlgorithm

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Initialize" << endmsg;

  // verify job options
  if ( m_inputLocations.empty() ) {
    error() << "No input locations specified." << endmsg;
    return StatusCode::FAILURE;
  }
  if ( m_outputLocation.empty() ) {
    error() << "No output locations specified." << endmsg;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode CreateFastTrackCollection::execute() {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;
  if ( !m_slowContainer.value() ) {
    // create output container and put it on TES
    SharedObjectsContainer<LHCb::Track>* out = new SharedObjectsContainer<LHCb::Track>;
    put( out, m_outputLocation );
    // get all input containers in turn and put track pointers into output
    for ( const std::string& src : m_inputLocations ) {
      LHCb::Tracks* input = get<LHCb::Tracks>( src );
      out->insert( input->begin(), input->end() );
    }
  } else {
    // count tracks so that we can make an output container of the
    // right size
    std::size_t ntracks = 0;
    for ( const std::string& src : m_inputLocations ) {
      const LHCb::Tracks* input = get<LHCb::Tracks>( src );
      ntracks += input->size();
    }
    // copy tracks the old way, using keyed containers
    LHCb::Tracks* out = new LHCb::Tracks();
    out->reserve( ntracks );
    put( out, m_outputLocation );
    for ( const std::string& src : m_inputLocations ) {
      const LHCb::Tracks* input = get<LHCb::Tracks>( src );
      for ( const LHCb::Track* tr : *input ) {
        // hopefully, I copy everything
        LHCb::Track* ntr = new LHCb::Track( *tr );
        ntr->setHistory( tr->history() );
        ntr->setType( tr->type() );
        out->insert( ntr );
      }
    }
  }

  return StatusCode::SUCCESS;
}
