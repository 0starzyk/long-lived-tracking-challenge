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
#include "LHCbAlgs/Transformer.h"

/**
 * This algorithm is a dummy std::vector<LHCb::Track> to keyed container converter, to allow the truth matching to work
 * for the upgrade
 *
 *  @author Renato Quagliani
 *  @date   25-01-2018
 */
class PrTrackConverter : public LHCb::Algorithm::Transformer<LHCb::Tracks( const std::vector<LHCb::Track>& )> {

public:
  PrTrackConverter( const std::string& name, ISvcLocator* pSvcLocator );
  LHCb::Tracks operator()( const std::vector<LHCb::Track>& inputTracks ) const override;
};

DECLARE_COMPONENT( PrTrackConverter )

PrTrackConverter::PrTrackConverter( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator, KeyValue{"InputTracksLocation", ""}, KeyValue{"OutKeyedTrackLocation", ""} ) {}

LHCb::Tracks PrTrackConverter::operator()( const std::vector<LHCb::Track>& inputTracks ) const {
  // Loop over the Tracks
  LHCb::Tracks OutputTracks;
  for ( auto& Source_Track : inputTracks ) {
    LHCb::Track* tr = new LHCb::Track;
    // Copy the track content in the new container
    tr->copy( Source_Track );
    OutputTracks.insert( tr );
  } // End loop over Tracks
  return OutputTracks;
}
