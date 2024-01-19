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
/** @file VPTrackSelector.cpp
 *
 *  Implementation file for reconstruction tool VPTrackSelector,
 *  based on VeloTrackSelector.
 *
 *  @author Christoph Hombach Christoph.Hombach@hep.manchester.ac.uk
 */

#include "TrackSelector.h"

class VPTrackSelector : public TrackSelector {

public:
  /// Constructor
  using TrackSelector::TrackSelector;

  /** Returns if the given track is selected or not
   *
   *  @param track Reference to the track to test
   *
   *  @return boolean indicating if the track is selected or not
   *  @retval true  Track is selected
   *  @retval false Track is rejected
   */
  bool accept( const LHCb::Track& track ) const override;

private:
  Gaudi::Property<size_t> m_requireTileOverlap{this, "RequireTileOverlap", false};
};

DECLARE_COMPONENT( VPTrackSelector )

//-----------------------------------------------------------------------------

bool VPTrackSelector::accept( const LHCb::Track& track ) const {
  bool acceptthis = true;
  if ( m_requireTileOverlap ) {
    acceptthis  = false;
    auto id     = track.lhcbIDs().begin();
    auto previd = id;
    for ( ++id; id < track.lhcbIDs().end() && !acceptthis; ++id )
      if ( id->isVP() ) {
        acceptthis = previd->isVP() && id->vpID().module() == previd->vpID().module() &&
                     id->vpID().sensor() != previd->vpID().sensor();
        previd = id;
      }
  }
  return acceptthis && TrackSelector::accept( track );
}
