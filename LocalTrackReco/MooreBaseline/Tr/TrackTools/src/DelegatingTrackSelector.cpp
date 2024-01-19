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
/** @file DelegatingTrackSelector.cpp
 *
 *  Implementation file for reconstruction tool : DelegatingTrackSelector
 *
 *  @author Chris Jones   Christopher.Rob.Jones@cern.ch
 *  @date   30/12/2005
 */

#include "TrackSelectorBase.h"

#include "Kernel/SynchronizedValue.h"

#include "GaudiKernel/SystemOfUnits.h"

#include <shared_mutex>
#include <sstream>

class DelegatingTrackSelector : public TrackSelectorBase {

public:
  using TrackSelectorBase::TrackSelectorBase;

  /** Returns if the given track is selected or not
   *
   *  @param aTrack Reference to the Track to test
   *
   *  @return boolean indicating if the track is selected or not
   *  @retval true  Track is selected
   *  @retval false Track is rejected
   */
  bool accept( const LHCb::Track& aTrack ) const override;

private:
  /// Access on demand the track selector for each track type
  ITrackSelector* trackSelector( const LHCb::Track& aTrack ) const;

  /// Track selector for each track type
  typedef GaudiUtils::HashMap<const LHCb::Track::Types, ITrackSelector*> TrackSelectors;

  /// Cache of TrackSelectors by type, with associated mutex for thread safety
  mutable LHCb::cxx::SynchronizedValue<TrackSelectors, std::shared_mutex> m_trSels;
};

DECLARE_COMPONENT( DelegatingTrackSelector )

bool DelegatingTrackSelector::accept( const LHCb::Track& aTrack ) const {
  if ( msgLevel( MSG::VERBOSE ) ) {
    verbose() << "Trying Track " << aTrack.key() << " " << aTrack.type() << " P=" << aTrack.p() << " Pt=" << aTrack.pt()
              << endmsg;
  }
  return ( checkTrackType( aTrack ) && trackSelector( aTrack )->accept( aTrack ) );
}

ITrackSelector* DelegatingTrackSelector::trackSelector( const LHCb::Track& aTrack ) const {
  // fast path
  ITrackSelector* selector = m_trSels.with_lock( [t = aTrack.type()]( const TrackSelectors& trSels ) {
    auto i = trSels.find( t );
    return i != trSels.end() ? i->second : nullptr;
  } );
  if ( !selector ) {
    // slow path
    selector = m_trSels.with_lock( [this, t = aTrack.type()]( TrackSelectors& trSels ) {
      // double check that we were not overtaken by some other thread
      auto i = trSels.find( t );
      if ( i != trSels.end() ) return i->second;
      // insert new selector in the m_trSels cache
      std::ostringstream name;
      name << toString( t );
      if ( msgLevel( MSG::DEBUG ) ) debug() << "Getting selector for " << name.str() << endmsg;
      if ( name.str().empty() ) Exception( "Empty track type" );
      return trSels[t] = tool<ITrackSelector>( "TrackSelector", name.str(), this );
      // FIXME Note that this way of creating tools without ToolHandles "hides" potential data
      // dependencies that the tool may have to the algorithm. See Issue Rec#171 and the more
      // detailed comment by Sascha in the Merge request Rec!2270
    } );
  }
  return selector;
}
