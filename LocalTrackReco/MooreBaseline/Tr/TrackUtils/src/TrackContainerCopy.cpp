/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

/** @class TrackContainerCopy TrackContainerCopy.h
 *
 *  Copy a container of tracks. By default do not copy tracks that failed the fit
 *
 *  Properties:
 *
 *  - inputLocations: Vector of input locations to copy.
 *  - outputLocation: Output location to copy the tracks to.
 *  - copyFailures: Also copy tracks that are flagged invalid?
 *  - Selector: The selector to select a subsample of tracks to copy (e.g.  TrackSelector )
 *
 *  @author M.Needham
 *  @date   30/05/2006
 */

#include "Event/Track.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/MergingTransformer.h"
#include "TrackInterfaces/ITrackSelector.h"
#include <string>

class TrackContainerCopy final : public LHCb::Algorithm::MergingTransformer<LHCb::Tracks(
                                     const Gaudi::Functional::vector_of_const_<LHCb::Track::Range>& )> {

public:
  TrackContainerCopy( const std::string& name, ISvcLocator* pSvcLocator )
      : MergingTransformer( name, pSvcLocator, {"inputLocations", {LHCb::TrackLocation::Velo}},
                            {"outputLocation", LHCb::TrackLocation::Default} ) {}

  LHCb::Tracks operator()( const Gaudi::Functional::vector_of_const_<LHCb::Track::Range>& trackLists ) const override {
    LHCb::Tracks outCont;
    for ( const auto& trackList : trackLists ) {
      for ( const auto* track : trackList ) {
        if ( ( !track->checkFlag( LHCb::Track::Flags::Invalid ) || m_copyFailures.value() ) &&
             ( !m_selector.isEnabled() || m_selector->accept( *track ) ) ) {
          outCont.insert( new LHCb::Track( *track ) );
        }
      }
    }
    return outCont;
  }

private:
  Gaudi::Property<bool>      m_copyFailures{this, "copyFailures", false}; ///< If true, copy also tracks that failed fit
  ToolHandle<ITrackSelector> m_selector{this, "Selector", ""};
};

DECLARE_COMPONENT( TrackContainerCopy )
