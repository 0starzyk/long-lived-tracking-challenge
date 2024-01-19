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
#pragma once

#include "Event/Track.h"

#include "GaudiKernel/IAlgTool.h"

#include <vector>

#include "GaudiKernel/IAlgTool.h"
#include "IHitExpectation.h"

/** @class IHitExpectation IDetailedHitExpectation.h TrackInterfaces/IDetailedHitExpectation
 *
 *  How many hits (of a given type) do we expect on a track ?
 *
 *  @author M. Schiller
 *  @date   2022-04-12
 */
struct IDetailedHitExpectation : extend_interfaces<IAlgTool> {
  DeclareInterfaceID( IDetailedHitExpectation, 1, 0 );

  /** small struct returning detailed hit info....
   */
  struct DetailedLayerInfo {
    LHCb::LHCbID id{};          //                   = 0;     ///< detector-specific ID identifying layer/sensor/...
    unsigned int global_id{};   // detector specific ID identifying a group of layer/sensor/...
    bool         found = false; ///< bool indicating if hit was found in layer/sensor/...
    constexpr DetailedLayerInfo( LHCb::LHCbID _id, unsigned int _global_id, bool _found )
        : id( _id ), global_id( _global_id ), found( _found ) {}
  };
  using DetailedInfo = std::vector<DetailedLayerInfo>;

  /** Returns number of hits expected
   *
   *  @param aTrack Reference to the Track to test
   *
   *  @return Info info
   */
  virtual DetailedInfo detailedExpectation( const LHCb::Track& aTrack, IGeometryInfo const& geometry ) const = 0;
};
