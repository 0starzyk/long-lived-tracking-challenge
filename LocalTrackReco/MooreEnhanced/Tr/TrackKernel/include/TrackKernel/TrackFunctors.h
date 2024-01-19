/*****************************************************************************\
* (c) Copyright 2018 CERN for the benefit of the LHCb Collaboration           *
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
#include "Event/TrackFitResult.h"
#include "Event/TrackFunctor.h"

namespace LHCb::Event {

  inline namespace v1 {

    inline const TrackFitResult* fitResult( const Track& t ) {
      return dynamic_cast<const TrackFitResult*>( t.fitResult() );
    }

    inline TrackFitResult* fitResult( Track& t ) { return dynamic_cast<TrackFitResult*>( t.fitResult() ); }

    inline auto measurements( const Track& t ) {
      const auto* fr = fitResult( t );
      return fr ? fr->measurements() : decltype( fr->measurements() ){};
    }

    inline auto nodes( const Track& t ) {
      const auto* fr = fitResult( t );
      return fr ? fr->nodes() : decltype( fr->nodes() ){};
    }

    inline auto nMeasurements( const Track& t ) {
      const auto* fr = fitResult( t );
      return fr ? fr->nActiveMeasurements() : 0;
    }

    inline auto nMeasurementsRemoved( const Track& t ) {
      const auto* fr = fitResult( t );
      return fr ? fr->nOutliers() : 0;
    }

    inline const LHCb::State& closestState( const Track& trk, double z ) {
      const auto* fr = fitResult( trk );
      if ( fr && !fr->nodes().empty() ) {
        auto iter = std::min_element( fr->nodes().begin(), fr->nodes().end(), TrackFunctor::distanceAlongZ( z ) );
        if ( iter == fr->nodes().end() ) throw GaudiException( "No state closest to z", __func__, StatusCode::FAILURE );
        return ( *iter )->state();
      } else {
        return trk.closestState( z );
      }
    }

    inline const LHCb::State& closestState( const Track& trk, const Gaudi::Plane3D& plane ) {
      const auto* fr = fitResult( trk );
      if ( fr && !fr->nodes().empty() ) {
        auto iter = std::min_element( fr->nodes().begin(), fr->nodes().end(), TrackFunctor::distanceToPlane( plane ) );
        if ( iter == fr->nodes().end() ) throw GaudiException( "No state closest to z", __func__, StatusCode::FAILURE );
        return ( *iter )->state();
      } else {
        return trk.closestState( plane );
      }
    }
  } // namespace v1
} // namespace LHCb::Event
