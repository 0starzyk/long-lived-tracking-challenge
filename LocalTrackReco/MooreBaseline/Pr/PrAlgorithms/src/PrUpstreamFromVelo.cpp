/*****************************************************************************\
* (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "Event/PrUpstreamTracks.h"
#include "Event/PrVeloTracks.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "LHCbAlgs/Transformer.h"

namespace {
  using TracksVP    = LHCb::Pr::Velo::Tracks;
  using TracksUT    = LHCb::Pr::Upstream::Tracks;
  using Transformer = LHCb::Algorithm::Transformer<TracksUT( TracksVP const& )>;
} // namespace

using TracksTag = LHCb::Pr::Upstream::Tag;
namespace Pr {
  /** @class UpstreamFromVelo PrUpstreamFromVelo.cpp
   *
   *  Converts a container of Velo tracks into a container of Upstream ones
   *  with some fixed pT value.
   */
  struct UpstreamFromVelo final : public Transformer {
    UpstreamFromVelo( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer( name, pSvcLocator, {"Input", ""}, {"Output", ""} ) {}

    TracksUT operator()( TracksVP const& inputTracks ) const override {
      using dType = SIMDWrapper::best::types;
      using F     = dType::float_v;
      TracksUT outputTracks{&inputTracks};
      outputTracks.reserve( inputTracks.size() );
      float invAssumedPT{1.f / m_assumedPT};
      for ( auto const& track : inputTracks.simd() ) {
        auto mask = track.loop_mask();
        auto cov  = track.StateCovX( 1 );
        auto pos  = track.StatePos( 1 );
        auto dir  = track.StateDir( 1 );

        // Assign q/p assuming q=+1 and pT is 'AssumedPT'
        auto    txy2 = dir.x * dir.x + dir.y * dir.y;
        F const qop  = invAssumedPT * sqrt( txy2 / ( 1 + txy2 ) );

        auto oTrack = outputTracks.compress_back( mask );
        oTrack.field<TracksTag::trackVP>().set( track.indices() );
        oTrack.field<TracksTag::VPHits>().resize( 0 );
        oTrack.field<TracksTag::UTHits>().resize( 0 );
        oTrack.field<TracksTag::State>().setQOverP( qop );
        oTrack.field<TracksTag::State>().setPosition( pos );
        oTrack.field<TracksTag::State>().setDirection( dir );
        oTrack.field<TracksTag::StateCovX>( LHCb::Pr::Upstream::CovXVector::x_x ).set( cov.x );
        oTrack.field<TracksTag::StateCovX>( LHCb::Pr::Upstream::CovXVector::x_tx ).set( cov.y );
        oTrack.field<TracksTag::StateCovX>( LHCb::Pr::Upstream::CovXVector::tx_tx ).set( cov.z );
      }
      return outputTracks;
    }

    Gaudi::Property<float> m_assumedPT{this, "AssumedPT", 4.5 * Gaudi::Units::GeV};
  };
} // namespace Pr

DECLARE_COMPONENT_WITH_ID( Pr::UpstreamFromVelo, "PrUpstreamFromVelo" )
