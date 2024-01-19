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

#include "Detector/FT/FTChannelID.h"
#include "Detector/UT/ChannelID.h"
#include "Detector/VP/VPChannelID.h"
#include "Event/PrLongTracks.h"
#include "Event/PrVeloTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "GaudiKernel/StdArrayAsProperty.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Transformer.h"
#include <vector>

/**
 * Converter between TracksFT SoA PoD and vector<Track_v2>
 *
 * @author Arthur Hennequin (CERN, LIP6)
 */

namespace {
  Gaudi::TrackSymMatrix covariance( float qOverP, std::array<float, 5> covarianceValues ) {
    Gaudi::TrackSymMatrix cov;
    cov( 0, 0 ) = covarianceValues[0];
    cov( 1, 1 ) = covarianceValues[1];
    cov( 2, 2 ) = covarianceValues[2];
    cov( 3, 3 ) = covarianceValues[3];
    cov( 4, 4 ) = covarianceValues[4] * qOverP * qOverP;
    return cov;
  }

} // namespace
class TracksFTConverter : public LHCb::Algorithm::Transformer<std::vector<LHCb::Event::v2::Track>(
                              const std::vector<LHCb::Event::v2::Track>&, const LHCb::Pr::Long::Tracks& )> {
  using Track  = LHCb::Event::v2::Track;
  using Tracks = LHCb::Pr::Long::Tracks;
  // From PrGeometryTool in PrAlgorithms

public:
  TracksFTConverter( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"TracksUTLocation", "Rec/Track/v2/UT"}, KeyValue{"TracksFTLocation", "Rec/Track/FT"}},
                     KeyValue{"OutputTracksLocation", "Rec/Track/v2/FT"} ) {}

  Gaudi::Property<std::array<float, 5>> m_covarianceValues{this, "covarianceValues", {4.0, 400.0, 4.e-6, 1.e-4, 0.1}};

  std::vector<Track> operator()( const std::vector<Track>& tracksUT, const Tracks& tracksFT ) const override {
    std::vector<Track> out;
    out.reserve( tracksFT.size() );
    m_nbTracksCounter += tracksFT.size();

    using dType = SIMDWrapper::scalar::types;
    using F     = dType::float_v;

    for ( auto const& track : tracksFT.scalar() ) {
      auto  hasUT    = track.trackUT().cast();
      auto& trackUT  = hasUT >= 0 ? tracksUT[hasUT] : tracksUT[track.trackVP().cast()];
      auto& newTrack = out.emplace_back( trackUT );
      newTrack.addToAncestors( trackUT );

      // set q/p in all of the existing states
      auto const qop     = track.qOverP().cast();
      auto const errQop2 = m_covarianceValues[4] * qop * qop;

      for ( auto& state : newTrack.states() ) {
        state.setQOverP( qop );
        state.setErrQOverP2( errQop2 );
      }

      // Add state end SciFi : 1
      LHCb::State       state;
      LHCb::StateVector s;
      Vec3<F>           pos = track.StatePos( 1 );
      Vec3<F>           dir = track.StateDir( 1 );
      s.setX( pos.x.cast() );
      s.setY( pos.y.cast() );
      s.setZ( pos.z.cast() );
      s.setTx( dir.x.cast() );
      s.setTy( dir.y.cast() );
      s.setQOverP( qop );
      state.setState( s );
      state.setCovariance( covariance( qop, m_covarianceValues ) );
      state.setLocation( LHCb::State::Location::AtT );

      newTrack.addToStates( state );

      // Add LHCbIds
      newTrack.setLhcbIDs( track.lhcbIDs(), LHCb::Tag::Unordered );

      newTrack.setType( Track::Type::Long );
      newTrack.setHistory( Track::History::PrForward );
      newTrack.setPatRecStatus( Track::PatRecStatus::PatRecIDs );
    }

    return out;
  };

private:
  mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
};

DECLARE_COMPONENT( TracksFTConverter )
