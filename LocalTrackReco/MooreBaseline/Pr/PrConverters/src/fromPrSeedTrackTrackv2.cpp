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
#include "Event/PrSeedTracks.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "GaudiKernel/StdArrayAsProperty.h"
#include "LHCbAlgs/ScalarTransformer.h"
#include "LHCbAlgs/Transformer.h"
#include "LHCbMath/MatVec.h"
#include "LHCbMath/bit_cast.h"
#include <vector>

/**
 * Converter between LHCb::Pr::Fitted::Forward::Tracks ( SoA PoD ) and vector<Track_v2>
 * @author Louis Henry
 */
namespace LHCb::Converters::Track::v2 {
  using SeedTag = Pr::Seeding::Tag;

  class fromPrSeedTrackTrackv2
      : public Algorithm::Transformer<std::vector<Event::v2::Track>( Pr::Seeding::Tracks const& )> {

  public:
    using KeyValue = Transformer::KeyValue;

    fromPrSeedTrackTrackv2( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer( name, pSvcLocator, {KeyValue{"InputTracks", ""}}, KeyValue{"OutputTracks", ""} ) {}

    std::vector<Event::v2::Track> operator()( Pr::Seeding::Tracks const& inputTracks ) const override {

      std::vector<Event::v2::Track> outputTracks;
      outputTracks.reserve( inputTracks.size() );
      for ( auto inputTrack : inputTracks.scalar() ) {
        auto& outTrack = outputTracks.emplace_back();
        outTrack.setType( Event::v2::Track::Type::Ttrack );
        outTrack.setHistory( Event::v2::Track::History::PrSeeding );
        outTrack.setPatRecStatus( Event::v2::Track::PatRecStatus::PatRecIDs );

        auto nHits      = ( inputTrack.field<SeedTag::FTHits>() ).size().cast();
        auto chi2PerDoF = ( inputTrack.get<SeedTag::Chi2PerDoF>() ).cast();
        auto getP       = outTrack.setChi2PerDoF( {chi2PerDoF, static_cast<int>( nHits ) - 5} );
        for ( int iID = 0; iID < nHits; iID++ ) outTrack.addToLhcbIDs( inputTrack.ft_lhcbID( iID ).LHCbID() );
        // Copying the states
        for ( int iState = 0; iState < 3; iState++ ) {
          auto                   state = inputTrack.getLHCbState( iState );
          Gaudi::TrackSymMatrix& cov   = state.covariance();
          cov( 0, 0 )                  = m_stateErrorX2;
          cov( 1, 1 )                  = m_stateErrorY2;
          cov( 2, 2 )                  = m_stateErrorTX2;
          cov( 3, 3 )                  = m_stateErrorTY2;
          cov( 4, 4 )                  = m_stateErrorP;
          outTrack.addToStates( state );
        }
      }
      m_nbTracksCounter += outputTracks.size();
      return outputTracks;
    }

  private:
    // - StateErrorX2: Error^2 on x-position (for making Track)
    Gaudi::Property<double> m_stateErrorX2{this, "StateErrorX2", 4.0};
    // - StateErrorY2: Error^2 on y-position (for making Track)
    Gaudi::Property<double> m_stateErrorY2{this, "StateErrorY2", 400.};
    // - StateErrorTX2: Error^2 on tx-slope (for making Track)
    Gaudi::Property<double> m_stateErrorTX2{this, "StateErrorTX2", 6.e-5};
    // - StateErrorTY2: Error^2 on ty-slope (for making Track)
    Gaudi::Property<double> m_stateErrorTY2{this, "StateErrorTY2", 1.e-4};
    // - StateErrorP:  Error^2 on momentum (for making Track)
    Gaudi::Property<double> m_stateErrorP{this, "StateErrorP", 0.15};

    mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of converted Tracks"};
  }; // namespace

  DECLARE_COMPONENT( fromPrSeedTrackTrackv2 )

} // namespace LHCb::Converters::Track::v2
