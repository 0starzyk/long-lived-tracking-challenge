/*****************************************************************************\
 * (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
 *                                                                             *
 * This software is distributed under the terms of the GNU General Public      *
 * Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
 *                                                                             *
 * In applying this licence, CERN does not waive the privileges and immunities *
 * granted to it by virtue of its status as an Intergovernmental Organization  *
 * or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/StateParameters.h"
#include "LHCbAlgs/Transformer.h"
#include "TrackInterfaces/ITrackFitter.h"

/** @class SeedMuonBuilder SeedMuonBuilder.h
 *
 * \brief  Make a SeedMuonTrack: match scifi and muon hits, afterwards apply kalman fit.
 *
 * Parameters:
 * -  MuonTracksLocation: where the muon standalone tracks come from
 * -  SciFiTracksLocation: where the scifi tracks come from
 * -  OutputSeedMuonTracks: where the scifimuontracks go to
 *
 */
//-----------------------------------------------------------------------------
class SeedMuonBuilder final
    : public LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Tracks&, const LHCb::Tracks&,
                                                        DetectorElement const& ),
                                          LHCb::DetDesc::usesConditions<DetectorElement>> {

public:
  // Standard constructor
  SeedMuonBuilder( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"MuonTracksLocation", ""}, KeyValue{"SciFiTracksLocation", ""},
                      KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}},
                     KeyValue{"OutputSeedMuonTracks", "Rec/Track/SeedMuon"} ){};

  LHCb::Tracks operator()( const LHCb::Tracks& muontracks, const LHCb::Tracks& scifitracks,
                           DetectorElement const& lhcb ) const override;

private:
  // properties
  Gaudi::Property<float> m_minSlopeDiff{this, "minSlopeDiff", 0.05};
  // counters
  mutable Gaudi::Accumulators::Counter<> m_countSeedMuons{this, "nSeedMuons"};
  // tools
  ToolHandle<ITrackFitter> m_tracksFitter{this, "Fitter", "TrackMasterFitter/Fitter"};
};

DECLARE_COMPONENT( SeedMuonBuilder )

LHCb::Tracks SeedMuonBuilder::operator()( const LHCb::Tracks& muontracks, const LHCb::Tracks& scifitracks,
                                          DetectorElement const& lhcb ) const {

  LHCb::Tracks outputTracks;
  outputTracks.reserve( muontracks.size() );

  // loop over muons
  for ( auto& muontrack : muontracks ) {
    const auto muonstate = muontrack->closestState( StateParameters::ZEndT );
    const auto muon_tx   = muonstate.tx();
    const auto muon_ty   = muonstate.ty();

    auto min_difference = std::numeric_limits<float>::infinity();
    auto goodCopy       = std::unique_ptr<LHCb::Track>{};

    // loop over scifi tracks
    for ( auto& scifitrack : scifitracks ) {
      const auto scifistate = scifitrack->closestState( StateParameters::ZEndT );
      const auto scifi_tx   = scifistate.tx();
      const auto scifi_ty   = scifistate.ty();
      const auto qp         = scifistate.qOverP();
      const auto cov        = scifistate.covariance();

      float slope_difference = std::sqrt( ( muon_tx - scifi_tx ) * ( muon_tx - scifi_tx ) +
                                          ( muon_ty - scifi_ty ) * ( muon_ty - scifi_ty ) );

      // make aCopy, containing the states from both
      auto aCopy = std::make_unique<LHCb::Track>();
      aCopy->addToAncestors( scifitrack );
      aCopy->addToAncestors( muontrack );
      for ( auto& scifistate : scifitrack->states() ) { aCopy->addToStates( *scifistate ); }
      aCopy->addToLhcbIDs( scifitrack->lhcbIDs() );
      for ( auto& muonstate : muontrack->states() ) { aCopy->addToStates( *muonstate ); }
      aCopy->addToLhcbIDs( muontrack->lhcbIDs() );
      aCopy->firstState().setQOverP( qp );
      for ( auto& state : aCopy->states() ) {
        state->setQOverP( qp );
        state->setCovariance( cov );
      }
      aCopy->setFitStatus( LHCb::Track::FitStatus::Unknown );
      aCopy->setFitHistory( LHCb::Track::FitHistory::Unknown );
      aCopy->setType( LHCb::Track::Types::SeedMuon );
      if ( slope_difference < min_difference ) {
        min_difference = slope_difference;
        goodCopy       = std::move( aCopy );
      }
    } // end of scifitracks loop

    if ( min_difference == std::numeric_limits<float>::infinity() || min_difference > m_minSlopeDiff )
      continue; // nothing reasonable found

    auto sc = m_tracksFitter->operator()( *goodCopy, *lhcb.geometry(), LHCb::Tr::PID::Muon() );
    if ( sc.isFailure() ) continue;

    ++m_countSeedMuons;
    outputTracks.insert( goodCopy.release() );
  } // end of muontrack loop

  return outputTracks;
}
