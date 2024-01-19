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
#include "Magnet/DeMagnet.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackInterfaces/ITrackFitter.h"
#include "TrackInterfaces/ITrackMomentumEstimate.h"

/** @class VeloMuonBuilder VeloMuonBuilder.h
 *
 * \brief  Make a ValoMuonTrack: match velo and muon tracks. afterwards apply kalman fit.
 *
 * Parameters:
 * -  MuonTracksLocation: where the muon standalone tracks come from
 * -  VeloTracksLocation: where the velo tracks come from
 * -  OutputVeloMuonTracks: where the velomuontracks go to
 * -  zmagnet: where the magnet's bending plane is (optimised for matching)
 * -  m_distcut: maximum distance between velo and muon tracks depending on the region
 *
 */
//-----------------------------------------------------------------------------
class VeloMuonBuilder final
    : public LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Tracks&, const LHCb::Tracks&,
                                                        DetectorElement const&, DeMagnet const& ),
                                          LHCb::DetDesc::usesConditions<DetectorElement, DeMagnet>> {

public:
  // Standard constructor
  VeloMuonBuilder( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"MuonTracksLocation", ""}, KeyValue{"VeloTracksLocation", ""},
                      KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top},
                      KeyValue{"MagnetLocation", LHCb::Det::Magnet::det_path}},
                     KeyValue{"OutputVeloMuonTracks", "Rec/Track/VeloMuon"} ){};

  LHCb::Tracks operator()( const LHCb::Tracks& muontracks, const LHCb::Tracks& velotracks, DetectorElement const& lhcb,
                           DeMagnet const& magnet ) const override;

private:
  mutable Gaudi::Accumulators::Counter<> m_countVeloMuons{this, "nVeloMuons"};

  Gaudi::Property<float>  m_distcutmultiplyer{this, "cutScale", 1};
  Gaudi::Property<double> m_zmagnet{this, "zmagnet", 5400.f * Gaudi::Units::mm};
  // determination on private ntuple from Run 2
  Gaudi::Property<std::vector<float>> m_distcut{this, "DistanceCut", {30 * 30, 60 * 60, 110 * 110, 200 * 200}};
  Gaudi::Property<std::vector<float>> s_xscale{this, "xScale", {0.06f, 0.1f, 0.15f, 0.15f}};

  Gaudi::Property<float> m_ptkickConstant{this, "pTkick", 1265.f * Gaudi::Units::MeV};
  Gaudi::Property<bool>  m_tCubicFit{this, "tCubicFit", true};

  ToolHandle<ITrackExtrapolator>     m_linearextrapolator{this, "LinearExtrapolator", "TrackLinearExtrapolator"};
  ToolHandle<ITrackFitter>           m_tracksFitter{this, "Fitter", "TrackMasterFitter/Fitter"};
  ToolHandle<ITrackMomentumEstimate> m_fastmomestimate{this, "FastMomEstimate", "FastMomentumEstimate"};
};

DECLARE_COMPONENT( VeloMuonBuilder )

LHCb::Tracks VeloMuonBuilder::operator()( const LHCb::Tracks& muontracks, const LHCb::Tracks& velotracks,
                                          DetectorElement const& lhcb, DeMagnet const& magnet ) const {

  LHCb::Tracks outputTracks;
  outputTracks.reserve( muontracks.size() );

  // loop over muons
  for ( auto& muontrack : muontracks ) {

    // extrapolate the muon track to the position in the magnet
    Gaudi::XYZPoint chamber = muontrack->position();
    Gaudi::XYZPoint muonpunktx, muonpunkty;
    muonpunkty = chamber;
    auto sc    = m_linearextrapolator->position( *( muontrack ), m_zmagnet, muonpunktx, *lhcb.geometry(),
                                              LHCb::Tr::PID::Muon() );
    if ( sc.isFailure() ) continue;

    auto minweight = std::numeric_limits<float>::infinity();
    // default using hit from M4; in case there is only 3 hits (SecondLoop = True), it can be a hit from M5
    const unsigned int MuonHitFromM4 = 2;
    auto               reg           = muontrack->lhcbIDs()[MuonHitFromM4].muonID().region();

    auto goodCopy = std::unique_ptr<LHCb::Track>();
    // loop over velos - chose the one with the closest distance at the extrapolated position in the magnet
    for ( auto& velotrack : velotracks ) {

      if ( velotrack->isVeloBackward() ) continue;

      // extrapolate the velo track to the position in the magnet
      Gaudi::XYZPoint velopunktx, velopunkty;
      sc = m_linearextrapolator->position( *velotrack, m_zmagnet, velopunktx, *lhcb.geometry(), LHCb::Tr::PID::Muon() );
      if ( sc.isFailure() ) continue;
      sc = m_linearextrapolator->position( *velotrack, chamber.z(), velopunkty, *lhcb.geometry(),
                                           LHCb::Tr::PID::Muon() );
      if ( sc.isFailure() ) continue;

      // now calculate distance
      float weighteddistance =
          float( ( velopunktx.x() - muonpunktx.x() ) * ( velopunktx.x() - muonpunktx.x() ) * s_xscale[reg] +
                 ( 1 - s_xscale[reg] ) * ( velopunkty.y() - muonpunkty.y() ) * ( velopunkty.y() - muonpunkty.y() ) );

      auto distancecut = m_distcut[reg] * m_distcutmultiplyer;
      if ( weighteddistance > distancecut ) { continue; }

      LHCb::State veloState = velotrack->firstState();
      double      xkick     = (double)( chamber.x() - veloState.x() ); // jstefaniak used interpolated value here

      double qp            = double( xkick / m_ptkickConstant / ( (double)chamber.z() - m_zmagnet ) );
      qp                   = qp * magnet.signedRelativeCurrent();
      double     sigmaqp   = qp * 0.15;
      const auto muonState = muontrack->closestState( StateParameters::ZEndT );
      sc                   = m_fastmomestimate->calculate( magnet, &veloState, &muonState, qp, sigmaqp, m_tCubicFit );
      Gaudi::TrackSymMatrix cov;
      cov( 0, 0 ) = 1.f;
      cov( 1, 1 ) = 1.f;
      cov( 2, 2 ) = 1.f;
      cov( 3, 3 ) = 1.f;
      cov( 4, 4 ) = sigmaqp * sigmaqp;

      auto aCopy = std::make_unique<LHCb::Track>();
      aCopy->addToAncestors( velotrack );
      aCopy->addToAncestors( muontrack );
      for ( auto& velostate : velotrack->states() ) { aCopy->addToStates( *velostate ); }
      aCopy->addToLhcbIDs( velotrack->lhcbIDs() );
      for ( auto& muonstate : muontrack->states() ) { aCopy->addToStates( *muonstate ); }
      aCopy->addToLhcbIDs( muontrack->lhcbIDs() );
      aCopy->firstState().setQOverP( qp );
      for ( auto& state : aCopy->states() ) {
        state->setQOverP( qp );
        state->setCovariance( cov );
      }
      aCopy->setFitStatus( LHCb::Track::FitStatus::Unknown );
      aCopy->setFitHistory( LHCb::Track::FitHistory::Unknown );
      aCopy->setType( LHCb::Track::Types::VeloMuon );

      float weight = weighteddistance;
      // safe the new track in goodCopy if the distance between velo and muon is so far the clostest
      if ( weight < minweight ) {
        minweight = weight;
        goodCopy  = std::move( aCopy );
      }
    }

    if ( minweight == std::numeric_limits<float>::infinity() ) continue; // -- nothing was found

    sc = m_tracksFitter->operator()( *goodCopy, *lhcb.geometry(), LHCb::Tr::PID::Muon() );
    if ( sc.isFailure() ) continue;
    ++m_countVeloMuons;
    outputTracks.insert( goodCopy.release() );
  }

  return outputTracks;
}
