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
#include "AIDA/IHistogram1D.h"
#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/Measurement.h"
#include "Event/State.h"
#include "Event/StateVector.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiHistoTool.h"
#include "GaudiAlg/IHistoTool.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiUtils/HistoStats.h"
#include "LHCbAlgs/Consumer.h"
#include "Linker/LinkedTo.h"
#include "TrackCheckerBase.h"
#include "TrackInterfaces/ITrackProjector.h"
#include "TrackInterfaces/ITrackProjectorSelector.h"
#include "TrackKernel/TrackFunctors.h"
#include <map>
#include <string>

/** @class TrackResChecker TrackResChecker.h
 *
 * Class for track monitoring
 *  @author M. Needham.
 *  @date   6-5-2007
 */

class TrackResChecker
    : public LHCb::Algorithm::Consumer<void( LHCb::Track::Range const&, LHCb::MCParticles const&,
                                             LHCb::LinksByKey const&, DetectorElement const& ),
                                       LHCb::DetDesc::usesBaseAndConditions<TrackCheckerBase, DetectorElement>> {

public:
  /** Standard constructor */
  TrackResChecker( const std::string& name, ISvcLocator* pSvcLocator );

  /** Algorithm initialize */
  StatusCode initialize() override;

  /** Algorithm execute */
  void operator()( LHCb::Track::Range const&, LHCb::MCParticles const&, LHCb::LinksByKey const&,
                   DetectorElement const& ) const override;

  /** Algorithm finalize */
  StatusCode finalize() override;

private:
  void resolutionHistos( IHistoTool const& histotool, LHCb::Track const& track, LHCb::MCParticle const& mcPart,
                         IGeometryInfo const& ) const;

  void pullplots( const IHistoTool& histotool, const LHCb::State& trueState, const LHCb::State& recState,
                  const std::string& location ) const;

  void plotsByMeasType( IHistoTool const& histotool, LHCb::Track const& track, LHCb::MCParticle const& mcPart,
                        IGeometryInfo const& ) const;

  const IHistoTool* createHistoTool( const std::string& name ) const;

private:
  Gaudi::Property<bool> m_plotsByMeasType{this, "PlotsByMeasType", false};

  ToolHandle<ITrackProjectorSelector>      m_projectorSelector{this, "ProjectorSelector",
                                                          "TrackProjectorSelector/Projector"};
  typedef std::map<int, const IHistoTool*> HistoToolMap;
  HistoToolMap                             m_histoTools;
};

DECLARE_COMPONENT( TrackResChecker )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TrackResChecker::TrackResChecker( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator,
                {KeyValue{"TracksInContainer", LHCb::TrackLocation::Default},
                 KeyValue{"MCParticleInContainer", LHCb::MCParticleLocation::Default},
                 KeyValue{"LinkerInTable", "Link/" + LHCb::TrackLocation::Default},
                 KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}} ) {}

StatusCode TrackResChecker::initialize() {
  return Consumer::initialize().andThen( [&] {
    m_histoTools[0] = createHistoTool( "ALL" );
    if ( splitByType() ) {
      using Type           = LHCb::Track::Types;
      constexpr auto types = std::array{Type::Velo,       Type::VeloBackward, Type::Long, Type::Upstream,
                                        Type::Downstream, Type::Ttrack,       Type::Muon};
      for ( auto type : types ) { m_histoTools[static_cast<int>( type )] = createHistoTool( toString( type ) ); }
    }
  } );
}

const IHistoTool* TrackResChecker::createHistoTool( const std::string& name ) const {
  IHistoTool*     htool  = tool<IHistoTool>( "HistoTool", name, this );
  GaudiHistoTool* ghtool = dynamic_cast<GaudiHistoTool*>( htool );
  ghtool->setHistoTopDir( histoPath() + "/" );
  std::string histodir = ghtool->histoDir();
  size_t      pos      = histodir.find( '.' );
  if ( pos != std::string::npos ) histodir.erase( 0, pos + 1 );
  ghtool->setHistoDir( histodir );
  return htool;
}

//=============================================================================
// Execute
//=============================================================================
void TrackResChecker::operator()( LHCb::Track::Range const& tracks, LHCb::MCParticles const& mcParts,
                                  LHCb::LinksByKey const& links, DetectorElement const& lhcb ) const {

  // loop over them
  for ( const auto* track : tracks ) {
    // Get the associated true particle
    const LHCb::MCParticle* mcparticle = mcTruth( *track, mcParts, links );

    if ( mcparticle
         // we actually just want to know if it passes the 'selector' inside the IMCReconstructible
         // && int(selector()->reconstructible(mcparticle)) > int(IMCReconstructible::NotReconstructible)
    ) {
      // split by type..
      const IHistoTool* histotool( 0 );
      if ( splitByType() ) {
        histotool = m_histoTools.at( static_cast<int>( track->type() ) );
      } else {
        histotool = m_histoTools.at( 0 );
      }

      // resolutions at predefined z.
      resolutionHistos( *histotool, *track, *mcparticle, *lhcb.geometry() );

      // prob chi^2
      histotool->plot1D( track->chi2PerDoF(), "chi2PerDof", "chi2PerDof", 0., 100., 1000 );
      histotool->plot1D( track->probChi2(), "probChi2", "probChi2", 0., 1., 50 );

      // fit status
      histotool->plot1D( static_cast<int>( track->fitStatus() ), "fitStatus", "fit status", -0.5, 4.5, 5 );

      histotool->plot1D( mcparticle->p() / Gaudi::Units::GeV, "truemom", "true p [GeV]", 0, 50, 100 );
      histotool->plot1D( mcparticle->pt() / Gaudi::Units::GeV, "truept", "true pT [GeV]", 0, 10, 100 );

      // Resolutions and pulls per Measurement type
      if ( m_plotsByMeasType.value() && nMeasurements( *track ) > 0 )
        plotsByMeasType( *histotool, *track, *mcparticle, *lhcb.geometry() );
    }
  }
}

//=============================================================================
//
//=============================================================================
void TrackResChecker::resolutionHistos( IHistoTool const& htool, LHCb::Track const& track,
                                        LHCb::MCParticle const& mcPart, IGeometryInfo const& geometry ) const {
  // pulls at vertex
  LHCb::State trueStateVertex;
  idealStateCreator()->createStateVertex( &mcPart, trueStateVertex ).ignore();
  LHCb::State vtxState;
  StatusCode  sc = extrapolator()->propagate( track, trueStateVertex.z(), vtxState, geometry );
  if ( sc.isSuccess() ) pullplots( htool, trueStateVertex, vtxState, "vertex" );

  // for vertex also make some 2-d plots
  if ( track.type() == LHCb::Track::Types::Long || track.type() == LHCb::Track::Types::Upstream ||
       track.type() == LHCb::Track::Types::Downstream || track.type() == LHCb::Track::Types::Ttrack ) {
    const double invp  = std::abs( track.firstState().qOverP() );
    const double ptrue = mcPart.p();
    const double eta   = mcPart.pseudoRapidity(); // track.pseudoRapidity();

    htool.plot2D( ptrue / Gaudi::Units::GeV, invp * ptrue - 1, "vertex/dpoverp_vs_p", "dp/p vs p", 0., 50., -0.1, 0.1,
                  25, 50 );
    htool.plot2D( eta, invp * ptrue - 1, "vertex/dpoverp_vs_eta", "dp/p vs eta", 2., 5., -0.05, 0.05, 20, 50 );

    const double invperr2 = track.firstState().covariance()( 4, 4 );
    if ( invperr2 > 0 ) {
      const double ppull = ( invp - 1 / ptrue ) / std::sqrt( invperr2 );
      htool.plot2D( ptrue / Gaudi::Units::GeV, ppull, "vertex/p_pull_vs_p", "p pull vs p", 0., 50., -10., 10., 25, 50 );
      htool.plot2D( eta, ppull, "vertex/p_pull_vs_eta", "p pull vs eta", 2., 5., -10., 10., 20, 50 );
    }
  }

  // fraction of tracks with correct charge
  bool correctcharge = track.firstState().qOverP() * mcPart.particleID().threeCharge() > 0;
  htool.plot1D( correctcharge, "correctcharge", "correct charge", -0.5, 1.5, 2 );

  if ( fullDetail() ) {
    for ( const LHCb::State* state : track.states() ) {
      // skip the closest to beam, since we already have it
      if ( state->location() == LHCb::State::Location::ClosestToBeam ) continue;
      double state_z = state->z();
      if ( state_z > mcPart.originVertex()->position().Z() && state_z < mcPart.endVertices().back()->position().Z() ) {
        LHCb::State trueState;
        StatusCode  sc = idealStateCreator()->createState( &mcPart, state_z, trueState, geometry );
        if ( sc.isSuccess() ) {
          std::string location = state->location() != LHCb::State::Location::LocationUnknown
                                     ? Gaudi::Utils::toString( state->location() )
                                     : format( "state_%d_mm", int( state_z ) );
          pullplots( htool, trueState, *state, location );
        }
      }
    }
  }
}

//=============================================================================
//
//=============================================================================
void TrackResChecker::pullplots( const IHistoTool& htool, const LHCb::State& trueState, const LHCb::State& recState,
                                 const std::string& location ) const {

  // save some typing
  const Gaudi::TrackVector&    vec     = recState.stateVector();
  const Gaudi::TrackVector&    trueVec = trueState.stateVector();
  const Gaudi::TrackSymMatrix& cov     = recState.covariance();
  const Gaudi::TrackSymMatrix& trueCov = trueState.covariance();
  const double                 dx      = vec( 0 ) - trueVec( 0 );
  const double                 dy      = vec( 1 ) - trueVec( 1 );
  const double                 dtx     = vec( 2 ) - trueVec( 2 );
  const double                 dty     = vec( 3 ) - trueVec( 3 );

  // fill the histograms
  htool.plot1D( dx, location + "/x_res", "x resolution / mm", -0.4, 0.4, 101 );
  htool.plot1D( dy, location + "/y_res", "y resolution / mm", -0.4, 0.4, 101 );
  htool.plot1D( dtx, location + "/tx_res", "tx resolution", -0.0025, 0.0025, 101 );
  htool.plot1D( dty, location + "/ty_res", "ty resolution", -0.0025, 0.0025, 101 );

  htool.plot1D( dx / sqrt( cov( 0, 0 ) + trueCov( 0, 0 ) ), location + "/xpull", "x pull", -5., 5., 101 );
  htool.plot1D( dy / sqrt( cov( 1, 1 ) + trueCov( 1, 1 ) ), location + "/ypull", "y pull", -5., 5., 101 );
  htool.plot1D( dtx / sqrt( cov( 2, 2 ) + trueCov( 2, 2 ) ), location + "/txpull", "tx pull", -5., 5., 101 );
  htool.plot1D( dty / sqrt( cov( 3, 3 ) + trueCov( 3, 3 ) ), location + "/typull", "ty pull", -5., 5., 101 );

  if ( std::abs( cov( 4, 4 ) ) > 1e-20 ) { // test that there was a momentum measurement
    const double qop      = vec( 4 );
    const double qoptrue  = trueVec( 4 );
    const double invp     = std::abs( qop );
    const double invptrue = std::abs( qoptrue );
    const double qoperr   = std::sqrt( cov( 4, 4 ) + trueCov( 4, 4 ) );
    // make two pulls, to be sensitive to both a curvature and a momentum bias
    htool.plot1D( ( qop - qoptrue ) / invptrue, location + "/qop_res", "qop ", -.02, .02, 101 );
    htool.plot1D( ( qop - qoptrue ) / qoperr, location + "/qoppull", "qop pull", -5., 5., 101 );
    htool.plot1D( ( invp - invptrue ) / qoperr, location + "/ppull", "p pull", -5., 5., 101 );
    htool.plot1D( invp / invptrue - 1, location + "/dpoverp", "dp/p", -0.05, 0.05, 101 );
    if ( invp > 0 )
      htool.plot1D( std::sqrt( cov( 4, 4 ) ) / invp, location + "/expecteddpoverp", "expected dp/p", 0., 0.01, 100 );
  }
}

//=============================================================================
//
//=============================================================================
void TrackResChecker::plotsByMeasType( IHistoTool const& htool, LHCb::Track const& track,
                                       LHCb::MCParticle const& mcPart, IGeometryInfo const& geometry ) const {
  for ( const auto& measure : measurements( track ) ) {

    LHCb::State trueStateAtMeas;
    StatusCode  sc = idealStateCreator()->createState( &mcPart, measure.z(), trueStateAtMeas, geometry );
    if ( sc.isSuccess() ) {

      const std::string dir = Gaudi::Utils::toString( measure.type() );
      LHCb::State       stateAtMeas;
      StatusCode        sc = extrapolator()->propagate( track, measure.z(), stateAtMeas, geometry );
      if ( sc.isSuccess() ) {
        // make pull plots as before
        pullplots( htool, trueStateAtMeas, stateAtMeas, dir );
      }

      // Monitor unbiased measurement resolutions
      ITrackProjector* proj = m_projectorSelector->projector( measure );
      if ( proj != 0 ) {

        auto [sc, projectResult] = proj->project( trueStateAtMeas, measure );

        if ( sc ) {

          struct Result {
            double res, errMeasure, chi2;
            int    ndf;
          };

          const auto [res, errMeasure, chi2, ndf] =
              std::visit( Gaudi::overload(
                              [&]( ITrackProjector::Project1DResult& projResult ) -> Result {
                                auto res        = projResult.residual[0];
                                auto errMeasure = projResult.errMeasure[0];
                                return {res, errMeasure, ( errMeasure > 0 ? std::pow( res / errMeasure, 2 ) : 0. ), 1};
                              },
                              [&]( ITrackProjector::Project2DResult& projResult ) -> Result {
                                auto residual   = projResult.residual;
                                auto errMeasure = projResult.errMeasure;
                                return {sqrt( pow( residual[0], 2 ) + pow( residual[1], 2 ) ),
                                        sqrt( pow( errMeasure[0], 2 ) + pow( errMeasure[1], 2 ) ),
                                        ( errMeasure > 0 ? std::pow( residual[0] / errMeasure[0], 2 ) +
                                                               std::pow( residual[1] / errMeasure[1], 2 )
                                                         : 0. ),
                                        2};
                              } ),
                          projectResult );

          htool.plot1D( res, dir + "/meas_res", " Measurement resolution", -0.5, 0.5, 100 );
          htool.plot1D( res / errMeasure, dir + "/meas_pull", " Measurement pull", -5., 5., 100 );
          htool.plot1D( chi2, dir + "/meas_chi2", " Measurement chi2", 0., 10., 200 );
          htool.plot2D( chi2, ndf, dir + "/meas_chi2ndf", "Measurement chi2 vs ndf", 0., 10., 0., 3., 200, 2 );
        }
      } else {
        Warning( "could not get projector for measurement", StatusCode::SUCCESS, 0 ).ignore();
      }
    }
  } // iterate measurements
}

//=============================================================================
//
//=============================================================================
StatusCode TrackResChecker::finalize() {
  info() << "     ************************************    " << endmsg;
  for ( const auto& ihtool : m_histoTools ) {
    const IHistoTool*     htool  = ihtool.second;
    const GaudiHistoTool* ghtool = dynamic_cast<const GaudiHistoTool*>( htool );
    for ( const auto& name :
          {"vertex/xpull", "vertex/ypull", "vertex/txpull", "vertex/typull", "vertex/ppull", "probChi2"} ) {
      const auto pull = htool->histo( HistoID( name ) );
      if ( pull )
        info() << ghtool->histoDir() << "/" << std::setiosflags( std::ios_base::left ) << std::setw( 10 )
               << pull->title() << " "
               << format( ":  mean =  %5.3f +/- %5.3f, RMS = %5.3f +/- %5.3f", pull->mean(),
                          Gaudi::Utils::HistoStats::meanErr( pull ), pull->rms(),
                          Gaudi::Utils::HistoStats::rmsErr( pull ) )
               << endmsg;
    }
    for ( const auto& name : {"vertex/x_res", "vertex/y_res"} ) {
      const auto res = htool->histo( HistoID( name ) );
      if ( res )
        info() << ghtool->histoDir() << "/" << res->title()
               << format( ":  RMS =  %5.3f +/- %5.3f micron", res->rms() * 1000,
                          Gaudi::Utils::HistoStats::rmsErr( res ) * 1000 )
               << endmsg;
    }
    const auto dpop = htool->histo( HistoID( "vertex/dpoverp" ) );
    if ( dpop )
      info() << ghtool->histoDir() << "/" << dpop->title()
             << format( ":  mean =  %6.4f +/- %6.4f, RMS =  %6.4f +/- %6.4f", dpop->mean(),
                        Gaudi::Utils::HistoStats::meanErr( dpop ), dpop->rms(),
                        Gaudi::Utils::HistoStats::rmsErr( dpop ) )
             << endmsg;
  }
  return TrackCheckerBase::finalize();
}
