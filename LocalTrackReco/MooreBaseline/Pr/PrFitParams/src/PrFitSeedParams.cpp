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
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/MCVertex.h"
#include "Event/StateParameters.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "PrFitParameters.h"
#include "PrFitParams/IPrFitTool.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrFitSeedParams
//
// 2006-12-08  : Olivier Callot
// 2013-01-23  : Yasmine Amhis
// Adapt to work with Fiber Tracker and UT
//-----------------------------------------------------------------------------

/** @class PrFitSeedParams PrFitSeedParams.h
 *
 *
 *  @author Olivier Callot
 *  @date   2012-07-03
 *  @modification on 2013-01-23  : Yasmine Amhis
 *  Adapt to work with Fiber Tracker and UT

 */
class PrFitSeedParams : public GaudiTupleAlg {
public:
  /// Standard constructor
  PrFitSeedParams( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution
  StatusCode finalize() override;   ///< Algorithm finalization

private:
  PublicToolHandle<IPrFitTool> m_fitTool{this, "FitTool", "PrFitTool"};
  std::string                  m_tupleName;
  double                       m_zRef;
  double                       m_zSeed;
  double                       m_zUT;
  int                          m_nEvent = 0;
  int                          m_nTrack = 0;

  std::vector<double> m_momentumScaleParams;
  std::vector<double> m_initialArrowParams;
  std::vector<double> m_zMagParams;

  PrFitParameters m_momentumScalePar;
  PrFitParameters m_initialArrowPar;
  PrFitParameters m_zMagPar;

  std::vector<double> m_dRatio;
  PrFitParameters     m_dRatioPar;

  std::vector<double> m_yCorrection;
  PrFitParameters     m_yCorrectionPar;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( PrFitSeedParams )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrFitSeedParams::PrFitSeedParams( const std::string& name, ISvcLocator* pSvcLocator )
    : GaudiTupleAlg( name, pSvcLocator ) {
  declareProperty( "TupleName", m_tupleName = "Track" );
  declareProperty( "ZRef", m_zRef = StateParameters::ZMidT );
  declareProperty( "ZSeed", m_zSeed = StateParameters::ZEndT );
  // THIS NEEDS TO BE UPDATED FOR THE UT
  declareProperty( "ZUT", m_zUT = 2485.0 * Gaudi::Units::mm );

  // these are the "tunables" for PatSeeding
  // for many of these parameters, PatSeeding only takes the first one
  // while we fit for several to disentangle dependencies which would
  // move the values we'd obtain away from their "best" value
  declareProperty( "InitialArrowParams", m_initialArrowParams = {4.21826e-09, -8.93796e-08, 0.372981} );
  declareProperty( "MomentumScaleParams", m_momentumScaleParams = {40.3751, 1163.24, -682850} );
  declareProperty( "zMagnetParams", m_zMagParams = {5372.27, -3111.41, 384.74} );
  declareProperty( "dRatioParams", m_dRatio = {-3.77e-4, 4.7} );
  declareProperty( "yCorrectionParams", m_yCorrection = {1.25e-14} );
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrFitSeedParams::initialize() {
  return GaudiTupleAlg::initialize().andThen( [&] {
    m_initialArrowPar.init( "InitialArrow", m_initialArrowParams );
    m_momentumScalePar.init( "MomentumScale", m_momentumScaleParams );
    m_zMagPar.init( "zMagnet", m_zMagParams );
    m_dRatioPar.init( "dRatio", m_dRatio );
    m_yCorrectionPar.init( "yCorrection", m_yCorrection );
  } );
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode PrFitSeedParams::execute() {

  debug() << "==> Execute" << endmsg;

  m_nEvent += 1;

  LHCb::MCParticles* partCtnr = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );

  auto trackInfo = MCTrackInfo{*get<LHCb::MCProperty>( LHCb::MCPropertyLocation::TrackInfo )};

  // Get the FT hits
  LHCb::MCHits* ftHits = get<LHCb::MCHits>( LHCb::MCHitLocation::FT );

  // Get the UT  hits
  LHCb::MCHits* UTHits = get<LHCb::MCHits>( LHCb::MCHitLocation::UT );

  for ( const LHCb::MCParticle* part : *partCtnr ) { // here loop over the MC particle container
    if ( 0 == trackInfo.fullInfo( part ) ) continue;
    if ( !trackInfo.hasT( part ) ) continue;
    if ( !trackInfo.hasUT( part ) && !trackInfo.hasVelo( part ) ) continue;
    const LHCb::MCVertex* vOrigin = part->originVertex();
    if ( !vOrigin ) continue;
    //== No velo -> request that the mother comes from near the beam line
    if ( !trackInfo.hasVelo( part ) ) {
      const LHCb::MCParticle* mother = vOrigin->mother();
      if ( !mother ) continue;
      if ( !mother->originVertex() ) continue;
      double rOrigin = mother->originVertex()->position().rho();
      if ( fabs( rOrigin ) > 4. ) continue;
    }
    bool                           hasInteractionVertex = false;
    SmartRefVector<LHCb::MCVertex> endV                 = part->endVertices();
    for ( SmartRefVector<LHCb::MCVertex>::const_iterator itV = endV.begin(); endV.end() != itV; itV++ ) {
      if ( ( *itV )->position().z() > StateParameters::ZEndT ) continue;
      hasInteractionVertex = true;
    }
    if ( hasInteractionVertex ) continue;

    double momentum = part->momentum().R();
    if ( 0 > part->particleID().threeCharge() ) { momentum = -1 * momentum; }

    // A container for used hits
    std::vector<Gaudi::XYZPoint> trHits;
    std::vector<Gaudi::XYZPoint> utHits;

    // Get the FT hits

    for ( const auto& fiberHit : *ftHits ) {
      if ( fiberHit->mcParticle() == part ) { trHits.push_back( fiberHit->midPoint() ); }
    }

    // Get the UT hits
    for ( const auto& iHitut : *UTHits ) {
      if ( iHitut->mcParticle() == part ) { utHits.push_back( iHitut->midPoint() ); }
    }

    m_nTrack += 1;
    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << format( "Track MC %4d z0 %7.2f p%8.2f Velo%2d UT%2d T%2d", part->key(), vOrigin->position().z(),
                         momentum / 1000, trackInfo.hasVelo( part ), trackInfo.hasUT( part ), trackInfo.hasT( part ) )
              << endmsg;
      for ( std::vector<Gaudi::XYZPoint>::iterator itP = trHits.begin(); trHits.end() != itP; ++itP ) {
        debug() << format( "  x%8.2f y%8.2f z%8.2f ", ( *itP ).x(), ( *itP ).y(), ( *itP ).z() ) << endmsg;
      }
    }

    Tuple tTrack = nTuple( m_tupleName, m_tupleName );
    tTrack->column( "moment", momentum ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    double pz = part->momentum().Z();
    tTrack->column( "slopeX", part->momentum().X() / pz ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "slopeY", part->momentum().Y() / pz ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    const auto xResult = m_fitTool->fitCubic( trHits, IPrFitTool::XY::X, m_zRef );
    const auto yResult = m_fitTool->fitLine( trHits, IPrFitTool::XY::Y, m_zRef );
    if ( !xResult || !yResult ) {
      err() << "Fit matrix is singular" << endmsg;
      continue;
    }

    const auto& [ax, bx, cx, dx] = *xResult;
    const auto& [ay, by]         = *yResult;

    m_momentumScalePar.setFun( 0, 1. );
    m_momentumScalePar.setFun( 1, by * by );
    m_momentumScalePar.setFun( 2, fabs( cx ) );

    double dp = 1. / ( cx * momentum ) - m_momentumScalePar.sum();
    m_momentumScalePar.addEvent( dp );

    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << format( "p %7.3f, N%4d, ax%8.2f bx%8.2f cx%8.2f dp%10.4f", momentum / 1000., trHits.size(), ax,
                         1.e3 * bx, 1.e6 * cx, dp )
              << endmsg;
    }

    // work out intercept point of line joining T1 and T3 with z = 0
    const double dzT1    = StateParameters::ZBegT - m_zRef;
    const double dzT3    = StateParameters::ZEndT - m_zRef;
    const double xT1     = ( ( dx * dzT1 + cx ) * dzT1 + bx ) * dzT1 + ax;
    const double xT3     = ( ( dx * dzT3 + cx ) * dzT3 + bx ) * dzT3 + ax;
    const double tx13    = ( xT3 - xT1 ) / ( dzT3 - dzT1 );
    const double xAtZero = ax - tx13 * m_zRef;
    // initial arrow depends in y and cx in T, so we determine it as
    // InitialArrow(by) = InitialArrow[0] + InitialArrow[1] * by * by
    // 		+ InitialArrow[2] * cx * cx
    m_initialArrowPar.setFun( 0, 1. );
    m_initialArrowPar.setFun( 1, by * by );
    m_initialArrowPar.setFun( 2, cx * cx );
    double dSag = 0.;
    if ( 1. < fabs( xAtZero ) ) {
      dSag = -cx / xAtZero - m_initialArrowPar.sum();
      m_initialArrowPar.addEvent( dSag );
    }

    double axt, bxt, cxt, ayt, byt, zMagnet, zEst;
    if ( 2 < utHits.size() ) {
      // m_fitTool->fitLine( ttHits, 0, m_zUT, axt, bxt );
      const auto xResult = m_fitTool->fitParabola( utHits, IPrFitTool::XY::X, m_zUT );
      const auto yResult = m_fitTool->fitLine( utHits, IPrFitTool::XY::Y, m_zUT );
      if ( !xResult || !yResult ) {
        err() << "UT fit matrix is singular" << endmsg;
        continue;
      }
      std::tie( axt, bxt, cxt ) = *xResult;
      std::tie( ayt, byt )      = *yResult;

      const double dz = m_zSeed - m_zRef;
      const double x0 = ax + dz * ( bx + dz * ( cx + dz * dx ) );
      const double tx = bx + dz * ( 2 * cx + dz * 3 * dx );

      zMagnet = ( axt - bxt * m_zUT - ( x0 - tx * m_zSeed ) ) / ( tx - bxt );

      m_zMagPar.setFun( 0, 1. );
      m_zMagPar.setFun( 1, by * by );
      m_zMagPar.setFun( 2, bx * bx );
      zEst = m_zMagPar.sum();
      if ( fabs( zMagnet - zEst ) < 100. ) m_zMagPar.addEvent( zMagnet - zEst );

    } else {
      zMagnet = 0.;
      zEst    = 0.;
      axt     = 0.;
      bxt     = 0.;
      cxt     = 0.;
      ayt     = 0.;
      byt     = 0.;
    }

    m_dRatioPar.setFun( 0, 1.0 );
    m_dRatioPar.setFun( 1, fabs( cx ) );
    if ( fabs( cx ) > 1e-10 && fabs( dx ) > 1e-15 ) m_dRatioPar.addEvent( dx / cx - m_dRatioPar.sum() );

    // check if we can improve pointing constraint in y
    // use only tracks which do not point to the origin too well, for those
    // with fabs(y(y=0)) < 40 cm, the correction is small anyway, so we avoid
    // to pollute our statistics with them (if the correction is slightly
    // off for the well pointing tracks, they will still point rather well if
    // we under- or overcorrect a little)
    const double y0 = ay - by * m_zRef;
    m_yCorrectionPar.setFun( 0, 1.0 );
    if ( 400. < fabs( y0 ) ) m_yCorrectionPar.addEvent( fabs( ( by * by * cx * cx ) / y0 ) - m_yCorrectionPar.sum() );

    tTrack->column( "ax", ax ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "bx", bx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "cx", cx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "dx", dx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "ay", ay ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "by", by ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    tTrack->column( "axt", axt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "bxt", bxt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "cxt", cxt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    tTrack->column( "ayt", ayt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "byt", byt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    tTrack->column( "x0", xAtZero ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "dp", dp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "dSag", dSag ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    tTrack->column( "zMagnet", zMagnet ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "zEst", zEst ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    tTrack->write().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

  } // end of the big loop

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrFitSeedParams::finalize() {

  debug() << "==> Finalize" << endmsg;
  MsgStream& msg = info();
  msg << "============================================" << endmsg;
  msg << "  Processed " << m_nEvent << " events and " << m_nTrack << " tracks. " << endmsg;
  msg << "============================================" << endmsg;

  m_initialArrowPar.updateParameters( msg );
  m_momentumScalePar.updateParameters( msg );
  m_zMagPar.updateParameters( msg );
  m_dRatioPar.updateParameters( msg );
  m_yCorrectionPar.updateParameters( msg );

  std::cout << std::endl;
  m_initialArrowPar.printPythonParams( name() );
  m_momentumScalePar.printPythonParams( name() );
  m_zMagPar.printPythonParams( name() );
  m_dRatioPar.printPythonParams( name() );
  m_yCorrectionPar.printPythonParams( name() );
  std::cout << std::endl;

  std::string toolName = "ToolSvc.PatSeedingTool";

  m_initialArrowPar.printPythonParams( toolName );
  m_momentumScalePar.printPythonParams( toolName );
  m_zMagPar.printPythonParams( toolName );
  m_dRatioPar.printPythonParams( toolName );
  m_yCorrectionPar.printPythonParams( toolName );
  std::cout << std::endl;

  return GaudiTupleAlg::finalize(); // must be called after all other actions
}

//=============================================================================