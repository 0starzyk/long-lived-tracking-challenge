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
#include "Event/State.h"
#include "Event/Track.h"
#include "Event/TrackParameters.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "Kernel/ILHCbMagnetSvc.h"
#include "Linker/LinkedFrom.h"
#include "Linker/LinkedTo.h"
#include "PrFitParameters.h"
#include "PrFitParams/IPrFitTool.h"
#include <vector>

//-----------------------------------------------------------------------------
// Implementation file for class : PrFitMatchParams
//
// 2017-02-22 : Sevda Esen
//-----------------------------------------------------------------------------

/** @class PrFitMatchParams PrFitMatchParams.h
 *  Parameterize the PrMatch tracks
 *
 *  Algorithm to calculate the different parameters in PatMatch, similar to FwdFitParams it tries to fit polynomial
 * dependencies. Uses purely MC information from MCHits ie this will only run when they are available (XDIGI, XDST) Only
 * Velo and T-station information is used. Some UT information is calculated, but not used at the moment.
 *
 *  - NTupleName: Name of the output nTuple
 *  - ZUT1: z-position of reference point in UT
 *  - ZRef: z-position of reference plane in T-stations
 *  - ZVelo: z-position of reference plane in Velo
 *  - zMagnetParams: Initial parameters for calculation of z of 'kink-position' of magnet
 *  - momentumParams: Initial paramters to calculate the momentum.
 *  - bendYParams: Initial paramters to calculate the bending in y.
 *  - MaxZVertex: Maximum z-position of PV to take tracks into account.
 *  - MinMomentum: Minimum momentum to take tracks into account.
 *  - zMatchY: z-position of matching in Y.
 *  - RequireUTHits: Require presence of UT hits for the calculation
 *  - MagnetScaleFactor: The scale factor for the magnet. Can be obtained with the field service, but hardcoded here
 * such that it can be run without. Note that only as many paramters are used for the caluclation as are given as input,
 * although internally more might be defined. The momentum paramters are not directly set in 'PatMatchTool', but in
 * 'FastMomentumEstimate', and only the cubic solution is calculated at the moment. Nothing else from the Brunel
 * sequence is needed when running this code, a possible way to run it would be:
 *
 *  @code
 *  def doIt():
 *      GaudiSequencer("BrunelSequencer").Members =  [ PrFitMatchParams("PrFitMatchParams", ZRef = 8520, zMagnetParams =
 * [ ... ], momentumParams = [ ... ], bendYParams = [ ... ]) ]
 *
 *  appendPostConfigAction( doIt )
 *  @endcode
 *
 *
 *  @author Michel De Cian
 *  @date   2014-12-21
 *  @copied by Sevda Esen for the upgrade
 */
class PrFitMatchParams : public GaudiTupleAlg {
public:
  /// Standard constructor
  using GaudiTupleAlg::GaudiTupleAlg;

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution
  StatusCode finalize() override;   ///< Algorithm finalization

private:
  ServiceHandle<ILHCbMagnetSvc> m_magFieldSvc{this, "MagneticField", "MagneticFieldSvc"};
  PublicToolHandle<IPrFitTool>  m_fitTool{this, "FitTool", "PrFitTool"};

  Gaudi::Property<std::string> m_tupleName{this, "NTupleName", "Track"};
  Gaudi::Property<std::string> m_inputLocationSeed{this, "SeedInput", LHCb::TrackLocation::Seed};

  Gaudi::Property<double> m_zUT1{this, "ZUT1", 2469.0 * Gaudi::Units::mm};
  Gaudi::Property<double> m_zRef{this, "ZRef", 9410.0 * Gaudi::Units::mm};
  Gaudi::Property<double> m_zVelo{this, "ZVelo", 0.0 * Gaudi::Units::mm};
  Gaudi::Property<double> m_maxZVertex{this, "MaxZVertex", 500 * Gaudi::Units::mm};
  Gaudi::Property<double> m_minMomentum{this, "MinMomentum", 2 * Gaudi::Units::GeV};
  Gaudi::Property<double> m_maxMomentum{this, "MaxMomentum", 1000 * Gaudi::Units::GeV};
  Gaudi::Property<double> m_zMatchY{this, "zMatchY", 8420 * Gaudi::Units::mm};
  Gaudi::Property<double> m_magnetScaleFactor{this, "MagnetScaleFactor", -1};

  Gaudi::Property<bool> m_requireUTHits{this, "RequireUTHits", false};
  Gaudi::Property<bool> m_resolution{this, "Resolution", true};

  Gaudi::Property<std::vector<double>> m_zMagParams{this, "zMagnetParams", {0., 0., 0., 0., 0.}};
  Gaudi::Property<std::vector<double>> m_momParams{this, "momentumParams", {0., 0., 0., 0., 0., 0.}};
  Gaudi::Property<std::vector<double>> m_bendYParams{this, "bendParamYParams", {0., 0.}};
  PrFitParameters                      m_zMagPar;
  PrFitParameters                      m_momPar;
  PrFitParameters                      m_bendParamY;

  float momentum( const LHCb::State* vState, const LHCb::State* tState );
  void  resolution();
  int   m_nEvent{0};
  int   m_nTrack{0};
};

DECLARE_COMPONENT( PrFitMatchParams )

//=============================================================================
// Initialisation. Check parameters
//=============================================================================
StatusCode PrFitMatchParams::initialize() {
  StatusCode sc = GaudiTupleAlg::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;             // error printed already by GaudiAlgorithm

  sc = m_magFieldSvc.retrieve();
  if ( sc.isFailure() ) return sc;

  m_magnetScaleFactor = m_magFieldSvc->signedRelativeCurrent(); // BUG/FIXME: this is the current for some a-priori
                                                                // arbitrary time, so basically 'undefined'
  debug() << "==> Initialize" << endmsg;

  if ( m_zMagParams.empty() )
    info() << "no starting values for magnet parameters provided. Will not calculate anything" << endmsg;
  if ( m_momParams.empty() )
    info() << "no starting values for momentum  parameters provided. Will not calculate anything" << endmsg;
  if ( m_bendYParams.empty() )
    info() << "no starting values for y bending  parameters provided. Will not calculate anything" << endmsg;

  m_zMagPar.init( "zMagnet", m_zMagParams );
  m_momPar.init( "momentum", m_momParams );
  m_bendParamY.init( "bendParamY", m_bendYParams );

  return sc;
}
//=============================================================================
// Main execution
//=============================================================================
StatusCode PrFitMatchParams::execute() {

  debug() << "==> Execute" << endmsg;

  m_nEvent++;

  // -- As this algorithm is the only one running, even the usual counter is not included.
  if ( m_nEvent % 100 == 0 ) info() << "Event " << m_nEvent << endmsg;

  debug() << "Processing event: " << m_nEvent << endmsg;

  if ( m_resolution ) {
    resolution();
    return StatusCode::SUCCESS;
  }

  LHCb::MCParticles* partCtnr = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );

  // Get the FT and UT hits
  LHCb::MCHits* ftHits = get<LHCb::MCHits>( LHCb::MCHitLocation::FT );
  LHCb::MCHits* utHits = get<LHCb::MCHits>( LHCb::MCHitLocation::UT );

  // Get the Velo hits
  LHCb::MCHits* veloHits = get<LHCb::MCHits>( LHCb::MCHitLocation::VP );

  const auto* flag = get<LHCb::MCProperty>( LHCb::MCPropertyLocation::TrackInfo );

  SmartRefVector<LHCb::MCVertex> vDecay;

  // A container for used hits
  std::vector<Gaudi::XYZPoint> trHits;
  std::vector<Gaudi::XYZPoint> UTHits;
  std::vector<Gaudi::XYZPoint> vHits;

  MCTrackInfo trackInfo{*flag};

  for ( auto pItr = partCtnr->begin(); partCtnr->end() != pItr; pItr++ ) {
    const LHCb::MCParticle* part = *pItr;

    // only long track reconstructible
    if ( !trackInfo.hasVeloAndT( part ) ) continue;

    // Select particles to work with:
    // == Origin vertex position
    const LHCb::MCVertex* vOrigin = part->originVertex();
    if ( nullptr == vOrigin ) continue;
    if ( m_maxZVertex < vOrigin->position().z() ) continue;
    // == Decay vertices

    // == Momentum cut
    double momentum = part->momentum().R();
    if ( m_minMomentum > momentum ) continue;
    if ( m_maxMomentum < momentum ) continue;
    if ( 11 == abs( part->particleID().pid() ) ) continue; // no electrons

    auto endV = part->endVertices();
    bool hasInteractionVertex =
        std::any_of( endV.begin(), endV.end(), []( const auto& itV ) { return itV->position().z() < 9500.; } );
    if ( hasInteractionVertex ) continue;

    debug() << "--- Found particle key " << part->key() << endmsg;

    UTHits.clear();
    trHits.clear();
    vHits.clear();

    // Get the Velo hits
    for ( auto vHitIt = veloHits->begin(); veloHits->end() != vHitIt; vHitIt++ ) {
      if ( ( *vHitIt )->mcParticle() == part ) { vHits.push_back( ( *vHitIt )->midPoint() ); }
    }

    // Get the FT hits
    for ( auto iHitFt = ftHits->begin(); ftHits->end() != iHitFt; iHitFt++ ) {
      if ( ( *iHitFt )->mcParticle() == part ) { trHits.push_back( ( *iHitFt )->midPoint() ); }
    }

    // Get the UT hits
    for ( auto iHitut = utHits->begin(); utHits->end() != iHitut; iHitut++ ) {
      if ( ( *iHitut )->mcParticle() == part ) { UTHits.push_back( ( *iHitut )->midPoint() ); }
    }

    if ( 12 > trHits.size() || 6 > vHits.size() ) continue;
    if ( m_requireUTHits ) {
      if ( 3 > UTHits.size() ) continue;
    }

    debug() << " particle momentum = " << part->momentum().R() / Gaudi::Units::GeV << " GeV" << endmsg;

    //== Fill ntuple
    double pz   = part->momentum().Z();
    double plXT = part->momentum().X() / pz;
    double plYT = part->momentum().Y() / pz;

    Tuple tTrack = nTuple( m_tupleName );

    tTrack->column( "pz", pz ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "plXT", plXT ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "plYT", plYT ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    m_nTrack++;

    //== Fit the UT area
    // -- This is not needed at the moment for the matching, but it was kept it
    // -- as it might be used at some point in the future
    debug() << "  UT: ";
    const auto utFitXLine  = m_fitTool->fitLine( UTHits, IPrFitTool::XY::X, m_zUT1 );
    const auto utFitY      = m_fitTool->fitLine( UTHits, IPrFitTool::XY::Y, m_zUT1 );
    const auto utFitXParab = m_fitTool->fitParabola( UTHits, IPrFitTool::XY::X, m_zUT1 );
    // m_fitTool->fitLine( UTHits, 1, m_zUT1, ayt, byt, cyt );
    if ( !utFitXLine || !utFitY || !utFitXParab ) {
      err() << "UT fit matrix is singular" << endmsg;
      continue;
    }

    const auto& [axt, bxt]         = *utFitXLine;
    const auto& [ayt, byt]         = *utFitY;
    const auto& [axt2, bxt2, cxt2] = *utFitXParab;

    debug() << format( " x %7.1f tx %7.4f   y %7.1f ty %7.4f ", axt, bxt, ayt, byt ) << endmsg;
    ;
    tTrack->column( "axt", axt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "bxt", bxt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "ayt", ayt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "byt", byt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "axt2", axt2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "bxt2", bxt2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "cxt2", cxt2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    std::vector<Gaudi::XYZPoint>::const_iterator itP;
    if ( msgLevel( MSG::DEBUG ) ) {
      for ( itP = UTHits.begin(); UTHits.end() > itP; itP++ ) {
        const double dz = ( *itP ).z() - m_zUT1;

        debug() << format( "    : %7.1f %7.1f %7.1f  dx %7.3f  dy %7.3f", ( *itP ).x(), ( *itP ).y(), ( *itP ).z(),
                           ( *itP ).x() - ( axt + bxt * dz ), ( *itP ).y() - ( ayt + byt * dz ) )
                << endmsg;
      }
    }

    // -- Fit the T-stations
    // -- x(z) = ax + bx*z + cx*z*z + dx*z*z*z
    // -- y(z) = ay + by*z
    const auto tFitX = m_fitTool->fitCubic( trHits, IPrFitTool::XY::X, m_zRef );
    const auto tFitY = m_fitTool->fitLine( trHits, IPrFitTool::XY::Y, m_zRef );
    if ( !tFitX || !tFitY ) {
      err() << "T-stations fit matrix is singular" << endmsg;
      continue;
    }

    const auto& [ax, bx, cx, dx] = *tFitX;
    const auto& [ay, by]         = *tFitY;

    tTrack->column( "ax", ax ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "bx", bx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "cx", cx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "dx", dx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "ay", ay ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "by", by ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << format( "  tr: x%7.1f %7.4f %7.3f %7.3f  y%7.1f %7.4f", ax, bx, 1.e6 * cx, 1.e9 * dx, ay, by )
              << endmsg;

      for ( itP = trHits.begin(); trHits.end() > itP; itP++ ) {
        const double dz = ( *itP ).z() - m_zRef;
        debug() << format( "    : %7.1f %7.1f %7.1f  dx %7.3f  dy %7.3f", ( *itP ).x(), ( *itP ).y(), ( *itP ).z(),
                           ( *itP ).x() - ( ax + bx * dz + cx * dz * dz + dx * dz * dz * dz ),
                           ( *itP ).y() - ( ay + by * dz ) )
                << endmsg;
      }
    }

    // -- Fit the velo area
    // -- x(z) = axv + bxv*z
    // -- y(z) = ayv + byv*z
    const auto veloFitX = m_fitTool->fitLine( vHits, IPrFitTool::XY::X, m_zVelo );
    const auto veloFitY = m_fitTool->fitLine( vHits, IPrFitTool::XY::Y, m_zVelo );
    if ( !veloFitX || !veloFitY ) {
      err() << "Velo fit matrix is singular" << endmsg;
      continue;
    }

    const auto& [axv, bxv] = *veloFitX;
    const auto& [ayv, byv] = *veloFitY;

    tTrack->column( "axv", axv ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "bxv", bxv ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "ayv", ayv ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "byv", byv ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << format( "  velo: x%7.1f %7.4f y%7.1f %7.4f", axv, bxv, ayv, byv ) << endmsg;

      for ( itP = vHits.begin(); vHits.end() > itP; itP++ ) {
        const double dz = ( *itP ).z() - m_zRef;
        debug() << format( "    : %7.1f %7.1f %7.1f  dx %7.3f  dy %7.3f", ( *itP ).x(), ( *itP ).y(), ( *itP ).z(),
                           ( *itP ).x() - ( ax + bx * dz + cx * dz * dz + dx * dz * dz * dz ),
                           ( *itP ).y() - ( ay + by * dz ) )
                << endmsg;
      }
    }

    // -- This is for finding the zMagnet position when using straight lines from the Velo and the T-stations
    // -- Only really makes sense in x, as for y the different y resolutions of Velo and T-stations matter a lot
    double zMagnet  = ( axv - bxv * m_zVelo - ( ax - bx * m_zRef ) ) / ( bx - bxv );
    double zMagnety = ( ayt - byt * m_zUT1 - ( ay - by * m_zRef ) ) / ( by - byt );
    double dSlope   = fabs( bx - bxv );
    double dSlopeY  = fabs( by - byv );

    m_zMagPar.setFun( 0, 1. );
    m_zMagPar.setFun( 1, fabs( dSlope ) );
    m_zMagPar.setFun( 2, dSlope * dSlope );
    m_zMagPar.setFun( 3, fabs( ax ) );
    m_zMagPar.setFun( 4, bxv * bxv );

    double zEst = m_zMagPar.sum();

    m_zMagPar.addEvent( zMagnet - zEst );

    tTrack->column( "zEst", zEst ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "zMagnet", zMagnet ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "zMagnety", zMagnety ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    // -- This is the parameter that defines the bending in y
    double bendParamY = ( ay - ayv ) + ( m_zMatchY - m_zRef ) * by - ( m_zMatchY - m_zVelo ) * byv;

    m_bendParamY.setFun( 0, dSlope * dSlope * byv );
    m_bendParamY.setFun( 1, dSlopeY * dSlopeY * byv );

    double bendParamYEst = m_bendParamY.sum();
    m_bendParamY.addEvent( bendParamY - bendParamYEst );
    // ------------------------------------------------------

    tTrack->column( "dSlope", dSlope ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    tTrack->column( "dSlope2", dSlope * dSlope ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    // -- Need to write something to calculate the momentum Params
    const double charge = part->particleID().threeCharge() / 3;
    const double qOverP = charge / momentum;

    // -- The magnet scale factor is hardcoded, as then one does not need to run the
    // -- field service
    const double proj = sqrt( ( 1. + bxv * bxv + byv * byv ) / ( 1. + bxv * bxv ) );
    const double coef = ( bxv - bx ) / ( proj * m_magnetScaleFactor * Gaudi::Units::GeV * qOverP );

    m_momPar.setFun( 0, 1. );
    m_momPar.setFun( 1, bx * bx );
    m_momPar.setFun( 2, bx * bx * bx * bx );
    m_momPar.setFun( 3, bx * bxv );
    m_momPar.setFun( 4, byv * byv );
    m_momPar.setFun( 5, byv * byv * byv * byv );

    double coefEst = m_momPar.sum();
    m_momPar.addEvent( coef - coefEst );

    tTrack->write().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Determine resolution
//=============================================================================
void PrFitMatchParams::resolution() {

  LHCb::LinksByKey const* velo_links =
      get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::TrackLocation::Velo ) );
  LHCb::LinksByKey const* seed_links =
      get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::TrackLocation::Seed ) );

  // Get the IT and TT hits
  LHCb::MCHits const* ftHits = get<LHCb::MCHits>( LHCb::MCHitLocation::FT );
  // LHCb::MCHits* utHits = get<LHCb::MCHits>( LHCb::MCHitLocation::UT );
  LHCb::MCHits const* vpHits = get<LHCb::MCHits>( LHCb::MCHitLocation::VP );

  // A container for used hits
  std::vector<Gaudi::XYZPoint> FTHits;
  std::vector<Gaudi::XYZPoint> VPHits;

  LHCb::MCParticles* mcParts = getIfExists<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );
  if ( msgLevel( MSG::ERROR ) && !mcParts )
    error() << "Could not find MCParticles at: " << LHCb::MCParticleLocation::Default << endmsg;

  const auto* flags = get<LHCb::MCProperty>( LHCb::MCPropertyLocation::TrackInfo );
  MCTrackInfo trackInfo{*flags};

  for ( const LHCb::MCParticle* mcPart : *mcParts ) {

    if ( mcPart == nullptr ) continue;
    if ( std::abs( mcPart->particleID().pid() ) == 11 ) continue;

    if ( !trackInfo.hasVeloAndT( mcPart ) ) continue;

    auto vtracks = LinkedFrom<LHCb::Track>{velo_links}.range( mcPart );
    auto ttracks = LinkedFrom<LHCb::Track>{seed_links}.range( mcPart );

    if ( vtracks.empty() ) continue;
    if ( ttracks.empty() ) continue;

    const LHCb::MCVertex* vOrigin = mcPart->originVertex();
    if ( !vOrigin ) continue;
    if ( 100. < vOrigin->position().R() ) {
      debug() << "Too far from beampipe" << endmsg;
      continue; // particles from close the beam line
    }

    const auto& endV = mcPart->endVertices();
    bool        hasInteractionVertex =
        std::any_of( endV.begin(), endV.end(), []( const auto& itV ) { return itV->position().z() < 9500.; } );
    if ( hasInteractionVertex ) {
      debug() << "Interaction vertex found. skipping" << endmsg;
      continue;
    }

    FTHits.clear();
    VPHits.clear();

    // Get the Velo hits
    for ( auto vHitIt = vpHits->begin(); vpHits->end() != vHitIt; vHitIt++ ) {
      if ( ( *vHitIt )->mcParticle() == mcPart ) { VPHits.push_back( ( *vHitIt )->midPoint() ); }
    }

    // Get the FT hits
    for ( auto iHitFt = ftHits->begin(); ftHits->end() != iHitFt; iHitFt++ ) {
      if ( ( *iHitFt )->mcParticle() == mcPart ) { FTHits.push_back( ( *iHitFt )->midPoint() ); }
    }

    if ( 12 > FTHits.size() || 6 > VPHits.size() ) continue;

    debug() << "got my hits  " << FTHits.size() + VPHits.size() << endmsg;

    //== Fit the TT area
    debug() << "  VP: ";
    const auto vpFitX = m_fitTool->fitLine( VPHits, IPrFitTool::XY::X, m_zVelo );
    const auto vpFitY = m_fitTool->fitLine( VPHits, IPrFitTool::XY::Y, m_zVelo );
    if ( !vpFitX || !vpFitY ) {
      err() << "VP fit matrix is singular" << endmsg;
      continue;
    }
    const auto& [axv, bxv] = *vpFitX;
    const auto& [ayv, byv] = *vpFitY;

    const auto ftFitX = m_fitTool->fitCubic( FTHits, IPrFitTool::XY::X, m_zRef );
    const auto ftFitY = m_fitTool->fitLine( FTHits, IPrFitTool::XY::Y, m_zRef );
    if ( !ftFitX || !ftFitY ) {
      err() << "FT fit matrix is singular" << endmsg;
      continue;
    }
    const double ax = std::get<0>( *ftFitX ), bx = std::get<1>( *ftFitX ),
                 // cx = std::get<2>(*ftFitX),
                 // dx = std::get<3>(*ftFitX),
                 // ay = std::get<0>(*ftFitY),
        by = std::get<1>( *ftFitY );

    const auto zMagnet = ( axv - bxv * m_zVelo - ( ax - bx * m_zRef ) ) / ( bx - bxv );
    const auto xMagnet = axv + ( zMagnet - m_zVelo ) * bxv;
    const auto yMagnet = ayv + ( zMagnet - m_zVelo ) * byv;

    LHCb::State const* tstate = &ttracks.front().closestState( 10000. );
    LHCb::State const* vstate = &vtracks.front().closestState( m_zVelo );
    // now other things
    const auto dSlopeExp   = std::abs( vstate->tx() - tstate->tx() );
    const auto dSlope2Exp  = dSlopeExp * dSlopeExp;
    const auto dSlopeYExp  = std::abs( vstate->ty() - tstate->ty() );
    const auto dSlopeY2Exp = dSlopeYExp * dSlopeYExp;

    const auto zMagnetExp = m_zMagParams[0] + m_zMagParams[1] * dSlopeExp + m_zMagParams[2] * dSlope2Exp +
                            m_zMagParams[3] * std::abs( tstate->x() ) + m_zMagParams[4] * vstate->tx() * vstate->tx();

    debug() << m_zMagParams[0] << "  " << m_zMagParams[1] << "  " << m_zMagParams[2] << "  " << m_zMagParams[3] << "  "
            << m_zMagParams[4] << endmsg;

    const float xVExp = vstate->x() + ( zMagnetExp - vstate->z() ) * vstate->tx();
    // -- This is the function that calculates the 'bending' in y-direction
    // -- The parametrisation can be derived with the MatchFitParams package
    const float yVExp = vstate->y() + ( m_zMatchY - vstate->z() ) * vstate->ty() +
                        vstate->ty() * ( m_bendYParams[0] * dSlope2Exp + m_bendYParams[1] * dSlopeY2Exp );

    const float xSExp = tstate->x() + ( zMagnetExp - tstate->z() ) * tstate->tx();
    const float ySExp = tstate->y() + ( m_zMatchY - tstate->z() ) * tstate->ty();

    const float distXExp = xSExp - xVExp;
    const float distYExp = ySExp - yVExp;

    const float dslXExp = vstate->tx() - tstate->tx();
    const float dslYExp = vstate->ty() - tstate->ty();
    const float teta2   = vstate->ty() * vstate->ty() + tstate->tx() * tstate->tx();

    // -- Check difference in momentum
    const auto p    = mcPart->momentum().R();
    const auto pExp = momentum( vstate, tstate );

    Tuple resoTuple = nTuple( "resoTuple", "resoTuple" );
    resoTuple->column( "zMagnet", zMagnet ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "yMagnet", yMagnet ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "xMagnet", xMagnet ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "zMagnetExp", zMagnetExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "yVExp", yVExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "ySExp", ySExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "distXExp", distXExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "distYExp", distYExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "xVExp", xVExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "xSExp", xSExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "dSlopeExp", dSlopeExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "pExp", pExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "pTrue", p ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "byvExp", vstate->ty() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "byv", byv ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "bxvExp", vstate->tx() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "bxv", bxv ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "byExp", tstate->ty() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "by", by ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "bxExp", tstate->tx() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "bx", bx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->column( "dslXExp", dslXExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "dslYExp", dslYExp ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    resoTuple->column( "teta2", teta2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    resoTuple->write().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  }
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrFitMatchParams::finalize() {

  debug() << "==> Finalize" << endmsg;

  MsgStream& msg = info() << "============================================" << endmsg;
  msg << "  Processed " << m_nEvent << " events and " << m_nTrack << " tracks. " << endmsg;
  msg << "============================================" << endmsg;
  m_zMagPar.updateParameters( msg );
  m_momPar.updateParameters( msg );
  m_bendParamY.updateParameters( msg );

  std::cout << std::endl;
  m_zMagPar.printPythonParams( name() );
  m_momPar.printPythonParams( name() );
  m_bendParamY.printPythonParams( name() );
  std::cout << std::endl;

  m_zMagPar.printParams( "PrMatch" );
  m_momPar.printParams( "PrMatch" );
  m_bendParamY.printParams( "PrMatch" );
  std::cout << std::endl;

  return GaudiTupleAlg::finalize(); // must be called after all other actions
}

//=============================================================================
//  Finalize
//=============================================================================
float PrFitMatchParams::momentum( const LHCb::State* vState, const LHCb::State* tState ) {
  const float txT     = tState->tx();
  const float txV     = vState->tx();
  const float tyV     = vState->ty();
  float       vars[5] = {txT * txT, txT * txT * txT * txT, txT * txV, tyV * tyV, tyV * tyV * tyV * tyV};

  const float coef = m_momParams[0] * vars[0] + m_momParams[1] * vars[1] + m_momParams[2] * vars[2] +
                     m_momParams[3] * vars[3] + m_momParams[4] * vars[4];

  debug() << m_momParams[0] << "  " << m_momParams[1] << "  " << m_momParams[2] << "  " << m_momParams[3] << "  "
          << m_momParams[4] << endmsg;

  debug() << vars[0] << "  " << vars[1] << "  " << vars[2] << "  " << vars[3] << "  " << vars[4] << endmsg;

  float proj        = sqrt( ( 1. + txV * txV + tyV * tyV ) / ( 1. + txV * txV ) );
  float scaleFactor = m_magFieldSvc->signedRelativeCurrent();

  float P = ( coef * Gaudi::Units::GeV * proj * scaleFactor ) / ( txV - txT );

  return P;
}
