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
#include "Event/MCVertex.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "PrFitParameters.h"
#include "PrFitParams/IPrFitTool.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrFitKsParams
//
// 2002-11-02 : Olivier Callot
// 2013-01-23  : Yasmine Amhis
// Adapt to work with Fiber Tracker and UT
//-----------------------------------------------------------------------------

/** @class PrFitKsParams PrFitKsParams.h
 *  Parameterize the KShort tracks
 *
 *  @author Olivier Callot
 *  @date   2012-07-03
 *  @modification on 2013-01-23  : Yasmine Amhis
 *  Adapt to work with Fiber Tracker and UT
 */
class PrFitKsParams : public GaudiTupleAlg {
public:
  /// Standard constructor
  using GaudiTupleAlg::GaudiTupleAlg;

  StatusCode execute() override;  ///< Algorithm execution
  StatusCode finalize() override; ///< Algorithm finalization

private:
  PublicToolHandle<IPrFitTool> m_fitTool{this, "FitTool", "PrFitTool"};
  Gaudi::Property<std::string> m_tupleName{this, "NTupleName", "Track"};
  Gaudi::Property<double>      m_zTT1{this, "ZTT1", 2469.0 * Gaudi::Units::mm}; // THIS NEEDS TO BE CHANGED FOR UT
  Gaudi::Property<double>      m_zRef{this, "ZRef", 9410.0 * Gaudi::Units::mm};

  Gaudi::Property<PrFitParameters> m_zMagPar{this, "zMagnetParams", {"zMagnet", {}}};
  Gaudi::Property<PrFitParameters> m_momPar{this, "momentumParams", {"momentum", {}}};

  int m_nEvent = 0;
  int m_nTrack = 0;
};

DECLARE_COMPONENT( PrFitKsParams )

//=============================================================================
// Main execution
//=============================================================================
StatusCode PrFitKsParams::execute() {
  debug() << "==> Execute" << endmsg;

  m_nEvent++;

  LHCb::MCParticles* partCtnr = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );

  // Get the UT  hits
  LHCb::MCHits* utHits = get<LHCb::MCHits>( LHCb::MCHitLocation::UT );
  // Get the FT hits
  LHCb::MCHits* ftHits = get<LHCb::MCHits>( LHCb::MCHitLocation::FT );

  LHCb::MCParticles::const_iterator pItr;
  const LHCb::MCParticle*           part;
  SmartRefVector<LHCb::MCVertex>    vDecay;

  const LHCb::MCParticle* kShort = 0;
  const LHCb::MCVertex*   kDecay = 0;

  // A container for used hits
  std::vector<Gaudi::XYZPoint> trHits;
  std::vector<Gaudi::XYZPoint> UTHits;

  for ( pItr = partCtnr->begin(); partCtnr->end() != pItr; pItr++ ) {
    part = *pItr;

    // Find the Pi from the kShort
    if ( 211 == abs( part->particleID().pid() ) ) {
      kDecay = part->originVertex();
      if ( 0 == kDecay ) continue;
      kShort = kDecay->mother();
      if ( 0 == kShort ) continue;
      if ( 310 != kShort->particleID().pid() ) continue;
      const LHCb::MCVertex* origin = kShort->originVertex();
      if ( 0 == origin ) continue;
      if ( 8. < origin->position().R() ) continue; // Kshorts from close the beam line

      bool                           hasInteractionVertex = false;
      SmartRefVector<LHCb::MCVertex> endV                 = part->endVertices();
      for ( SmartRefVector<LHCb::MCVertex>::const_iterator itV = endV.begin(); endV.end() != itV; itV++ ) {
        if ( ( *itV )->position().z() < 9500. ) hasInteractionVertex = true;
      }
      if ( hasInteractionVertex ) continue;

      debug() << "--- Found pi key " << part->key() << endmsg;

      UTHits.clear();
      trHits.clear();

      // Get the UT hits
      for ( LHCb::MCHits::const_iterator iHitut = utHits->begin(); utHits->end() != iHitut; iHitut++ ) {
        if ( ( *iHitut )->mcParticle() == part ) { UTHits.push_back( ( *iHitut )->midPoint() ); }
      }

      // Get the FT hits
      for ( LHCb::MCHits::const_iterator fHitIt = ftHits->begin(); ftHits->end() != fHitIt; fHitIt++ ) {
        if ( ( *fHitIt )->mcParticle() == part ) { trHits.push_back( ( *fHitIt )->midPoint() ); }
      }
      if ( 3 > UTHits.size() || 12 > trHits.size() ) continue;

      debug() << "=== Found a good K0S Decay : " << kShort->key() << " decay at "
              << format( "%7.1f %7.1f %7.1f", kDecay->position().x(), kDecay->position().y(), kDecay->position().z() )
              << endmsg;
      debug() << " pion momentum = " << part->momentum().R() / Gaudi::Units::GeV << " GeV" << endmsg;

      //== Fill ntuple
      double pz    = kShort->momentum().Z();
      double kSlXT = kShort->momentum().X() / pz;
      double kSlYT = kShort->momentum().Y() / pz;

      Tuple tTrack = nTuple( m_tupleName.value(), m_tupleName.value() );

      m_nTrack++;

      tTrack->column( "KVertX", kDecay->position().X() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "KVertY", kDecay->position().Y() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "KVertZ", kDecay->position().Z() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "KMoment", kShort->momentum().R() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "KSlopeX", kSlXT ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "KSlopeY", kSlYT ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      //== Fit the Pi+
      pz = part->momentum().Z();
      tTrack->column( "moment", part->momentum().R() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "slopeX", part->momentum().X() / pz ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "slopeY", part->momentum().Y() / pz ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      //== Fit the TT area
      if ( msgLevel( MSG::DEBUG ) ) debug() << "  UT: ";
      const auto fitUTX = m_fitTool->fitLine( UTHits, IPrFitTool::XY::X, m_zTT1 );
      const auto fitUTY = m_fitTool->fitLine( UTHits, IPrFitTool::XY::Y, m_zTT1 );
      if ( !fitUTX || !fitUTY ) {
        err() << "UT fit matrix is singular" << endmsg;
        continue;
      }

      const auto& [axt, bxt] = *fitUTX;
      const auto& [ayt, byt] = *fitUTY;

      if ( msgLevel( MSG::DEBUG ) )
        debug() << format( " x %7.1f tx %7.4f   y %7.1f ty %7.4f ", axt, bxt, ayt, byt ) << endmsg;
      ;
      tTrack->column( "axt", axt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "bxt", bxt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "ayt", ayt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "byt", byt ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      std::vector<Gaudi::XYZPoint>::const_iterator itP;
      if ( msgLevel( MSG::DEBUG ) ) {
        for ( itP = UTHits.begin(); UTHits.end() > itP; itP++ ) {
          const double dz = ( *itP ).z() - m_zTT1;
          debug() << format( "    : %7.1f %7.1f %7.1f  dx %7.3f  dy %7.3f", ( *itP ).x(), ( *itP ).y(), ( *itP ).z(),
                             ( *itP ).x() - ( axt + bxt * dz ), ( *itP ).y() - ( ayt + byt * dz ) )
                  << endmsg;
        }
      }

      const auto fitTrX = m_fitTool->fitCubic( trHits, IPrFitTool::XY::X, m_zRef );
      const auto fitTrY = m_fitTool->fitLine( trHits, IPrFitTool::XY::Y, m_zRef );
      if ( !fitTrX || !fitTrY ) {
        err() << "Track fit matrix is singular" << endmsg;
        continue;
      }

      const auto& [ax, bx, cx, dx] = *fitTrX;
      const auto& [ay, by]         = *fitTrY;

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

      double zMagnet = ( axt - bxt * m_zTT1 - ( ax - bx * m_zRef ) ) / ( bx - bxt );
      tTrack->column( "zMagnet", zMagnet ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      m_zMagPar.value().setFun( 0, 1. );
      m_zMagPar.value().setFun( 1, by * by );
      m_zMagPar.value().setFun( 2, bx * bx );

      double zEst = m_zMagPar.value().sum();
      tTrack->column( "zEst", zEst ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      m_zMagPar.value().addEvent( zMagnet - zEst );

      double bytPred = by + 5. * by * fabs( by ) * ( bx - bxt ) * ( bx - bxt );
      double aytPred = ay + ( zMagnet - m_zRef ) * by + ( m_zTT1 - zMagnet ) * bytPred;
      aytPred -= 2000 * by * ( bx - bxt ) * ( bx - bxt );

      tTrack->column( "bytPred", bytPred ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "aytPred", aytPred ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      tTrack->column( "xMagnet", axt + ( zMagnet - m_zTT1 ) * bxt )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->column( "yMagnet", ayt + ( zMagnet - m_zTT1 ) * byt )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      m_momPar.value().setFun( 0, 1. );
      m_momPar.value().setFun( 1, bx * bx );
      m_momPar.value().setFun( 2, by * by );

      double dSlope = fabs( bx - bxt );
      double dp     = dSlope * part->momentum().R() - m_momPar.value().sum();

      m_momPar.value().addEvent( dp );
      tTrack->column( "PredMom", m_momPar.value().sum() / dSlope )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      double xMag = ax + ( zEst - m_zRef ) * bx;
      tTrack->column( "xPred", xMag * m_zTT1 / zEst ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      tTrack->write().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      debug() << format( "  zMag%7.1f xMag%7.1f yMag%7.1f back%7.1f", zMagnet, axt + ( zMagnet - m_zTT1 ) * bxt,
                         ayt + ( zMagnet - m_zTT1 ) * byt, ay + ( zMagnet - m_zRef ) * by )
              << endmsg;
    }
  }
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrFitKsParams::finalize() {

  debug() << "==> Finalize" << endmsg;

  MsgStream& msg = info() << "============================================" << endmsg;
  msg << "  Processed " << m_nEvent << " events and " << m_nTrack << " tracks. " << endmsg;
  msg << "============================================" << endmsg;
  m_zMagPar.value().updateParameters( msg );
  m_momPar.value().updateParameters( msg );

  std::cout << std::endl;
  m_zMagPar.value().printPythonParams( name() );
  m_momPar.value().printPythonParams( name() );
  std::cout << std::endl;

  m_zMagPar.value().printParams( "PatKShort" );
  m_momPar.value().printParams( "PatKShort" );
  std::cout << std::endl;

  return GaudiTupleAlg::finalize(); // must be called after all other actions
}
//=============================================================================
