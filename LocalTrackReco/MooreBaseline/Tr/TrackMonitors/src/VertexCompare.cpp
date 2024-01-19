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
#include "AIDA/IHistogram1D.h"
#include "Event/MCVertex.h"
#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiAlg/Tuples.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiUtils/HistoStats.h"
#include "MCInterfaces/IForcedBDecayTool.h"
#include <Event/MCTrackInfo.h>
#include <Linker/LinkedTo.h>
//-----------------------------------------------------------------------------
// Implementation file for class : VertexCompare
//-----------------------------------------------------------------------------

struct MCPVInfo {
  LHCb::MCVertex*                pMCPV;             // pointer to MC PV
  int                            nRecTracks;        // number of reconstructed tracks from this MCPV
  int                            nRecBackTracks;    // number of reconstructed backward tracks
  int                            indexRecPVInfo;    // index to reconstructed PVInfo (-1 if not reco)
  int                            nCorrectTracks;    // correct tracks belonging to reconstructed PV
  int                            multClosestMCPV;   // multiplicity of closest reconstructable MCPV
  double                         distToClosestMCPV; // distance to closest reconstructible MCPV
  std::vector<LHCb::MCParticle*> m_mcPartInMCPV;
  std::vector<LHCb::Track*>      m_recTracksInMCPV;
};

struct RecPVInfo {
public:
  int              nTracks;       // number of tracks in a vertex
  int              nBackTracks;   // number of backward tracks in a vertex
  Gaudi::XYZPoint  position;      // position
  Gaudi::XYZPoint  positionSigma; // position sigmas
  int              indexMCPVInfo; // index to MCPVInfo
  LHCb::RecVertex* pRECPV;        // pointer to REC PV
  bool             rec1, rec2;
};

class VertexCompare : public GaudiTupleAlg {
public:
  /// Standard constructor
  VertexCompare( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution
  StatusCode finalize() override;   ///< Algorithm finalization

private:
  bool debugLevel() const { return msgLevel( MSG::DEBUG ) || msgLevel( MSG::VERBOSE ); }

  ToolHandle<IForcedBDecayTool> m_forcedBtool{this, "ForcedBDecayTool", "ForcedBDecayTool"};

  //   Gaudi::Property<int> m_nTracksToBeRecble{this,"nTracksToBeRecble",  5};
  Gaudi::Property<bool> m_produceHistogram{this, "produceHistogram", true};
  Gaudi::Property<bool> m_produceNtuple{this, "produceNtuple", true};
  //   Gaudi::Property<double> m_dzIsolated{this,"dzIsolated",         10.0 * Gaudi::Units::mm};
  //   Gaudi::Property<std::string> m_inputTracksName{this,"inputTracksName",     LHCb::TrackLocation::Default};
  Gaudi::Property<std::string> m_inputVerticesName1{this, "inputVerticesName1", LHCb::RecVertexLocation::Primary};
  Gaudi::Property<std::string> m_inputVerticesName2{this, "inputVerticesName2", LHCb::RecVertexLocation::Primary};
  //   std::string m_pvSeedingName;
  Gaudi::Accumulators::Counter<> m_nVtx{this, "Number of pairs of vertices in processed events"};

  void printRat( std::string mes, int a, int b );
  //  bool getInputTracks( std::vector<LHCb::Track*>& vecOfTracks,  std::string& trackLoc);
  bool getInputVertices( std::vector<LHCb::RecVertex*>& vecOfVertices1, std::vector<LHCb::RecVertex*>& vecOfVertices2 );
  void matchByDistance( std::vector<LHCb::RecVertex>& verthalf, std::vector<LHCb::RecVertex>& vertfull,
                        std::vector<int>& link );
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( VertexCompare )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
VertexCompare::VertexCompare( const std::string& name, ISvcLocator* pSvcLocator )
    : GaudiTupleAlg( name, pSvcLocator ) {}
//=============================================================================
// Initialization
//=============================================================================
StatusCode VertexCompare::initialize() {
  return GaudiTupleAlg::initialize().andThen( [&] { m_nVtx.reset(); } );
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode VertexCompare::execute() {
  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  std::vector<LHCb::RecVertex*> vecOfVertices1;
  std::vector<LHCb::RecVertex*> vecOfVertices2;

  bool vertices_ok = getInputVertices( vecOfVertices1, vecOfVertices2 );
  if ( !vertices_ok ) return StatusCode::SUCCESS; // return SUCCESS anyway

  std::vector<LHCb::RecVertex> fullVrt1, fullVrt2;

  if ( debugLevel() ) debug() << "  Vtx Properities       x       y       z      chi2/ndof     ntracks" << endmsg;
  // Fill reconstructed PV info
  for ( LHCb::RecVertex* pv : vecOfVertices1 ) {
    fullVrt1.push_back( *pv );
    //     int nTracks      = pv->tracks().size();
    if ( debugLevel() ) debug() << m_inputVerticesName1 << endmsg;
    if ( debugLevel() )
      debug() << "              " << pv->position().x() << "   " << pv->position().y() << "   " << pv->position().z()
              << "   " << pv->chi2PerDoF() << "   " << pv->tracks().size() << endmsg;
  } // end of loop over vertices1

  for ( LHCb::RecVertex* pv : vecOfVertices2 ) {
    fullVrt2.push_back( *pv );
    //     int nTracks      = pv->tracks().size();
    if ( debugLevel() ) debug() << m_inputVerticesName2 << endmsg;
    if ( debugLevel() )
      debug() << "              " << pv->position().x() << "   " << pv->position().y() << "   " << pv->position().z()
              << "   " << pv->chi2PerDoF() << "   " << pv->tracks().size() << endmsg;
  } // end of loop over vertices2
  if ( debugLevel() ) debug() << "fullVrt1 size   " << fullVrt1.size() << endmsg;
  if ( debugLevel() ) debug() << "fullVrt2 size   " << fullVrt2.size() << endmsg;
  int size_diff = fullVrt1.size() - fullVrt2.size();
  if ( m_produceHistogram.value() ) { plot( double( size_diff ), 1001, "size_diff", -5.5, 5.5, 11 ); }
  if ( m_produceNtuple.value() ) {
    Tuple myTuple_evt = nTuple( 101, "PV_nTuple_evt", CLID_ColumnWiseTuple );
    myTuple_evt->column( "size_diff", double( size_diff ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    myTuple_evt->column( "size_1", double( fullVrt1.size() ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    myTuple_evt->column( "size_2", double( fullVrt2.size() ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    myTuple_evt->write().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  }

  int ntracks_part  = 0;
  int ntracks_part2 = 0;
  int dtracks       = 0;

  int              oIt = 0;
  std::vector<int> link;
  if ( fullVrt1.size() != 0 && fullVrt2.size() != 0 )
    matchByDistance( fullVrt1, fullVrt2, link );
  else
    return StatusCode::SUCCESS;
  for ( const LHCb::RecVertex& vrtf : fullVrt1 ) {
    ++m_nVtx;
    if ( link.at( oIt ) == -1 ) continue;
    Gaudi::SymMatrix3x3 covPV_part1 = vrtf.covMatrix();
    double              sigx_part1  = ( covPV_part1( 0, 0 ) );
    double              sigy_part1  = ( covPV_part1( 1, 1 ) );
    double              sigz_part1  = ( covPV_part1( 2, 2 ) );
    double              covxx1      = sigx_part1;
    double              covyy1      = sigy_part1;
    double              covzz1      = sigz_part1;
    double              covxy1      = ( covPV_part1( 0, 1 ) );
    double              covxz1      = ( covPV_part1( 0, 2 ) );
    double              covyz1      = ( covPV_part1( 1, 2 ) );

    Gaudi::SymMatrix3x3 covPV_part2 = fullVrt2.at( link.at( oIt ) ).covMatrix();
    double              sigx_part2  = ( covPV_part2( 0, 0 ) );
    double              sigy_part2  = ( covPV_part2( 1, 1 ) );
    double              sigz_part2  = ( covPV_part2( 2, 2 ) );

    double covxx2 = sigx_part2;
    double covyy2 = sigy_part2;
    double covzz2 = sigz_part2;
    double covxy2 = ( covPV_part2( 0, 1 ) );
    double covxz2 = ( covPV_part2( 0, 2 ) );
    double covyz2 = ( covPV_part2( 1, 2 ) );

    double x1 = vrtf.position().x();
    double y1 = vrtf.position().y();
    double z1 = vrtf.position().z();

    double x2 = fullVrt2.at( link.at( oIt ) ).position().x();
    double y2 = fullVrt2.at( link.at( oIt ) ).position().y();
    double z2 = fullVrt2.at( link.at( oIt ) ).position().z();

    double dx = vrtf.position().x() - fullVrt2.at( link.at( oIt ) ).position().x();
    double dy = vrtf.position().y() - fullVrt2.at( link.at( oIt ) ).position().y();
    double dz = vrtf.position().z() - fullVrt2.at( link.at( oIt ) ).position().z();

    double errx   = sqrt( sigx_part1 + sigx_part2 );
    double erry   = sqrt( sigy_part1 + sigy_part2 );
    double errz   = sqrt( sigz_part1 + sigz_part2 );
    ntracks_part  = vrtf.tracks().size();
    ntracks_part2 = fullVrt2.at( link.at( oIt ) ).tracks().size();
    dtracks       = ntracks_part - ntracks_part2;

    if ( m_produceHistogram.value() ) {
      plot( dx, 1021, "dx", -0.25, 0.25, 50 );
      plot( dy, 1022, "dy", -0.25, 0.25, 50 );
      plot( dz, 1023, "dz", -1.5, 1.5, 60 );
      plot( x1, 2021, "x1", -0.25, 0.25, 50 );
      plot( y1, 2022, "y1", -0.25, 0.25, 50 );
      plot( z1, 2023, "z1", -100., 100., 50 );
      plot( x2, 3021, "x2", -0.25, 0.25, 50 );
      plot( y2, 3022, "y2", -0.25, 0.25, 50 );
      plot( z2, 3023, "z2", -100., 100., 50 );
      plot( std::sqrt( sigx_part1 ), 4011, "x err 1", 0., .1, 50 );
      plot( std::sqrt( sigy_part1 ), 4012, "y err 1", 0., .1, 50 );
      plot( std::sqrt( sigz_part1 ), 4013, "z err 1", 0., .5, 50 );
      plot( std::sqrt( sigx_part2 ), 4021, "x err 2", 0., .1, 50 );
      plot( std::sqrt( sigy_part2 ), 4022, "y err 2", 0., .1, 50 );
      plot( std::sqrt( sigz_part2 ), 4023, "z err 2", 0., .5, 50 );
      plot( dx / errx, 1031, "pullx", -5., 5., 50 );
      plot( dy / erry, 1032, "pully", -5., 5., 50 );
      plot( dz / errz, 1033, "pullz", -5., 5., 50 );
      plot( double( ntracks_part ), 1041, "ntracks_part1", 0., 150., 50 );
      plot( double( ntracks_part2 ), 1042, "ntracks_part2", 0., 150., 50 );
      plot( double( dtracks ), 1043, "dtracks", 0., 150., 50 );
      // dz to get false vertex rate from data
      if ( vecOfVertices1.size() > 1 ) {
        for ( unsigned int i1 = 0; i1 < vecOfVertices1.size(); i1++ ) {
          for ( unsigned int i2 = 0; i2 < vecOfVertices1.size(); i2++ ) {
            if ( i2 != i1 ) {
              double vdz = vecOfVertices1[i1]->position().z() - vecOfVertices1[i2]->position().z();
              plot( vdz, 1201, "dz vertices 1", -150., 150., 100 );
            }
          }
        }
      }
      if ( vecOfVertices2.size() > 1 ) {
        for ( unsigned int i1 = 0; i1 < vecOfVertices2.size(); i1++ ) {
          for ( unsigned int i2 = 0; i2 < vecOfVertices2.size(); i2++ ) {
            if ( i2 != i1 ) {
              double vdz = vecOfVertices2[i1]->position().z() - vecOfVertices2[i2]->position().z();
              plot( vdz, 1202, "dz vertices 2", -150., 150., 100 );
            }
          }
        }
      }
    }
    if ( m_produceNtuple.value() ) {
      Tuple myTuple = nTuple( 102, "PV_nTuple", CLID_ColumnWiseTuple );
      myTuple->column( "ntracks_part", double( ntracks_part ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "ntracks_part2", double( ntracks_part2 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "dtracks", double( dtracks ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "dx", dx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "dy", dy ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "dz", dz ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "x1", x1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "y1", y1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "z1", z1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "x2", x2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "y2", y2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "z2", z2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "errx", errx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "erry", erry ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "errz", errz ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covxx1", covxx1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covyy1", covyy1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covzz1", covzz1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covxy1", covxy1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covxz1", covxz1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covyz1", covyz1 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covxx2", covxx2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covyy2", covyy2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covzz2", covzz2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covxy2", covxy2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covxz2", covxz2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      myTuple->column( "covyz2", covyz2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      myTuple->write().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    }
    oIt++;
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode VertexCompare::finalize() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Finalize" << endmsg;

  info() << " ============================================" << endmsg;
  info() << " Efficiencies for reconstructed vertices:    " << endmsg;
  info() << " ============================================" << endmsg;
  info() << " " << endmsg;

  //   info() << " PV is isolated if dz to closest reconstructible MC PV >  "
  //          << m_dzIsolated << " mm" << endmsg;
  //   std::string ff = "by counting tracks";
  // /*  if ( !m_matchByTracks )*/ ff = "by dz distance";
  //   info() << " Two splited vertices matched:  "
  //          <<  ff << endmsg;

  //   info() << " " << endmsg;
  //     printRat("All",       m_nPartVtx,       m_nVtx );

  const AIDA::IHistogram1D* dx    = histo( HistoID( 1021 ) );
  const AIDA::IHistogram1D* pullx = histo( HistoID( 1031 ) );
  const AIDA::IHistogram1D* dy    = histo( HistoID( 1022 ) );
  const AIDA::IHistogram1D* pully = histo( HistoID( 1032 ) );
  const AIDA::IHistogram1D* dz    = histo( HistoID( 1023 ) );
  const AIDA::IHistogram1D* pullz = histo( HistoID( 1033 ) );
  if ( dx ) {
    info() << "      ---------------------------------------" << endmsg;
    info() << "dx:    "
           << format( "mean =  %5.3f +/- %5.3f, RMS =  %5.3f +/- %5.3f", dx->mean(),
                      Gaudi::Utils::HistoStats::meanErr( dx ), dx->rms(), Gaudi::Utils::HistoStats::rmsErr( dx ) )
           << endmsg;
  }
  if ( dy ) {
    info() << "dy:    "
           << format( "mean =  %5.3f +/- %5.3f, RMS =  %5.3f +/- %5.3f", dy->mean(),
                      Gaudi::Utils::HistoStats::meanErr( dy ), dy->rms(), Gaudi::Utils::HistoStats::rmsErr( dy ) )
           << endmsg;
  }
  if ( dz ) {
    info() << "dz:    "
           << format( "mean =  %5.3f +/- %5.3f, RMS =  %5.3f +/- %5.3f", dz->mean(),
                      Gaudi::Utils::HistoStats::meanErr( dz ), dz->rms(), Gaudi::Utils::HistoStats::rmsErr( dz ) )
           << endmsg;
  }
  info() << "      ---------------------------------------" << endmsg;
  if ( pullx ) {
    info() << "pullx: "
           << format( "mean =  %5.3f +/- %5.3f, RMS =  %5.3f +/- %5.3f", pullx->mean(),
                      Gaudi::Utils::HistoStats::meanErr( pullx ), pullx->rms(),
                      Gaudi::Utils::HistoStats::rmsErr( pullx ) )
           << endmsg;
  }
  if ( pully ) {
    info() << "pully: "
           << format( "mean =  %5.3f +/- %5.3f, RMS =  %5.3f +/- %5.3f", pully->mean(),
                      Gaudi::Utils::HistoStats::meanErr( pully ), pully->rms(),
                      Gaudi::Utils::HistoStats::rmsErr( pully ) )
           << endmsg;
  }
  if ( pullz ) {
    info() << "pullz: "
           << format( "mean =  %5.3f +/- %5.3f, RMS =  %5.3f +/- %5.3f", pullz->mean(),
                      Gaudi::Utils::HistoStats::meanErr( pullz ), pullz->rms(),
                      Gaudi::Utils::HistoStats::rmsErr( pullz ) )
           << endmsg;
  }
  info() << " ============================================" << endmsg;
  //
  return GaudiTupleAlg::finalize(); // Must be called after all other actions
}

//=============================================================================
//  Match vertices by distance
//=============================================================================
void VertexCompare::matchByDistance( std::vector<LHCb::RecVertex>& vertfull, std::vector<LHCb::RecVertex>& verthalf,
                                     std::vector<int>& link ) {

  if ( verthalf.size() > vertfull.size() && debugLevel() ) debug() << "half.size > full.size" << endmsg;
  for ( int imc = 0; imc < (int)vertfull.size(); imc++ ) {
    //     if ( mcpvvec[imc].indexRecPVInfo  > -1) continue;
    double mindist  = 999999.;
    int    indexrec = -1;
    for ( int irec = 0; irec < (int)verthalf.size(); irec++ ) {
      if ( std::count( link.begin(), link.end(), irec ) != 0 ) continue;
      double dist = fabs( verthalf.at( irec ).position().z() - vertfull.at( imc ).position().z() );
      if ( dist < mindist ) {
        mindist  = dist;
        indexrec = irec;
      }
    }
    if ( debugLevel() ) debug() << "original vertex " << imc << " linked to " << indexrec << " half vertex." << endmsg;
    link.push_back( indexrec );
  }

  for ( int imc = 0; imc < (int)vertfull.size(); imc++ ) {
    int count = std::count( link.begin(), link.end(), imc );
    if ( count > 1 && debugLevel() ) debug() << "linked twice to vertex " << imc << endmsg;
  }
}

//=============================================================================
//  printRat
//=============================================================================
void VertexCompare::printRat( std::string mes, int a, int b ) {

  double rat = 0.;
  if ( b > 0 ) rat = 1.0 * a / b;

  // reformat message
  unsigned int len  = 20;
  std::string  pmes = mes;
  while ( pmes.length() < len ) { pmes += " "; }
  pmes += " : ";

  info() << pmes << format( " %6.3f ( %7d / %8d )", rat, a, b ) << endmsg;
}

//=============================================================================
//  Get input vertices
//=============================================================================
bool VertexCompare::getInputVertices( std::vector<LHCb::RecVertex*>& vecOfVertices1,
                                      std::vector<LHCb::RecVertex*>& vecOfVertices2 ) {

  std::string verticesName1;
  if ( m_inputVerticesName1 == "Offline" ) {
    verticesName1 = LHCb::RecVertexLocation::Primary;
  } else if ( m_inputVerticesName1 == "3D" ) {
    verticesName1 = LHCb::RecVertexLocation::Velo3D;
  } else if ( m_inputVerticesName1 == "2D" ) {
    verticesName1 = LHCb::RecVertexLocation::Velo2D;
  } else {
    verticesName1 = m_inputVerticesName1;
  }

  std::string verticesName2;
  if ( m_inputVerticesName2 == "Offline" ) {
    verticesName2 = LHCb::RecVertexLocation::Primary;
  } else if ( m_inputVerticesName2 == "3D" ) {
    verticesName2 = LHCb::RecVertexLocation::Velo3D;
  } else if ( m_inputVerticesName2 == "2D" ) {
    verticesName2 = LHCb::RecVertexLocation::Velo2D;
  } else {
    verticesName2 = m_inputVerticesName2;
  }

  if ( verticesName1 == "none" ) {
    debug() << " Vertices 1 not specified " << verticesName1 << endmsg;
    return false;
  }
  if ( verticesName2 == "none" ) {
    debug() << " Vertices 2 not specified " << verticesName2 << endmsg;
    return false;
  }

  LHCb::RecVertices* recoVertices1 = getIfExists<LHCb::RecVertices>( verticesName1 );
  if ( !recoVertices1 ) {
    debug() << " No vertices at " << verticesName1 << endmsg;
    return false;
  }

  LHCb::RecVertices* recoVertices2 = getIfExists<LHCb::RecVertices>( verticesName2 );
  if ( !recoVertices2 ) {
    debug() << " No vertices at " << verticesName2 << endmsg;
    return false;
  }

  std::copy( recoVertices1->begin(), recoVertices1->end(), std::back_inserter( vecOfVertices1 ) );
  std::copy( recoVertices2->begin(), recoVertices2->end(), std::back_inserter( vecOfVertices2 ) );

  return true;
}
