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
#include "Event/MCVertex.h"
#include "Event/RecVertex.h"
#include "LHCbAlgs/Consumer.h"

/**
 *  @author Olivier Callot
 *  @date   2011-11-16
 */
class VertexChecker
    : public LHCb::Algorithm::Consumer<void( LHCb::MCVertices const&, LHCb::RecVertices const&, LHCb::Tracks const& )> {
public:
  VertexChecker( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override; ///< Algorithm initialization
  void       operator()( LHCb::MCVertices const&, LHCb::RecVertices const&,
                   LHCb::Tracks const& ) const override; ///< Algorithm execution
  StatusCode finalize() override;                              ///< Algorithm finalization

  int mcBin( LHCb::MCVertex* vert ) const {
    int nPart  = vert->products().size();
    int binNum = nPart / m_binSizeMC;
    if ( binNum >= m_nbBin ) binNum = m_nbBin - 1;
    return binNum;
  }

  int recBin( LHCb::RecVertex* vert ) const {
    int nPart  = vert->tracks().size();
    int binNum = nPart / m_binSizeRec;
    if ( binNum >= m_nbBin ) binNum = m_nbBin - 1;
    return binNum;
  }

private:
  Gaudi::Property<std::string> m_inputLocation{this, "InputLocation", LHCb::RecVertexLocation::Primary};
  Gaudi::Property<double>      m_deltaZForMatch{this, "DeltaZForMatch", 1.000 * Gaudi::Units::mm};
  Gaudi::Property<double>      m_minIPForTrack{this, "MinIPForTrack", 0.150 * Gaudi::Units::mm};
  Gaudi::Property<double>      m_maxIPForTrack{this, "MaxIPForTrack", 3.000 * Gaudi::Units::mm};
  Gaudi::Property<double>      m_maxRadius{this, "MaxRadius", 3.000 * Gaudi::Units::mm};
  Gaudi::Property<int>         m_nbBin{this, "NbBin", 10};
  Gaudi::Property<int>         m_binSizeMC{this, "BinSizeMC", 20};
  Gaudi::Property<int>         m_binSizeRec{this, "BinSizeRec", 5};

  // TODO: replace with Gaudi::Accumulator
  mutable std::vector<int> m_mcVertices;
  mutable std::vector<int> m_mcFound;
  mutable std::vector<int> m_recVertices;
  mutable std::vector<int> m_recFake;

  mutable double m_s0;

  mutable double m_sx;
  mutable double m_sx2;

  mutable double m_sy;
  mutable double m_sy2;

  mutable double m_sz;
  mutable double m_sz2;

  mutable int m_nbLargeIP;
  mutable int m_nEvent;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( VertexChecker )
namespace {
  struct LowerByZ {
    bool operator()( const LHCb::MCVertex* lhs, const LHCb::MCVertex* rhs ) const {
      return lhs->position().z() < rhs->position().z();
    }
  };
} // namespace

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
VertexChecker::VertexChecker( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer{name,
               pSvcLocator,
               {KeyValue{"MCVertexLocation", LHCb::MCVertexLocation::Default},
                KeyValue{"RecVertexLocation", LHCb::RecVertexLocation::Primary},
                KeyValue{"TrackLocation", LHCb::TrackLocation::Velo}}} {}
//=============================================================================
// Initialization
//=============================================================================
StatusCode VertexChecker::initialize() {
  return Consumer::initialize().andThen( [&] {
    m_mcVertices.resize( m_nbBin, 0 );
    m_mcFound.resize( m_nbBin, 0 );
    m_recVertices.resize( m_nbBin, 0 );
    m_recFake.resize( m_nbBin, 0 );
    m_s0        = 0.;
    m_sx        = 0.;
    m_sx2       = 0.;
    m_sy        = 0.;
    m_sy2       = 0.;
    m_sz        = 0.;
    m_sz2       = 0.;
    m_nbLargeIP = 0;
    m_nEvent    = 0;
  } );
}

//=============================================================================
// Main execution
//=============================================================================
void VertexChecker::operator()( LHCb::MCVertices const& mcVertices, LHCb::RecVertices const& myVertices,
                                LHCb::Tracks const& veloTracks ) const {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;
  bool debug = msgLevel( MSG::DEBUG );

  std::vector<LHCb::MCVertex*> mcPvs;
  std::copy_if( mcVertices.begin(), mcVertices.end(), std::back_inserter( mcPvs ), [&]( const auto& v ) {
    return v->isPrimary() && ( v->products().size() > 2 ) && ( v->position().rho() <= m_maxRadius );
  } );

  std::sort( mcPvs.begin(), mcPvs.end(), LowerByZ() );

  auto itMc = mcPvs.begin();
  auto itPv = myVertices.begin();
  while ( itMc != mcPvs.end() && itPv != myVertices.end() ) {
    double dx = ( *itMc )->position().x() - ( *itPv )->position().x();
    double dy = ( *itMc )->position().y() - ( *itPv )->position().y();
    double dz = ( *itMc )->position().z() - ( *itPv )->position().z();
    if ( fabs( dz ) < m_deltaZForMatch ) {
      if ( debug )
        info() << format( " === Matched   dx%7.3f dy%7.3f  z%7.3f  dz%7.3f", dx, dy, ( *itMc )->position().z(), dz )
               << endmsg;
      int n = mcBin( *itMc );
      m_mcVertices[n]++;
      m_mcFound[n]++;
      int m = recBin( *itPv );
      m_recVertices[m]++;
      ++itMc;
      ++itPv;
      m_s0 += 1.;
      m_sx += dx;
      m_sx2 += dx * dx;
      m_sy += dy;
      m_sy2 += dy * dy;
      m_sz += dz;
      m_sz2 += dz * dz;

    } else if ( dz > 0 ) {
      if ( ( *itPv )->position().rho() < m_maxRadius ) {
        if ( debug )
          info() << format( " +++ fake reco  x%7.3f  y%7.3f  z%7.3f   n %3d", ( *itPv )->position().x(),
                            ( *itPv )->position().y(), ( *itPv )->position().z(), ( *itPv )->tracks().size() )
                 << endmsg;
        int m = recBin( *itPv );
        m_recVertices[m]++;
        m_recFake[m]++;
      }
      ++itPv;
    } else {
      if ( debug )
        info() << format( " --- miss reco  x%7.3f  y%7.3f  z%7.3f   n %3d", ( *itMc )->position().x(),
                          ( *itMc )->position().y(), ( *itMc )->position().z(), ( *itMc )->products().size() )
               << endmsg;
      int n = mcBin( *itMc );
      m_mcVertices[n]++;
      ++itMc;
    }
  }

  for ( ; itMc != mcPvs.end(); ++itMc ) {
    if ( debug )
      info() << format( " --- miss reco  x%7.3f  y%7.3f  z%7.3f   n %3d", ( *itMc )->position().x(),
                        ( *itMc )->position().y(), ( *itMc )->position().z(), ( *itMc )->products().size() )
             << endmsg;
    int n = mcBin( *itMc );
    m_mcVertices[n]++;
  }

  for ( ; itPv != myVertices.end(); ++itPv ) {
    if ( ( *itPv )->position().rho() < m_maxRadius ) {
      if ( debug )
        info() << format( " +++ fake reco  x%7.3f  y%7.3f  z%7.3f   n %3d", ( *itPv )->position().x(),
                          ( *itPv )->position().y(), ( *itPv )->position().z(), ( *itPv )->tracks().size() )
               << endmsg;
      int m = recBin( *itPv );
      m_recVertices[m]++;
      m_recFake[m]++;
    } else {
      info() << format( " ??? vertex at large rho   x%7.3f  y%7.3f  z%7.3f   n %3d", ( *itPv )->position().x(),
                        ( *itPv )->position().y(), ( *itPv )->position().z(), ( *itPv )->tracks().size() )
             << endmsg;
    }
  }

  // Number of tracks with IP between min and max values
  for ( const auto& itT : veloTracks ) {
    if ( itT->isVeloBackward() ) continue;
    Gaudi::XYZPoint  point = itT->position();
    Gaudi::XYZVector dir   = itT->slopes();
    auto [lowestIP, best]  = std::accumulate(
        myVertices.begin(), myVertices.end(), std::pair{1.e9, static_cast<const LHCb::RecVertex*>( nullptr )},
        [&]( auto current, const auto* v ) {
          double dx  = point.x() + ( v->position().z() - point.z() ) * dir.x() - v->position().x();
          double dy  = point.y() + ( v->position().z() - point.z() ) * dir.y() - v->position().y();
          double ip2 = dx * dx + dy * dy;
          return ( !current.second || ip2 < current.first ) ? std::pair{ip2, v} : current;
        } );
    lowestIP = sqrt( lowestIP );
    if ( lowestIP > m_minIPForTrack && lowestIP < m_maxIPForTrack ) {
      if ( debug )
        info() << format( "Large IP %6.3f zTr %7.3f zBest %7.3f", lowestIP, point.z(), best->position().z() ) << endmsg;
      m_nbLargeIP++;
    }
  }
  m_nEvent++;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode VertexChecker::finalize() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Finalize" << endmsg;

  info() << "=== MC Vertex reconstruction efficiency ===" << endmsg;
  int    totVertices = 0;
  int    totFound    = 0;
  double eff         = 0.;

  for ( unsigned int kk = 0; m_mcVertices.size() > kk; ++kk ) {
    totVertices += m_mcVertices[kk];
    totFound += m_mcFound[kk];
    eff = 0.;
    if ( 0 < m_mcVertices[kk] ) { eff = 100. * m_mcFound[kk] / double( m_mcVertices[kk] ); }
    if ( kk != m_mcVertices.size() - 1 ) {
      info() << format( "   MC %3d to %3d parts:  N %6d found %6d  eff %6.2f %%", m_binSizeMC * kk,
                        m_binSizeMC * ( kk + 1 ), m_mcVertices[kk], m_mcFound[kk], eff )
             << endmsg;
    } else {
      info() << format( "   MC   over %3d parts:  N %6d found %6d  eff %6.2f %%", m_binSizeMC * kk, m_mcVertices[kk],
                        m_mcFound[kk], eff )
             << endmsg;
    }
  }
  eff = 100. * totFound / double( totVertices );
  info() << format( "   MC total              N %6d found %6d  eff %6.2f %%", totVertices, totFound, eff ) << endmsg
         << "=== Fake rate of reconstructed vertices ===" << endmsg;
  int totRec  = 0;
  int totFake = 0;
  for ( unsigned int kk = 0; m_recVertices.size() > kk; ++kk ) {
    totRec += m_recVertices[kk];
    totFake += m_recFake[kk];
    eff = 0.;
    if ( 0 < m_recVertices[kk] ) { eff = 100. * m_recFake[kk] / double( m_recVertices[kk] ); }
    if ( kk == m_recVertices.size() - 1 ) {
      info() << format( "  REC   over %3d tracks: N %6d fake %6d  ghost %6.2f %%", 5 * kk, m_recVertices[kk],
                        m_recFake[kk], eff )
             << endmsg;
    } else {
      info() << format( "  REC %3d to %3d tracks: N %6d fake %6d  ghost %6.2f %%", 5 * kk, 5 * ( kk + 1 ),
                        m_recVertices[kk], m_recFake[kk], eff )
             << endmsg;
    }
  }
  eff = 100. * totFake / double( totRec );
  info() << format( "  REC total              N %6d fake %6d  ghost %6.2f %%", totRec, totFake, eff ) << endmsg;

  double meanX = m_sx / m_s0;
  info() << format( "Distance to MC: <x> %7.3f sx %7.3f", meanX, sqrt( m_sx2 / m_s0 - meanX * meanX ) ) << endmsg;
  double meanY = m_sy / m_s0;
  info() << format( "Distance to MC: <y> %7.3f sy %7.3f", meanY, sqrt( m_sy2 / m_s0 - meanY * meanY ) ) << endmsg;
  double meanZ = m_sz / m_s0;
  info() << format( "Distance to MC: <z> %7.3f sz %7.3f", meanZ, sqrt( m_sz2 / m_s0 - meanZ * meanZ ) ) << endmsg;

  double trackRate = double( m_nbLargeIP ) / double( m_nEvent );
  info() << format( "Number of large ( %5.3f to %5.3f mm) IP tracks/event = %7.2f", m_minIPForTrack.value(),
                    m_maxIPForTrack.value(), trackRate )
         << endmsg;

  return Consumer::finalize(); // must be called after all other actions
}
