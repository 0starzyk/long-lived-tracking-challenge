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

#include "MuonDet/DeMuonDetector.h"
#include "MuonDet/MuonNamespace.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"

#include "Event/MuonCoord.h"
#include "Event/State.h"
#include "Event/Track.h"

#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "LHCbAlgs/Consumer.h"

#include "AIDA/IHistogram1D.h"

namespace {
  struct Cache {
    double zM1, MAXsizeX, MAXsizeY;
  };
  struct MuonHit {
    const LHCb::MuonCoord* coord;
    double                 x, y, z;
    double                 dx, dy, dz;
  };
} // namespace

namespace LHCb {
  /**
   *  @author Stefania Vecchi
   *  @date   2010-01-22
   */
  class TrackMuonMatchMonitor
      : public LHCb::Algorithm::Consumer<
            void( Track::Range const&, MuonCoords const&, DetectorElement const&, DeMuonDetector const&, Cache const& ),
            DetDesc::usesBaseAndConditions<GaudiHistoAlg, DetectorElement, DeMuonDetector, Cache>> {
  public:
    TrackMuonMatchMonitor( const std::string& name, ISvcLocator* pSvcLoc )
        : Consumer{name,
                   pSvcLoc,
                   {KeyValue{"TracksLocation", TrackLocation::Default}, KeyValue{"MuonCoords", "Raw/Muon/Coords"},
                    KeyValue{"StandardGeometryTop", standard_geometry_top},
                    KeyValue{"MuonDetectorPath", DeMuonLocation::Default}, KeyValue{"CachePath", name + "_Cache"}}} {}

    StatusCode initialize() override; ///< Algorithm initialization
    void       operator()( Track::Range const&, MuonCoords const&, DetectorElement const&, DeMuonDetector const&,
                     Cache const& ) const override;

  private:
    ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackLinearExtrapolator"};

    Gaudi::Property<int>    m_iMS{this, "WhichStation", 0};
    Gaudi::Property<double> m_maxErrX{this, "MaxErrX", 5 * Gaudi::Units::mm};
    Gaudi::Property<double> m_maxErrY{this, "MaxErrY", 20 * Gaudi::Units::mm};

    static constexpr int nREGIONS = 4;
    AIDA::IHistogram1D * m_resx_a[nREGIONS], *m_resy_a[nREGIONS], *m_resx_c[nREGIONS], *m_resy_c[nREGIONS];
    double               m_hisxmax[nREGIONS];
  };

  // Declaration of the Algorithm Factory
  DECLARE_COMPONENT_WITH_ID( TrackMuonMatchMonitor, "TrackMuonMatchMonitor" )

} // namespace LHCb

StatusCode LHCb::TrackMuonMatchMonitor::initialize() {
  return Consumer::initialize().andThen( [&]() {
    addConditionDerivation(
        {inputLocation<DeMuonDetector>()}, inputLocation<Cache>(), [this]( DeMuonDetector const& muonDet ) -> Cache {
          return {muonDet.getStationZ( m_iMS ), muonDet.getOuterX( m_iMS ), muonDet.getOuterY( m_iMS )};
        } );
    std::string name;
    setHistoTopDir( "Track/" );
    for ( int iR = 0; iR < nREGIONS; ++iR ) {
      // in the "signal" region +/- 6 error units respect the track extrapolation point
      unsigned int nbin = 100;
      double       max  = 80. + 40. * float( iR );
      double       min  = -max;
      name              = "resX_ASide_M1R" + std::to_string( iR + 1 );
      m_resx_a[iR]      = book1D( name, name, min, max, nbin );
      name              = "resY_ASide_M1R" + std::to_string( iR + 1 );
      m_resy_a[iR]      = book1D( name, name, min, max, nbin );
      name              = "resX_CSide_M1R" + std::to_string( iR + 1 );
      m_resx_c[iR]      = book1D( name, name, min, max, nbin );
      name              = "resY_CSide_M1R" + std::to_string( iR + 1 );
      m_resy_c[iR]      = book1D( name, name, min, max, nbin );
      m_hisxmax[iR]     = max;
    }

    return StatusCode::SUCCESS;
  } );
}

void LHCb::TrackMuonMatchMonitor::operator()( LHCb::Track::Range const& tTracks, LHCb::MuonCoords const& coords,
                                              DetectorElement const& lhcb, DeMuonDetector const& muonDet,
                                              Cache const& cache ) const {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;
  if ( tTracks.empty() ) return;

  if ( coords.empty() ) {
    if ( msgLevel( MSG::DEBUG ) ) debug() << " No hits retrieved , skip event" << endmsg;
    return;
  }

  // cache the position of the hits. that saves a lot of time.
  std::vector<MuonHit> muonhits;
  muonhits.reserve( coords.size() );
  for ( const auto& coord : coords ) {
    if ( m_iMS == int( coord->key().station() ) ) { // only the Chosen station
      auto pos = muonDet.position( coord->key() );
      if ( pos ) {
        auto& hit = muonhits.emplace_back();
        hit.coord = coord;
        hit.x     = pos->x();
        hit.y     = pos->y();
        hit.z     = pos->z();
        hit.dx    = pos->dX();
        hit.dy    = pos->dY();
        hit.dz    = pos->dZ();
      }
    }
  }

  if ( msgLevel( MSG::DEBUG ) ) debug() << " Found " << tTracks.size() << " tracks in the container " << endmsg;
  for ( const Track* track : tTracks ) {
    if ( track->hasT() && track->chi2PerDoF() < 5 && track->p() > 1 * Gaudi::Units::GeV ) {

      State      stateAtM1;
      StatusCode sc = m_extrapolator->propagate( *track, cache.zM1, stateAtM1, *lhcb.geometry() );

      if ( sc.isSuccess() && std::abs( stateAtM1.x() ) < cache.MAXsizeX && std::abs( stateAtM1.y() ) < cache.MAXsizeY &&
           std::sqrt( stateAtM1.errX2() ) < m_maxErrX && std::sqrt( stateAtM1.errY2() ) < m_maxErrY ) {

        for ( const MuonHit& hit : muonhits ) {

          int    region = hit.coord->key().region();
          double deltaZ = hit.z - stateAtM1.z();
          double deltaX = hit.x - ( stateAtM1.x() + stateAtM1.tx() * deltaZ );
          double deltaY = hit.y - ( stateAtM1.y() + stateAtM1.ty() * deltaZ );

          if ( std::abs( deltaX ) < m_hisxmax[region] && std::abs( deltaY ) < m_hisxmax[region] ) {

            AIDA::IHistogram1D *tempx, *tempy;

            tempx = hit.x > 0 ? m_resx_a[region] : m_resx_c[region];
            tempy = hit.x > 0 ? m_resy_a[region] : m_resy_c[region];

            tempx->fill( deltaX ); // X residuals on the same Z as the hit
            tempy->fill( deltaY ); // Y residuals on the same Z as the hit
          }
        }
      }
    }
  }
}
