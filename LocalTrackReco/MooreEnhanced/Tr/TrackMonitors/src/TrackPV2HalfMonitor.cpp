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

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/ODIN.h"
#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "TrackInterfaces/IPVOfflineTool.h"
#include "TrackInterfaces/ITrackVertexer.h"
#include "TrackKernel/TrackStateVertex.h"

#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Consumer.h"
#include <Gaudi/Accumulators/Histogram.h>

#include <algorithm>
namespace GA  = Gaudi::Accumulators;
using GAFTYPE = double;
using GAAxis  = GA::Axis<GAFTYPE>;
using GAH1    = Gaudi::Accumulators::Histogram<1, GA::atomicity::full, GAFTYPE>;
using GAPr    = Gaudi::Accumulators::ProfileHistogram<1>;

class TrackPV2HalfMonitor
    : public LHCb::Algorithm::Consumer<void( LHCb::ODIN const&, LHCb::Track::Range const&, DetectorElement const& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DetectorElement>> {
public:
  /** Standard construtor */
  TrackPV2HalfMonitor( const std::string& name, ISvcLocator* pSvcLoc )
      : Consumer{name,
                 pSvcLoc,
                 {KeyValue{"ODINLocation", LHCb::ODINLocation::Default},
                  KeyValue{"TrackContainer", LHCb::TrackLocation::Default},
                  KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}}} {}

  /** Algorithm initialize */
  StatusCode initialize() override;

  /** Algorithm execute */
  void operator()( LHCb::ODIN const&, LHCb::Track::Range const&, DetectorElement const& ) const override;

private:
  PublicToolHandle<ITrackVertexer> m_vertexer{this, "TrackVertexer", "TrackVertexer"};
  ToolHandle<IPVOfflineTool>       m_toolpv{this, "PVOfflineTool", "PVOfflineTool"};
  Gaudi::Property<GAFTYPE>         m_zpvmin{this, "MinZPV", -20 * Gaudi::Units::cm};
  Gaudi::Property<GAFTYPE>         m_zpvmax{this, "MaxZPV", 20 * Gaudi::Units::cm};
  Gaudi::Property<GAFTYPE>         m_limpvx{this, "limPx", 2. * Gaudi::Units::mm};
  Gaudi::Property<GAFTYPE>         m_limpvy{this, "limPy", 1. * Gaudi::Units::mm};
  Gaudi::Property<GAFTYPE>         m_limpvz{this, "limPz", 150. * Gaudi::Units::mm};
  Gaudi::Property<GAFTYPE>         m_meanpvx{this, "meanPx", 0. * Gaudi::Units::mm};
  Gaudi::Property<GAFTYPE>         m_meanpvy{this, "meanPy", 0. * Gaudi::Units::mm};
  Gaudi::Property<GAFTYPE>         m_limdpvx{this, "limDPx", 0.5 * Gaudi::Units::mm};
  Gaudi::Property<GAFTYPE>         m_limdpvy{this, "limDPy", 0.5 * Gaudi::Units::mm};
  Gaudi::Property<GAFTYPE>         m_limdpvz{this, "limDPz", 1. * Gaudi::Units::mm};
  Gaudi::Property<GAFTYPE>         m_limchi2{this, "limChi2", 10.};
  Gaudi::Property<unsigned int>    m_nprbins{this, "NumProfileBins", 20};
  Gaudi::Property<unsigned int>    m_PV_trackmin{this, "MinNumTrPerPV", 5};

  // Left side histograms
  mutable std::optional<GAH1> m_numTracksPerPVLeft;
  mutable std::optional<GAH1> m_xPosPVLeft;
  mutable std::optional<GAH1> m_yPosPVLeft;
  mutable std::optional<GAH1> m_zPosPVLeft;
  mutable std::optional<GAH1> m_PVChi2Left;
  mutable std::optional<GAH1> m_refittedPVChi2Left;
  // right side histograms
  mutable std::optional<GAH1> m_numTracksPerPVRight;
  mutable std::optional<GAH1> m_xPosPVRight;
  mutable std::optional<GAH1> m_yPosPVRight;
  mutable std::optional<GAH1> m_zPosPVRight;
  mutable std::optional<GAH1> m_PVChi2Right;
  mutable std::optional<GAH1> m_refittedPVChi2Right;
  // Left-right histograms
  mutable std::optional<GAH1> m_deltaxPVLeftRight;
  mutable std::optional<GAH1> m_deltayPVLeftRight;
  mutable std::optional<GAH1> m_deltazPVLeftRight;
  mutable std::optional<GAH1> m_deltaxPVLeftRightPull;
  mutable std::optional<GAH1> m_deltayPVLeftRightPull;
  mutable std::optional<GAH1> m_deltazPVLeftRightPull;
  // EventTime histo
  mutable std::optional<GAH1> m_eventTimeHisto;
  // Left-right profiles
  mutable std::optional<GAPr> m_deltayzPVLeftRight;
  mutable std::optional<GAPr> m_deltaxzPVLeftRight;
  mutable std::optional<GAPr> m_deltazzPVLeftRight;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TrackPV2HalfMonitor )

//=============================================================================
// Initialization
//=============================================================================
StatusCode TrackPV2HalfMonitor::initialize() {
  return Consumer::initialize().andThen( [&] {
    if ( histoTopDir().empty() ) setHistoTopDir( "Track/" );
    m_numTracksPerPVLeft.emplace( this, "Left PV Num of track", "Left PV Num of track", GAAxis{100, 0.5, 100.5} );
    m_xPosPVLeft.emplace( this, "Left PV x position", "Left PV x position",
                          GAAxis{200, -m_limpvx.value() + m_meanpvx, m_limpvx.value() + m_meanpvx} );
    m_yPosPVLeft.emplace( this, "Left PV y position", "Left PV y position",
                          GAAxis{200, -m_limpvy.value() + m_meanpvy, m_limpvy.value() + m_meanpvy} );
    m_zPosPVLeft.emplace( this, "Left PV z position", "Left PV z position", GAAxis{200, -m_limpvz, m_limpvz} );
    m_PVChi2Left.emplace( this, "Left PV Chi2 per dof", "Left PV Chi2 per dof", GAAxis{100, 0, m_limchi2} );
    m_refittedPVChi2Left.emplace( this, "Left re-fitted PV Chi2 per dof re-fitted", "Left re-fitted PV Chi2 per dof",
                                  GAAxis{100, 0, m_limchi2} );
    m_numTracksPerPVRight.emplace( this, "Right PV Num of track", "Right PV Num of track", GAAxis{100, 0.5, 100.5} );
    m_xPosPVRight.emplace( this, "Right PV x position", "Right PV x position",
                           GAAxis{200, -m_limpvx.value() + m_meanpvx, m_limpvx.value() + m_meanpvx} );
    m_yPosPVRight.emplace( this, "Right PV y position", "Right PV y position",
                           GAAxis{200, -m_limpvy.value() + m_meanpvy, m_limpvy.value() + m_meanpvy} );
    m_zPosPVRight.emplace( this, "Right PV z position", "Right PV z position", GAAxis{200, -m_limpvz, m_limpvz} );
    m_PVChi2Right.emplace( this, "Right PV Chi2 per dof", "Right PV Chi2 per dof", GAAxis{100, 0, m_limchi2} );
    m_refittedPVChi2Right.emplace( this, "Right re-fitted PV Chi2 per dof re-fitted", "Right re-fitted PV Chi2 per dof",
                                   GAAxis{100, 0, m_limchi2} );
    // Left-right histograms
    m_deltaxPVLeftRight.emplace( this, "Left-Right PV delta x", "Left-Right PV delta x",
                                 GAAxis{100, -m_limdpvx, m_limdpvx} );
    m_deltayPVLeftRight.emplace( this, "Left-Right PV delta y", "Left-Right PV delta y",
                                 GAAxis{100, -m_limdpvy, m_limdpvy} );
    m_deltazPVLeftRight.emplace( this, "Left-Right PV delta z", "Left-Right PV delta z",
                                 GAAxis{100, -m_limdpvz, m_limdpvz} );
    m_deltaxPVLeftRightPull.emplace( this, "Left-Right PV delta x pull", "Left-Right PV delta x pull",
                                     GAAxis{100, -5, 5} );
    m_deltayPVLeftRightPull.emplace( this, "Left-Right PV delta y pull", "Left-Right PV delta y pull",
                                     GAAxis{100, -5, 5} );
    m_deltazPVLeftRightPull.emplace( this, "Left-Right PV delta z pull", "Left-Right PV delta z pull",
                                     GAAxis{100, -5, 5} );
    // EventTime histo
    m_eventTimeHisto.emplace( this, "EventTime TrackPV2HalfMonitor", "TimeMinute", GAAxis{1000, 0, 30000} );
    // Left-right profiles
    m_deltaxzPVLeftRight.emplace( this, "PV left-right delta x versus z", "PV left-right delta x versus z",
                                  GAAxis{m_nprbins, m_zpvmin, m_zpvmax} );
    m_deltayzPVLeftRight.emplace( this, "PV left-right delta y versus z", "PV left-right delta y versus z",
                                  GAAxis{m_nprbins, m_zpvmin, m_zpvmax} );
    m_deltazzPVLeftRight.emplace( this, "PV left-right delta z versus z", "PV left-right delta z versus z",
                                  GAAxis{m_nprbins, m_zpvmin, m_zpvmax} );
  } );
}

//=============================================================================
// Structure
//=============================================================================

namespace {
  std::vector<const LHCb::Track*> myconvert( const SmartRefVector<LHCb::Track>& tracks ) {
    return {tracks.begin(), tracks.end()};
  }

  template <class TrackContainer, class Predicate>
  std::vector<const LHCb::Track*> myselect( const TrackContainer& tracks, Predicate selector ) {
    std::vector<const LHCb::Track*> rc;
    std::copy_if( tracks.begin(), tracks.end(), std::back_inserter( rc ),
                  [&]( const auto& t ) { return selector( t ); } );
    return rc;
  }

  struct TrackVeloSidePredicate {
    int m_sign;
    TrackVeloSidePredicate( int asign ) : m_sign( asign ) {}
    bool operator()( const LHCb::Track* track ) const {
      int  side         = 0;
      bool allhitsleft  = true;
      bool allhitsright = true;

      const std::vector<LHCb::LHCbID>& track_ids = track->lhcbIDs();
      for ( const auto& track_id : track_ids ) {
        if ( !track_id.isVP() ) continue;

        auto vp_id = track_id.vpID();
        // 0 should be right, 1 left side
        bool sidepos = vp_id.sidepos();

        allhitsleft  = allhitsleft && sidepos;
        allhitsright = allhitsright && !sidepos;
      }
      if ( allhitsleft ) side = +1;
      if ( allhitsright ) side = -1;
      return side * m_sign > 0;
    }
  };
} // namespace

//=============================================================================
// Execute
//=============================================================================

void TrackPV2HalfMonitor::operator()( LHCb::ODIN const& odin, LHCb::Track::Range const& alltracks,
                                      DetectorElement const& lhcb ) const {

  ulonglong evTimeGps = odin.gpsTime();

  long long int tzero = 1270064494071424ll; // there is an ll at the end, so that C++ knows this has to be a long long
  float         eventTimeGpsMinute = (float)( ( ( evTimeGps - tzero ) / 1000000. ) / 3600. );

  // get the input data
  std::vector<const LHCb::Track*> selectedtracks;
  std::copy_if( alltracks.begin(), alltracks.end(), std::back_inserter( selectedtracks ),
                []( const LHCb::Track* tr ) { return tr->hasVelo() && tr->chi2PerDoF() < 10; } );

  // split the track in right and left, to evalute PV by only right or left tracks
  auto lefttracks  = myselect( selectedtracks, TrackVeloSidePredicate( +1 ) );
  auto righttracks = myselect( selectedtracks, TrackVeloSidePredicate( -1 ) );
  if ( lefttracks.size() >= 2 && righttracks.size() >= 2 ) {
    std::vector<LHCb::RecVertex> leftoutvtxvec;
    std::vector<LHCb::RecVertex> rightoutvtxvec;
    m_toolpv->reconstructMultiPVFromTracks( righttracks, rightoutvtxvec, *lhcb.geometry() ).ignore();
    m_toolpv->reconstructMultiPVFromTracks( lefttracks, leftoutvtxvec, *lhcb.geometry() ).ignore();

    const LHCb::RecVertex* leftvertex    = nullptr;
    const LHCb::RecVertex* rightvertex   = nullptr;
    int                    n_goodleftPV  = 0;
    int                    n_goodrightPV = 0;

    for ( const LHCb::RecVertex& pv : leftoutvtxvec ) {
      ++( *m_numTracksPerPVLeft )[pv.tracks().size()];
      if ( pv.tracks().size() >= m_PV_trackmin ) {
        ++n_goodleftPV;
        leftvertex = &pv;
        ++( *m_xPosPVLeft )[pv.position().x()];
        ++( *m_yPosPVLeft )[pv.position().y()];
        ++( *m_zPosPVLeft )[pv.position().z()];
        ++( *m_PVChi2Left )[pv.chi2() / pv.nDoF()];
        std::vector<const LHCb::Track*> pvtracks       = myconvert( pv.tracks() );
        auto                            refittedvertex = m_vertexer->fit( pvtracks, *lhcb.geometry() );
        if ( refittedvertex ) ++( *m_refittedPVChi2Left )[refittedvertex->chi2() / refittedvertex->nDoF()];
      }
    }

    for ( const LHCb::RecVertex& pv : rightoutvtxvec ) {
      ++( *m_numTracksPerPVRight )[pv.tracks().size()];
      if ( pv.tracks().size() >= m_PV_trackmin ) {
        ++n_goodrightPV;
        rightvertex = &pv;
        ++( *m_xPosPVRight )[pv.position().x()];
        ++( *m_yPosPVRight )[pv.position().y()];
        ++( *m_zPosPVRight )[pv.position().z()];
        ++( *m_PVChi2Right )[pv.chi2() / pv.nDoF()];
        std::vector<const LHCb::Track*> pvtracks       = myconvert( pv.tracks() );
        auto                            refittedvertex = m_vertexer->fit( pvtracks, *lhcb.geometry() );
        if ( refittedvertex ) ++( *m_refittedPVChi2Right )[refittedvertex->chi2() / refittedvertex->nDoF()];
      }
    }

    if ( leftoutvtxvec.size() == 1 && rightoutvtxvec.size() == 1 && rightvertex && leftvertex ) {

      if ( msgLevel( MSG::DEBUG ) )
        debug() << "Found " << n_goodrightPV << " Right PV and " << n_goodleftPV << " Left PV" << endmsg;

      const Gaudi::XYZVector dx  = leftvertex->position() - rightvertex->position();
      const auto             cov = leftvertex->covMatrix() + rightvertex->covMatrix();
      ++( *m_deltaxPVLeftRight )[dx.x()];
      ++( *m_deltayPVLeftRight )[dx.y()];
      ++( *m_deltazPVLeftRight )[dx.z()];
      ++( *m_deltaxPVLeftRightPull )[dx.x() / std::sqrt( cov( 0, 0 ) )];
      ++( *m_deltayPVLeftRightPull )[dx.y() / std::sqrt( cov( 1, 1 ) )];
      ++( *m_deltazPVLeftRightPull )[dx.z() / std::sqrt( cov( 2, 2 ) )];
      ++( *m_eventTimeHisto )[eventTimeGpsMinute];

      if ( std::abs( dx.z() ) < m_limdpvz ) {
        double z = 0.5 * ( leftvertex->position().z() + rightvertex->position().z() );
        if ( std::abs( dx.x() ) < m_limdpvx ) ( *m_deltaxzPVLeftRight )[z] += dx.x();
        if ( std::abs( dx.y() ) < m_limdpvy ) ( *m_deltayzPVLeftRight )[z] += dx.y();
        if ( std::abs( dx.z() ) < m_limdpvz ) ( *m_deltazzPVLeftRight )[z] += dx.z();
      }
    }
  }
}
