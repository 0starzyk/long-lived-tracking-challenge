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
#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/ToolHandle.h"
#include "Kernel/HitPattern.h"
#include "LHCbAlgs/Transformer.h"
#include "TrackInterfaces/IPVOfflineTool.h"

#include <algorithm>

using RecVertices = LHCb::RecVertices;
using RecVertex   = LHCb::RecVertex;
using Output      = std::tuple<RecVertices, RecVertices, RecVertices>;

class PatPVLeftRightFinderMerger
    : public LHCb::Algorithm::MultiTransformer<Output( LHCb::Track::Range const&, DetectorElement const& ),
                                               LHCb::DetDesc::usesConditions<DetectorElement>> {
public:
  /** Standard construtor */
  PatPVLeftRightFinderMerger( const std::string& name, ISvcLocator* pSvcLoc )
      : MultiTransformer{name,
                         pSvcLoc,
                         {KeyValue{"InputTracks", LHCb::TrackLocation::Default},
                          KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}},
                         {KeyValue{"OutputLeftVertices", ""}, KeyValue{"OutputRightVertices", ""},
                          KeyValue{"OutputVertices", ""}}} {}

  /** Algorithm execute */
  Output operator()( LHCb::Track::Range const&, DetectorElement const& lhcbtop ) const override;

private:
  Gaudi::Property<double>    m_maxDeltaZForMerge{this, "MaxDeltaZForMerge", 1 * Gaudi::Units::mm};
  Gaudi::Property<double>    m_maxDeltaZChi2ForMerge{this, "MaxDeltaZChi2ForMerge", 9};
  Gaudi::Property<bool>      m_addUnmergedVertices{this, "AddUnmergedVertices", false};
  ToolHandle<IPVOfflineTool> m_toolpv{this, "PVOfflineTool", "PVOfflineTool"};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( PatPVLeftRightFinderMerger )

//=============================================================================
// Structure
//=============================================================================

namespace {

  RecVertices convert( const std::vector<LHCb::RecVertex>& vertices ) {
    RecVertices rc;
    rc.reserve( vertices.size() );
    for ( const auto& v : vertices ) rc.insert( new RecVertex{v} );
    return rc;
  }

} // namespace

//=============================================================================
// Execute
//=============================================================================

Output PatPVLeftRightFinderMerger::operator()( LHCb::Track::Range const& alltracks,
                                               const DetectorElement&    lhcbtop ) const {

  // divide the tracks in to 'valid' left and right tracks. we omit 'mixed' tracks.
  std::vector<const LHCb::Track*> lefttracks, righttracks;
  for ( const auto& tr : alltracks ) {
    LHCb::HitPattern hitpattern{tr->lhcbIDs()};
    if ( hitpattern.numVeloA() > hitpattern.numVeloC() )
      lefttracks.push_back( &( *tr ) );
    else if ( hitpattern.numVeloC() > hitpattern.numVeloA() )
      righttracks.push_back( &( *tr ) );
  }

  debug() << "Number of input tracks: " << alltracks.size() << ", " << lefttracks.size() << ", " << righttracks.size()
          << endmsg;

  // construct left and right vertices
  std::vector<LHCb::RecVertex> leftvertices, rightvertices;
  if ( lefttracks.size() >= 2 )
    m_toolpv->reconstructMultiPVFromTracks( lefttracks, leftvertices, *lhcbtop.geometry() ).ignore();
  if ( righttracks.size() >= 2 )
    m_toolpv->reconstructMultiPVFromTracks( righttracks, rightvertices, *lhcbtop.geometry() ).ignore();

  // info() << "Number of left, right vertices: " << leftvertices.size() << ", " << rightvertices.size() << endmsg ;
  // Now match them
  // 1. combine (keeping track of side) and sort
  std::vector<std::tuple<int, const LHCb::RecVertex*>> allvertices;
  for ( const auto& v : leftvertices ) allvertices.emplace_back( 0, &v );
  for ( const auto& v : rightvertices ) allvertices.emplace_back( 1, &v );
  std::sort( allvertices.begin(), allvertices.end(), []( const auto& lhs, const auto& rhs ) {
    return std::get<1>( lhs )->position().z() < std::get<1>( rhs )->position().z();
  } );

  // 2. cluster
  LHCb::RecVertices mergedvertices;
  const auto        N = allvertices.size();
  if ( N > 0 ) {
    for ( size_t i = 0; i < N - 1; ++i ) {

      const auto& [iside, ivtx] = allvertices[i];
      const auto& [jside, jvtx] = allvertices[i + 1];
      bool merged               = false;
      if ( iside != jside ) {
        // consider if we can merge these vertices
        const auto dz   = jvtx->position().z() - ivtx->position().z();
        const auto zcov = jvtx->covMatrix()( 2, 2 ) + ivtx->covMatrix()( 2, 2 );
        // info() << "dz, dzchi2: " << dz << " " << dz*dz/zcov << endmsg ;
        if ( dz < m_maxDeltaZForMerge && dz * dz < m_maxDeltaZChi2ForMerge * zcov ) {
          // check that there isn't a better candidate to merge the next vertex with
          bool veto = false;
          if ( i + 2 < N ) {
            const auto& [kside, kvtx] = allvertices[i + 2];
            veto                      = ( kside != jside ) && ( kvtx->position().z() - jvtx->position().z() < dz );
          }
          if ( !veto ) {
            // merge these two vertices
            merged         = true;
            auto mergedvtx = new LHCb::RecVertex{};
            // we should perhaps make this a proper weighted average.
            Gaudi::XYZPoint position{0.5 * ( Gaudi::XYZVector{ivtx->position()} + Gaudi::XYZVector{jvtx->position()} )};
            mergedvtx->setPosition( position );
            mergedvtx->setCovMatrix( 0.5 * ( ivtx->covMatrix() + jvtx->covMatrix() ) );
            mergedvtx->setChi2AndDoF( ivtx->chi2() + jvtx->chi2(), ivtx->nDoF() + jvtx->nDoF() );
            // add all the tracks
            for ( const auto& itrk : ivtx->tracksWithWeights() ) mergedvtx->addToTracks( itrk.first, itrk.second );
            for ( const auto& jtrk : jvtx->tracksWithWeights() ) mergedvtx->addToTracks( jtrk.first, jtrk.second );
            mergedvertices.insert( mergedvtx );
            // make sure to increase the pointer
            ++i;
          }
        }
      }
      if ( !merged && m_addUnmergedVertices ) mergedvertices.insert( new RecVertex{*ivtx} );
    }
  }
  debug() << "Number of merged vertices: " << mergedvertices.size() << endmsg;

  // Also create output containers for the left and right vertices, just in case anyone is interested.
  RecVertices keyedleftvertices{convert( leftvertices )};
  RecVertices keyedrightvertices{convert( rightvertices )};
  return Output{std::move( keyedleftvertices ), std::move( keyedrightvertices ), std::move( mergedvertices )};
}
