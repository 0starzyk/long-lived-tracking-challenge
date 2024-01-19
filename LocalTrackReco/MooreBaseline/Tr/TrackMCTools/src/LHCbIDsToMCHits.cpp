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
#include "Event/FTLiteCluster.h"
#include "Event/MCHit.h"
#include "Event/MCProperty.h"
#include "Event/MuonCoord.h"
#include "Event/Track.h"
#include "Event/UTCluster.h"
#include "Event/VPCluster.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "Linker/LinkedTo.h"
#include "MCInterfaces/ILHCbIDsToMCHits.h"

using namespace LHCb;

/** @class LHCbIDsToMCHits LHCbIDsToMCHits.h
 *
 *  Link ids to MCHits
 *
 *  WARNING: this code is NOT threadsafe, and should only be used in single-threaded jobs...
 *
 *  @author M.Needham
 *  @date   31/04/2006
 */

class LHCbIDsToMCHits : public extends<GaudiTool, ILHCbIDsToMCHits, IIncidentListener> {

public:
  /// constructer
  using extends::extends;

  /** initialize */
  StatusCode initialize() override;

  /**
     Trivial link from list of IDs to all MCHits contributing
     @param start  iterator to first id
     @param stop   iterator to last id
     @param output vector by reference
     @return StatusCode
  */
  StatusCode link( LHCbIDs::const_iterator start, LHCbIDs::const_iterator stop, LinkMap& output ) const override;

  /**
     Trivial link from list of ALL ids in track to MCHits contributing
     @param aTrack track
     @param output vector by reference
     @return StatusCode
  */
  StatusCode link( const LHCb::Track& aTrack, LinkMap& output ) const override;

  /**
     Trivial link from single id to MCHits contributing
     @param id
     @param output vector by reference
     @return StatusCode
  */
  StatusCode link( LHCb::LHCbID id, LinkMap& output ) const override;

  /** Implement the handle method for the Incident service.
   *  This is used to inform the tool of software incidents.
   *
   *  @param incident The incident identifier
   */
  void handle( const Incident& incident ) override;

private:
  template <typename ID, typename LINKER>
  void linkToDetTruth( const ID& id, LINKER&& aLinker, LinkMap& output ) const;
  void linkUT( const LHCb::LHCbID& id, LinkMap& output ) const;
  void linkVP( const LHCb::LHCbID& id, LinkMap& output ) const;
  void linkMuon( const LHCb::LHCbID& id, LinkMap& output ) const;
  void linkFT( const LHCb::LHCbID& id, LinkMap& output ) const;

  mutable LHCb::LinksByKey const* m_utLinks   = nullptr;
  mutable LHCb::LinksByKey const* m_vpLinks   = nullptr;
  mutable LHCb::LinksByKey const* m_muonLinks = nullptr;
  mutable LHCb::LinksByKey const* m_ftLinks   = nullptr;

  static constexpr auto m_endString = "2MCHits";
};

DECLARE_COMPONENT( LHCbIDsToMCHits )

#include "Event/MCHit.h"

/// Link LHCbID to MCHits in detector in question
template <typename ID, typename LINKER>
void LHCbIDsToMCHits::linkToDetTruth( const ID& id, LINKER&& aLinker, LinkMap& output ) const {
  auto r = aLinker.range( id );
  if ( !r.empty() ) {
    for ( const auto& aHit : r ) output[&aHit] += 1;
  } else {
    output[nullptr] += 1;
  }
}

StatusCode LHCbIDsToMCHits::initialize() {
  return extends::initialize().andThen( [&] { incSvc()->addListener( this, IncidentType::BeginEvent ); } );
}

StatusCode LHCbIDsToMCHits::link( LHCbIDs::const_iterator start, LHCbIDs::const_iterator stop, LinkMap& output ) const {

  for ( auto iter = start; iter != stop; ++iter ) {
    link( *iter, output ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  } // iter

  return StatusCode::SUCCESS;
}

StatusCode LHCbIDsToMCHits::link( const Track& aTrack, LinkMap& output ) const {
  auto const& ids = aTrack.lhcbIDs();
  return link( ids.begin(), ids.end(), output );
}

StatusCode LHCbIDsToMCHits::link( LHCbID id, LinkMap& output ) const {

  // switch statement from hell
  switch ( id.detectorType() ) {
  case LHCbID::channelIDtype::UT:
    linkUT( id, output );
    break;
  case LHCbID::channelIDtype::VP:
    linkVP( id, output );
    break;
  case LHCbID::channelIDtype::Muon:
    linkMuon( id, output );
    break;
  case LHCbID::channelIDtype::FT:
    linkFT( id, output );
    break;

  default:
    Warning( "Unknown type !", StatusCode::SUCCESS, 10 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    break;
  }

  return StatusCode::SUCCESS;
}

void LHCbIDsToMCHits::handle( const Incident& incident ) {
  if ( IncidentType::BeginEvent == incident.type() ) {
    m_utLinks   = nullptr;
    m_vpLinks   = nullptr;
    m_muonLinks = nullptr;
    m_ftLinks   = nullptr;
    ;
  }
}

void LHCbIDsToMCHits::linkUT( const LHCbID& lhcbid, LinkMap& output ) const {
  if ( !m_utLinks ) {
    m_utLinks = SmartDataPtr<LHCb::LinksByKey>{
        evtSvc(), LHCb::LinksByKey::linkerName( LHCb::UTClusterLocation::UTClusters + m_endString )};
    if ( !m_utLinks ) { throw GaudiException( "no UTLinker", "LHCbIDsToMCHits", StatusCode::FAILURE ); }
  }
  linkToDetTruth( lhcbid.utID(), LinkedTo<LHCb::MCHit>{m_utLinks}, output );
}

void LHCbIDsToMCHits::linkVP( const LHCbID& lhcbid, LinkMap& output ) const {

  if ( !m_vpLinks ) {
    m_vpLinks = SmartDataPtr<LHCb::LinksByKey>{
        evtSvc(), LHCb::LinksByKey::linkerName( LHCb::VPClusterLocation::Default + m_endString )};
    if ( !m_vpLinks ) { throw GaudiException( "no vPLinker", "LHCbIDsToMCHits", StatusCode::FAILURE ); }
  }
  linkToDetTruth( lhcbid.vpID(), LinkedTo<LHCb::MCHit>{m_vpLinks}, output );
}

void LHCbIDsToMCHits::linkFT( const LHCbID& lhcbid, LinkMap& output ) const {
  if ( !m_ftLinks ) {
    m_ftLinks = SmartDataPtr<LHCb::LinksByKey>{
        evtSvc(), LHCb::LinksByKey::linkerName( LHCb::FTLiteClusterLocation::Default + m_endString )};
    if ( !m_ftLinks ) { throw GaudiException( "no FTLinker", "LHCbIDsToMCHits", StatusCode::FAILURE ); }
  }
  linkToDetTruth( lhcbid.ftID(), LinkedTo<LHCb::MCHit>{m_ftLinks}, output );
}

void LHCbIDsToMCHits::linkMuon( const LHCbID& lhcbid, LinkMap& output ) const {

  if ( !m_muonLinks ) {
    m_muonLinks = SmartDataPtr<LHCb::LinksByKey>{
        evtSvc(), LHCb::LinksByKey::linkerName( LHCb::MuonCoordLocation::MuonCoords + m_endString )};
    if ( !m_muonLinks ) { throw GaudiException( "no MuonLinker", "LHCbIDsToMCHits", StatusCode::FAILURE ); }
  }
  linkToDetTruth( lhcbid.muonID(), LinkedTo<LHCb::MCHit>{m_muonLinks}, output );
}
