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
#include "Event/MCParticle.h"
#include "Event/MuonCoord.h"
#include "Event/Track.h"
#include "Event/UTCluster.h"
#include "Event/VPCluster.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "Linker/LinkedTo.h"
#include "MCInterfaces/ILHCbIDsToMCParticles.h"

/** @class LHCbIDsToMCParticles LHCbIDsToMCParticles.h
 *
 *  Link ids to MCParticles
 *
 *  @author M.Needham
 *  @date   31/04/2006
 */

class LHCbIDsToMCParticles : public extends<GaudiTool, ILHCbIDsToMCParticles, IIncidentListener> {

public:
  /// constructer
  using extends::extends;

  /** initialize */
  StatusCode initialize() override;

  /**
    Trivial link from list of IDs to all particles contributing
    @param start  iterator to first id
    @param stop   iterator to last id
    @param output vector by reference
    @return StatusCode
  */
  StatusCode link( LHCbIDs::const_iterator start, LHCbIDs::const_iterator stop, LinkMap& output ) const override;

  /**
    Trivial link from list of ALL ids in track to particles contributing
    @param aTrack track
    @param output vector by reference
    @return StatusCode
  */
  StatusCode link( const LHCb::Track& aTrack, LinkMap& output ) const override;

  /**
    Trivial link from single id to particles contributing
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
  void linkToDetTruth( const ID& id, const LINKER& aLinker, LinkMap& output ) const;

  StatusCode linkMuon( const LHCb::LHCbID& id, LinkMap& output ) const;
  StatusCode linkVP( const LHCb::LHCbID& id, LinkMap& output ) const;
  StatusCode linkUT( const LHCb::LHCbID& id, LinkMap& output ) const;
  StatusCode linkFT( const LHCb::LHCbID& id, LinkMap& output ) const;

  mutable LHCb::LinksByKey const* m_muonLinks = nullptr;
  mutable LHCb::LinksByKey const* m_vpLinks   = nullptr;
  mutable LHCb::LinksByKey const* m_utLinks   = nullptr;
  mutable LHCb::LinksByKey const* m_ftLinks   = nullptr;
};

#include "Event/MCParticle.h"
template <typename ID, typename LINKER>
void LHCbIDsToMCParticles::linkToDetTruth( const ID& id, const LINKER& aLinker, LinkMap& output ) const {

  auto r = aLinker.range( id );
  if ( !r.empty() ) {
    for ( const auto& p : r ) output[&p] += 1;
  } else {
    output[nullptr] += 1;
  }
}

DECLARE_COMPONENT( LHCbIDsToMCParticles )
using namespace LHCb;

StatusCode LHCbIDsToMCParticles::initialize() {
  return extends::initialize().andThen( [&] {
    // warning() << "This code is not thread-safe. Do NOT use in multi-threaded jobs" << endmsg;
    incSvc()->addListener( this, IncidentType::BeginEvent );
  } );
}

StatusCode LHCbIDsToMCParticles::link( LHCbIDs::const_iterator start, LHCbIDs::const_iterator stop,
                                       LinkMap& output ) const {
  return std::accumulate( start, stop, StatusCode{StatusCode::SUCCESS}, [&]( StatusCode sc, const auto& id ) {
    return sc.andThen( [&] { return link( id, output ); } );
  } );
}

StatusCode LHCbIDsToMCParticles::link( const Track& aTrack, LinkMap& output ) const {
  auto const& ids = aTrack.lhcbIDs();
  return link( ids.begin(), ids.end(), output );
}

StatusCode LHCbIDsToMCParticles::link( LHCbID id, LinkMap& output ) const {
  switch ( id.detectorType() ) {
  case LHCbID::channelIDtype::Muon:
    return linkMuon( id, output );
  case LHCbID::channelIDtype::VP:
    return linkVP( id, output );
  case LHCbID::channelIDtype::UT:
    return linkUT( id, output );
  case LHCbID::channelIDtype::FT:
    return linkFT( id, output );
  default:
    return Warning( "Unknown type !", StatusCode::SUCCESS, 10 );
  }
}

void LHCbIDsToMCParticles::handle( const Incident& incident ) {
  if ( IncidentType::BeginEvent == incident.type() ) {
    m_muonLinks = nullptr;
    m_vpLinks   = nullptr;
    m_utLinks   = nullptr;
    m_ftLinks   = nullptr;
  }
}

StatusCode LHCbIDsToMCParticles::linkMuon( const LHCbID& lhcbid, LinkMap& output ) const {
  if ( !m_muonLinks ) {
    m_muonLinks =
        SmartDataPtr<LHCb::LinksByKey>( evtSvc(), LHCb::LinksByKey::linkerName( LHCb::MuonCoordLocation::MuonCoords ) );
    if ( !m_muonLinks ) { return Error( "no MuonLinker", StatusCode::FAILURE, 10 ); }
  }
  linkToDetTruth( lhcbid.muonID(), LinkedTo<LHCb::MCParticle>( m_muonLinks ), output );
  return StatusCode::SUCCESS;
}

StatusCode LHCbIDsToMCParticles::linkVP( const LHCbID& lhcbid, LinkMap& output ) const {
  if ( !m_vpLinks ) {
    m_vpLinks =
        SmartDataPtr<LHCb::LinksByKey>( evtSvc(), LHCb::LinksByKey::linkerName( LHCb::VPClusterLocation::Default ) );
    if ( !m_vpLinks ) { return Error( "no VPLinker", StatusCode::FAILURE, 10 ); }
  }
  linkToDetTruth( lhcbid.vpID(), LinkedTo<LHCb::MCParticle>{m_vpLinks}, output );
  return StatusCode::SUCCESS;
}

StatusCode LHCbIDsToMCParticles::linkUT( const LHCbID& lhcbid, LinkMap& output ) const {
  if ( !m_utLinks ) {
    m_utLinks =
        SmartDataPtr<LHCb::LinksByKey>( evtSvc(), LHCb::LinksByKey::linkerName( LHCb::UTClusterLocation::UTClusters ) );
    if ( !m_utLinks ) { return Error( "no UTLinker", StatusCode::FAILURE, 10 ); }
  }
  linkToDetTruth( lhcbid.utID(), LinkedTo<LHCb::MCParticle>{m_utLinks}, output );
  return StatusCode::SUCCESS;
}

StatusCode LHCbIDsToMCParticles::linkFT( const LHCbID& lhcbid, LinkMap& output ) const {
  if ( !m_ftLinks ) {
    m_ftLinks = SmartDataPtr<LHCb::LinksByKey>( evtSvc(),
                                                LHCb::LinksByKey::linkerName( LHCb::FTLiteClusterLocation::Default ) );
    if ( !m_ftLinks ) { return Error( "no FTLinker", StatusCode::FAILURE, 10 ); }
  }
  linkToDetTruth( lhcbid.ftID(), LinkedTo<LHCb::MCParticle>{m_ftLinks}, output );
  return StatusCode::SUCCESS;
}
