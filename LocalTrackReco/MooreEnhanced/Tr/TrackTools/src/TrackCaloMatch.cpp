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

/** @class TrackCaloMatch TrackCaloMatch.h
 *
 * Implementation of TrackCaloMatch tool
 * see interface header for description
 *
 *  @author M.Needham
 *  @date   30/12/2005
 */

#include "CaloFutureUtils/CaloFuture2Track.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "Relations/IRelation.h"
#include "TrackInterfaces/ITrackCaloMatch.h"
#include <string>

using namespace LHCb;

class TrackCaloMatch : public extends<GaudiTool, ITrackCaloMatch, IIncidentListener> {

public:
  /// constructor
  using extends::extends;

  StatusCode initialize() override;

  /// the method
  double energy( const LHCb::Track& aTrack ) const override;

  /** Implement the handle method for the Incident service.
   *  This is used to inform the tool of software incidents.
   *
   *  @param incident The incident identifier
   */
  void handle( const Incident& incident ) override;

private:
  void initEvent() const;

  typedef IRelation<LHCb::Track, float> Table;
  double                                energy( const LHCb::Track& aTrack, const TrackCaloMatch::Table* table ) const;

  mutable Table* m_ecalE = nullptr;
  mutable Table* m_hcalE = nullptr;
  mutable Table* m_psE   = nullptr;

  std::string m_ecalLocation;
  std::string m_hcalLocation;
  std::string m_prsLocation;

  Gaudi::Property<double> m_alpha{this, "alpha", 8.};
  Gaudi::Property<double> m_beta{this, "beta", 1.};
  Gaudi::Property<double> m_gamma{this, "gamma", 1.};

  mutable bool m_configured = false;
};

DECLARE_COMPONENT( TrackCaloMatch )

StatusCode TrackCaloMatch::initialize() {

  StatusCode sc = GaudiTool::initialize();
  if ( sc.isFailure() ) { return Error( "Failed to initialize", sc ); }

  m_ecalLocation = rootInTES() + CaloFutureIdLocation::EcalE;
  m_hcalLocation = rootInTES() + CaloFutureIdLocation::HcalE;
  m_prsLocation  = rootInTES() + CaloFutureIdLocation::PrsE;

  incSvc()->addListener( this, IncidentType::BeginEvent );

  return StatusCode::SUCCESS;
}

double TrackCaloMatch::energy( const Track& aTrack ) const {

  // get the input - seeds
  if ( !m_configured ) initEvent();

  double eEcal = energy( aTrack, m_ecalE );
  double eHcal = energy( aTrack, m_hcalE );
  double ePrs  = energy( aTrack, m_psE );

  // known bug - sometimes ps gives -ve energy
  if ( ePrs < 0 ) {
    eEcal = 0;
    ePrs  = 0;
  }

  return ( m_alpha * ePrs ) + ( m_beta * eEcal ) + ( m_gamma * eHcal );
}

double TrackCaloMatch::energy( const Track& aTrack, const TrackCaloMatch::Table* table ) const {

  Table::Range aRange = table->relations( &aTrack );
  return ( !aRange.empty() ? aRange.front().to() : 0 );
}

void TrackCaloMatch::handle( const Incident& incident ) {
  if ( IncidentType::BeginEvent == incident.type() ) { m_configured = false; }
}

void TrackCaloMatch::initEvent() const {

  m_configured = true;

  using namespace LHCb::CaloFutureIdLocation;
  m_ecalE = get<Table>( m_ecalLocation );
  m_hcalE = get<Table>( m_hcalLocation );
  m_psE   = get<Table>( m_prsLocation );
}
