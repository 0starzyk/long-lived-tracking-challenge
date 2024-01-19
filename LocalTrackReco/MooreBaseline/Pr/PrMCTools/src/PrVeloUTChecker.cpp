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
#include "Event/MCParticle.h"
#include "Event/Track.h"
#include "Event/UTCluster.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Linker/LinkedFrom.h"
#include "Linker/LinkerTool.h"
#include "MCInterfaces/IMCReconstructible.h"

/** @class PrVeloUTChecker PrVeloUTChecker.h "PrMCTools/PrVeloUTChecker"
 *
 * Class for VeloTT studies
 *  @author A. Di Canto
 *  @date   2012-05-18
 */

class PrVeloUTChecker : public GaudiHistoAlg {
public:
  /** Standard construtor */
  PrVeloUTChecker( const std::string& name, ISvcLocator* pSvcLocator );

  /** Algorithm execute */
  StatusCode execute() override;

  /** Algorithm initialize */
  StatusCode initialize() override;

  /** Algorithm finalize */
  StatusCode finalize() override;

private:
  void printStatistics();
  void getTTtruth( const LHCb::MCParticle*, int& hits, int& layers );
  bool bAncestor( const LHCb::MCParticle* ) const;
  void fillEffHistos( double var, std::string title, double min, double max, bool found,
                      std::vector<bool> const& flags );
  void fillGhostHistos( double var, std::string title, double min, double max, bool found,
                        std::vector<bool> const& flags );

  /** Data memebers */
  std::string                    m_ForwardContainer; ///< Input Forward tracks container location
  std::string                    m_VeloContainer;    ///< Input Velo tracks container location
  std::string                    m_VeloTTContainer;  ///< Input VeloTT tracks container location
  ToolHandle<IMCReconstructible> m_selector{this, "MCReconstrictible",
                                            "MCReconstructible/Selector"}; ///< Pointer to selector

  class MyCounter {
  public:
    std::vector<std::string> name;
    std::vector<int>         num;
    std::vector<int>         den;

    void addCategory( std::string );
    void count( bool, std::vector<bool> const& );
  };

  MyCounter m_eff_counter;
  MyCounter m_ghosts_counter;
};

DECLARE_COMPONENT( PrVeloUTChecker )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrVeloUTChecker::PrVeloUTChecker( const std::string& name, ISvcLocator* pSvcLocator )
    : GaudiHistoAlg( name, pSvcLocator ) {
  declareProperty( "VeloTTContainer", m_VeloTTContainer = LHCb::TrackLocation::VeloTT );
  declareProperty( "VeloContainer", m_VeloContainer = LHCb::TrackLocation::Velo );
  declareProperty( "ForwardContainer", m_ForwardContainer = LHCb::TrackLocation::Forward );
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrVeloUTChecker::initialize() {
  // Mandatory initialization of GaudiAlgorithm
  StatusCode sc = GaudiHistoAlg::initialize();
  if ( sc.isFailure() ) { return sc; }

  m_eff_counter.addCategory( "reco'ble as long                     " );
  m_eff_counter.addCategory( "reco'ble as long + reco'ed in Velo   " );
  m_eff_counter.addCategory( "reco'ble as long + reco'ed by Forward" );
  m_eff_counter.addCategory( "reco'ble as long + in TT acceptance  " );
  m_eff_counter.addCategory( "loose HLT                            " );
  m_eff_counter.addCategory( "loose HLT + from B                   " );
  m_eff_counter.addCategory( "loose HLT + in TT acceptance         " );
  m_eff_counter.addCategory( "loose HLT + in TT acceptance + from B" );
  m_eff_counter.addCategory( "tight HLT                            " );
  m_eff_counter.addCategory( "tight HLT + from B                   " );
  m_eff_counter.addCategory( "tight HLT + in TT acceptance         " );
  m_eff_counter.addCategory( "tight HLT + in TT acceptance + from B" );

  m_ghosts_counter.addCategory( "all         " );
  m_ghosts_counter.addCategory( " + good Velo" );
  m_ghosts_counter.addCategory( " + loose HLT" );
  m_ghosts_counter.addCategory( " + tight HLT" );

  return sc;
}

//=============================================================================
// Execute
//=============================================================================
StatusCode PrVeloUTChecker::execute() {

  // check input data
  if ( !exist<LHCb::Tracks>( m_VeloTTContainer ) )
    return Warning( m_VeloTTContainer + " not found", StatusCode::SUCCESS, 0 );
  if ( !exist<LHCb::Tracks>( m_VeloContainer ) )
    return Warning( m_VeloContainer + " not found", StatusCode::SUCCESS, 0 );
  if ( !exist<LHCb::Tracks>( m_ForwardContainer ) )
    return Warning( m_ForwardContainer + " not found", StatusCode::SUCCESS, 0 );

  using AsctTool = LinkerTool<LHCb::Track, LHCb::MCParticle>;
  // get the association tables
  auto       associator  = AsctTool( evtSvc(), m_VeloTTContainer );
  const auto directTable = associator.direct();
  if ( !directTable ) return Error( "Failed to find direct table for VeloTT tracks", StatusCode::FAILURE );
  const auto inverseTable = associator.inverse();
  if ( !inverseTable ) return Error( "Failed to find inverse table for VeloTT tracks", StatusCode::FAILURE );

  auto       associator_velo  = AsctTool( evtSvc(), m_VeloContainer );
  const auto directTable_velo = associator_velo.direct();
  if ( !directTable_velo ) return Error( "Failed to find direct table for Velo tracks", StatusCode::FAILURE );
  const auto inverseTable_velo = associator_velo.inverse();
  if ( !inverseTable_velo ) return Error( "Failed to find inverse table for Velo tracks", StatusCode::FAILURE );

  auto       associator_forward   = AsctTool( evtSvc(), m_ForwardContainer );
  const auto inverseTable_forward = associator_forward.inverse();
  if ( !inverseTable_forward ) return Error( "Failed to find inverse table for Forward tracks", StatusCode::FAILURE );

  // compute performances
  // computeEfficiency;
  // loop over MC particles
  //=============================================================================
  // Loop over MC particles and look for reconstructed tracks
  //=============================================================================
  const LHCb::MCParticles* particles = get<LHCb::MCParticles>( LHCb::MCParticleLocation::Default );
  for ( auto const& ip : *particles ) {
    if ( !m_selector->isReconstructibleAs( IMCReconstructible::ChargedLong, ip ) ) continue;
    if ( std::abs( ip->particleID().pid() ) == 11 ) continue;

    std::vector<bool> flags;
    flags.push_back( true ); // reco'ble as long

    flags.push_back( !inverseTable_velo.relations( ip ).empty() ); // reco'ed in Velo

    bool reconstructed_forward = !inverseTable_forward.relations( ip ).empty();
    flags.push_back( reconstructed_forward ); // reco'ed by Forward

    int nTTlayers, nTThits;
    getTTtruth( ip, nTThits, nTTlayers );
    flags.push_back( reconstructed_forward && ( nTTlayers > 2 ) ); // reco'ed by Forward and in TT acceptance

    bool loose = reconstructed_forward && ip->momentum().P() > 3000. && ip->momentum().Pt() > 500.;
    bool tight = reconstructed_forward && ip->momentum().P() > 3000. && ip->momentum().Pt() > 1250.;
    bool fromB = bAncestor( ip );

    flags.push_back( loose );                               // loose HLT
    flags.push_back( loose && fromB );                      // loose HLT from B
    flags.push_back( loose && ( nTTlayers > 2 ) );          // loose HLT and in TT acceptance
    flags.push_back( loose && fromB && ( nTTlayers > 2 ) ); // loose HLT from B and in TT acceptance

    flags.push_back( tight );                               // tight HLT
    flags.push_back( tight && fromB );                      // tight HLT from B
    flags.push_back( tight && ( nTTlayers > 2 ) );          // tight HLT and in TT acceptance
    flags.push_back( tight && fromB && ( nTTlayers > 2 ) ); // tight HLT from B and in TT acceptance

    // count how many are reconstructed
    bool reconstructed = !inverseTable.relations( ip ).empty();

    if ( flags.size() != m_eff_counter.name.size() )
      error() << "Number of flags (" << flags.size() << ") does not coincide with number of categories ("
              << m_eff_counter.name.size() << ") in the efficiency counter" << endmsg;

    m_eff_counter.count( reconstructed, flags );

    fillEffHistos( ip->momentum().P() / 1e3, "p [GeV/c]", 0., 50., reconstructed, flags );
    fillEffHistos( ip->momentum().Pt() / 1e3, "p_{T} [GeV/c]", 0., 3., reconstructed, flags );
    fillEffHistos( ip->momentum().Eta(), "#eta", 1., 6., reconstructed, flags );
    fillEffHistos( ip->momentum().Phi(), "#phi [rad]", -4., 4., reconstructed, flags );
    fillEffHistos( nTTlayers, "MC TT layers", 0., 5., reconstructed, flags );
  }

  // computeGhostRate
  // Loop over tracks
  //=============================================================================
  // Loop over tracks and look for matched MC particles
  //=============================================================================
  for ( auto const& it : *get<LHCb::Tracks>( m_VeloTTContainer ) ) {

    LHCb::Track* veloTr     = *( ( it->ancestors() ).begin() );
    bool         velo_ghost = directTable_velo.relations( veloTr ).empty();
    bool         loose      = !velo_ghost && it->p() > 3000. && it->pt() > 500.;
    bool         tight      = !velo_ghost && it->p() > 3000. && it->pt() > 1250.;
    // all tracks,  good velo, loose HLT, tight HLT
    const std::vector<bool> flags = {true, !velo_ghost, loose, tight};

    // count how many are ghosts
    bool ghost = directTable.relations( it ).empty();

    if ( flags.size() != m_ghosts_counter.name.size() )
      error() << "Number of flags (" << flags.size() << ") does not coincide with number of categories ("
              << m_eff_counter.name.size() << ") in the ghosts counter" << endmsg;

    m_ghosts_counter.count( ghost, flags );

    fillGhostHistos( it->p() / 1e3, "p [GeV/c]", 0., 50., ghost, flags );
    fillGhostHistos( it->pt() / 1e3, "p_{T} [GeV/c]", 0., 3., ghost, flags );
    fillGhostHistos( it->pseudoRapidity(), "#eta", 1., 6., ghost, flags );
    fillGhostHistos( it->phi(), "#phi [rad]", -4., 4., ghost, flags );
    fillGhostHistos( it->chi2PerDoF(), "#chi^{2}/ndf", 0., 5., ghost, flags );
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
// Finalization
//=============================================================================
StatusCode PrVeloUTChecker::finalize() {
  printStatistics();
  return GaudiHistoAlg::finalize();
}

//=============================================================================
// Get MC TT hist for particle p
//=============================================================================
void PrVeloUTChecker::getTTtruth( const LHCb::MCParticle* p, int& nhits, int& nlayers ) {
  const LHCb::LinksByKey* tt_links =
      get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( LHCb::UTClusterLocation::UTClusters ) );

  nhits = 0;
  std::vector<int> fired;
  for ( const auto& TTCluster : LinkedFrom<LHCb::UTCluster>{tt_links}.range( p ) ) {
    if ( !TTCluster.isUT() ) continue;
    nhits++;
    int layer = TTCluster.layer() + 10 * TTCluster.station();
    if ( std::find( fired.begin(), fired.end(), layer ) == fired.end() ) fired.push_back( layer );
  }
  nlayers = fired.size();
}

//=============================================================================
// Look if particle mcPart comes from a b-hadron decay
//=============================================================================
bool PrVeloUTChecker::bAncestor( const LHCb::MCParticle* mcPart ) const {
  bool                    fromB  = false;
  const LHCb::MCParticle* mother = mcPart->mother();
  while ( mother != 0 && fromB == false ) {
    fromB  = mother->particleID().hasBottom() && ( mother->particleID().isMeson() || mother->particleID().isBaryon() );
    mother = mother->mother();
  }
  return fromB;
}

//=============================================================================
// Print final summary
//=============================================================================
void PrVeloUTChecker::printStatistics() {
  always() << "**** Velo TT efficiencies:" << endmsg;
  for ( unsigned int kk = 0; m_eff_counter.name.size() > kk; ++kk ) {
    if ( 0 == m_eff_counter.den[kk] ) {
      always() << "  " << m_eff_counter.name[kk] << " -- no particles found" << endmsg;
      continue;
    }
    double eff = double( m_eff_counter.num[kk] ) / double( m_eff_counter.den[kk] );
    double err = sqrt( eff * ( 1. - eff ) / double( m_eff_counter.den[kk] ) );
    always() << "  " << m_eff_counter.name[kk] << format( " (%4.2f +/- %.2f)%%", 100. * eff, 100. * err ) << endmsg;
  }

  always() << "**** Velo TT ghost rates:" << endmsg;
  for ( unsigned int kk = 0; m_ghosts_counter.name.size() > kk; ++kk ) {
    if ( 0 == m_ghosts_counter.den[kk] ) {
      always() << "  " << m_ghosts_counter.name[kk] << " -- no tracks found" << endmsg;
      continue;
    }
    double eff = double( m_ghosts_counter.num[kk] ) / double( m_ghosts_counter.den[kk] );
    double err = sqrt( eff * ( 1. - eff ) / double( m_ghosts_counter.den[kk] ) );
    always() << "  " << m_ghosts_counter.name[kk] << format( " (%4.2f +/- %.2f)%%", 100. * eff, 100. * err ) << endmsg;
  }
}

//=============================================================================
// Filling histos
//=============================================================================
void PrVeloUTChecker::fillEffHistos( double var, std::string title, double min, double max, bool found,
                                     std::vector<bool> const& flags ) {
  for ( unsigned int kk = 0; flags.size() > kk; ++kk ) {
    if ( !flags[kk] ) continue;
    std::string histoname = "MC particle " + title + " - " + m_eff_counter.name[kk];
    plot( var, histoname, min, max );
    if ( !found ) continue;
    plot( var, "reco'ed " + histoname, min, max );
  }
}

void PrVeloUTChecker::fillGhostHistos( double var, std::string title, double min, double max, bool found,
                                       std::vector<bool> const& flags ) {
  for ( unsigned int kk = 0; flags.size() > kk; ++kk ) {
    if ( !flags[kk] ) continue;
    std::string histoname = "track " + title + " - " + m_ghosts_counter.name[kk];
    plot( var, histoname, min, max );
    if ( !found ) continue;
    plot( var, "ghost " + histoname, min, max );
  }
}

//=============================================================================
// Methods of MyCounter
//=============================================================================
void PrVeloUTChecker::MyCounter::addCategory( std::string s ) {
  name.push_back( std::move( s ) );
  num.push_back( 0 );
  den.push_back( 0 );
}

void PrVeloUTChecker::MyCounter::count( bool found, std::vector<bool> const& flags ) {
  for ( unsigned int kk = 0; flags.size() > kk; ++kk ) {
    if ( !flags[kk] ) continue;
    ++den[kk];
    if ( !found ) continue;
    ++num[kk];
  }
}
