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
#include "Event/TwoProngVertex.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ITrajPoca.h"
#include "Kernel/ParticleProperty.h"
#include "Magnet/DeMagnet.h"
#include "TrackInterfaces/ITrackVertexer.h"
#include "TrackKernel/TrackTraj.h"

#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "LHCbAlgs/Transformer.h"

#include "gsl/gsl_cdf.h"

namespace {
  //=============================================================================
  // Check of two tracks do not have two much overlap to be combined as TwoProng
  //=============================================================================

  bool overlap( const LHCb::Track& trackA, const LHCb::Track& trackB ) {
    // for now, just look at common ancestors
    return std::any_of( trackA.ancestors().begin(), trackA.ancestors().end(), [&]( const auto& ancA ) {
      return std::any_of( trackB.ancestors().begin(), trackB.ancestors().end(),
                          [&]( const auto& ancB ) { return ancA == ancB; } );
    } );
  }

} // namespace
using TwoProngVertices = KeyedContainer<LHCb::TwoProngVertex, Containers::HashMap>;

class TrackDiMuonMonitor final : public LHCb::Algorithm::Transformer<
                                     TwoProngVertices( LHCb::Tracks const&, DetectorElement const&, DeMagnet const& ),
                                     LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DetectorElement, DeMagnet>> {
public:
  /** Standard construtor */
  TrackDiMuonMonitor( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer{name,
                    pSvcLocator,
                    {KeyValue{"TrackLocation", LHCb::TrackLocation::Muon},
                     KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top},
                     KeyValue{"Magnet", LHCb::Det::Magnet::det_path}},
                    {"DiMuonLocation", "Rec/Vertex/DiMuon"}} {}

  /** Algorithm initialize */
  StatusCode initialize() override;

  /** Algorithm execute */
  TwoProngVertices operator()( LHCb::Tracks const&, DetectorElement const&, DeMagnet const& ) const override;

private:
  bool accept( const LHCb::TwoProngVertex& vertex ) const;

private:
  PublicToolHandle<ITrajPoca>      m_pocatool{this, "PocaTool", "TrajPoca"};
  PublicToolHandle<ITrackVertexer> m_vertexer{this, "TrackVertexer", "TrackVertexer"};
  double                           m_maxDistance = 5 * Gaudi::Units::mm;
  Gaudi::Property<double>          m_minMass{this, "MinMass", 2.6 * Gaudi::Units::GeV};
  Gaudi::Property<double>          m_maxMass{this, "MaxMass", 4.0 * Gaudi::Units::GeV};
  double                           m_minJPsiMass = 3.065 * Gaudi::Units::GeV;
  double                           m_maxJPsiMass = 3.125 * Gaudi::Units::GeV;
  Gaudi::Property<double>          m_maxChi2TwoProngVertex{this, "MaxChi2TwoProngVertex", 25};
  double                           m_muonmass = 0.;

  mutable Gaudi::Accumulators::AveragingCounter<> m_num_combinations{this, "numcombinations"};
  mutable Gaudi::Accumulators::AveragingCounter<> m_num_postracks{this, "numpostracks"};
  mutable Gaudi::Accumulators::AveragingCounter<> m_num_negtracks{this, "numnegtracks"};
  mutable Gaudi::Accumulators::AveragingCounter<> m_num_selected{this, "numselected"};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TrackDiMuonMonitor )

StatusCode TrackDiMuonMonitor::initialize() {
  return Transformer::initialize().andThen( [&] {
    auto propertysvc = service<LHCb::IParticlePropertySvc>( "LHCb::ParticlePropertySvc", true );
    if ( !propertysvc ) return StatusCode::FAILURE;
    const auto* muon = propertysvc->find( "mu+" );
    if ( !muon ) return StatusCode::FAILURE;
    m_muonmass = muon->mass();
    setHistoTopDir( "Track/" );
    return StatusCode::SUCCESS;
  } );
}

TwoProngVertices TrackDiMuonMonitor::operator()( LHCb::Tracks const& muontracks, DetectorElement const& lhcb,
                                                 DeMagnet const& magneticfield ) const {

  // Sort them by charge, make some cuts
  std::vector<const LHCb::Track*> postracks, negtracks;
  for ( auto track : muontracks ) {
    // require tracks with T and (TT or Velo)
    bool accepted = track->hasT() && ( track->hasVelo() || track->hasUT() );
    if ( accepted && !track->ancestors().empty() ) {
      // muon track: take the first ancestor
      track = track->ancestors().front();
    } else {
      track = nullptr;
    }
    if ( track ) { ( track->firstState().qOverP() > 0 ? postracks : negtracks ).push_back( track ); }
  }

  // turn them into trajectories
  std::vector<LHCb::TrackTraj> postrajs;
  std::vector<LHCb::TrackTraj> negtrajs;
  postrajs.reserve( postracks.size() );
  negtrajs.reserve( negtracks.size() );
  for ( const auto& pos : postracks ) postrajs.emplace_back( *pos, &magneticfield );
  for ( const auto& neg : negtracks ) negtrajs.emplace_back( *neg, &magneticfield );

  m_num_combinations += postracks.size() * negtracks.size();
  m_num_postracks += postracks.size();
  m_num_negtracks += negtracks.size();

  // Create the output container
  TwoProngVertices dimuoncontainer;

  for ( auto ipos = postracks.begin(); ipos != postracks.end(); ++ipos )
    for ( auto ineg = negtracks.begin(); ineg != negtracks.end(); ++ineg ) {

      if ( overlap( **ipos, **ineg ) ) continue;
      const LHCb::TrackTraj& postraj = postrajs[std::distance( postracks.begin(), ipos )];
      const LHCb::TrackTraj& negtraj = negtrajs[std::distance( negtracks.begin(), ineg )];

      double z_seed = postraj.beginRange() + 1 * Gaudi::Units::mm;
      double mupos( z_seed ), muneg( z_seed );

      // Calls pocatool
      Gaudi::XYZVector deltaX;
      StatusCode       sc = m_pocatool->minimize( postraj, mupos, negtraj, muneg, deltaX, 0.001 * Gaudi::Units::mm );

      if ( !sc.isSuccess() ) {
        Warning( "TrajPoca Failure", StatusCode::SUCCESS, 0 ).ignore();
        continue;
      }
      double distance = deltaX.R();
      double z        = 0.5 * ( mupos + muneg );
      if ( distance > m_maxDistance ) continue;

      // now make an invariant mass cut
      Gaudi::XYZVector     mompos  = postraj.momentum( mupos );
      Gaudi::XYZVector     momneg  = negtraj.momentum( muneg );
      double               mumass2 = m_muonmass * m_muonmass;
      Gaudi::LorentzVector p4pos( mompos.X(), mompos.Y(), mompos.Z(), std::sqrt( mompos.Mag2() + mumass2 ) );
      Gaudi::LorentzVector p4neg( momneg.X(), momneg.y(), momneg.Z(), std::sqrt( momneg.Mag2() + mumass2 ) );
      double               dimuonmass = ( p4pos + p4neg ).M();

      if ( dimuonmass < m_minMass || m_maxMass < dimuonmass ) continue;

      // finally, create the vertex and cut on the chisquare
      LHCb::State posstate = postraj.state( z );
      LHCb::State negstate = negtraj.state( z );
      auto        vertex   = m_vertexer->fit( posstate, negstate, *lhcb.geometry() );
      if ( !vertex ) continue;
      vertex->addToTracks( *ipos );
      vertex->addToTracks( *ineg );
      if ( !accept( *vertex ) ) continue;
      dimuoncontainer.add( vertex.release() );
    }

  m_num_selected += dimuoncontainer.size();

  plot( dimuoncontainer.size(), "multiplicity", "J/psi candidate multiplicity", -0.5, 20.5, 21 );

  return dimuoncontainer;
}

bool TrackDiMuonMonitor::accept( const LHCb::TwoProngVertex& vertex ) const {
  bool          rc( false );
  static double chi2max  = gsl_cdf_chisq_Qinv( 1e-6, 1 );
  double        chi2     = vertex.chi2();
  double        chi2prob = chi2 < chi2max ? gsl_cdf_chisq_Q( chi2, 1 ) : 0;
  plot( chi2prob, "chi2prob", "chi2prob", 0, 1 );

  if ( vertex.chi2() < m_maxChi2TwoProngVertex ) {
    Gaudi::LorentzVector p4A  = vertex.p4A( m_muonmass );
    Gaudi::LorentzVector p4B  = vertex.p4B( m_muonmass );
    Gaudi::LorentzVector p4   = p4A + p4B;
    double               mass = p4.M();
    plot( mass, "mass", "dimuon mass", m_minMass, m_maxMass );

    if ( m_minJPsiMass < mass && mass < m_maxJPsiMass ) {

      plot( vertex.position().x(), "vxJPsi", "J/psi candidate vertex x", -2, 2 );
      plot( vertex.position().y(), "vyJPsi", "J/psi candidate vertex y", -2, 2 );
      plot( vertex.position().z(), "vzJPsi", "J/psi candidate vertex z", -200, 200 );
      plot( chi2prob, "chi2probJPsi", "J/psi candidate chi2prob", 0, 1 );
      plot( mass, "massJPsi", "J/psi candidate mass", 3.05 * Gaudi::Units::GeV, 3.15 * Gaudi::Units::GeV );
      const LHCb::Track* postrack = vertex.trackA();
      const LHCb::Track* negtrack = vertex.trackB();
      if ( postrack->firstState().qOverP() < 0 ) std::swap( postrack, negtrack );
      double ppos = postrack->firstState().momentum().R();
      double pneg = negtrack->firstState().momentum().R();
      double pdif = ppos - pneg;
      profile1D( pdif, mass, "massVersusMomDif", "dimuon mass versus p_{pos} - p_{neg}", -50 * Gaudi::Units::GeV,
                 50 * Gaudi::Units::GeV, 20, "", 3.065 * Gaudi::Units::GeV, 3.125 * Gaudi::Units::GeV );
      plot( p4.P(), "momJPsi", "JPsi candidate momentum", 0, 100 * Gaudi::Units::GeV );
      plot( pdif, "momdifJPsi", "p_{pos} - p_{neg} for JPsi candidates", -50 * Gaudi::Units::GeV,
            50 * Gaudi::Units::GeV );
      plot( p4.Pt(), "ptJPsi", "JPsi candidate Pt", 0, 10 * Gaudi::Units::GeV );
      plot( p4.Eta(), "etaJPsi", "JPsi candidate eta", 2, 5 );

      rc = true;
    }
  }
  return rc;
}
