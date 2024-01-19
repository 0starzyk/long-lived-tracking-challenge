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
#include "Event/TwoProngVertex.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/ToolHandle.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ParticleProperty.h"
#include "LHCbAlgs/Consumer.h"
#include "TrackInterfaces/ITrackVertexer.h"
#include "gsl/gsl_cdf.h"

class TrackV0Monitor
    : public LHCb::Algorithm::Consumer<void( const LHCb::TwoProngVertices&, const LHCb::RecVertex::Range& ),
                                       Gaudi::Functional::Traits::BaseClass_t<GaudiHistoAlg>> {
public:
  /** Standard construtor */
  TrackV0Monitor( const std::string& name, ISvcLocator* pSvcLocator );

  /** Algorithm initialize */
  StatusCode initialize() override;

  /** Algorithm execute */
  void operator()( const LHCb::TwoProngVertices&, const LHCb::RecVertex::Range& ) const override;

private:
  void process( const LHCb::TwoProngVertex& vertex, const LHCb::RecVertex::Range& pvs, const std::string& dir ) const;
  void computeDecayLength( const LHCb::TwoProngVertex& vertex, const LHCb::RecVertex::Range& pvs, double& ipchi2,
                           double& decaylength, double& decaylengtherr ) const;

private:
  ToolHandle<ITrackVertexer> m_vertexer{"TrackVertexer"};
  Gaudi::Property<double>    m_maxIPChi2{this, "MaxIPChi2", 25.0};
  double                     m_pionmass{0};
  double                     m_protonmass{0};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TrackV0Monitor )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TrackV0Monitor::TrackV0Monitor( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer{name,
               pSvcLocator,
               {KeyValue{"V0Location", LHCb::TwoProngVertexLocation::Default},
                KeyValue{"PrimaryLocation", LHCb::RecVertexLocation::Primary}}} {}

StatusCode TrackV0Monitor::initialize() {
  setHistoTopDir( "Track/" );
  StatusCode                  sc          = Consumer::initialize();
  LHCb::IParticlePropertySvc* propertysvc = svc<LHCb::IParticlePropertySvc>( "LHCb::ParticlePropertySvc", true );
  if ( !propertysvc ) return StatusCode::FAILURE;

  const LHCb::ParticleProperty* pion = propertysvc->find( "pi+" );
  if ( !pion ) return StatusCode::FAILURE;
  m_pionmass                           = pion->mass();
  const LHCb::ParticleProperty* proton = propertysvc->find( "p+" );
  if ( !proton ) return StatusCode::FAILURE;
  m_protonmass = proton->mass();
  // const LHCb::ParticleProperty* kshort = propertysvc->find( "KS0" ) ;
  //     m_k0pid = kshort ? LHCb::ParticleID(kshort->pdgID()) : 0 ;
  //     const LHCb::ParticleProperty* lambda = propertysvc->find( "Lambda0" );
  //     m_lambdapid     = lambda ? LHCb::ParticleID(lambda->pdgID()) : 0 ;
  //     m_lambdabarpid = lambda ? LHCb::ParticleID(lambda->antiParticle()->pdgID()) : 0 ;

  return sc;
}

void TrackV0Monitor::process( const LHCb::TwoProngVertex& vertex, const LHCb::RecVertex::Range& pvs,
                              const std::string& dir ) const {
  static double chi2max  = gsl_cdf_chisq_Qinv( 1e-6, 1 );
  double        chi2     = vertex.chi2();
  double        chi2prob = chi2 < chi2max ? gsl_cdf_chisq_Q( chi2, 1 ) : 0;
  double        pipimass = vertex.mass( m_pionmass, m_pionmass );
  double        ppimass  = vertex.mass( m_protonmass, m_pionmass );
  double        pipmass  = vertex.mass( m_pionmass, m_protonmass );
  if ( vertex.trackA()->firstState().qOverP() < 0 ) std::swap( pipmass, ppimass );

  // compute the ipchi2 and the decaylength
  double ipchi2( -1 ), decaylength( 0 ), decaylengtherr( 1 );
  computeDecayLength( vertex, pvs, ipchi2, decaylength, decaylengtherr );
  plot( ipchi2, dir + "/ipchi2", "ip chi2", 0, 100 );
  if ( ipchi2 < m_maxIPChi2 ) {
    plot( chi2prob, dir + "/chi2prob", "chi2prob", 0, 1 );
    plot( vertex.position().z(), dir + "/z", "z", -100, 2000 );
    plot( pipimass, dir + "/pipimass", "pi+ pi- mass", 450, 550 );
    plot( ppimass, dir + "/ppimass", "p+ pi- mass", 1090, 1140 );
    plot( pipmass, dir + "/pipmass", "p- pi+ mass", 1090, 1140 );
    plot( decaylength / decaylengtherr, dir + "/dls", "decay length significance", -10, 100 );
  }

  if ( fullDetail() ) {
    plot( pipimass, dir + "/pipimassall", "pi+ pi- mass", 450, 550 );
    plot( ppimass, dir + "/ppimassall", "p+ pi- mass", 1090, 1140 );
    plot( pipmass, dir + "/pipmassall", "p- pi+ mass", 1090, 1140 );
    plot( decaylength / decaylengtherr, dir + "/dlsall", "decay length significance", -10, 100 );
    plot( vertex.position().x(), dir + "/x", "x", -100, 100 );
    plot( vertex.position().y(), dir + "/y", "y", -100, 100 );
  }
}

void TrackV0Monitor::operator()( const LHCb::TwoProngVertices& v0container, const LHCb::RecVertex::Range& pvs ) const {
  plot( v0container.size(), "v0multiplicity", "v0multiplicity", -0.5, 20.5, 21 );
  for ( const LHCb::TwoProngVertex* vertex : v0container ) {
    // determine the type: Long-Long, DS-DS, or DS-Long
    const LHCb::Track::Types trackA = vertex->trackA()->type();
    const LHCb::Track::Types trackB = vertex->trackB()->type();
    auto                     type   = trackA == LHCb::Track::Types::Long && trackB == LHCb::Track::Types::Long
                    ? "LongLong"
                    : trackA == LHCb::Track::Types::Downstream && trackB == LHCb::Track::Types::Downstream
                          ? "DownstreamDownstream"
                          : "Mixed";
    process( *vertex, pvs, type );
  }
}

void TrackV0Monitor::computeDecayLength( const LHCb::TwoProngVertex& vertex, const LHCb::RecVertex::Range& pvs,
                                         double& ipchi2, double& decaylength, double& decaylengtherr ) const {
  // some default values
  ipchi2         = -1;
  decaylength    = 0;
  decaylengtherr = 1;

  for ( const auto& pv : pvs ) {
    double thisipchi2( -1 ), thisdecaylength( 0 ), thisdecaylengtherr( 1 );
    bool   success = m_vertexer->computeDecayLength( vertex, *pv, thisipchi2, thisdecaylength, thisdecaylengtherr );
    if ( success && ( ipchi2 < 0 || thisipchi2 < ipchi2 ) ) {
      ipchi2         = thisipchi2;
      decaylength    = thisdecaylength;
      decaylengtherr = thisdecaylengtherr;
    }
  }
}
