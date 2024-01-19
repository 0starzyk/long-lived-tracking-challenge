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
#include "Event/Particle.h"
#include "Event/ProtoParticle.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ParticleProperty.h"
#include "LHCbAlgs/Consumer.h"

class TrackTune : public LHCb::Algorithm::Consumer<void( const LHCb::Track::Range&, const LHCb::Particle::Range& ),
                                                   Gaudi::Functional::Traits::BaseClass_t<GaudiTupleAlg>> {

public:
  TrackTune( const std::string& name, ISvcLocator* pSvc );

  StatusCode initialize() override;

  void operator()( const LHCb::Track::Range&, const LHCb::Particle::Range& ) const override;

private:
  const LHCb::Track* track( const LHCb::Particle& part ) const;

  bool isFound( const LHCb::Track::Range& tracks, const LHCb::Particle& part ) const;

  std::vector<const LHCb::Particle*> select( const LHCb::Particle::Range& input ) const;

  bool inMassRange( const LHCb::Particle& particle ) const;

  Gaudi::Property<std::string> m_resonanceName{this, "resonanceName", "J/psi(1S)"};

  Gaudi::Property<double> m_deltaMass{this, "resonance", 110.};
  double                  m_minMass = 0;
  double                  m_maxMass = 0;
  Gaudi::Property<bool>   m_selectBest{this, "selectBest", true};
  Gaudi::Property<double> m_minPurityCut{this, "minPurityCut", 0.7};
};

DECLARE_COMPONENT( TrackTune )

TrackTune::TrackTune( const std::string& name, ISvcLocator* pSvc )
    : Consumer{name,
               pSvc,
               {KeyValue{"ParticleLocation", "/Event/Dimuon/Phys/SelDiMuonInciLoose/Particles"},
                KeyValue{"TrackLocation", LHCb::TrackLocation::Default}}} {}

StatusCode TrackTune::initialize() {

  static const std::string histoDir = "Track/";
  if ( "" == histoTopDir() ) setHistoTopDir( histoDir );

  // Mandatory initialization
  StatusCode sc = Consumer::initialize();
  if ( sc.isFailure() ) { return sc; }

  LHCb::IParticlePropertySvc*   propertysvc = svc<LHCb::IParticlePropertySvc>( "LHCb::ParticlePropertySvc", true );
  const LHCb::ParticleProperty* prop        = propertysvc->find( m_resonanceName );
  if ( !prop ) { return Error( "Failed to find resonance", StatusCode::SUCCESS ); }

  m_minMass = prop->mass() - m_deltaMass;
  m_maxMass = prop->mass() + m_deltaMass;
  propertysvc->release();

  info() << "MinMass " << m_minMass << " MaxMass " << m_maxMass << endmsg;

  return StatusCode::SUCCESS;
}

void TrackTune::operator()( const LHCb::Track::Range& tracks, const LHCb::Particle::Range& particles ) const {
  Tuple myTuple = nTuple( "Candidates" );
  for ( const auto* t : select( particles ) ) {

    myTuple << Tuples::Column( "M", t->measuredMass() ) << Tuples::Column( "found", isFound( tracks, *t ) )
            << Tuples::Column( "PT", t->pt() ) << Tuples::Column( "Candidates", particles.size() );

    myTuple->write().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  }
}

const LHCb::Track* TrackTune::track( const LHCb::Particle& part ) const {
  const LHCb::ProtoParticle* proto = part.proto();
  return ( !proto || proto->charge() == 0 ) ? nullptr : proto->track();
}

bool TrackTune::isFound( const LHCb::Track::Range& tracks, const LHCb::Particle& part ) const {
  bool        ok        = true;
  const auto& daughters = part.daughters();
  for ( auto iter = daughters.begin(); iter != daughters.end() && ok; ++iter ) {
    const LHCb::Track* aTrack = track( **iter );
    if ( !aTrack ) info() << "Failed to find track " << endmsg;
    const double nHits        = aTrack->nLHCbIDs();
    bool         matchedTrack = false;
    for ( auto iterT = tracks.begin(); iterT != tracks.end() && !matchedTrack; ++iterT ) {
      const double fracCommon = aTrack->nCommonLhcbIDs( **iterT ) / double( nHits );
      plot( fracCommon, "purity", "purity", 0., 2., 100 );
      if ( fracCommon > m_minPurityCut ) matchedTrack = true;
    } // tracks
    if ( !matchedTrack ) ok = false;
  } // particles

  return ok;
}

std::vector<const LHCb::Particle*> TrackTune::select( const LHCb::Particle::Range& input ) const {
  if ( !m_selectBest ) return {input.begin(), input.end()};

  double                bestChi2 = 9999.;
  const LHCb::Particle* bestPart = nullptr;
  for ( const auto& i : input ) {
    if ( !inMassRange( *i ) ) continue;
    auto vert = i->endVertex();
    if ( vert->chi2PerDoF() < bestChi2 ) {
      bestChi2 = vert->chi2PerDoF();
      bestPart = i;
    }
  }
  if ( bestPart ) return {bestPart};
  return {};
}

bool TrackTune::inMassRange( const LHCb::Particle& particle ) const {
  double m = particle.measuredMass();
  return ( m > m_minMass && m < m_maxMass );
}
