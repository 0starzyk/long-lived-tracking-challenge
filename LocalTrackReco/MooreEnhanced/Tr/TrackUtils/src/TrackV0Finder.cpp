/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "DetDesc/DetectorElement.h"
#include "Event/PrimaryVertex.h"
#include "Event/StateParameters.h"
#include "Event/TwoProngVertex.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ITrajPoca.h"
#include "Kernel/ParticleID.h"
#include "Kernel/ParticleProperty.h"
#include "LHCbAlgs/Transformer.h"
#include "TrackInterfaces/ITrackStateProvider.h"
#include "TrackInterfaces/ITrackVertexer.h"
#include "TrackKernel/TrackTraj.h"

#include "GaudiKernel/Vector4DTypes.h"

#include <iterator>
#include <numeric>

namespace {
  // output type of the algorithm
  using OutContainer = KeyedContainer<LHCb::TwoProngVertex, Containers::HashMap>;

  // set of constants
  constexpr double maxTrackChi2PerDoF{10};
  constexpr double zmin{-100 * Gaudi::Units::cm};
  constexpr double zmax{StateParameters::ZEndTT};
  constexpr double minVertexPVChi2{25};                   // 3 dofs
  constexpr double sigmaBFlightX{1.0 * Gaudi::Units::mm}; // approx RMS of flight length of B in X
  constexpr double sigmaBFlightY{1.0 * Gaudi::Units::mm}; // approx RMS of flight length of B in Y
  constexpr double sigmaBFlightZ{20 * Gaudi::Units::mm};  // approx RMS of flight length of B in Z
  constexpr double stateZTolerance{TrackParameters::propagationTolerance};
} // namespace

namespace LHCb {

  /**
   *  @author Wouter HULSBERGEN
   *  @date   2007-10-08
   */
  class TrackV0Finder : public Algorithm::Transformer<OutContainer( RecVertex::Range const&, Track::Range const&,
                                                                    DetectorElement const& ),
                                                      DetDesc::usesConditions<DetectorElement>> {
  public:
    TrackV0Finder( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer( name, pSvcLocator,
                       {{"PVContainer", RecVertexLocation::Primary},
                        {"TrackContainers", TrackLocation::Default},
                        {"StandardGeometryTop", standard_geometry_top}},
                       {"V0Container", RecVertexLocation::V0} ) {}
    StatusCode   initialize() override;
    OutContainer operator()( RecVertex::Range const&, Track::Range const&, DetectorElement const& ) const override;

  private:
    void constrainToVertex( const Gaudi::XYZPoint& pos, const Gaudi::LorentzVector& p4,
                            const Gaudi::SymMatrix7x7& cov7x7, const VertexBase& pv, double& chi2, double& decaylength,
                            double& decaylengtherr ) const;

    template <class PVContainer>
    bool hasV0Topology( TwoProngVertex& vertex, const PVContainer& pvs, double& decaylength ) const;

    ToolHandle<ITrajPoca>               m_pocatool{"TrajPoca"};
    ToolHandle<ITrackVertexer>          m_vertexer{"TrackVertexer"};
    ToolHandle<ITrackStateProvider>     m_stateprovider{this, "StateProvider", "TrackStateProvider"};
    ServiceHandle<IParticlePropertySvc> m_propertysvc{this, "ParticlePropertySvc", "ParticlePropertySvc"};
    const ParticleProperty*             m_ksProperty{nullptr};
    const ParticleProperty*             m_lambdaProperty{nullptr};
    const ParticleProperty*             m_pionProperty{nullptr};
    const ParticleProperty*             m_protonProperty{nullptr};

    Gaudi::Property<int>    m_maxNumCommonHits{this, "MaxNumCommonHits", -1};
    Gaudi::Property<double> m_minDeltaZ{this, "MinDeltaZ", 1 * Gaudi::Units::cm};
    Gaudi::Property<double> m_maxDocaDD{this, "MaxDocaDD", 1 * Gaudi::Units::mm};
    Gaudi::Property<double> m_maxDocaLL{this, "MaxDocaLL", 5 * Gaudi::Units::mm};
    Gaudi::Property<double> m_maxChi2V0Vertex{this, "MaxChi2V0Vertex", 16};         // 1 dof
    Gaudi::Property<double> m_maxChi2PVConstraint{this, "MaxChi2PVConstraint", 25}; // 2 dofs
    Gaudi::Property<double> m_minDecayLengthSignificance{this, "MinDecayLengthSignificance", 5};
    Gaudi::Property<bool>   m_correctBFlight{this, "CorrectBFlight", true};
    Gaudi::Property<double> m_ksmasscutLL{this, "KsMassCutLL", 40 * Gaudi::Units::MeV};
    Gaudi::Property<double> m_ksmasscutDD{this, "KsMassCutDD", 60 * Gaudi::Units::MeV};
    Gaudi::Property<double> m_lambdamasscut{this, "LambdaMassCut", 20 * Gaudi::Units::MeV};
    Gaudi::Property<double> m_minCTauKs{this, "MinCTauKs", -1};         // ctau =27 mm
    Gaudi::Property<double> m_minCTauLambda{this, "MinCTauLambda", -1}; // ctau =79 mm
    Gaudi::Property<bool>   m_useExtrapolator{this, "UseExtrapolator", true};
    Gaudi::Property<bool>   m_excludePVTracks{this, "ExcludePVTracks", true};
    Gaudi::Property<bool>   m_rejectUpstreamHits{this, "RejectUpstreamHits", true};
    Gaudi::Property<bool>   m_addStateAtVertex{this, "AddStateAtVertex", true};

    mutable Gaudi::Accumulators::StatCounter<>            m_nSelected{this, "numselected"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_trajPocaFailure{this, "TrajPoca Failure"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_vertexerFailure{this, "TrackVertexer Failure"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::WARNING> m_extrapolationFailure{
        this, "Extrapolation failed. Rely on trajectory interpolation"};
  };

  // Declaration of the Algorithm Factory
  DECLARE_COMPONENT_WITH_ID( TrackV0Finder, "TrackV0Finder" )

} // namespace LHCb

StatusCode LHCb::TrackV0Finder::initialize() {
  return Transformer::initialize().andThen( [&] {
    m_ksProperty     = m_propertysvc->find( "KS0" );
    m_lambdaProperty = m_propertysvc->find( "Lambda0" );
    m_pionProperty   = m_propertysvc->find( "pi+" );
    m_protonProperty = m_propertysvc->find( "p+" );
    if ( !m_ksProperty || !m_lambdaProperty || !m_pionProperty || !m_protonProperty ) {
      throw GaudiException( "Did not find all properties.", "TrackV0Finder::initialize", StatusCode::FAILURE );
    }
    return StatusCode::SUCCESS;
  } );
}

namespace {

  template <class PVContainer>
  bool inAnyVertex( const LHCb::Track& track, const PVContainer& pvs ) {
    return std::any_of( pvs.begin(), pvs.end(),
                        [id = LHCb::PrimaryVertex::uniqueVeloSegmentID( track )](
                            typename PVContainer::const_reference pv ) { return pv->contains( id ); } );
  }

  Gaudi::Vector3 transform( const Gaudi::XYZVector& vec ) { return {vec.X(), vec.Y(), vec.Z()}; }

  struct TrackTrajPair final {
    TrackTrajPair( const LHCb::Track* atrack, const LHCb::TrackTraj* atraj )
        : track( atrack ), traj( atraj ) {} // FIXME: C++20 -- delete this line...
    const LHCb::Track*     track;
    const LHCb::TrackTraj* traj;
  };

} // namespace

OutContainer LHCb::TrackV0Finder::operator()( RecVertex::Range const& oldpvcontainer, Track::Range const& tracks,
                                              DetectorElement const& lhcb ) const {
  // Get the defaulf geometry
  auto& geometry = *lhcb.geometry();
  // Create the output container
  OutContainer v0container;
  // Get the primary vertices.
  std::vector<std::unique_ptr<PrimaryVertex>> pvcontainer;
  pvcontainer.reserve( oldpvcontainer.size() );
  std::transform( oldpvcontainer.begin(), oldpvcontainer.end(), std::back_inserter( pvcontainer ),
                  []( const RecVertex* pv ) { return std::make_unique<PrimaryVertex>( *pv ); } );
  // Locate the one that's most downstream.
  double zprimary = ( !pvcontainer.empty() ? std::accumulate( pvcontainer.begin(), pvcontainer.end(), 9999.0,
                                                              []( double z, const std::unique_ptr<PrimaryVertex>& pv ) {
                                                                return std::min( z, pv->position().z() );
                                                              } )
                                           : -1000.0 );

  // Sort them by charge, make some cuts
  std::vector<TrackTrajPair> postracks, negtracks;
  // Get the Tracks
  for ( const auto& track : tracks ) {
    // require tracks with T and (TT or Velo)
    if ( ( track->type() == Track::Types::Long || track->type() == Track::Types::Downstream ) &&
         track->chi2PerDoF() < maxTrackChi2PerDoF &&
         // remove tracks from PVs, if required
         ( !m_excludePVTracks.value() || !inAnyVertex( *track, pvcontainer ) ) ) {
      const TrackTraj* traj = m_stateprovider->trajectory( *track, geometry );
      if ( !traj ) continue;
      if ( track->firstState().qOverP() > 0 )
        postracks.emplace_back( track, traj );
      else
        negtracks.emplace_back( track, traj );
    }
  }

  const double pimass     = m_pionProperty->mass();
  const double pmass      = m_protonProperty->mass();
  const double ksmass     = m_ksProperty->mass();
  const double lambdamass = m_lambdaProperty->mass();

  for ( const auto& ipos : postracks )
    for ( const auto& ineg : negtracks )
      if ( ( m_maxNumCommonHits < 0 ) || ( int( ipos.track->nCommonLhcbIDs( *ineg.track ) ) <= m_maxNumCommonHits ) ) {

        const TrackTraj& postraj = *ipos.traj;
        const TrackTraj& negtraj = *ineg.traj;

        // Calls pocatool
        double           mupos( 0 ), muneg( 0 );
        Gaudi::XYZVector deltaX;
        StatusCode       sc = m_pocatool->minimize( postraj, mupos, negtraj, muneg, deltaX, 0.001 * Gaudi::Units::mm );
        if ( sc.isFailure() ) {
          ++m_trajPocaFailure;
          continue;
        }
        double distance = deltaX.R();
        double z        = 0.5 * ( mupos + muneg );
        // make the cut on the distance and the z-position
        bool isLL = ipos.track->hasVelo() && ineg.track->hasVelo();
        if ( !( ( distance < m_maxDocaLL || ( distance < m_maxDocaDD && !isLL ) ) && zmin < z && z < zmax &&
                ( m_excludePVTracks.value() || zprimary + m_minDeltaZ < z ) ) ) {
          continue;
        }

        // now make an invariant mass cut. resolution is only a
        // bit worse, and it really beats down combinatorics at
        // little cost.

        // Gaudi::XYZVector mompos = vertex->p3A() ;
        Gaudi::XYZVector mompos = postraj.momentum( mupos );
        // Gaudi::XYZVector momneg = vertex->p3B() ;
        Gaudi::XYZVector     momneg = negtraj.momentum( muneg );
        Gaudi::LorentzVector p4pos( mompos.X(), mompos.Y(), mompos.Z(), std::sqrt( mompos.Mag2() + pimass * pimass ) );
        Gaudi::LorentzVector p4neg( momneg.X(), momneg.y(), momneg.Z(), std::sqrt( momneg.Mag2() + pimass * pimass ) );
        Gaudi::LorentzVector p4pipi   = p4pos + p4neg;
        double               pipimass = p4pipi.M();
        p4pos.SetE( std::sqrt( mompos.Mag2() + pmass * pmass ) );
        Gaudi::LorentzVector p4ppi   = p4pos + p4neg;
        double               ppimass = p4ppi.M();
        p4pos.SetE( std::sqrt( mompos.Mag2() + pimass * pimass ) );
        p4neg.SetE( std::sqrt( momneg.Mag2() + pmass * pmass ) );
        Gaudi::LorentzVector p4pip         = p4pos + p4neg;
        double               pipmass       = p4pip.M();
        double               mom           = p4pipi.P();
        bool                 iskscandidate = std::abs( pipimass - ksmass ) < ( isLL ? m_ksmasscutLL : m_ksmasscutDD );
        bool                 islambdacandidate     = std::abs( ppimass - lambdamass ) < m_lambdamasscut;
        bool                 isantilambdacandidate = std::abs( pipmass - lambdamass ) < m_lambdamasscut;

        if ( !( iskscandidate || islambdacandidate || isantilambdacandidate ) ) continue;

        State                           posstate = postraj.state( z );
        State                           negstate = negtraj.state( z );
        std::unique_ptr<TwoProngVertex> vertex{m_vertexer->fit( posstate, negstate, geometry )};

        double decaylength;
        if ( !vertex ) {
          ++m_vertexerFailure;
          continue;
        }

        if ( !( vertex->chi2() < m_maxChi2V0Vertex && hasV0Topology( *vertex, pvcontainer, decaylength ) ) ) continue;

        // cut on ctau
        if ( m_minCTauKs > 0 || m_minCTauLambda > 0 ) {
          iskscandidate         = iskscandidate && decaylength * pipimass / mom > m_minCTauKs;
          islambdacandidate     = islambdacandidate && decaylength * ppimass / mom > m_minCTauLambda;
          isantilambdacandidate = isantilambdacandidate && decaylength * pipmass / mom > m_minCTauLambda;
        }

        // one last check: test that there are no hits upstream of the vertex on either track
        bool hasUpstreamHits = false;
        if ( m_rejectUpstreamHits.value() ) {
          const State* mstatepos = ipos.track->stateAt( State::Location::FirstMeasurement );
          const State* mstateneg = ineg.track->stateAt( State::Location::FirstMeasurement );
          double       zveto     = vertex->position().z() - 5 * std::sqrt( vertex->covMatrix()( 2, 2 ) );
          hasUpstreamHits        = ( mstatepos && mstatepos->z() < zveto ) || ( mstateneg && mstateneg->z() < zveto );
        }
        if ( !( !hasUpstreamHits && ( iskscandidate || islambdacandidate || isantilambdacandidate ) ) ) continue;

        // The states above were obtained from the
        // trajectory approximation. We actually prefer to
        // use real extrapolators. We do that after the
        // final selection to save time.
        if ( m_useExtrapolator.value() ) {
          sc = m_stateprovider->state( posstate, *( ipos.track ), z, geometry, stateZTolerance );
          if ( sc.isSuccess() )
            posstate.setLocation( State::Location::V0Vertex );
          else {
            ++m_extrapolationFailure;
            posstate = postraj.state( z );
          }
          sc = m_stateprovider->state( negstate, *( ineg.track ), z, geometry, stateZTolerance );
          if ( sc.isSuccess() )
            negstate.setLocation( State::Location::V0Vertex );
          else {
            ++m_extrapolationFailure;
            negstate = negtraj.state( z );
          }
          std::unique_ptr<TwoProngVertex> newvertex{m_vertexer->fit( posstate, negstate, geometry )};
          if ( newvertex ) vertex = std::move( newvertex );
        }

        vertex->addToTracks( ipos.track );
        vertex->addToTracks( ineg.track );
        if ( iskscandidate ) {
          auto pid = ParticleID( m_ksProperty->pdgID() );
          vertex->addPID( pid );
        }
        if ( islambdacandidate ) {
          auto pid = ParticleID( m_lambdaProperty->pdgID() );
          vertex->addPID( pid );
        }
        if ( isantilambdacandidate ) {
          auto pid = ParticleID( m_lambdaProperty->antiParticle()->pdgID() );
          vertex->addPID( pid );
        }
        // transfer ownership to the output container
        v0container.add( vertex.release() );
        // finally add the states to the tracks, provided the extrapolation was done properly
        // (the 'Vertex' location was only set for properly interpolated tracks.
        if ( m_addStateAtVertex.value() ) {
          if ( posstate.location() == State::Location::V0Vertex )
            ( const_cast<Track*>( ipos.track ) )->addToStates( posstate );
          if ( negstate.location() == State::Location::V0Vertex )
            ( const_cast<Track*>( ineg.track ) )->addToStates( negstate );
        }
      }

  m_nSelected += v0container.size();
  return v0container;
}

void LHCb::TrackV0Finder::constrainToVertex( const Gaudi::XYZPoint& pos, const Gaudi::LorentzVector& p4,
                                             const Gaudi::SymMatrix7x7& cov7, const LHCb::VertexBase& pv, double& chi2,
                                             double& decaylength, double& decaylengtherr ) const {
  // This calculation is basically a 1-iteration beamspot fit. The
  // constraint is
  //
  //    r = x - lambda p/|p| - xbs
  //
  // where x and p are the position of the decay vertex of the
  // candidate and its momentum, lambda is the decaylength and xbs
  // the position of the beamspot. The covariance in the constraint
  // is
  //
  //    V = Vbs + Vxx - a * Vxp - a Vxp^T + a^2 * Vpp
  //
  // where a=lambda/|p|^2. It needs an initial estimate for the
  // flightlength, for which we simply take the projection of deltaX
  // on the direction. We now minimize  the chisquare contribution
  //
  //     chi^2 = r^T V^{-1} r
  //
  // for lambda.

  Gaudi::Vector3   dx    = transform( pos - pv.position() );
  Gaudi::XYZVector p3    = p4.Vect();
  double           p3mag = p3.R();
  Gaudi::Vector3   dir   = transform( p3.Unit() );

  Gaudi::SymMatrix3x3 W = pv.covMatrix(); // we'll repace this with a constant error that contains B motion

  if ( m_correctBFlight.value() ) {
    // For determining whether this candidate is compatible with the
    // PV, we want to take into account the the Ks may come from a
    // B. So, we need to add something to the PV vertex. However, we
    // want to count the B flight length in the decay length, so we
    // don't want to add the B flight length to the PV error in z. The
    // trick is to remove the contribution of B-flight along the V0
    // direction.

    Gaudi::SymMatrix3x3 covBFlight;
    covBFlight( 0, 0 ) = sigmaBFlightX * sigmaBFlightX;
    covBFlight( 1, 1 ) = sigmaBFlightY * sigmaBFlightY;
    covBFlight( 2, 2 ) = sigmaBFlightZ * sigmaBFlightZ;

    // now project out the part perpendicular to the direction.
    //   W +=  (1-P) * covB * (1-P)
    // where P is the projection matrix
    //   P_ij = dir_i * dir_j
    //
    // I am sure that there is something left to optimize here ...
    Gaudi::SymMatrix3x3 UnitMinusP;
    UnitMinusP( 0, 0 ) = UnitMinusP( 1, 1 ) = UnitMinusP( 2, 2 ) = 1;
    for ( size_t irow = 0; irow < 3; ++irow )
      for ( size_t jrow = 0; jrow <= irow; ++jrow ) UnitMinusP( irow, jrow ) -= dir( irow ) * dir( jrow );

    // here we could use that W is diagonal. that saves a lot of time
    // W +=  ROOT::Math::Similarity(covBFlight,UnitMinusP) ;
    for ( size_t irow = 0; irow < 3; ++irow )
      for ( size_t jrow = 0; jrow <= irow; ++jrow )
        for ( size_t krow = 0; krow < 3; ++krow )
          W( irow, jrow ) += UnitMinusP( irow, krow ) * covBFlight( krow, krow ) * UnitMinusP( krow, jrow );
  }

  // double a = (ROOT::Math::Transpose(dir)*dx)/p3mag  ;
  double a = ROOT::Math::Dot( dir, dx ) / p3mag;
  for ( size_t row = 0; row < 3; ++row )
    for ( size_t col = 0; col <= row; ++col )
      W( row, col ) +=
          cov7( row, col ) + a * a * cov7( row + 3, col + 3 ) - a * ( cov7( row + 3, col ) + cov7( col + 3, row ) );

  int OK = W.InvertChol();
  if ( !OK ) error() << "inversion error in constrainToVertex" << endmsg;

  double halfdChi2dLam2 = ROOT::Math::Similarity( W, dir );
  decaylength           = ROOT::Math::Dot( dir, W * dx ) / halfdChi2dLam2;
  decaylengtherr        = std::sqrt( 1 / halfdChi2dLam2 );

  Gaudi::Vector3 res = dx - decaylength * dir;

  chi2 = ROOT::Math::Similarity( W, res );
}

template <class PVContainer>
bool LHCb::TrackV0Finder::hasV0Topology( LHCb::TwoProngVertex& vertex, const PVContainer& pvs,
                                         double& decaylength ) const {
  // returns true if V0 candidate accepted. we veto two types of background:
  // * V0 candidates that do not point to any PV ( chi2 > m_maxChi2PVConstraint)
  // * V0 candidates that point to a PV but with too small decay length

  bool accept = true;

  if ( !pvs.empty() ) {
    double bestipchi2( -1 ); // ip chi2
    double bestvvchi2( -1 ); // vertex-vertex chi2
    double decaylengtherr = 1;
    decaylength           = 0;

    // first find the best fitting primary vertex
    const Gaudi::XYZPoint& pos    = vertex.position();
    Gaudi::SymMatrix7x7    cov7x7 = vertex.covMatrix7x7( 0, 0 );
    Gaudi::LorentzVector   p4     = vertex.momentum( 0, 0 );
    for ( const auto& pv : pvs ) {
      double thisipchi2, thisdecaylength, thisdecaylengtherr;
      constrainToVertex( pos, p4, cov7x7, *pv, thisipchi2, thisdecaylength, thisdecaylengtherr );
      // is this the best fitting PV?
      if ( bestipchi2 < 0 || thisipchi2 < bestipchi2 ) {
        bestipchi2     = thisipchi2;
        decaylength    = thisdecaylength;
        decaylengtherr = thisdecaylengtherr;
      }
      double thisvvchi2 = thisipchi2 + std::pow( thisdecaylength / thisdecaylengtherr, 2 );
      if ( bestvvchi2 < 0 || thisvvchi2 < bestvvchi2 ) bestvvchi2 = thisvvchi2;
    }
    accept = bestvvchi2 > minVertexPVChi2 && bestipchi2 < m_maxChi2PVConstraint &&
             decaylength / decaylengtherr > m_minDecayLengthSignificance;
  }

  return accept;
}
