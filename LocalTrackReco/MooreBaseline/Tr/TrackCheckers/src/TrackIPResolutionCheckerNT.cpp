/*****************************************************************************\
* (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "Event/LinksByKey.h"
#include "Event/MCHeader.h"
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Event/MCVertex.h"
#include "Event/PrimaryVertices.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "Event/TrackVertexUtils.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "Kernel/HitPattern.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Consumer.h"
#include "TrackKernel/CubicStateInterpolationTraj.h"
#include "TrackKernel/TrackFunctors.h"
#include "boost/functional/hash.hpp"
#include <map>
/** @class TrackIPResolutionCheckerBase
 *
 *  This base class provides the logic to create nTuples containing
 *  MC and reconstructed information related to the impact parameter
 *  and its calculation. Two versions exist of this algorithm,
 *  one which uses MCHits (and thus requires xdigi/xdst files),
 *  named TrackIPResolutionCheckerNTMCHits, and one which only
 *  requires links between the MCParticles and the tracks,
 *  called TrackIPResolutionCheckerNT.
 *
 *  @author W. Hulsbergen, P. Tsopelas, L.Dufour (modernisation only)
 */

namespace {
  using MCHitVector = std::vector<const LHCb::MCHit*>;
  using MCHitMap    = std::map<const LHCb::MCParticle*, MCHitVector>;
  // FIXME: we could template this on the PV container, if needed.
  using Vertices = LHCb::Event::PV::PrimaryVertexContainer;

  const LHCb::MCVertex* findMCOriginVertex( const LHCb::MCParticle& particle, double decaylengthtolerance ) {
    // take this particle, walk up the decay tree and determine its originvertex
    const LHCb::MCVertex* originvertex = particle.originVertex();
    if ( originvertex ) {
      const LHCb::MCParticle* mother = originvertex->mother();
      if ( mother && mother != &particle ) {
        const LHCb::MCVertex* motheroriginvertex = mother->originVertex();
        if ( motheroriginvertex &&
             ( motheroriginvertex == originvertex ||
               ( ( motheroriginvertex->position() - originvertex->position() ).R() ) < decaylengthtolerance ) )
          originvertex = findMCOriginVertex( *mother, decaylengthtolerance );
      }
    }
    return originvertex;
  }

  int mcVertexType( const LHCb::MCVertex& vertex ) {
    int rc( -1 );
    if ( vertex.isPrimary() )
      rc = 0;
    else if ( vertex.mother() ) {
      const LHCb::MCParticle* mother = vertex.mother();
      if ( mother->particleID().hasBottom() && ( mother->particleID().isMeson() || mother->particleID().isBaryon() ) )
        rc = 3;
      else if ( mother->particleID().hasCharm() &&
                ( mother->particleID().isMeson() || mother->particleID().isBaryon() ) )
        rc = 2;
      else if ( mother->particleID().hasStrange() )
        rc = 1;
      else
        rc = 4;
    }
    return rc;
  }

  struct ForwardMCHitSorter {
    bool operator()( const LHCb::MCHit* lhs, const LHCb::MCHit* rhs ) { return lhs->entry().z() < rhs->entry().z(); }
  };
  struct BackwardMCHitSorter {
    bool operator()( const LHCb::MCHit* lhs, const LHCb::MCHit* rhs ) { return rhs->entry().z() < lhs->entry().z(); }
  };

  void fillTuple( Tuples::Tuple theTuple, const LHCb::Track::Range& tracks, const LHCb::MCParticles& mcparticles,
                  const LHCb::MCHeader& mcheader, const LHCb::LinksByKey& linker, const Vertices& pvs,
                  const LHCb::MCHits* hits ) {
    // create the list of true PVs. count how many reconstructed tracks in each PV.
    std::map<const LHCb::MCVertex*, int> truepvs;
    // create a map from all MCParticles to MCHits
    MCHitMap mchitmap;

    if ( hits ) {
      // first collect
      for ( const LHCb::MCHit* mchit : *hits ) {
        if ( mchit->mcParticle() ) { mchitmap[mchit->mcParticle()].push_back( mchit ); }
      }
      // now sort them
      for ( auto& [vtx, hits] : mchitmap ) {
        if ( vtx->momentum().Pz() > 0 )
          std::sort( hits.begin(), hits.end(), ForwardMCHitSorter() );
        else
          std::sort( hits.begin(), hits.end(), BackwardMCHitSorter() );
      }
    }

    int itrack = 0;

    for ( const LHCb::Track* track : tracks ) {
      // keep track of track multiplicity
      theTuple->column( "itrack", int( itrack ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      ++itrack;
      theTuple->column( "ntrack", int( tracks.size() ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "numtruePV", int( mcheader.numOfPrimaryVertices() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      if ( !track->hasVelo() ) continue;

      double x  = track->firstState().x();
      double y  = track->firstState().y();
      double z  = track->firstState().z();
      double tx = track->firstState().tx();
      double ty = track->firstState().ty();

      theTuple->column( "probChi2", track->probChi2() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "chi2", track->chi2() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "ndof", track->nDoF() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "type", static_cast<int>( track->type() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "x", x ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "y", y ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "z", z ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "tx", tx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "ty", ty ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "qop", track->firstState().qOverP() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "veloChi2", track->info( LHCb::Track::AdditionalInfo::FitVeloChi2, 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "veloNdof", track->info( LHCb::Track::AdditionalInfo::FitVeloNDoF, 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "TChi2", track->info( LHCb::Track::AdditionalInfo::FitTChi2, 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "TNdof", track->info( LHCb::Track::AdditionalInfo::FitTNDoF, 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "backward", track->isVeloBackward() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      const LHCb::State* stateAtFirstHit = track->stateAt( LHCb::State::Location::FirstMeasurement );
      theTuple->column( "firsthittx", double( stateAtFirstHit ? stateAtFirstHit->tx() : 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "firsthitty", double( stateAtFirstHit ? stateAtFirstHit->ty() : 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "firsthitx", double( stateAtFirstHit ? stateAtFirstHit->x() : 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "firsthity", double( stateAtFirstHit ? stateAtFirstHit->y() : 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "firsthitz", double( stateAtFirstHit ? stateAtFirstHit->z() : 0 ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      LHCb::HitPattern hitpattern( track->lhcbIDs() );
      theTuple->column( "numVeloStations", int( hitpattern.numVeloStations() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "numVeloStationsOverlap", int( hitpattern.numVeloStationsOverlap() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "numVeloHoles", int( hitpattern.numVeloHoles() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "numTLayers", int( hitpattern.numTLayers() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "numVeloA", int( hitpattern.numVeloA() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "numVeloC", int( hitpattern.numVeloC() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->column( "numhits", int( track->lhcbIDs().size() ) )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      const LHCb::State&                state = track->firstState();
      LHCb::CubicStateInterpolationTraj tracktraj( state, Gaudi::XYZVector() );
      Gaudi::XYZPoint                   trkpos( state.position() );
      Gaudi::XYZVector                  trkdir( state.slopes().Unit() );

      auto fit            = fitResult( *track );
      int  trackWasFitted = fit && !fit->nodes().empty();
      theTuple->column( "trackWasFitted", trackWasFitted ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      // We also want to monitor the reconstructed IP, so the IP with respect to the reconstructed PVs.
      theTuple->column( "numrecPV", int( pvs.size() ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      const LHCb::VertexBase* recpv( nullptr );
      double                  bestip2( 0 );

      for ( const auto& thispv : pvs ) {
        Gaudi::XYZVector dx    = trkpos - thispv.position();
        Gaudi::XYZVector delta = dx - trkdir * dx.Dot( trkdir );
        double           ip2   = delta.Mag2();
        if ( recpv == 0 || ip2 < bestip2 ) {
          bestip2 = ip2;
          recpv   = &thispv;
        }
      }

      theTuple->column( "cov00", state.covariance()( 0, 0 ) ).ignore();
      theTuple->column( "cov11", state.covariance()( 1, 1 ) ).ignore();
      theTuple->column( "cov22", state.covariance()( 2, 2 ) ).ignore();
      theTuple->column( "cov33", state.covariance()( 3, 3 ) ).ignore();
      theTuple->column( "cov44", state.covariance()( 4, 4 ) ).ignore();

      theTuple->column( "cov01", state.covariance()( 0, 1 ) ).ignore();
      theTuple->column( "cov02", state.covariance()( 0, 2 ) ).ignore();
      theTuple->column( "cov03", state.covariance()( 0, 3 ) ).ignore();
      theTuple->column( "cov04", state.covariance()( 0, 4 ) ).ignore();

      theTuple->column( "cov12", state.covariance()( 2, 1 ) ).ignore();
      theTuple->column( "cov13", state.covariance()( 1, 3 ) ).ignore();
      theTuple->column( "cov14", state.covariance()( 1, 4 ) ).ignore();

      theTuple->column( "cov23", state.covariance()( 2, 3 ) ).ignore();
      theTuple->column( "cov24", state.covariance()( 2, 4 ) ).ignore();

      theTuple->column( "cov34", state.covariance()( 4, 3 ) ).ignore();

      if ( recpv ) {
        LHCb::State stateAtVtx = tracktraj.state( recpv->position().z() );
        // now compute the errors. this isn't quite right because:
        // - PV is biased
        // - correction for Z error is approximate
        double tx        = stateAtVtx.tx();
        double recipxerr = std::sqrt( state.covariance()( 0, 0 ) + recpv->covMatrix()( 0, 0 ) +
                                      2 * tx * recpv->covMatrix()( 0, 2 ) + tx * tx * recpv->covMatrix()( 2, 2 ) );
        double ty        = stateAtVtx.ty();
        double recipyerr = std::sqrt( state.covariance()( 1, 1 ) + recpv->covMatrix()( 1, 1 ) +
                                      2 * ty * recpv->covMatrix()( 1, 2 ) + ty * ty * recpv->covMatrix()( 2, 2 ) );

        theTuple->column( "recIP3D", std::sqrt( bestip2 ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recIPx", stateAtVtx.x() - recpv->position().x() )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recIPy", stateAtVtx.y() - recpv->position().y() )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recPVx", recpv->position().x() ).ignore();
        theTuple->column( "recPVy", recpv->position().y() ).ignore();
        theTuple->column( "recPVz", recpv->position().z() ).ignore();
        theTuple->column( "stateAtVertexX", stateAtVtx.position().x() ).ignore();
        theTuple->column( "stateAtVertexY", stateAtVtx.position().y() ).ignore();
        theTuple->column( "stateAtVertexZ", stateAtVtx.position().z() ).ignore();
        theTuple->column( "stateAtVertexTx", stateAtVtx.tx() ).ignore();
        theTuple->column( "stateAtVertexTy", stateAtVtx.ty() ).ignore();

        theTuple->column( "recpv_c00", recpv->covMatrix()( 0, 0 ) ).ignore();
        theTuple->column( "recpv_c01", recpv->covMatrix()( 0, 1 ) ).ignore();
        theTuple->column( "recpv_c02", recpv->covMatrix()( 0, 2 ) ).ignore();
        theTuple->column( "recpv_c03", recpv->covMatrix()( 0, 3 ) ).ignore();

        theTuple->column( "recpv_c11", recpv->covMatrix()( 1, 1 ) ).ignore();
        theTuple->column( "recpv_c12", recpv->covMatrix()( 1, 2 ) ).ignore();
        theTuple->column( "recpv_c13", recpv->covMatrix()( 1, 3 ) ).ignore();

        theTuple->column( "recpv_c22", recpv->covMatrix()( 2, 2 ) ).ignore();
        theTuple->column( "recpv_c23", recpv->covMatrix()( 2, 3 ) ).ignore();

        theTuple->column( "recpv_c33", recpv->covMatrix()( 3, 3 ) ).ignore();

        theTuple->column( "stateAtVertex_c00", stateAtVtx.covariance()( 0, 0 ) ).ignore();
        theTuple->column( "stateAtVertex_c01", stateAtVtx.covariance()( 0, 1 ) ).ignore();
        theTuple->column( "stateAtVertex_c02", stateAtVtx.covariance()( 0, 2 ) ).ignore();
        theTuple->column( "stateAtVertex_c03", stateAtVtx.covariance()( 0, 3 ) ).ignore();

        theTuple->column( "stateAtVertex_c11", stateAtVtx.covariance()( 1, 1 ) ).ignore();
        theTuple->column( "stateAtVertex_c12", stateAtVtx.covariance()( 1, 2 ) ).ignore();
        theTuple->column( "stateAtVertex_c13", stateAtVtx.covariance()( 1, 3 ) ).ignore();

        theTuple->column( "stateAtVertex_c22", stateAtVtx.covariance()( 2, 2 ) ).ignore();
        theTuple->column( "stateAtVertex_c23", stateAtVtx.covariance()( 2, 3 ) ).ignore();

        theTuple->column( "stateAtVertex_c33", stateAtVtx.covariance()( 3, 3 ) ).ignore();

        theTuple->column( "recIPxerr", recipxerr ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recIPyerr", recipyerr ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

        auto bestipchi2 = LHCb::TrackVertexUtils::vertexChi2( stateAtVtx, recpv->position(), recpv->covMatrix() );
        theTuple->column( "recIPChi2", bestipchi2 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      } else {
        theTuple->column( "recIP3D", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recIPx", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recIPy", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recIPxerr", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recIPyerr", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recIPChi2", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "recPVx", -99999999999999. ).ignore();
        theTuple->column( "recPVy", -99999999999999. ).ignore();
        theTuple->column( "recPVz", -99999999999999. ).ignore();

        theTuple->column( "stateAtVertexX", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertexY", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertexZ", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertexTx", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertexTy", -99999999999999. ).ignore();

        theTuple->column( "recpv_c00", -99999999999999. ).ignore();
        theTuple->column( "recpv_c01", -99999999999999. ).ignore();
        theTuple->column( "recpv_c02", -99999999999999. ).ignore();
        theTuple->column( "recpv_c03", -99999999999999. ).ignore();

        theTuple->column( "recpv_c11", -99999999999999. ).ignore();
        theTuple->column( "recpv_c12", -99999999999999. ).ignore();
        theTuple->column( "recpv_c13", -99999999999999. ).ignore();

        theTuple->column( "recpv_c22", -99999999999999. ).ignore();
        theTuple->column( "recpv_c23", -99999999999999. ).ignore();

        theTuple->column( "recpv_c33", -99999999999999. ).ignore();

        theTuple->column( "stateAtVertex_c00", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertex_c01", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertex_c02", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertex_c03", -99999999999999. ).ignore();

        theTuple->column( "stateAtVertex_c11", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertex_c12", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertex_c13", -99999999999999. ).ignore();

        theTuple->column( "stateAtVertex_c22", -99999999999999. ).ignore();
        theTuple->column( "stateAtVertex_c23", -99999999999999. ).ignore();

        theTuple->column( "stateAtVertex_c33", -99999999999999. ).ignore();

      } // end of recpv check

      // now do things linked to MC truth
      const LHCb::MCParticle* mcparticle{nullptr};
      double                  maxWeight{0};
      linker.applyToLinks(
          track->key(), [&maxWeight, &mcparticle, &mcparticles]( unsigned int, unsigned int mcPartKey, float weight ) {
            if ( weight > maxWeight ) {
              maxWeight  = weight;
              mcparticle = static_cast<const LHCb::MCParticle*>( mcparticles.containedObject( mcPartKey ) );
            }
          } );

      bool hasMCMatch = mcparticle && mcparticle->originVertex();
      theTuple->column( "waslinked", hasMCMatch ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      // The "typeofprefix" tag is used to see which particles come from a True Primary Vertex,
      // not a True Primary vertex or whether they are ghosts.
      bool   isFromPV = hasMCMatch && mcparticle->originVertex()->isPrimary();
      double typeofprefix;
      typeofprefix = isFromPV ? 0 : ( hasMCMatch ? 1 : 2 );
      theTuple->column( "typeofprefix", typeofprefix ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      bool hits_ctrl            = false;
      bool mchitmapit_ctrl      = false;
      bool stateAtFirstHit_ctrl = false;

      int mcOVT = -1; // set to -1 for ghosts
      if ( hasMCMatch ) {
        Gaudi::XYZPoint         trueorigin    = mcparticle->originVertex()->position();
        const LHCb::MCParticle* mother        = mcparticle->originVertex()->mother();
        const LHCb::MCVertex*   UltOrigVertex = findMCOriginVertex( *mcparticle, 1e-3 );
        const LHCb::MCParticle* ulmother      = UltOrigVertex ? UltOrigVertex->mother() : 0;
        mcOVT                                 = mcVertexType( UltOrigVertex );
        theTuple->column( "truepid", mcparticle->particleID().pid() )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        bool hasBottom = ulmother && ulmother->particleID().hasBottom() &&
                         ( ulmother->particleID().isMeson() || ulmother->particleID().isBaryon() );
        theTuple->column( "ulmotherHasBottom", hasBottom ? 1 : 0 )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        bool hasCharm = ulmother && ulmother->particleID().hasCharm() &&
                        ( ulmother->particleID().isMeson() || ulmother->particleID().isBaryon() );
        theTuple->column( "ulmotherHasCharm", hasCharm ? 1 : 0 )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "ulmotherHasStrange", ulmother && ulmother->particleID().hasStrange() ? 1 : 0 )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truemotherpid", mother ? mother->particleID().pid() : 0 )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truevertexpid", ulmother ? ulmother->particleID().pid() : 0 )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truex", trueorigin.x() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truey", trueorigin.y() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truez", trueorigin.z() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truemom", mcparticle->p() / Gaudi::Units::GeV )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truepx", mcparticle->momentum().Px() / Gaudi::Units::GeV )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truepy", mcparticle->momentum().Py() / Gaudi::Units::GeV )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truepz", mcparticle->momentum().Pz() / Gaudi::Units::GeV )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truept", mcparticle->momentum().Pt() / Gaudi::Units::GeV )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "trueeta", mcparticle->momentum().eta() )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truephi", mcparticle->momentum().phi() )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

        // we can use the CubicStateInterpolationTraj (or later tracktraj)
        // to do the error propagation. it is just a line!
        LHCb::State state = tracktraj.state( trueorigin.z() );
        double      tx    = state.tx();
        double      ty    = state.ty();
        double      IPx   = state.x() - trueorigin.x();
        double      IPy   = state.y() - trueorigin.y();
        double      IP3D  = std::sqrt( ( IPx * IPx + IPy * IPy ) / ( 1 + tx * tx + ty * ty ) );
        theTuple->column( "IPx", IPx ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPy", IPy ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IP3D", IP3D ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPxerr", std::sqrt( state.covariance()( 0, 0 ) ) )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPyerr", std::sqrt( state.covariance()( 1, 1 ) ) )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

        if ( hits ) {
          hits_ctrl                           = true;
          MCHitMap::const_iterator mchitmapit = mchitmap.find( mcparticle );
          if ( mchitmapit == mchitmap.end() ) {
            theTuple->column( "nummchits", int( 0 ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
          } else {
            mchitmapit_ctrl = true;
            // store the z-position of the first MC hit
            const MCHitVector& mchits = mchitmapit->second;
            theTuple->column( "nummchits", int( mchits.size() ) )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            // std::cout<< "mchits: " << mchits.size() << std::endl ;

            const LHCb::MCHit* mchit  = mchits.front();
            const LHCb::MCHit* mchitL = mchits.back();
            Gaudi::XYZPoint    poshit = mchit->entry();
            theTuple->column( "edep", mchit->energy() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "firstmchitx", poshit.x() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "firstmchity", poshit.y() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "firstmchitz", poshit.z() ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "firstmchitdz", mchit->displacement().z() )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "lastmchitz", mchitL->entry().z() )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "truetxfirstmchit", mchit->dxdz() )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "truetyfirstmchit", mchit->dydz() )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            double dz = poshit.z() - z;
            theTuple->column( "IPxfirstmchit", ( x + dz * tx ) - poshit.x() )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "IPyfirstmchit", ( y + dz * ty ) - poshit.y() )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

            if ( stateAtFirstHit ) {
              stateAtFirstHit_ctrl = true;
              // locate closest MCHit
              const LHCb::MCHit* closestmchit = mchit;
              for ( const LHCb::MCHit* anmchit : mchits ) {
                if ( std::abs( stateAtFirstHit->z() - anmchit->entry().z() ) <
                     std::abs( stateAtFirstHit->z() - mchit->entry().z() ) )
                  mchit = anmchit;
              }
              theTuple->column( "truetxfirsthit", closestmchit->dxdz() )
                  .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
              theTuple->column( "truetyfirsthit", closestmchit->dydz() )
                  .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
              Gaudi::XYZPoint posmchit = closestmchit->entry();
              double          dz       = posmchit.z() - stateAtFirstHit->z();
              theTuple->column( "IPxfirsthit", ( stateAtFirstHit->x() + dz * stateAtFirstHit->tx() ) - posmchit.x() )
                  .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
              theTuple->column( "IPyfirsthit", ( stateAtFirstHit->y() + dz * stateAtFirstHit->ty() ) - posmchit.y() )
                  .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            }

            // let's now extrapolate the mchit of the first hit to the z position of the vertex,
            // as if there were no scattering
            dz                        = trueorigin.z() - poshit.z();
            double extrapolatedmchitx = poshit.x() + dz * mchit->dxdz();
            double extrapolatedmchity = poshit.y() + dz * mchit->dydz();

            dz = trueorigin.z() - state.z();
            theTuple->column( "IPxfirsthitatvertex", ( state.x() + dz * state.tx() ) - extrapolatedmchitx )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
            theTuple->column( "IPyfirsthitatvertex", ( state.y() + dz * state.ty() ) - extrapolatedmchity )
                .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
          } // numhits check
        }   // mchits exist check
      } else {
        theTuple->column( "truepid", -999999 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "ulmotherHasBottom", 0 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "ulmotherHasCharm", 0 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "ulmotherHasStrange", 0 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truemotherpid", -999999 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truevertexpid", 0 ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truex", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truey", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truez", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truemom", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truepx", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truepy", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truepz", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truept", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "trueeta", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truephi", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPx", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPy", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IP3D", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPxerr", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPyerr", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      } // end of MCMatch check

      if ( !hits_ctrl ) {
        theTuple->column( "nummchits", int( 0 ) ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      }

      if ( !mchitmapit_ctrl ) {
        theTuple->column( "edep", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "firstmchitx", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "firstmchity", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "firstmchitz", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "firstmchitdz", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "lastmchitz", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truetxfirstmchit", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truetyfirstmchit", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPxfirstmchit", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPyfirstmchit", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPxfirsthitatvertex", -9999999999999. )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPyfirsthitatvertex", -9999999999999. )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      }

      if ( !stateAtFirstHit_ctrl ) {
        theTuple->column( "truetxfirsthit", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "truetyfirsthit", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPxfirsthit", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        theTuple->column( "IPyfirsthit", -9999999999999. ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      }

      theTuple->column( "mcOriginVertexType", mcOVT ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      theTuple->write().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
    }
  }
} // namespace

class TrackIPResolutionCheckerNT final
    : public LHCb::Algorithm::Consumer<void( const LHCb::Track::Range&, const LHCb::MCParticles&, const LHCb::MCHeader&,
                                             const LHCb::LinksByKey&, const Vertices& ),
                                       Gaudi::Functional::Traits::BaseClass_t<GaudiTupleAlg>> {
public:
  TrackIPResolutionCheckerNT( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"TrackContainer", LHCb::TrackLocation::Velo},
                   KeyValue{"MCParticleInput", LHCb::MCParticleLocation::Default},
                   KeyValue{"MCHeaderLocation", LHCb::MCHeaderLocation::Default},
                   KeyValue{"LinkerLocation", "Link/Pr/LHCbID"},
                   KeyValue{"PVContainer", LHCb::Event::PV::DefaultLocation}} ) {}
  void operator()( const LHCb::Track::Range& tracks, const LHCb::MCParticles& mcparticles,
                   const LHCb::MCHeader& mcheader, const LHCb::LinksByKey& linker,
                   const Vertices& pvs ) const override {
    return fillTuple( nTuple( "tracks", "", CLID_ColumnWiseTuple ), tracks, mcparticles, mcheader, linker, pvs,
                      nullptr );
  }
};
DECLARE_COMPONENT( TrackIPResolutionCheckerNT )

class TrackIPResolutionCheckerNTMCHits final
    : public LHCb::Algorithm::Consumer<void( const LHCb::Track::Range&, const LHCb::MCParticles&, const LHCb::MCHeader&,
                                             const LHCb::LinksByKey&, const Vertices&, const LHCb::MCHits& ),
                                       Gaudi::Functional::Traits::BaseClass_t<GaudiTupleAlg>> {
public:
  TrackIPResolutionCheckerNTMCHits( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {
                      KeyValue{"TrackContainer", LHCb::TrackLocation::Velo},
                      KeyValue{"MCParticleInput", LHCb::MCParticleLocation::Default},
                      KeyValue{"MCHeaderLocation", LHCb::MCHeaderLocation::Default},
                      KeyValue{"LinkerLocation", "Link/Pr/LHCbID"},
                      KeyValue{"PVContainer", LHCb::Event::PV::DefaultLocation},
                      KeyValue{"MCHitsLocation", "/Event/MC/VP/Hits"},
                  } ) {}

  void operator()( const LHCb::Track::Range& tracks, const LHCb::MCParticles& mcparticles,
                   const LHCb::MCHeader& mcheader, const LHCb::LinksByKey& linker, const Vertices& pvs,
                   const LHCb::MCHits& mchits ) const override {
    return fillTuple( nTuple( "tracks", "", CLID_ColumnWiseTuple ), tracks, mcparticles, mcheader, linker, pvs,
                      &mchits );
  }
};
DECLARE_COMPONENT( TrackIPResolutionCheckerNTMCHits )
