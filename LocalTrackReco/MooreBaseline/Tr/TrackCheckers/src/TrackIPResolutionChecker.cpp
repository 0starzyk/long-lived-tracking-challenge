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
#include "Event/MCHeader.h"
#include "Event/MCParticle.h"
#include "Event/RecVertex.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "Event/TrackFitResult.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiUtils/Aida2ROOT.h"
#include "LHCbMath/ValueWithError.h"
#include "Linker/LinkedTo.h"
#include "TH1.h"
#include "TH2.h"
#include "TrackKernel/CubicStateInterpolationTraj.h"
#include "TrackKernel/TrackFunctors.h"

/** @class TrackIPResolutionChecker TrackIPResolutionChecker.cpp "TrackCheckers/TrackIPResolutionChecker"
 *
 * Class for track monitoring
 *  @author W. Hulsbergen and P. Tsopelas
 *  @date   21-01-2013
 */

class TrackIPResolutionChecker : public GaudiHistoAlg {

public:
  /** Standard construtor */
  using GaudiHistoAlg::GaudiHistoAlg;

  /** Algorithm execute */
  StatusCode execute() override;

  /** Algorithm finalize */
  StatusCode finalize() override;

private:
  void createResolutionProfile( const HistoID& inputname, const HistoID& outputname );

private:
  Gaudi::Property<std::string> m_tracklocation{this, "TrackContainer",
                                               LHCb::TrackLocation::Velo}; ///< Input Tracks container location
  Gaudi::Property<std::string> m_linkerlocation{this, "LinkerLocation", LHCb::TrackLocation::Velo};
  Gaudi::Property<std::string> m_pvlocation{this, "PVContainer", LHCb::RecVertexLocation::Primary};
  Gaudi::Property<int>         m_minNumTracksReconstructablePV{this, "MinNumTracksReconstructablePV", 5};
};

DECLARE_COMPONENT( TrackIPResolutionChecker )

namespace {
  Gaudi::Math::ValueWithError calc3sigmarms( const TH1& h1 ) {
    // finally, the 3sigma RMS. that one is most hard. we want the limits
    // that are 3*RMS such that the RMS within these limits is
    // consistent.

    // lets first map the histogram into something more symmetric. the
    // lazy way.
    const TH1* h1pos = &h1;
    TH1F*      ownedh1pos( 0 );
    if ( h1.GetXaxis()->GetXmin() < 0 ) {
      // create a new histogram that has exactly half the bins and just fill it with abs(x)
      ownedh1pos = new TH1F( "h1pos", "", h1.GetNbinsX() / 2, 0, h1.GetXaxis()->GetXmax() );
      for ( int ibin = 1; ibin <= h1.GetNbinsX(); ++ibin )
        ownedh1pos->Fill( std::abs( h1.GetBinCenter( ibin ) ), h1.GetBinContent( ibin ) );
      h1pos = ownedh1pos;
    }

    // now just start counting
    double sumw( 0 ), sumx2w( 0 );
    double xtrunc( 0 );
    int    ibin = 1;
    for ( ; ibin <= h1pos->GetNbinsX(); ++ibin ) {
      double x  = h1pos->GetBinCenter( ibin );
      double c  = h1pos->GetBinContent( ibin );
      double up = h1pos->GetXaxis()->GetBinUpEdge( ibin );

      double newsumw   = sumw + c;
      double newsumx2w = sumx2w + c * x * x;
      if ( sumw > 0 && x > h1pos->GetMean() ) {
        double newrms = sqrt( newsumx2w / newsumw );
        if ( 3 * newrms < up ) {
          double drms = newrms - sqrt( sumx2w / sumw );
          double frac = ( 3 * drms ) / h1pos->GetXaxis()->GetBinWidth( ibin );
          // std::cout << frac << std::endl ;
          xtrunc = ( up - ( 1 - frac ) * h1pos->GetXaxis()->GetBinWidth( ibin ) ) / 3;
          break;
        }
      }
      sumw   = newsumw;
      sumx2w = newsumx2w;
    }
    if ( ibin > h1pos->GetNbinsX() && sumw > 0 ) // just return the rms
      xtrunc = sqrt( sumx2w / sumw );
    if ( ownedh1pos ) delete ownedh1pos;
    return sumw > 0 ? Gaudi::Math::ValueWithError( xtrunc, xtrunc * xtrunc / sumw )
                    : Gaudi::Math::ValueWithError( 0, 0 );
  }

  // count the number of charged stable daughters in the LHCb acceptance for this particle
  size_t numChargedDaughtersInAcceptance( const LHCb::MCVertex& vertex ) {
    size_t rc( 0 );
    for ( const LHCb::MCParticle* particle : vertex.products() ) {
      // if 'resonance' then just add it
      const LHCb::MCVertex* endvertex = particle->endVertices().size() > 0 ? &( *( particle->endVertices()[0] ) ) : 0;
      if ( endvertex && std::abs( endvertex->position().z() - vertex.position().z() ) < 0.001 * Gaudi::Units::mm )
        rc += numChargedDaughtersInAcceptance( *endvertex );
      // add charged daughters that live long enough
      else if (                                     // require that it is charged
          particle->particleID().threeCharge() != 0 // charged
          // require a minimum momentum to traverse the velo
          && particle->p() > 300 * Gaudi::Units::MeV
          // require that its eta is within the acceptance
          && 2 < particle->pseudoRapidity() &&
          particle->pseudoRapidity() < 5
          // require that there is no endvertex, or that it is sort-of outside the velo
          && ( endvertex == 0 || ( endvertex->position() - vertex.position() ).Rho() > 10 * Gaudi::Units::cm ) ) {
        ++rc;
      }
    }
    return rc;
  }

} // namespace

// routine that creates a resolution profile from a 2D histogram
void TrackIPResolutionChecker::createResolutionProfile( const HistoID& inputname, const HistoID& outputname ) {
  // get the histogram
  const TH2* h2 = Gaudi::Utils::Aida2ROOT::aida2root( histo2D( inputname ) );
  // create the 1D histogram
  if ( !h2 ) {
    debug() << "Cannot find histogram with name: " << inputname << endmsg;
  } else {
    std::string title = h2->GetTitle();
    size_t      ipos  = title.find( " versus" );
    if ( ipos == std::string::npos ) ipos = title.find( " vs" );
    if ( ipos != std::string::npos ) title.insert( ipos, " resolution" );

    TH1* h1 = Gaudi::Utils::Aida2ROOT::aida2root(
        book1D( outputname, title.c_str(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(), h2->GetNbinsX() ) );

    // info() << outputname << " = " << title.c_str() << " xmin: " << h2->GetXaxis()->GetXmin() << " xmax: " <<
    // h2->GetXaxis()->GetXmax() << " nbins: " << h2->GetNbinsX() << endmsg;

    // now fill it
    for ( int ibin = 1; ibin <= h2->GetNbinsX(); ++ibin ) {
      TH1* h1tmp = h2->ProjectionY( "tmp", ibin, ibin );
      if ( h1tmp->Integral() > 5 ) {
        Gaudi::Math::ValueWithError sigma = calc3sigmarms( *h1tmp );
        h1->SetBinContent( ibin, sigma.value() );
        h1->SetBinError( ibin, sigma.error() );
        // info() << sigma.value() << " " << sigma.error() << endmsg;
      }
      delete h1tmp;
    }
  }
}

//=============================================================================
// Initialization. Check parameters
//=============================================================================
StatusCode TrackIPResolutionChecker::finalize() {
  // fill some resolution profiles
  const char names[3][256] = {"IP/Velo/", "IP/Long/", "IP/Backward/"};
  for ( int i = 0; i < 3; ++i ) {
    std::string prefix( names[i] );
    createResolutionProfile( prefix + "IP3DVsInvTruePtH2", prefix + "IP3DResolutionVsInvTruePt" );
    createResolutionProfile( prefix + "IPxVsInvTruePtH2", prefix + "IPxResolutionVsInvTruePt" );
    createResolutionProfile( prefix + "IPyVsInvTruePtH2", prefix + "IPyResolutionVsInvTruePt" );
  }
  createResolutionProfile( "PV/dxVersusNTrk", "PV/PVXResolutionVsNTrk" );
  createResolutionProfile( "PV/dyVersusNTrk", "PV/PVYResolutionVsNTrk" );
  createResolutionProfile( "PV/dzVersusNTrk", "PV/PVZResolutionVsNTrk" );

  // Efficiencies profiles
  const TH1* h1_ip3d_b = Gaudi::Utils::Aida2ROOT::aida2root( histo1D( ( HistoID ) "Eff/IP3D_b" ) );
  const TH1* h1_ip3d   = Gaudi::Utils::Aida2ROOT::aida2root( histo1D( ( HistoID ) "Eff/IP3D" ) );

  int nbins = h1_ip3d->GetNbinsX();
  for ( int ibin = 1; ibin <= nbins; ++ibin ) {
    double eff_b = h1_ip3d_b->Integral( ibin, nbins ) / h1_ip3d_b->Integral( 1, nbins );
    profile1D( ibin * 0.005 - 0.0025, eff_b, "Eff/eff_b", "e(IP3D_b > x)", 0, 0.5 );

    double eff = h1_ip3d->Integral( ibin, nbins ) / h1_ip3d->Integral( 1, nbins );
    profile1D( ibin * 0.005 - 0.0025, eff, "Eff/eff", "e(IP3D > x)", 0, 0.5 );
  }

  return GaudiHistoAlg::finalize();
}

bool hasB( const LHCb::MCParticle* mcp ) {
  for ( ; mcp; mcp = mcp->originVertex()->mother() ) {
    if ( mcp->particleID().hasBottom() && ( mcp->particleID().isMeson() || mcp->particleID().isBaryon() ) ) return true;
  }
  return false;
}

bool hasD( const LHCb::MCParticle* mcp ) {
  for ( ; mcp; mcp = mcp->originVertex()->mother() ) {
    if ( mcp->particleID().hasCharm() && ( mcp->particleID().isMeson() || mcp->particleID().isBaryon() ) ) return true;
  }
  return false;
}

//=============================================================================
// Execute
//=============================================================================
StatusCode TrackIPResolutionChecker::execute() {
  // some constant that should probably become properties
  const double maxinvpt   = 3.0;
  const double maxip      = 0.5;
  const int    nbinsinvpt = 15;

  // get the list of tracks
  const LHCb::Track::Range tracks = get<LHCb::Track::Range>( m_tracklocation );

  // get the list of PVs as well
  const LHCb::RecVertex::Range pvs = get<LHCb::RecVertex::Range>( m_pvlocation );

  // get the linker table track -> mcparticle
  const auto* links = get<LHCb::LinksByKey>( LHCb::LinksByKey::linkerName( m_linkerlocation ) );

  // efficiency check
  for ( const LHCb::Track* track : tracks ) {
    const LHCb::MCParticle* mcparticle = LinkedTo<LHCb::MCParticle>{links}.range( track ).try_front();
    if ( mcparticle && mcparticle->originVertex() ) {
      // make some selection for the IP study
      if ( mcparticle->momentum().Pz() > 0 && 2 < mcparticle->pseudoRapidity() && mcparticle->pseudoRapidity() < 5 ) {

        std::string prefix = std::string( "Eff/" );

        double trueptGeV = mcparticle->momentum().Pt() / Gaudi::Units::GeV;
        if ( trueptGeV > 0.2 ) {
          // Gaudi::XYZPoint trueorigin = mcparticle->originVertex()->position();
          Gaudi::XYZPoint trueorigin = mcparticle->primaryVertex()->position();

          // for now, the track state is simply firststate
          const LHCb::State& firststate = track->firstState();
          // we can use the CubicStateInterpolationTraj (or later tracktraj) to do the error propagation
          LHCb::CubicStateInterpolationTraj tracktraj( firststate, Gaudi::XYZVector() );
          LHCb::State                       state = tracktraj.state( trueorigin.z() );
          double                            tx    = state.tx();
          double                            ty    = state.ty();
          double                            IPx   = state.x() - trueorigin.x();
          double                            IPy   = state.y() - trueorigin.y();
          double                            IP3D  = std::sqrt( ( IPx * IPx + IPy * IPy ) / ( 1 + tx * tx + ty * ty ) );

          if ( hasB( mcparticle ) )
            plot1D( IP3D, prefix + "IP3D_b", "IP 3D from B", 0, maxip );
          else if ( hasD( mcparticle ) )
            plot1D( IP3D, prefix + "IP3D_d", "IP 3D from D", 0, maxip );
          else
            plot1D( IP3D, prefix + "IP3D", "IP 3D background", 0, maxip );
        }
      }
    }
  }

  // loop over the tracks
  for ( const LHCb::Track* track : tracks ) {
    const LHCb::MCParticle* mcparticle = LinkedTo<LHCb::MCParticle>{links}.range( track ).try_front();
    if ( mcparticle && mcparticle->originVertex() ) {
      // make some selection for the IP study
      if ( mcparticle->momentum().Pz() > 0 && 2 < mcparticle->pseudoRapidity() && mcparticle->pseudoRapidity() < 5 ) {

        std::string prefix = std::string( "IP/" ) + Gaudi::Utils::toString( track->type() ) + "/";

        double trueptGeV = mcparticle->momentum().Pt() / Gaudi::Units::GeV;
        plot1D( trueptGeV, prefix + "TruePtH1", "true Pt in GeV", 0, 5 );
        plot1D( mcparticle->p() / Gaudi::Units::GeV, prefix + "TruePH1", "true momentum in GeV", 0, 100 );
        const auto* fit            = fitResult( *track );
        int         trackWasFitted = fit && !fit->nodes().empty();
        plot1D( trackWasFitted, prefix + "TrackWasFittedH1", "track was fitted by K-filter", -0.5, 1.5, 2 );

        if ( trueptGeV > 0.2 ) {
          double truephi      = mcparticle->momentum().Phi();
          double invtrueptGeV = 1 / trueptGeV;

          plot( track->lhcbIDs().size(), "NumLHCBIDs", "Number of LHCbIDs", -0.5, 40.5, 41 );
          // Gaudi::XYZPoint trueorigin = mcparticle->originVertex()->position();
          Gaudi::XYZPoint trueorigin = mcparticle->primaryVertex()->position();
          plot( trueorigin.z(), prefix + "TrueOriginZ", "True origin z", -100, 100 );
          plot( trueorigin.x(), prefix + "TrueOriginX", "True origin x", -0.5, 0.5 );
          plot( trueorigin.y(), prefix + "TrueOriginY", "True origin y", -0.5, 0.5 );

          // for now, the track state is simply firststate
          const LHCb::State& firststate = track->firstState();
          // we can use the CubicStateInterpolationTraj (or later tracktraj) to do the error propagation
          LHCb::CubicStateInterpolationTraj tracktraj( firststate, Gaudi::XYZVector() );
          LHCb::State                       state = tracktraj.state( trueorigin.z() );
          double                            tx    = state.tx();
          double                            ty    = state.ty();
          double                            IPx   = state.x() - trueorigin.x();
          double                            IPy   = state.y() - trueorigin.y();
          // double                            IP3D  = std::sqrt( ( IPx * IPx + IPy * IPy ) / ( 1 + tx * tx + ty * ty )
          // );
          double IP3D = std::sqrt( ( IPx * IPx + IPy * IPy + ( IPx * ty - IPy * tx ) * ( IPx * ty - IPy * tx ) ) /
                                   ( 1 + tx * tx + ty * ty ) );

          /*Vec3<float> A = Vec3<float>(trueorigin.x(), trueorigin.y(), trueorigin.z());
          Vec3<float> B = Vec3<float>(state.x(), state.y(), trueorigin.z());
          Vec3<float> u = Vec3<float>(state.tx(), state.ty(), 1.f);

          double                            IP3D = (B - A).cross(u).mag() / u.mag();*/

          plot2D( invtrueptGeV, IPx, prefix + "IPxVsInvTruePtH2", "IPx versus 1/pt_true", 0, maxinvpt, -maxip, maxip,
                  nbinsinvpt, 200 );
          plot2D( invtrueptGeV, IPy, prefix + "IPyVsInvTruePtH2", "IPy versus 1/pt_true", 0, maxinvpt, -maxip, maxip,
                  nbinsinvpt, 200 );
          plot2D( invtrueptGeV, IP3D, prefix + "IP3DVsInvTruePtH2", "IP versus 1/pt_true", 0, maxinvpt, 0, maxip,
                  nbinsinvpt, 200 );

          // True IP : True PV & True MCParticle
          Gaudi::XYZPoint mcpos    = mcparticle->originVertex()->position();
          double          trueIPx  = mcpos.x() - trueorigin.x();
          double          trueIPy  = mcpos.y() - trueorigin.y();
          double          trueIPz  = mcpos.z() - trueorigin.z();
          double          trueIP3D = std::sqrt( trueIPx * trueIPx + trueIPy * trueIPy + trueIPz * trueIPz );

          plot1D( trueIPx, prefix + "trueIPxH1", "True IP x", -maxip, maxip );
          plot1D( trueIPy, prefix + "trueIPyH1", "True IP y", -maxip, maxip );
          plot1D( trueIP3D, prefix + "trueIP3DH1", "True IP 3D", 0, maxip );

          plot1D( std::abs( IP3D - trueIP3D ), prefix + "IP3Derr", "abs(IP3D - IP3D_{True})", 0, 0.1 );

          auto   mom    = mcparticle->momentum();
          double trueTx = mom.x() / mom.z();
          double trueTy = mom.y() / mom.z();
          double maxErr = 0.006;
          plot1D( ( state.tx() - trueTx ), prefix + "errTxH1", "tx - trueTx", -maxErr, maxErr );
          plot1D( ( state.ty() - trueTy ), prefix + "errTyH1", "ty - trueTy", -maxErr, maxErr );

          if ( trueptGeV > 0.5 ) {
            const double pi = M_PI;
            plot2D( truephi, IPx, prefix + "IPxVsTruePhiH2", "IPx versus phi_{true}", -pi, pi, -maxip, maxip, 24, 100 );
            plot2D( truephi, IPy, prefix + "IPyVsTruePhiH2", "IPy versus phi_{true}", -pi, pi, -maxip, maxip, 24, 100 );
          }

          plot1D( IPx, prefix + "IPxH1", "IP x", -maxip, maxip );
          plot1D( IPy, prefix + "IPyH1", "IP y", -maxip, maxip );
          plot1D( IP3D, prefix + "IP3DH1", "IP 3D", 0, maxip );
          plot1D( IPx / std::sqrt( state.covariance()( 0, 0 ) ), prefix + "IPxPullH1", "IP x pull", -5, 5 );
          plot1D( IPy / std::sqrt( state.covariance()( 1, 1 ) ), prefix + "IPyPullH1", "IP y pull", -5, 5 );
        }
      }
    }
  }

  // keep track of which true PVs were found
  const LHCb::MCHeader*           mcheader = get<const LHCb::MCHeader*>( LHCb::MCHeaderLocation::Default );
  std::set<const LHCb::MCVertex*> foundpvs;
  for ( const LHCb::RecVertex* pv : pvs ) {
    // plot the number of tracks per PV
    plot1D( pv->tracks().size(), "PV/RecoPVNumTracksH1", "Number of reconstructed tracks in reconstructed PVs", 0.5,
            50.5, 50 );

    // find the closest PV, in absolute distance
    const LHCb::MCVertex* truepv( 0 );
    double                bestdistance2 = 0;
    for ( const LHCb::MCVertex* thistruepv : mcheader->primaryVertices() ) {
      double distance2 = ( pv->position() - thistruepv->position() ).Mag2();
      if ( truepv == 0 || distance2 < bestdistance2 ) {
        bestdistance2 = distance2;
        truepv        = thistruepv;
      }
    }
    if ( truepv ) {
      foundpvs.insert( truepv );
      Gaudi::XYZVector delta = pv->position() - truepv->position();
      plot1D( delta.x(), "PV/dxH1", "x_{PV} - x_{PV,true}", -0.08, 0.08 );
      plot1D( delta.y(), "PV/dyH1", "y_{PV} - y_{PV,true}", -0.08, 0.08 );
      plot1D( delta.z(), "PV/dzH1", "z_{PV} - z_{PV,true}", -0.4, 0.4 );
      plot2D( pv->tracks().size(), delta.x(), "PV/dxVersusNTrk", "x_{PV} - x_{PV,true} vs number of tracks", 5, 105,
              -0.08, 0.08, 20, 80 );
      plot2D( pv->tracks().size(), delta.y(), "PV/dyVersusNTrk", "y_{PV} - y_{PV,true} vs number of tracks", 5, 105,
              -0.08, 0.08, 20, 80 );
      plot2D( pv->tracks().size(), delta.z(), "PV/dzVersusNTrk", "z_{PV} - z_{PV,true} vs number of tracks", 5, 105,
              -0.4, 0.4, 20, 80 );
      plot1D( delta.x() / std::sqrt( pv->covMatrix()( 0, 0 ) ), "PV/dxpullH1", "x_{PV} pull", -5, 5 );
      plot1D( delta.y() / std::sqrt( pv->covMatrix()( 1, 1 ) ), "PV/dypullH1", "y_{PV} pull", -5, 5 );
      plot1D( delta.z() / std::sqrt( pv->covMatrix()( 2, 2 ) ), "PV/dzpullH1", "z_{PV} pull", -5, 5 );
    }
  }

  // plot the total number of true PVs
  size_t numtruepvs = mcheader->numOfPrimaryVertices();
  plot1D( numtruepvs, "PV/NumTruePVs", "Number of true PVs", -0.5, 20.5, 21 );

  // make a plot of PV efficiency versus Z and versus number of tracks
  size_t numreconstructablepvs( 0 );
  for ( const LHCb::MCVertex* pv : mcheader->primaryVertices() ) {
    size_t numreconstructable = numChargedDaughtersInAcceptance( *pv );
    bool   found              = foundpvs.find( pv ) != foundpvs.end();
    profile1D( double( numreconstructable ), found, "PV/PVEfficiencyVsNumReconstructableTracks",
               "PV efficiency versus num tracks in acceptance", 0.5, 50.5, 50 );
    if ( numreconstructable >= 5 ) {
      profile1D( pv->position().z(), found, "PV/PVEfficiencyVsZ", "PV efficiency versus Z", -100, 100, 40 );
      profile1D( numtruepvs, found, "PV/PVEfficiencyVsNumTruePVs", "PV efficiency versus num true PVs", 0.5, 20.5, 20 );
      ++numreconstructablepvs;
    }
  }
  plot2D( numreconstructablepvs, pvs.size(), "PV/NumPVsVsNumTrueReconstructablePVs",
          "Number of PVs vs number of true reconstructable PVs", -0.5, 20.5, -0.5, 20.5, 21, 21 );

  // Finally, we also want to monitor the reconstructed IP, so the IP with respect to the reconstructed PVs.
  if ( !pvs.empty() ) {
    for ( const LHCb::Track* track : tracks )
      if ( !track->isVeloBackward() ) {
        const LHCb::MCParticle* mcparticle = LinkedTo<LHCb::MCParticle>{links}.range( track ).try_front();
        // distinghuish secondaries, from primaries, from ghosts
        bool        hasMCMatch = mcparticle && mcparticle->originVertex();
        bool        isFromPV   = hasMCMatch && mcparticle->originVertex()->isPrimary();
        std::string prefix =
            isFromPV ? "IPRecPV/TruePrimary/" : ( hasMCMatch ? "IPRecPV/TrueSecondary/" : "IPRecPV/Ghost/" );

        // find the closest PV, again using the minimal distance
        const LHCb::RecVertex*            pv( 0 );
        double                            bestip2( 0 );
        const LHCb::State&                state = track->firstState();
        LHCb::CubicStateInterpolationTraj tracktraj( state, Gaudi::XYZVector() );
        Gaudi::XYZPoint                   trkpos( state.position() );
        Gaudi::XYZVector                  trkdir( state.slopes().Unit() );

        for ( const LHCb::RecVertex* thispv : pvs ) {
          Gaudi::XYZVector dx    = trkpos - thispv->position();
          Gaudi::XYZVector delta = dx - trkdir * dx.Dot( trkdir );
          double           ip2   = delta.Mag2();
          if ( pv == 0 || ip2 < bestip2 ) {
            bestip2 = ip2;
            pv      = thispv;
          }
        }

        plot1D( std::sqrt( bestip2 ), prefix + "IP3DH1", "IP 3D wrt to reconstructed vertex", 0, maxip );

        LHCb::State stateAtVtx = tracktraj.state( pv->position().z() );
        double      dx         = stateAtVtx.x() - pv->position().x();
        double      dy         = stateAtVtx.y() - pv->position().y();
        plot1D( dx, prefix + "IPxH1", "IP X wrt reconstructed vertex", -maxip, maxip );
        plot1D( dy, prefix + "IPyH1", "IP Y wrt reconstructed vertex", -maxip, maxip );

        // now compute the errors. this isn't quite right because
        // - PV is biased
        // - correction for Z error is approximate
        double tx      = stateAtVtx.tx();
        double sigmaX2 = state.covariance()( 0, 0 ) + pv->covMatrix()( 0, 0 ) + 2 * tx * pv->covMatrix()( 0, 2 ) +
                         tx * tx * pv->covMatrix()( 2, 2 );
        plot1D( dx / std::sqrt( sigmaX2 ), prefix + "IPxPullH1", "IP X pull wrt reconstructed vertex", -5, 5 );
        double ty      = stateAtVtx.ty();
        double sigmaY2 = state.covariance()( 1, 1 ) + pv->covMatrix()( 1, 1 ) + 2 * ty * pv->covMatrix()( 1, 2 ) +
                         ty * ty * pv->covMatrix()( 2, 2 );
        plot1D( dx / std::sqrt( sigmaY2 ), prefix + "IPyPullH1", "IP Y pull wrt reconstructed vertex", -5, 5 );
      }
  } // end of !pvs.empty()

  return StatusCode::SUCCESS;
}
