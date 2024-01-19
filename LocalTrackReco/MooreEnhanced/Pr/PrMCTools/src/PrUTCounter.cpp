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

#include "PrUTCounter.h"

#include "GaudiAlg/IHistoTool.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrUTCounter
//
// 2006-06-28 : Olivier Callot
// 2015-01-17 : modified by Svende Braun, Michel de Cian to work with PrChecker2.cpp
//-----------------------------------------------------------------------------

void PrUTCounter::addSelection( std::string name, bool writeHisto, bool ) {
  if ( name.size() + 5 > m_titleSize ) { m_titleSize = name.size() + 5; }
  m_name.push_back( name );
  m_writeHisto.push_back( writeHisto );
  m_nbTrack.push_back( 0. );
  m_mcHits.push_back( 0. );
  m_foundOK.push_back( 0. );
  m_wrong.push_back( 0. );

  m_nbTrack3.push_back( 0. );
  m_mcHits3.push_back( 0. );
  m_foundOK3.push_back( 0. );
  m_wrong3.push_back( 0. );
}

void PrUTCounter::initEvent( IHistoTool const* htool, ITrackExtrapolator const* extrapolator, int,
                             LHCb::Track::Range const& tracks, LHCb::LinksByKey const& links,
                             IGeometryInfo const& geometry ) {
  // Process the ghost tracks
  for ( const auto* tr : tracks ) {
    StatusCode sc = StatusCode::FAILURE;
    sc.ignore();

    // -- Needed for x/y plots at a given z
    LHCb::State state;
    LHCb::State state2;
    bool        doXYZPlots = m_xyPlots;
    if ( m_xyPlots ) {
      state         = tr->closestState( 9000 );
      StatusCode sc = extrapolator->propagate( state, 9000, geometry );
      if ( sc.isFailure() ) { doXYZPlots = false; }

      state2 = tr->closestState( 2485 );
      sc     = extrapolator->propagate( state2, 2485, geometry );
      if ( sc.isFailure() ) { doXYZPlots = false; }

      if ( std::isnan( state.x() ) || std::isnan( state.y() ) || std::isnan( state2.x() ) || std::isnan( state2.y() ) )
        doXYZPlots = false;
    }

    if ( !links.hasEntry( *tr ) ) {
      double nbInUT = 0;
      for ( std::vector<LHCb::LHCbID>::const_iterator itId = tr->lhcbIDs().begin(); tr->lhcbIDs().end() != itId;
            ++itId ) {
        if ( ( *itId ).isUT() ) { nbInUT += 1.; }
      }
      m_nbGhost += 1.;
      if ( htool && m_writeHistos > 0 ) {
        htool->plot1D( tr->pseudoRapidity(), m_title + "/Eta_Ghosts", "Eta_Ghosts", 0., 7., 50 );
        if ( tr->type() != LHCb::Track::Types::Velo ) {
          htool->plot1D( tr->momentum().Phi(), m_title + "/Phi_Ghosts", "Phi_Ghosts", -3.142, 3.142, 25 );
          htool->plot1D( tr->pt(), m_title + "/Pt_Ghosts", "Pt_Ghosts", 0., 10000., 100 );
          htool->plot1D( tr->p(), m_title + "/P_Ghosts", "P_Ghosts", 0., 100000., 100 );
          htool->plot2D( tr->pseudoRapidity(), tr->p(), m_title + "/EtaP_Ghosts", "EtaP_Ghosts", 0., 7., 0., 100000.,
                         20, 20 );
          htool->plot2D( tr->pseudoRapidity(), tr->momentum().Phi(), m_title + "/EtaPhi_Ghosts", "EtaPhi_Ghosts", 0.,
                         7., -3.142, 3.142, 20, 20 );
        }
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/XYZ9000_Ghosts", "XYZ9000_Ghosts", -3000, 3000., -3000.,
                         3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/XYZ2485_Ghosts", "XYZ2485_Ghosts", -1000, 1000., -1000.,
                         1000.0, 100, 100 );
        }
      }
      m_nbGhostHit += nbInUT;
    }
    if ( htool && m_writeHistos > 0 ) {
      htool->plot1D( tr->pseudoRapidity(), m_title + "/Eta_Total", "Eta_Total", 0., 7., 50 );
      if ( tr->type() != LHCb::Track::Types::Velo ) {
        htool->plot1D( tr->momentum().Phi(), m_title + "/Phi_Total", "Phi_Total", -3.142, 3.142, 25 );
        htool->plot1D( tr->pt(), m_title + "/Pt_Total", "Pt_Total", 0., 10000., 100 );
        htool->plot1D( tr->p(), m_title + "/P_Total", "P_Total", 0., 100000., 100 );
        htool->plot2D( tr->pseudoRapidity(), tr->p(), m_title + "/EtaP_Total", "EtaP_Total", 0., 7., 0., 100000., 20,
                       20 );
        htool->plot2D( tr->pseudoRapidity(), tr->momentum().Phi(), m_title + "/EtaPhi_Total", "EtaPhi_Total", 0., 7.,
                       -3.142, 3.142, 20, 20 );
      }
      if ( doXYZPlots ) {
        htool->plot2D( state.x(), state.y(), m_title + "/XYZ9000_Total", "XYZ9000_Total", -3000, 3000., -3000., 3000.0,
                       100, 100 );
        htool->plot2D( state2.x(), state2.y(), m_title + "/XYZ2485_Ghosts", "XYZ2485_Ghosts", -1000, 1000., -1000.,
                       1000.0, 100, 100 );
      }
    }
    m_totTrack++;
    if ( m_triggerNumbers && ( tr->type() != LHCb::Track::Types::Velo ) && tr->p() > 3000. && tr->pt() > 500. ) {
      if ( !links.hasEntry( *tr ) ) {
        m_totGhostTrigger++;
        if ( htool && m_writeHistos > 0 ) {
          htool->plot1D( tr->pseudoRapidity(), m_title + "/Eta_Ghosts_P>3GeV_Pt>0.5GeV", "Eta_Ghosts_P>3GeV_Pt>0.5GeV",
                         0., 7., 50 );
          htool->plot1D( tr->momentum().Phi(), m_title + "/Phi_Ghosts_P>3GeV_Pt>0.5GeV", "Phi_Ghosts_P>3GeV_Pt>0.5GeV",
                         -3.142, 3.142, 25 );
          htool->plot1D( tr->pt(), m_title + "/Pt_Ghosts_P>3GeV_Pt>0.5GeV", "Pt_Ghosts_P>3GeV_Pt>0.5GeV", 0., 10000.,
                         100 );
          htool->plot1D( tr->p(), m_title + "/P_Ghosts_P>3GeV_Pt>0.5GeV", "P_Ghosts_P>3GeV_Pt>0.5GeV", 0., 100000.,
                         100 );
          htool->plot2D( tr->pseudoRapidity(), tr->p(), m_title + "/EtaP_Ghosts_P>3GeV_Pt>0.5GeV",
                         "EtaP_Ghosts_P>3GeV_Pt>0.5GeV", 0., 7., 0., 100000., 20, 20 );
          htool->plot2D( tr->pseudoRapidity(), tr->momentum().Phi(), m_title + "/EtaPhi_Ghosts_P>3GeV_Pt>0.5GeV",
                         "EtaPhi_Ghosts_P>3GeV_Pt>0.5GeV", 0., 7., -3.142, 3.142, 20, 20 );
          if ( doXYZPlots ) {
            htool->plot2D( state.x(), state.y(), m_title + "/XYZ9000_Ghosts_P>3GeV_Pt>0.5GeV",
                           "XYEff_Ghosts_P>3GeV_Pt>0.5GeV", -3000, 3000., -3000., 3000.0, 100, 100 );
            htool->plot2D( state2.x(), state2.y(), m_title + "/XYZ2485_Ghosts", "XYZ2485_Ghosts", -1000, 1000., -1000.,
                           1000.0, 100, 100 );
          }
        }
      }
      if ( htool && m_writeHistos > 0 ) {
        htool->plot1D( tr->pseudoRapidity(), m_title + "/Eta_Total_P>3GeV_Pt>0.5GeV", "Eta_Total_P>3GeV_Pt>0.5GeV", 0.,
                       7., 50 );
        htool->plot1D( tr->momentum().Phi(), m_title + "/Phi_Total_P>3GeV_Pt>0.5GeV", "Phi_Total_P>3GeV_Pt>0.5GeV",
                       -3.142, 3.142, 25 );
        htool->plot1D( tr->pt(), m_title + "/Pt_Total_P>3GeV_Pt>0.5GeV", "Pt_Total_P>3GeV_Pt>0.5GeV", 0., 10000., 100 );
        htool->plot1D( tr->p(), m_title + "/P_Total_P>3GeV_Pt>0.5GeV", "P_Total_P>3GeV_Pt>0.5GeV", 0., 100000., 100 );
        htool->plot2D( tr->pseudoRapidity(), tr->p(), m_title + "/EtaP_Total_P>3GeV_Pt>0.5GeV",
                       "EtaP_Ghosts_P>3GeV_Pt>0.5GeV", 0., 7., 0., 100000., 20, 20 );
        htool->plot2D( tr->pseudoRapidity(), tr->momentum().Phi(), m_title + "/EtaPhi_Total_P>3GeV_Pt>0.5GeV",
                       "EtaPhi_Total_P>3GeV_Pt>0.5GeV", 0., 7., -3.142, 3.142, 20, 20 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/XYZ9000_Total_P>3GeV_Pt>0.5GeV",
                         "XYZ9000_Total_P>3GeV_Pt>0.5GeV", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/XYZ2485_Ghosts", "XYZ2485_Ghosts", -1000, 1000., -1000.,
                         1000.0, 100, 100 );
        }
      }
      m_totTrackTrigger++;
    }
  } // end track loop
} // end initialize event

void PrUTCounter::countAndPlot( IHistoTool const* htool, ITrackExtrapolator const* extrapolator,
                                LHCb::MCParticle const* part, std::vector<bool> flags, std::vector<LHCb::LHCbID>& ids,
                                int, std::map<const LHCb::Track*, double> const& trackList,
                                IGeometryInfo const& geometry ) {
  if ( flags.size() > m_name.size() ) { throw std::string( "Flag size mismatch" ); }

  unsigned int nbTrack    = 0;
  unsigned int maxRecHits = 0.;

  std::vector<LHCb::LHCbID> ttIds;

  for ( std::vector<LHCb::LHCbID>::const_iterator itId = ids.begin(); ids.end() != itId; ++itId ) {
    if ( ( *itId ).isUT() ) { ttIds.push_back( *itId ); }
  }
  std::vector<bool> shallIPlotTheHistograms( flags.size(), false );

  for ( unsigned int kk = 0; flags.size() > kk; ++kk ) {
    if ( flags[kk] ) {
      for ( auto [tr, weight] : trackList ) {
        unsigned int nbOK    = 0;
        unsigned int nbWrong = 0;

        for ( std::vector<LHCb::LHCbID>::const_iterator itId = tr->lhcbIDs().begin(); tr->lhcbIDs().end() != itId;
              ++itId ) {
          if ( ( *itId ).isUT() ) {
            LHCb::LHCbID t     = ( *itId );
            bool         found = false;
            for ( std::vector<LHCb::LHCbID>::const_iterator itMc = ttIds.begin(); ttIds.end() != itMc; ++itMc ) {
              if ( t == ( *itMc ) ) found = true;
            }
            if ( found ) {
              shallIPlotTheHistograms[kk] = true;
              nbOK++;
            } else {
              nbWrong++;
            }
          }
        }
        nbTrack++;
        m_nbTrack[kk] += 1.;
        m_mcHits[kk] += ttIds.size();
        m_foundOK[kk] += nbOK;
        m_wrong[kk] += nbWrong;
        if ( 2 < ttIds.size() ) {
          m_nbTrack3[kk] += 1.;
          m_mcHits3[kk] += ttIds.size();
          m_foundOK3[kk] += nbOK;
          m_wrong3[kk] += nbWrong;
        }
        if ( 0 < nbTrack ) { maxRecHits = std::max( maxRecHits, nbOK ); }
      }
    }
  }

  const double mcq = part->particleID().threeCharge() > 0 ? 1. : -1.;

  // -- Needed for x/y plots at a given z
  LHCb::State state, state2;
  bool        doXYZPlots = m_xyPlots;
  if ( m_xyPlots ) {
    const double vX     = part->originVertex()->position().X();
    const double vY     = part->originVertex()->position().Y();
    const double vZ     = part->originVertex()->position().Z();
    double       slopex = -10000., slopey = -10000.;
    if ( std::fabs( part->momentum().Pz() ) > 0.00001 ) {
      slopex = ( part->momentum().Px() ) / ( part->momentum().Pz() );
      slopey = ( part->momentum().Py() ) / ( part->momentum().Pz() );
    }
    const double q_over_p = ( mcq / 3 ) / ( part->momentum().P() );

    state.setState( vX, vY, vZ, slopex, slopey, q_over_p );
    StatusCode sc = extrapolator->propagate( state, 9000, geometry );
    if ( sc.isFailure() ) { doXYZPlots = false; }

    state2.setState( vX, vY, vZ, slopex, slopey, q_over_p );
    sc = extrapolator->propagate( state2, 2485, geometry );
    if ( sc.isFailure() ) { doXYZPlots = false; }

    if ( std::isnan( state.x() ) || std::isnan( state.y() ) || std::isnan( state2.x() ) || std::isnan( state2.y() ) )
      doXYZPlots = false;
  }

  for ( unsigned int k = 0; flags.size() > k; ++k ) {
    // -- Protect against nonphysical states
    if ( std::isnan( state.x() ) || std::isnan( state.y() ) ) continue;
    if ( !htool ) break;
    if ( m_writeHistos < ( m_writeHisto[k] ? 1 : 2 ) ) continue;
    if ( !flags[k] ) continue;
    htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_reconstructible",
                   m_name[k] + "_Eta_reconstructible", 0., 7., 50 );
    htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_reconstructible",
                   m_name[k] + "_Phi_reconstructible", -3.142, 3.142, 25 );
    htool->plot1D( part->momentum().Pt(), m_title + "/" + m_name[k] + "_Pt_reconstructible",
                   m_name[k] + "_Pt_reconstructible", 0., 10000., 100 );
    htool->plot1D( part->momentum().P(), m_title + "/" + m_name[k] + "_P_reconstructible",
                   m_name[k] + "_P_reconstructible", 0., 100000., 100 );
    if ( m_writeHistos > 1 ) {
      htool->plot1D( nbTrack, m_title + "/" + m_name[k] + "_expectedHits_reconstructible",
                     m_name[k] + "_expectedHits_reconstructible", -0.5, 20.5, 21 );
      htool->plot2D( part->momentum().Eta(), part->momentum().P(), m_title + "/" + m_name[k] + "_EtaP_reconstructible",
                     m_name[k] + "_EtaP_reconstructible", 0., 7., 0., 100000., 20, 20 );
      htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                     m_title + "/" + m_name[k] + "_EtaPhi_reconstructible", m_name[k] + "_EtaPhi_reconstructible", 0.,
                     7., -3.241, 3.142, 20, 20 );
      if ( doXYZPlots ) {
        htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_reconstructible",
                       "_XYZ9000_reconstructible", -3000, 3000., -3000., 3000.0, 100, 100 );
        htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_reconstructible",
                       "_XYZ2485_reconstructible", -1000, 1000., -1000., 1000.0, 100, 100 );
      }
      if ( mcq > 0 ) {
        htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_pos_reconstructible",
                       m_name[k] + "_Eta_pos_reconstructible", 0., 7., 50 );
        htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_pos_reconstructible",
                       m_name[k] + "_Phi_pos_reconstructible", -3.142, 3.142, 25 );
        htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                       m_title + "/" + m_name[k] + "_EtaP_pos_reconstructible", m_name[k] + "_EtaP_pos_reconstructible",
                       0., 7., 0., 100000., 20, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       m_title + "/" + m_name[k] + "_EtaPhi_pos_reconstructible",
                       m_name[k] + "_EtaPhi_pos_reconstructible", 0., 7., -3.241, 3.142, 20, 20 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_pos_reconstructible",
                         "_XYZ9000_pos_reconstructible", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_pos_reconstructible",
                         "_XYZ2485_pos_reconstructible", -1000, 1000., -1000., 1000.0, 100, 100 );
        }
      } else {
        htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_neg_reconstructible",
                       m_name[k] + "_Eta_neg_reconstructible", 0., 7., 50 );
        htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_neg_reconstructible",
                       m_name[k] + "_Phi_neg_reconstructible", -3.142, 3.142, 25 );
        htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                       m_title + "/" + m_name[k] + "_EtaP_neg_reconstructible", m_name[k] + "_EtaP_neg_reconstructible",
                       0., 7., 0., 100000., 20, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       m_title + "/" + m_name[k] + "_EtaPhi_neg_reconstructible",
                       m_name[k] + "_EtaPhi_neg_reconstructible", 0., 7., -3.241, 3.142, 20, 20 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_neg_reconstructible",
                         "_XYZ9000_neg_reconstructible", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_neg_reconstructible",
                         "_XYZ2485_neg_reconstructible", -1000, 1000., -1000., 1000.0, 100, 100 );
        }
      }
    }
    if ( !shallIPlotTheHistograms[k] ) continue;
    htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_reconstructed",
                   m_name[k] + "_Eta_reconstructed", 0., 7., 50 );
    htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_reconstructed",
                   m_name[k] + "_Phi_reconstructed", -3.142, 3.142, 25 );
    htool->plot1D( part->momentum().Pt(), m_title + "/" + m_name[k] + "_Pt_reconstructed",
                   m_name[k] + "_Pt_reconstructed", 0., 10000., 100 );
    htool->plot1D( part->momentum().P(), m_title + "/" + m_name[k] + "_P_reconstructed", m_name[k] + "_P_reconstructed",
                   0., 100000., 100 );
    //
    // htool->plot1D(pmeas,m_title+"/"+m_name[k]+"_P_reconstructedmeasured",m_name[k]+"_P_reconstructedmeasured",0.,100000.,50);
    // htool->plot1D(ptmeas,m_title+"/"+m_name[k]+"_Pt_reconstructedmeasured",m_name[k]+"_Pt_reconstructedmeasured",0.,10000.,50);
    // htool->plot1D(etameas,m_title+"/"+m_name[k]+"_Eta_reconstructedmeasured",m_name[k]+"_Eta_reconstructedmeasured",0.,7.,50);
    // htool->plot1D(phimeas,m_title+"/"+m_name[k]+"_Phi_reconstructedmeasured",m_name[k]+"_Phi_reconstructedmeasured",-3.142,3.142,25);

    if ( m_writeHistos > 1 ) {
      htool->plot1D( nbTrack, m_title + "/" + m_name[k] + "_expectedHits_reconstructed",
                     m_name[k] + "_expectedHits_reconstructed", -0.5, 20.5, 21 );
      htool->plot1D( maxRecHits, m_title + "/" + m_name[k] + "_reconstructedHits", m_name[k] + "_reconstructedHits",
                     -0.5, 20.5, 21 );
      htool->plot1D( maxRecHits / nbTrack, m_title + "/" + m_name[k] + "_HitEff", m_name[k] + "_HitEff", 0.0, 1.1, 50 );
      htool->plot2D( part->momentum().Eta(), part->momentum().P(), m_title + "/" + m_name[k] + "_EtaP_reconstructed",
                     m_name[k] + "_EtaP_reconstructed", 0., 7., 0., 100000., 20, 20 );
      htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                     m_title + "/" + m_name[k] + "_EtaPhi_reconstructed", m_name[k] + "_EtaPhi_reconstructed", 0., 7.,
                     -3.241, 3.142, 20, 20 );
      if ( doXYZPlots ) {
        htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_reconstructed",
                       "_XYZ9000_reconstructed", -3000, 3000., -3000., 3000.0, 100, 100 );
        htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_reconstructed",
                       "_XYZ2485_reconstructed", -1000, 1000., -1000., 1000.0, 100, 100 );
      }

      if ( mcq > 0 ) {
        htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_pos_reconstructed",
                       m_name[k] + "_Eta_pos_reconstructed", 0., 7., 50 );
        htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_pos_reconstructed",
                       m_name[k] + "_Phi_pos_reconstructed", -3.142, 3.142, 25 );
        htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                       m_title + "/" + m_name[k] + "_EtaP_pos_reconstructed", m_name[k] + "_EtaP_pos_reconstructed", 0.,
                       7., 0., 100000., 20, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       m_title + "/" + m_name[k] + "_EtaPhi_pos_reconstructed", m_name[k] + "_EtaPhi_pos_reconstructed",
                       0., 7., -3.241, 3.142, 20, 20 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_pos_reconstructed",
                         "_XYZ9000_pos_reconstructed", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_pos_reconstructed",
                         "_XYZ2485_pos_reconstructed", -1000, 1000., -1000., 1000.0, 100, 100 );
        }
      } else {
        htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_neg_reconstructed",
                       m_name[k] + "_Eta_neg_reconstructed", 0., 7., 50 );
        htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_neg_reconstructed",
                       m_name[k] + "_Phi_neg_reconstructed", -3.142, 3.142, 25 );
        htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                       m_title + "/" + m_name[k] + "_EtaP_neg_reconstructed", m_name[k] + "_EtaP_neg_reconstructed", 0.,
                       7., 0., 100000., 20, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       m_title + "/" + m_name[k] + "_EtaPhi_neg_reconstructed", m_name[k] + "_EtaPhi_neg_reconstructed",
                       0., 7., -3.241, 3.142, 20, 20 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_neg_reconstructed",
                         "_XYZ9000_neg_reconstructed", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_neg_reconstructed",
                         "_XYZ2485_neg_reconstructed", -1000, 1000., -1000., 1000.0, 100, 100 );
        }
      }
    }
  }
}

void PrUTCounter::printStatistics( MsgStream& info, std::string location ) {
  if ( m_totTrack == 0 ) return;
  m_title.resize( m_titleSize, ' ' );
  std::string strigger = "for P>3GeV,Pt>0.5GeV";
  strigger.resize( m_titleSize, ' ' );

  info << "**** UT Efficiency for " << location << " ****  ";
  if ( 0 != m_nbGhost ) {
    double bad = m_nbGhostHit / m_nbGhost;
    info << format( "%6.0f ghost, %5.2f UT per track", m_nbGhost, bad ) << endmsg;
  }
  if ( m_triggerNumbers ) {
    double gosttrig = 0;
    if ( m_totTrackTrigger != 0 ) gosttrig = 100. * m_totGhostTrigger / m_totTrackTrigger;
    info << "**** " << strigger
         << format( "%7d tracks including       %7d ghosts [%4.1f %%]  ****", m_totTrackTrigger, m_totGhostTrigger,
                    gosttrig )
         << endmsg;
  }

  for ( unsigned int kk = 0; m_name.size() > kk; ++kk ) {
    if ( 0.5 > m_nbTrack[kk] ) continue;
    double eff      = 0.;
    double fraceff  = 0.;
    double bad      = 0.;
    double fracbad  = 0.;
    double meanHits = 0.;
    if ( 0.5 < m_nbTrack[kk] ) {
      meanHits = m_mcHits[kk] / m_nbTrack[kk];
      eff      = m_foundOK[kk] / m_nbTrack[kk];
      fraceff  = 100. * eff / meanHits;
      bad      = m_wrong[kk] / m_nbTrack[kk];
      fracbad  = 100. * bad / ( eff + bad );
    }

    double eff3      = 0.;
    double fraceff3  = 0.;
    double bad3      = 0.;
    double fracbad3  = 0.;
    double meanHits3 = 0.;
    if ( 0.5 < m_nbTrack3[kk] ) {
      meanHits3 = m_mcHits3[kk] / m_nbTrack3[kk];
      eff3      = m_foundOK3[kk] / m_nbTrack3[kk];
      fraceff3  = 100. * eff3 / meanHits3;
      bad3      = m_wrong3[kk] / m_nbTrack3[kk];
      fracbad3  = 100. * bad3 / ( eff3 + bad3 );
    }

    std::string nameformat  = m_name[kk];
    std::string nameformat2 = m_name[kk];
    std::string blank( m_titleSize - ( nameformat.size() ), ' ' );
    std::string blank2( m_titleSize - ( nameformat.size() + 5.0 ), ' ' );
    nameformat2 = blank2 + nameformat;
    nameformat  = blank + nameformat;
    info << "  " << nameformat
         << format( " :%6.0f tr %5.2f from %5.2f mcUT [%5.1f %%] %5.2f ghost hits on real tracks [%4.1f %%]",
                    m_nbTrack[kk], eff, meanHits, fraceff, bad, fracbad )
         << endmsg;
    info << "  " << nameformat2 + " >3UT"
         << format( " :%6.0f tr %5.2f from %5.2f mcUT [%5.1f %%] %5.2f ghost hits on real tracks [%4.1f %%]",
                    m_nbTrack3[kk], eff3, meanHits3, fraceff3, bad3, fracbad3 )
         << endmsg;
  }
  info << endmsg;
}
