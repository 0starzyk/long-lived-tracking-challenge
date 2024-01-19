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

#include "PrTrackCounter.h"

#include "GaudiAlg/IHistoTool.h"

#include <boost/algorithm/string/replace.hpp>

//-----------------------------------------------------------------------------
// Implementation file for class : PrTrackCounter
//
// 2005-06-02 : Olivier Callot
// Modified by Wenbin Qian for VP Efficiency
//-----------------------------------------------------------------------------

void PrTrackCounter::addSelection( std::string name, bool writeHisto, bool plotNegEta ) {
  if ( name.size() > m_titleSize ) { m_titleSize = name.size(); }
  m_name.push_back( name );
  m_writeHisto.push_back( writeHisto );
  m_plotNegEta.push_back( plotNegEta );
  m_wanted.push_back( 0 );
  m_counted.push_back( 0 );
  m_velofirstcounter.push_back( 0 );
  m_velolastcounter.push_back( 0 );
  m_clone.push_back( 0 );
  m_purity.push_back( 0. );
  m_hitEff.push_back( 0. );
  m_hitEffFirstVeloHits.push_back( 0. );
  m_hitEffLastVeloHits.push_back( 0. );
}

void PrTrackCounter::initEvent( IHistoTool const* htool, ITrackExtrapolator const* extrapolator, int nPV,
                                LHCb::Track::Range const& tracks, LHCb::LinksByKey const& links,
                                IGeometryInfo const& geometry ) {

  double nbTracks = 0;
  double nbGhost  = 0;

  if ( htool && m_writeHistos > 0 ) htool->plot1D( nPV, m_title + "/nPV_Events", "nPV_Events", -0.5, 20.5, 21 );

  for ( const LHCb::Track* track : tracks ) {
    if ( track->checkFlag( LHCb::Track::Flags::Invalid ) ) continue;
    if ( ( m_trackType != LHCb::Track::Types::Unknown ) && ( track->type() != m_trackType ) ) continue;
    if ( track->ghostProbability() > m_ghostProbCut && track->ghostProbability() < 999 ) continue;
    bool eta25 = !m_eta25cut || ( track->pseudoRapidity() > 2. && track->pseudoRapidity() < 5. );
    if ( !eta25 ) continue;
    const auto& lhcbIDs = track->lhcbIDs();
    // Count the number of hits in each sub-detectors
    double n_hits_Velo  = 0.;
    double n_hits_UT    = 0.;
    double n_hits_Scifi = 0.;
    for ( std::vector<LHCb::LHCbID>::const_iterator lhcbID = lhcbIDs.begin(); lhcbIDs.end() != lhcbID; ++lhcbID ) {
      if ( ( *lhcbID ).isVP() ) {
        n_hits_Velo += 1.;
      } else if ( ( *lhcbID ).isUT() ) {
        n_hits_UT += 1.;
      } else if ( ( *lhcbID ).isFT() ) {
        n_hits_Scifi += 1.;
      }
    }

    // -- Needed for x/y plots at a given z
    LHCb::State state;
    LHCb::State state2;
    bool        doXYZPlots = m_xyPlots;
    if ( m_xyPlots ) {
      state         = track->closestState( 9000 );
      StatusCode sc = extrapolator->propagate( state, 9000, geometry );
      if ( sc.isFailure() ) { doXYZPlots = false; }

      state2 = track->closestState( 2485 );
      sc     = extrapolator->propagate( state2, 2485, geometry );
      if ( sc.isFailure() ) { doXYZPlots = false; }

      if ( std::isnan( state.x() ) || std::isnan( state.y() ) || std::isnan( state2.x() ) || std::isnan( state2.y() ) )
        doXYZPlots = false;
    }

    if ( !links.hasEntry( *track ) ) {
      nbGhost++;
      if ( htool && m_writeHistos > 0 ) {
        htool->plot1D( nPV, m_title + "/nPV_Ghosts", "nPV_Ghosts", -0.5, 20.5, 21 );
        htool->plot1D( track->pseudoRapidity(), m_title + "/Eta_Ghosts", "Eta_Ghosts", 0., 7., 50 );
        htool->plot1D( lhcbIDs.size(), m_title + "/nHits_all_Ghosts", "nHits_all_Ghosts", -0.5, 35.5, 35 );
        htool->plot1D( n_hits_Velo, m_title + "/nHits_Velo_Ghosts", "nHits_Velo_Ghosts", -0.5, 25.5, 26 );
        htool->plot1D( n_hits_UT, m_title + "/nHits_UT_Ghosts", "nHits_UT_Ghosts", -0.5, 10.5, 11 );
        htool->plot1D( n_hits_Scifi, m_title + "/nHits_Scifi_Ghosts", "nHits_Scifi_Ghosts", -0.5, 15.5, 16 );
        htool->plot2D( n_hits_Scifi, n_hits_Velo, m_title + "/nHits_Scifi_Velo_Ghosts", "nHits_Scifi_Velo_Ghosts", -0.5,
                       15.5, -0.5, 25.5, 16, 26 );
        if ( track->type() != LHCb::Track::Types::Velo ) {
          htool->plot1D( track->momentum().Phi(), m_title + "/Phi_Ghosts", "Phi_Ghosts", -3.142, 3.142, 25 );
          htool->plot1D( track->pt(), m_title + "/Pt_Ghosts", "Pt_Ghosts", 0., 10000., 100 );
          htool->plot1D( track->p(), m_title + "/P_Ghosts", "P_Ghosts", 0., 100000., 100 );
          htool->plot2D( track->pseudoRapidity(), track->p(), m_title + "/EtaP_Ghosts", "EtaP_Ghosts", 0., 7., 0.,
                         100000., 20, 20 );
          htool->plot2D( track->pseudoRapidity(), track->momentum().Phi(), m_title + "/EtaPhi_Ghosts", "EtaPhi_Ghosts",
                         0., 7., -3.142, 3.142, 20, 20 );
        }
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/XYZ9000_Ghosts", "XYZ9000_Ghosts", -3000, 3000., -3000.,
                         3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/XYZ2485_Ghosts", "XYZ2485_Ghosts", -1000, 1000., -1000.,
                         1000.0, 100, 100 );
        }
      }
    }
    if ( htool && m_writeHistos > 0 ) {
      htool->plot1D( nPV, m_title + "/nPV_Total", "nPV_Total", -0.5, 20.5, 21 );
      htool->plot1D( track->pseudoRapidity(), m_title + "/Eta_Total", "Eta_Total", 0., 7., 50 );
      htool->plot1D( lhcbIDs.size(), m_title + "/nHits_all_Total", "nHits_all_Total", -0.5, 35.5, 35 );
      htool->plot1D( n_hits_Velo, m_title + "/nHits_Velo_Total", "nHits_Velo_Total", -0.5, 25.5, 26 );
      htool->plot1D( n_hits_UT, m_title + "/nHits_UT_Total", "nHits_UT_Total", -0.5, 10.5, 11 );
      htool->plot1D( n_hits_Scifi, m_title + "/nHits_Scifi_Total", "nHits_Scifi_Total", -0.5, 15.5, 16 );
      htool->plot2D( n_hits_Scifi, n_hits_Velo, m_title + "/nHits_Scifi_Velo_Ghosts", "nHits_Scifi_Velo_Ghosts", -0.5,
                     15.5, -0.5, 25.5, 16, 26 );

      if ( track->type() != LHCb::Track::Types::Velo ) {
        htool->plot1D( track->momentum().Phi(), m_title + "/Phi_Total", "Phi_Total", -3.142, 3.142, 25 );
        htool->plot1D( track->pt(), m_title + "/Pt_Total", "Pt_Total", 0., 10000., 100 );
        htool->plot1D( track->p(), m_title + "/P_Total", "P_Total", 0., 100000., 100 );
        htool->plot2D( track->pseudoRapidity(), track->p(), m_title + "/EtaP_Total", "EtaP_Total", 0., 7., 0., 100000.,
                       20, 20 );
        htool->plot2D( track->pseudoRapidity(), track->momentum().Phi(), m_title + "/EtaPhi_Total", "EtaPhi_Total", 0.,
                       7., -3.142, 3.142, 20, 20 );
      }
      if ( doXYZPlots ) {
        htool->plot2D( state.x(), state.y(), m_title + "/XYZ9000_Total", "XYZ9000_Total", -3000, 3000., -3000., 3000.0,
                       100, 100 );
        htool->plot2D( state2.x(), state2.y(), m_title + "/XYZ2485_Total", "XYZ2485_Total", -1000, 1000., -1000.,
                       1000.0, 100, 100 );
      }
    }
    m_totTrack++;
    nbTracks++;
    if ( m_triggerNumbers && ( track->type() != LHCb::Track::Types::Velo ) && track->p() > 3000. &&
         track->pt() > 500. ) {
      if ( !links.hasEntry( *track ) ) {
        m_totGhostTrigger++;
        if ( htool && m_writeHistos > 0 ) {
          htool->plot1D( nPV, m_title + "/nPV_Ghosts_P>3GeV_Pt>0.5GeV", "nPV_Ghosts_P>3GeV_Pt>0.5GeV", -0.5, 20.5, 21 );
          htool->plot1D( track->pseudoRapidity(), m_title + "/Eta_Ghosts_P>3GeV_Pt>0.5GeV",
                         "Eta_Ghosts_P>3GeV_Pt>0.5GeV", 0., 7., 100 );
          htool->plot1D( track->momentum().Phi(), m_title + "/Phi_Ghosts_P>3GeV_Pt>0.5GeV",
                         "Phi_Ghosts_P>3GeV_Pt>0.5GeV", -3.142, 3.142, 25 );
          htool->plot1D( track->pt(), m_title + "/Pt_Ghosts_P>3GeV_Pt>0.5GeV", "Pt_Ghosts_P>3GeV_Pt>0.5GeV", 0., 10000.,
                         100 );
          htool->plot1D( track->p(), m_title + "/P_Ghosts_P>3GeV_Pt>0.5GeV", "P_Ghosts_P>3GeV_Pt>0.5GeV", 0., 100000.,
                         100 );
          htool->plot2D( track->pseudoRapidity(), track->p(), m_title + "/EtaP_Ghosts_P>3GeV_Pt>0.5GeV",
                         "EtaP_Ghosts_P>3GeV_Pt>0.5GeV", 0., 7., 0., 100000., 20, 20 );
          htool->plot2D( track->pseudoRapidity(), track->momentum().Phi(), m_title + "/EtaPhi_Ghosts_P>3GeV_Pt>0.5GeV",
                         "EtaPhi_Ghosts_P>3GeV_Pt>0.5GeV", 0., 7., -3.142, 3.142, 20, 20 );
          if ( doXYZPlots ) {
            htool->plot2D( state.x(), state.y(), m_title + "/XYZ9000_Ghosts_P>3GeV_Pt>0.5GeV",
                           "XYZ9000_Ghosts_P>3GeV_Pt>0.5GeV", -3000, 3000., -3000., 3000.0, 100, 100 );
            htool->plot2D( state2.x(), state2.y(), m_title + "/XYZ2485_Ghosts_P>3GeV_Pt>0.5GeV",
                           "XYZ2485_Ghosts_P>3GeV_Pt>0.5GeV", -1000, 1000., -1000., 1000.0, 100, 100 );
          }
        }
      }
      if ( htool && m_writeHistos > 0 ) {
        htool->plot1D( nPV, m_title + "/nPV_Total_P>3GeV_Pt>0.5GeV", "nPV_Total_P>3GeV_Pt>0.5GeV", -0.5, 20.5, 21 );
        htool->plot1D( track->pseudoRapidity(), m_title + "/Eta_Total_P>3GeV_Pt>0.5GeV", "Eta_Total_P>3GeV_Pt>0.5GeV",
                       0., 7., 100 );
        htool->plot1D( track->momentum().Phi(), m_title + "/Phi_Total_P>3GeV_Pt>0.5GeV", "Phi_Total_P>3GeV_Pt>0.5GeV",
                       -3.142, 3.142, 25 );
        htool->plot1D( track->pt(), m_title + "/Pt_Total_P>3GeV_Pt>0.5GeV", "Pt_Total_P>3GeV_Pt>0.5GeV", 0., 10000.,
                       100 );
        htool->plot1D( track->p(), m_title + "/P_Total_P>3GeV_Pt>0.5GeV", "P_Total_P>3GeV_Pt>0.5GeV", 0., 100000.,
                       100 );
        htool->plot2D( track->pseudoRapidity(), track->p(), m_title + "/EtaP_Total_P>3GeV_Pt>0.5GeV",
                       "EtaP_Ghosts_P>3GeV_Pt>0.5GeV", 0., 7., 0., 100000., 20, 20 );
        htool->plot2D( track->pseudoRapidity(), track->momentum().Phi(), m_title + "/EtaPhi_Total_P>3GeV_Pt>0.5GeV",
                       "EtaPhi_Total_P>3GeV_Pt>0.5GeV", 0., 7., -3.142, 3.142, 20, 20 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/XYZ9000_Total_P>3GeV_Pt>0.5GeV",
                         "XYZ9000_Total_P>3GeV_Pt>0.5GeV", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/XYZ2485_Total_P>3GeV_Pt>0.5GeV",
                         "XYZ2485_Total_P>3GeV_Pt>0.5GeV", -1000, 1000., -1000., 1000.0, 100, 100 );
        }
      }
      m_totTrackTrigger++;
    }
  }
  m_totGhost += nbGhost;
  double fracGhost = 0.;
  if ( 0 < nbTracks ) fracGhost = double( nbGhost ) / nbTracks;
  m_fracGhost += fracGhost;
  m_nEvent += 1.;
}

void PrTrackCounter::countAndPlot( IHistoTool const* htool, ITrackExtrapolator const* extrapolator,
                                   LHCb::MCParticle const* part, std::vector<bool> flags,
                                   std::vector<LHCb::LHCbID>& ids, int nPV,
                                   std::map<const LHCb::Track*, double> const& trackList,
                                   IGeometryInfo const&                        geometry ) {

  // LHCbIDs should be ordered
  assert( std::is_sorted( ids.begin(), ids.end() ) );

  if ( flags.size() > m_name.size() ) { throw std::string( "Flag size mismatch" ); }

  bool found = false;
  int  clone = 0;

  // trackList link tracks related to the MCparticle passed
  if ( LHCb::Track::Types::Unknown == m_trackType ) {
    // TODO loop

    if ( trackList.size() != 0 ) {
      found = true;
      clone = trackList.size() - 1;
    }
  } else {
    for ( auto [tr, weight] : trackList ) {
      if ( tr->ghostProbability() > m_ghostProbCut && tr->ghostProbability() < 999 ) continue;
      if ( tr->type() == m_trackType ) {
        if ( !found ) {
          found = true;
        } else {
          clone++;
        }
      }
    }
  }

  //== Count how many of the proper type...
  double nTrue = 0.;

  if ( HitType::Unspecified == m_hitTypesToCheck ) {
    nTrue = double( ids.size() );
  } else {
    for ( std::vector<LHCb::LHCbID>::const_iterator itId = ids.begin(); ids.end() != itId; ++itId ) {
      if ( ( *itId ).isVP() ) {
        if ( 0 != ( m_hitTypesToCheck & HitType::VP ) ) nTrue += 1.;
      } else if ( ( *itId ).isUT() ) {
        if ( 0 != ( m_hitTypesToCheck & HitType::UT ) ) nTrue += 1.;
      } else if ( ( *itId ).isFT() ) {
        if ( 0 != ( m_hitTypesToCheck & HitType::FT ) ) nTrue += 1.;
      }
    }
  }
  unsigned int maxRecHits = 0.;
  for ( unsigned int kk = 0; flags.size() > kk; ++kk ) {
    if ( flags[kk] ) {
      m_wanted[kk]++;
      if ( found ) {
        m_counted[kk]++;
        m_clone[kk] += clone;
        for ( auto [tr, weight] : trackList ) {
          if ( ( m_trackType != LHCb::Track::Types::Unknown ) && ( tr->type() != m_trackType ) ) continue;
          if ( tr->ghostProbability() > m_ghostProbCut && tr->ghostProbability() < 999 ) continue;
          m_purity[kk] += weight;
          unsigned int nbMeas              = 0;
          unsigned int nbMeasFirstVeloHits = 0;
          unsigned int nbMeaslastVeloHits  = 0;
          for ( std::vector<LHCb::LHCbID>::const_iterator itId = tr->lhcbIDs().begin(); tr->lhcbIDs().end() != itId;
                ++itId ) {
            auto const foundID = std::find( ids.begin(), ids.end(), *itId );
            if ( foundID == ids.end() ) continue;
            if ( ( *itId ).isVP() ) {
              if ( 0 != ( m_hitTypesToCheck & HitType::VP ) ) {
                if ( std::distance( ids.begin(), foundID ) < m_firstNVeloHits ) {
                  nbMeasFirstVeloHits += 1;
                } else {
                  nbMeaslastVeloHits += 1;
                }
                nbMeas += 1;
              }
            } else if ( ( *itId ).isFT() ) {
              if ( 0 != ( m_hitTypesToCheck & HitType::FT ) ) nbMeas += 1;
            } else if ( ( *itId ).isUT() ) {
              if ( 0 != ( m_hitTypesToCheck & HitType::UT ) ) nbMeas += 1;
            }
          }
          if ( 0 < nTrue ) {
            maxRecHits = std::max( maxRecHits, nbMeas );
            if ( m_hitTypesToCheck == HitType::VP ) {
              if ( nTrue > m_firstNVeloHits ) {
                double eff = double( nbMeaslastVeloHits ) / ( nTrue - m_firstNVeloHits );
                m_hitEffLastVeloHits[kk] += eff;
                ++m_velolastcounter[kk];
              }
              m_hitEffFirstVeloHits[kk] += double( nbMeasFirstVeloHits ) / double( m_firstNVeloHits );
              ++m_velofirstcounter[kk];
            }
            double eff = double( nbMeas ) / nTrue;
            m_hitEff[kk] += eff;
          }
          htool->plot1D( nbMeas, m_title + "/" + m_name[kk] + "_nbMeasHits_reconstructed",
                         m_name[kk] + "_nbMeasHits_reconstructed", -0.5, 20.5, 21 );
        }
      } // end found loop
    }
  }

  double prodx = part->originVertex()->position().X();
  double prody = part->originVertex()->position().Y();
  double docaz = -100.0;
  if ( part->momentum().Pt() > 0.00001 ) {
    docaz = std::abs( 1. / part->momentum().Pt() * ( prodx * part->momentum().Py() - prody * part->momentum().Px() ) );
  }

  const LHCb::MCParticle* mother = part;
  while ( mother->mother() != NULL ) mother = mother->mother();
  double PVz = mother->originVertex()->position().Z();
  double mcq = part->particleID().threeCharge() > 0 ? 1. : -1.;

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
    htool->plot1D( nPV, m_title + "/" + m_name[k] + "_nPV_reconstructible", m_name[k] + "_nPV_reconstructible", -0.5,
                   20.5, 21 );
    if ( m_plotNegEta[k] )
      htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_reconstructible",
                     m_name[k] + "_Eta_reconstructible", -7., 7., 100 );
    else
      htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_reconstructible",
                     m_name[k] + "_Eta_reconstructible", 0., 7., 50 );
    htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_reconstructible",
                   m_name[k] + "_Phi_reconstructible", -3.142, 3.142, 25 );
    htool->plot1D( part->momentum().Pt(), m_title + "/" + m_name[k] + "_Pt_reconstructible",
                   m_name[k] + "_Pt_reconstructible", 0., 10000., 100 );
    htool->plot1D( part->momentum().P(), m_title + "/" + m_name[k] + "_P_reconstructible",
                   m_name[k] + "_P_reconstructible", 0., 100000., 100 );
    htool->plot1D( part->originVertex()->position().Z(), m_title + "/" + m_name[k] + "_z_reconstructible",
                   m_name[k] + "_z_reconstructible", -550., 200., 200 );

    if ( m_writeHistos > 1 ) {
      htool->plot1D( nTrue, m_title + "/" + m_name[k] + "_expectedHits_reconstructible",
                     m_name[k] + "_expectedHits_reconstructible", -0.5, 20.5, 21 );
      htool->plot1D( docaz, m_title + "/" + m_name[k] + "_docaz_reconstructible", m_name[k] + "_docaz_reconstructible",
                     0., 10., 50 );
      htool->plot1D( PVz, m_title + "/" + m_name[k] + "_PVz_reconstructible", m_name[k] + "_PVz_reconstructible", -200.,
                     200., 50 );
      if ( m_plotNegEta[k] ) {
        htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                       m_title + "/" + m_name[k] + "_EtaP_reconstructible", m_name[k] + "_EtaP_reconstructible", -7.,
                       7., 0., 100000., 40, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       m_title + "/" + m_name[k] + "_EtaPhi_reconstructible", m_name[k] + "_EtaPhi_reconstructible",
                       -7., 7., -3.142, 3.142, 40, 20 );
      } else {
        htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                       m_title + "/" + m_name[k] + "_EtaP_reconstructible", m_name[k] + "_EtaP_reconstructible", 0., 7.,
                       0., 100000., 20, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       m_title + "/" + m_name[k] + "_EtaPhi_reconstructible", m_name[k] + "_EtaPhi_reconstructible", 0.,
                       7., -3.142, 3.142, 20, 20 );
      }
      if ( doXYZPlots ) {
        htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_reconstructible",
                       "_XYZ9000_reconstructible", -3000, 3000., -3000., 3000.0, 100, 100 );
        htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_reconstructible",
                       "_XYZ2485_reconstructible", -1000, 1000., -1000., 1000.0, 100, 100 );
      }
      if ( mcq > 0 ) {
        if ( m_plotNegEta[k] ) {
          htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_pos_reconstructible",
                         m_name[k] + "_Eta_pos_reconstructible", -7., 7., 100 );
          htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                         m_title + "/" + m_name[k] + "_EtaP_pos_reconstructible",
                         m_name[k] + "_EtaP_pos_reconstructible", -7., 7., 0., 100000., 40, 20 );
          htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                         m_title + "/" + m_name[k] + "_EtaPhi_pos_reconstructible",
                         m_name[k] + "_EtaPhi_pos_reconstructible", -7., 7., -3.142, 3.142, 40, 20 );
        } else {
          htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_pos_reconstructible",
                         m_name[k] + "_Eta_pos_reconstructible", 0., 7., 50 );
          htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                         m_title + "/" + m_name[k] + "_EtaP_pos_reconstructible",
                         m_name[k] + "_EtaP_pos_reconstructible", 0., 7., 0., 100000., 20, 20 );
          htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                         m_title + "/" + m_name[k] + "_EtaPhi_pos_reconstructible",
                         m_name[k] + "_EtaPhi_pos_reconstructible", 0., 7., -3.142, 3.142, 20, 20 );
        }
        htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_pos_reconstructible",
                       m_name[k] + "_Phi_pos_reconstructible", -3.142, 3.142, 25 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_pos_reconstructible",
                         "_XYZ9000_pos_reconstructible", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_pos_reconstructible",
                         "_XYZ2485_pos_reconstructible", -1000, 1000., -1000., 1000.0, 100, 100 );
        }
      } else {
        if ( m_plotNegEta[k] ) {
          htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_neg_reconstructible",
                         m_name[k] + "_Eta_neg_reconstructible", -7., 7., 100 );
          htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                         m_title + "/" + m_name[k] + "_EtaP_neg_reconstructible",
                         m_name[k] + "_EtaP_neg_reconstructible", -7., 7., 0., 100000., 40, 20 );
          htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                         m_title + "/" + m_name[k] + "_EtaPhi_neg_reconstructible",
                         m_name[k] + "_EtaPhi_neg_reconstructible", -7., 7., -3.142, 3.142, 40, 20 );
        } else {
          htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_neg_reconstructible",
                         m_name[k] + "_Eta_neg_reconstructible", 0., 7., 50 );
          htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                         m_title + "/" + m_name[k] + "_EtaP_neg_reconstructible",
                         m_name[k] + "_EtaP_neg_reconstructible", 0., 7., 0., 100000., 20, 20 );
          htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                         m_title + "/" + m_name[k] + "_EtaPhi_neg_reconstructible",
                         m_name[k] + "_EtaPhi_neg_reconstructible", 0., 7., -3.142, 3.142, 20, 20 );
        }
        htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_neg_reconstructible",
                       m_name[k] + "_Phi_neg_reconstructible", -3.142, 3.142, 25 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_neg_reconstructible",
                         "_XYZ900_neg_reconstructible", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_neg_reconstructible",
                         "_XYZ2485_neg_reconstructible", -1000, 1000., -1000., 1000.0, 100, 100 );
        }
      }
    }
    if ( !found ) continue;
    htool->plot1D( nPV, m_title + "/" + m_name[k] + "_nPV_reconstructed", m_name[k] + "_nPV_reconstructed", -0.5, 20.5,
                   21 );
    if ( m_plotNegEta[k] )
      htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_reconstructed",
                     m_name[k] + "_Eta_reconstructed", -7., 7., 100 );
    else
      htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_reconstructed",
                     m_name[k] + "_Eta_reconstructed", 0., 7., 50 );
    htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_reconstructed",
                   m_name[k] + "_Phi_reconstructed", -3.142, 3.142, 25 );
    htool->plot1D( part->momentum().Pt(), m_title + "/" + m_name[k] + "_Pt_reconstructed",
                   m_name[k] + "_Pt_reconstructed", 0., 10000., 100 );
    htool->plot1D( part->momentum().P(), m_title + "/" + m_name[k] + "_P_reconstructed", m_name[k] + "_P_reconstructed",
                   0., 100000., 100 );
    htool->plot1D( part->originVertex()->position().Z(), m_title + "/" + m_name[k] + "_z_reconstructed",
                   m_name[k] + "_z_reconstructed", -550., 200., 200 );
    //
    // htool->plot1D(pmeas,m_title+"/"+m_name[k]+"_P_reconstructedmeasured",m_name[k]+"_P_reconstructedmeasured",0.,100000.,100);
    // htool->plot1D(ptmeas,m_title+"/"+m_name[k]+"_Pt_reconstructedmeasured",m_name[k]+"_Pt_reconstructedmeasured",0.,10000.,100);
    // htool->plot1D(etameas,m_title+"/"+m_name[k]+"_Eta_reconstructedmeasured",m_name[k]+"_Eta_reconstructedmeasured",0.,7.,50);
    // htool->plot1D(phimeas,m_title+"/"+m_name[k]+"_Phi_reconstructedmeasured",m_name[k]+"_Phi_reconstructedmeasured",-3.142,3.142,25);

    if ( m_writeHistos > 1 ) {
      htool->plot1D( nTrue, m_title + "/" + m_name[k] + "_expectedHits_reconstructed",
                     m_name[k] + "_expectedHits_reconstructed", -0.5, 20.5, 21 );
      htool->plot1D( maxRecHits, m_title + "/" + m_name[k] + "_reconstructedHits", m_name[k] + "_reconstructedHits",
                     -0.5, 20.5, 21 );
      htool->plot1D( maxRecHits / nTrue, m_title + "/" + m_name[k] + "_HitEff", m_name[k] + "_HitEff", 0.0, 1.1, 50 );
      htool->plot1D( docaz, m_title + "/" + m_name[k] + "_docaz_reconstructed", m_name[k] + "_docaz_reconstructed", 0.,
                     10., 50 );
      htool->plot1D( PVz, m_title + "/" + m_name[k] + "_PVz_reconstructed", m_name[k] + "_PVz_reconstructed", -200.,
                     200., 50 );
      if ( m_plotNegEta[k] ) {
        htool->plot2D( part->momentum().Eta(), part->momentum().P(), m_title + "/" + m_name[k] + "_EtaP_reconstructed",
                       m_name[k] + "_EtaP_reconstructed", -7., 7., 0., 100000., 40, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       m_title + "/" + m_name[k] + "_EtaPhi_reconstructed", m_name[k] + "_EtaPhi_reconstructed", -7.,
                       7., -3.142, 3.142, 40, 20 );
      } else {
        htool->plot2D( part->momentum().Eta(), part->momentum().P(), m_title + "/" + m_name[k] + "_EtaP_reconstructed",
                       m_name[k] + "_EtaP_reconstructed", 0., 7., 0., 100000., 20, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       m_title + "/" + m_name[k] + "_EtaPhi_reconstructed", m_name[k] + "_EtaPhi_reconstructed", 0., 7.,
                       -3.142, 3.142, 20, 20 );
      }
      if ( doXYZPlots ) {
        htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_reconstructed",
                       "_XYZ9000_reconstructed", -3000, 3000., -3000., 3000.0, 100, 100 );
        htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_reconstructed",
                       "_XYZ2485_reconstructed", -1000, 1000., -1000., 1000.0, 100, 100 );
      }
      if ( mcq > 0 ) {
        if ( m_plotNegEta[k] ) {
          htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_pos_reconstructed",
                         m_name[k] + "_Eta_pos_reconstructed", -7., 7., 100 );
          htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                         m_title + "/" + m_name[k] + "_EtaP_pos_reconstructed", m_name[k] + "_EtaP_pos_reconstructed",
                         -7., 7., 0., 100000., 40, 20 );
          htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                         m_title + "/" + m_name[k] + "_EtaPhi_pos_reconstructed",
                         m_name[k] + "_EtaPhi_pos_reconstructed", -7., 7., -3.142, 3.142, 40, 20 );
        } else {
          htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_pos_reconstructed",
                         m_name[k] + "_Eta_pos_reconstructed", 0., 7., 50 );
          htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                         m_title + "/" + m_name[k] + "_EtaP_pos_reconstructed", m_name[k] + "_EtaP_pos_reconstructed",
                         0., 7., 0., 100000., 20, 20 );
          htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                         m_title + "/" + m_name[k] + "_EtaPhi_pos_reconstructed",
                         m_name[k] + "_EtaPhi_pos_reconstructed", 0., 7., -3.142, 3.142, 20, 20 );
        }
        htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_pos_reconstructed",
                       m_name[k] + "_Phi_pos_reconstructed", -3.142, 3.142, 25 );
        if ( doXYZPlots ) {
          htool->plot2D( state.x(), state.y(), m_title + "/" + m_name[k] + "_XYZ9000_pos_reconstructed",
                         "_XYZ9000_pos_reconstructed", -3000, 3000., -3000., 3000.0, 100, 100 );
          htool->plot2D( state2.x(), state2.y(), m_title + "/" + m_name[k] + "_XYZ2485_pos_reconstructed",
                         "_XYZ2485_pos_reconstructed", -1000, 1000., -1000., 1000.0, 100, 100 );
        }
      } else {
        if ( m_plotNegEta[k] ) {
          htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_neg_reconstructed",
                         m_name[k] + "_Eta_neg_reconstructed", -7., 7., 100 );
          htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                         m_title + "/" + m_name[k] + "_EtaP_neg_reconstructed", m_name[k] + "_EtaP_neg_reconstructed",
                         -7., 7., 0., 100000., 40, 20 );
          htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                         m_title + "/" + m_name[k] + "_EtaPhi_neg_reconstructed",
                         m_name[k] + "_EtaPhi_neg_reconstructed", -7., 7., -3.142, 3.142, 40, 20 );
        } else {
          htool->plot1D( part->momentum().Eta(), m_title + "/" + m_name[k] + "_Eta_neg_reconstructed",
                         m_name[k] + "_Eta_neg_reconstructed", 0., 7., 50 );
          htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                         m_title + "/" + m_name[k] + "_EtaP_neg_reconstructed", m_name[k] + "_EtaP_neg_reconstructed",
                         0., 7., 0., 100000., 20, 20 );
          htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                         m_title + "/" + m_name[k] + "_EtaPhi_neg_reconstructed",
                         m_name[k] + "_EtaPhi_neg_reconstructed", 0., 7., -3.142, 3.142, 20, 20 );
        }
        htool->plot1D( part->momentum().Phi(), m_title + "/" + m_name[k] + "_Phi_neg_reconstructed",
                       m_name[k] + "_Phi_neg_reconstructed", -3.142, 3.142, 25 );
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

void PrTrackCounter::printStatistics( MsgStream& info, std::string ) {
  if ( 0 == m_nEvent ) return;
  std::string short_title = m_title;

  double totT = m_totTrack + 0.00000000001;
  double frac = 100. * double( m_totGhost ) / totT;
  m_title.resize( m_titleSize, ' ' );
  std::string strigger      = "for P>3GeV,Pt>0.5GeV";
  std::string short_trigger = strigger;
  boost::replace_all( short_trigger, ">", "$>$" );
  strigger.resize( m_titleSize, ' ' );

  info << "**** " << m_title
       << format( "%7d tracks including        %7d ghosts [%5.2f %%], Event average %5.2f %% ****", m_totTrack,
                  m_totGhost, frac, 100. * m_fracGhost / m_nEvent )
       << endmsg;
  if ( m_triggerNumbers && m_totTrackTrigger > 0 )
    info << "**** " << strigger
         << format( "%7d tracks including       %7d ghosts [%4.2f %%]  ****", m_totTrackTrigger, m_totGhostTrigger,
                    100. * m_totGhostTrigger / m_totTrackTrigger )
         << endmsg;

  std::FILE* table = nullptr;
  if ( m_writetex ) {
    std::string outfile = m_texoutdir + m_texoutname + "_" + short_title + ".tex";
    info << "writing LaTeX table to " << outfile << endmsg;
    table = std::fopen( outfile.c_str(), "w" );
    std::fprintf( table, "\\begin{table}[]\n" );
    std::fprintf( table, "\t\\begin{center}\n" );
    std::fprintf( table, "\t\\resizebox*{\\textwidth}{!}{\n" );
    std::fprintf(
        table,
        "\t\t\\begin{tabular}{rr@{ }lr@{ $[$}r@{ $\\%%]$}r@{ }l@{ $[$}r@{ $\\%%]$}l@{ }r@{ $\\%%$, hitEff: }r}\n" );
    std::fprintf( table,
                  "\t\t\t\\multicolumn{1}{@{}l}{\\textbf{%s}} & %7d & \\multicolumn{3}{@{}l}{tracks including} & %7d & "
                  "ghosts & %4.2f & \\multicolumn{3}{@{}l}{, Event average %5.2f \\%%}\\\\ \n",
                  short_title.c_str(), m_totTrack, m_totGhost, frac, 100. * m_fracGhost / m_nEvent );
    if ( m_triggerNumbers && m_totTrackTrigger > 0 ) {
      std::fprintf( table,
                    "\t\t\t\\multicolumn{1}{@{}l}{ %s } & %7d & \\multicolumn{3}{@{}l}{tracks including} & %7d & "
                    "ghosts & %4.2f & \\multicolumn{3}{l}{ }\\\\ \n",
                    short_trigger.c_str(), m_totTrackTrigger, m_totGhostTrigger,
                    100. * m_totGhostTrigger / m_totTrackTrigger );
    }
  }

  for ( unsigned int kk = 0; m_name.size() > kk; ++kk ) {
    if ( 0 == m_wanted[kk] ) continue;
    frac                   = 100. * double( m_counted[kk] ) / double( m_wanted[kk] );
    double      nTot       = double( m_counted[kk] + m_clone[kk] ) + 0.00000001;
    double      clo        = 100. * double( m_clone[kk] ) / nTot;
    double      purity     = 100. * m_purity[kk] / nTot;
    double      hitEff     = 100. * m_hitEff[kk] / nTot;
    std::string nameformat = m_name[kk];
    std::string blank( m_titleSize - ( nameformat.size() ), ' ' );
    nameformat = nameformat + blank;
    info << "  " << nameformat
         << format( " :%8d from %8d [%6.2f %%] %6d clones [%5.2f %%]", m_counted[kk], m_wanted[kk], frac, m_clone[kk],
                    clo );

    if ( m_hitTypesToCheck == HitType::VP ) {
      double VhitEffFirst = m_velofirstcounter[kk] > 0 ? 100. * m_hitEffFirstVeloHits[kk] / m_velofirstcounter[kk] : -1;
      double VhitEffLast  = m_velolastcounter[kk] > 0 ? 100. * m_hitEffLastVeloHits[kk] / m_velolastcounter[kk] : -1;
      info << format( ", purity:%6.2f %%, hitEff:%6.2f %%, hitEffFirst%s:%6.2f %%, hitEffLast:%6.2f %%", purity, hitEff,
                      std::to_string( m_firstNVeloHits ).c_str(), VhitEffFirst, VhitEffLast )
           << endmsg;
    } else {
      info << format( ", purity:%6.2f %%, hitEff:%6.2f %%", purity, hitEff ) << endmsg;
    }

    if ( m_writetex ) {
      boost::replace_all( nameformat, "_", "\\_" );
      boost::replace_all( nameformat, ">", "$>$" );
      std::fprintf( table,
                    "%s : & %8d & from & %8d & %5.2f & %6d & clones & %4.2f &, purity: & %6.2f & %6.2f \\%% \\\\ \n",
                    nameformat.c_str(), m_counted[kk], m_wanted[kk], frac, m_clone[kk], clo, purity, hitEff );
    }
  }
  info << endmsg;

  if ( m_writetex ) {
    std::fprintf( table, "\t\t\\end{tabular}\n" );
    std::fprintf( table, "\t }\n" );
    std::fprintf( table, "\t \\caption[%s %s]{%s %s}\\label{tab:%s_%s}\n", m_texoutname.c_str(), short_title.c_str(),
                  m_texoutname.c_str(), short_title.c_str(), m_texoutname.c_str(), short_title.c_str() );
    std::fprintf( table, "\t\\end{center}\n" );
    std::fprintf( table, "\\end{table}\n" );
    std::fclose( table );
  }
}
