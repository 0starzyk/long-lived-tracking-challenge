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

#include "Event/Track.h"
#include "Kernel/LHCbID.h"
#include "LHCbAlgs/Transformer.h"

#include "Gaudi/Accumulators.h"
#include "GaudiKernel/SmartIF.h"
#include "GaudiKernel/SystemOfUnits.h"

#include <map>
#include <memory>
#include <string>

namespace LHCb {

  /**
   *  A small algorithm that filters tracks of a certain type.
   *
   *  @author Jan Amoraal
   *  @date   2007-07-11
   */
  class MuonTrackSelector : public LHCb::Algorithm::Transformer<LHCb::Tracks( const LHCb::Tracks& )> {

  public:
    typedef std::vector<LHCb::LHCbID>                               LHCBIDS;
    typedef std::map<std::string, bool ( LHCb::LHCbID::* )() const> LHCbDetChecks;

    MuonTrackSelector( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer( name, pSvcLocator, {KeyValue{"TracksInputContainer", TrackLocation::Default}},
                       KeyValue{"TracksOutputContainer", "Alignment/FilteredTracks"} ) {}

    LHCb::Tracks operator()( const LHCb::Tracks& inputTracks ) const override;

  private:
    Gaudi::Property<int>         m_nStation{this, "minHitStation", 2};
    Gaudi::Property<int>         m_theR{this, "TheRegion", -1};
    Gaudi::Property<double>      m_pcut{this, "MuonPcut", 0. * Gaudi::Units::GeV};
    Gaudi::Property<double>      m_muonChisquareCut{this, "MuonChisquareCut", 0.};
    Gaudi::Property<bool>        m_calo{this, "inCaloAcceptance", false};
    Gaudi::Property<bool>        m_noOverlap{this, "noOverlap", false};
    Gaudi::Property<std::string> m_trackType{this, "TrackType", "Long"};

    mutable Gaudi::Accumulators::Counter<>          m_lowMomentum{this, "Not selected (low momentum)"};
    mutable Gaudi::Accumulators::Counter<>          m_unwantedRegion{this, "Not selected (unwanted region)"};
    mutable Gaudi::Accumulators::Counter<>          m_doesnothitRegion{this, "Not selected (does not hit region)"};
    mutable Gaudi::Accumulators::Counter<>          m_moreThanOneRegion{this, "Not selected (more than one region)"};
    mutable Gaudi::Accumulators::Counter<>          m_lowNbHits{this, "Not selected (low number of hit station)"};
    mutable Gaudi::Accumulators::Counter<>          m_overlapHitStation{this, "Not selected (overlap hit stations)"};
    mutable Gaudi::Accumulators::Counter<>          m_offCaloAcceptance{this, "Not selected (off the CALO acceptance)"};
    mutable Gaudi::Accumulators::Counter<>          m_badChisquare{this, "Not selected (bad chisquare)"};
    mutable Gaudi::Accumulators::AveragingCounter<> m_filteredTracks{this, "Nb selected tracks"};
    mutable Gaudi::Accumulators::AveragingCounter<> m_nbInputTracks{this, "Nb input tracks"};

    bool printDebug() const { return msgLevel( MSG::DEBUG ); };

    void filterMuonTrack( LHCb::Track* track, LHCb::Tracks& outputContainer ) const;
  };

  DECLARE_COMPONENT_WITH_ID( MuonTrackSelector, "MuonTrackSelector" )

  void MuonTrackSelector::filterMuonTrack( Track* track, Tracks& outputContainer ) const {
    if ( m_pcut && track->p() < m_pcut ) {
      ++m_lowMomentum;
      debug() << " Discard the track due to the low momentum" << track->p() << endmsg;
      return; // select high momentum tracks
    }

    bool ASide = false;
    bool CSide = false;
    int  MA[5] = {0, 0, 0, 0, 0};
    int  MC[5] = {0, 0, 0, 0, 0};
    int  MR[4] = {0, 0, 0, 0};

    const std::vector<LHCb::LHCbID>& lhcbids = track->lhcbIDs();
    for ( std::vector<LHCb::LHCbID>::const_iterator id = lhcbids.begin(); id != lhcbids.end(); ++id ) {
      if ( id->isMuon() ) {
        int iS = id->muonID().station();
        int iR = id->muonID().region();
        int iQ = id->muonID().quarter();
        if ( iQ < 2 )
          MA[iS] = 1; // A-side
        else
          MC[iS] = 1;
        MR[iR]++;
      }
    }

    int MAside = MA[0] + MA[1] + MA[2] + MA[3] + MA[4];
    int MCside = MC[0] + MC[1] + MC[2] + MC[3] + MC[4];

    if ( m_theR > -1 && m_theR < 4 ) {
      for ( int iR = 0; iR < 4; iR++ ) {
        if ( iR != m_theR && MR[iR] != 0 ) {
          debug() << " Discard the track since it hits unwanted Region" << iR << endmsg;
          ++m_unwantedRegion;
          return;
        }
      }
      if ( MR[m_theR] == 0 ) {
        debug() << " Discard the track since it doesn't hit Region" << m_theR << endmsg;
        ++m_doesnothitRegion;
        return;
      }
    } else if ( m_theR > 9 && m_theR < 40 ) {
      for ( int iR = 0; iR < 4; iR++ ) {
        if ( iR > int( m_theR / 10 ) && MR[iR] != 0 ) {
          debug() << " Discard the track since it hits unwanted Region" << iR << endmsg;
          ++m_unwantedRegion;
          return;
        }
      }
    } else if ( m_theR == 10 ) {
      int iMR = 0;
      for ( int iR = 0; iR < 4; iR++ ) {
        if ( MR[iR] != 0 ) iMR++;
      }
      if ( iMR > 1 ) {
        debug() << " Discard the track since it hits more than one region" << iMR << endmsg;
        ++m_moreThanOneRegion;
        return;
      }
    }

    if ( MAside != 0 && MCside == 0 )
      ASide = true;
    else if ( MAside == 0 && MCside != 0 )
      CSide = true;

    if ( MAside + MCside < m_nStation ) {
      debug() << " Discard the track due to the low number of hit station " << MAside + MCside << endmsg;
      ++m_lowNbHits;
      return; /// requires at least some hits somewhere
    }

    if ( m_noOverlap && !( ( ASide || CSide ) && ( MAside > m_nStation || MCside > m_nStation ) ) ) {
      debug() << " Discard the track since overlaps hit station Cside " << MCside << " Aside " << MAside << endmsg;
      ++m_overlapHitStation;
      return;
    } else {
      debug() << " Track hit station Cside " << MCside << " Aside " << MAside << endmsg;
    }

    if ( m_calo ) {
      State stateAtCALO;
      if ( fabs( stateAtCALO.x() ) > 3900 && fabs( stateAtCALO.y() ) > 3150 ) {
        debug() << " Discard the track since falls off the CALO acceptance x" << fabs( stateAtCALO.x() ) << " y "
                << stateAtCALO.y() << endmsg;
        ++m_offCaloAcceptance;
        return; // out the calo acceptance
      }
    }

    if ( m_muonChisquareCut > 0 && track->chi2PerDoF() > m_muonChisquareCut ) {
      debug() << " Discard the track due to the chisquare " << track->chi2PerDoF() << endmsg;
      ++m_badChisquare;
      return;
    }

    debug() << " Track selected " << endmsg;
    auto    clonedTrack = std::make_unique<Track>( *track, track->key() );
    LHCBIDS ids         = clonedTrack->lhcbIDs();
    outputContainer.add( clonedTrack.release() );
  }

  Tracks MuonTrackSelector::operator()( const Tracks& inputTracks ) const {
    m_nbInputTracks += inputTracks.size();
    Tracks output;
    for ( Track* t : inputTracks ) { this->filterMuonTrack( t, output ); }
    m_filteredTracks += output.size();
    return output;
  }

} // namespace LHCb
