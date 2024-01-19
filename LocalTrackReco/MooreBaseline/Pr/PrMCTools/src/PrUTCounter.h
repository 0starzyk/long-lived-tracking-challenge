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
#ifndef PRUTCOUNTER_H
#define PRUTCOUNTER_H 1

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/Track.h"

#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/IHistoTool.h"

#include "TrackInterfaces/ITrackExtrapolator.h"

#include <memory>

/** @class PrUTCounter PrUTCounter.h
 *
 *  @author Olivier Callot
 *  @date   2006-06-28
 *  2015-01-17 : modified by Svende Braun, Michel de Cian to work with PrChecker2.cpp
 */
class PrUTCounter final {
public:
  void initEvent( IHistoTool const* htool, ITrackExtrapolator const* extrapolator, int nPV,
                  LHCb::Track::Range const& tracks, LHCb::LinksByKey const& links, IGeometryInfo const& geometry );

  void countAndPlot( IHistoTool const* htool, ITrackExtrapolator const* extrapolator, LHCb::MCParticle const* part,
                     std::vector<bool> flags, std::vector<LHCb::LHCbID>& ids, int nPV,
                     std::map<const LHCb::Track*, double> const& trackList, IGeometryInfo const& geometry );

  void addSelection( std::string name, bool writeHisto, bool plotNegEta = false );

  void printStatistics( MsgStream& info, std::string location );

  std::string title() { return m_title; }
  void        setTitle( std::string& title ) { m_title = title; }
  void        setFirstNVeloHits( unsigned int ) {}
  void        setWriteHistos( int write ) { m_writeHistos = write; };
  void        setUseEta25Cut( bool cut ) { m_eta25cut = cut; };
  void        setTriggerNumbers( bool numbers ) { m_triggerNumbers = numbers; };
  void        setXYPlots( bool xyPlots ) { m_xyPlots = xyPlots; };
  void        setHitTypesToCheck( int ){};
  void        setTrackType( LHCb::Track::Types ) {}
  void        setGhostProbCut( double ) {}
  void        setTeXName( const std::string&, const std::string& ) {}

private:
  int  m_writeHistos{-1};
  bool m_eta25cut{false};
  bool m_triggerNumbers{false};
  bool m_xyPlots{false};

  std::string  m_title;
  unsigned int m_titleSize{0};

  int    m_totTrack{0};
  double m_nbGhost{0.};
  double m_nbGhostHit{0.};
  int    m_totTrackTrigger{0}; ///< Total number of tracks processed
  int    m_totGhostTrigger{0};

  std::vector<std::string> m_name;       ///< Name of the sub-counters
  std::vector<bool>        m_writeHisto; ///< Make histograms for this container
  std::vector<double>      m_nbTrack;
  std::vector<double>      m_mcHits;  ///< Nb of MC hits on tracks
  std::vector<double>      m_foundOK; ///< Nb of correct hits
  std::vector<double>      m_wrong;   ///< Nb of wrong ones
  std::vector<double>      m_nbTrack3;
  std::vector<double>      m_mcHits3;  ///< Nb of MC hits on tracks with >= 3 UT hits
  std::vector<double>      m_foundOK3; ///< Nb of correct hits 3 UT hits
  std::vector<double>      m_wrong3;   ///< Nb of wrong ones 3 UT hits
};
#endif // PRUTCOUNTER_H
