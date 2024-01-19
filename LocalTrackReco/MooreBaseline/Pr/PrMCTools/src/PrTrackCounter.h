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
#ifndef PRTRACKCOUNTER_H
#define PRTRACKCOUNTER_H 1

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/Track.h"

#include "GaudiAlg/GaudiTool.h"
#include "GaudiAlg/IHistoTool.h"

#include "TrackInterfaces/ITrackExtrapolator.h"

#include <memory>

/** @class PrTrackCounter PrTrackCounter.h
 *  This is a counter for track efficiency measurement.
 *
 *  @author Olivier Callot
 *  @date   2005-06-10
 *  Modified by Wenbin Qian for the VP Pat Efficiency
 *  @adpate to A-Team framework 2007-08-20 SHM
 */
class PrTrackCounter final {
public:
  void initEvent( IHistoTool const* htool, ITrackExtrapolator const* extrapolator, int nPV,
                  LHCb::Track::Range const& tracks, LHCb::LinksByKey const& links, IGeometryInfo const& geometry );

  void countAndPlot( IHistoTool const* htool, ITrackExtrapolator const* extrapolator, LHCb::MCParticle const* part,
                     std::vector<bool> flags, std::vector<LHCb::LHCbID>& ids, int nPV,
                     std::map<const LHCb::Track*, double> const& trackList, IGeometryInfo const& geometry );

  void addSelection( std::string name, bool writeHisto, bool plotNegEta = false );

  void printStatistics( MsgStream&, std::string location );

  std::string title() { return m_title; }
  void        setTitle( std::string& title ) { m_title = title; }
  void        setFirstNVeloHits( unsigned int firstNVeloHits ) { m_firstNVeloHits = firstNVeloHits; }
  void        setWriteHistos( int write ) { m_writeHistos = write; }
  void        setUseEta25Cut( bool cut ) {
    m_eta25cut = cut;
  } // to be deprecated, cuts on track properties should be done in filters before the algorithm.
  void setTriggerNumbers( bool numbers ) { m_triggerNumbers = numbers; }
  void setXYPlots( bool xyPlots ) { m_xyPlots = xyPlots; };
  void setHitTypesToCheck( int data ) { m_hitTypesToCheck = data; }
  void setTrackType( LHCb::Track::Types type ) { m_trackType = type; }
  void setGhostProbCut( double cutVal ) {
    m_ghostProbCut = cutVal;
  } // to be deprecated, cuts on track properties should be done in filters before the algorithm.
  void setTeXName( const std::string& directory, const std::string& name ) {
    m_writetex   = true;
    m_texoutdir  = directory;
    m_texoutname = name;
  }

  // Has to mirror the the dictionary HitType in PrUpgradeChecking.py
  enum HitType { Unspecified = 0, VP = 3, UT = 4, FT = 8 };

private:
  int         m_writeHistos{-1};
  bool        m_eta25cut{false};
  bool        m_triggerNumbers{false};
  bool        m_xyPlots{false};
  int         m_hitTypesToCheck{HitType::Unspecified};
  double      m_ghostProbCut{999};
  bool        m_writetex{false};
  std::string m_texoutname;
  std::string m_texoutdir;

  LHCb::Track::Types m_trackType{LHCb::Track::Types::Unknown};

  std::string  m_title;
  unsigned int m_titleSize{0};
  unsigned int m_firstNVeloHits{3};

  // total variables
  int    m_totTrack{0};        ///< Total number of tracks processed
  int    m_totGhost{0};        ///< Total number of ghosts
  int    m_totTrackTrigger{0}; ///< Total number of tracks processed
  int    m_totGhostTrigger{0}; ///< Total number of ghosts
  double m_fracGhost{0.};
  double m_nEvent{0.};

  std::vector<std::string> m_name;       ///< Name of the sub-counters
  std::vector<bool>        m_writeHisto; ///< Make histograms for this container
  std::vector<bool>        m_plotNegEta; ///< plot negative eta values (for example for Velo tracks)
  std::vector<int>         m_wanted;     ///< Nb MC tracks measurable.
  std::vector<int>         m_counted;    ///< counters for statistics
  std::vector<unsigned>    m_velofirstcounter;
  std::vector<unsigned>    m_velolastcounter;
  std::vector<int>         m_clone;               ///< counters for clones
  std::vector<double>      m_purity;              ///< Sum of purity (linker weight)
  std::vector<double>      m_hitEff;              ///< Sum of hitEfficiency
  std::vector<double>      m_hitEffFirstVeloHits; ///< Sum of hitEfficiency
  std::vector<double>      m_hitEffLastVeloHits;  ///< Sum of hitEfficiency
};

#endif // PRTRACKCOUNTER_H
