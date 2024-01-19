
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

/** @class PrMultiplicityChecker PrMultiplicityChecker.cpp
 *  Check the quality of the pattern recognition, by comparing to MC information
 *  Produces efficiency and ghost rate vs multiplicity.
 * This class is available in python with PrMultiplicityChecker as template arguments
 * The respective name to import is PrMultiplicityChecker
 *
 *  As a default selection cuts of old PrChecker are used. The following cuts are predefined:
 *  - isVelo, isUp, isLong, isSeed,
 */

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/Track.h"

#include "Event/MCTrackInfo.h"

#include "Event/PrHits.h"
#include "Event/PrVeloTracks.h"
#include "PrKernel/PrFTHitHandler.h"
#include "PrKernel/UTHitHandler.h"

#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiAlg/GaudiHistoTool.h"
#include "GaudiAlg/IHistoTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Consumer.h"

#include "LoKi/IMCHybridFactory.h"
#include "LoKi/MCParticles.h"
#include "LoKi/Primitives.h"

#include "boost/algorithm/string/replace.hpp"
#include "fmt/format.h"

namespace {
  using KeyValue = std::pair<std::string, std::string>;
  namespace VP   = LHCb::Pr::VP;
} // namespace

class PrMultiplicityChecker
    : public LHCb::Algorithm::Consumer<void( const LHCb::Track::Range&, const LHCb::LinksByKey&,
                                             const LHCb::MCParticles&, const LHCb::MCVertices&, const LHCb::MCProperty&,
                                             const VP::Hits&, const LHCb::Pr::Velo::Tracks&, const UT::HitHandler&,
                                             const PrFTHitHandler<PrHit>& ),
                                       Gaudi::Functional::Traits::BaseClass_t<GaudiHistoAlg>> {
public:
  PrMultiplicityChecker( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"Tracks", ""}, KeyValue{"Links", ""},
                   KeyValue{"MCParticleInput", LHCb::MCParticleLocation::Default},
                   KeyValue{"MCVerticesInput", LHCb::MCVertexLocation::Default},
                   KeyValue{"MCPropertyInput", LHCb::MCPropertyLocation::TrackInfo}, KeyValue{"Velo_Hits", ""},
                   KeyValue{"Velo_Tracks", ""}, KeyValue{"UT_Hits", UTInfo::HitLocation},
                   KeyValue{"FT_Hits", PrFTInfo::FTHitsLocation}} ){};

  void operator()( const LHCb::Track::Range& tracks, const LHCb::LinksByKey& links, const LHCb::MCParticles& mcParts,
                   const LHCb::MCVertices& mcVert, const LHCb::MCProperty& flags, const VP::Hits& Velo_Hits,
                   const LHCb::Pr::Velo::Tracks& Velo_Tracks, const UT::HitHandler& UT_Hits,
                   const PrFTHitHandler<PrHit>& FT_Hits ) const override {

    // in debug mode, check consistency of the inputs
    assert( links.sourceClassID() == LHCb::Track::classID() &&
            "Incompatible link table in PrMultiplicityChecker. Source should be Track" );
    assert( links.targetClassID() == LHCb::MCParticle::classID() &&
            "Incompatible link table in PrMultiplicityChecker. Target should be McParticle" );

    auto         trackInfo = MCTrackInfo{flags};
    unsigned int nPrim     = std::count_if( mcVert.begin(), mcVert.end(), [&]( const auto& vertex ) {
      if ( !vertex->isPrimary() ) return false;
      int nbVisible = std::count_if( mcParts.begin(), mcParts.end(), [&]( const auto& part ) {
        return part->primaryVertex() == vertex && trackInfo.hasVelo( part );
      } );
      return nbVisible > 4;
    } );

    this->initEvent( m_histoTool.get(), nPrim, tracks, links, Velo_Hits, Velo_Tracks, UT_Hits, FT_Hits );

    //== Build a table (vector of map) of Track -> weigth per MCParticle, indexed by MCParticle key.
    std::vector<std::map<const LHCb::Track*, double>> tracksForParticle;
    links.applyToAllLinks( [&tracksForParticle, &tracks]( int trackKey, unsigned int mcPartKey, float weight ) {
      auto track = std::find_if( tracks.begin(), tracks.end(), [&trackKey]( auto t ) { return t->key() == trackKey; } );
      if ( track != tracks.end() ) {
        if ( tracksForParticle.size() <= mcPartKey ) { tracksForParticle.resize( mcPartKey + 1 ); }
        tracksForParticle[mcPartKey][*track] += weight;
      }
    } );

    std::vector<LHCb::LHCbID>            noids;
    std::map<const LHCb::Track*, double> noTracksList;

    for ( const auto part : mcParts ) {
      if ( 0 == trackInfo.fullInfo( part ) ) continue;

      // make list of tags of reconstructible MC particle
      std::bitset<4> type_tags{};                    //  list of tags of reconstructible MC particle
      type_tags.set( 0, trackInfo.hasVelo( part ) ); // velo
      type_tags.set( 1, trackInfo.hasVelo( part ) && trackInfo.hasUT( part ) ); // upstream
      type_tags.set( 2, trackInfo.hasVeloAndT( part ) );                        // long
      type_tags.set( 3, trackInfo.hasT( part ) );                               // seed
      if ( type_tags.none() ) continue;

      try {
        this->countAndPlot( m_histoTool.get(), part, nPrim, type_tags,
                            tracksForParticle.size() > (unsigned int)part->key() ? tracksForParticle[part->key()]
                                                                                 : noTracksList,
                            Velo_Hits, Velo_Tracks, UT_Hits, FT_Hits );
      } catch ( const std::string& msg ) { Warning( msg ).ignore(); }
    }
  };

  StatusCode finalize() override {
    info() << "Results" << endmsg;
    printStatistics( info(), inputLocation<0>() );
    return Consumer::finalize();
  };

private:
  ToolHandle<IHistoTool>       m_histoTool{this, "HistoTool", "HistoTool/PrMultiplicityChecker"};
  Gaudi::Property<std::string> m_title{this, "Title", ""};

  mutable int    m_totTrack{0}; ///< Total number of tracks processed
  mutable int    m_totGhost{0}; ///< Total number of ghosts
  mutable double m_fracGhost{0};
  mutable double m_nEvent{0};

  static constexpr auto tag_name = std::array{"Velo", "Upstream", "Forward", "Seed"};

  void initEvent( const IHistoTool* htool, const int nPV, const LHCb::Track::Range& tracks,
                  const LHCb::LinksByKey& links, const VP::Hits& Velo_Hits, const LHCb::Pr::Velo::Tracks& Velo_Tracks,
                  const UT::HitHandler& UT_Hits, const PrFTHitHandler<PrHit>& FT_Hits ) const {

    double nbTracks = 0;
    double nbGhost  = 0;

    for ( const LHCb::Track* track : tracks ) {

      const auto& lhcbIDs = track->lhcbIDs();
      // Count the number of hits in each sub-detectors
      double n_hits_Velo  = 0.;
      double n_hits_UT    = 0.;
      double n_hits_Scifi = 0.;
      for ( const auto id : lhcbIDs ) {
        if ( id.isVP() )
          n_hits_Velo += 1;
        else if ( id.isUT() )
          n_hits_UT += 1;
        else if ( id.isFT() )
          n_hits_Scifi += 1;
      }

      std::vector<std::string> track_types{"_Total"};
      if ( !links.hasEntry( *track ) ) {
        nbGhost++;
        track_types.push_back( "_Ghosts" );
      }

      for ( auto const& tr_type : track_types ) {
        htool->plot1D( nPV, m_title + "/nPV" + tr_type, "nPV" + tr_type, -0.5, 20.5, 21 );
        htool->plot1D( track->momentum().Phi(), m_title + "/Phi" + tr_type, "Phi" + tr_type, -3.142, 3.142, 25 );
        htool->plot1D( track->pt(), m_title + "/Pt" + tr_type, "Pt" + tr_type, 0., 10000., 100 );
        htool->plot1D( track->p(), m_title + "/P" + tr_type, "P" + tr_type, 0., 100000., 100 );
        htool->plot1D( track->pseudoRapidity(), m_title + "/Eta" + tr_type, "Eta" + tr_type, 0., 7., 50 );
        htool->plot1D( lhcbIDs.size(), m_title + "/nHits_all" + tr_type, "nHits_all" + tr_type, -0.5, 35.5, 35 );
        htool->plot1D( n_hits_Velo, m_title + "/nHits_Velo" + tr_type, "nHits_Velo" + tr_type, -0.5, 25.5, 26 );
        htool->plot1D( n_hits_UT, m_title + "/nHits_UT" + tr_type, "nHits_UT" + tr_type, -0.5, 10.5, 11 );
        htool->plot1D( n_hits_Scifi, m_title + "/nHits_Scifi" + tr_type, "nHits_Scifi" + tr_type, -0.5, 15.5, 16 );
        htool->plot1D( FT_Hits.hits().size(), m_title + "/FTHits" + tr_type, "FTHits" + tr_type, 0., 30000., 300 );
        htool->plot1D( UT_Hits.nbHits(), m_title + "/UTHits" + tr_type, "UTHits" + tr_type, 0., 30000., 300 );
        htool->plot1D( Velo_Hits.size(), m_title + "/VeloHits" + tr_type, "VeloHits" + tr_type, 0., 30000., 300 );
        htool->plot1D( Velo_Tracks.size(), m_title + "/VeloTracks" + tr_type, "VeloTracks" + tr_type, 0., 3000., 300 );
      }
      nbTracks++;
    }
    m_totTrack += nbTracks;
    m_totGhost += nbGhost;

    double fracGhost = 0.;
    if ( 0 < nbTracks ) fracGhost = double( nbGhost ) / nbTracks;
    m_fracGhost += fracGhost;
    m_nEvent += 1.;

    htool->plot1D( nbTracks, m_title + "/nbTracks", "nbTracks", 0., 5000., 500 );
  }

  void countAndPlot( const IHistoTool* htool, const LHCb::MCParticle* part, const int nPV,
                     const std::bitset<4> type_tags, const std::map<const LHCb::Track*, double>& trackList,
                     const VP::Hits& Velo_Hits, const LHCb::Pr::Velo::Tracks& Velo_Tracks,
                     const UT::HitHandler& UT_Hits, const PrFTHitHandler<PrHit>& FT_Hits ) const {

    bool found = false;
    found      = !trackList.empty();

    double prodx = part->originVertex()->position().X();
    double prody = part->originVertex()->position().Y();
    double docaz = -100.0;
    if ( part->momentum().Pt() > 0.00001 ) {
      docaz =
          std::abs( 1. / part->momentum().Pt() * ( prodx * part->momentum().Py() - prody * part->momentum().Px() ) );
    }

    const LHCb::MCParticle* mother = part;
    while ( mother->mother() != NULL ) mother = mother->mother();
    double PVz = mother->originVertex()->position().Z();

    std::vector<std::string> particle_tags{"_reconstructible"};
    if ( found ) particle_tags.push_back( "_reconstructed" );
    for ( unsigned int k = 0; type_tags.size() > k; ++k ) {
      if ( type_tags[k] == 0 ) continue;
      for ( auto const& r_tag : particle_tags ) {
        htool->plot1D( nPV, fmt::format( "{}/{}_nPV{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_nPV{}", tag_name[k], r_tag ), -0.5, 20.5, 21 );
        htool->plot1D( part->momentum().Eta(), fmt::format( "{}/{}_Eta{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_Eta{}", tag_name[k], r_tag ), 0., 7., 50 );
        htool->plot1D( part->momentum().Phi(), fmt::format( "{}/{}_Phi{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_Phi{}", tag_name[k], r_tag ), -3.142, 3.142, 25 );
        htool->plot1D( part->momentum().Pt(), fmt::format( "{}/{}_Pt{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_Pt{}", tag_name[k], r_tag ), 0., 10000., 100 );
        htool->plot1D( part->momentum().P(), fmt::format( "{}/{}_P{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_P{}", tag_name[k], r_tag ), 0., 100000., 100 );
        htool->plot1D( part->originVertex()->position().Z(), fmt::format( "{}/{}_z{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_z{}", tag_name[k], r_tag ), -550., 200., 200 );
        htool->plot1D( docaz, fmt::format( "{}/{}_docaz{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_docaz{}", tag_name[k], r_tag ), 0., 10., 50 );
        htool->plot1D( PVz, fmt::format( "{}/{}_PVz{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_PVz{}", tag_name[k], r_tag ), -200., 200., 50 );
        htool->plot1D( FT_Hits.hits().size(), fmt::format( "{}/{}_FTHits{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_FTHits{}", tag_name[k], r_tag ), 0., 30000., 300 );
        htool->plot2D( FT_Hits.hits().size(), part->momentum().P(),
                       fmt::format( "{}/{}_FTHitsP{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_FTHitsP{}", tag_name[k], r_tag ), 0., 30000., 0., 100000., 20, 20 );
        htool->plot1D( UT_Hits.nbHits(), fmt::format( "{}/{}_UTHits{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_UTHits{}", tag_name[k], r_tag ), 0., 30000., 300 );
        htool->plot1D( Velo_Hits.size(), fmt::format( "{}/{}_VeloHits{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_VeloHits{}", tag_name[k], r_tag ), 0., 30000., 300 );
        htool->plot1D( Velo_Tracks.size(), fmt::format( "{}/{}_VeloTracks{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_VeloTracks{}", tag_name[k], r_tag ), 0., 3000., 30 );
        htool->plot2D( part->momentum().Eta(), part->momentum().P(),
                       fmt::format( "{}/{}_EtaP{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_EtaP{}", tag_name[k], r_tag ), 0., 7., 0., 100000., 20, 20 );
        htool->plot2D( part->momentum().Eta(), part->momentum().Phi(),
                       fmt::format( "{}/{}_EtaPhi{}", m_title, tag_name[k], r_tag ),
                       fmt::format( "{}_EtaPhi{}", tag_name[k], r_tag ), 0., 7., -3.142, 3.142, 20, 20 );
      }
    }
  }

  void printStatistics( MsgStream& info, std::string ) const {
    if ( 0 == m_nEvent ) return;
    double totT = m_totTrack + 0.00000000001;
    double frac = 100. * double( m_totGhost ) / totT;
    info << "**** " << m_title
         << format( "%7d tracks including        %7d ghosts [%5.2f %%], Event average %5.2f %% ****", m_totTrack,
                    m_totGhost, frac, 100. * m_fracGhost / m_nEvent )
         << endmsg;
  }
};

DECLARE_COMPONENT_WITH_ID( PrMultiplicityChecker, "PrMultiplicityChecker" )
