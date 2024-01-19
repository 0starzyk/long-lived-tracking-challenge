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
#include "PrUTCounter.h"

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/LinksByKey.h"
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Event/MCProperty.h"
#include "Event/MCTrackInfo.h"
#include "Event/MCVertex.h"
#include "Event/Track.h"
#include "LoKi/IMCHybridFactory.h"
#include "LoKi/MCParticles.h"
#include "LoKi/Primitives.h"
#include "MCInterfaces/IMCReconstructible.h"
#include "MuonDet/DeMuonDetector.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include "GaudiAlg/GaudiHistoAlg.h"
#include "GaudiAlg/GaudiHistoTool.h"
#include "GaudiAlg/IHistoTool.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Consumer.h"

#include "boost/algorithm/string/replace.hpp"

/**
 *  Check the quality of the pattern recognition, by comparing to MC information
 *  Produces efficiency, ghost rate and clone rate numbers.
 *  Parameters:
 *   - [deprecated] Eta25Cut: Only consider particles with 2 < eta < 5? (default: false)
 *   - CheckNegEtaPlot: Check eta plotting range, plot negative values only if no eta25 cut was applied (default: false)
 *   - TriggerNumbers: Give numbers for p > 3GeV, pT > 500 MeV? (default: false)
 *     if selected long_fromB_P>3GeV_Pt>0.5GeV cut is added to each track container
 *   - VetoElectrons: Take electrons into account in numbers? (default: true)
 *   - WriteTexOutput: Writes the statistics table to a TeX file (default: false)
 *     which is dumped to the location specified in TexOutputFolder
 *   - MyCuts: selection cuts to be applied (default: empty)
 *   - WriteHistos: whether to plot histograms via IHistoTool. Values are -1 (default, no histograms), 1 (histograms), 2
 * (mote histograms : expectedHits, docaz, PVz, EtaP, EtaPhi, efficiency maps @z=9000mm XYZ9000 and @z=2485mm XYZ2485)
 *
 * This class is templated and is available in python with PrCounter and PrUTCounter as template arguments
 * The respective names to import are PrChecker and PrUTHitChecker
 *
 * Typical usage :
 *
 * @code
 *   from Configurables import PrChecker
 *   mychecker = PrChecker("PrCheckerVelo",
 *                          Title="Velo",
 *                          Tracks = "Rec/Track/Velo",
 *                          TriggerNumbers=True,
 *                          MyCuts = { "01_velo" : "isVelo",
 *                                     "02_long" : "isLong",
 *                                     "03_long>5GeV" : "isLong & over5" } ))
 *
 *   from Configurables import LoKi__Hybrid__MCTool
 *   myFactory = LoKi__Hybrid__MCTool("MCHybridFactory")
 *   myFactory.Modules = [ "LoKiMC.decorators" ]
 *   mychecker.addTool( myFactory )
 *  @endcode
 *
 *  As a default selection cuts of old PrChecker are used. The following cuts are predefined:
 *  - is(Not)Long, is(Not)Velo, is(Not)Down, is(Not)Up, is(Not)UT, is(Not)Seed,
 *  - fromB, fromD, BOrDMother, fromKsFromB, strange
 *  - is(Not)Electron, eta25, over5, trigger
 *
 *  and LoKi syntax (LoKi::MCParticles) can be used for kinematical cuts, e.g. (MCPT> 2300), here the '()' are
 * essential.
 *
 *  NB: If you care about the implementation: The cut-strings are converted into two types of functors:
 *  - LoKi-type functors (hence all LoKi::MCParticles cuts work)
 *  - and custom-defined ones, mostly for type of reconstructibility and daughter-criteria (like 'isNotLong')
 *  where in the end all functors are evaluated on each MCParticle for each track container to define the
 * reconstructibility. If a MCParticle is actually reconstructed is checked. A large part of the code just deals with
 * the conversion of the strings into functors.
 *
 */

template <typename InternalCounter>
class PrCheckerAlgorithm
    : public LHCb::Algorithm::Consumer<void( LHCb::Track::Range const&, LHCb::MCParticles const&,
                                             LHCb::MCVertices const&, LHCb::MCProperty const&, LHCb::LinksByKey const&,
                                             LHCb::LinksByKey const&, DetectorElement const& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DetectorElement>> {
public:
  PrCheckerAlgorithm( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"Tracks", ""}, KeyValue{"MCParticleInput", LHCb::MCParticleLocation::Default},
                   KeyValue{"MCVerticesInput", LHCb::MCVertexLocation::Default},
                   KeyValue{"MCPropertyInput", LHCb::MCPropertyLocation::TrackInfo}, KeyValue{"Links", ""},
                   KeyValue{"LinkTableLocation", "Link/Pr/LHCbID"},
                   KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}} ) {}

  /// Algorithm initialization
  StatusCode initialize() override;
  /// Algorithm execution
  void operator()( LHCb::Track::Range const&, LHCb::MCParticles const&, LHCb::MCVertices const&,
                   LHCb::MCProperty const&, LHCb::LinksByKey const&, LHCb::LinksByKey const&,
                   DetectorElement const& ) const override;
  /// Algorithm finalization
  StatusCode finalize() override;

private:
  Gaudi::Property<std::string>        m_title{this, "Title", ""};
  Gaudi::Property<unsigned int>       m_firstNVeloHits{this, "FirstNVeloHits", 3};
  Gaudi::Property<unsigned int>       m_hitTypesToCheck{this, "HitTypesToCheck", 0};
  Gaudi::Property<LHCb::Track::Types> m_trackType{this, "TrackType", LHCb::Track::Types::Unknown};

  // The counters are not at all thread safe for the moment
  // So we protect their use by a mutex. Not optimal but MC checking
  // does not need to be absolutely fast
  mutable InternalCounter m_counter;
  mutable std::mutex      m_counterMutex;

  // -- histograming options
  Gaudi::Property<int> m_writeHistos{this, "WriteHistos", -1};

  Gaudi::Property<bool> m_eta25cut{this, "Eta25Cut", false}; // to be deprecated, cuts on track properties should be
                                                             // done in filters before the algorithm.
  Gaudi::Property<bool> m_checkNegEtaPlot{
      this, "CheckNegEtaPlot", false,
      "eta plotting range check : plot negative values only if no eta25 cut was applied"};
  Gaudi::Property<bool>  m_triggerNumbers{this, "TriggerNumbers", false};
  Gaudi::Property<bool>  m_vetoElectrons{this, "VetoElectrons", true};
  Gaudi::Property<bool>  m_xyPlots{this, "XYPlots", false};
  Gaudi::Property<float> m_ghostProbCut{this, "GhostProbCut", 1.0}; // to be deprecated, cuts on track properties should
                                                                    // be done in filters before the algorithm.

  Gaudi::Property<bool>        m_writetexfile{this, "WriteTexOutput", false};
  Gaudi::Property<std::string> m_texfilename{this, "TexOutputName", "efficiencies"};
  Gaudi::Property<std::string> m_texfolder{this, "TexOutputFolder", ""};

  enum recAs {
    isLong,
    isNotLong,
    isDown,
    isNotDown,
    isUp,
    isNotUp,
    isVelo,
    isNotVelo,
    isUT,
    isNotUT,
    isSeed,
    isNotSeed,
    strange,
    fromB,
    fromD,
    fromKsFromB,
    isElectron,
    isNotElectron,
    BOrDMother,
    PairProd,
    isDecay,
    fromHI,
    fromPV,
    hasVeloOverlap,
    hasVeloCrossingSide,
    muonHitsInAllStations,
    muonHitsInAtLeastTwoStations,
    isMuon,
    isPion,
    fromSignal
  };
  static constexpr int size_recAs = fromSignal + 1;
  std::string_view     toString( recAs r ) {
    constexpr std::array<std::pair<recAs, std::string_view>, size_recAs> recAs_labels = {
        {{isLong, "isLong"},
         {isNotLong, "isNotLong"},
         {isDown, "isDown"},
         {isNotDown, "isNotDown"},
         {isUp, "isUp"},
         {isNotUp, "isNotUp"},
         {isVelo, "isVelo"},
         {isNotVelo, "isNotVelo"},
         {isUT, "isUT"},
         {isNotUT, "isNotUT"},
         {isSeed, "isSeed"},
         {isNotSeed, "isNotSeed"},
         {strange, "strange"},
         {fromB, "fromB"},
         {fromD, "fromD"},
         {fromKsFromB, "fromKsFromB"},
         {isElectron, "isElectron"},
         {isNotElectron, "isNotElectron"},
         {BOrDMother, "BOrDMother"},
         {PairProd, "PairProd"},
         {isDecay, "isDecay"},
         {fromHI, "fromHI"},
         {fromPV, "fromPV"},
         {hasVeloOverlap, "hasVeloOverlap"},
         {hasVeloCrossingSide, "hasVeloCrossingSide"},
         {muonHitsInAllStations, "muonHitsInAllStations"},
         {muonHitsInAtLeastTwoStations, "muonHitsInAtLeastTwoStations"},
         {isMuon, "isMuon"},
         {isPion, "isPion"},
         {fromSignal, "fromSignal"}}};
    auto i = std::find_if( recAs_labels.begin(), recAs_labels.end(), [r]( const auto& p ) { return p.first == r; } );
    if ( i != recAs_labels.end() ) return i->second;
    throw std::runtime_error( "bad recAs value" );
  }

  // convert strings to normal cuts ==> called m_otherCuts (without LoKi Hybrid factory)
  /**
   *  Predefined selection cuts: it converts strings to normal cuts, used by addOtherCuts
   */
  class isTrack {
    recAs m_kind;

  public:
    isTrack( recAs kind ) { m_kind = kind; };
    /// Functor that checks if the MCParticle fulfills certain criteria, e.g. reco'ble as long track, B daughter, ...
    bool operator()( LHCb::MCParticle* mcp, MCTrackInfo* mcInfo, std::vector<LHCb::LHCbID> const& lhcbIds ) const {
      switch ( m_kind ) {
      case isLong:
        return mcInfo->hasVeloAndT( mcp );
      case isNotLong:
        return !mcInfo->hasVeloAndT( mcp );
      case isDown:
        return mcInfo->hasT( mcp ) && mcInfo->hasUT( mcp );
      case isNotDown:
        return !( mcInfo->hasT( mcp ) && mcInfo->hasUT( mcp ) );
      case isUp:
        return mcInfo->hasVelo( mcp ) && mcInfo->hasUT( mcp );
      case isNotUp:
        return !( mcInfo->hasVelo( mcp ) && mcInfo->hasUT( mcp ) );
      case isVelo:
        return mcInfo->hasVelo( mcp );
      case isNotVelo:
        return !mcInfo->hasVelo( mcp );
      case isSeed:
        return mcInfo->hasT( mcp );
      case isNotSeed:
        return !mcInfo->hasT( mcp );
      case isUT:
        return mcInfo->hasUT( mcp );
      case isNotUT:
        return !mcInfo->hasUT( mcp );
      case isElectron:
        return std::abs( mcp->particleID().pid() ) == 11;
      case isNotElectron:
        return std::abs( mcp->particleID().pid() ) != 11;
      case hasVeloOverlap: {
        std::vector<LHCb::LHCbID> vpids;
        for ( auto const& lhcbid : lhcbIds ) {
          if ( lhcbid.isVP() ) {
            auto vp_sensor = lhcbid.vpID().sensor();
            auto vp_module = lhcbid.vpID().module();
            for ( auto const& id2 : vpids ) {
              // if 2 ids are in the same module but different sensors, there is an overlap
              if ( vp_module == id2.vpID().module() && vp_sensor != id2.vpID().sensor() ) return true;
            }
            vpids.push_back( lhcbid );
          }
        }
        return false;
      }
      case hasVeloCrossingSide: {
        std::vector<LHCb::LHCbID> vpids;
        for ( auto const& lhcbid : lhcbIds ) {
          if ( lhcbid.isVP() ) {
            auto vp_station = lhcbid.vpID().station();
            auto vp_module  = lhcbid.vpID().module();
            for ( auto const& id2 : vpids ) {
              // if 2 ids are in the same station but different modules, the particle is crossing side
              if ( vp_station == id2.vpID().station() && vp_module != id2.vpID().module() ) return true;
            }
            vpids.push_back( lhcbid );
          }
        }
        return false;
      }
      case muonHitsInAllStations: {
        // Check that lhcbIds contains at least one muon hit in each muon station
        constexpr auto          n_stations = 4;
        std::bitset<n_stations> seen;
        for ( auto const& lhcbid : lhcbIds ) {
          if ( lhcbid.isMuon() ) { seen[lhcbid.muonID().station()] = true; }
        }
        return seen.all();
      }
      case muonHitsInAtLeastTwoStations: {
        // Check that lhcbIds contains at least one muon hit in at least two muon stations
        constexpr auto          n_stations = 4;
        std::bitset<n_stations> seen;
        for ( auto const& lhcbid : lhcbIds ) {
          if ( lhcbid.isMuon() ) { seen[lhcbid.muonID().station()] = true; }
        }
        return ( seen.count() >= 2 );
      }
      case isMuon:
        return std::abs( mcp->particleID().pid() ) == 13;
      case isPion:
        return std::abs( mcp->particleID().pid() ) == 211;
      case fromSignal:
        return mcp->fromSignal();
      default:;
      }

      if ( !mcp->originVertex() ) return false;

      const LHCb::MCParticle* mother = mcp->originVertex()->mother();
      if ( mother ) {
        if ( mother->originVertex() ) {
          double rOrigin = mother->originVertex()->position().rho();
          if ( std::abs( rOrigin ) < 5. ) {
            int pid = std::abs( mother->particleID().pid() );
            // -- MCParticle is coming from a strang particle
            if ( m_kind == strange ) {
              return ( 130 == pid ||  // K0L
                       310 == pid ||  // K0S
                       321 == pid ||  // K+
                       3122 == pid || // Lambda
                       3222 == pid || // Sigma+
                       3212 == pid || // Sigma0
                       3112 == pid || // Sigma-
                       3322 == pid || // Xsi0
                       3312 == pid || // Xsi-
                       3334 == pid    // Omega-
              );
            }
            // -- It's a Kshort from a b Hadron
            if ( m_kind == fromKsFromB ) {
              auto gmom = mother->originVertex()->mother();
              return gmom && 310 == pid && 2 == mcp->originVertex()->products().size() &&
                     gmom->particleID().hasBottom() &&
                     ( gmom->particleID().isMeson() || gmom->particleID().isBaryon() );
            }
          }
        }
      }
      // -- It's a daughter of a B or D hadron
      bool motherB = false;
      bool motherD = false;
      for ( ; mother; mother = mother->originVertex()->mother() ) {
        if ( mother->particleID().isMeson() || mother->particleID().isBaryon() ) {
          if ( mother->particleID().hasBottom() ) motherB = true;
          if ( mother->particleID().hasCharm() ) motherD = true;
        }
      }
      switch ( m_kind ) {
      case fromD:
        return motherD;
      case fromB:
        return motherB;
      case BOrDMother:
        return motherD || motherB;
      default:;
      }
      // -- It's from a decay, from gamma->ee pair production or from Hadronic Interaction.
      // -- isDecay includes both DecayVertex and OscillatedAndDecay
      auto t = mcp->originVertex()->type();
      switch ( m_kind ) {
      case isDecay:
        return t == LHCb::MCVertex::MCVertexType::DecayVertex || t == LHCb::MCVertex::MCVertexType::OscillatedAndDecay;
      case PairProd:
        return t == LHCb::MCVertex::MCVertexType::PairProduction;
      case fromHI:
        return t == LHCb::MCVertex::MCVertexType::HadronicInteraction;
      case fromPV:
        return t == LHCb::MCVertex::MCVertexType::ppCollision;
      default:;
      }

      return false;
    }
  };

  /**
   *  Class that adds selection cuts defined in isTrack to cuts
   */
  class addOtherCuts {
    std::vector<isTrack> m_cuts;

  public:
    void addCut( recAs cat ) { m_cuts.emplace_back( cat ); }

    /// Functor that evaluates all 'isTrack' cuts
    bool operator()( LHCb::MCParticle* mcp, MCTrackInfo* mcInfo, std::vector<LHCb::LHCbID> const& lhcbIds ) const {
      return std::all_of( m_cuts.begin(), m_cuts.end(),
                          [&]( const auto& cut ) { return cut( mcp, mcInfo, lhcbIds ); } );
    }
  };

  // maps for each track container with {cut name,selection cut}
  Gaudi::Property<std::map<std::string, std::string>> m_map{this, "MyCuts", {}};

  /** @brief makes vector of second elements of DefaultCutMap --> needed as input for m_Cuts */
  std::vector<std::string> getMyCut( std::map<std::string, std::string> myCutMap ) {
    std::vector<std::string> dummy;
    std::transform( myCutMap.begin(), myCutMap.end(), std::back_inserter( dummy ),
                    []( const auto& p ) { return p.second; } );
    return dummy;
  }

  ToolHandle<IHistoTool>             m_histoTool{this, "HistoTool", "HistoTool/PrCheckerHistos"};
  ToolHandle<LoKi::IMCHybridFactory> m_factory{this, "LoKiFactory", "LoKi::Hybrid::MCTool/MCHybridFactory:PUBLIC"};
  ToolHandle<ITrackExtrapolator>     m_extrapolator{this, "TrackMasterExtrapolator", "TrackMasterExtrapolator"};
  std::vector<LoKi::Types::MCCut>    m_MCCuts;
  std::vector<addOtherCuts>          m_MCCuts2;
};

template <typename InternalCounter>
StatusCode PrCheckerAlgorithm<InternalCounter>::initialize() {
  StatusCode sc = Consumer::initialize();
  if ( sc.isFailure() ) { return sc; }

  static const std::string histoDir = "Track/";
  if ( "" == histoTopDir() ) setHistoTopDir( histoDir );
  m_histoTool.retrieve().ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ ); // needs to be done for next line to
                                                                                  // work
  GaudiHistoTool* ghtool = dynamic_cast<GaudiHistoTool*>( m_histoTool.get() );

  // -- catch the possible failure of the dynamic cast
  if ( !ghtool ) {
    error() << "Dynamic cast of Gaudi Histogramming Tool failed!" << endmsg;
    return StatusCode::FAILURE;
  }

  ghtool->setHistoDir( histoDir + name() );

  m_counter.setTitle( m_title.value() );
  m_counter.setFirstNVeloHits( m_firstNVeloHits );
  m_counter.setWriteHistos( m_writeHistos );
  m_counter.setHitTypesToCheck( m_hitTypesToCheck );
  m_counter.setTrackType( m_trackType );
  m_counter.setXYPlots( m_xyPlots );

  for ( auto pair : m_map ) {
    if ( m_checkNegEtaPlot ) {
      // define eta plotting range: plot negative values only if no eta25 cut was applied
      const std::string etaString( "eta25" );
      std::size_t       found = pair.second.find( etaString );
      if ( !m_eta25cut && found == std::string::npos ) {
        m_counter.addSelection( pair.first, true, true );
      } else {
        m_counter.addSelection( pair.first, true );
      }
    } else {
      m_counter.addSelection( pair.first, true );
    }
  }

  m_counter.setUseEta25Cut( m_eta25cut );
  if ( m_eta25cut ) {
    warning() << "Property Eta25Cut to be deprecated. Please use dedicated algorithm to filter on track properties. "
                 "For example see in PrUpgradeChecking.py"
              << endmsg;
  }
  m_counter.setGhostProbCut( m_ghostProbCut );
  if ( m_ghostProbCut != 1.0 ) {
    warning() << "Property GhostProbCut to be deprecated. Please use dedicated algorithm to filter on track "
                 "properties. For example see in PrUpgradeChecking.py"
              << endmsg;
  }

  m_counter.setTriggerNumbers( m_triggerNumbers );
  if ( m_writetexfile.value() ) m_counter.setTeXName( m_texfolder, m_texfilename );

  // -- convert all strings into functors
  for ( std::string cutString : getMyCut( m_map ) ) { // loop over 2nd element of Cuts = strings of cuts

    m_MCCuts2.emplace_back();

    // flag to circumvent veto of electrons for selections that only look at elctrons
    bool ExplicitlyKeepElectrons = cutString.find( "isElectron" ) != std::string::npos;
    // -- extract aliases from cutString and replace with 'MCTRUE'
    for ( int i = 0; i != size_recAs; ++i ) {
      recAs alias = static_cast<recAs>( i );

      std::size_t found = cutString.find( toString( alias ) );

      if ( found != std::string::npos ) { // if found then
        m_MCCuts2.back().addCut( alias ); // other components are already filled //add this
                                          // category of cuts to addOtherCuts()
        cutString.replace( found, toString( alias ).length(),
                           "MCTRUE" ); // replace found at position found, with length of string
                                       // to replace, replace it with string "" (Loki Cut)
      }
    }

    // -- Veto electrons or not
    if ( m_vetoElectrons && !ExplicitlyKeepElectrons ) {
      m_MCCuts2.back().addCut( isNotElectron );
      boost::replace_first( cutString, toString( isNotElectron ), "MCTRUE" );
    }

    // -- LoKi cuts: define aliases for better use
    const std::string etaString( "eta25" );
    std::size_t       found = cutString.find( etaString );
    if ( m_eta25cut && found == std::string::npos ) {
      cutString.append( " & eta25" );
      found = cutString.find( etaString );
    }
    if ( found != std::string::npos ) {
      cutString.replace( found, etaString.length(), "(MCETA > 2.0) & (MCETA < 5.0)" );
    }

    boost::replace_first( cutString, "over5", "(MCP > 5000)" );
    boost::replace_first( cutString, "trigger", "(MCP>3000) & (MCPT>500)" );
    // ---------------------------------------------------------------------------------

    LoKi::Types::MCCut tmp = LoKi::BasicFunctors<const LHCb::MCParticle*>::BooleanConstant( false ); //
    sc                     = m_factory->get( cutString, tmp );
    m_MCCuts.push_back( tmp );
    if ( sc.isFailure() ) { return Error( "Error from IMCHybridFactory", sc ); } // RETURN
  }

  return StatusCode::SUCCESS;
}

template <typename InternalCounter>
void PrCheckerAlgorithm<InternalCounter>::
     operator()( LHCb::Track::Range const& tracks, LHCb::MCParticles const& mcParts, LHCb::MCVertices const& mcVert,
            LHCb::MCProperty const& flags, LHCb::LinksByKey const& tr2McLink, LHCb::LinksByKey const& mc2IdLink,
            DetectorElement const& geometry ) const {
  // in debug mode, check consistency of the inputs
  assert( tr2McLink.sourceClassID() == LHCb::Track::classID() &&
          "Incompatible link table in PrCheckerAlgorithm. Source should be Track" );
  assert( tr2McLink.targetClassID() == LHCb::MCParticle::classID() &&
          "Incompatible link table in PrCheckerAlgorithm. Target should be McParticle" );

  MCTrackInfo  trackInfo = {flags};
  unsigned int nPrim     = std::count_if( mcVert.begin(), mcVert.end(), [&]( const auto& vertex ) {
    if ( !vertex->isPrimary() ) return false;
    int nbVisible = std::count_if( mcParts.begin(), mcParts.end(), [&]( const auto& part ) {
      return part->primaryVertex() == vertex && trackInfo.hasVelo( part );
    } );
    return nbVisible > 4;
  } );

  m_counter.initEvent( m_histoTool.get(), m_extrapolator.get(), nPrim, tracks, tr2McLink, geometry );

  //== Build a table (vector of map) of Track -> weigth per MCParticle, indexed by MCParticle key.
  std::vector<std::map<const LHCb::Track*, double>> tracksForParticle;
  tr2McLink.applyToAllLinks( [&tracksForParticle, &tracks]( int trackKey, unsigned int mcPartKey, float weight ) {
    auto track = std::find_if( tracks.begin(), tracks.end(), [&trackKey]( auto t ) { return t->key() == trackKey; } );
    if ( track != tracks.end() ) {
      if ( tracksForParticle.size() <= mcPartKey ) { tracksForParticle.resize( mcPartKey + 1 ); }
      tracksForParticle[mcPartKey][*track] += weight;
    }
  } );

  //== Build a table (vector of vectors) of LHCbID per MCParticle, indexed by MCParticle key.
  std::vector<std::vector<LHCb::LHCbID>> idsForParticle;
  mc2IdLink.applyToAllLinks( [&idsForParticle]( unsigned int id, unsigned int mcPartKey, float ) {
    if ( idsForParticle.size() <= mcPartKey ) { idsForParticle.resize( mcPartKey + 1 ); }
    idsForParticle[mcPartKey].emplace_back( id );
  } );

  std::vector<LHCb::LHCbID>            noids;
  std::map<const LHCb::Track*, double> noTracksList;
  for ( const auto part : mcParts ) {
    if ( 0 == trackInfo.fullInfo( part ) ) continue;
    // get the LHCbIDs that go with the MCParticle
    auto& ids = idsForParticle.size() > (unsigned int)part->key() ? idsForParticle[part->key()] : noids;
    // cuts
    std::vector<bool> flags;
    for ( unsigned int i = 0; i < m_MCCuts.size(); ++i ) {
      flags.push_back( m_MCCuts[i]( part ) && m_MCCuts2[i]( part, &trackInfo, ids ) );
    }
    try {
      m_counter.countAndPlot( m_histoTool.get(), m_extrapolator.get(), part, flags, ids, nPrim,
                              tracksForParticle.size() > (unsigned int)part->key() ? tracksForParticle[part->key()]
                                                                                   : noTracksList,
                              geometry );
    } catch ( const std::string& msg ) {
      Warning( msg ).ignore();
      if ( msgLevel( MSG::DEBUG ) ) {
        debug() << "... Flag size " << flags.size() << " >  " << m_counter.title().size() << " declared selections"
                << endmsg;
      }
    }
  }
}

template <typename InternalCounter>
StatusCode PrCheckerAlgorithm<InternalCounter>::finalize() {
  info() << "Results" << endmsg;
  m_counter.printStatistics( info(), inputLocation<0>() );
  return Consumer::finalize(); // must be called after all other actions
}

// Easier with typedefs to avoid strange syntax in Macros
typedef PrCheckerAlgorithm<PrTrackCounter> PrTrackChecker;
typedef PrCheckerAlgorithm<PrUTCounter>    PrUTHitChecker;

DECLARE_COMPONENT_WITH_ID( PrTrackChecker, "PrTrackChecker" )
DECLARE_COMPONENT_WITH_ID( PrUTHitChecker, "PrUTHitChecker" )
