/*****************************************************************************\
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "AIDA/IHistogram1D.h"
#include "DetDesc/Condition.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/MCHeader.h"
#include "Event/MCParticle.h"
#include "Event/MCProperty.h"
#include "Event/MCTrackInfo.h"
#include "Event/MCVertex.h"
#include "Event/PrimaryVertices.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "Event/Track_v2.h"
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiAlg/Tuples.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiUtils/HistoStats.h"
#include "Kernel/STLExtensions.h"
#include "LHCbAlgs/Consumer.h"
#include "VPDet/DeVP.h"
#include "fmt/format.h"
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

//-----------------------------------------------------------------------------
// Implementation file for class : PrimaryVertexChecker
//-----------------------------------------------------------------------------

namespace {
  using Vertices = LHCb::Event::PV::PrimaryVertexContainer;

  std::string fmtRat( std::string_view name, Gaudi::Accumulators::BinomialCounter<> const& mcpv,
                      Gaudi::Accumulators::BinomialCounter<> const& eff,
                      Gaudi::Accumulators::BinomialCounter<> const& m_false ) {
    return fmt::format(
        "{:<25} : {:8} from {:8} ({:8}-{:<8}) [ {:5.2f} %], false {:4} from reco. {:8} ({:8}+{:<4}) [ {:5.2f} %] ",
        name, eff.nTrueEntries(), eff.nEntries(), mcpv.nEntries(), mcpv.nFalseEntries(), eff.efficiency() * 100,
        m_false.nTrueEntries(), m_false.nEntries(), m_false.nFalseEntries(), m_false.nTrueEntries(),
        m_false.efficiency() * 100 );
  }
  std::string fmtAvTracks( std::string_view name, Gaudi::Accumulators::AveragingCounter<> const& av_tracks,
                           Gaudi::Accumulators::AveragingCounter<> const& av_mcp ) {
    return fmt::format( "{:<25} : av. PV tracks: {:6.2f} [MC: {:6.2f}]", name, av_tracks.mean(), av_mcp.mean() );
  }

  double check_histogram( AIDA::IHistogram1D const* h, bool rms ) { return h ? ( rms ? h->rms() : h->mean() ) : 0.; }

  std::string fmtRes( std::string_view name, double x, double y, double z ) {
    return fmt::format( "{:<25} :  x: {:+5.3f}, y: {:+5.3f}, z: {:+5.3f}", name, x, y, z );
  }

  struct MCPVInfo {
    LHCb::MCVertex const* pMCPV             = {nullptr}; // pointer to MC PV
    int                   nRecTracks        = {0};       // number of reconstructed tracks from this MCPV
    int                   nRecBackTracks    = {0};       // number of reconstructed backward tracks
    int                   indexRecPVInfo    = {-1};      // index to reconstructed PVInfo (-1 if not reco)
    int                   nCorrectTracks    = {0};       // correct tracks belonging to reconstructed PV
    int                   multClosestMCPV   = {0};       // multiplicity of closest reconstructable MCPV
    double                distToClosestMCPV = {999999.}; // distance to closest reconstructible MCPV
    bool                  decayCharm        = {false};   // type of mother particle
    bool                  decayBeauty       = {false};
    bool                  decayStrange      = {false};
  };

  template <typename VERTEXTYPE>
  struct RecPVInfo {
    int               nTracks; // number of tracks
    double            chi2;
    double            nDoF;
    int               mother;
    Gaudi::XYZPoint   position;      // position
    Gaudi::XYZPoint   positionSigma; // position sigmas
    int               indexMCPVInfo; // index to MCPVInfo
    VERTEXTYPE const* pRECPV{nullptr};
  };

  MCPVInfo const* closestMCPV( LHCb::span<MCPVInfo const> rblemcpv, MCPVInfo const& mc ) {
    if ( rblemcpv.size() < 2 ) return nullptr;
    struct Result {
      MCPVInfo const* ptr  = nullptr;
      double          dist = std::numeric_limits<double>::max();
    };
    return std::accumulate( rblemcpv.begin(), rblemcpv.end(), Result{},
                            [&mc]( Result r, auto const& i ) {
                              if ( i.pMCPV == mc.pMCPV ) return r;
                              double dist = ( i.pMCPV->position() - mc.pMCPV->position() ).R();
                              return r.dist < dist ? r : Result{&i, dist};
                            } )
        .ptr;
  }

  template <typename RecPVInfoType>
  MCPVInfo const* closestMCPV( LHCb::span<MCPVInfo const> rblemcpv, RecPVInfoType const& rec ) {
    if ( rblemcpv.size() < 2 ) return nullptr;
    struct Result {
      MCPVInfo const* ptr  = nullptr;
      double          dist = std::numeric_limits<double>::max();
    };
    return std::accumulate( rblemcpv.begin(), rblemcpv.end(), Result{},
                            [&rec]( Result r, auto const& i ) {
                              double dist = ( i.pMCPV->position() - rec.pRECPV->position() ).R();
                              return r.dist < dist ? r : Result{&i, dist};
                            } )
        .ptr;
  }

  template <typename RecPVInfoType>
  void match_mc_vertex_by_distance( int ipv, LHCb::span<RecPVInfoType> rinfo, LHCb::span<MCPVInfo> mcpvvec ) {

    struct Result {
      MCPVInfo const* ptr  = nullptr;
      double          dist = std::numeric_limits<double>::max();
    };
    auto [ptr, mindist] = std::accumulate( mcpvvec.begin(), mcpvvec.end(), Result{},
                                           [z = rinfo[ipv].position.z()]( auto r, auto const& mcpv ) {
                                             if ( mcpv.indexRecPVInfo > -1 ) return r;
                                             double dist = std::abs( mcpv.pMCPV->position().z() - z );
                                             return r.dist < dist ? r : Result{&mcpv, dist};
                                           } );
    if ( ptr && mindist < 5.0 * rinfo[ipv].positionSigma.z() ) {
      auto indexmc                    = ptr - mcpvvec.data();
      rinfo[ipv].indexMCPVInfo        = indexmc;
      mcpvvec[indexmc].indexRecPVInfo = ipv;
    }
  }

  enum class recoAs {
    all,
    isolated,
    close,
    ntracks_low,
    ntracks_high,
    z_low,
    z_middle,
    z_high,
    beauty,
    charm,
    strange,
    other,
    first,
    second,
    third,
    fourth,
    fifth
  };

  constexpr auto All = std::array{
      recoAs::all,      recoAs::isolated, recoAs::close,  recoAs::ntracks_low, recoAs::ntracks_high, recoAs::z_low,
      recoAs::z_middle, recoAs::z_high,   recoAs::beauty, recoAs::charm,       recoAs::strange,      recoAs::other,
      recoAs::first,    recoAs::second,   recoAs::third,  recoAs::fourth,      recoAs::fifth};
  constexpr auto Part  = std::array{recoAs::all, recoAs::beauty, recoAs::charm, recoAs::strange, recoAs::other};
  constexpr auto Basic = std::array{recoAs::all,          recoAs::isolated, recoAs::close,    recoAs::ntracks_low,
                                    recoAs::ntracks_high, recoAs::z_low,    recoAs::z_middle, recoAs::z_high,
                                    recoAs::beauty,       recoAs::charm,    recoAs::strange,  recoAs::other};

  constexpr int size_basic  = static_cast<int>( recoAs::other ) + 1;
  constexpr int size_recoAs = static_cast<int>( recoAs::fifth ) + 1;
  constexpr int size_multi  = size_recoAs - size_basic;
  constexpr int begin_multi = static_cast<int>( recoAs::first );

  void collectProductss( LHCb::MCVertex const& mcpv, LHCb::MCVertex const& mcvtx,
                         std::vector<LHCb::MCParticle const*>& allprods ) {
    for ( auto const& idau : mcvtx.products() ) {
      double dv2 = ( mcpv.position() - idau->originVertex()->position() ).Mag2();
      if ( dv2 > ( 100. * Gaudi::Units::mm ) * ( 100. * Gaudi::Units::mm ) ) continue;
      allprods.push_back( idau.target() );
      for ( auto const& ivtx : idau->endVertices() ) collectProductss( mcpv, *ivtx, allprods );
    }
  }

  void count_reconstructible_mc_particles( LHCb::span<MCPVInfo> mcpvvec, LHCb::MCProperty const& flags,
                                           bool requireVelo ) {
    const MCTrackInfo trInfo = {flags};
    for ( auto& infomc : mcpvvec ) {
      LHCb::MCVertex const*                avtx         = infomc.pMCPV;
      unsigned                             mcPartInMCPV = 0;
      std::vector<LHCb::MCParticle const*> allproducts;
      collectProductss( *avtx, *avtx, allproducts );

      for ( auto const pmcp : allproducts ) {
        if ( pmcp->particleID().isMeson() || pmcp->particleID().isBaryon() ) {
          if ( pmcp->particleID().hasBottom() ) infomc.decayBeauty = true;
          if ( pmcp->particleID().hasCharm() ) infomc.decayCharm = true;
          if ( pmcp->particleID().hasStrange() ) infomc.decayStrange = true;
          if ( !requireVelo || trInfo.hasVelo( pmcp ) ) {
            double dv2 = ( avtx->position() - pmcp->originVertex()->position() ).Mag2();
            if ( dv2 < 0.0000001 && pmcp->p() > 100. * Gaudi::Units::MeV ) ++mcPartInMCPV;
          }
        }
      }
      infomc.nRecTracks = mcPartInMCPV;
    }
  }

  // int count_velo_tracks( const std::vector<LHCb::Event::v2::WeightedTrack>& tracksPV ) {
  //  return std::count_if( tracksPV.begin(), tracksPV.end(), []( const auto& t ) { return t.track->hasVelo(); } );
  //}
} // namespace

/// Helper to format recoAs enum class as the underlying number
template <>
struct fmt::formatter<recoAs> : fmt::formatter<std::underlying_type_t<recoAs>> {
  template <typename FormatContext>
  auto format( recoAs v, FormatContext& ctx ) {
    using und_t = std::underlying_type_t<recoAs>;
    return formatter<und_t>::format( static_cast<und_t>( v ), ctx );
  }
};

class PrimaryVertexChecker
    : public LHCb::Algorithm::Consumer<void( LHCb::Tracks const&, Vertices const&, LHCb::MCVertices const&,
                                             LHCb::MCParticles const&, LHCb::MCHeader const&, LHCb::MCProperty const&,
                                             DeVP const& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiTupleAlg, DeVP>> {
public:
  /// Standard constructor
  PrimaryVertexChecker( std::string const& name, ISvcLocator* pSvcLocator )
      : Consumer{name,
                 pSvcLocator,
                 {KeyValue{"inputTracksName", LHCb::TrackLocation::Default},
                  KeyValue{"inputVerticesName", LHCb::Event::PV::DefaultLocation},
                  KeyValue{"MCVertexInput", LHCb::MCVertexLocation::Default},
                  KeyValue{"MCParticleInput", LHCb::MCParticleLocation::Default},
                  KeyValue{"MCHeaderLocation", LHCb::MCHeaderLocation::Default},
                  KeyValue{"MCPropertyInput", LHCb::MCPropertyLocation::TrackInfo},
                  KeyValue{"DeVP", LHCb::Det::VP::det_path}}} {}

  StatusCode initialize() override; ///< Algorithm initialization
  void       operator()( LHCb::Tracks const& tracks, Vertices const& vertices, LHCb::MCVertices const& mcvtx,
                   LHCb::MCParticles const& mcps, LHCb::MCHeader const& mcheader, LHCb::MCProperty const& flags,
                   DeVP const& ) const override; ///< Algorithm execution
  StatusCode finalize() override;                      ///< Algorithm finalization

  using VertexType = Vertices::value_type;

private:
  Gaudi::Property<int>    m_nTracksToBeRecble{this, "nTracksToBeRecble", 4};  // min number of tracks in PV
  Gaudi::Property<bool>   m_produceHistogram{this, "produceHistogram", true}; // producing histograms (light version)
  Gaudi::Property<bool>   m_produceNtuple{this, "produceNtuple", false};      // producing NTuples (full version)
  Gaudi::Property<bool>   m_requireVelo{this, "RequireVelo", true};           // requiring VELO for tracks
  Gaudi::Property<double> m_dzIsolated{this, "dzIsolated", 10.0 * Gaudi::Units::mm}; // split close/isolated PVs
  Gaudi::Property<int>    m_nTracksToPrint{this, "nTracksToPrint", 10};              // split low/high multiplicity PVs
  Gaudi::Property<double> m_zToPrint{this, "zToPrint", 50.0 * Gaudi::Units::mm};     // split in z

  mutable Gaudi::Accumulators::Counter<>                            m_nevt{this, "nEvents"};
  mutable std::map<recoAs, Gaudi::Accumulators::BinomialCounter<>>  m_false;     // False PVs vs Reco PVs
  mutable std::map<recoAs, Gaudi::Accumulators::BinomialCounter<>>  m_eff;       // MC reconstructible vs Reconstructed
  mutable std::map<recoAs, Gaudi::Accumulators::BinomialCounter<>>  m_mcpv;      // MC vs MC reconstrucible
  mutable std::map<recoAs, Gaudi::Accumulators::AveragingCounter<>> m_av_mcp;    // average mc particles in MCPV
  mutable std::map<recoAs, Gaudi::Accumulators::AveragingCounter<>> m_av_tracks; // average tracks in RecoPV

  std::string toString( recoAs n ) const;
  bool        checkCondition( MCPVInfo const& MCPV, recoAs n ) const;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT_WITH_ID( PrimaryVertexChecker, "PrimaryVertexChecker" )

StatusCode PrimaryVertexChecker::initialize() { return Consumer::initialize(); }

//=============================================================================
// Main execution
//=============================================================================
void PrimaryVertexChecker::operator()( LHCb::Tracks const& tracks, Vertices const& recoVtx,
                                       LHCb::MCVertices const& mcvtx, LHCb::MCParticles const& mcps,
                                       LHCb::MCHeader const& mcheader, LHCb::MCProperty const& flags,
                                       DeVP const& vpdet ) const {
  ++m_nevt;

  const auto beamspot = vpdet.beamSpot();

  if ( msgLevel( MSG::DEBUG ) )
    debug() << " # tracks: " << tracks.size() << "  # vertices: " << recoVtx.size() << endmsg;

  // Fill reconstucted PV info
  std::vector<RecPVInfo<VertexType>> recpvvec;
  recpvvec.reserve( recoVtx.size() );

  for ( auto const& pv : recoVtx ) {
    RecPVInfo<VertexType> recinfo;
    recinfo.pRECPV        = &pv;
    recinfo.position      = pv.position();
    auto const& covPV     = pv.covMatrix();
    recinfo.positionSigma = Gaudi::XYZPoint{sqrt( covPV( 0, 0 ) ), sqrt( covPV( 1, 1 ) ), sqrt( covPV( 2, 2 ) )};
    recinfo.nTracks       = pv.nTracks();
    recinfo.chi2          = pv.chi2();
    recinfo.nDoF          = pv.nDoF();
    recinfo.indexMCPVInfo = -1;
    recinfo.mother        = 0;
    recinfo.indexMCPVInfo = -1;
    recpvvec.push_back( recinfo );
  }

  // Fill MC PV info
  std::vector<MCPVInfo> mcpvvec;
  mcpvvec.reserve( mcvtx.size() );

  for ( auto const& mcpv : mcvtx ) {
    LHCb::MCParticle const* motherPart = mcpv->mother();
    if ( !motherPart && mcpv->type() == LHCb::MCVertex::MCVertexType::ppCollision ) {
      MCPVInfo mcprimvert;
      mcprimvert.pMCPV = mcpv;
      mcpvvec.emplace_back( mcprimvert );
    }
  }

  count_reconstructible_mc_particles( mcpvvec, flags, m_requireVelo );

  std::sort( mcpvvec.begin(), mcpvvec.end(),
             []( MCPVInfo const& lhs, MCPVInfo const& rhs ) { return lhs.nRecTracks > rhs.nRecTracks; } );

  std::vector<MCPVInfo> rblemcpv;
  std::vector<MCPVInfo> not_rble_but_visible;
  std::vector<MCPVInfo> not_rble;

  for ( auto const& itmc : mcpvvec ) {
    if ( itmc.nRecTracks >= m_nTracksToBeRecble ) {
      rblemcpv.push_back( itmc );
    } else if ( itmc.nRecTracks > 1 ) {
      not_rble_but_visible.push_back( itmc );
    } else {
      not_rble.push_back( itmc );
    }
  }

  for ( int ipv = 0; ipv < static_cast<int>( recpvvec.size() ); ipv++ ) {
    match_mc_vertex_by_distance( ipv, LHCb::span<RecPVInfo<VertexType>>{recpvvec}, rblemcpv );
  }

  rblemcpv.insert( rblemcpv.end(), not_rble_but_visible.begin(), not_rble_but_visible.end() );
  rblemcpv.insert( rblemcpv.end(), not_rble.begin(), not_rble.end() );

  debug() << endmsg << " MC vertices " << endmsg;
  debug() << " ===================================" << endmsg;
  for ( auto const& [ipv, rmcpv] : LHCb::range::enumerate( rblemcpv ) ) {
    LHCb::MCVertex const* mcpv = rmcpv.pMCPV;
    debug() << fmt::format( " {:3} {:3}  xyz ( {:7.4f} {:7.4f} {:8.3f} )   nrec = {:4}", ipv, rmcpv.indexRecPVInfo,
                            mcpv->position().x(), mcpv->position().y(), mcpv->position().z(), rmcpv.nRecTracks )
            << ( rmcpv.indexRecPVInfo < 0 ? "  NOTRECO" : "" ) << endmsg;
  }

  debug() << " -----------------------------------" << endmsg << endmsg;

  debug() << endmsg << " REC vertices " << endmsg;
  debug() << " ===================================" << endmsg;
  for ( auto const& [ipv, recpv] : LHCb::range::enumerate( recpvvec ) ) {
    debug() << fmt::format(
                   " {:3} {:3}  xyz ( {:7.4f} {:7.4f} {:8.3f} )  ntra = {:4}  sigxyz ( {:7.4f} {:7.4f} {:8.4f} )   "
                   "chi2/NDF = {:7.2f}",
                   ipv, recpv.indexMCPVInfo, recpv.position.x(), recpv.position.y(), recpv.position.z(), recpv.nTracks,
                   recpv.positionSigma.x(), recpv.positionSigma.y(), recpv.positionSigma.z(),
                   recpv.pRECPV->chi2PerDoF() )
            << ( recpvvec[ipv].indexMCPVInfo < 0 ? "  FALSE" : "" ) << endmsg;
  }
  debug() << " -----------------------------------" << endmsg;

  // find nr of false PV
  int nFalsePV_real = 0;
  for ( auto const& [ipv, recpv] : LHCb::range::enumerate( recpvvec ) ) {
    int    fake   = 0;
    double x      = recpv.position.x();
    double y      = recpv.position.y();
    double z      = recpv.position.z();
    double r      = std::sqrt( x * x + y * y );
    double errx   = recpv.positionSigma.x();
    double erry   = recpv.positionSigma.y();
    double errz   = recpv.positionSigma.z();
    double errr   = std::sqrt( ( ( x * errx ) * ( x * errx ) + ( y * erry ) * ( y * erry ) ) / ( x * x + y * y ) );
    int    mother = recpv.mother;
    double chi2   = recpv.chi2;
    double nDoF   = recpv.nDoF;

    if ( recpv.indexMCPVInfo < 0 ) {
      fake = 1;

      if ( auto cmc = closestMCPV( rblemcpv, recpv ); cmc ) {
        for ( auto const& n : Part ) m_false[n] += checkCondition( *cmc, n );

        double dist = ( cmc->pMCPV->position() - recpv.pRECPV->position() ).R();
        m_false[dist > m_dzIsolated.value() ? recoAs::isolated : recoAs::close] += true;
        m_false[recpv.nTracks >= m_nTracksToPrint.value() ? recoAs::ntracks_high : recoAs::ntracks_low] += true;
        m_false[recpv.position.z() < -m_zToPrint.value()
                    ? recoAs::z_low
                    : ( recpv.position.z() < m_zToPrint.value() ? recoAs::z_middle : recoAs::z_high )] += true;

        int idx = cmc - rblemcpv.data();
        if ( idx < size_multi ) { m_false[All[begin_multi + idx]] += true; }
      }

      bool vis_found = false;
      for ( unsigned int imc = 0; imc < not_rble_but_visible.size(); imc++ ) {
        if ( not_rble_but_visible[imc].indexRecPVInfo > -1 ) continue;
        double dist = fabs( mcpvvec[imc].pMCPV->position().z() - recpv.position.z() );
        if ( dist < 5.0 * recpv.positionSigma.z() ) {
          vis_found                                = true;
          not_rble_but_visible[imc].indexRecPVInfo = 10;
          break;
        }
      } // imc

      if ( !vis_found ) nFalsePV_real++;
    }
    if ( m_produceNtuple ) {
      Tuple myTuple2 = nTuple( 102, "PV_nTuple2", CLID_ColumnWiseTuple );
      myTuple2->column( "fake", double( fake ) ).ignore();
      myTuple2->column( "r", double( r ) ).ignore();
      myTuple2->column( "x", double( x ) ).ignore();
      myTuple2->column( "y", double( y ) ).ignore();
      myTuple2->column( "z", double( z ) ).ignore();
      myTuple2->column( "errr", double( errr ) ).ignore();
      myTuple2->column( "errz", double( errz ) ).ignore();
      myTuple2->column( "errx", double( errx ) ).ignore();
      myTuple2->column( "erry", double( erry ) ).ignore();
      myTuple2->column( "mother", double( mother ) ).ignore();
      myTuple2->column( "chi2", double( chi2 ) ).ignore();
      myTuple2->column( "nDoF", double( nDoF ) ).ignore();
      myTuple2->write().ignore();
    }
  }
  // Fill distance to closest recble MC PV and its multiplicity
  for ( auto& mcpv : rblemcpv ) {
    auto cmc               = closestMCPV( rblemcpv, mcpv );
    mcpv.distToClosestMCPV = ( cmc ? ( cmc->pMCPV->position() - mcpv.pMCPV->position() ).R() : 999999. );
    mcpv.multClosestMCPV   = ( cmc ? cmc->nRecTracks : 0 );
  }

  for ( auto const& itmc : rblemcpv ) {
    for ( auto const& n : Basic ) {
      bool cut = checkCondition( itmc, n );
      if ( cut ) {
        if ( itmc.nRecTracks < m_nTracksToBeRecble ) {
          m_mcpv[n] += false;
        } else {
          m_mcpv[n] += true;
          m_av_mcp[n] += itmc.nRecTracks;
          if ( itmc.indexRecPVInfo < 0 ) {
            m_eff[n] += false;
          } else {
            m_eff[n] += true;
            m_false[n] += false;
          }
        }
        if ( itmc.indexRecPVInfo > -1 ) { m_av_tracks[n] += recpvvec[itmc.indexRecPVInfo].nTracks; }
      }
    }
  }

  int mcpv = 0;
  int high = 0;

  for ( auto& itmc : rblemcpv ) {
    double x             = -99999.;
    double y             = -99999.;
    double dx            = -99999.;
    double dy            = -99999.;
    double dz            = -99999.;
    double r             = -99999.;
    double z             = -99999.;
    double errx          = -99999.;
    double erry          = -99999.;
    double errz          = -99999.;
    double errr          = -99999.;
    double chi2          = -999999.;
    double nDoF          = -999999.;
    int    reconstructed = 0;
    int    ntracks_pvrec = 0;
    int    ntracks_pvmc  = 0;
    int    dtrcks        = 0;
    int    pvevt         = 0;
    int    mother        = 0;
    int    assoctrks     = 0;
    int    nassoctrks    = 0;

    double zMC = itmc.pMCPV->position().z();
    double yMC = itmc.pMCPV->position().y();
    double xMC = itmc.pMCPV->position().x();
    double rMC =
        std::sqrt( ( xMC - beamspot.x() ) * ( xMC - beamspot.x() ) + ( yMC - beamspot.y() ) * ( yMC - beamspot.y() ) );

    if ( mcpv < size_multi ) {
      if ( itmc.nRecTracks < m_nTracksToBeRecble ) {
        m_mcpv[All[begin_multi + mcpv]] += false;
      } else {
        m_mcpv[All[begin_multi + mcpv]] += true;
        m_av_mcp[All[begin_multi + mcpv]] += itmc.nRecTracks;
        if ( itmc.indexRecPVInfo < 0 ) {
          m_eff[All[begin_multi + mcpv]] += false;
        } else {
          m_eff[All[begin_multi + mcpv]] += true;
          m_false[All[begin_multi + mcpv]] += false;
        }
      }
    }

    if ( int indRec = itmc.indexRecPVInfo; indRec > -1 ) {
      auto const& recpv = recpvvec[indRec];
      ++high;
      ++pvevt;
      reconstructed = 1;
      dx            = recpv.position.x() - itmc.pMCPV->position().x();
      dy            = recpv.position.y() - itmc.pMCPV->position().y();
      dz            = recpv.position.z() - itmc.pMCPV->position().z();
      x             = recpv.position.x();
      y             = recpv.position.y();
      z             = recpv.position.z();
      r    = std::sqrt( ( x - beamspot.x() ) * ( x - beamspot.x() ) + ( y - beamspot.y() ) * ( y - beamspot.y() ) );
      errx = recpv.positionSigma.x();
      erry = recpv.positionSigma.y();
      errz = recpv.positionSigma.z();
      errr = std::sqrt( ( ( x * errx ) * ( x * errx ) + ( y * erry ) * ( y * erry ) ) / ( x * x + y * y ) );
      ntracks_pvrec = recpv.nTracks;
      ntracks_pvmc  = itmc.pMCPV->products().size();
      dtrcks        = ntracks_pvmc - ntracks_pvrec;
      mother        = recpv.mother;
      chi2          = recpv.chi2;
      nDoF          = recpv.nDoF;

      if ( mcpv < size_multi ) { m_av_tracks[All[begin_multi + mcpv]] += recpv.nTracks; }
      // Filling histograms
      if ( m_produceHistogram ) {
        plot( itmc.pMCPV->position().x(), 1001, "xmc", -0.25, 0.25, 50 );
        plot( itmc.pMCPV->position().y(), 1002, "ymc", -0.25, 0.25, 50 );
        plot( itmc.pMCPV->position().z(), 1003, "zmc", -20, 20, 50 );
        plot( recpv.position.x(), 1011, "xrd", -0.25, 0.25, 50 );
        plot( recpv.position.y(), 1012, "yrd", -0.25, 0.25, 50 );
        plot( recpv.position.z(), 1013, "zrd", -20, 20, 50 );
        plot( dx, 1021, "dx", -0.25, 0.25, 50 );
        plot( dy, 1022, "dy", -0.25, 0.25, 50 );
        plot( dz, 1023, "dz", -0.5, 0.5, 50 );
        plot( dx / errx, 1031, "pullx", -5., 5., 50 );
        plot( dy / erry, 1032, "pully", -5., 5., 50 );
        plot( dz / errz, 1033, "pullz", -5., 5., 50 );
        if ( itmc.nRecTracks < 10 ) {
          plot( dx, 1101, "dx", -0.25, 0.25, 50 );
          plot( dy, 1102, "dy", -0.25, 0.25, 50 );
          plot( dz, 1103, "dz", -0.5, 0.5, 50 );
        } else if ( itmc.nRecTracks >= 10 && itmc.nRecTracks < 30 ) {
          plot( dx, 1111, "dx", -0.25, 0.25, 50 );
          plot( dy, 1112, "dy", -0.25, 0.25, 50 );
          plot( dz, 1113, "dz", -0.5, 0.5, 50 );
        } else {
          plot( dx, 1121, "dx", -0.25, 0.25, 50 );
          plot( dy, 1122, "dy", -0.25, 0.25, 50 );
          plot( dz, 1123, "dz", -0.5, 0.5, 50 );
        }

        if ( itmc.pMCPV->position().z() < -m_zToPrint ) {
          plot( dx, 1201, "dx", -0.25, 0.25, 50 );
          plot( dy, 1202, "dy", -0.25, 0.25, 50 );
          plot( dz, 1203, "dz", -0.5, 0.5, 50 );
        } else if ( itmc.pMCPV->position().z() < m_zToPrint ) {
          plot( dx, 1211, "dx", -0.25, 0.25, 50 );
          plot( dy, 1212, "dy", -0.25, 0.25, 50 );
          plot( dz, 1213, "dz", -0.5, 0.5, 50 );
        } else {
          plot( dx, 1221, "dx", -0.25, 0.25, 50 );
          plot( dy, 1222, "dy", -0.25, 0.25, 50 );
          plot( dz, 1223, "dz", -0.5, 0.5, 50 );
        }

        plot( double( ntracks_pvrec ), 1041, "ntracks_pvrec", 0., 150., 50 );
        plot( double( dtrcks ), 1042, "mcrdtracks", 0., 150., 50 );
        if ( pvevt == 1 ) {
          plot( double( recpvvec.size() ), 1051, "nPVperEvt", -0.5, 5.5, 6 );
          for ( auto const& rp : recpvvec ) assoctrks += rp.nTracks;
          nassoctrks = tracks.size() - assoctrks;
          plot( double( nassoctrks ), 1052, "nassoctrks", 0., 150., 50 );
        }
      }
    }
    mcpv++;

    // Filling ntuple
    if ( m_produceNtuple ) {
      int    nmrc         = 0;
      double dist_closest = itmc.distToClosestMCPV;
      Tuple  myTuple      = nTuple( 101, "PV_nTuple", CLID_ColumnWiseTuple );
      myTuple->column( "reco", double( reconstructed ) ).ignore();
      myTuple->column( "isol", dist_closest > m_dzIsolated.value() ? 1. : 0. ).ignore();
      myTuple->column( "ntracks", double( ntracks_pvrec ) ).ignore();
      myTuple->column( "nrectrmc", double( itmc.nRecTracks ) ).ignore();
      myTuple->column( "dzclose", dist_closest ).ignore();
      myTuple->column( "nmcpv", double( rblemcpv.size() ) ).ignore();
      myTuple->column( "mtruemcpv", double( mcheader.numOfPrimaryVertices() ) ).ignore();
      myTuple->column( "nmcallpv", double( mcpvvec.size() ) ).ignore();
      myTuple->column( "nrecpv", double( recpvvec.size() ) ).ignore();
      myTuple->column( "decayCharm", double( itmc.decayCharm ) ).ignore();
      myTuple->column( "decayBeauty", double( itmc.decayBeauty ) ).ignore();
      myTuple->column( "decayStrange", double( itmc.decayStrange ) ).ignore();
      myTuple->column( "multirec", double( high ) ).ignore();
      myTuple->column( "multimc", double( mcpv ) ).ignore();
      myTuple->column( "dx", dx ).ignore();
      myTuple->column( "dy", dy ).ignore();
      myTuple->column( "dz", dz ).ignore();
      myTuple->column( "x", x ).ignore();
      myTuple->column( "y", y ).ignore();
      myTuple->column( "r", r ).ignore();
      myTuple->column( "zMC", zMC ).ignore();
      myTuple->column( "yMC", yMC ).ignore();
      myTuple->column( "xMC", xMC ).ignore();
      myTuple->column( "rMC", rMC ).ignore();
      myTuple->column( "z", z ).ignore();
      myTuple->column( "xBeam", beamspot.x() ).ignore();
      myTuple->column( "yBeam", beamspot.y() ).ignore();
      myTuple->column( "errx", errx ).ignore();
      myTuple->column( "erry", erry ).ignore();
      myTuple->column( "errz", errz ).ignore();
      myTuple->column( "errr", errr ).ignore();
      myTuple->column( "mother", double( mother ) ).ignore();
      myTuple->column( "evtnr", double( m_nevt.nEntries() ) ).ignore();
      myTuple->column( "chi2", double( chi2 ) ).ignore();
      myTuple->column( "nDoF", double( nDoF ) ).ignore();
      myTuple->column( "size_tracks", double( tracks.size() ) ).ignore();
      myTuple->column( "size_mcp", double( mcps.size() ) ).ignore();
      myTuple->column( "mcpvrec", double( nmrc ) ).ignore();
      myTuple->write().ignore();
    }
  }
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrimaryVertexChecker::finalize() {
  info() << "     ************************************ " << endmsg << " MC PV is reconstructible if at least "
         << m_nTracksToBeRecble.value() << "  tracks are reconstructed" << endmsg
         << " MC PV is isolated if dz to closest reconstructible MC PV >  " << m_dzIsolated.value() << " mm" << endmsg
         << " MC efficiency split by tracks with threshold: (" << m_nTracksToBeRecble.value() << ","
         << m_nTracksToPrint.value() << "), >= " << m_nTracksToPrint.value() << endmsg
         << " MC efficiency split by z position: <-" << m_zToPrint.value() << ", (-" << m_zToPrint.value() << ","
         << m_zToPrint.value() << "), >" << m_zToPrint.value() << endmsg << " REC and MC vertices matched:  by distance"
         << endmsg << endmsg;
  for ( auto const& n : All ) { info() << fmtRat( toString( n ), m_mcpv[n], m_eff[n], m_false[n] ) << endmsg; }
  info() << endmsg;
  for ( auto const& n : All ) { info() << fmtAvTracks( toString( n ), m_av_tracks[n], m_av_mcp[n] ) << endmsg; }
  info()
      << endmsg
      << fmtRes( "1_res_all", check_histogram( histo( HistoID( 1021 ) ), true ),
                 check_histogram( histo( HistoID( 1022 ) ), true ), check_histogram( histo( HistoID( 1023 ) ), true ) )
      << endmsg
      << fmtRes( "2_res_ntracks<10", check_histogram( histo( HistoID( 1101 ) ), true ),
                 check_histogram( histo( HistoID( 1102 ) ), true ), check_histogram( histo( HistoID( 1103 ) ), true ) )
      << endmsg
      << fmtRes( "3_res_ntracks(10,30)", check_histogram( histo( HistoID( 1111 ) ), true ),
                 check_histogram( histo( HistoID( 1112 ) ), true ), check_histogram( histo( HistoID( 1113 ) ), true ) )
      << endmsg
      << fmtRes( "4_res_ntracks>30", check_histogram( histo( HistoID( 1121 ) ), true ),
                 check_histogram( histo( HistoID( 1122 ) ), true ), check_histogram( histo( HistoID( 1123 ) ), true ) )
      << endmsg
      << fmtRes( "5_res_z<-50", check_histogram( histo( HistoID( 1201 ) ), true ),
                 check_histogram( histo( HistoID( 1202 ) ), true ), check_histogram( histo( HistoID( 1203 ) ), true ) )
      << endmsg
      << fmtRes( "6_res_z(-50,50)", check_histogram( histo( HistoID( 1211 ) ), true ),
                 check_histogram( histo( HistoID( 1212 ) ), true ), check_histogram( histo( HistoID( 1213 ) ), true ) )
      << endmsg
      << fmtRes( "7_res_z>50", check_histogram( histo( HistoID( 1221 ) ), true ),
                 check_histogram( histo( HistoID( 1222 ) ), true ), check_histogram( histo( HistoID( 1223 ) ), true ) )
      << endmsg << endmsg
      << fmtRes( "1_pull_width_all", check_histogram( histo( HistoID( 1031 ) ), true ),
                 check_histogram( histo( HistoID( 1032 ) ), true ), check_histogram( histo( HistoID( 1033 ) ), true ) )
      << endmsg
      << fmtRes( "1_pull_mean_all", check_histogram( histo( HistoID( 1031 ) ), false ),
                 check_histogram( histo( HistoID( 1032 ) ), false ),
                 check_histogram( histo( HistoID( 1033 ) ), false ) )
      << endmsg << endmsg;

  return Consumer::finalize(); // Must be called after all other actions
}

bool PrimaryVertexChecker::checkCondition( MCPVInfo const& MCPV, recoAs n ) const {
  switch ( n ) {
  case recoAs::all:
    return true;
  case recoAs::isolated:
    return ( MCPV.distToClosestMCPV > m_dzIsolated );
  case recoAs::close:
    return ( MCPV.distToClosestMCPV <= m_dzIsolated );
  case recoAs::ntracks_low:
    return ( MCPV.nRecTracks >= m_nTracksToBeRecble && MCPV.nRecTracks < m_nTracksToPrint );
  case recoAs::ntracks_high:
    return ( MCPV.nRecTracks >= m_nTracksToPrint );
  case recoAs::z_low:
    return ( MCPV.pMCPV->position().z() < -m_zToPrint );
  case recoAs::z_middle:
    return ( MCPV.pMCPV->position().z() >= -m_zToPrint && MCPV.pMCPV->position().z() < m_zToPrint );
  case recoAs::z_high:
    return ( MCPV.pMCPV->position().z() >= m_zToPrint );
  case recoAs::beauty:
    return ( MCPV.decayBeauty );
  case recoAs::charm:
    return ( MCPV.decayCharm );
  case recoAs::strange:
    return ( MCPV.decayStrange );
  case recoAs::other:
    return ( !( MCPV.decayBeauty ) && !( MCPV.decayCharm ) && !( MCPV.decayStrange ) );
  default:
    return false;
  }
}

std::string PrimaryVertexChecker::toString( recoAs n ) const {
  switch ( n ) {
  case recoAs::all:
    return fmt::format( "0{} all", recoAs::all );
  case recoAs::isolated:
    return fmt::format( "0{} isolated", recoAs::isolated );
  case recoAs::close:
    return fmt::format( "0{} close", recoAs::close );
  case recoAs::ntracks_low:
    return fmt::format( "0{} ntracks<{}", recoAs::ntracks_low, m_nTracksToPrint.value() );
  case recoAs::ntracks_high:
    return fmt::format( "0{} ntracks>={}", recoAs::ntracks_high, m_nTracksToPrint.value() );
  case recoAs::z_low:
    return fmt::format( "0{} z<{:2.1f}", recoAs::z_low, -m_zToPrint.value() );
  case recoAs::z_middle:
    return fmt::format( "0{} z in ({:2.1f}, {:2.1f})", recoAs::z_middle, -m_zToPrint.value(), +m_zToPrint.value() );
  case recoAs::z_high:
    return fmt::format( "0{} z >={:2.1f}", recoAs::z_high, +m_zToPrint.value() );
  case recoAs::beauty:
    return fmt::format( "0{} decayBeauty", recoAs::beauty );
  case recoAs::charm:
    return fmt::format( "0{} decayCharm", recoAs::charm );
  case recoAs::strange:
    return fmt::format( "{} decayStrange", recoAs::strange );
  case recoAs::other:
    return fmt::format( "{} other", recoAs::other );
  case recoAs::first:
    return fmt::format( "{} 1MCPV", recoAs::first );
  case recoAs::second:
    return fmt::format( "{} 2MCPV", recoAs::second );
  case recoAs::third:
    return fmt::format( "{} 3MCPV", recoAs::third );
  case recoAs::fourth:
    return fmt::format( "{} 4MCPV", recoAs::fourth );
  case recoAs::fifth:
    return fmt::format( "{} 5MCPV", recoAs::fifth );
  default:
    return "not defined";
  }
}
