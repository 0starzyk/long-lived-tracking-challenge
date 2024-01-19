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
// Include files
#include "Event/MCHeader.h"
#include "Event/MCVertex.h"
#include "Event/RecVertex_v2.h"
//#include "Event/RecVertex.h"
#include "Event/TrackVertexUtils.h"
#include "Event/Track_v2.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "LHCbAlgs/Consumer.h"
#include "LHCbMath/BloomFilter.h"
#include <limits>

/** @class TrackPVChecker
 *
 *
 *  @author Wouter Hulsbergen
 *  @date   2018
 */

using namespace LHCb;
using namespace LHCb::Event::v2;
using PrimaryVertexContainer = LHCb::Event::v2::RecVertices;
using Track                  = LHCb::Event::v2::Track;
// using PrimaryVertexContainer = LHCb::RecVertices ;
using TrackContainer = std::vector<Track>;
// using TrackContainer =  KeyedContainer<LHCb::Event::v1::Track,Containers::KeyedObjectManager<Containers::hashmap> > ;

class TrackPVChecker : public LHCb::Algorithm::Consumer<void( PrimaryVertexContainer const&, TrackContainer const&,
                                                              const LHCb::MCParticles& /*, const LHCb::LinksByKey&*/ ),
                                                        Gaudi::Functional::Traits::useGaudiHistoAlg> {
public:
  /// Standard constructor
  TrackPVChecker( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"VertexLocation", "" /* LHCb::RecVertexLocation::Primary*/},
                   KeyValue{"TrackLocation", "" /* LHCb::TrackLocation::Default*/},
                   KeyValue{"MCParticleLocation", LHCb::MCParticleLocation::Default}
                   /*,KeyValue{"LinkerInTable", "Link/" + LHCb::TrackLocation::Default}*/} ) {}

  void operator()( PrimaryVertexContainer const&, ///< Algorithm execution
                   TrackContainer const&, const LHCb::MCParticles& /*, const LHCb::LinksByKey&*/ ) const override;

private:
  Gaudi::Property<double> m_maxIPChi2ForTrackMatch{this, "MaxIPChi2ForTrackMatch", 10};
  Gaudi::Property<double> m_deltaZForMatch{this, "DeltaZForMatch", 1 * Gaudi::Units::mm};
  Gaudi::Property<int>    m_minNumParticles{this, "MinNumParticles", 3};
  Gaudi::Property<double> m_maxRadius{this, "MaxRadius", 3 * Gaudi::Units::mm};
};

DECLARE_COMPONENT( TrackPVChecker )

namespace {
  template <class T>
  std::vector<Track const*> const _extract_tracks( const T& tracks ) {
    std::vector<Track const*> ts( tracks.size(), nullptr );
    std::transform( tracks.begin(), tracks.end(), ts.begin(), []( auto const& wt ) { return wt.track; } );
    return ts;
  }

  template <class T>
  std::vector<float> const _extract_weights( const T& tracks ) {
    std::vector<float> ts( tracks.size() );
    std::transform( tracks.begin(), tracks.end(), ts.begin(), []( auto const& wt ) { return wt.weight; } );
    return ts;
  }

  template <class TrackContainer, class MCParticleContainer>
  std::vector<std::pair<const Track*, const LHCb::MCParticle*>>
  trackToMCParticle( const TrackContainer& tracks, const MCParticleContainer& mcparticles ) {
    // first collect all MC particles that we could potentially match to tracks
    std::vector<const LHCb::MCParticle*> mctracks;
    constexpr std::array                 pids{11, 13, 211, 321, 2212};
    for ( const auto& particle : mcparticles ) {
      if ( std::find( pids.begin(), pids.end(), std::abs( particle->particleID().pid() ) ) != pids.end() &&
           particle->originVertex() && std::abs( particle->originVertex()->position().z() ) < 700 * Gaudi::Units::mm ) {
        float L2{0};
        for ( const auto& decayvtx : particle->endVertices() ) {
          auto thisL2 = ( decayvtx->position() - particle->originVertex()->position() ).Mag2();
          if ( L2 < thisL2 ) L2 = thisL2;
        }
        if ( std::sqrt( L2 ) > 100 * Gaudi::Units::mm ) mctracks.push_back( particle );
      }
    }
    // now match every track to the best mcparticle, assuming that there is one
    std::vector<std::pair<const Track*, const LHCb::MCParticle*>> matches;
    for ( const auto& trk : tracks ) {
      const LHCb::State& s            = trk.firstState();
      auto               weightmatrix = s.covariance().Sub<Gaudi::SymMatrix4x4>( 0, 0 );
      weightmatrix.Invert();
      double                  bestchi2 = 100;
      const LHCb::MCParticle* bestmcp{0};
      for ( const auto& mcp : mctracks ) {
        const float mctx = mcp->momentum().Px() / mcp->momentum().Pz();
        const float mcty = mcp->momentum().Py() / mcp->momentum().Pz();
        float       dtx  = s.tx() - mctx;
        float       dty  = s.ty() - mcty;
        if ( std::abs( dtx ) < 0.01 && std::abs( dty ) < 0.01 ) {
          const auto mcpos = mcp->originVertex()->position();
          float      dz    = s.z() - mcpos.z();
          float      dx    = s.x() - ( mcpos.x() + mctx * dz );
          float      dy    = s.y() - ( mcpos.y() + mcty * dz );
          if ( std::abs( dx ) < 10 * Gaudi::Units::mm && std::abs( dy ) < 10 * Gaudi::Units::mm ) {
            // compute the chi2 for the match
            Gaudi::Vector4 delta{dx, dy, dtx, dty};
            double         chi2 = ROOT::Math::Similarity( delta, weightmatrix );
            if ( chi2 < bestchi2 ) {
              bestmcp  = mcp;
              bestchi2 = chi2;
            }
          }
        }
      }
      // std::cout << "Best match: " << bestchi2 << std::endl ;
      matches.emplace_back( &trk, bestmcp );
    }
    return matches;
  }

  const LHCb::MCVertex* mcoriginvertex( const LHCb::MCParticle& mcp ) {
    const LHCb::MCVertex*   mcv    = mcp.originVertex();
    const LHCb::MCParticle* mother = mcv->mother();
    if ( mcv ) {
      // check if the particle has a mother
      if ( mother && mother->originVertex() ) {
        if ( std::abs( mcv->position().z() - mother->originVertex()->position().z() ) < 0.002 * Gaudi::Units::mm )
          mcv = mcoriginvertex( *mother );
      }
    } else if ( mother )
      mcv = mcoriginvertex( *mother );
    return mcv;
  }

  enum ParticleOrigin { Unmatched = 0, Primary, Beauty, Charm, Strange, Conversion, Unknown };
  ParticleOrigin vertexOrigin( const LHCb::MCVertex& mcvertex ) {
    ParticleOrigin origin = ParticleOrigin::Unknown;
    if ( mcvertex.isPrimary() )
      origin = ParticleOrigin::Primary;
    else {
      // get the mother
      auto mother = mcvertex.mother();
      if ( mother ) {
        if ( mother->particleID().hasBottom() )
          origin = ParticleOrigin::Beauty;
        else if ( mother->particleID().hasCharm() )
          origin = ParticleOrigin::Charm;
        else if ( mother->particleID().hasStrange() )
          origin = ParticleOrigin::Strange;
        else if ( mother->particleID().pid() == 22 )
          origin = ParticleOrigin::Conversion;
        else
          origin = ParticleOrigin::Unknown;
      }
    }
    return origin;
  }

  // template<typename HitIDContainer>
  // uint32_t uniqueVeloSegmentID( const HitIDContainer& lhcbids) {
  //   const uint32_t* begin   = reinterpret_cast<const uint32_t*>( &( *( lhcbids.begin() ) ) );
  //   const uint32_t* end     = reinterpret_cast<const uint32_t*>( &( *( lhcbids.end() ) ) );
  //   BloomFilterImpl::HashFNV1a<uint32_t, uint32_t> hashobj;
  //   // const uint32_t hash1 = hashobj(lhcbid);
  //   // hashing a second object on top of it
  //   // uint32_t hash2 = hashobj(lhcbid, hash1);

  //   // hashing a vector of objects - slightly ugly, since we need to get the
  //   // initial value of an "empty" hash from deep within the implementation...
  //   uint32_t        hash   = BloomFilterImpl::__doFNV1a<uint32_t>::hashinit;
  //   const uint32_t* it     = begin;
  //   const uint32_t  maxval = LHCbID{LHCb::LHCbID::channelIDtype::UT, 0}.lhcbID(); // PK-R3C: was TT
  //   while ( it < end && *it < maxval ) {
  //     hash = hashobj( *it, hash );
  //     ++it;
  //   }
  //   return hash;
  // }

  // struct TrackWithVeloId
  // {
  //   const Track* track{null_ptr} ;
  //   uint32_t uniqueID{0} ;
  //   const Track& operator*() const { return *track ; }
  //   TrackWithVeloId( const LHCb::Track& atrack ) :
  //     track{&atrack},uniqueID{uniqueVeloSegmentID(atrack.lhcbIDs())}{}
  // } ;

  struct MCVertexWithTracks {
    MCVertexWithTracks( const LHCb::MCVertex* vtx, const std::vector<const Track*>& tracks )
        : mcvertex{vtx}, recotracks{tracks} {
      std::sort( std::begin( recotracks ), std::end( recotracks ) );
    }
    bool                          isPrimary() const { return mcvertex->isPrimary(); }
    const LHCb::MCVertex*         mcvertex{0};
    std::vector<const Track*>     recotracks;
    std::vector<const RecVertex*> matches;
  };

  struct RecVertexInfo {
    const RecVertex*                       recvertex;
    std::vector<const Track*>              sortedtracks;
    std::vector<const MCVertexWithTracks*> matches;
    RecVertexInfo( const RecVertex* v ) : recvertex{v}, sortedtracks{_extract_tracks( v->tracks() )} {
      std::sort( sortedtracks.begin(), sortedtracks.end() );
    }
  };

  double sqr( double x ) { return x * x; }

  // adapted from std::set_intersection
  template <class InputIterator1, class InputIterator2>
  size_t count_intersection( InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
                             InputIterator2 last2 ) {
    size_t rc( 0 );
    while ( first1 != last1 && first2 != last2 ) {
      if ( *first1 < *first2 ) {
        ++first1;
      } else if ( *first2 < *first1 ) {
        ++first2;
      } else {
        ++first1;
        ++first2;
        ++rc;
      }
    }
    return rc;
  }
} // namespace

void TrackPVChecker::operator()( PrimaryVertexContainer const& vertices, ///< Algorithm execution
                                 TrackContainer const&         tracks,
                                 const LHCb::MCParticles&      mcparticles /*, const LHCb::LinksByKey&*/ ) const {
  // convert the keyed container to a container of pointers, or even of objects, for simplicity.
  const int MaxNumPVs = 15;
  // get the beamline
  const auto beamline = Gaudi::XYZVector{0, 0, 0};

  // see if we can match the tracks
  auto tracktomcparticlemap = trackToMCParticle( tracks, mcparticles );
  // now match every track to an MC origin-vertex, walking up as far as we can on the decay tree
  std::map<const LHCb::MCVertex*, std::vector<const Track*>> mcVertexToTracks;
  for ( const auto& m : tracktomcparticlemap ) {
    const auto&           trk      = m.first;
    const LHCb::MCVertex* mcvertex = m.second ? mcoriginvertex( *m.second ) : 0;
    // fill some info on the track, provided it actually has an mc match
    if ( mcvertex ) {
      // fill the map from vertex to tracks
      mcVertexToTracks[mcvertex].push_back( m.first );

      // for the rest we only care about tracks that actually come from a primary
      if ( mcvertex->isPrimary() ) {
        // extrapolate
        LHCb::State state = trk->firstState();
        state.linearTransportTo( mcvertex->position().z() );
        // double doca = LHCb::TrackVertexUtils::doca( state, mcvertex->position() ) ;
        const double ipchi2 = LHCb::TrackVertexUtils::vertexChi2( state, mcvertex->position(), Gaudi::SymMatrix3x3{} );
        plot1D( ipchi2, "trkipchi2tomcvertex", 0.0, 25.0 );
        Gaudi::XYZVector dxtrue = state.position() - mcvertex->position();
        plot( dxtrue.x(), "trkdxtrue", -0.4, 0.4 );
        plot( dxtrue.y(), "trkdytrue", -0.4, 0.4 );
        plot( dxtrue.x() / std::sqrt( state.covariance()( 0, 0 ) ), "trkdxtruepull", "track x pull", -5, 5 );
        plot( dxtrue.y() / std::sqrt( state.covariance()( 1, 1 ) ), "trkdytruepull", "track y pull", -5, 5 );
        plot( trk->chi2PerDoF(), "trkchi2dof", "track chi2/dof", 0, 5 );
        // let's see if we get better pulls if we scale the errors
        const double covscale = trk->chi2PerDoF() < 1 ? 1 : trk->chi2PerDoF();
        plot( dxtrue.x() / std::sqrt( covscale * state.covariance()( 0, 0 ) ), "trkdxpullscaled", "track x pull scaled",
              -5, 5 );
        plot( dxtrue.y() / std::sqrt( covscale * state.covariance()( 1, 1 ) ), "trkdypullscaled", "track y pull scaled",
              -5, 5 );

        // extrapolate the track to the beam line
        const auto   tx = state.tx();
        const auto   ty = state.ty();
        const double dz =
            ( tx * ( beamline.x() - state.x() ) + ty * ( beamline.y() - state.y() ) ) / ( tx * tx + ty * ty );
        state.linearTransportTo( state.z() + dz );
        // compute the error on z0
        Gaudi::SymMatrix2x2 W = state.covariance().Sub<Gaudi::SymMatrix2x2>( 0, 0 );
        W.Invert();
        Gaudi::Vector2 txy( tx, ty );
        double         z0weight = ROOT::Math::Similarity( Gaudi::Vector2{tx, ty}, W );
        double         z0err    = std::sqrt( 1 / z0weight );
        plot( z0err, "trkz0err", 0, 2 );
        plot( z0weight, "trkz0weight", 0, 400 );
        profile1D( mcvertex->position().z(), log( z0err ) / M_LN10, "trklogz0errvsz", "log(z0 error) versus z", -150,
                   150 );
        double blchi2 = ROOT::Math::Similarity( Gaudi::Vector2{beamline.x() - state.x(), beamline.y() - state.y()}, W );
        plot( blchi2, "trkblchi2", 0, 50 );

        // look at the distance to that vertex
        Gaudi::XYZVector dx0 = state.position() - mcvertex->position();
        plot( dx0.z(), "trkdz0", -5.0, 5.0 );
        plot( dx0.y(), "trkdx0", -0.4, 0.4 );
        plot( dx0.x(), "trkdy0", -0.4, 0.4 );
        const double t2 = sqr( state.tx() ) + sqr( state.ty() );
        const double t  = std::sqrt( t2 );
        profile1D( t, z0err, "trkz0errvst", 0, 0.4 );
        profile1D( t, std::abs( dx0.z() ), "trkabsdzvst", 0, 0.4 );
        profile1D( t, std::abs( dx0.y() ), "trkabsdyvst", 0, 0.4 );
        profile1D( t, std::abs( dx0.x() ), "trkabsdxvst", 0, 0.4 );
      }
    }
  }

  // get the full list of MC PVs
  // const auto mcvertices = get<LHCb::MCVertices>(LHCb::MCVertexLocation::Default );
  const LHCb::MCHeader* mcheader = get<LHCb::MCHeader>( LHCb::MCHeaderLocation::Default );
  // plot the number of tracks per primary vertex
  const auto numTruePVs = mcheader->primaryVertices().size();
  plot( numTruePVs, "nTruePVs", -0.5, 15.5, 16 );
  for ( const auto& mcvertex : mcheader->primaryVertices() )
    plot( mcVertexToTracks[mcvertex].size(), "trackspermcpv", "tracks in mc pv", -0.5, 100.5, 101 );

  // Choose the subset with at least two reconstructed tracks
  std::vector<MCVertexWithTracks> mcvertices;
  for ( const auto& it : mcVertexToTracks )
    if ( it.second.size() > 1 ) mcvertices.emplace_back( it.first, it.second );

  // order these vertices in z
  std::sort( mcvertices.begin(), mcvertices.end(), []( const auto& a, const auto& b ) -> bool {
    return a.mcvertex->position().z() < b.mcvertex->position().z();
  } );

  // let's check how many tracks we assign with the new method
  const auto numReconstructableTruePVs =
      std::count_if( mcvertices.begin(), mcvertices.end(), []( const auto& v ) -> bool { return v.isPrimary(); } );
  plot( numReconstructableTruePVs, "nReconstructableTruePVs", -0.5, 15.5, 16 );

  for ( const auto& mcv : mcvertices ) {
    const auto N = mcv.recotracks.size();
    if ( mcv.isPrimary() )
      plot( N, "trackspermcprimary", "tracks in mc primary", -0.5, 100.5, 101 );
    else
      plot( N, "trackspermcsecondary", "tracks in mc secondary", -0.5, 10.5, 11 );
    plot( vertexOrigin( *( mcv.mcvertex ) ), "mcvertexorigin", 0.5, 6.5, 6 );

    // For the new seeding method, we would like to know the maximum opening angle between the tracks.
    if ( mcv.isPrimary() ) {
      double maxdt2{0};
      for ( size_t itrk = 0; itrk < N; ++itrk )
        for ( size_t jtrk = 0; jtrk < itrk; ++jtrk ) {
          const auto dtx = mcv.recotracks[itrk]->firstState().tx() - mcv.recotracks[jtrk]->firstState().tx();
          const auto dty = mcv.recotracks[itrk]->firstState().ty() - mcv.recotracks[jtrk]->firstState().ty();
          const auto dt2 = dtx * dtx + dty * dty;
          if ( maxdt2 < dt2 ) maxdt2 = dt2;
        }
      plot1D( std::sqrt( maxdt2 ), "maxopeningangle", "max opening angle", 0, 1.0, 100 );
    }
  }

  // get the reconstructed vertices and match them to the most overlapping MC vertex
  std::vector<const RecVertex*> unmatchedvertices;
  for ( const auto& pv : vertices ) {
    RecVertexInfo recvertexinfo{&pv};
    const auto&   tracks = recvertexinfo.sortedtracks;

    MCVertexWithTracks* bestmcvertex{nullptr};
    double              bestabsdeltaz{9999};
    size_t              bestoverlap{0};
    size_t              totaloverlap{0};

    // std::cout << "Reco PV tracks: " << tracks.size() << std::endl ;
    // for(const auto& tr : tracks) std::cout << tr << std::endl ;
    for ( auto& mcpv : mcvertices ) {
      // I don't understand how this even compiles ...
      size_t overlap = count_intersection( std::begin( tracks ), std::end( tracks ), std::begin( mcpv.recotracks ),
                                           std::end( mcpv.recotracks ) );
      // std::cout << "Overlap: " << mcpv.recotracks.size() << " " << overlap << std::endl ;
      totaloverlap += overlap;

      // Now that the vertex->track matching is broken for the default
      // alg, we'll also need an alternative matching.

      // We will also match if it is close enough. We need this now
      // because in some versions we have the wrong tracks matched to
      // the vertex.
      const double absdz    = std::abs( pv.position().z() - mcpv.mcvertex->position().z() );
      const double dzpull   = absdz / std::sqrt( pv.covMatrix()( 2, 2 ) );
      bool         goodpull = dzpull < 5.0;
      if ( goodpull ) {
        // valid match: let's just match it: this is what we'll use for efficiency, fakerate etc
        mcpv.matches.emplace_back( &pv );
        recvertexinfo.matches.emplace_back( &mcpv );
      }
      if ( overlap >= 2 || goodpull ) {
        if ( !bestmcvertex || ( overlap > bestoverlap ) || ( overlap == bestoverlap && ( absdz ) < bestabsdeltaz ) ) {
          bestmcvertex  = &mcpv;
          bestoverlap   = overlap;
          bestabsdeltaz = absdz;
        }
      }
    }
    plot1D( recvertexinfo.matches.size(), "nummc2recmatches", -0.5, 4.5, 5 );
    const double totalweight{std::accumulate( pv.tracks().begin(), pv.tracks().end(), double( 0.0 ),
                                              []( auto&& init, auto const& rhs ) { return init + rhs.weight; } )};
    // info() << "Matched? : " << bestmcvertex << endmsg ;

    // bool matched = false;
    if ( bestmcvertex ) {
      // this is what we'll use for resolution
      const float purity = tracks.size() > 0 ? float( bestoverlap ) / tracks.size() : 0;
      if ( bestmcvertex->isPrimary() )
        plot( purity, "matchpurityprimary", "purity for primary vertex", 0.0, 1.02, 51 );
      else
        plot( purity, "matchpuritysecondary", "purity for secondary vertex", 0.0, 1.02, 51 );
      // if(  (float(bestoverlap)/totaloverlap > minPurityForMatch) || (bestabsdeltaz < 1*Gaudi::Units::mm) ) {
      // bestmcvertex->matches.emplace_back(pv) ;
      // fwewstd::cout << "Adding match:" << pv.position().z() << " " << bestmcvertex->mcvertex->position().z() <<
      // std::endl ;
      // matched = true;
      /*
  } else {
  std::cout << "Cannot match: " << float(bestoverlap)/totaloverlap << " " << bestabsdeltaz << std::endl ;
  }
      */
      // make reso/pull plots
      if ( bestmcvertex->isPrimary() ) {
        auto delta = pv.position() - bestmcvertex->mcvertex->position();
        plot1D( delta.x(), "dx", -0.1, 0.1 );
        plot1D( delta.y(), "dy", -0.1, 0.1 );
        plot1D( delta.z(), "dz", -0.5, 0.5 );
        const double errx{std::sqrt( pv.covMatrix()( 0, 0 ) )};
        const double erry{std::sqrt( pv.covMatrix()( 1, 1 ) )};
        const double errz{std::sqrt( pv.covMatrix()( 2, 2 ) )};
        plot1D( errx, "errx", 0, 0.1 );
        plot1D( erry, "erry", 0, 0.1 );
        plot1D( errz, "errz", 0, 0.5 );
        plot1D( delta.x() / errx, "dxpull", -5, 5 );
        plot1D( delta.y() / erry, "dypull", -5, 5 );
        plot1D( delta.z() / errz, "dzpull", -5, 5 );
        if ( delta.z() / errz < 5 ) {
          plot1D( totalweight, "totalweightmatchedpvs", "total track weight for matched pvs", 0, 50 );
          plot1D( tracks.size(), "numtracksmatchedpvs", "number of tracks for matched pvs", 0, 50 );
        }
      } else {
        auto delta = pv.position() - bestmcvertex->mcvertex->position();
        plot1D( delta.x(), "secdx", -0.1, 0.1 );
        plot1D( delta.y(), "secdy", -0.1, 0.1 );
        plot1D( delta.z(), "secdz", -0.5, 0.5 );
        plot1D( tracks.size(), "secntracks", 0.5, 100.5, 100 );
      }
    }

    profile1D( numTruePVs, recvertexinfo.matches.empty(), "fakeratevsntruepv", "fake rate vs num true PVs", 0.5,
               MaxNumPVs + 0.5, MaxNumPVs );
    profile1D( tracks.size(), recvertexinfo.matches.empty(), "fakeratevsntracks", "fake rate vs num tracks in vertex",
               1.5, 41.5, 40 );
    if ( recvertexinfo.matches.size() > 0 ) {
      const float mergerate = recvertexinfo.matches.size() - 1;
      profile1D( numTruePVs, mergerate, "mergeratevsntruepv", "fake rate vs num true PVs", 0.5, MaxNumPVs + 0.5,
                 MaxNumPVs );
    } else {
      plot1D( totalweight, "totalweightfakes", "total track weight for fakes", 0, 50 );
      plot1D( tracks.size(), "numtracksfakes", "number of tracks for fakes", 0, 50 );
      plot1D( std::sqrt( pv.covMatrix()( 2, 2 ) ), "errzfakes", "sigma(z) for fakes", 0.0, 1.0 );
    }

    if ( recvertexinfo.matches.empty() ) {
      unmatchedvertices.push_back( &pv );
      if ( tracks.size() > 5 ) {
        std::stringstream os;
        os << "How could a PV with more than 5 tracks not be matched to MC? " << pv.position().z() << " "
           << tracks.size() << " " << bestoverlap << " " << totaloverlap << std::endl;
        os << "MC vertices: " << std::endl;
        for ( auto& mcpv : mcvertices ) {
          os << mcpv.recotracks.size() << " "
             << count_intersection( std::begin( tracks ), std::end( tracks ), std::begin( mcpv.recotracks ),
                                    std::end( mcpv.recotracks ) )
             << " " << mcpv.isPrimary() << " " << mcpv.mcvertex->position().z() << std::endl;
        }
        os << "Rec vertices: " << std::endl;
        for ( const auto& jpv : vertices ) os << jpv.position().z() << " " << jpv.tracks().size() << std::endl;
        verbose() << os.str() << endmsg;
      }
    }
  }

  // plot chi2 difference between vertices
  {
    size_t N = vertices.size();
    for ( size_t i = 0; i < N; ++i )
      for ( size_t j = 0; j < i; ++j ) {
        const auto&         ivertex  = vertices[i];
        const auto&         jvertex  = vertices[j];
        auto                deltavec = ivertex.position() - jvertex.position();
        Gaudi::Vector3      delta{deltavec.x(), deltavec.y(), deltavec.z()};
        Gaudi::SymMatrix3x3 W = ivertex.covMatrix() + jvertex.covMatrix();
        W.InvertChol();
        double dchi2 = ROOT::Math::Similarity( delta, W );
        plot1D( dchi2, "deltachi2", "vertex-vertex delta chi2", 0, 50 );
      }
  }

  // make some plots
  std::vector<MCVertexWithTracks> mcpvs;
  std::copy_if( mcvertices.begin(), mcvertices.end(), std::back_inserter( mcpvs ),
                []( const auto& v ) { return v.isPrimary(); } );
  plot1D( mcpvs.size(), "numreconstructibleprimaries", "num reconstructible primary vertices", -0.5, MaxNumPVs + 0.5,
          MaxNumPVs + 1 );
  plot1D( mcvertices.size() - mcpvs.size(), "numreconstructiblesecondaries", "num reconstructible secondary vertices",
          -0.5, MaxNumPVs + 0.5, MaxNumPVs + 1 );
  plot2D( mcpvs.size(), vertices.size(), "numpvsvsnummcpvs", "num reco pvs vs num reconstructible pvs", -0.5,
          MaxNumPVs + 0.5, -0.5, MaxNumPVs + 0.5, MaxNumPVs + 1, MaxNumPVs + 1 );

  bool                      somethingstrange{false};
  const MCVertexWithTracks* prev{0};
  for ( const auto& mcpv : mcvertices ) {
    const double dx = mcpv.mcvertex->position().x() - beamline.x();
    const double dy = mcpv.mcvertex->position().y() - beamline.y();
    const double R  = std::sqrt( dx * dx + dy * dy );

    if ( mcpv.isPrimary() ) {
      plot1D( mcpv.mcvertex->position().x(), "mcpvx", -1.0, 1.0 );
      plot1D( mcpv.mcvertex->position().y(), "mcpvy", -1.0, 1.0 );
      plot1D( mcpv.mcvertex->position().z(), "mcpvz", -300.0, 300.0 );
      plot1D( dx, "mcpvdx", -0.4, 0.4 );
      plot1D( dy, "mcpvdy", -0.4, 0.4 );
      plot1D( mcpv.matches.size(), "numrec2mcmatches", -0.5, 4.5, 5 );
      profile1D( mcpv.recotracks.size(), !mcpv.matches.empty(), "effvsntracks",
                 "efficiency vs number of tracks in true PV", -0.5, 40.5, 41 );
      if ( !mcpv.matches.empty() ) {
        const size_t nclones = mcpv.matches.size() - 1;
        profile1D( mcpv.recotracks.size(), nclones, "cloneratevsntracks", "clone rate vs num reco tracks in true PV",
                   0.5, 100.5 );
        profile1D( numTruePVs, nclones, "cloneratevsnprimaries", "clone rate vs number of true PVs", 0.5,
                   MaxNumPVs + 0.5, MaxNumPVs );
      }

      if ( mcpv.recotracks.size() >= 5 ) {
        profile1D( R, !mcpv.matches.empty(), "effvsR", "efficiency vs Rxy", 0, 0.2, 20 );
        profile1D( mcpv.mcvertex->position().z(), !mcpv.matches.empty(), "effvsZ", "efficiency vs z", -150., 150., 30 );
        profile1D( numTruePVs, !mcpv.matches.empty(), "effvsnprimaries", "efficiency vs number of true PVs", 0.5,
                   MaxNumPVs + 0.5, MaxNumPVs );
        if ( prev && prev->recotracks.size() >= 5 ) {
          double deltaz = mcpv.mcvertex->position().z() - prev->mcvertex->position().z();
          profile1D( deltaz, ( !mcpv.matches.empty() || !prev->matches.empty() ), "effonevsdeltaz",
                     "efficiency to find one vs #Delta z", 0, 10, 40 );
          profile1D( deltaz, !mcpv.matches.empty() && !prev->matches.empty(), "effbothvsdeltaz",
                     "efficiency to find both vs #Delta z", 0, 10, 40 );
          profile1D( deltaz, !mcpv.matches.empty() && !prev->matches.empty() && mcpv.matches != prev->matches,
                     "effuniquevsdeltaz", "efficiency to find both vs #Delta z", 0, 10, 40 );

          plot1D( deltaz, "distancetonextmcvertex", 0.0, 50.0 );
          // see if there is an overlapping reco vertex
          int noverlap{0};
          for ( const auto& recpv1 : mcpv.matches )
            for ( const auto& recpv2 : prev->matches )
              if ( recpv1 == recpv2 ) ++noverlap;
          profile1D( deltaz, noverlap, "mergeratevsdeltaz", "merge rate vs #Delta z", 0, 2, 20 );
        }
      }
      if ( mcpv.matches.size() > 1 ) {
        auto pv1 = mcpv.matches[0];
        auto pv2 = mcpv.matches[1];
        if ( pv2->tracks().size() > pv1->tracks().size() ) std::swap( pv1, pv2 );
        plot1D( pv1->tracks().size(), "multiplematchntrk1", 0.5, 50.5 );
        plot1D( pv2->tracks().size(), "multiplematchntrk2", 0.5, 50.5 );
        plot1D( pv2->position().z() - pv1->position().z(), "multiplematchdz", -5.0, 5.0 );
      }

      // how can we miss a vertex with 8 tracks?
      if ( mcpv.recotracks.size() >= 8 && mcpv.matches.empty() ) somethingstrange = true;
      prev = &mcpv;
    } else {
      profile1D( R, !mcpv.matches.empty(), "effvsRsec", "efficiency vs Rxy for secondaries", 0, 1.0 );
      profile1D( mcpv.recotracks.size(), !mcpv.matches.empty(), "effvsntrackssec",
                 "efficiency vs num reco tracks for secondaries", -0.5, 10.5, 11 );
    }
  }

  if ( somethingstrange ) {
    std::stringstream os;
    os << "MC PVs: " << std::endl;
    for ( const auto& mcpv : mcvertices )
      os << mcpv.isPrimary() << std::setprecision( 3 ) << std::setw( 10 ) << mcpv.mcvertex->position().z()
         << std::setw( 10 ) << ( mcpv.mcvertex->position() - beamline ).Rho() << " " << mcpv.recotracks.size() << " "
         << mcpv.matches.size() << " " << std::endl;
    os << "REC PVs: " << std::endl;
    for ( const auto& recpv : vertices ) os << recpv.position().z() << " " << recpv.tracks().size() << std::endl;
    verbose() << os.str() << endmsg;
  }

  std::set<const Track*> usedtracks;
  plot1D( vertices.size(), "numpvs", "pv multiplicity", -0.5, MaxNumPVs + 0.5, MaxNumPVs + 1 );
  for ( const auto& recpv : vertices ) {
    plot1D( recpv.tracks().size(), "pvnumtracks", "Number of tracks", 0.5, 100.5, 100 );
    plot1D( recpv.chi2PerDoF(), "pvchi2dof", "Chi2 per dof", 0., 10., 50 );
    plot1D( recpv.position().x(), "pvx", -1.0, 1.0 );
    plot1D( recpv.position().y(), "pvy", -1.0, 1.0 );
    plot1D( recpv.position().z(), "pvz", -300.0, 300.0 );
    plot1D( std::sqrt( recpv.covMatrix()( 0, 0 ) ), "pvsigx", 0, 1. );
    plot1D( std::sqrt( recpv.covMatrix()( 1, 1 ) ), "pvsigy", 0, 1. );
    plot1D( std::sqrt( recpv.covMatrix()( 2, 2 ) ), "pvsigz", 0, 1. );

    size_t numbackward{0};
    // now the trouble is that v1 track and v2 recvertex have a
    // different track with weights vector. definitely like the
    // v2::recvertex solution better.
    for ( const auto& wtrk : recpv.tracks() ) {
      plot1D( wtrk.weight, "pvtrkweight", "Track weight", 0, 1 );
      if ( wtrk.track ) plot1D( int( wtrk.track->type() ), "pvtrktype", "Track type", -0.5, 20.5, 21 );
      usedtracks.insert( wtrk.track );
      if ( !wtrk.track || wtrk.track->isVeloBackward() ) ++numbackward;
    }
    plot1D( numbackward, "pvnumbackwardtracks", "Number of backward tracks", -0.5, 50.5, 51 );
  }

  for ( const auto* recpv : unmatchedvertices ) {
    plot1D( recpv->tracks().size(), "unmatchednumtracks", 0.5, 100.5, 50 );
    const double errz{std::sqrt( recpv->covMatrix()( 2, 2 ) )};
    plot1D( errz, "unmatchederrz", 0, 0.5 );
    for ( const auto& trk : _extract_weights( recpv->tracks() ) ) plot1D( trk, "unmatchedtrkweight", 0, 1 );
  }

  if ( tracks.size() > 0 ) {
    const float fractionused = usedtracks.size() / float( tracks.size() );
    plot( fractionused, "fractionusedtracks", "Fraction of used tracks", 0.0, 1.0 );
    profile1D( numTruePVs, fractionused, "fractionusedvsntruepv", "Fraction of used tracks vs num PVs", 0.5,
               MaxNumPVs + 0.5, MaxNumPVs );

    // also count the number of tracks with pv association and the number of used tracks with pv association
    size_t numpvtracks{0}, numusedpvtracks{0};
    for ( const auto& mcpv : mcvertices ) {
      if ( mcpv.isPrimary() ) {
        numpvtracks += mcpv.recotracks.size();
        for ( const auto& trk : mcpv.recotracks )
          if ( std::find( usedtracks.begin(), usedtracks.end(), trk ) != usedtracks.end() )
            ++numusedpvtracks;
          else {
            // plot the ipchi2 to the closest vertex
            double ipchi2{999.};
            for ( const auto& recpv : vertices ) {
              double thisipchi2 =
                  LHCb::TrackVertexUtils::vertexChi2( trk->firstState(), recpv.position(), recpv.covMatrix() );
              ipchi2 = std::min( ipchi2, thisipchi2 );
            }
            plot1D( ipchi2, "ipchi2formissedtracks", "IP chi2 for missed tracks", 0.0, 25.0 );
          }
      }
    }

    if ( !vertices.empty() ) {
      for ( const auto& trk : tracks ) {
        // find out if this track is from a secondary or a primary or unknown
        ParticleOrigin origin = ParticleOrigin::Unmatched;
        for ( const auto& vtotrks : mcVertexToTracks ) {
          auto it = std::find( vtotrks.second.begin(), vtotrks.second.end(), &trk );
          if ( it != vtotrks.second.end() ) origin = vertexOrigin( *vtotrks.first );
        }
        plot1D( int( origin ), "trackorigin", "track origin", -0.5, 6.5, 7 );

        // compute the ip chi2 to the closest reconstructed vertex
        double ipchi2{std::numeric_limits<double>::max()};
        for ( const auto& recpv : vertices ) {
          double thisipchi2 =
              LHCb::TrackVertexUtils::vertexChi2( trk.firstState(), recpv.position(), recpv.covMatrix() );
          ipchi2 = std::min( ipchi2, thisipchi2 );
        }
        bool used = std::find( usedtracks.begin(), usedtracks.end(), &trk ) != usedtracks.end();
        profile1D( int( origin ), float( used ), "fractionoftrksusedvstrkorigin",
                   "fraction of tracks used vs track origin", -0.5, 5.5, 6 );

        if ( used ) {
          plot1D( ipchi2, "ipchi2forusedtracks", "IP chi2 for used tracks", 0.0, 25.0 );
        } else {
          plot1D( ipchi2, "ipchi2forunusedtracks", "IP chi2 for unused tracks", 0.0, 25.0 );
        }
        const double covscale     = trk.chi2PerDoF() < 1 ? 1 : trk.chi2PerDoF();
        const double scaledipchi2 = ipchi2 / covscale;
        switch ( origin ) {
        case Primary:
          plot1D( ipchi2, "ipchi2fortruepvtracks", "IP chi2 for true PV tracks", 0.0, 25.0 );
          plot1D( scaledipchi2, "scaledipchi2fortruepvtracks", "scaled IP chi2 for true PV tracks", 0.0, 25.0 );
          break;
        case Beauty:
          plot1D( ipchi2, "ipchi2fortruebeautytracks", "IP chi2 for true beauty tracks", 0.0, 25.0 );
          plot1D( scaledipchi2, "scaledipchi2fortruebeautytracks", "scaled IP chi2 for true beauty tracks", 0.0, 25.0 );
          break;
        case Charm:
          plot1D( ipchi2, "ipchi2fortruecharmtracks", "IP chi2 for true charm tracks", 0.0, 25.0 );
          plot1D( scaledipchi2, "scaledipchi2fortruecharmtracks", "scaled IP chi2 for true charm tracks", 0.0, 25.0 );
          break;
        case Strange:
          plot1D( ipchi2, "ipchi2fortruestrangetracks", "IP chi2 for true strange tracks", 0.0, 25.0 );
          plot1D( scaledipchi2, "scaledipchi2fortruestrangetracks", "scaled IP chi2 for true strange tracks", 0.0,
                  25.0 );
          break;
        default:
          break;
        }
      }
    }
    if ( numpvtracks > 0 ) plot( numusedpvtracks / float( numpvtracks ), "fractionusedpvtracks", 0.0, 1.0 );
  }
}
