/*****************************************************************************\
# (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      #
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "DetDesc/Condition.h"
#include "Event/RecVertex_v2.h"
#include "Event/Track.h"
#include "LHCbAlgs/Transformer.h"
#include "PrKernel/PrSelection.h"
#include "Timer.h"
#include "VPDet/DeVP.h"

#include "GaudiKernel/ToolHandle.h"

#include <boost/container/static_vector.hpp>

#include <array>
#include <vector>

//#define TIMINGHISTOGRAMMING 1

#ifdef TIMINGHISTOGRAMMING
#  include "AIDA/IProfile1D.h"
#  include "GaudiKernel/IHistogramSvc.h"
#endif

using namespace LHCb::Event::v2;
using TrackBeamLineVertexFinderOutput = std::tuple<std::vector<RecVertex>, Pr::Selection<Track>>;

/**
 * PV finding strategy:
 * step 1: select tracks with velo info and cache some information useful for PV finding
 * step 2: fill a histogram with the z of the poca to the beamline
 * step 3: do a peak search in that histogram ('vertex seeds')
 * step 4: assign tracks to the closest seed ('partitioning')
 * step 5: fit the vertices with an adapative vertex fit
 *
 *  @author Wouter Hulsbergen (Nikhef, 2018)
 **/
class TrackBeamLineVertexFinder : public LHCb::Algorithm::MultiTransformer<TrackBeamLineVertexFinderOutput(
                                                                               std::vector<Track> const&, DeVP const& ),
                                                                           LHCb::DetDesc::usesConditions<DeVP>> {
public:
  TrackBeamLineVertexFinder( const std::string& name, ISvcLocator* pSvcLocator );
  TrackBeamLineVertexFinderOutput operator()( std::vector<LHCb::Event::v2::Track> const&, DeVP const& ) const override;
  StatusCode                      initialize() override;

private:
  Gaudi::Property<uint32_t> m_minNumTracksPerVertex{this, "MinNumTracksPerVertex", 4};
  Gaudi::Property<float>    m_zmin{this, "MinZ", -300 * Gaudi::Units::mm, "Min z position of vertex seed"};
  Gaudi::Property<float>    m_zmax{this, "MaxZ", +300 * Gaudi::Units::mm, "Max z position of vertex seed"};
  Gaudi::Property<float>    m_dz{this, "ZBinSize", 0.25 * Gaudi::Units::mm, "Z histogram bin size"};
  Gaudi::Property<float>    m_maxTrackZ0Err{this, "MaxTrackZ0Err", 1.5 * Gaudi::Units::mm,
                                         "Maximum z0-error for adding track to histo"};
  Gaudi::Property<float>    m_minDensity{this, "MinDensity", 0. / Gaudi::Units::mm,
                                      "Minimal density at cluster peak  (inverse resolution)"};
  Gaudi::Property<float>    m_minDipDensity{this, "MinDipDensity", 3.0 / Gaudi::Units::mm,
                                         "Minimal depth of a dip to split cluster (inverse resolution)"};
  Gaudi::Property<float>    m_minTracksInSeed{this, "MinTrackIntegralInSeed", 2.5};
  Gaudi::Property<float>    m_maxVertexRho{this, "BeamSpotRCut", 0.3 * Gaudi::Units::mm,
                                        "Maximum distance of vertex to beam line"};
  Gaudi::Property<uint32_t> m_maxFitIter{this, "MaxFitIter", 5, "Maximum number of iterations for vertex fit"};
  Gaudi::Property<float> m_maxDeltaChi2{this, "MaxDeltaChi2", 12, "Maximum chi2 contribution of track to vertex fit"};
  Gaudi::Property<float> m_maxTrackBLChi2{this, "MaxTrackBLChi2", 10,
                                          "Maximum chi2 of track to beam line contributing to seeds"};
#ifdef TIMINGHISTOGRAMMING
  AIDA::IProfile1D* m_timeperstepPr{nullptr};
  AIDA::IProfile1D* m_timevsntrksPr{nullptr};
  AIDA::IProfile1D* m_timevsnvtxPr{nullptr};
#endif

  static const uint16_t Nztemplatebins        = 16;
  static const uint16_t Nztemplates           = 32;
  static const uint16_t TemplateNormalization = 128;
  uint16_t              m_ztemplates[2 * Nztemplates][Nztemplatebins]; // odd and even, see note
  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_nbPVsCounter{this, "Nb PVs"};
};

DECLARE_COMPONENT( TrackBeamLineVertexFinder )

TrackBeamLineVertexFinder::TrackBeamLineVertexFinder( const std::string& name, ISvcLocator* pSvcLocator )
    : MultiTransformer(
          name, pSvcLocator,
          {KeyValue{"InputTracks", LHCb::Event::v2::TrackLocation::Default}, KeyValue{"DEVP", LHCb::Det::VP::det_path}},
          {KeyValue{"OutputVertices", LHCb::Event::v2::RecVertexLocation::Primary},
           KeyValue{"UnusedTracks", "Rec/Track/NonPVTracks"}} ) {}

StatusCode TrackBeamLineVertexFinder::initialize() {
  return MultiTransformer::initialize().andThen( [&] {
#ifdef TIMINGHISTOGRAMMING
    auto hsvc       = service<IHistogramSvc>( "HistogramDataSvc", true );
    m_timeperstepPr = hsvc->bookProf( name() + "/timeperstep", "time per step", 20, -0.5, 19.5 );
    m_timevsntrksPr = hsvc->bookProf( name() + "/timevsntrks", "time vs number of tracks", 50, -0.5, 249.5 );
    m_timevsnvtxPr  = hsvc->bookProf( name() + "/timevsnvtx", "time vs number of vertices", 12, -0.5, 11.5 );
#endif

    // Fill the odd templates first
    const double sqrthalf = std::sqrt( 0.5 );
    {
      // of course, since the thing is symmetric we can do this more efficiently, but that's not quite worth it now
      const double zmaxodd = m_dz * ( Nztemplatebins / 2 - 1 + 0.5 );
      const double zminodd = -zmaxodd;
      for ( int itemplate = 0; itemplate < Nztemplates; ++itemplate ) {
        const double sigmaz   = m_maxTrackZ0Err * double( itemplate + 1 ) / Nztemplates;
        double       integral = 0.5 * std::erf( sqrthalf * zminodd / sigmaz );
        for ( int ibin = 0; ibin < Nztemplatebins - 1; ++ibin ) {
          double thisintegral = 0.5 * std::erf( sqrthalf * ( zminodd + ( ibin + 1 ) * m_dz ) / sigmaz );
          double bincontent   = thisintegral - integral;
          m_ztemplates[2 * itemplate + 1][ibin] = int( bincontent * TemplateNormalization );
          integral                              = thisintegral;
        }
        m_ztemplates[2 * itemplate + 1][Nztemplatebins - 1] = 0;
      }
    }

    // even templates
    {
      // of course, since the thing is symmetric we can do this more efficiently, but that's not quite worth it now
      const double zmaxeven = m_dz * Nztemplatebins / 2;
      const double zmineven = -zmaxeven;
      for ( int itemplate = 0; itemplate < Nztemplates; ++itemplate ) {
        const double sigmaz   = m_maxTrackZ0Err * double( itemplate + 1 ) / Nztemplates;
        double       integral = 0.5 * std::erf( sqrthalf * zmineven / sigmaz );
        for ( int ibin = 0; ibin < Nztemplatebins; ++ibin ) {
          double thisintegral               = 0.5 * std::erf( sqrthalf * ( zmineven + ( ibin + 1 ) * m_dz ) / sigmaz );
          double bincontent                 = thisintegral - integral;
          m_ztemplates[2 * itemplate][ibin] = int( bincontent * TemplateNormalization );
          integral                          = thisintegral;
        }
      }
    }
    return StatusCode::SUCCESS;
  } );
}

namespace {

  // structure with minimal track info needed for PV search
  struct PVTrack {
    PVTrack() {}
    PVTrack( const LHCb::State& state, double dz, uint16_t _index )
        : z{float( state.z() + dz )}
        , x{float( state.x() + dz * state.tx() ), float( state.y() + dz * state.ty() )}
        , tx{float( state.tx() ), float( state.ty() )}
        , index{_index} {
      // perhaps we should invert it /before/ switching to single FPP?
      // it doesn't seem to make much difference.
      const auto& V   = state.covariance();
      auto        dz2 = dz * dz;
      W( 0, 0 )       = V( 0, 0 ) + 2 * dz * V( 2, 0 ) + dz2 * V( 2, 2 );
      W( 1, 0 )       = V( 1, 0 ) + dz * ( V( 3, 0 ) + V( 2, 1 ) ) + dz2 * V( 3, 2 );
      W( 1, 1 )       = V( 1, 1 ) + 2 * dz * V( 3, 1 ) + dz2 * V( 3, 3 );
      W.Invert();
    }
    float                z{0};
    Gaudi::Vector2F      x;        /// position (x,y)
    Gaudi::Vector2F      tx;       /// direction (tx,ty)
    Gaudi::SymMatrix2x2F W;        /// weightmatrix
    uint16_t             index{0}; /// index in the list with tracks
  };

  template <typename FTYPE>
  auto sqr( FTYPE x ) {
    return x * x;
  }

  struct Extremum {
    Extremum( uint16_t _index, uint16_t _value, uint32_t _integral )
        : index{_index}, value{_value}, integral{_integral} {}
    uint16_t index;
    uint16_t value;
    uint32_t integral;
  };

  struct Cluster {
    Cluster( uint16_t _izfirst, uint16_t _izlast, uint16_t _izmax )
        : izfirst{_izfirst}, izlast{_izlast}, izmax{_izmax} {}
    uint16_t izfirst;
    uint16_t izlast;
    uint16_t izmax;
  };

  struct SeedZWithIteratorPair {
    using iterator = std::vector<PVTrack>::iterator;
    float    z;
    iterator begin;
    iterator end;
    SeedZWithIteratorPair( float _z, iterator _begin, iterator _end ) : z{_z}, begin{_begin}, end{_end} {}
  };

  // Need a small extension to the track when fitting the
  // vertex. Caching this information doesn't seem to help much
  // though.
  struct PVTrackInVertex : PVTrack {
    PVTrackInVertex( const PVTrack& trk ) : PVTrack{trk} {
      ROOT::Math::SMatrix<float, 3, 2> H;
      H( 0, 0 ) = H( 1, 1 ) = 1;
      H( 2, 0 )             = -trk.tx( 0 );
      H( 2, 1 )             = -trk.tx( 1 );
      HW                    = H * W;
      HWH                   = ROOT::Math::Similarity( H, W );
    }
    ROOT::Math::SMatrix<float, 3, 2> HW;
    Gaudi::SymMatrix3x3F             HWH;
    float                            weight{1};
  };

  struct Vertex {
    Gaudi::XYZPoint                         position;
    Gaudi::SymMatrix3x3                     poscov;
    std::vector<std::pair<unsigned, float>> tracks; // index to track + weight in vertex fit
    double                                  chi2;
  };

  // This implements the adapative vertex fit with Tukey's weights.
  Vertex fitAdaptive( const std::vector<PVTrack>::iterator& tracksbegin,
                      const std::vector<PVTrack>::iterator& tracksend, const Gaudi::XYZPoint& seedposition,
                      uint16_t maxNumIter, float chi2max ) {
    // make vector of TrackInVertex objects
    std::vector<PVTrackInVertex> tracks( tracksbegin, tracksend );
    bool                         converged = false;
    Vertex                       vertex;
    auto&                        vtxpos = vertex.position;
    Gaudi::SymMatrix3x3F         vtxcov;
    vtxpos = seedposition;
    const float maxDeltaZConverged{0.001};
    float       chi2tot{0};
    uint16_t    nselectedtracks{0};
    uint16_t    iter{0};
    for ( ; iter < maxNumIter && !converged; ++iter ) {
      Gaudi::SymMatrix3x3F halfD2Chi2DX2;
      Gaudi::Vector3F      halfDChi2DX;
      chi2tot         = 0;
      nselectedtracks = 0;
      Gaudi::Vector2F vtxposvec{float( vtxpos.x() ), float( vtxpos.y() )};
      for ( auto& trk : tracks ) {
        // compute the chi2
        const float           dz   = vtxpos.z() - trk.z;
        const Gaudi::Vector2F res  = vtxposvec - ( trk.x + dz * trk.tx );
        float                 chi2 = ROOT::Math::Similarity( res, trk.W );
        // compute the weight.
        trk.weight = 0.0f;
        if ( chi2 < chi2max ) { // to branch or not, that is the question!
          ++nselectedtracks;
          // Tukey's weight
          trk.weight = sqr( 1.f - chi2 / chi2max );
          // trk.weight = chi2 < 1 ? 1 : sqr( 1. - (chi2-1) / (chi2max-1) ) ;
          // += operator does not work for mixed FP types
          halfD2Chi2DX2 += trk.weight * trk.HWH;
          halfDChi2DX += trk.weight * trk.HW * res;
          // if I use expressions, it crashes!
          // const Gaudi::SymMatrix3x3F thisHalfD2Chi2DX2 = weight * ROOT::Math::Similarity(H, trk.W ) ;
          // const Gaudi::Vector3F HWr = trk.HW * res ;
          // for(int irow=0; irow<3; ++irow) {
          //   halfDChi2DX(irow) += trk.weight * HWr(irow) ;
          //   for(int icol=0; icol<=irow; ++icol)
          //     halfD2Chi2DX2(irow,icol) += trk.weight * trk.HWH(irow,icol) ;
          // }
          chi2tot += trk.weight * chi2;
        }
      }
      if ( nselectedtracks >= 2 ) {
        // compute the new vertex covariance
        vtxcov = halfD2Chi2DX2;
        /*int OK =*/vtxcov.InvertChol();

        // compute the delta w.r.t. the reference
        Gaudi::Vector3F delta = -1.0F * vtxcov * halfDChi2DX;

        // note: this is only correct if chi2 was chi2 of reference!
        chi2tot += ROOT::Math::Dot( delta, halfDChi2DX );

        // update the position
        vtxpos.SetX( vtxpos.x() + delta( 0 ) );
        vtxpos.SetY( vtxpos.y() + delta( 1 ) );
        vtxpos.SetZ( vtxpos.z() + delta( 2 ) );
        converged = std::abs( delta( 2 ) ) < maxDeltaZConverged;
      } else {
        break;
      }
    } // end iteration loop
    // std::cout << "Number of iterations: " << iter << " " << nselectedtracks << std::endl ;
    for ( int irow = 0; irow < 3; ++irow )
      for ( int icol = 0; icol <= irow; ++icol ) vertex.poscov( irow, icol ) = vtxcov( irow, icol );
    vertex.chi2 = chi2tot;
    vertex.tracks.reserve( tracks.size() );
    for ( const auto& trk : tracks ) {
      if ( trk.weight > 0 ) vertex.tracks.emplace_back( trk.index, trk.weight );
    }
    return vertex;
  }
} // namespace

TrackBeamLineVertexFinderOutput TrackBeamLineVertexFinder::
                                operator()( std::vector<LHCb::Event::v2::Track> const& tracks, DeVP const& deVP ) const {
  // Get the beamline. this only accounts for position, not
  // rotation. that's something to improve!
  Gaudi::XYZPoint beamline = deVP.beamSpot();

  // get the tracks
#ifdef TIMINGHISTOGRAMMING
  Timer timer[20];
  timer[9].start();
  timer[1].start();
#endif
  // Step 1: select tracks with velo info, compute the poca to the
  // beamline. cache the covariance matrix at this position. I'd
  // rather us a combination of copy_if and transform, but don't know
  // how to do that efficiently.
  const auto           Ntrk = tracks.size();
  std::vector<PVTrack> pvtracks( Ntrk ); // allocate everything upfront. don't use push_back/emplace_back
  {
    auto        it         = pvtracks.begin();
    const float halfwindow = ( Nztemplatebins / 2 + 1 ) * m_dz;
    const float zmin       = m_zmin + halfwindow;
    const float zmax       = m_zmax - halfwindow;
    for ( uint16_t index{0}; index < Ntrk; ++index ) {
      const auto& trk = tracks[index];
      if ( trk.hasVelo() ) {
        // compute the (chance in) z of the poca to the beam axis
        const LHCb::State& s  = trk.firstState();
        const auto         tx = s.tx();
        const auto         ty = s.ty();
        const double dz   = ( tx * ( beamline.x() - s.x() ) + ty * ( beamline.y() - s.y() ) ) / ( tx * tx + ty * ty );
        const double newz = s.z() + dz;
        if ( zmin < newz && newz < zmax ) {
          *it = PVTrack{s, dz, index};
          ++it;
        }
      }
    }
    pvtracks.erase( it, pvtracks.end() );
  }

#ifdef TIMINGHISTOGRAMMING
  timer[1].stop();
  timer[2].start();
#endif

  // Step 2: fill a histogram with the z position of the poca. Use the
  // projected vertex error on that position as the width of a
  // gauss. Divide the gauss properly over the bins. This is quite
  // slow: some simplification may help here.

  // we need to define what a bin is: integral between
  //   zmin + ibin*dz and zmin + (ibin+1)*dz
  // we'll have lot's of '0.5' in the code below. at some point we may
  // just want to shift the bins.

  // this can be changed into an std::accumulate
  const int             Nbins = ( m_zmax - m_zmin ) / m_dz;
  std::vector<uint16_t> zhisto( Nbins, 0 );
  for ( const auto& trk : pvtracks ) {
    // bin in which z0 is, in floating point
    const int jbin        = int( 4 * ( trk.z - m_zmin ) / m_dz ); // we need factor 4 to determine odd or even
    const int dbin        = jbin % 4;
    const int minbin      = jbin / 4 - Nztemplatebins / 2 + ( ( dbin == 0 ) ? 0 : 1 );
    const int oddtemplate = ( dbin == 0 || dbin == 3 ) ? 0 : 1;
    // make sure the template fits. make sure first and last bin
    // remain empty. there will be lot's of funny effects at edge of
    // histogram but we do not care.

    // about half the time is spent in logic here and computation of
    // zerr; and the other half in the additions.
    const float zweight = ROOT::Math::Similarity( trk.W, trk.tx );
    const float zerr    = 1 / std::sqrt( zweight ); // Q_rsqrt(zweight) ;
    // compute the bealine chi2. we make a cut for tracks contributing to seeding
    const float blchi2 = ROOT::Math::Similarity(
        trk.W, Gaudi::Vector2F{float( beamline.x() ) - trk.x( 0 ), float( beamline.y() ) - trk.x( 1 )} );
    const int zerrbin = Nztemplates * zerr / m_maxTrackZ0Err;
    if ( zerrbin < Nztemplates && blchi2 < m_maxTrackBLChi2 ) {
      // now perform a quick addition. this could easily be
      // parallelized. let's hope the compilor unroles the loop.
      for ( uint16_t i = 0; i < Nztemplatebins; ++i ) zhisto[i + minbin] += m_ztemplates[2 * zerrbin + oddtemplate][i];
    }
  }

#ifdef TIMINGHISTOGRAMMING
  timer[2].stop();
  timer[3].start();
#endif
  // Step 3: perform a peak search in the histogram. This used to be
  // very simple but the logic needed to find 'significant dips' made
  // it a bit more complicated. In the end it doesn't matter so much
  // because it takes relatively little time.

  std::vector<Cluster> clusters;
  {
    // step A: make 'ProtoClusters'
    // Step B: for each such ProtoClusters
    //    - find the significant extrema (an odd number, start with a minimum. you can always achieve this by adding a
    //    zero bin at the beginning)
    //      an extremum is a bin-index, plus the integral till that point, plus the content of the bin
    //    - find the highest extremum and
    //       - try and partition at the lowest minimum besides it
    //       - if that doesn't work, try the other extremum
    //       - if that doesn't work, accept as cluster

    // Step A: make 'proto-clusters': these are subsequent bins with non-zero content and an integral above the
    // threshold.
    const uint32_t minTracksInSeed = m_minTracksInSeed * TemplateNormalization;
    const uint32_t mindip          = m_minDipDensity * m_dz * TemplateNormalization; // need to invent something
    const uint32_t minpeak         = m_minDensity * m_dz * TemplateNormalization;

    using BinIndex = uint16_t;
    // FIXME: how dangerous is it to use a fixed capacity vector?
    boost::container::static_vector<BinIndex, 64> clusteredges;
    // std::vector<BinIndex> clusteredges ;
    {
      bool     prevempty = true;
      uint32_t integral  = zhisto[0];
      for ( BinIndex i = 1; i < Nbins; ++i ) {
        integral += zhisto[i];
        bool empty = zhisto[i] < 1;
        if ( empty != prevempty ) {
          if ( empty ) {
            if ( integral >= minTracksInSeed )
              clusteredges.emplace_back( i );
            else
              clusteredges.pop_back();
            integral = 0;
          } else {
            clusteredges.emplace_back( i - 1 );
          }
          prevempty = empty;
        }
      }
    }
    // Step B: turn these into clusters. There can be more than one cluster per proto-cluster.
    const size_t Nproto = clusteredges.size() / 2;
    for ( uint16_t i = 0; i < Nproto; ++i ) {
      const BinIndex ibegin = clusteredges[i * 2];
      const BinIndex iend   = clusteredges[i * 2 + 1];
      // find the extrema. all this complicated logic is to be able to
      // split close seeds.
      boost::container::static_vector<Extremum, 32> extrema;
      {
        bool     rising   = true;
        uint32_t integral = zhisto[ibegin];
        extrema.emplace_back( ibegin, zhisto[ibegin], integral );
        for ( uint16_t i = ibegin; i < iend; ++i ) {
          const auto value = zhisto[i];
          if ( value != zhisto[i + 1] ) {
            bool stillrising = zhisto[i + 1] > value;
            if ( stillrising != rising ) { // found a local extremum
              const uint16_t n          = extrema.size();
              const int16_t  deltavalue = value - extrema.back().value;
              if ( rising ) {       // found a local maximum
                if ( n % 2 == 1 ) { // the last one was a minimum. check that this maximum is high enough. we always
                                    // accept the first maximum.
                  if ( ( n == 1 && value >= minpeak ) || deltavalue >= int( mindip ) )
                    extrema.emplace_back( i, value, integral + value / 2 );
                } else { // the last one was a maximum, but apparently there were no good minima in between
                  if ( deltavalue > 0 ) {
                    extrema.pop_back();
                    extrema.emplace_back( i, value, integral + value / 2 );
                  }
                }
              } else {              // found a local minimum
                if ( n % 2 == 0 ) { // the last one was a maximum. check that this minimum is small enough
                  if ( -1 * deltavalue >= int( mindip ) ) extrema.emplace_back( i, value, integral + 0.5f * value );
                } else { // the last one was a minimum, but apparently there were no good maxima in between
                  if ( deltavalue < 0 ) {
                    extrema.pop_back();
                    extrema.emplace_back( i, value, integral + value / 2 );
                  }
                }
              }
            }
            rising = stillrising;
          }
          integral += value;
        }
        if ( extrema.size() % 2 == 1 ) { // last was minimum. this one should replace it
          extrema.pop_back();
        }
        extrema.emplace_back( iend, zhisto[iend], integral );
      }

      // FIXME: temporary logic check
      // if( extrema.size()%2==0 ) {
      // warning() << "ERROR: even number of extrema found." << extrema.size() << endmsg ;
      //}

      if ( extrema.size() >= 3 ) {
        // now partition on  extrema
        const auto                                   N = extrema.size();
        boost::container::static_vector<Cluster, 16> subclusters;
        // std::vector<Cluster> subclusters ;
        if ( N > 3 ) {
          for ( uint32_t i = 1; i < N / 2 + 1; ++i ) {
            if ( extrema[2 * i].integral - extrema[2 * i - 2].integral > minTracksInSeed ) {
              subclusters.emplace_back( extrema[2 * i - 2].index, extrema[2 * i].index, extrema[2 * i - 1].index );
            }
          }
        }

        if ( subclusters.empty() ) {
          clusters.emplace_back( extrema.front().index, extrema.back().index, extrema[1].index );
        } else {
          // adjust the limit of the first and last to extend to the entire protocluster
          subclusters.front().izfirst = ibegin;
          subclusters.back().izlast   = iend;
          clusters.insert( std::end( clusters ), std::begin( subclusters ), std::end( subclusters ) );
        }
      }
    }
  }

#ifdef TIMINGHISTOGRAMMING
  timer[3].stop();
  timer[4].start();
#endif

  // Step 4: partition the set of tracks by vertex seed: just
  // choose the closest one. The easiest is to loop over tracks and
  // assign to closest vertex by looping over all vertices. However,
  // that becomes very slow as time is proportional to both tracks and
  // vertices. A better method is to rely on the fact that vertices
  // are sorted in z, and then use std::partition, to partition the
  // track list on the midpoint between two vertices. The logic is
  // slightly complicated to deal with partitions that have too few
  // tracks. I checked it by comparing to the 'slow' method.

  // I found that this funny weighted 'maximum' is better than most other inexpensive solutions.
  auto zClusterMean = [this, &zhisto]( auto izmax ) -> float {
    const uint16_t* b   = zhisto.data() + izmax;
    int             d1  = *b - *( b - 1 );
    int             d2  = *b - *( b + 1 );
    float           idz = d1 + d2 > 0 ? ( 0.5f * ( d1 - d2 ) ) / ( d1 + d2 ) : 0.0f;
    return m_zmin + m_dz * ( izmax + idz + 0.5f );
  };

  std::vector<SeedZWithIteratorPair> seedsZWithIteratorPair;
  seedsZWithIteratorPair.reserve( clusters.size() );

  if ( !clusters.empty() ) {
    std::vector<PVTrack>::iterator it    = pvtracks.begin();
    int                            iprev = 0;
    for ( int i = 0; i < int( clusters.size() ) - 1; ++i ) {
      // const float zmid = 0.5*(zseeds[i+1].z+zseeds[i].z) ;
      const float                    zmid = m_zmin + m_dz * 0.5f * ( clusters[i].izlast + clusters[i + 1].izfirst + 1 );
      std::vector<PVTrack>::iterator newit =
          std::partition( it, pvtracks.end(), [zmid]( const auto& trk ) { return trk.z < zmid; } );
      // complicated logic to get rid of partitions that are too small, doign the least amount of work
      if ( std::distance( it, newit ) >= m_minNumTracksPerVertex ) {
        seedsZWithIteratorPair.emplace_back( zClusterMean( clusters[i].izmax ), it, newit );
        iprev = i;
      } else {
        // if the partition is too small, then repartition the stuff we
        // have just isolated and assign to the previous and next. You
        // could also 'skip' this partition, but then you do too much
        // work for the next.
        if ( !seedsZWithIteratorPair.empty() && newit != it ) {
          const float zmid = m_zmin + m_dz * ( clusters[iprev].izlast + clusters[i + 1].izfirst + 0.5f );
          newit            = std::partition( it, newit, [zmid]( const auto& trk ) { return trk.z < zmid; } );
          // update the last one
          seedsZWithIteratorPair.back().end = newit;
        }
      }
      it = newit;
    }
    // Make sure to add the last partition
    if ( std::distance( it, pvtracks.end() ) >= m_minNumTracksPerVertex ) {
      seedsZWithIteratorPair.emplace_back( zClusterMean( clusters.back().izmax ), it, pvtracks.end() );
    } else if ( !seedsZWithIteratorPair.empty() ) {
      seedsZWithIteratorPair.back().end = pvtracks.end();
    }
  }

#ifdef TIMINGHISTOGRAMMING
  timer[4].stop();
  timer[5].start();
#endif

  // Step 5: perform the adaptive vertex fit for each seed.
  std::vector<Vertex> vertices;
  std::transform( seedsZWithIteratorPair.begin(), seedsZWithIteratorPair.end(), std::back_inserter( vertices ),
                  [&]( const auto& seed ) {
                    return fitAdaptive( seed.begin, seed.end, Gaudi::XYZPoint{beamline.x(), beamline.y(), seed.z},
                                        m_maxFitIter, m_maxDeltaChi2 );
                  } );

#ifdef TIMINGHISTOGRAMMING
  timer[5].stop();
  timer[6].start();
#endif

  // Steps that we could still take:
  // * remove vertices with too little tracks
  // * assign unused tracks to other vertices
  // * merge vertices that are close

  // create the output containers
  TrackBeamLineVertexFinderOutput output{std::vector<RecVertex>(), {tracks, Pr::details::alwaysFalse{}}};
  auto&                           recvertexcontainer = std::get<0>( output );
  //  std::vector<LHCb::Event::v2::RecVertex> recvertexcontainer;
  recvertexcontainer.reserve( vertices.size() );
  std::vector<bool> unusedtrackflags( Ntrk, true );
  const auto        maxVertexRho2 = sqr( m_maxVertexRho.value() );
  uint32_t          used_tracks{0};
  for ( const auto& vertex : vertices ) {
    const auto beamlinedx   = vertex.position.x() - beamline.x();
    const auto beamlinedy   = vertex.position.y() - beamline.y();
    const auto beamlinerho2 = sqr( beamlinedx ) + sqr( beamlinedy );
    if ( vertex.tracks.size() >= m_minNumTracksPerVertex && beamlinerho2 < maxVertexRho2 ) {
      int   nDoF      = 2 * vertex.tracks.size() - 3;
      auto& recvertex = recvertexcontainer.emplace_back( vertex.position, vertex.poscov,
                                                         LHCb::Event::v2::Track::Chi2PerDoF{vertex.chi2 / nDoF, nDoF} );
      recvertex.setTechnique( LHCb::Event::v2::RecVertex::RecVertexType::Primary );
      for ( const auto& dau : vertex.tracks ) {
        recvertex.addToTracks( &( tracks[dau.first] ), dau.second );
        unusedtrackflags[dau.first] = false;
        ++used_tracks;
      }
    }
  }
  m_nbPVsCounter += recvertexcontainer.size();
#ifdef TIMINGHISTOGRAMMING
  timer[6].stop();
  timer[7].start();
#endif
  auto& unusedtrackscontainer = std::get<1>( output );
  // the indices in Pr::Selection need to be sorted.
  // std::sort is expensive, so we use counting sort instead.
  // unusedtrackscontainer.m_container = tracks ;
  // unusedtrackscontainer.m_indices = std::move( unusedtracks ) ;
  // std::sort( begin(unusedtrackscontainer.m_indices), end(unusedtrackscontainer.m_indices) ) ;

  unusedtrackscontainer.m_indices.reserve( Ntrk - used_tracks );
  for ( uint32_t i = 0; i < Ntrk; ++i )
    if ( unusedtrackflags[i] ) unusedtrackscontainer.m_indices.push_back( i );

#ifdef TIMINGHISTOGRAMMING
  timer[7].stop();
  timer[9].stop();
  for ( int i = 0; i < 20; ++i ) m_timeperstepPr->fill( float( i ), timer[i].total() );
  m_timevsntrksPr->fill( pvtracks.size(), timer[9].total() );
  m_timevsnvtxPr->fill( vertices.size(), timer[9].total() );
#endif
  return output;
}
