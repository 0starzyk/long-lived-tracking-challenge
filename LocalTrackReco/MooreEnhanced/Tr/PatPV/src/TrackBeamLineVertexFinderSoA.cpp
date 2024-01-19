/*****************************************************************************\
# (c) Copyright 2000-2021 CERN for the benefit of the LHCb Collaboration      #
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/PrVeloTracks.h"
#include "Event/RecVertex_v2.h"
#include "Event/Track.h"
#include "GaudiKernel/ToolHandle.h"
#include "Kernel/EventLocalAllocator.h"
#include "Kernel/STLExtensions.h"
#include "LHCbAlgs/Transformer.h"
#include "VPDet/DeVP.h"

#ifdef TIMINGHISTOGRAMMING
#  include "AIDA/IProfile1D.h"
#  include "GaudiKernel/IHistogramSvc.h"
#  include "Timer.h"
#endif

#include "Event/PrimaryVertices.h"

// boost includes
#include <boost/container/static_vector.hpp>

// std includes
#include <array>
#include <limits>
#include <vector>

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
namespace {

  // c++20's remove_cvref
  template <typename T>
  struct remove_cvref {
    using type = std::remove_cv_t<std::remove_reference_t<T>>;
  };

  template <typename T>
  using remove_cvref_t = typename remove_cvref<T>::type;

  using simd    = SIMDWrapper::best::types;
  using float_v = simd::float_v;
  using int_v   = simd::int_v;

  template <typename T>
  auto to_std_array( T&& some_v ) {
    if constexpr ( std::is_same_v<remove_cvref_t<T>, float_v> ) {
      std::array<float, simd::size> tmp;
      some_v.store( tmp.data() );
      return tmp;
    } else if ( std::is_same_v<remove_cvref_t<T>, int_v> ) {
      std::array<int, simd::size> tmp;
      some_v.store( tmp.data() );
      return tmp;
    }
  }

  constexpr float defaultMinZ     = -300 * Gaudi::Units::mm;
  constexpr float defaultMaxZ     = +300 * Gaudi::Units::mm;
  constexpr float defaultZBinSize = 0.25 * Gaudi::Units::mm;
  constexpr int   defaultNZBins   = int( ( defaultMaxZ - defaultMinZ ) / defaultZBinSize );
} // namespace

using namespace LHCb::Event::PV;

class TrackBeamLineVertexFinderSoA
    : public LHCb::Algorithm::Transformer<PrimaryVertexContainer( const EventContext&, const LHCb::Pr::Velo::Tracks&,
                                                                  const LHCb::Pr::Velo::Tracks&, const DeVP& ),
                                          LHCb::DetDesc::usesConditions<DeVP>> {
public:
  /// Standard constructor
  TrackBeamLineVertexFinderSoA( const std::string& name, ISvcLocator* pSvcLocator );
  /// Initialization
  StatusCode initialize() override;
  /// Execution
  PrimaryVertexContainer operator()( const EventContext&, const LHCb::Pr::Velo::Tracks&, const LHCb::Pr::Velo::Tracks&,
                                     const DeVP& ) const override;

private:
  Gaudi::Property<uint32_t> m_minNumTracksPerVertex{this, "MinNumTracksPerVertex", 4};
  Gaudi::Property<float>    m_zmin{this, "MinZ", defaultMinZ, "Min z position of vertex seed"};
  Gaudi::Property<float>    m_zmax{this, "MaxZ", defaultMaxZ, "Max z position of vertex seed"};
  Gaudi::Property<float>    m_dz{this, "ZBinSize", defaultZBinSize, "Z histogram bin size"};
  Gaudi::Property<float>    m_maxTrackZ0Err{this, "MaxTrackZ0Err", 1.5 * Gaudi::Units::mm,
                                         "Maximum z0-error for adding track to histo"};
  Gaudi::Property<float>    m_minDensity{this, "MinDensity", 0. / Gaudi::Units::mm,
                                      "Minimal density at cluster peak  (inverse resolution)"};
  Gaudi::Property<float>    m_minDipDensity{this, "MinDipDensity", 3.0 / Gaudi::Units::mm,
                                         "Minimal depth of a dip to split cluster (inverse resolution)"};
  Gaudi::Property<float>    m_minTracksInSeed{this, "MinTrackIntegralInSeed", 2.5};
  Gaudi::Property<float>    m_maxVertexRho{this, "BeamSpotRCut", 0.3 * Gaudi::Units::mm,
                                        "Maximum distance of vertex to beam line"};
  Gaudi::Property<uint32_t> m_maxFitIter{this, "MaxFitIter", 10, "Maximum number of iterations for vertex fit"};
  Gaudi::Property<float> m_maxDeltaChi2{this, "MaxDeltaChi2", 12, "Maximum chi2 contribution of track to vertex fit"};
  Gaudi::Property<float> m_maxDeltaZConverged{this, "MaxDeltaZConverged", 0.001,
                                              "Limit on change in z to determine if vertex fit has converged"};
  Gaudi::Property<float> m_maxDeltaChi2Converged{this, "MaxDeltaChi2Converged", 0.01,
                                                 "Limit on change in chi2 to determine if vertex fit has converged"};
  Gaudi::Property<float> m_minVertexZSeparationChi2{this, "MinVertexZSeparationChi2", 16};
  Gaudi::Property<float> m_minVertexZSeparation{this, "MinVertexZSeparation", 0.5f * Gaudi::Units::mm};
  Gaudi::Property<float> m_maxTrackBLChi2{this, "MaxTrackBLChi2", 10,
                                          "Maximum chi2 of track to beam line contributing to seeds"};
  Gaudi::Property<float> m_beamLineOffsetX{this, "BeamLineOffsetX", 0.0f, "X correction applied to beamline position"};
  Gaudi::Property<float> m_beamLineOffsetY{this, "BeamLineOffsetY", 0.0f, "Y correction applied to beamline position"};
  static constexpr uint16_t Nztemplatebins        = 16;
  static constexpr uint16_t Nztemplates           = 32;
  static constexpr uint16_t TemplateNormalization = 128;
  uint16_t                  m_ztemplates[2 * Nztemplates][Nztemplatebins]; // odd and even, see note
  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_nbPVsCounter{this, "Nb PVs"};

#ifdef TIMINGHISTOGRAMMING
  AIDA::IProfile1D* m_timeperstepPr{nullptr};
  AIDA::IProfile1D* m_timevsntrksPr{nullptr};
  AIDA::IProfile1D* m_timevsnvtxPr{nullptr};
#endif
};

DECLARE_COMPONENT( TrackBeamLineVertexFinderSoA )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TrackBeamLineVertexFinderSoA::TrackBeamLineVertexFinderSoA( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator,
                   {KeyValue{"TracksBackwardLocation", "Rec/Track/VeloBackward"},
                    KeyValue{"TracksLocation", "Rec/Track/Velo"}, KeyValue{"DeVP", LHCb::Det::VP::det_path}},
                   KeyValue{"OutputVertices", LHCb::Event::v2::RecVertexLocation::Primary} ) {}

//=============================================================================
// ::initialize()
//=============================================================================
StatusCode TrackBeamLineVertexFinderSoA::initialize() {
  auto sc = Transformer::initialize();

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
        double thisintegral                   = 0.5 * std::erf( sqrthalf * ( zminodd + ( ibin + 1 ) * m_dz ) / sigmaz );
        double bincontent                     = thisintegral - integral;
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

#ifdef TIMINGHISTOGRAMMING
  auto hsvc       = service<IHistogramSvc>( "HistogramDataSvc", true );
  m_timeperstepPr = hsvc->bookProf( name() + "/timeperstep", "time per step", 20, -0.5, 19.5 );
  m_timevsntrksPr = hsvc->bookProf( name() + "/timevsntrks", "time vs number of tracks", 50, -0.5, 249.5 );
  m_timevsnvtxPr  = hsvc->bookProf( name() + "/timevsnvtx", "time vs number of vertices", 12, -0.5, 11.5 );
#endif

  return sc;
}

//=============================================================================
// ::execute()
//=============================================================================

namespace {

  using namespace LHCb::Event;

  template <typename FTYPE>
  inline auto sqr( FTYPE x ) {
    return x * x;
  }
  //
  struct Extremum {
    Extremum( uint16_t _index, uint16_t _value, uint32_t _integral )
        : index{_index}, value{_value}, integral{_integral} {}
    uint16_t index;
    uint16_t value;
    uint32_t integral;
  };
  //
  struct Cluster {
    Cluster( uint16_t _izfirst, uint16_t _izlast, uint16_t _izmax )
        : izfirst{_izfirst}, izlast{_izlast}, izmax{_izmax} {}
    uint16_t izfirst;
    uint16_t izlast;
    uint16_t izmax;
    uint16_t ntracks{0};
  };

  // inline std::ostream& operator<<( std::ostream& os, Cluster const& c ) {
  //   os << "[" << c.izfirst << ", " << c.izlast << ", " << c.izmax << "]";
  //   return os;
  // }

  // we put this in local scope such that we can access it from various standalone routines
#ifdef TIMINGHISTOGRAMMING
  std::array<LHCb::TrackKernel::Timer, 20> timer;
  void                                     resettimers() {
    for ( auto& t : timer ) t.reset();
  }
#endif
} // namespace

namespace LHCb::Event::PV {

  namespace PVTrackTag {
    struct zerr : float_field {};   // error in coordinate at beam line
    struct blchi2 : float_field {}; // beamline chi2
    template <typename T>
    using pvseedingtrack_t = SOACollection<T, PVTrackTag::veloindex, PVTrackTag::pvindex, z, zerr, blchi2>;
  } // namespace PVTrackTag

  struct PVSeedingTracks : PVTrackTag::pvseedingtrack_t<PVSeedingTracks> {
    using base_t = typename PVTrackTag::pvseedingtrack_t<PVSeedingTracks>;
    using base_t::base_t;
  };

} // namespace LHCb::Event::PV

PrimaryVertexContainer TrackBeamLineVertexFinderSoA::operator()( const EventContext&           evtCtx,
                                                                 const LHCb::Pr::Velo::Tracks& tracksBackward,
                                                                 const LHCb::Pr::Velo::Tracks& tracksForward,
                                                                 const DeVP&                   vpdet ) const {
  // Get the beamline. this only accounts for position, not
  // rotation. that's something to improve! I have considered caching
  // this (with a handle for changes in the geometry, but the
  // computation is so fast that it isn't worth it.)
  // const auto beamline = Gaudi::XYZVector{beamline.X, beamline.Y, 0};
  // const Vec3<float_v> BL = Vec3<float_v>( beamline.X, beamline.Y, 0 );

  // Get the memory resource
  auto memResource = LHCb::getMemResource( evtCtx );

  // Step 1: select tracks with velo info, compute the poca to the
  // beamline. cache the covariance matrix at this position. I'd
  // rather us a combination of copy_if and transform, but don't know
  // how to do that efficiently.
#ifdef TIMINGHISTOGRAMMING
  resettimers();
  timer[9].start();
  timer[1].start();
#endif

  // Allow for correction for beamline position derived from conditions
  const auto beamspot  = vpdet.beamSpot();
  const auto beamlineX = beamspot.x() + m_beamLineOffsetX;
  const auto beamlineY = beamspot.y() + m_beamLineOffsetY;

  // actually we only need to store the (zbeam,veloindex) and fill a histogram. the rest we do not really need.
  PVSeedingTracks pvseedtracks;
  auto            pvseedtracks_simd = pvseedtracks.simd();
  pvseedtracks.reserve( tracksForward.size() + tracksBackward.size() );
  int icontainer{0};
  for ( const auto& tracks : {&tracksForward, &tracksBackward} ) {
    for ( auto const& track : tracks->simd() ) {
      auto loop_mask = track.loop_mask();
      auto index     = pvseedtracks.size();
      pvseedtracks.resize( index + popcount( loop_mask ) );

      // there must be a more efficient way to copy from one container to the other
      // zero indicates the state at the beamline
      Vec3<float_v> pos  = track.StatePos( 0 );
      Vec3<float_v> dir  = track.StateDir( 0 );
      Vec3<float_v> covX = track.StateCovX( 0 ); // covXX, covXTx, covTxTx (beats me why this is stored as a Vec3)
      Vec3<float_v> covY = track.StateCovY( 0 );

      // compute the z coordinate closest to the beamline and store it
      auto const dx    = beamlineX - pos.x;
      auto const dy    = beamlineY - pos.y;
      auto const tx    = dir.x;
      auto const ty    = dir.y;
      auto const dz    = select( loop_mask, ( tx * dx + ty * dy ) / ( tx * tx + ty * ty ), 0.f );
      auto const zbeam = pos.z + dz;

      auto const    Vx          = covX.x + 2 * dz * covX.y + dz * dz * covX.z;
      auto const    Vy          = covY.x + 2 * dz * covY.y + dz * dz * covY.z;
      float_v const Wx          = select( loop_mask, 1 / Vx, 1.f );
      float_v const Wy          = select( loop_mask, 1 / Vy, 1.f );
      float_v const zweight     = tx * Wx * tx + ty * Wy * ty;
      float_v const zerr        = select( loop_mask, 1 / sqrt( zweight ), 1.f ); // Q_rsqrt(zweight) ;
      float_v const blchi2      = dx * Wx * dx + dy * Wy * dy;
      auto          pvseedtrack = pvseedtracks_simd[index];
      pvseedtrack.field<PVTrackTag::z>().set( zbeam );
      pvseedtrack.field<PVTrackTag::zerr>().set( zerr );
      pvseedtrack.field<PVTrackTag::blchi2>().set( blchi2 );
      pvseedtrack.field<PVTrackTag::veloindex>().set( icontainer, track.indices() );
    }
    ++icontainer;
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
  // NOTE/TODO: As discussed in Rec#122, this could (and perhaps should) be drawn from an algorithm-local pool instead
  //            of the event-local pool used here. Alternatively, Wouter suggests that specific example could just be
  //            std::array, but this doesn't invalidate the algorithm-local idea!
  const int                                               nZBins = ( m_zmax - m_zmin ) / m_dz;
  boost::container::small_vector<uint16_t, defaultNZBins> zhisto( nZBins, 0 );
  const float                                             halfwindow = ( Nztemplatebins / 2 + 1 ) * m_dz;
  const float                                             zmin       = m_zmin + halfwindow;
  const float                                             zmax       = m_zmax - halfwindow;
  {
    for ( const auto pvtrack : pvseedtracks.simd() ) {
      auto const loop_mask = pvtrack.loop_mask();
      const auto zbeam     = pvtrack.get<PVTrackTag::z>();
      const auto zerr      = pvtrack.get<PVTrackTag::zerr>();
      const auto blchi2    = pvtrack.get<PVTrackTag::blchi2>();

      // bin in which z0 is, in floating point
      auto const jbin =
          int_v{4 * ( zbeam - m_zmin.value() ) / m_dz.value()}; // we need factor 4 to determine odd or even
      int_v const dbin = jbin & 3;                              // & 3 is the same as % 4
      auto const  minbin =
          to_std_array( ( jbin >> 2 ) - Nztemplatebins / 2 + select( dbin == 0, int_v( 0 ), int_v( 1 ) ) );
      auto const oddtemplate = to_std_array( select( ( dbin == 0 ) || ( dbin == 3 ), int_v( 0 ), int_v( 1 ) ) );
      //// make sure the template fits. make sure first and last bin
      //// remain empty. there will be lot's of funny effects at edge of
      //// histogram but we do not care.

      //// compute the bealine chi2. we make a cut for tracks contributing to seeding
      // (bl_x, bl_y)^T * W * (bl_x, bl_y) for diagonal W
      auto const zerrbin = int_v{Nztemplates / m_maxTrackZ0Err.value() * zerr};
      auto const mask    = to_std_array( int_v( float_v( zerrbin ) < Nztemplates && blchi2 < m_maxTrackBLChi2.value() &&
                                             zmin < zbeam && zbeam < zmax && loop_mask ) );
      auto const zerrbin_arr = to_std_array( zerrbin );
      /// add selected entries to the histogram.
      for ( std::size_t simd_idx = 0; simd_idx < simd::size; ++simd_idx ) {
        if ( mask[simd_idx] ) {
          for ( int j = 0; j < Nztemplatebins; ++j ) {
            zhisto[j + minbin[simd_idx]] += m_ztemplates[2 * zerrbin_arr[simd_idx] + oddtemplate[simd_idx]][j];
          }
        }
      }
    }
  }

#ifdef TIMINGHISTOGRAMMING
  timer[2].stop();
  timer[3].start();
#endif
  ////
  ////// Step 3: perform a peak search in the histogram. This used to be
  ////// very simple but the logic needed to find 'significant dips' made
  ////// it a bit more complicated. In the end it doesn't matter so much
  ////// because it takes relatively little time.
  ////
  boost::container::small_vector<Cluster, 32> clusters;
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
    boost::container::small_vector<BinIndex, 64> clusteredges;
    {
      bool     prevempty = true;
      uint32_t integral  = zhisto[0];
      for ( BinIndex i = 1; i < zhisto.size(); ++i ) {
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
      boost::container::small_vector<Extremum, 32> extrema;
      {
        bool     rising   = true;
        uint32_t integral = zhisto[ibegin];
        extrema.emplace_back( ibegin, zhisto[ibegin], integral );
        for ( uint16_t i = ibegin; i < iend; ++i ) {
          const auto value = zhisto[i];
          if ( value != zhisto[i + 1] ) {
            bool stillrising = zhisto[i + 1] > value;
            if ( stillrising != rising ) { // found a local extremum
              const uint8_t n          = extrema.size();
              const int16_t deltavalue = value - extrema.back().value;
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
        const auto                                  N = extrema.size();
        boost::container::small_vector<Cluster, 16> subclusters;
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

  // std::cout << "Number of clusters: " << clusters.size() << std::endl ;
#ifdef TIMINGHISTOGRAMMING
  timer[3].stop();
#endif

  // FIXME: we set up a lot of navigation below which is difficult to
  // fix if a PV is removed in the final selection step. Luckily this
  // happens so rarely that we can choose a lazy solution: if a PV is
  // removed, we just kill the corresponding seed and start over. We
  // can improve this in the future.
  PrimaryVertexContainer output{memResource};
  auto&                  vertices = output.vertices;
  auto&                  pvtracks = output.tracks;
  pvtracks.prvelocontainers[0]    = &tracksForward;
  pvtracks.prvelocontainers[1]    = &tracksBackward;

  if ( !clusters.empty() ) {

    // Step 4: partition the set of tracks by vertex seed: just choose
    // the closest one. we used to do this with std::partition, but that
    // requires that the tracks be reordered etc. instead, we now do
    // something else: each bin in the histogram is assigned to a PV. we
    // then simply recompute for each track in which bin it belongs, and
    // then choose the PV corresponding to that bin.

    // what makes this complicated is that we need to remove clusters that are too small.
#ifdef TIMINGHISTOGRAMMING
    timer[4].start();
#endif
    bool clusteringready = false;
    while ( !clusteringready && clusters.size() > 0 ) {
      boost::container::small_vector<int, defaultNZBins> pvindexhisto( nZBins, 0 );
      const auto                                         N          = clusters.size();
      const int                                          minpvindex = 0;
      const int                                          maxpvindex = N - 1;
      int                                                ibin       = 0;
      for ( int ipv = minpvindex; ipv < maxpvindex; ++ipv ) {
        // this will be the source of a small bias for nearby clusters
        const int zmid = ( clusters[ipv].izlast + clusters[ipv + 1].izfirst ) / 2;
        while ( ibin <= zmid ) pvindexhisto[ibin++] = ipv;
      }
      while ( ibin < nZBins ) pvindexhisto[ibin++] = maxpvindex;
      // assign all the tracks
      for ( const auto& pvtrack : pvseedtracks.simd() ) {
        auto ibin     = int_v{( pvtrack.get<PVTrackTag::z>() - m_zmin.value() ) / m_dz.value()};
        auto truncbin = max( 0, min( ibin, nZBins - 1 ) );
        auto pvindex  = gather( pvindexhisto.data(), truncbin );
        pvseedtracks.store<PVTrackTag::pvindex>( pvtrack.offset(), pvindex );
      }
      // count the number of tracks assigned to each vertex.
      for ( auto& clus : clusters ) clus.ntracks = 0;
      for ( const auto& pvtrack : pvseedtracks.scalar() )
        ++( clusters[pvtrack.get<PVTrackTag::pvindex>().cast()].ntracks );
      // remove clusters that are too small. this should happen rarely,
      // so it doesn't make sense to reuse the histogram that we already
      // had. we will first remove the smallest cluster, then restart
      // the loop.
      auto smallestcluster = std::min_element( clusters.begin(), clusters.end(),
                                               []( const auto& a, const auto& b ) { return a.ntracks < b.ntracks; } );
      clusteringready      = smallestcluster->ntracks >= m_minNumTracksPerVertex;
      if ( !clusteringready ) clusters.erase( smallestcluster );
    }
  }

  //// Step 5: Seed the primary vertices and assign tracks
  // Now we need to construct the output
#ifdef TIMINGHISTOGRAMMING
  timer[4].stop();
#endif
  if ( !clusters.empty() ) {
#ifdef TIMINGHISTOGRAMMING
    timer[5].start();
#endif
    // std::vector<Vertex, LHCb::Allocators::EventLocal<Vertex>> vertices{memResource};
    // First initialize the vertex parameters. I found that this funny
    // weighted 'maximum' is better than most other inexpensive
    // solutions.
    auto zClusterMean = [this, &zhisto]( auto izmax ) -> float {
      const uint16_t* b   = zhisto.data() + izmax;
      int             d1  = *b - *( b - 1 );
      int             d2  = *b - *( b + 1 );
      float           idz = d1 + d2 > 0 ? ( 0.5f * ( d1 - d2 ) ) / ( d1 + d2 ) : 0.0f;
      return m_zmin + m_dz * ( izmax + idz + 0.5f );
    };

    vertices.clear();
    vertices.reserve( clusters.size() );
    uint32_t trackoffset{0};
    for ( const auto& clus : clusters ) {
      auto& vtx = vertices.emplace_back( Gaudi::XYZPoint{beamlineX, beamlineY, zClusterMean( clus.izmax )} );
      vtx.setRange( trackoffset, trackoffset );
      trackoffset += clus.ntracks;
    }

    // Sort the tracks by PV. Cannot really vectorise this, but since it's got minimal instructions it is probably
    // fine
    pvtracks.reserve( pvseedtracks.size() + simd::size ); // make sure we have some extra
    pvtracks.resize( pvseedtracks.size() );
    // To prevent some trouble in the future, we'll explicitly initialize the PV index for the padding
    for ( size_t offset = 0; offset < pvtracks.size(); offset += simd::size )
      pvtracks.store<PVTrackTag::pvindex>( offset, int_v( 0 ) );

    auto& velotrackmap = output.velotrackmap;
    velotrackmap.resize( pvseedtracks.size(), -1 );
    auto pvtracksscalar = pvtracks.scalar();
    for ( const auto& pvseedtrack : pvseedtracks.scalar() ) {
      // read the position from the vertex
      const auto pvindex = pvseedtrack.get<PVTrackTag::pvindex>().cast();
      auto&      vertex  = vertices[pvindex];
      auto       offset  = vertex.end();
      // point the pv track to the right velo segment
      pvtracksscalar[offset].field<PVTrackTag::veloindex>().set(
          std::get<0>( pvseedtrack.field<PVTrackTag::veloindex>().index() ),
          std::get<1>( pvseedtrack.field<PVTrackTag::veloindex>().index() ) );
      pvtracksscalar[offset].field<PVTrackTag::pvindex>().set( pvindex );
      // pvtracksscalar[pvseedtrack.offset()].field<PVTrackTag::pvtrackindex>().set( offset ) ;
      velotrackmap[pvseedtrack.offset()] = offset;
      // update the vertex
      vertex.setSize( vertex.size() + 1 );
    }

#ifdef TIMINGHISTOGRAMMING
    timer[15].start();
#endif
    populateFromVeloTracks( output, tracksForward, tracksBackward );
    // populateFromVeloTracks( output );
#ifdef TIMINGHISTOGRAMMING
    timer[15].stop();
#endif

#ifdef TIMINGHISTOGRAMMING
    timer[5].stop();
    timer[6].start();
#endif
    // Call the fit. This can now be done separately for every
    // vertex. However, because a fit 'corrupts' the data of the
    // next vertex, we better always fit them in order.
    fitAdaptive( output.vertices, pvtracks, m_maxDeltaChi2, m_maxDeltaZConverged, m_maxDeltaChi2Converged,
                 m_maxFitIter );

#ifdef TIMINGHISTOGRAMMING
    timer[6].stop();
    timer[7].start();
#endif

    // Now perform two additional steps:
    // * merging (this seem especially needed on very busy events, like lead-lead)
    // * remove vertices with too little tracks or too large distance to beamline
    // If this leads to a change in the number of vertices, we will also perform another vertex fit

    // Merge vertices that are compatible with the previous vertex
    boost::container::small_vector<PrimaryVertex, 16> mergedvertices;
    mergedvertices.reserve( vertices.size() );
    mergedvertices.push_back( vertices.front() );
    for ( size_t i = 1; i < vertices.size(); ++i ) {
      const auto& ivertex = vertices[i];
      auto&       jvertex = mergedvertices.back();
      const auto  deltaz  = ivertex.position().z() - jvertex.position().z();
      const auto  icovz   = ivertex.covMatrix()( 2, 2 );
      const auto  jcovz   = jvertex.covMatrix()( 2, 2 );
      const auto  covz    = icovz + jcovz;
      if ( ( std::abs( deltaz ) < m_minVertexZSeparation ) ||
           ( deltaz * deltaz < covz * m_minVertexZSeparationChi2 ) ) {
        const auto begin = jvertex.begin();
        if ( jcovz > icovz ) jvertex = ivertex;
        jvertex.setRange( begin, ivertex.end() );
      } else {
        mergedvertices.push_back( ivertex );
      }
    }

    // Remove any vertices that do not pass the selection. We assign their tracks to the next or previous vertex,
    // depending on which one is closer.
    boost::container::small_vector<PrimaryVertex, 16> selectedvertices;
    selectedvertices.reserve( mergedvertices.size() );
    const auto maxVertexRho2 = sqr( m_maxVertexRho.value() );
    for ( size_t i = 0; i < mergedvertices.size(); ++i ) {
      const auto& vertex       = mergedvertices[i];
      const auto  beamlinedx   = vertex.position().x() - beamlineX;
      const auto  beamlinedy   = vertex.position().y() - beamlineY;
      const auto  beamlinerho2 = sqr( beamlinedx ) + sqr( beamlinedy );
      if ( vertex.nTracks() >= int( m_minNumTracksPerVertex ) && beamlinerho2 < maxVertexRho2 ) {
        selectedvertices.push_back( vertex );
      } else {
        const bool prevexists = selectedvertices.size() > 0;
        const bool nextexists = i + 1 < mergedvertices.size();
        if ( prevexists && nextexists ) {
          // add to the closest one
          const auto prevdz = selectedvertices.back().position().z() - vertex.position().z();
          const auto nextdz = mergedvertices[i + 1].position().z() - vertex.position().z();
          if ( std::abs( prevdz ) < std::abs( nextdz ) ) {
            selectedvertices.back().setRange( selectedvertices.back().begin(), vertex.end() );
          } else {
            mergedvertices[i + 1].setRange( vertex.begin(), mergedvertices[i + 1].end() );
          }
        } else if ( prevexists ) {
          selectedvertices.back().setRange( selectedvertices.back().begin(), vertex.end() );
        } else if ( nextexists ) {
          mergedvertices[i + 1].setRange( vertex.begin(), mergedvertices[i + 1].end() );
        }
      }
    }

    if ( selectedvertices.size() < vertices.size() ) {
      // Clear the output vertices and copy back the goodvertices
      vertices.clear();
      vertices.insert( vertices.begin(), selectedvertices.begin(), selectedvertices.end() );

      // Reassign the track->PV pointers
      auto pvtrackssimd = pvtracks.simd();
      for ( size_t i = 0; i < vertices.size(); ++i ) {
        const auto& vertex = vertices[i];
        const int_v pvindex_v{int( i )};
        for ( auto offset = vertex.begin(); offset < vertex.end(); offset += simd::size )
          pvtrackssimd[offset].field<PVTrackTag::pvindex>().set( pvindex_v );
      }

      // Refit
      for ( auto& vertex : vertices )
        fitAdaptive( vertex, pvtracks, m_maxDeltaChi2, m_maxDeltaZConverged, m_maxDeltaChi2Converged, m_maxFitIter );
    }
#ifdef TIMINGHISTOGRAMMING
    timer[7].stop();
#endif
  }

  // At some point we need to set the indices of the PVs.
  output.setIndices();

  // LHCb::Event::v2::RecVertices recvertexcontainer{memResource};
  m_nbPVsCounter += vertices.size();
  // std::cout << "TrackBeamLineVertexFinderSoA::operator end" << std::endl ;

  // PrimaryVertexContainer output{std::move( primaryvertices ), std::move( pvtracks )};
#ifdef TIMINGHISTOGRAMMING
  // timer[8].stop();
  timer[9].stop();
  for ( int i = 0; i < 20; ++i ) m_timeperstepPr->fill( float( i ), timer[i].total() );
  m_timevsntrksPr->fill( pvseedtracks.size(), timer[9].total() );
  m_timevsnvtxPr->fill( vertices.size(), timer[9].total() );
#endif
  return output;
}

// The converter:

#include "Event/RecVertex.h"
#include "Event/RecVertex_v2.h"
#include "Event/Track_v2.h"

namespace LHCb::Event::v2 {
  using Tracks = std::vector<Track>;
}

class PVToRecConverterV1 : public LHCb::Algorithm::Transformer<LHCb::RecVertices(
                               const LHCb::Event::PV::PrimaryVertexContainer&, const LHCb::Track::Range& )> {
public:
  PVToRecConverterV1( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"InputVertices", LHCb::Event::PV::DefaultLocation}, KeyValue{"InputTracks", ""}},
                     KeyValue{"OutputVertices", LHCb::RecVertexLocation::Primary} ) {}
  LHCb::RecVertices operator()( const LHCb::Event::PV::PrimaryVertexContainer& vertices,
                                const LHCb::Track::Range&                      keyed_tracks ) const override {
    // create the output container
    LHCb::RecVertices recvertexcontainer;
    recvertexcontainer.reserve( vertices.size() );
    const auto pvtracks = vertices.tracks.scalar();
    for ( const auto& vertex : vertices ) {
      auto recvertex = std::make_unique<LHCb::RecVertex>( vertex.position() );
      recvertex->setChi2( vertex.chi2() );
      recvertex->setNDoF( vertex.nDoF() );
      recvertex->setCovMatrix( vertex.covMatrix() );
      recvertex->setTechnique( LHCb::RecVertex::RecVertexType::Primary );
      // The following relies on the Velo tracks being in the right order ([forward,backward])
      const auto fwdsize = vertices.tracks.prvelocontainers[0]->size();
      for ( int trkindex = vertex.begin(); trkindex != vertex.end(); ++trkindex ) {
        const auto pvtrack = pvtracks[trkindex];
        const auto weight  = pvtrack.weight().cast();
        if ( weight > 0 ) {
          auto [containerIdx, veloIdx] = pvtrack.field<PVTrackTag::veloindex>().index();
          const int containerindex     = containerIdx.cast();
          const int velotrackindex     = veloIdx.cast() + ( containerindex == 0 ? 0 : fwdsize );
          if ( velotrackindex >= 0 && velotrackindex < int( keyed_tracks.size() ) )
            recvertex->addToTracks( keyed_tracks[velotrackindex], weight );
        }
      }
      recvertexcontainer.add( recvertex.release() );
    }
    return recvertexcontainer;
  }
};

class PVToRecConverterV2 : public LHCb::Algorithm::Transformer<LHCb::Event::v2::RecVertices(
                               const EventContext&,
                               // const LHCb::Pr::Velo::Tracks&,
                               // const LHCb::Pr::Velo::Tracks&,
                               const LHCb::Event::PV::PrimaryVertexContainer&, const LHCb::Event::v2::Tracks& )> {
public:
  PVToRecConverterV2( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {// KeyValue{"TracksBackwardLocation", "Rec/Track/VeloBackward"},
                      // KeyValue{"TracksLocation", "Rec/Track/Velo"},
                      KeyValue{"InputVertices", LHCb::Event::PV::DefaultLocation}, KeyValue{"InputTracks", ""}},
                     KeyValue{"OutputVertices", LHCb::Event::v2::RecVertexLocation::Primary} ) {}

  LHCb::Event::v2::RecVertices operator()( const EventContext&                            evtCtx,
                                           const LHCb::Event::PV::PrimaryVertexContainer& pvdata,
                                           // const LHCb::Pr::Velo::Tracks& backwardvelotracks,
                                           // const LHCb::Pr::Velo::Tracks& forwardvelotracks,
                                           const LHCb::Event::v2::Tracks& v2tracks ) const override {
    auto memResource = LHCb::getMemResource( evtCtx );
    // create the output container and copy the vertex data
    LHCb::Event::v2::RecVertices recvertexcontainer{memResource};
    recvertexcontainer.reserve( pvdata.vertices.size() );
    const auto pvtracks = pvdata.tracks.scalar();
    const auto fwdsize  = pvdata.tracks.prvelocontainers[0]->size();
    for ( const auto& vertex : pvdata.vertices ) {
      auto& recvertex = recvertexcontainer.emplace_back(
          vertex.position(), vertex.covMatrix(),
          LHCb::Event::v2::Track::Chi2PerDoF{vertex.chi2() / vertex.nDoF(), vertex.nDoF()} );
      recvertex.setTechnique( LHCb::Event::v2::RecVertex::RecVertexType::Primary );
      recvertex.reserve( vertex.nTracks() );
      for ( int trkindex = vertex.begin(); trkindex != vertex.end(); ++trkindex ) {
        const auto pvtrack = pvtracks[trkindex];
        const auto weight  = pvtrack.weight().cast();
        if ( weight > 0 ) {
          auto [containerIdx, veloIdx] = pvtrack.field<PVTrackTag::veloindex>().index();
          const int containerindex     = containerIdx.cast();
          const int velotrackindex     = veloIdx.cast() + ( containerindex == 0 ? 0 : fwdsize );
          recvertex.addToTracks( &( v2tracks[velotrackindex] ), weight );
        }
      }
    }
    return recvertexcontainer;
  }
};

DECLARE_COMPONENT( PVToRecConverterV1 )
DECLARE_COMPONENT( PVToRecConverterV2 )

class RecV1ToPVConverter
    : public LHCb::Algorithm::Transformer<LHCb::Event::PV::PrimaryVertexContainer( const LHCb::RecVertices& )> {
public:
  RecV1ToPVConverter( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator, {KeyValue{"InputVertices", LHCb::RecVertexLocation::Primary}},
                     KeyValue{"OutputVertices", LHCb::Event::PV::DefaultLocation} ) {}
  LHCb::Event::PV::PrimaryVertexContainer operator()( const LHCb::RecVertices& vertices ) const override {
    // create the output container
    LHCb::Event::PV::PrimaryVertexContainer pvcontainer;
    pvcontainer.vertices.reserve( vertices.size() );
    for ( const auto& vertex : vertices ) {
      auto& newpv = pvcontainer.vertices.emplace_back( vertex->position() );
      newpv.setNDoF( vertex->nDoF() );
      newpv.setChi2( vertex->chi2() );
      newpv.setCovMatrix( /* LHCb::LinAlg::convert(*/ vertex->covMatrix() );
    }
    pvcontainer.setIndices();
    return pvcontainer;
  }
};

DECLARE_COMPONENT( RecV1ToPVConverter )
