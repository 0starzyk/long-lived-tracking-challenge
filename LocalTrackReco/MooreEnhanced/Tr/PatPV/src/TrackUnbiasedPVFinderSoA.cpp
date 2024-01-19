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
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/PrVeloTracks.h"
#include "Event/RecVertex_v2.h"
#include "Event/Track.h"
#include "GaudiKernel/ToolHandle.h"
#include "Kernel/EventLocalAllocator.h"
#include "Kernel/STLExtensions.h"
#include "LHCbAlgs/Transformer.h"
#include "VPDet/DeVP.h"

#include "Event/TrackVertexUtils.h"
#include "LHCbMath/StateVertexUtils.h"

#ifdef TIMINGHISTOGRAMMING
#  include "AIDA/IProfile1D.h"
#  include "GaudiKernel/IHistogramSvc.h"
#  include "TrackKernel/Timer.h"
#endif

#include "Event/PrimaryVertices.h"

// boost includes
#include <boost/container/static_vector.hpp>

// std includes
#include <array>
#include <limits>
#include <vector>

#include <yaml-cpp/yaml.h>

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

  using simd    = SIMDWrapper::best::types;
  using float_v = simd::float_v;
  using int_v   = simd::int_v;

} // namespace

using namespace LHCb::Event::PV;

class TrackUnbiasedPVFinderSoA
    : public LHCb::Algorithm::Transformer<PrimaryVertexContainer( const EventContext&, const LHCb::Pr::Velo::Tracks&,
                                                                  const LHCb::Pr::Velo::Tracks&, const DeVP& ),
                                          LHCb::DetDesc::usesConditions<DeVP>> {
public:
  /// Standard constructor
  TrackUnbiasedPVFinderSoA( const std::string& name, ISvcLocator* pSvcLocator );
  /// Initialization
  StatusCode initialize() override;
  /// Execution
  /// Execution
  PrimaryVertexContainer operator()( const EventContext&, const LHCb::Pr::Velo::Tracks&, const LHCb::Pr::Velo::Tracks&,
                                     const DeVP& ) const override;

private:
  Gaudi::Property<uint32_t> m_minNumTracksPerVertex{this, "MinNumTracksPerVertex", 5};
  Gaudi::Property<float>    m_zmin{this, "MinZ", -300 * Gaudi::Units::mm, "Min z position of vertex seed"};
  Gaudi::Property<float>    m_zmax{this, "MaxZ", +300 * Gaudi::Units::mm, "Max z position of vertex seed"};
  Gaudi::Property<float>    m_maxTrackZ0Err{this, "MaxTrackZ0Err", 1.5 * Gaudi::Units::mm,
                                         "Maximum z0-error for adding track to histo"};
  Gaudi::Property<float>    m_maxVertexRho{this, "BeamSpotRCut", 0.3 * Gaudi::Units::mm,
                                        "Maximum distance of vertex to beam line"};
  Gaudi::Property<uint32_t> m_maxFitIter{this, "MaxFitIter", 10, "Maximum number of iterations for vertex fit"};
  Gaudi::Property<float> m_maxDeltaChi2{this, "MaxDeltaChi2", 12, "Maximum chi2 contribution of track to vertex fit"};
  Gaudi::Property<float> m_maxDeltaZConverged{this, "MaxDeltaZConverged", 0.001,
                                              "Limit on change in z to determine if vertex fit has converged"};
  Gaudi::Property<float> m_maxDeltaChi2Converged{this, "MaxDeltaChi2Converged", 0.01,
                                                 "Limit on change in chi2 to determine if vertex fit has converged"};
  Gaudi::Property<float> m_minVertexZSeparationChi2{this, "MinVertexZSeparationChi2", 25};
  Gaudi::Property<float> m_minVertexZSeparation{this, "MinVertexZSeparation", 1.0f * Gaudi::Units::mm};
  Gaudi::Property<float> m_beamLineOffsetX{this, "BeamLineOffsetX", 0.0f, "X correction applied to beamline position"};
  Gaudi::Property<float> m_beamLineOffsetY{this, "BeamLineOffsetY", 0.0f, "Y correction applied to beamline position"};
  mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_nbPVsCounter{this, "Nb PVs"};

#ifdef TIMINGHISTOGRAMMING
  AIDA::IProfile1D* m_timeperstepPr{nullptr};
  AIDA::IProfile1D* m_timevsntrksPr{nullptr};
  AIDA::IProfile1D* m_timevsnvtxPr{nullptr};
#endif
};

DECLARE_COMPONENT( TrackUnbiasedPVFinderSoA )

TrackUnbiasedPVFinderSoA::TrackUnbiasedPVFinderSoA( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator,
                   {KeyValue{"TracksBackwardLocation", "Rec/Track/VeloBackward"},
                    KeyValue{"TracksLocation", "Rec/Track/Velo"},
                    KeyValue{"BeamSpotLocation", "AlgorithmSpecific-" + name + "-beamspot"}},
                   KeyValue{"OutputVertices", LHCb::Event::v2::RecVertexLocation::Primary} ) {}

StatusCode TrackUnbiasedPVFinderSoA::initialize() {
  return Transformer::initialize().andThen( [&] {
#ifdef TIMINGHISTOGRAMMING
    auto hsvc       = service<IHistogramSvc>( "HistogramDataSvc", true );
    m_timeperstepPr = hsvc->bookProf( name() + "/timeperstep", "time per step", 20, -0.5, 19.5 );
    m_timevsntrksPr = hsvc->bookProf( name() + "/timevsntrks", "time vs number of tracks", 50, -0.5, 249.5 );
    m_timevsnvtxPr  = hsvc->bookProf( name() + "/timevsnvtx", "time vs number of vertices", 12, -0.5, 11.5 );
#endif
    return StatusCode::SUCCESS;
  } );
}

namespace {

  using namespace LHCb::Event;

  template <typename FTYPE>
  inline auto sqr( FTYPE x ) {
    return x * x;
  }

  // we put this in local scope such that we can access it from various standalone routines
#ifdef TIMINGHISTOGRAMMING
  std::array<Timer, 20> timer;
  void                  resettimers() {
    for ( auto& t : timer ) t.reset();
  }
#endif
} // namespace

namespace LHCb::Event::PV {
  namespace PVTrackTag {

    // struct veloindex : int_field {}; // index to the corresponding velo track
    // struct pvindex : int_field {};   // index to the corresponding pv

    struct status : int_field {};
    enum struct Status { Unused, Used };
    struct t2 : float_field {}; // error in coordinate at beam line
    struct ip2 : float_field {};
    struct zweight : float_field {};
    template <typename T>
    using pvseedingtrack_t = SOACollection<T, z, x, y, tx, ty, Vx, Vy, ip2, ipchi2, pvindex, veloindex, t2, zweight>;

  } // namespace PVTrackTag

  struct PVSeedingTracks : PVTrackTag::pvseedingtrack_t<PVSeedingTracks> {
    using base_t = typename PVTrackTag::pvseedingtrack_t<PVSeedingTracks>;
    using base_t::base_t;
  };
} // namespace LHCb::Event::PV

namespace {

  struct LightStateVector {
    float       m_z;
    float       m_x;
    float       m_y;
    float       m_tx;
    float       m_ty;
    const auto& z() const { return m_z; }
    const auto& x() const { return m_x; }
    const auto& y() const { return m_y; }
    const auto& tx() const { return m_tx; }
    const auto& ty() const { return m_ty; }
    template <typename PVTrack>
    LightStateVector( const PVTrack& track )
        : m_z{track.template field<PVTrackTag::z>().get().cast()}
        , m_x{track.template field<PVTrackTag::x>().get().cast()}
        , m_y{track.template field<PVTrackTag::y>().get().cast()}
        , m_tx{track.template field<PVTrackTag::tx>().get().cast()}
        , m_ty{track.template field<PVTrackTag::ty>().get().cast()} {}
  };

  struct LightState : public LightStateVector {
    Gaudi::SymMatrix4x4F m_cov;
    const auto&          covariance() const { return m_cov; }
    template <typename PVTrack>
    LightState( const PVTrack& track ) : LightStateVector{track} {
      m_cov( 0, 0 ) = track.template field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::xx ).get().cast();
      m_cov( 2, 0 ) = track.template field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::xtx ).get().cast();
      m_cov( 2, 2 ) = track.template field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::txtx ).get().cast();
      m_cov( 1, 1 ) = track.template field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::xx ).get().cast();
      m_cov( 3, 1 ) = track.template field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::xtx ).get().cast();
      m_cov( 3, 3 ) = track.template field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::txtx ).get().cast();
    }
  };

  LHCb::State toState( const LightState& in ) {
    LHCb::State out;
    out.setX( in.x() );
    out.setY( in.y() );
    out.setZ( in.z() );
    out.setTx( in.tx() );
    out.setTy( in.ty() );
    for ( int irow = 0; irow < 4; ++irow )
      for ( int icol = 0; icol <= irow; ++icol ) out.covariance()( irow, icol ) = in.covariance()( irow, icol );
    return out;
  }

  template <typename PVTrack>
  LightStateVector makestate( const PVTrack& track ) {
    return LightStateVector{track};
  }

  struct VertexSeed {
    Gaudi::XYZPoint     position;
    Gaudi::SymMatrix3x3 covariance;
    int                 index{-1};
    size_t              ntracks{0};
    size_t              nclosetracks{0};
    float               chi2{0};
  };
} // namespace

PrimaryVertexContainer TrackUnbiasedPVFinderSoA::operator()( const EventContext&           evtCtx,
                                                             const LHCb::Pr::Velo::Tracks& tracksBackward,
                                                             const LHCb::Pr::Velo::Tracks& tracksForward,
                                                             const DeVP&                   deVP ) const {
  /*

    Some observations:
    - two-track seeding requires two-track combinatorics
    - ways to reduce the number of combinations
        - limiting the number of tracks that participate (cut on angle, or sigma(z); cannot make it very tight)
        - flagging used tracks
        - do some sort of sorting to get to NlogN. just don't know how to do that without beamline.
    - PVs are only separated in z. SVs are separated in z from their PV. so actually, vertices are always separated in
    z.
    - sorting tracks was about as slow as the vertex fit. sorting is NlogN
    - the standard alg assigns about 90% of tracks to a vertex: the remaining 10% is not compatible with 4-prong
    vertices.
    - the really slow events will be those in which there are no seeds: many incompatible tracks

    Requirements:
    - a seed has at least two tracks that are incompatible with all other seeds (dchi2>10)
    - a seed has a maximum chi2 of 4 (1?!)
    - a seed has at least N tracks with delta-chi2 < 4
    - a seed has a small error in z. tracks mst make large angle with each other. something will small eigenvalues of
    cov matrix


    Lot's of random toughts:
    - we can create seeds and partition simultaneously. just not sure that that's much gain.
    - we can give a track 3 states: 0) unused ; 1) compatible with seed ; 2) close to seed ; 3) in seed
       2) or 3) probably does not matter, or does it?
    - only category 0 is considered for seeding.
    - category 0) and 1) are open for repartitioning: cat 0) is expected to be about 10%.
    - the repartitioning can be done dring the combinatorics loop if we also store the delta-chi2

    - we'll flag tracks as used quicker, if we can more quickly find the good seeds. so perhaps sorting tracks in 'z'
    still makes sense.

    - a seed with large errors will absorbe poor tracks. so we need to find the small error seeds first

    - parallel track make fore poor seeds.

    - we can do a histogramming method with seed z. however, that would still require full track-seed combinatorics
    afterwards. so, not sure we would gain.

    - let's just sort tracks in tx^2+ty^2. then require for seeds a maximum value of (tx1*tx2 + ty1*ty2 + 1), or a
    minimum valueof dt^2 = dtx^2 + dty^2.

    - make a histogram of
       - dt for every pair of tracks in a PV
       - max dt in a PV. I hope this has a minimum at about 0.1 rad.


   */

  // Get the beamline. this only accounts for position, not
  // rotation. that's something to improve!

  // Get the memory resource
  auto memResource = LHCb::getMemResource( evtCtx );

  // Allow for correction for beamline position derived from conditions
  Gaudi::XYZPoint beamline  = deVP.beamSpot();
  const auto      beamlineX = beamline.x() + m_beamLineOffsetX;
  const auto      beamlineY = beamline.y() + m_beamLineOffsetY;

  // Step 1: select tracks with velo info, compute the poca to the
  // beamline. cache the covariance matrix at this position. I'd
  // rather us a combination of copy_if and transform, but don't know
  // how to do that efficiently.
#ifdef TIMINGHISTOGRAMMING
  resettimers();
  timer[9].start();
  timer[1].start();
#endif

  // actually we only need to store the (zbeam,veloindex) and fill a histogram. the rest we do not really need.
  PVSeedingTracks pvseedtracks;
  const size_t    Nvelo = tracksForward.size() + tracksBackward.size();
  pvseedtracks.reserve( Nvelo + simd::size ); // reserve one extra for overflow of the padding when filling?
  // make sure to set the padding of the pv index since we need that to be valid for all simd values
  pvseedtracks.simd()[simd::size * ( Nvelo / simd::size )].field<PVTrackTag::pvindex>().set( int_v{0} );

  int icontainer{0};
  for ( const auto& tracks : {&tracksForward, &tracksBackward} ) {
    // const size_t cursize = pvseedtracks.size();
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

      // unlike TBLFV ignore the extrapolation to the beam line. however, we need the track parameters to determine the
      // vtx position

      // compute the z coordinate closest to the beamline and store it
      auto const tx = dir.x;
      auto const ty = dir.y;
      auto const Vx = covX.x; // + 2 * dz * covX.y + dz * dz * covX.z; // we could also choose the 'minimal' value here
      auto const Vy = covY.x; // + 2 * dz * covY.y + dz * dz * covY.z;
      auto const Wx = select( loop_mask, 1 / Vx, 0.001f );
      auto const Wy = select( loop_mask, 1 / Vy, 0.001f );
      auto const zweight = tx * Wx * tx + ty * Wy * ty;
      // auto const detperpcov = ( Wx * (1+tx*tx) * Wy *  (1+ty*ty) ) ;

      // auto const zerr2    = select( loop_mask, 1 / zweight , 1.f ); // Q_rsqrt(zweight) ;
      auto const t2          = tx * tx + ty * ty;
      auto       pvseedtrack = pvseedtracks.simd()[index];
      pvseedtrack.field<PVTrackTag::t2>().set( t2 );
      pvseedtrack.field<PVTrackTag::zweight>().set( zweight );
      pvseedtrack.field<PVTrackTag::pvindex>().set( int_v{0} );
      pvseedtrack.field<PVTrackTag::ipchi2>().set( float_v{9999.} );
      pvseedtrack.field<PVTrackTag::ip2>().set( float_v{9999.} );
      pvseedtrack.field<PVTrackTag::z>().set( pos.z );
      pvseedtrack.field<PVTrackTag::x>().set( pos.x );
      pvseedtrack.field<PVTrackTag::y>().set( pos.y );
      pvseedtrack.field<PVTrackTag::tx>().set( dir.x );
      pvseedtrack.field<PVTrackTag::ty>().set( dir.y );
      pvseedtrack.field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::xx ).set( covX.x );
      pvseedtrack.field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::xtx ).set( covX.y );
      pvseedtrack.field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::txtx ).set( covX.z );
      pvseedtrack.field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::xx ).set( covY.x );
      pvseedtrack.field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::xtx ).set( covY.y );
      pvseedtrack.field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::txtx ).set( covY.z );
      pvseedtrack.field<PVTrackTag::veloindex>().set( icontainer, track.indices() );
    }
    ++icontainer;
  }

#ifdef TIMINGHISTOGRAMMING
  timer[1].stop();
  timer[2].start();
#endif

  // Step 2: order the tracks in t2
  std::vector<std::pair<int, float>> orderedtracks;
  auto                               pvseedtracks_scalar = pvseedtracks.scalar();
  for ( auto tr : pvseedtracks_scalar ) orderedtracks.emplace_back( tr.offset(), tr.get<PVTrackTag::zweight>().cast() );
  std::sort( orderedtracks.begin(), orderedtracks.end(),
             []( const auto& lhs, const auto& rhs ) { return lhs.second > rhs.second; } );

#ifdef TIMINGHISTOGRAMMING
  timer[2].stop();
  timer[3].start();
#endif
  // Step 3:
  std::vector<VertexSeed> seeds;
  const auto              N                 = orderedtracks.size();
  auto                    pvseedtracks_simd = pvseedtracks.simd();

  // const float m_minSeedTrackDistToOtherSeeds = 0.4 * Gaudi::Units::mm;
  const float m_maxSeedDoca = 0.2 * Gaudi::Units::mm;
  // const float mindoca2                       = m_minSeedTrackDistToOtherSeeds * m_minSeedTrackDistToOtherSeeds;
  const float m_minSeedTrackChi2ToOtherSeeds = 25.0; // m_maxDeltaChi2 ;
  const float mindt2                         = 0.05 * 0.05;
  const float m_maxSeedChi2                  = 1;
  const float m_maxSeedZError                = 0.5 * Gaudi::Units::mm;
  const float m_minTrackZWeight              = 0; // 1/m_maxTrackZ0Err ;

  for ( unsigned int i = 0; i < N && orderedtracks[i].second > m_minTrackZWeight; ++i ) {
    const auto trackA = pvseedtracks_scalar[orderedtracks[i].first];
    if ( /*trackA.field<PVTrackTag::ip2>().get() > mindoca2*/
         trackA.field<PVTrackTag::ipchi2>().get() > m_minSeedTrackChi2ToOtherSeeds ) {
      const auto stateA = LightState{trackA};

      // this loop can in theory be vectorized
      for ( unsigned int j = 0; j < i; ++j ) {
        const auto trackB = pvseedtracks_scalar[orderedtracks[j].first];

        if ( /*trackB.field<PVTrackTag::ip2>().get() > mindoca2*/
             trackB.field<PVTrackTag::ipchi2>().get() > m_minSeedTrackChi2ToOtherSeeds ) {
          // compute a vertex. this needs to be fast. we would like to make some cuts on vertex quality, etc.
          // we could also first cut on the doca, to make this faster.
          const auto stateB = LightState{trackB};

          // cut on the opening angle of the two tracks
          const auto dtx = trackA.field<PVTrackTag::tx>().get() - trackB.field<PVTrackTag::tx>().get();
          const auto dty = trackA.field<PVTrackTag::ty>().get() - trackB.field<PVTrackTag::ty>().get();
          const auto dt2 = dtx * dtx + dty * dty;
          if ( dt2 > mindt2 ) {
            // if we just need the position, then a vertex fit is actually not needed

            Gaudi::XYZPoint     seedpos;
            Gaudi::XYZPoint     vertexpos;
            Gaudi::SymMatrix3x3 vertexweight;
            Gaudi::SymMatrix3x3 vertexcov;
            // const auto chi2 = StateVertexUtils::vertex(stateA, stateB, seedpos, vertexweight, vertexcov ) ;
            const auto doca = LHCb::StateVertexUtils::doca( stateA, stateB );
            LHCb::StateVertexUtils::poca( stateA, stateB, seedpos );

            // make a requirement on the vertex z error. would perhaps like to use 'density' here.
            // minuslogdensity = 0.5*chi2 + 0.5*log( det( vertexcov ) )

            if ( std::abs( doca ) < m_maxSeedDoca ) {

              const auto chi2 = LHCb::StateVertexUtils::vertexChi2( stateA, stateB );
              if ( chi2 < m_maxSeedChi2 ) {

                const auto chi2alt = LHCb::TrackVertexUtils::vertex( toState( stateA ), toState( stateB ), vertexpos,
                                                                     vertexweight, vertexcov );

                if ( vertexcov( 2, 2 ) < m_maxSeedZError * m_maxSeedZError ) {
                  /*
                    std::cout << "Accepted seed: "
                              << seeds.size() << " " << i << " " << j << " "
                              << trackA.field<PVTrackTag::ip2>().get() << " "
                              << trackB.field<PVTrackTag::ip2>().get() << " : "
                              << doca << " " << chi2 << " " << chi2alt << " : "
                              << seedpos << " " << vertexpos << " " << std::sqrt( vertexcov(2,2))
                              << std::endl ;
                  */

                  // flag these two tracks
                  const auto iseed = seeds.size();
                  trackA.field<PVTrackTag::pvindex>().set( iseed );
                  trackB.field<PVTrackTag::pvindex>().set( iseed );
                  // trackA.field<PVSeedingTrackTag::ip2>().set( 0.5*chi2 ) ; // this is not quite right, right?
                  // trackB.field<PVSeedingTrackTag::ip2>().set( 0.5*chi2 ) ;

                  //
                  VertexSeed seed;
                  seed.position   = vertexpos;
                  seed.covariance = vertexcov;
                  seed.index      = iseed;
                  seed.chi2       = chi2alt;
                  seeds.emplace_back( seed );

                  // this loop can be vectorized!
                  for ( auto ktrack : pvseedtracks_simd ) {
                    // const auto unused = ktrack.field<PVSeedingTrackTag::pvindex>.get() < 0 ;
                    // compute the ip2 to the seed. I think that doca is actually better, but harder to decide on a
                    // cut.

                    // very expensive to compute an LHCb::State here. There must be a vectorized form already.
                    // const auto state  = ktrack.state() ;
                    // what is even cheaper is to compute the ip2 ignoring the vertex error. why don't we just choose
                    // the closest?!

                    // const auto ip2 = TrackVertexUtils::vertexChi2( state, vertexpos, vertexcov ) ;

                    // temporarily use ip2 = docaXY. we really don't want to use errors.
                    const auto dz = seedpos.z() - ktrack.field<PVTrackTag::z>().get();
                    const auto dx =
                        ktrack.field<PVTrackTag::x>().get() + dz * ktrack.field<PVTrackTag::tx>().get() - seedpos.x();
                    const auto dy =
                        ktrack.field<PVTrackTag::y>().get() + dz * ktrack.field<PVTrackTag::ty>().get() - seedpos.y();
                    const auto ip2         = dx * dx + dy * dy;
                    const auto ip2best     = ktrack.field<PVTrackTag::ip2>().get();
                    const auto pvindexbest = ktrack.field<PVTrackTag::pvindex>().get();
                    const auto closer      = ip2 < ip2best;
                    if ( popcount( closer ) ) {
                      ktrack.field<PVTrackTag::ip2>().set( select( closer, ip2, ip2best ) );
                      ktrack.field<PVTrackTag::pvindex>().set( select( closer, iseed, pvindexbest ) );
                      // let's also compute the IPchi2 to the selected PV. at this point we should wonder why we don't
                      // just run the vertex fit.
                      const auto Vx =
                          ktrack.field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::xx ).get() +
                          2 * dz * ktrack.field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::xtx ).get() +
                          dz * dz * ktrack.field<PVTrackTag::Vx>( PVTrackTag::XTxCovMatrixElement::txtx ).get();
                      const auto Vy =
                          ktrack.field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::xx ).get() +
                          2 * dz * ktrack.field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::xtx ).get() +
                          dz * dz * ktrack.field<PVTrackTag::Vy>( PVTrackTag::XTxCovMatrixElement::txtx ).get();
                      const auto ipchi2     = dx * dx / Vx + dy * dy / Vy;
                      const auto ipchi2best = ktrack.field<PVTrackTag::ipchi2>().get();
                      ktrack.field<PVTrackTag::ipchi2>().set( select( closer, ipchi2, ipchi2best ) );
                    }
                  }
                  break;
                }
              }
            }
          }
        }
      }
    }
  }

  // step 3b: prune the seeds to remove duplicates? or just fit and remove afterwards?

  /*
    - remove any seed that has less than N tracks assigned

  */
  if ( !seeds.empty() ) {
    // First do some counting, for which we still need the seeds in the original order
    for ( const auto ktrack : pvseedtracks_scalar ) {
      auto& seed = seeds[ktrack.field<PVTrackTag::pvindex>().get().cast()];
      seed.ntracks += 1;
      if ( ktrack.field<PVTrackTag::ipchi2>().get().cast() < m_maxDeltaChi2 ) seed.nclosetracks += 1;
    }

    // Next sort the seeds in z
    std::sort( seeds.begin(), seeds.end(),
               []( const auto& lhs, const auto& rhs ) { return lhs.position.z() < rhs.position.z(); } );

    // now delete all the seeds that are not acceptable because they have too few tracks left
    ////std::cout << "Seeds: " << seeds.size() << std::endl ;
    // for( const auto& seed : seeds )
    //  std::cout << seed.position <<  " " << seed.ntracks << " " << seed.nclosetracks << " " <<
    //  std::sqrt(seed.covariance(2,2)) << std::endl ;

    // std::cout << "Creating map" << std::endl ;
    // Map from the original index to the index in the sorted PV list

    // Map from index in the sorted list to an index in the sorted list for an accepted vertex
    const size_t     N = seeds.size();
    std::vector<int> indexmap( N, 0 );
    for ( size_t i = 0; i < N; ++i ) indexmap[i] = i;

    // Map from sorted index to the index in the selected PV list
    std::vector<int> selindexmap( seeds.size(), -1 );

    // Apply some sort of selection: to be implemented
    std::vector<VertexSeed> selectedseeds;
    selectedseeds.reserve( N );
    for ( size_t i = 0; i < N; ++i ) {
      const auto& iseed = seeds[i];
      // compute the distance to the last selected vertex
      const bool prevexists = selectedseeds.size() > 0;
      // make sure that this one is not too close to the previous seed. if so, assign it to that.
      if ( prevexists ) {
        auto&      jseed  = selectedseeds.back();
        const auto deltaz = jseed.position.z() - iseed.position.z();
        const auto covz   = jseed.covariance( 2, 2 ) + iseed.covariance( 2, 2 );
        if ( std::abs( deltaz ) < m_minVertexZSeparation || deltaz * deltaz < covz * m_minVertexZSeparationChi2 ) {
          if ( jseed.nclosetracks < iseed.nclosetracks ) jseed = iseed;
          selindexmap[i] = selectedseeds.size() - 1;
        }
      }
      if ( selindexmap[i] < 0 ) {
        if ( iseed.nclosetracks >= m_minNumTracksPerVertex ) {
          // add it to the list of accepted seeds
          selectedseeds.push_back( iseed );
          selindexmap[i] = selectedseeds.size() - 1;
        } else {
          // merge with the closest one
          const bool nextexists = i + 1 < N;
          if ( prevexists && nextexists ) {
            // add to the closest one
            auto&      jseed  = selectedseeds.back();
            const auto prevdz = jseed.position.z() - iseed.position.z();
            const auto nextdz = seeds[i + 1].position.z() - iseed.position.z();
            if ( std::abs( prevdz ) < std::abs( nextdz ) ) {
              selindexmap[i] = selectedseeds.size() - 1;
            } else {
              indexmap[i] = i + 1;
            }
          } else if ( prevexists ) {
            selindexmap[i] = selectedseeds.size() - 1;
          } else if ( nextexists ) {
            indexmap[i] = i + 1;
          }
        }
      }
    }

    // update the indexmap to the selected seeds
    if ( selectedseeds.size() > 0 ) {
      // std::cout << "Ready to do something with the seeds maps: " << seeds.size() << " " << selectedseeds.size() <<
      // std::endl ;

      // If a vertex was not selected, it should point to a vertex
      // in the selected list. Sometimes, it needs to be
      // redirected. That's what this loop does:
      for ( int i = 0; i < int( N ); ++i ) {
        if ( selindexmap[i] < 0 ) {
          if ( indexmap[i] == i ) {
            std::cout << "TrackUnbiasedPVFinder Serious problem: " << i << " " << selindexmap[i] << std::endl;
            break;
          }
          auto j = indexmap[i];
          while ( selindexmap[j] < 0 ) j = indexmap[j];
          selindexmap[i] = selindexmap[j];
        }
      }
      // std::cout << "ready redirecting" << std::endl ;
      // Now everything is indexed in terms of the position in the
      // sorted list. However, we need it for the original
      // index. That's what we do here:
      std::vector<int> sortedselindexmap( seeds.size(), 0 );
      for ( size_t i = 0; i < N; ++i ) sortedselindexmap[seeds[i].index] = selindexmap[i];
      // std::cout << "final check of the map" << std::endl ;
      // Check that the selindexmap is fine
      for ( size_t i = 0; i < N; ++i ) {
        if ( sortedselindexmap[i] < 0 || sortedselindexmap[i] >= int( selectedseeds.size() ) )
          std::cout << "Something is wrong with the map: " << i << " " << sortedselindexmap[i] << std::endl;
      }
      // Make sure to reassign the tracks
      // std::cout << "Reassign tracks" << std::endl ;
      // for(int i=0; i<seeds.size(); ++i) std::cout << i << " --> " << indexmap[ i ] << std::endl ;
      for ( const auto ktrack : pvseedtracks_simd ) {
        // std::cout << ktrack.field<PVTrackTag::pvindex>().get() << std::endl ;
        auto newpvindex = gather( sortedselindexmap.data(), ktrack.field<PVTrackTag::pvindex>().get() );
        ktrack.field<PVTrackTag::pvindex>().set( newpvindex );
      }
      // Make sure to update the pv index
      for ( size_t i = 0; i < selectedseeds.size(); ++i ) selectedseeds[i].index = i;

      // Check that every track points to a vertex
      for ( const auto ktrack : pvseedtracks_scalar ) {
        const auto pvindex = ktrack.field<PVTrackTag::pvindex>().get().cast();
        if ( pvindex < 0 || pvindex >= int( selectedseeds.size() ) ) {
          std::cout << "Track points nowhere: " << ktrack.offset() << " " << pvindex << std::endl;
        }
      }

      // update the counters
      for ( size_t i = 0; i < selectedseeds.size(); ++i ) selectedseeds[i].ntracks = 0;
      for ( const auto ktrack : pvseedtracks_scalar ) {
        auto& seed = selectedseeds[ktrack.field<PVTrackTag::pvindex>().get().cast()];
        seed.ntracks += 1;
      }
    }
    seeds = selectedseeds;
    // std::cout << "Selected seeds: " << seeds.size() << std::endl ;
    // for( const auto& seed : seeds )
    //  std::cout << seed.position <<  " " << seed.ntracks << " " << seed.nclosetracks << " " <<
    //  std::sqrt(seed.covariance(2,2)) << std::endl ;
  }

#ifdef TIMINGHISTOGRAMMING
  timer[3].stop();
#endif

  // Create the output
  PrimaryVertexContainer output{memResource};
  auto&                  vertices = output.vertices;
  auto&                  pvtracks = output.tracks;
  pvtracks.prvelocontainers[0]    = &tracksForward;
  pvtracks.prvelocontainers[1]    = &tracksBackward;

  if ( !seeds.empty() ) {
#ifdef TIMINGHISTOGRAMMING
    timer[5].start();
#endif
    // std::vector<Vertex, LHCb::Allocators::EventLocal<Vertex>> vertices{memResource};
    // First initialize the vertex parameters. I found that this funny
    // weighted 'maximum' is better than most other inexpensive
    // solutions.

    vertices.reserve( seeds.size() );
    uint32_t trackoffset{0};
    for ( const auto& clus : seeds ) {
      auto& vtx = vertices.emplace_back( clus.position );
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

    // Check that every track points to a velo track
    for ( const auto& pvtrack : pvtracks.scalar() ) {
      auto [containerIdx, veloIdx] = pvtrack.field<PVTrackTag::veloindex>().index();
      if ( containerIdx.cast() < 0 || containerIdx.cast() > 1 || veloIdx.cast() < 0 ||
           veloIdx.cast() >= int( pvtracks.prvelocontainers[containerIdx.cast()]->size() ) )
        std::cout << "Pv track not assigned properly: " << containerIdx << " " << veloIdx << std::endl;
    }

    populateFromVeloTracks( output, tracksForward, tracksBackward );

#ifdef TIMINGHISTOGRAMMING
    timer[5].stop();
    timer[6].start();
#endif
    // Call the fit. This can now be done separately for every
    // vertex. However, because a fit may 'corrupt' the data of the
    // next vertex, we better always fit them in order.
    fitAdaptive( output.vertices, pvtracks, m_maxDeltaChi2, m_maxDeltaZConverged, m_maxDeltaChi2Converged, m_maxFitIter,
                 1 );

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
  // std::cout << "TrackUnbiasedPVFinderSoA::operator end" << std::endl ;

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
