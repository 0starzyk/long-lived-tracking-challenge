/***************************************************************************** \
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "PrUTMagnetTool.h"
#include <DetDesc/GenericConditionAccessorHolder.h>
#include <Event/PrHits.h>
#include <Event/PrUpstreamTracks.h>
#include <Event/PrVeloTracks.h>
#include <Event/UTSectorHelper.h>
#include <Event/UTTrackUtils.h>
#include <GaudiAlg/ISequencerTimerTool.h>
#include <LHCbAlgs/Transformer.h>
#include <LHCbMath/GeomFun.h>
#include <Magnet/DeMagnet.h>
#include <PrKernel/PrMutUTHits.h>
#include <UTDAQ/UTDAQHelper.h>
#include <UTDAQ/UTInfo.h>
#include <vdt/sqrt.h>

/** @class PrVeloUT PrVeloUT.h
 *
 *  PrVeloUT algorithm. This is just a wrapper,
 *  the actual pattern recognition is done in the 'PrVeloUTTool'.
 *
 *  - InputTracksName: Input location for Velo tracks
 *  - OutputTracksName: Output location for VeloTT tracks
 *  - TimingMeasurement: Do a timing measurement?
 *
 *  @author Mariusz Witek
 *  @date   2007-05-08
 *  @update for A-Team framework 2007-08-20 SHM
 *
 *  2017-03-01: Christoph Hasse (adapt to future framework)
 *  2019-04-26: Arthur Hennequin (change data Input/Output)
 *  2020-08-26: Peilian Li (change data Input/Output to SOACollection)
 */

namespace LHCb::Pr {

  using simd   = SIMDWrapper::best::types;
  using scalar = SIMDWrapper::scalar::types;

  constexpr auto nanMomentum = std::numeric_limits<float>::quiet_NaN();

  constexpr static int                  batchSize  = 2 * simd::size;
  [[maybe_unused]] constexpr static int maxNumCols = 3; // if needed, algo can be templated with this
  [[maybe_unused]] constexpr static int maxNumRows = 3; // if needed, algo can be templated with this
  [[maybe_unused]] constexpr static int maxNumSectors =
      maxNumCols * maxNumRows; // if needed, algo can be templated with this
  const int totalUTLayers = static_cast<int>( UTInfo::DetectorNumbers::TotalLayers );

  struct ProtoTracks final {

    std::array<float, simd::size> wbs;
    std::array<float, simd::size> xMidFields;
    std::array<float, simd::size> invKinkVeloDists;

    // -- this is for the hits
    // -- this does _not_ include overlap hits, so only 4 per track
    std::array<float, 4 * batchSize> xs;
    std::array<float, 4 * batchSize> zs;
    std::array<float, 4 * batchSize> weightss{}; // this needs to be zero-initialized
    std::array<float, 4 * batchSize> sins;
    std::array<int, 4 * batchSize>   ids;
    std::array<int, 4 * batchSize>   hitIndexs;

    // -- this is the output of the fit
    std::array<float, batchSize> qps;
    std::array<float, batchSize> chi2TTs;
    std::array<float, batchSize> xTTs;
    std::array<float, batchSize> xSlopeTTs;
    std::array<float, batchSize> ys;

    // -- and this the original state (in the Velo)
    std::array<float, 3 * batchSize> statePoss;
    std::array<float, 2 * batchSize> stateDirs;
    std::array<int, batchSize>       ancestorIndexs;

    // -- and this an index to find the hit containers
    std::array<int, batchSize> hitContIndexs;

    std::size_t size{0};
    SOA_ACCESSOR_VAR( x, &( xs[pos * batchSize] ), int pos )
    SOA_ACCESSOR_VAR( z, &( zs[pos * batchSize] ), int pos )
    SOA_ACCESSOR_VAR( weight, &( weightss[pos * batchSize] ), int pos )
    SOA_ACCESSOR_VAR( sin, &( sins[pos * batchSize] ), int pos )
    SOA_ACCESSOR_VAR( id, &( ids[pos * batchSize] ), int pos )
    SOA_ACCESSOR_VAR( hitIndex, &( hitIndexs[pos * batchSize] ), int pos )

    SOA_ACCESSOR( qp, qps.data() )
    SOA_ACCESSOR( chi2TT, chi2TTs.data() )
    SOA_ACCESSOR( xTT, xTTs.data() )
    SOA_ACCESSOR( xSlopeTT, xSlopeTTs.data() )
    SOA_ACCESSOR( y, ys.data() )

    SOA_ACCESSOR( ancestorIndex, ancestorIndexs.data() )
    SOA_ACCESSOR( hitContIndex, hitContIndexs.data() )
    VEC3_SOA_ACCESSOR( pos, (float*)&( statePoss[0] ), (float*)&( statePoss[batchSize] ),
                       (float*)&( statePoss[2 * batchSize] ) )
    VEC3_XY_SOA_ACCESSOR( dir, (float*)&( stateDirs[0] ), (float*)&( stateDirs[batchSize] ), 1.0f )

    SOA_ACCESSOR( wb, wbs.data() )
    SOA_ACCESSOR( xMidField, xMidFields.data() )
    SOA_ACCESSOR( invKinkVeloDist, invKinkVeloDists.data() )

    void initTracks( int indexVal, float maxPseudoChi2 ) {
      hitIndexs.fill( indexVal );
      chi2TTs.fill( maxPseudoChi2 );
      size = 0;
    }

    template <typename dType>
    void fillHelperParams( Vec3<typename dType::float_v> pos, Vec3<typename dType::float_v> dir, const float zKink,
                           const float sigmaVeloSlope ) {

      using F = typename dType::float_v;

      F( pos.x + dir.x * ( zKink - pos.z ) ).store( xMidFields.data() );
      F a = sigmaVeloSlope * ( zKink - pos.z );
      F( 1.0f / ( a * a ) ).store( wbs.data() );
      F( 1.0f / ( zKink - pos.z ) ).store( invKinkVeloDists.data() );
    }

    template <typename dType>
    void storeSimpleFitInfo( typename dType::float_v chi2TT, typename dType::float_v qp, typename dType::float_v xTT,
                             typename dType::float_v xSlopeTT, unsigned int trackIndex ) {

      using F = typename dType::float_v;

      F( chi2TT ).store( &chi2TTs[trackIndex] );
      F( qp ).store( &qps[trackIndex] );
      F( xTT ).store( &xTTs[trackIndex] );
      F( xSlopeTT ).store( &xSlopeTTs[trackIndex] );
    }

    // -- this runs over all 4 layers, even if no hit was found
    // -- but it fills a weight of 0
    // -- Note: These are not "physical" layers, as the hits are ordered such that only
    // -- the last one can be not filled.
    template <typename dType>
    void fillHitInfo( const LHCb::Event::Zip<SIMDWrapper::Scalar, LHCb::Pr::UT::Mut::Hits>& hitsInL,
                      unsigned int                                                          trackIndex ) {

      using F = typename dType::float_v;
      using I = typename dType::int_v;
      auto nValid{0};
      if ( hitsInL.size() != 0 ) {
        for ( int i = 0; i < totalUTLayers; ++i ) {
          int        hitI   = hitIndexs[trackIndex + i * batchSize];
          const auto valid  = ( hitI != -1 );
          const auto weight = valid ? hitsInL[hitI].weight().cast() : 0.0f;
          nValid += valid;
          hitI = std::max( 0, hitI );
          LHCb::LHCbID id( LHCb::Detector::UT::ChannelID( hitsInL[hitI].channelID().cast() ) );

          F( weight ).store( &weightss[trackIndex + i * batchSize] );
          F( hitsInL[hitI].x() ).store( &xs[trackIndex + i * batchSize] );
          F( hitsInL[hitI].z() ).store( &zs[trackIndex + i * batchSize] );
          F( hitsInL[hitI].sin() ).store( &sins[trackIndex + i * batchSize] );
          I( hitsInL[hitI].index() ).store( &hitIndexs[trackIndex + i * batchSize] );
          I( id.lhcbID() ).store( &ids[trackIndex + i * batchSize] );
        }
      }
      // this is useful in filterMode
      if ( !nValid ) { F( nanMomentum ).store( &qps[trackIndex] ); }
    }
  };

  namespace {
    struct VeloUTGeomCache final {
      VeloUTGeomCache( const DeUTDetector& det, const UTMagnetTool::Cache& magtoolcache ) : common( det ) {
        // m_zMidUT is a position of normalization plane which should to be close to z middle of UT ( +- 5 cm ).
        // Cached once in VeloUTTool at initialization. No need to update with small UT movement.
        zMidUT = magtoolcache.zCenterUT;
        // zMidField and distToMomentum is properly recalculated in PrUTMagnetTool when B field changes
        distToMomentum = 1.0f / magtoolcache.dist2mom;
      }
      VeloUTGeomCache() = default;
      UTDAQ::GeomCache common;
      float            zMidUT{0.}, distToMomentum{0.};
    };
  } // namespace

  class VeloUT : public LHCb::Algorithm::Transformer<Upstream::Tracks( const EventContext&, const Velo::Tracks&,
                                                                       const LHCb::Pr::UT::Hits&,
                                                                       const VeloUTGeomCache&, const DeMagnet& ),
                                                     LHCb::DetDesc::usesConditions<VeloUTGeomCache, DeMagnet>> {

  public:
    /// Standard constructor
    VeloUT( const std::string& name, ISvcLocator* pSvcLocator );

    StatusCode initialize() override;

    Upstream::Tracks operator()( const EventContext&, const Velo::Tracks&, const LHCb::Pr::UT::Hits&,
                                 const VeloUTGeomCache&, const DeMagnet& ) const override final;

  private:
    Gaudi::Property<float> m_minMomentum{this, "MinMomentum", 1500.f * Gaudi::Units::MeV};
    Gaudi::Property<float> m_minPT{this, "MinPT", 300.f * Gaudi::Units::MeV};
    Gaudi::Property<float> m_minMomentumFinal{this, "MinMomentumFinal", 2500.f * Gaudi::Units::MeV};
    Gaudi::Property<float> m_minPTFinal{this, "MinPTFinal", 425.f * Gaudi::Units::MeV};
    Gaudi::Property<float> m_maxPseudoChi2{this, "MaxPseudoChi2", 1280.};
    Gaudi::Property<float> m_yTol{this, "YTolerance", 0.5 * Gaudi::Units::mm}; // 0.8
    Gaudi::Property<float> m_yTolSlope{this, "YTolSlope", 0.08};               // 0.2
    Gaudi::Property<float> m_hitTol{this, "HitTol", 0.8 * Gaudi::Units::mm};   // 0.8
    Gaudi::Property<float> m_deltaTx{this, "DeltaTx", 0.018};                  // 0.02
    Gaudi::Property<float> m_maxXSlope{this, "MaxXSlope", 0.350};
    Gaudi::Property<float> m_maxYSlope{this, "MaxYSlope", 0.300};
    Gaudi::Property<float> m_centralHoleSize{this, "CentralHoleSize", 33. * Gaudi::Units::mm};
    Gaudi::Property<float> m_intraLayerDist{this, "IntraLayerDist", 15.0 * Gaudi::Units::mm};
    Gaudi::Property<float> m_overlapTol{this, "OverlapTol", 0.5 * Gaudi::Units::mm};
    Gaudi::Property<float> m_passHoleSize{this, "PassHoleSize", 40. * Gaudi::Units::mm};
    Gaudi::Property<float> m_LD3Hits{this, "LD3HitsMin", -0.5};
    Gaudi::Property<float> m_LD4Hits{this, "LD4HitsMin", -0.5};
    Gaudi::Property<int>   m_minLayers{this, "MinLayers", 3};

    Gaudi::Property<bool> m_printVariables{this, "PrintVariables", false};
    Gaudi::Property<bool> m_passTracks{this, "PassTracks", false};
    Gaudi::Property<bool> m_finalFit{this, "FinalFit", true};
    Gaudi::Property<bool> m_fiducialCuts{this, "FiducialCuts", true};
    Gaudi::Property<bool> m_filterMode{this, "FilterMode", false};

    mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_seedsCounter{this, "#seeds"};
    mutable Gaudi::Accumulators::SummingCounter<unsigned int> m_tracksCounter{this, "#tracks"};

    LHCb::UT::TrackUtils::MiniStates getStates( const Velo::Tracks& inputTracks, Upstream::Tracks& outputTracks,
                                                float zMidUT ) const;

    bool getHitsScalar( const LHCb::Pr::UT::Hits& hh, const LHCb::UT::TrackUtils::MiniStates& filteredStates,
                        const std::array<LHCb::UT::TrackUtils::BoundariesNominal, 4>& compBoundsArray,
                        LHCb::Pr::UT::Mut::Hits& hitsInLayers, const std::size_t t, float zMidUT,
                        int minLayers = totalUTLayers - 1 ) const;

    inline void findHits( const LHCb::Pr::UT::Hits& hh, const simd::float_v& yProto, const simd::float_v& ty,
                          const simd::float_v& tx, const simd::float_v& xOnTrackProto, const simd::float_v& tolProto,
                          const simd::float_v& xTolNormFact, LHCb::Pr::UT::Mut::Hits& mutHits,
                          const simd::float_v& yTol, const int firstIndex, const int lastIndex ) const;

    template <bool forward>
    bool formClusters( const LHCb::Pr::UT::Mut::Hits& hitsInLayers, ProtoTracks& pTracks, const int trackIndex,
                       float zMidUT ) const;

    template <typename BdlTable>
    void prepareOutputTrackSIMD( const ProtoTracks&                                    protoTracks,
                                 const std::array<LHCb::Pr::UT::Mut::Hits, batchSize>& hitsInLayers,
                                 Upstream::Tracks& outputTracks, const Velo::Tracks& inputTracks,
                                 const BdlTable& bdlTable, const DeMagnet& magnet, float zMidUT ) const;

    /// Multipurpose tool for Bdl and deflection
    ToolHandle<UTMagnetTool> m_PrUTMagnetTool{this, "PrUTMagnetTool", "PrUTMagnetTool"};

    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_too_much_in_filtered{
        this, "Reached the maximum number of tracks in filteredStates!!"};
    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_too_much_in_boundaries{
        this, "Reached the maximum number of tracks in Boundaries!!"};

    constexpr static float c_zKink{1780.0};
    constexpr static float c_sigmaVeloSlope{0.10 * Gaudi::Units::mrad};
    constexpr static float c_invSigmaVeloSlope{10.0 / Gaudi::Units::mrad};
  };
} // namespace LHCb::Pr

//-----------------------------------------------------------------------------
// Implementation file for class : PrVeloUT
//
// 2007-05-08: Mariusz Witek
// 2017-03-01: Christoph Hasse (adapt to future framework)
// 2019-04-26: Arthur Hennequin (change data Input/Output)
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_COMPONENT_WITH_ID( LHCb::Pr::VeloUT, "PrVeloUT" )

using TracksTag  = LHCb::Pr::Upstream::Tag;
namespace HitTag = LHCb::Pr::UT::Mut::HitTag;
namespace LHCb::Pr {
  namespace {

    simd::mask_v CholeskyDecomposition3( const std::array<simd::float_v, 6>& mat, std::array<simd::float_v, 3>& rhs ) {
      // -- copied from Root::Math::CholeskyDecomp
      // -- first decompose
      std::array<simd::float_v, 6> dst;
      simd::mask_v                 mask = !( mat[0] > simd::float_v{0.0f} );
      dst[0]                            = max( 1e-6f, dst[0] ); // that's only needed if you care about FPE
      dst[0]                            = 1.0f / sqrt( mat[0] );
      dst[1]                            = mat[1] * dst[0];
      dst[2]                            = mat[2] - dst[1] * dst[1];
      mask                              = mask || ( dst[2] < simd::float_v{0.0f} );
      dst[2]                            = max( 1e-6f, dst[2] ); // that's only needed if you care about FPE
      dst[2]                            = 1.0f / sqrt( dst[2] );
      dst[3]                            = mat[3] * dst[0];
      dst[4]                            = ( mat[4] - dst[1] * dst[3] ) * dst[2];
      dst[5]                            = mat[5] - ( dst[3] * dst[3] + dst[4] * dst[4] );
      mask                              = mask || ( dst[5] < simd::float_v{0.0f} );
      dst[5]                            = max( 1e-6f, dst[5] ); // that's only needed if you care about FPE
      dst[5]                            = 1.0f / sqrt( dst[5] );

      // -- then solve
      // -- solve Ly = rhs
      const simd::float_v y0 = rhs[0] * dst[0];
      const simd::float_v y1 = ( rhs[1] - dst[1] * y0 ) * dst[2];
      const simd::float_v y2 = ( rhs[2] - ( dst[3] * y0 + dst[4] * y1 ) ) * dst[5];
      // solve L^Tx = y, and put x into rhs
      rhs[2] = (y2)*dst[5];
      rhs[1] = ( y1 - ( dst[4] * rhs[2] ) ) * dst[2];
      rhs[0] = ( y0 - ( dst[3] * rhs[2] + dst[1] * rhs[1] ) ) * dst[0];

      return mask;
    }

    // -- parameters that describe the z position of the kink point as a function of ty in a 4th order polynomial (even
    // terms only)
    constexpr auto magFieldParams = std::array{2010.0f, -2240.0f, -71330.f};

    // perform a fit using trackhelper's best hits with y correction, improve qop estimate
    simd::float_v fastfitterSIMD( std::array<simd::float_v, 4>& improvedParams, const ProtoTracks& protoTracks,
                                  const float zMidUT, const simd::float_v qpxz2p, const int t,
                                  simd::mask_v& goodFitMask ) {

      const Vec3<simd::float_v> pos = protoTracks.pos<simd::float_v>( t );
      const Vec3<simd::float_v> dir = protoTracks.dir<simd::float_v>( t );

      const simd::float_v x  = pos.x;
      const simd::float_v y  = pos.y;
      const simd::float_v z  = pos.z;
      const simd::float_v tx = dir.x;
      const simd::float_v ty = dir.y;
      const simd::float_v zKink =
          magFieldParams[0] - ty * ty * magFieldParams[1] - ty * ty * ty * ty * magFieldParams[2];
      const simd::float_v xMidField = x + tx * ( zKink - z );

      const simd::float_v zDiff = 0.001f * ( zKink - zMidUT );

      // -- This is to avoid division by zero...
      const simd::float_v pHelper = max( abs( protoTracks.qp<simd::float_v>( t ) * qpxz2p ), 1e-9f );
      const simd::float_v invP    = pHelper * 1.f / sqrt( 1.0f + ty * ty );

      // these resolution are semi-empirical, could be tuned and might not be correct for low momentum.
      const simd::float_v error1 =
          0.14f + 10000.0f * invP; // this is the resolution due to multiple scattering between Velo and UT
      const simd::float_v error2 = 0.12f + 3000.0f * invP; // this is the resolution due to the finite Velo resolution
      const simd::float_v error  = error1 * error1 + error2 * error2;
      const simd::float_v weight = 1.0f / error;

      std::array<simd::float_v, 6> mat = {weight, weight * zDiff, weight * zDiff * zDiff, 0.0f, 0.0f, 0.0f};
      std::array<simd::float_v, 3> rhs = {weight * xMidField, weight * xMidField * zDiff, 0.0f};

      for ( int i = 0; i < 4; ++i ) {
        // -- there are 3-hit candidates, but we'll
        // -- just treat them like 4-hit candidates
        // -- with 0 weight for the last hit
        const simd::float_v ui = protoTracks.x<simd::float_v>( t, i );
        const simd::float_v dz = 0.001f * ( protoTracks.z<simd::float_v>( t, i ) - zMidUT );
        const simd::float_v w  = protoTracks.weight<simd::float_v>( t, i );
        const simd::float_v ta = protoTracks.sin<simd::float_v>( t, i );
        mat[0] += w;
        mat[1] += w * dz;
        mat[2] += w * dz * dz;
        mat[3] += w * ta;
        mat[4] += w * dz * ta;
        mat[5] += w * ta * ta;
        rhs[0] += w * ui;
        rhs[1] += w * ui * dz;
        rhs[2] += w * ui * ta;
      }

      goodFitMask = !CholeskyDecomposition3( mat, rhs );

      const simd::float_v xUTFit      = rhs[0];
      const simd::float_v xSlopeUTFit = 0.001f * rhs[1];
      const simd::float_v offsetY     = rhs[2];

      const simd::float_v distX = ( xMidField - xUTFit - xSlopeUTFit * ( zKink - zMidUT ) );
      // -- This takes into account that the distance between a point and track is smaller than the distance on the
      // x-axis
      const simd::float_v distCorrectionX2 = 1.0f / ( 1 + xSlopeUTFit * xSlopeUTFit );
      simd::float_v       chi2 = weight * ( distX * distX * distCorrectionX2 + offsetY * offsetY / ( 1.0f + ty * ty ) );

      for ( int i = 0; i < 4; ++i ) {

        const simd::float_v dz   = protoTracks.z<simd::float_v>( t, i ) - zMidUT;
        const simd::float_v w    = protoTracks.weight<simd::float_v>( t, i );
        const simd::float_v dist = ( protoTracks.x<simd::float_v>( t, i ) - xUTFit - xSlopeUTFit * dz -
                                     offsetY * protoTracks.sin<simd::float_v>( t, i ) );

        chi2 += w * dist * dist * distCorrectionX2;
      }

      // new VELO slope x
      const simd::float_v xb =
          0.5f * ( ( xUTFit + xSlopeUTFit * ( zKink - zMidUT ) ) + xMidField ); // the 0.5 is empirical
      const simd::float_v xSlopeVeloFit = ( xb - x ) / ( zKink - z );

      improvedParams = {xUTFit, xSlopeUTFit, y + ty * ( zMidUT - z ) + offsetY, chi2};

      // calculate q/p
      const simd::float_v sinInX  = xSlopeVeloFit / sqrt( 1.0f + xSlopeVeloFit * xSlopeVeloFit + ty * ty );
      const simd::float_v sinOutX = xSlopeUTFit / sqrt( 1.0f + xSlopeUTFit * xSlopeUTFit + ty * ty );
      return ( sinInX - sinOutX );
    }

    // -- Evaluate the linear discriminant
    // -- Coefficients derived with LD method for p, pT and chi2 with TMVA
    template <int nHits>
    simd::float_v evaluateLinearDiscriminantSIMD( const std::array<simd::float_v, 3>& inputValues ) {

      constexpr auto coeffs =
          ( nHits == 3 ? std::array{0.162880166064f, -0.107081172665f, 0.134153123662f, -0.137764853657f}
                       : std::array{0.235010729187f, -0.0938323617311f, 0.110823681145f, -0.170467109599f} );

      assert( coeffs.size() == inputValues.size() + 1 );

      return simd::float_v{coeffs[0]} + coeffs[1] * log<simd::float_v>( inputValues[0] ) +
             coeffs[2] * log<simd::float_v>( inputValues[1] ) + coeffs[3] * log<simd::float_v>( inputValues[2] );
    }

    /*
    simd::float_v calcXTol( const simd::float_v minMom, const simd::float_v ty ) {
      return ( 38000.0f / minMom + 0.25f ) * ( 1.0f + ty * ty * 0.8f );
    }
    */
    // --------------------------------------------------------------------
    // -- Helper function to calculate the planeCode: 0 - 1 - 2 - 3
    // --------------------------------------------------------------------
    int planeCode( unsigned int id ) {
      const int station = ( (unsigned int)id & static_cast<unsigned int>( UTInfo::MasksBits::StationMask ) ) >>
                          static_cast<int>( UTInfo::MasksBits::StationBits );
      const int layer = ( (unsigned int)id & static_cast<unsigned int>( UTInfo::MasksBits::LayerMask ) ) >>
                        static_cast<int>( UTInfo::MasksBits::LayerBits );
      return 2 * ( station - 1 ) + ( layer - 1 );
    }

    // --------------------------------------------------------------------
    // -- Helper function to find duplicates in hits in the output
    // --------------------------------------------------------------------
    [[maybe_unused]] bool findDuplicates( const Upstream::Tracks& outputTracks ) {
      for ( auto const& track : outputTracks.scalar() ) {
        std::vector<LHCbID> IDs;
        IDs.reserve( track.nUTHits().cast() );
        for ( int h = 0; h < track.nUTHits().cast(); h++ ) { IDs.emplace_back( track.ut_lhcbID( h ).LHCbID() ); }

        std::sort( IDs.begin(), IDs.end() );
        if ( std::adjacent_find( IDs.begin(), IDs.end() ) != IDs.end() ) return false;
      }
      return true;
    }
    // --------------------------------------------------------------------

    // -- These things are all hardcopied from the PrTableForFunction
    // -- and PrUTMagnetTool
    // -- If the granularity or whatever changes, this will give wrong results
    simd::int_v masterIndexSIMD( const simd::int_v index1, const simd::int_v index2, const simd::int_v index3 ) {
      return ( index3 * 11 + index2 ) * 31 + index1;
    }

    constexpr auto minValsBdl = std::array{-0.3f, -250.0f, 0.0f};
    constexpr auto maxValsBdl = std::array{0.3f, 250.0f, 800.0f};
    constexpr auto deltaBdl   = std::array{0.02f, 50.0f, 80.0f};
    // constexpr auto dxDyHelper = std::array{0.0f, 1.0f, -1.0f, 0.0f};
    // ===========================================================================================
    // -- 2 helper functions for fit
    // -- Pseudo chi2 fit, templated for 3 or 4 hits
    // ===========================================================================================
    inline __attribute__( ( always_inline ) ) void
    addHit( span<float, 3> mat, span<float, 2> rhs, const LHCb::Pr::UT::Mut::Hits& hits, int index, float zMidUT ) {
      const auto& hit = hits.scalar()[index];
      const float ui  = hit.x().cast();
      const float ci  = hit.cos().cast();
      const float dz  = 0.001f * ( hit.z().cast() - zMidUT );
      const float wi  = hit.weight().cast();
      mat[0] += wi * ci;
      mat[1] += wi * ci * dz;
      mat[2] += wi * ci * dz * dz;
      rhs[0] += wi * ui;
      rhs[1] += wi * ui * dz;
    }
    template <std::size_t N>
    inline __attribute__( ( always_inline ) ) void
    simpleFit( const std::array<int, N>& indices, const LHCb::Pr::UT::Mut::Hits& hits, ProtoTracks& pTracks,
               const int trackIndex, float zMidUT, float zKink, float invSigmaVeloSlope ) {
      static_assert( N == 3 || N == 4 );

      // -- Scale the z-component, to not run into numerical problems
      // -- with floats
      const float wb              = pTracks.wb<scalar::float_v>( 0 ).cast();
      const float xMidField       = pTracks.xMidField<scalar::float_v>( 0 ).cast();
      const float invKinkVeloDist = pTracks.invKinkVeloDist<scalar::float_v>( 0 ).cast();
      const float stateX          = pTracks.pos<scalar::float_v>( trackIndex ).x.cast();
      const float stateTx         = pTracks.dir<scalar::float_v>( trackIndex ).x.cast();

      const float zDiff = 0.001f * ( zKink - zMidUT );
      auto        mat   = std::array{wb, wb * zDiff, wb * zDiff * zDiff};
      auto        rhs   = std::array{wb * xMidField, wb * xMidField * zDiff};

      auto const muthit = hits.scalar();
      std::for_each( indices.begin(), indices.end(),
                     [&]( const auto index ) { addHit( mat, rhs, hits, index, zMidUT ); } );

      ROOT::Math::CholeskyDecomp<float, 2> decomp( mat.data() );
      if ( !decomp ) return;

      decomp.Solve( rhs );

      const float xSlopeTTFit = 0.001f * rhs[1];
      const float xTTFit      = rhs[0];

      // new VELO slope x
      const float xb            = xTTFit + xSlopeTTFit * ( zKink - zMidUT );
      const float xSlopeVeloFit = ( xb - stateX ) * invKinkVeloDist;
      const float chi2VeloSlope = ( stateTx - xSlopeVeloFit ) * invSigmaVeloSlope;

      const float chi2TT = std::accumulate( indices.begin(), indices.end(), chi2VeloSlope * chi2VeloSlope,
                                            [&]( float chi2, const int index ) {
                                              const float du =
                                                  ( xTTFit + xSlopeTTFit * ( muthit[index].z().cast() - zMidUT ) ) -
                                                  muthit[index].x().cast();
                                              return chi2 + muthit[index].weight().cast() * ( du * du );
                                            } ) /
                           ( N + 1 - 2 );

      if ( chi2TT < pTracks.chi2TT<scalar::float_v>( trackIndex ).cast() ) {

        // calculate q/p
        const float sinInX  = xSlopeVeloFit * vdt::fast_isqrtf( 1.0f + xSlopeVeloFit * xSlopeVeloFit );
        const float sinOutX = xSlopeTTFit * vdt::fast_isqrtf( 1.0f + xSlopeTTFit * xSlopeTTFit );
        const float qp      = ( sinInX - sinOutX );

        pTracks.storeSimpleFitInfo<scalar>( chi2TT, qp, xTTFit, xSlopeTTFit, trackIndex );
        for ( std::size_t i = 0; i < N; i++ ) { pTracks.store_hitIndex<scalar::int_v>( trackIndex, i, indices[i] ); }
        if constexpr ( N == 3 ) { pTracks.store_hitIndex<scalar::int_v>( trackIndex, 3, -1 ); }
      }
    }
  } // namespace

  //=============================================================================
  // Standard constructor, initializes variables
  //=============================================================================
  VeloUT::VeloUT( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"InputTracksName", "Rec/Track/Velo"}, KeyValue{"UTHits", UTInfo::HitLocation},
                      KeyValue{"GeometryInfo", "AlgorithmSpecific-" + name + "-UTGeometryInfo"},
                      KeyValue{"Magnet", LHCb::Det::Magnet::det_path}},
                     KeyValue{"OutputTracksName", "Rec/Track/UT"} ) {}

  /// Initialization
  StatusCode VeloUT::initialize() {
    return Transformer::initialize().andThen( [&] {
      addConditionDerivation( {DeUTDetLocation::location(), m_PrUTMagnetTool->cacheLocation()},
                              inputLocation<VeloUTGeomCache>(),
                              []( const DeUTDetector& det, const UTMagnetTool::Cache& magtoolcache ) {
                                return VeloUTGeomCache( det, magtoolcache );
                              } );
    } );
  }

  //=============================================================================
  // Main execution
  //=============================================================================
  Upstream::Tracks VeloUT::operator()( const EventContext& evtCtx, const Velo::Tracks& inputTracks,
                                       const LHCb::Pr::UT::Hits& hh, const VeloUTGeomCache& velogeom,
                                       const DeMagnet& magnet ) const {

    Upstream::Tracks outputTracks{&inputTracks, Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx )};
    outputTracks.reserve( inputTracks.size() );
    m_seedsCounter += inputTracks.size();

    const auto& geometry = velogeom.common;
    const auto& bdlTable = m_PrUTMagnetTool->BdlTable();

    LHCb::UT::TrackUtils::MiniStates filteredStates = getStates( inputTracks, outputTracks, velogeom.zMidUT );

    simd::float_v invMinPt       = 1.0f / m_minPT.value();
    simd::float_v invMinMomentum = 1.0f / m_minMomentum.value();

    // -- Used for the calculation of the size of the search windows
    constexpr const std::array<float, totalUTLayers> normFact{0.95f, 1.0f, 1.36f, 1.41f};

    auto compBoundsArray = LHCb::UTDAQ::findAllSectorsExtrap<
        LHCb::UT::TrackUtils::BoundariesNominal, LHCb::UT::TrackUtils::BoundariesNominalTag::types,
        LHCb::UT::TrackUtils::maxNumRowsBoundariesNominal, LHCb::UT::TrackUtils::maxNumColsBoundariesNominal>(
        filteredStates, geometry,
        [&]( int layerIndex, simd::float_v x, simd::float_v y, simd::float_v z, simd::float_v tx, simd::float_v ty,
             simd::float_v ) {
          // -- this 0.002 seems a little odd...
          const simd::float_v theta     = max( 0.002f, sqrt( tx * tx + ty * ty ) );
          const simd::float_v invMinMom = min( invMinPt * theta, invMinMomentum );

          const simd::float_v xTol =
              abs( velogeom.distToMomentum * invMinMom * normFact[layerIndex] ) - abs( tx ) * m_intraLayerDist.value();
          const simd::float_v yTol = m_yTol.value() + m_yTolSlope.value() * xTol;

          const simd::float_v zGeo{geometry.layers[layerIndex].z};
          const simd::float_v dxDy{geometry.layers[layerIndex].dxDy};

          const simd::float_v yAtZ   = y + ty * ( zGeo - z );
          const simd::float_v xLayer = x + tx * ( zGeo - z );
          const simd::float_v yLayer = yAtZ + yTol * dxDy;
          return std::make_tuple( xLayer, yLayer, xTol, yTol );
        },
        m_minLayers );

    auto hitsInLayers = LHCb::make_object_array<LHCb::Pr::UT::Mut::Hits, batchSize>( Zipping::generateZipIdentifier(),
                                                                                     LHCb::getMemResource( evtCtx ) );
    for ( auto& hits : hitsInLayers ) { hits.reserve( 32 ); }
    ProtoTracks pTracks;

    // -- We cannot put all found hits in an array, as otherwise the stack overflows
    // -- so we just do the whole thing in batches
    const std::size_t filteredStatesSize = filteredStates.size();
    for ( std::size_t t = 0; t < filteredStatesSize; t += batchSize ) {

      // -- This is scalar, as the hits are found in a scalar way
      filteredStates.clear();
      for ( std::size_t t2 = 0; t2 < batchSize && t2 + t < filteredStatesSize; ++t2 ) {
        const auto fSize = filteredStates.size();
        hitsInLayers[fSize].clear();
        hitsInLayers[fSize].layerIndices.fill( -1 );
        const bool foundHits = getHitsScalar( hh, filteredStates, compBoundsArray, hitsInLayers[fSize], t + t2,
                                              velogeom.zMidUT, m_minLayers );
        filteredStates.copy_back<scalar>( filteredStates, t + t2, foundHits );
      }

      pTracks.initTracks( -1, m_maxPseudoChi2.value() );
      for ( const auto& fState : filteredStates.scalar() ) {
        const auto tEff = fState.offset();

        Vec3<scalar::float_v> pos{fState.get<LHCb::UT::TrackUtils::MiniStateTag::State>().x(),
                                  fState.get<LHCb::UT::TrackUtils::MiniStateTag::State>().y(),
                                  fState.get<LHCb::UT::TrackUtils::MiniStateTag::State>().z()};
        Vec3<scalar::float_v> dir{fState.get<LHCb::UT::TrackUtils::MiniStateTag::State>().tx(),
                                  fState.get<LHCb::UT::TrackUtils::MiniStateTag::State>().ty(), 1.f};

        const int trackIndex = pTracks.size;
        pTracks.fillHelperParams<scalar>( pos, dir, c_zKink, c_sigmaVeloSlope );
        pTracks.store_pos<scalar::float_v>( trackIndex, pos );
        pTracks.store_dir<scalar::float_v>( trackIndex, dir );

        if ( !formClusters<true>( hitsInLayers[tEff], pTracks, trackIndex, velogeom.zMidUT ) ) {
          formClusters<false>( hitsInLayers[tEff], pTracks, trackIndex, velogeom.zMidUT );
        }

        if ( !m_filterMode && pTracks.hitIndex<scalar::int_v>( trackIndex, 0 ).cast() == -1 ) continue;

        scalar::int_v ancestorIndex = fState.get<LHCb::UT::TrackUtils::MiniStateTag::index>();
        pTracks.store_ancestorIndex<scalar::int_v>( trackIndex, ancestorIndex );
        pTracks.store_hitContIndex<scalar::int_v>( trackIndex, tEff );
        // -- this runs over all 4 layers, even if no hit was found
        // -- but it fills a weight of 0
        // -- Note: These are not "physical" layers, as the hits are ordered such that only
        // -- the last one can be not filled.
        auto hitsInL = hitsInLayers[tEff].scalar();
        pTracks.fillHitInfo<scalar>( hitsInL, trackIndex );
        pTracks.size++;
      }

      // padding to avoid FPEs
      if ( ( pTracks.size + simd::size ) < batchSize ) {
        pTracks.store_pos<simd::float_v>( pTracks.size, Vec3<simd::float_v>( 1.f, 1.f, 1.f ) );
        pTracks.store_dir<simd::float_v>( pTracks.size, Vec3<simd::float_v>( 1.f, 1.f, 1.f ) );
      }
      prepareOutputTrackSIMD( pTracks, hitsInLayers, outputTracks, inputTracks, bdlTable, magnet, velogeom.zMidUT );
    }

    // -- The algorithm should not store duplicated hits...
    assert( findDuplicates( outputTracks ) && "Hit duplicates found" );

    m_tracksCounter += outputTracks.size();
    return outputTracks;
  }
  //=============================================================================
  // Get the state, do some cuts
  //=============================================================================
  __attribute__( ( flatten ) ) LHCb::UT::TrackUtils::MiniStates
  VeloUT::getStates( const Velo::Tracks& inputTracks, Upstream::Tracks& outputTracks, float zMidUT ) const {

    const int                        EndVelo = 1;
    LHCb::UT::TrackUtils::MiniStates filteredStates{Zipping::generateZipIdentifier(),
                                                    {inputTracks.get_allocator().resource()}};
    filteredStates.reserve( inputTracks.size() );

    const auto centralHoleR2 = simd::float_v{m_centralHoleSize * m_centralHoleSize};

    for ( auto const& velotrack : inputTracks.simd() ) {
      auto const loopMask = velotrack.loop_mask();
      auto const trackVP  = velotrack.indices();
      auto       pos      = velotrack.StatePos( EndVelo );
      auto       dir      = velotrack.StateDir( EndVelo );
      auto       covX     = velotrack.StateCovX( EndVelo );

      simd::float_v xMidUT = pos.x + dir.x * ( zMidUT - pos.z );
      simd::float_v yMidUT = pos.y + dir.y * ( zMidUT - pos.z );

      simd::mask_v centralHoleMask = xMidUT * xMidUT + yMidUT * yMidUT < centralHoleR2;
      simd::mask_v slopesMask   = ( ( abs( dir.x ) > m_maxXSlope.value() ) || ( abs( dir.y ) > m_maxYSlope.value() ) );
      simd::mask_v passHoleMask = abs( xMidUT ) < m_passHoleSize.value() && abs( yMidUT ) < m_passHoleSize.value();
      simd::mask_v mask         = centralHoleMask || slopesMask;
      simd::mask_v csMask       = loopMask && !mask && ( !simd::mask_v{m_passTracks.value()} || !passHoleMask );

      auto fState = filteredStates.compress_back<SIMDWrapper::InstructionSet::Best>( csMask );
      fState.field<LHCb::UT::TrackUtils::MiniStateTag::State>().setPosition( pos.x, pos.y, pos.z );
      fState.field<LHCb::UT::TrackUtils::MiniStateTag::State>().setDirection( dir.x, dir.y );
      fState.field<LHCb::UT::TrackUtils::MiniStateTag::State>().setQOverP( 0.f );
      fState.field<LHCb::UT::TrackUtils::MiniStateTag::index>().set( trackVP );

      if ( m_passTracks ) {
        auto outMask = loopMask && passHoleMask; // not sure if correct...

        auto const track = outputTracks.compress_back( outMask );
        track.field<TracksTag::trackVP>().set( velotrack.indices() );
        track.field<TracksTag::State>().setQOverP( 0.f );
        track.field<TracksTag::State>().setPosition( pos );
        track.field<TracksTag::State>().setDirection( dir );
        track.field<TracksTag::UTHits>().resize( 0 );
        track.field<TracksTag::VPHits>().resize( velotrack.nHits() );
        track.field<TracksTag::StateCovX>( 0 ).set( covX.x );
        track.field<TracksTag::StateCovX>( 1 ).set( covX.y );
        track.field<TracksTag::StateCovX>( 2 ).set( covX.z );
        for ( int idx = 0; idx < velotrack.nHits().hmax( outMask ); ++idx ) {
          track.field<TracksTag::VPHits>()[idx].template field<TracksTag::Index>().set( velotrack.vp_index( idx ) );
          track.field<TracksTag::VPHits>()[idx].template field<TracksTag::LHCbID>().set( velotrack.vp_lhcbID( idx ) );
        }
      }
    }
    return filteredStates;
  }

  //=============================================================================
  // Find the hits
  //=============================================================================
  inline __attribute__( ( always_inline ) ) bool
  VeloUT::getHitsScalar( const LHCb::Pr::UT::Hits& hh, const LHCb::UT::TrackUtils::MiniStates& filteredStates,
                         const std::array<LHCb::UT::TrackUtils::BoundariesNominal, totalUTLayers>& compBoundsArray,
                         LHCb::Pr::UT::Mut::Hits& hitsInLayers, const std::size_t t, float zMidUT,
                         int minLayers ) const {

    // -- This is for some sanity checks later
    constexpr const int maxSectorsPerRegion = static_cast<int>( UTInfo::SectorNumbers::MaxSectorsPerRegion );
    constexpr const int maxLayer            = totalUTLayers;
    constexpr const int maxRegion           = static_cast<int>( UTInfo::DetectorNumbers::Regions );
    [[maybe_unused]] constexpr const int maxSectorNumber =
        maxSectorsPerRegion + ( ( maxLayer - 1 ) * maxRegion + ( maxRegion - 1 ) ) * maxSectorsPerRegion;

    const simd::float_v tolProto{m_yTol.value()};

    const auto fState  = filteredStates.scalar();
    const auto xState  = fState[t].get<LHCb::UT::TrackUtils::MiniStateTag::State>().x().cast();
    const auto yState  = fState[t].get<LHCb::UT::TrackUtils::MiniStateTag::State>().y().cast();
    const auto zState  = fState[t].get<LHCb::UT::TrackUtils::MiniStateTag::State>().z().cast();
    const auto txState = fState[t].get<LHCb::UT::TrackUtils::MiniStateTag::State>().tx().cast();
    const auto tyState = fState[t].get<LHCb::UT::TrackUtils::MiniStateTag::State>().ty().cast();

    // in filter mode tracks close to the hole in the centre of the UT may have no hits
    if ( m_filterMode ) {
      const auto xMidUT = xState + txState * ( zMidUT - zState );
      const auto yMidUT = yState + tyState * ( zMidUT - zState );
      const auto rMidUT = xMidUT * xMidUT + yMidUT * yMidUT;
      minLayers         = rMidUT < m_passHoleSize * m_passHoleSize ? 0 : minLayers;
    }

    std::size_t nSize   = 0;
    int         nLayers = 0;

    // -- the protos could be precomputed
    const simd::float_v yProto{yState - tyState * zState};
    const simd::float_v xOnTrackProto{xState - txState * zState};
    const simd::float_v ty{tyState};
    const simd::float_v tx{txState};

    // -- the second condition is to ensure at least 3 layers with hits
    for ( int layerIndex = 0; layerIndex < totalUTLayers && layerIndex - nLayers <= totalUTLayers - minLayers;
          ++layerIndex ) {

      const auto          compBoundsArr = compBoundsArray[layerIndex].scalar();
      const auto          xTolS = compBoundsArr[t].get<LHCb::UT::TrackUtils::BoundariesNominalTag::xTol>().cast();
      const auto          nPos  = compBoundsArr[t].get<LHCb::UT::TrackUtils::BoundariesNominalTag::nPos>().cast();
      const simd::float_v yTol  = m_yTol.value() + m_yTolSlope.value() * xTolS;
      const simd::float_v xTol  = xTolS + abs( tx ) * m_intraLayerDist.value();

      assert( nPos < maxNumSectors && "nPos out of bound" );

      for ( int j = 0; j < nPos; j++ ) {

        const int sectA = compBoundsArr[t].get<LHCb::UT::TrackUtils::BoundariesNominalTag::sects>( j ).cast();
        const int sectB = ( j == nPos - 1 )
                              ? sectA
                              : compBoundsArr[t].get<LHCb::UT::TrackUtils::BoundariesNominalTag::sects>( j + 1 ).cast();

        assert( sectA != LHCb::UTDAQ::paddingSectorNumber && "sectA points to padding element" );
        assert( sectB != LHCb::UTDAQ::paddingSectorNumber && "sectB points to padding element" );
        assert( ( sectA > -1 ) && ( sectA < maxSectorNumber ) && "sector number out of bound" );
        assert( ( sectB > -1 ) && ( sectB < maxSectorNumber ) && "sector number out of bound" );

        // -- Sector is allowed to be a duplicate if it is the last element (as it has no consequence then)
        assert( ( ( sectA != sectB ) || ( j == nPos - 1 ) ) && "duplicated sectors" );

        // -- The idea is to merge adjacent ranges of indices, such that collecting hits is more efficient
        // -- let's try to make it branchless
        const std::pair<int, int>& temp       = hh.indices( sectA );
        const std::pair<int, int>& temp2      = hh.indices( sectB );
        const int                  firstIndex = temp.first;
        // -- We put the lastIndex to the end of the next container if they join up
        // -- Note that this is _not_ fulfilled if the sector has elements and is a duplicate,
        // -- but this only happens if it is the padding element, in which case we are already at the last
        // -- element of the loop. If it happens in any other case, the asssert fires.
        const auto shift     = ( temp2.first == temp.second );
        const auto lastIndex = shift ? temp2.second : temp.second;
        j += shift;

        findHits( hh, yProto, ty, tx, xOnTrackProto, tolProto, xTol, hitsInLayers, yTol, firstIndex, lastIndex );
      }

      nLayers += ( nSize != hitsInLayers.size() );

      hitsInLayers.layerIndices[layerIndex] = nSize;
      nSize                                 = hitsInLayers.size();
    }

    // -- only use these hits, if we have at least 3 layers
    return nLayers >= minLayers;
  }
  // ==============================================================================
  // -- Method that finds the hits in a given layer within a certain range
  // ==============================================================================
  inline __attribute__( ( always_inline ) ) void
  VeloUT::findHits( const LHCb::Pr::UT::Hits& hh, const simd::float_v& yProto, const simd::float_v& ty,
                    const simd::float_v& tx, const simd::float_v& xOnTrackProto, const simd::float_v& tolProto,
                    const simd::float_v& xTolNormFact, LHCb::Pr::UT::Mut::Hits& mutHits, const simd::float_v& yTol,
                    const int firstIndex, const int lastIndex ) const {

    const auto& myHits = hh;
    const auto  myHs   = myHits.simd();

    for ( int i = firstIndex; i < lastIndex; i += simd::size ) {

      const auto mH = myHs[i];

      // -- Calculate distance between straight line extrapolation from Velo and hit position
      const simd::float_v yy = yProto + ty * mH.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>();
      const simd::float_v xx =
          mH.get<LHCb::Pr::UT::UTHitsTag::xAtYEq0>() + yy * mH.get<LHCb::Pr::UT::UTHitsTag::dxDy>();
      const simd::float_v xOnTrack = xOnTrackProto + tx * mH.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>();
      const simd::float_v absdx    = abs( xx - xOnTrack );

      if ( none( absdx < xTolNormFact ) ) continue;
      auto loopMask = simd::loop_mask( i, lastIndex );

      // is there anything like minmax?
      const simd::float_v yMin =
          min( mH.get<LHCb::Pr::UT::UTHitsTag::yBegin>(), mH.get<LHCb::Pr::UT::UTHitsTag::yEnd>() );
      const simd::float_v yMax =
          max( mH.get<LHCb::Pr::UT::UTHitsTag::yBegin>(), mH.get<LHCb::Pr::UT::UTHitsTag::yEnd>() );

      const simd::float_v tol  = yTol + absdx * tolProto;
      auto                mask = ( yMin - tol < yy && yy < yMax + tol ) && ( absdx < xTolNormFact ) && loopMask;

      if ( none( mask ) ) continue;
      auto muthit = mutHits.compress_back<SIMDWrapper::InstructionSet::Best>( mask );
      muthit.field<HitTag::xs>().set( xx );
      muthit.field<HitTag::zs>().set( mH.get<LHCb::Pr::UT::UTHitsTag::zAtYEq0>() );
      muthit.field<HitTag::coss>().set( mH.get<LHCb::Pr::UT::UTHitsTag::cos>() );
      muthit.field<HitTag::sins>().set( mH.get<LHCb::Pr::UT::UTHitsTag::cos>() * -1.0f *
                                        mH.get<LHCb::Pr::UT::UTHitsTag::dxDy>() );
      muthit.field<HitTag::weights>().set( mH.get<LHCb::Pr::UT::UTHitsTag::weight>() );
      muthit.field<HitTag::channelIDs>().set( mH.get<LHCb::Pr::UT::UTHitsTag::channelID>() );
      muthit.field<HitTag::indexs>().set( simd::indices( i ) ); // fill the index in the original hit container
    }
  }
  //=========================================================================
  // Form clusters
  //=========================================================================
  template <bool forward>
  inline __attribute__( ( always_inline ) ) bool VeloUT::formClusters( const LHCb::Pr::UT::Mut::Hits& hitsInLayers,
                                                                       ProtoTracks& pTracks, const int trackIndex,
                                                                       float zMidUT ) const {

    const int begin0 = forward ? hitsInLayers.layerIndices[0] : hitsInLayers.layerIndices[3];
    const int end0   = forward ? hitsInLayers.layerIndices[1] : hitsInLayers.size();

    const int begin1 = forward ? hitsInLayers.layerIndices[1] : hitsInLayers.layerIndices[2];
    const int end1   = forward ? hitsInLayers.layerIndices[2] : hitsInLayers.layerIndices[3];

    const int begin2 = forward ? hitsInLayers.layerIndices[2] : hitsInLayers.layerIndices[1];
    const int end2   = forward ? hitsInLayers.layerIndices[3] : hitsInLayers.layerIndices[2];

    const int begin3            = forward ? hitsInLayers.layerIndices[3] : hitsInLayers.layerIndices[0];
    const int end3              = forward ? hitsInLayers.size() : hitsInLayers.layerIndices[1];
    bool      fourLayerSolution = false;

    const float stateTx = pTracks.dir<scalar::float_v>( trackIndex ).x.cast();
    const auto& hitsInL = hitsInLayers.scalar();

    // -- this is scalar for the moment
    for ( int i0 = begin0; i0 < end0; ++i0 ) {

      const float xhitLayer0 = hitsInL[i0].x().cast();
      const float zhitLayer0 = hitsInL[i0].z().cast();

      // Loop over Second Layer
      for ( int i2 = begin2; i2 < end2; ++i2 ) {

        const float xhitLayer2 = hitsInL[i2].x().cast();
        const float zhitLayer2 = hitsInL[i2].z().cast();

        const float tx = ( xhitLayer2 - xhitLayer0 ) / ( zhitLayer2 - zhitLayer0 );

        if ( std::abs( tx - stateTx ) > m_deltaTx ) continue;

        int   bestHit1Index = -1;
        float hitTol        = m_hitTol;

        for ( int i1 = begin1; i1 < end1; ++i1 ) {

          const float xhitLayer1 = hitsInL[i1].x().cast();
          const float zhitLayer1 = hitsInL[i1].z().cast();

          const float xextrapLayer1 = xhitLayer0 + tx * ( zhitLayer1 - zhitLayer0 );
          if ( std::abs( xhitLayer1 - xextrapLayer1 ) < hitTol ) {
            hitTol        = std::abs( xhitLayer1 - xextrapLayer1 );
            bestHit1Index = i1;
          }
        }

        if ( fourLayerSolution && bestHit1Index == -1 ) continue;

        int bestHit3Index = -1;
        hitTol            = m_hitTol;
        for ( int i3 = begin3; i3 < end3; ++i3 ) {

          const float xhitLayer3 = hitsInL[i3].x().cast();
          const float zhitLayer3 = hitsInL[i3].z().cast();

          const float xextrapLayer3 = xhitLayer2 + tx * ( zhitLayer3 - zhitLayer2 );

          if ( std::abs( xhitLayer3 - xextrapLayer3 ) < hitTol ) {
            hitTol        = std::abs( xhitLayer3 - xextrapLayer3 );
            bestHit3Index = i3;
          }
        }
        // -- All hits found
        if ( bestHit1Index != -1 && bestHit3Index != -1 ) {
          simpleFit( std::array{i0, bestHit1Index, i2, bestHit3Index}, hitsInLayers, pTracks, trackIndex, zMidUT,
                     c_zKink, c_invSigmaVeloSlope );

          if ( !fourLayerSolution && pTracks.hitIndex<scalar::int_v>( trackIndex, 0 ).cast() != -1 ) {
            fourLayerSolution = true;
          }
          continue;
        }

        // -- Nothing found in layer 3
        if ( !fourLayerSolution && bestHit1Index != -1 ) {
          simpleFit( std::array{i0, bestHit1Index, i2}, hitsInLayers, pTracks, trackIndex, zMidUT, c_zKink,
                     c_invSigmaVeloSlope );
          continue;
        }
        // -- Noting found in layer 1
        if ( !fourLayerSolution && bestHit3Index != -1 ) {
          simpleFit( std::array{i0, bestHit3Index, i2}, hitsInLayers, pTracks, trackIndex, zMidUT, c_zKink,
                     c_invSigmaVeloSlope );
          continue;
        }
      }
    }
    return fourLayerSolution;
  }
  //=========================================================================
  // Create the Velo-UT tracks
  //=========================================================================
  template <typename BdlTable>
  inline __attribute__( ( always_inline ) ) __attribute__( ( flatten ) ) void
  VeloUT::prepareOutputTrackSIMD( const ProtoTracks&                                    protoTracks,
                                  const std::array<LHCb::Pr::UT::Mut::Hits, batchSize>& hitsInLayers,
                                  Upstream::Tracks& outputTracks, const Velo::Tracks& inputTracks,
                                  const BdlTable& bdlTable, const DeMagnet& magnet, float zMidUT ) const {

    auto const velozipped = inputTracks.simd();
    const auto pSize      = protoTracks.size;
    for ( std::size_t t = 0; t < pSize; t += simd::size ) {
      //== Handle states. copy Velo one, add TT.
      const simd::float_v zOrigin =
          select( protoTracks.dir<simd::float_v>( t ).y > 0.001f,
                  protoTracks.pos<simd::float_v>( t ).z -
                      protoTracks.pos<simd::float_v>( t ).y / protoTracks.dir<simd::float_v>( t ).y,
                  protoTracks.pos<simd::float_v>( t ).z -
                      protoTracks.pos<simd::float_v>( t ).x / protoTracks.dir<simd::float_v>( t ).x );

      // -- this is to filter tracks where the fit had a too large chi2
      simd::mask_v fourHitTrack = protoTracks.weight<simd::float_v>( t, 3 ) > 0.0001f;

      // const float bdl1    = m_PrUTMagnetTool->bdlIntegral(helper.state.ty,zOrigin,helper.state.z);

      // -- These are calculations, copied and simplified from PrTableForFunction
      // -- FIXME: these rely on the internal details of PrTableForFunction!!!
      //           and should at least be put back in there, and used from here
      //           to make sure everything _stays_ consistent...
      auto var = std::array{protoTracks.dir<simd::float_v>( t ).y, zOrigin, protoTracks.pos<simd::float_v>( t ).z};

      simd::int_v index1 = min( max( simd::int_v{( var[0] + 0.3f ) / 0.6f * 30}, 0 ), 30 );
      simd::int_v index2 = min( max( simd::int_v{( var[1] + 250 ) / 500 * 10}, 0 ), 10 );
      simd::int_v index3 = min( max( simd::int_v{var[2] / 800 * 10}, 0 ), 10 );

      simd::float_v bdl = gather( bdlTable.table().data(), masterIndexSIMD( index1, index2, index3 ) );

      // -- TODO: check if we can go outside this table...
      const std::array<simd::float_v, 3> bdls =
          std::array{gather( bdlTable.table().data(), masterIndexSIMD( index1 + 1, index2, index3 ) ),
                     gather( bdlTable.table().data(), masterIndexSIMD( index1, index2 + 1, index3 ) ),
                     gather( bdlTable.table().data(), masterIndexSIMD( index1, index2, index3 + 1 ) )};

      const std::array<simd::float_v, 3> boundaries = {-0.3f + simd::float_v{index1} * deltaBdl[0],
                                                       -250.0f + simd::float_v{index2} * deltaBdl[1],
                                                       0.0f + simd::float_v{index3} * deltaBdl[2]};

      // -- This is an interpolation, to get a bit more precision
      simd::float_v addBdlVal{0.0f};
      for ( int i = 0; i < 3; ++i ) {

        // -- this should make sure that values outside the range add nothing to the sum
        var[i] = select( minValsBdl[i] > var[i], boundaries[i], var[i] );
        var[i] = select( maxValsBdl[i] < var[i], boundaries[i], var[i] );

        const simd::float_v dTab_dVar = ( bdls[i] - bdl ) / deltaBdl[i];
        const simd::float_v dVar      = ( var[i] - boundaries[i] );
        addBdlVal += dTab_dVar * dVar;
      }
      bdl += addBdlVal;
      // ----

      // -- order is: x, tx, y, chi2
      std::array<simd::float_v, 4> finalParams = {
          protoTracks.xTT<simd::float_v>( t ), protoTracks.xSlopeTT<simd::float_v>( t ),
          protoTracks.pos<simd::float_v>( t ).y +
              protoTracks.dir<simd::float_v>( t ).y * ( zMidUT - protoTracks.pos<simd::float_v>( t ).z ),
          protoTracks.chi2TT<simd::float_v>( t )};

      const simd::float_v qpxz2p  = -1.0f / bdl * 3.3356f / Gaudi::Units::GeV;
      simd::mask_v        fitMask = simd::mask_true();
      simd::float_v       qp      = m_finalFit ? fastfitterSIMD( finalParams, protoTracks, zMidUT, qpxz2p, t, fitMask )
                                    : protoTracks.qp<simd::float_v>( t ) /
                                          sqrt( 1.0f + protoTracks.dir<simd::float_v>( t ).y *
                                                           protoTracks.dir<simd::float_v>( t ).y ); // is this correct?

      qp                      = select( fitMask, qp, protoTracks.qp<simd::float_v>( t ) );
      const simd::float_v qop = select( abs( bdl ) < 1.e-8f, simd::float_v{1000.0f}, qp * qpxz2p );

      // -- Don't make tracks that have grossly too low momentum
      // -- Beware of the momentum resolution!
      const simd::float_v p = abs( 1.0f / qop );
      const simd::float_v pt =
          p * sqrt( protoTracks.dir<simd::float_v>( t ).x * protoTracks.dir<simd::float_v>( t ).x +
                    protoTracks.dir<simd::float_v>( t ).y * protoTracks.dir<simd::float_v>( t ).y );

      const simd::float_v xUT  = finalParams[0];
      const simd::float_v txUT = finalParams[1];
      const simd::float_v yUT  = finalParams[2];

      // -- apply some fiducial cuts
      // -- they are optimised for high pT tracks (> 500 MeV)
      simd::mask_v fiducialMask = simd::mask_false();

      if ( m_fiducialCuts ) {
        const float magSign = magnet.signedRelativeCurrent();

        fiducialMask = ( magSign * qop < 0.0f && xUT > -48.0f && xUT < 0.0f && abs( yUT ) < 33.0f );
        fiducialMask = fiducialMask || ( magSign * qop > 0.0f && xUT < 48.0f && xUT > 0.0f && abs( yUT ) < 33.0f );

        fiducialMask = fiducialMask || ( magSign * qop < 0.0f && txUT > 0.09f + 0.0003f * pt );
        fiducialMask = fiducialMask || ( magSign * qop > 0.0f && txUT < -0.09f - 0.0003f * pt );
      }

      // -- evaluate the linear discriminant and reject ghosts
      // -- the values only make sense if the final fit is performed
      simd::mask_v mvaMask = simd::mask_true();

      if ( m_finalFit ) {

        const simd::float_v fourHitDisc  = evaluateLinearDiscriminantSIMD<4>( {p, pt, finalParams[3]} );
        const simd::float_v threeHitDisc = evaluateLinearDiscriminantSIMD<3>( {p, pt, finalParams[3]} );

        simd::mask_v fourHitMask  = fourHitDisc > m_LD4Hits.value();
        simd::mask_v threeHitMask = threeHitDisc > m_LD3Hits.value();

        // -- only have 3 or 4 hit tracks
        mvaMask = ( fourHitTrack && fourHitMask ) || ( !fourHitTrack && threeHitMask );
      }

      const auto pPTMask        = p > m_minMomentumFinal.value() && pt > m_minPTFinal.value();
      const auto loopMask       = simd::loop_mask( t, pSize );
      const auto validTrackMask = pPTMask && !fiducialMask && mvaMask && loopMask;

      const auto finalMask     = m_filterMode ? loopMask : validTrackMask;
      const auto finalQoP      = select( validTrackMask, qop, nanMomentum );
      const auto ancestor      = protoTracks.ancestorIndex<simd::int_v>( t );
      const auto velo_ancestor = velozipped.gather( ancestor, finalMask );
      const auto EndVelo       = 1;
      const auto currentsize   = outputTracks.size();

      const auto oTrack = outputTracks.compress_back( finalMask );
      oTrack.field<TracksTag::trackVP>().set( ancestor );
      oTrack.field<TracksTag::State>().setQOverP( finalQoP );
      oTrack.field<TracksTag::State>().setPosition( velo_ancestor.StatePos( EndVelo ) );
      oTrack.field<TracksTag::State>().setDirection( velo_ancestor.StateDir( EndVelo ) );
      oTrack.field<TracksTag::VPHits>().resize( velo_ancestor.nHits() );
      auto covX = velo_ancestor.StateCovX( EndVelo );
      oTrack.field<TracksTag::StateCovX>( 0 ).set( covX.x );
      oTrack.field<TracksTag::StateCovX>( 1 ).set( covX.y );
      oTrack.field<TracksTag::StateCovX>( 2 ).set( covX.z );
      for ( int idx = 0; idx < velo_ancestor.nHits().hmax( finalMask ); ++idx ) {
        oTrack.field<TracksTag::VPHits>()[idx].field<TracksTag::Index>().set( velo_ancestor.vp_index( idx ) );
        oTrack.field<TracksTag::VPHits>()[idx].field<TracksTag::LHCbID>().set( velo_ancestor.vp_lhcbID( idx ) );
      }

      if ( m_filterMode ) {
        oTrack.field<TracksTag::UTHits>().resize( 0 );
        continue;
      }

      const auto txArray = SIMDWrapper::to_array( txUT );
      const auto xArray  = SIMDWrapper::to_array( xUT );

      std::array<int, simd::size> nUTHits{};
      // -- This is needed to find the planeCode of the layer with the missing hit
      std::array<int, simd::size> sumLayArray{};

      // -- from here on, go over each track individually to find and add the overlap hits
      // -- this is not particularly elegant...
      // -- As before, these are "pseudo layers", i.e. it is not guaranteed that if i > j, z[i] > z[j]

      for ( int iLayer = 0; iLayer < totalUTLayers; ++iLayer ) {
        int trackIndex2 = 0;
        for ( unsigned int t2 = 0; t2 < simd::size; ++t2 ) {
          if ( !testbit( finalMask, t2 ) ) continue;
          const auto tscalar = t + t2;

          const bool goodHit = ( protoTracks.weight<scalar::float_v>( tscalar, iLayer ).cast() > 0.0001f );
          const auto hitIdx  = protoTracks.hitIndex<scalar::int_v>( tscalar, iLayer );
          const auto id      = protoTracks.id<scalar::int_v>( tscalar, iLayer );

          // -- Only add the hit, if it is not in an empty layer (that sounds like a tautology,
          // -- but given that one always has 4 hits, even if only 3 make sense, it is needed)
          // -- Only the last pseudo-layer can be an empty layer
          if ( goodHit ) {
            auto hits = outputTracks.scalar()[currentsize + trackIndex2].field<TracksTag::UTHits>();
            hits.resize( nUTHits[t2] + 1 );
            hits[nUTHits[t2]].template field<TracksTag::Index>().set( hitIdx );
            hits[nUTHits[t2]].template field<TracksTag::LHCbID>().set( LHCb::Event::lhcbid_v<scalar>( id ) );
            nUTHits[t2] += 1;
          }
          // --
          // -----------------------------------------------------------------------------------
          // -- The idea of the following code is: In layers where we have found a hit, we search for
          // -- overlap hits.
          // -- In layers where no hit was found initially, we use the better parametrization of the final
          // -- track fit to pick up hits that were lost in the initial search
          // -----------------------------------------------------------------------------------
          const float zhit         = goodHit ? protoTracks.z<scalar::float_v>( tscalar, iLayer ).cast() : zMidUT;
          const float xhit         = goodHit ? protoTracks.x<scalar::float_v>( tscalar, iLayer ).cast() : xArray[t2];
          const int   hitContIndex = protoTracks.hitContIndex<scalar::int_v>( tscalar ).cast();

          // -- The total sum of all plane codes is: 0 + 1 + 2 + 3 = 6
          // -- We can therefore get the plane code of the last pseudo-layer
          // -- as: 6 - sumOfAllOtherPlaneCodes
          const auto pC = goodHit ? planeCode( id.cast() ) : 6 - sumLayArray[t2];
          sumLayArray[t2] += pC;

          assert( pC > -1 && pC < 4 && "plane code for overlap out of bound" );

          const float txUTS = txArray[t2];

          const int begin = hitsInLayers[hitContIndex].layerIndices[pC];
          const int end =
              ( pC == 3 ) ? hitsInLayers[hitContIndex].size() : hitsInLayers[hitContIndex].layerIndices[pC + 1];
          const auto& hitsInL = hitsInLayers[hitContIndex].scalar();
          for ( int index2 = begin; index2 < end; ++index2 ) {
            const float zohit = hitsInL[index2].z().cast();
            if ( zohit == zhit ) continue;

            const float xohit   = hitsInL[index2].x().cast();
            const float xextrap = xhit + txUTS * ( zohit - zhit );
            if ( xohit - xextrap < -m_overlapTol ) continue;
            if ( xohit - xextrap > m_overlapTol ) break;

            if ( nUTHits[t2] >= int( LHCb::Pr::Upstream::Tracks::MaxUTHits ) )
              continue; // get this number from PrUpstreamTracks!!!
            const scalar::int_v utidx = hitsInL[index2].index();
            LHCb::LHCbID        oid( LHCb::Detector::UT::ChannelID( hitsInL[index2].channelID().cast() ) );
            const scalar::int_v lhcbid = bit_cast<int>( oid.lhcbID() );

            auto hits = outputTracks.scalar()[currentsize + trackIndex2].field<TracksTag::UTHits>();
            hits.resize( nUTHits[t2] + 1 );
            hits[nUTHits[t2]].template field<TracksTag::Index>().set( utidx );
            hits[nUTHits[t2]].template field<TracksTag::LHCbID>().set( LHCb::Event::lhcbid_v<scalar>( lhcbid ) );

            nUTHits[t2] += 1;
            // only one overlap hit
            break; // this should ensure there are never more than 8 hits on the track
          }
          trackIndex2++;
        }
      }
    }
  }
} // namespace LHCb::Pr
