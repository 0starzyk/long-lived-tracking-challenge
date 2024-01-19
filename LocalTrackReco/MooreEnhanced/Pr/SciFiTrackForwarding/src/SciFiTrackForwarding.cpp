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
#include "DetDesc/Condition.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "LHCbAlgs/Transformer.h"

#include "Event/PrLongTracks.h"
#include "Event/PrUpstreamTracks.h"
#include "Event/StateParameters.h"
#include "FTDet/DeFTDetector.h"

#include "Magnet/DeMagnet.h"
#include "PrKernel/PrFTZoneHandler.h"
#include "PrKernel/PrSelection.h"
#include "SciFiTrackForwardingHits.h"
#include "vdt/log.h"

#include <limits>
#include <string>
#include <tuple>

#include <boost/container/static_vector.hpp>

/** @class SciFiTrackForwarding SciFiTrackForwarding.cpp
 *
 *  \brief Alternative algorithm for the HLT1 Forward tracking of the upgrade
 *
 *   The general approach is to find a hit in each x-layer in station 3
 *   which is used to form a track seed from which one can forward this track
 *   into the other stations.
 *
 *   For further information see presentations in
 *   - https://indico.cern.ch/event/810764/
 *   - https://indico.cern.ch/event/786084/
 *
 *   In the current implementation, the seeding doublets are searched over two iterations.
 *   First iteration, performed over all VeloUT reconstructed tracks:
 *     doublets of hits are searched in the S3L0 --> S3L3 layers,
 *     they are then prolonged to S2 and S1.
 *   Second iteration, performed over the VeloUT tracks
 *     which have not been prolonged to the SciFi in the first iteration:
 *     doublets of hits are searched in the S2L0 --> S2L3 layers,
 *     they are then prolonged to S3 and S1.
 *
 *
 * FIXME
 * - The currently employed linear fit only uses x-hits.
 *   This leads to incresing momentum resolution for higher momentum tracks.
 */

///////////////////////////////////////////////////////////////////////////////
// anonymous namespace for helper functions local to this compilation unit
///////////////////////////////////////////////////////////////////////////////
using TracksTag = LHCb::Pr::Long::Tag;
namespace {

  using scalar = SIMDWrapper::scalar::types;
  using sF     = scalar::float_v;
  using sI     = scalar::int_v;

  using vec = SIMDWrapper::best::types;
  using vF  = vec::float_v;
  using vI  = vec::int_v;

  // stolen from Niklas
  // c++20's remove_cvref
  template <typename T>
  struct remove_cvref {
    using type = std::remove_cv_t<std::remove_reference_t<T>>;
  };

  template <typename T>
  using remove_cvref_t = typename remove_cvref<T>::type;

  template <typename T>
  auto to_std_array( T&& some_v ) {
    if constexpr ( std::is_same_v<remove_cvref_t<T>, vF> ) {
      std::array<float, vec::size> tmp;
      some_v.store( tmp.data() );
      return tmp;
    } else if ( std::is_same_v<remove_cvref_t<T>, vI> ) {
      std::array<int, vec::size> tmp;
      some_v.store( tmp.data() );
      return tmp;
    }
  }

  using TracksUT = LHCb::Pr::Upstream::Tracks;
  using TracksFT = LHCb::Pr::Long::Tracks;

  // constants for extrapolation polynomials from x hit in S3L0
  // to the corresponding x hit in other stations and layers
  double constexpr ExtFacS3_alpha = 1.4706824654297932e-5;
  double constexpr ExtFacS3_beta  = -3.152273942420754e-09;
  double constexpr ExtFacS3_gamma = -0.0003351012155394656;

  // the same as above, but starting from S2L0 extrapolating to other layers
  double constexpr ExtFacS2_alpha = 2.6591242069515776e-5;
  double constexpr ExtFacS2_beta  = -7.989227792416518e-09;
  double constexpr ExtFacS2_gamma = -0.000461618451581078;

  struct GeomCache {
    std::array<float, 24> LayerZPos{std::numeric_limits<float>::signaling_NaN()};
    std::array<float, 24> dZdYFiber{std::numeric_limits<float>::signaling_NaN()};
    std::array<float, 12> dXdYFiber{std::numeric_limits<float>::signaling_NaN()};

    // polynomial coefficients for extrapolations
    std::array<float, 48> ext_coeff{std::numeric_limits<float>::signaling_NaN()};

    FTZoneCache::PrFTZoneHandler zoneHandler;

    GeomCache() = default;

    GeomCache( DeFT const& ftDet ) : zoneHandler( ftDet ) {
      // get the layer z-positions and slopes to cache them for later usage
      for ( int i{0}; i < 12; ++i ) {
        // the first 12 values are the bottom layers
        LayerZPos[i] = zoneHandler.zone( 2 * i ).z();
        dZdYFiber[i] = zoneHandler.zone( 2 * i ).dzDy();

        // the last 12 values are the top layers
        LayerZPos[12 + i] = zoneHandler.zone( 2 * i + 1 ).z();
        dZdYFiber[12 + i] = zoneHandler.zone( 2 * i + 1 ).dzDy();
      }

      // get the uv layer fiber slopes to cache them for later usage
      for ( int i{0}; i < 6; ++i ) {
        // the first 6 values are the bottom uv layers
        dXdYFiber[i] = -zoneHandler.zone( LHCb::Detector::FT::stereoZones[2 * i] ).dxDy();

        // the last 6 values are the top uv layers
        dXdYFiber[6 + i] = -zoneHandler.zone( LHCb::Detector::FT::stereoZones[2 * i + 1] ).dxDy();
      }

      // calculate some coefficients for the extrapolation into other stations and put them in a cache
      int cnt = 0;

      // used in the first iteration, seeding from S3
      for ( int i : {0, 3, 4, 7} ) {
        ext_coeff[cnt]        = ExtFacS3_alpha * std::pow( LayerZPos[i] - LayerZPos[8], 2 );
        ext_coeff[12 + cnt++] = ExtFacS3_alpha * std::pow( LayerZPos[12 + i] - LayerZPos[20], 2 );

        ext_coeff[cnt]        = ExtFacS3_beta * std::pow( LayerZPos[i] - LayerZPos[8], 3 );
        ext_coeff[12 + cnt++] = ExtFacS3_beta * std::pow( LayerZPos[12 + i] - LayerZPos[20], 3 );

        ext_coeff[cnt]        = ExtFacS3_gamma * std::pow( LayerZPos[i] - LayerZPos[8], 2 );
        ext_coeff[12 + cnt++] = ExtFacS3_gamma * std::pow( LayerZPos[12 + i] - LayerZPos[20], 2 );
      }

      cnt = 24;
      // used in the second iteration, seeding from S2
      for ( int i : {0, 3, 8, 11} ) {
        int const from = i < 4 ? 4 : 7;

        ext_coeff[cnt]        = ExtFacS2_alpha * std::pow( LayerZPos[i] - LayerZPos[from], 2 );
        ext_coeff[12 + cnt++] = ExtFacS2_alpha * std::pow( LayerZPos[12 + i] - LayerZPos[from + 12], 2 );

        ext_coeff[cnt]        = ExtFacS2_beta * std::pow( LayerZPos[i] - LayerZPos[from], 3 );
        ext_coeff[12 + cnt++] = ExtFacS2_beta * std::pow( LayerZPos[12 + i] - LayerZPos[from + 12], 3 );

        ext_coeff[cnt]        = ExtFacS2_gamma * std::pow( LayerZPos[i] - LayerZPos[from], 2 );
        ext_coeff[12 + cnt++] = ExtFacS2_gamma * std::pow( LayerZPos[12 + i] - LayerZPos[from + 12], 2 );
      }
    }
  };

  // constants for the extrapolation polynomial from Velo to SciFi, used to determine searchwindows and minPT cut border
  std::array<float, 8> constexpr toSciFiExtParams{4824.31956565f,  426.26974766f,   7071.08408876f, 12080.38364257f,
                                                  14077.79607408f, 13909.31561208f, 9315.34184959f, 3209.49021545f};

  // constants for the y-correction polynomial
  std::array<float, 8> constexpr deltaYParams{3.78837f, 73.1636f, 7353.89f,  -6347.68f,
                                              20270.3f, 3721.02f, -46038.2f, 230943.f};

  // parameters used for the momentum determination at the very end
  std::array<float, 8> constexpr MomentumParams{1239.4073749458162, 486.05664058906814, 6.7158701518424815,
                                                632.7283787142547,  2358.5758035677504, -9256.27946160669,
                                                241.4601040854867,  42.04859549174048};

  // Parameters for the linear discriminant which we use to reject ghosts
  std::array<float, 8> LDAParams{5.35770606, 1.22786675, -1.17615119, 2.41396581,
                                 2.09624748, 1.7029565,  -2.69818528, 1.88471171};

  // internal helper class to keep track of our best candidate
  struct SciFiTrackForwardingCand {
    // array used to store indices of hits
    boost::container::static_vector<int, 12> ids;
    float                                    PQ{0.f};    // Momentum times charge of candidate
    float                                    finaltx{0}; // x-slope at the first x-hit in station 3 after linear fit
    float                                    newx0{0};   // x position of x-hit in S3 after linear fit
    // float               quality{50.f};  // the default value of quality acts as cut in the later selection of
    // candidates
    float quality{-4.5f}; // the default value of quality acts as cut in the later selection of candidates
    int   numHits{0};     // number of hits this candidate holds
  };

  // Custom span that allow for negative indices
  // needed for the upper_bound below
  template <typename T>
  class span {
  public:
    span() : m_data( nullptr ), m_size( 0 ) {}
    span( T* data, int size ) : m_data( data ), m_size( size ) {}
    T&  operator[]( int i ) const { return m_data[i]; }
    int size() const { return m_size; }

  private:
    T*  m_data;
    int m_size;
  };

  // modified fast_upper_bound from Arthur
  template <class T>
  int get_closest_hit_idx_binary( span<const T> const& hitvec, T const value ) {
    int size = hitvec.size();
    int low  = 0;

    int half1 = size / 2;
    low += ( hitvec[low + half1] <= value ) * ( size - half1 );
    size      = half1;
    int half2 = size / 2;
    low += ( hitvec[low + half2] <= value ) * ( size - half2 );
    size      = half2;
    int half3 = size / 2;
    low += ( hitvec[low + half3] <= value ) * ( size - half3 );
    size      = half3;
    int half4 = size / 2;
    low += ( hitvec[low + half4] <= value ) * ( size - half4 );
    size      = half4;
    int half5 = size / 2;
    low += ( hitvec[low + half5] <= value ) * ( size - half5 );
    size = half5;

    do {
      int half = size / 2;
      low += ( hitvec[low + half] <= value ) * ( size - half );
      size = half;
    } while ( size > 0 );

    return low - ( std::abs( hitvec[low] - value ) >= std::abs( hitvec[low - 1] - value ) );
  }

  // Simple linear search
  // since my vector container is padded by sentinels I can get rid of a check in the loop and simply advance until
  // condition is met
  [[using gnu: hot, const]] inline int get_closest_hit_idx( SciFiTrackForwardingHits::hits_t const& hitvec, int start,
                                                            float const val ) noexcept {

    // vector of comparison value
    vF vval{val};
    for ( ;; ) {
      // load vector of data
      vF vvec{hitvec.data() + start};
      // comparison mask with true/false values
      auto const mask = vval > vvec;

      // if less than simd::size whe can stop
      if ( !all( mask ) ) {
        start += popcount( mask );
        break;
      }
      // all values where smaller thus do next iteration
      start += vec::size;
    }
    return start - ( hitvec[start] - val >= val - hitvec[start - 1] );
  }

  // helper function to fill a vector with some information later used in the linear fit
  // only to reduce code duplication and readability later on
  [[using gnu: always_inline]] void inline fillFitVec( std::array<float, 24>& fitvec, int& fitcnt, float const x,
                                                       float const zdelta, float const curve, float const residual ) {
    fitvec[fitcnt++] = x;
    fitvec[fitcnt++] = zdelta;
    fitvec[fitcnt++] = curve;
    fitvec[fitcnt++] = residual;
  }

  std::tuple<float, float, float> simple_fit( std::array<float, 24> const& fitvec, int const fitcnt, float const newtx,
                                              float const firstXhit ) {
    // quick linear fit in x0 and tx, curvature terms remain fixed
    float s0  = 0;
    float sz  = 0;
    float sz2 = 0;
    float sd  = 0;
    float sdz = 0;

    for ( int idx{0}; idx < fitcnt; idx += 4 ) {
      // for now the weight is just 1
      s0 += 1;
      sz += fitvec[idx + 1];
      sz2 += fitvec[idx + 1] * fitvec[idx + 1];
      sd += fitvec[idx + 3];
      sdz += fitvec[idx + 3] * fitvec[idx + 1];
    }

    float const den     = 1.f / ( sz * sz - 6 * sz2 );
    float const xoff    = ( sdz * sz - sd * sz2 ) * den;
    float const txoff   = ( sd * sz - 6 * sdz ) * den;
    float const finaltx = newtx + txoff;
    float const newx0   = firstXhit + xoff;

    // calculate a chi2 after the fit
    float chi2 = ( fitvec[0] - newx0 ) * ( fitvec[0] - newx0 );
    for ( int idx{4}; idx < fitcnt; idx += 4 ) {
      // for now the weight is just 1
      float const diff = fitvec[idx] - ( newx0 + fitvec[idx + 1] * finaltx + fitvec[idx + 2] );
      chi2 += diff * diff;
    }
    return {finaltx, newx0, chi2};
  }

  template <bool is_first_loop>
  float calc_mom( float const finaltx, float const endv_tx, float const endv_ty2, float const mag_scale,
                  float const factor ) {

    // determine momentum * charge estimate, polynomial was fitted on MC
    const float finaltx2 = finaltx * finaltx;
    const float coef =
        ( MomentumParams[0] + finaltx2 * ( MomentumParams[1] + MomentumParams[2] * finaltx2 ) +
          MomentumParams[3] * finaltx * endv_tx + endv_ty2 * ( MomentumParams[4] + MomentumParams[5] * endv_ty2 ) +
          MomentumParams[6] * endv_tx * endv_tx );

    if constexpr ( is_first_loop ) {
      return mag_scale * coef / ( finaltx - endv_tx ) + factor * MomentumParams[7];
    } else {
      return 0.9818f * ( mag_scale * coef / ( finaltx - endv_tx ) + factor * MomentumParams[7] );
    }
  }

  float calc_quality( float const endv_pt, float const factor, float const firstXhit, float const straighExt,
                      float const endv_pq, float const endv_qp, float const candPQ, float const wrong_sign_pt,
                      float const chi2, float const endv_txy, float const S3UVDelta, float const S2UVDelta,
                      float const S1UVDelta, int const numuvhits, int const numxhits ) {

    // the condition evaluates to true if the track is a charge flipped one
    float deltaMom = ( endv_pt > wrong_sign_pt and ( factor * ( firstXhit - straighExt ) < 0 ) )
                         ? 3.f * std::abs( endv_qp * ( endv_pq + candPQ ) )
                         : std::abs( endv_qp * ( endv_pq - candPQ ) );

    // counteract worse momentum resolution at higher momenutm by scaling down
    // the observed delta
    if ( std::abs( candPQ ) > 40000 ) deltaMom *= 40000.f / std::abs( candPQ );

    // For now a linear discriminant seems to do well enough. But if needed this
    // could be substituted by a neural net for better rejection of ghosts
    return LDAParams[0] * approx_log( sF( 1.f + deltaMom ) ).cast() +
           LDAParams[1] * approx_log( sF( 1.f + chi2 ) ).cast() +
           LDAParams[2] *
               approx_log( sF( std::abs( candPQ ) * endv_txy / std::sqrt( 1 + endv_txy * endv_txy ) ) ).cast() +
           LDAParams[3] * approx_log( sF( S3UVDelta ) ).cast() + LDAParams[4] * approx_log( sF( S2UVDelta ) ).cast() +
           LDAParams[5] * approx_log( sF( S1UVDelta ) ).cast() + LDAParams[6] * static_cast<float>( numuvhits ) +
           LDAParams[7] * static_cast<float>( numxhits );
  }

  struct config_loop_1 {
    // average z-position of the kink position inside the magnet (tuned on MC),
    float static constexpr zPosKinkMagnet = 5282.5f;
    // idx of first layer of seed station in Z-position array
    int static constexpr seed_first_layer_z_idx{8};
    int static constexpr next_first_layer_z_idx{4};
    int static constexpr last_first_layer_z_idx{0};
    // the array of extrapolation coefficients is 48 entries long
    // first 24 are for the first loop 24-47 for the second loop
    int static constexpr ext_coeff_start{0};
    // Config 1 starts seeding in last station
    // indices to select layers of third station
    int static constexpr seed_station_idx_first{4};
    int static constexpr seed_station_idx_last{5};
    // station 2 is the first station to be extrapolated into
    int static constexpr next_station_idx_first{2};
    int static constexpr next_station_idx_last{3};
    float static constexpr doublet_curve_c0{-0.231615f};
    float static constexpr doublet_curve_c1{33.1256f};
    float static constexpr doublet_curve_c2{10.4693f};
    float static constexpr seed_scale_ycorr{1.005f};
    float static constexpr next_scale_ycorr{0.83f};
    float static constexpr last_scale_ycorr{0.615f};
    // cut window parameters
    float static constexpr doublet_win_slope = 4.25f;
    float static constexpr doublet_win_off   = 0.48f;
    float static constexpr doublet_win_max   = 3.f;
    float static constexpr seed_uv_win_slope = 10.f;
    float static constexpr seed_uv_win_off   = .75f;

    float static constexpr next_l0_win_slope = 2.f;
    float static constexpr next_l0_win_off   = 2.f;
    float static constexpr next_l3_win_slope = 1.f;
    float static constexpr next_l3_win_off   = 1.8f;
    float static constexpr next_uv_win_slope = 10.f;
    float static constexpr next_uv_win_off   = .75f;

    float static constexpr last_l0_win_slope = 2.6f;
    float static constexpr last_l0_win_off   = 3.6f;
    float static constexpr last_l3_win_slope = 2.6f;
    float static constexpr last_l3_win_off   = 3.5f;
    float static constexpr last_uv_win_slope = 7.f;
    float static constexpr last_uv_win_off   = 1.f;

    int static constexpr last_station_idx_first{0};
    int static constexpr last_station_idx_last{1};
  };

  struct config_loop_2 {
    float static constexpr zPosKinkMagnet = 5195.f;
    int static constexpr seed_first_layer_z_idx{4};
    int static constexpr next_first_layer_z_idx{8};
    int static constexpr last_first_layer_z_idx{0};
    // the array of extrapolation coefficients is 48 entries long
    // first 24 are for the first loop 24-47 for the second loop
    int static constexpr ext_coeff_start{24};
    // Config 1 starts seeding in last station
    // indices to select layers of third station
    int static constexpr seed_station_idx_first{2};
    int static constexpr seed_station_idx_last{3};
    // station 2 is the first station to be extrapolated into
    int static constexpr next_station_idx_first{4};
    int static constexpr next_station_idx_last{5};
    float static constexpr doublet_curve_c0{-1.050f};
    float static constexpr doublet_curve_c1{17.74f};
    float static constexpr doublet_curve_c2{20.62};
    float static constexpr seed_scale_ycorr{0.83f};
    float static constexpr next_scale_ycorr{1.005f};
    float static constexpr last_scale_ycorr{0.615f};
    // cut window parameters
    float static constexpr doublet_win_slope = 4.3f;
    float static constexpr doublet_win_off   = 0.6f;
    float static constexpr doublet_win_max   = 3.2f;
    float static constexpr seed_uv_win_slope = 10.f;
    float static constexpr seed_uv_win_off   = 1.f;

    float static constexpr next_l0_win_slope = 1.5f;
    float static constexpr next_l0_win_off   = 2.3f;
    float static constexpr next_l3_win_slope = 1.5f;
    float static constexpr next_l3_win_off   = 1.5f;
    float static constexpr next_uv_win_slope = 10.f;
    float static constexpr next_uv_win_off   = 1.f;

    float static constexpr last_l0_win_slope = 2.5f;
    float static constexpr last_l0_win_off   = 2.f;
    float static constexpr last_l3_win_slope = 1.5f;
    float static constexpr last_l3_win_off   = 2.3f;
    float static constexpr last_uv_win_slope = 7.5f;
    float static constexpr last_uv_win_off   = 1.5f;
    int static constexpr last_station_idx_first{0};
    int static constexpr last_station_idx_last{1};
  };

} // namespace

///////////////////////////////////////////////////////////////////////////////
//
//  Fast Upgrade Forward Tracking Algorithm
//
//  The general approach is to find a hit in each x-layer in station 3
//  which is used to form a track seed from which one can forward this track
//  into the other stations.
//
//  For further information see presentations in
//  - https://indico.cern.ch/event/810764/
//  - https://indico.cern.ch/event/786084/
//
//
///////////////////////////////////////////////////////////////////////////////
class SciFiTrackForwarding
    : public LHCb::Algorithm::Transformer<TracksFT( EventContext const&, SciFiTrackForwardingHits const&,
                                                    TracksUT const&, GeomCache const&, DeMagnet const& ),
                                          LHCb::DetDesc::usesConditions<GeomCache, DeMagnet>> {

public:
  SciFiTrackForwarding( std::string const& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {KeyValue{"HitsLocation", "Rec/SciFiTrackForwarding/Hits"},
                      KeyValue{"InputTracks", "Rec/Track/UT"},
                      KeyValue{"GeometryCache", "FT:" + name + "-ZoneHandlerCache"},
                      KeyValue{"Magnet", LHCb::Det::Magnet::det_path}},
                     KeyValue{"Output", "Rec/Track/FT"} ) {}

  StatusCode initialize() override {
    return Transformer::initialize().andThen( [&] {
      addConditionDerivation<GeomCache( DeFT const& )>( {DeFTDetectorLocation::Default}, inputLocation<GeomCache>() );
    } );
  }

  TracksFT operator()( EventContext const&, SciFiTrackForwardingHits const&, TracksUT const&, GeomCache const&,
                       DeMagnet const& ) const override;

  template <typename config>
  [[using gnu: const, hot]] SciFiTrackForwardingCand const
  find_track( GeomCache const&, SciFiTrackForwardingHits const&, int const, Vec3<SIMDWrapper::scalar::float_v> const&,
              Vec3<SIMDWrapper::scalar::float_v> const&, float const, float const, float const, float const,
              float const, float const, std::array<float, 12> const& ) const;

  mutable Gaudi::Accumulators::SummingCounter<unsigned> m_counter_created_tracks{this, "Created long tracks"};
  mutable Gaudi::Accumulators::Counter<>                m_counter_2_loop{this, "2nd Loop"};

  Gaudi::Property<float> p_minPt{this, "MinPt", 490.f};

  // if  a track exhibits a higher PT we open a symmetric search window
  // this is to counter the possibility of a wrong charge estimate from the UT
  Gaudi::Property<float> p_wrongSignPT{this, "WrongSignPT", 2000.f};

  // parameters for the momentum dependent upper error limit on the x-search window estimation
  Gaudi::Property<float> p_UpperLimitOffset{this, "UpperLimit_offset", 100.f};
  Gaudi::Property<float> p_UpperLimitSlope{this, "UpperLimit_slope", 2800.f};
  Gaudi::Property<float> p_UpperLimitMax{this, "UpperLimit_max", 600.f};
  Gaudi::Property<float> p_UpperLimitMin{this, "UpperLimit_min", 150.f};

  // same as above for the lower limit
  Gaudi::Property<float> p_LowerLimitOffset{this, "LowerLimit_offset", 50.f};
  Gaudi::Property<float> p_LowerLimitSlope{this, "LowerLimit_slope", 1400.f};
  Gaudi::Property<float> p_LowerLimitMax{this, "LowerLimit_max", 600.f};

  // second loop for searching doublets
  Gaudi::Property<bool> p_SecondLoop{this, "SecondLoop", true};

private:
  // this enables me to print things automagically space separated (less typing)
  // also when enable is set to false, all of this code is removed which is nice for me to test what this actually
  // costs
  template <bool enable = false, typename... Args>
  void mydebug( Args&&... args ) const {
    if constexpr ( enable ) {
      if ( msgLevel( MSG::DEBUG ) ) { ( ( debug() << std::forward<Args>( args ) << " " ), ... ) << endmsg; }
    }
  }
};

DECLARE_COMPONENT( SciFiTrackForwarding )

TracksFT SciFiTrackForwarding::operator()( EventContext const& evtCtx, SciFiTrackForwardingHits const& hithandler,
                                           TracksUT const& tracks, GeomCache const& cache,
                                           DeMagnet const& magnet ) const {

  // factor of -1 because everything was trained with magdown, so for magdown we want factor = 1, magup = -1
  const float magscalefactor = -1 * magnet.signedRelativeCurrent();

  TracksFT Output{tracks.getVeloAncestors(),
                  &tracks,
                  nullptr,
                  LHCb::Event::Enum::Track::History::SciFiTrackForwarding,
                  Zipping::generateZipIdentifier(),
                  LHCb::getMemResource( evtCtx )};
  Output.reserve( tracks.size() );

  mydebug( "LayerZPos", cache.LayerZPos );

  // I want this in GeV but default in LHCb is MeV so this is a compromise to have the property in MeV
  float const minPTCutValue = p_minPt / float{Gaudi::Units::GeV};

  auto buffer_2_loop = m_counter_2_loop.buffer();
  // start our loop over UT tracks

  auto const utzipped          = tracks.simd();
  auto const vec_mag_scale     = vF{magscalefactor};
  auto const neg_vec_mag_scale = vF{-1.f * magscalefactor};
  for ( auto const& track : utzipped ) {

    auto loop_mask = track.loop_mask();
    // placeholder for the best candidate
    mydebug( "++++++++++ Start new UT track ++++++++++++++++" );

    auto const& endv_pos = track.StatePos();
    auto const& endv_dir = track.StateDir();

    // everything I need from the end velo state
    auto const endv_tx   = endv_dir.x;
    auto const endv_ty   = endv_dir.y;
    auto const endv_qp   = track.qOverP();
    auto const endv_y    = endv_pos.y;
    auto const endv_z    = endv_pos.z;
    auto const endv_ty2  = endv_ty * endv_ty;
    auto const endv_txy2 = endv_tx * endv_tx + endv_ty2;

    // determine if a track is expected to intersect the seeding station from the top or the bottom
    // check is performed at last layer of second station.
    auto const ones                = vI{1};
    auto const zeroes              = vI{0};
    auto const track_in_upper_half = signselect( ( endv_y + ( cache.LayerZPos[8] - endv_z ) * endv_ty ), ones, zeroes );
    // offset for LayerZPos & ext_coeff and dxdy array
    auto const up_down_offset = to_std_array( track_in_upper_half * ( cache.LayerZPos.size() / 2 ) );

    // charge * magnet polarity factor needed in various polynomials below
    auto const factor = signselect( endv_qp, vec_mag_scale, neg_vec_mag_scale );

    // charge over momentum is usally in MeV so make it GeV
    auto const qOpGeV = endv_qp * float{Gaudi::Units::GeV} * magscalefactor;
    // calculate an upper error boundary on the tracks x-prediction,
    auto const upper_error =
        min( max( p_UpperLimitOffset.value() + p_UpperLimitSlope.value() * abs( qOpGeV ), p_UpperLimitMin.value() ),
             p_UpperLimitMax.value() );

    // calculate a lower error boundary on the tracks x-prediction,
    auto const lower_error =
        min( p_LowerLimitOffset.value() + p_LowerLimitSlope.value() * abs( qOpGeV ), p_LowerLimitMax.value() );

    // extrapolate the track x-position into the scifi to determine search windows
    // common term for both extrapolations below
    auto const term1 =
        toSciFiExtParams[0] + endv_tx * ( -toSciFiExtParams[1] * factor + toSciFiExtParams[2] * endv_tx ) +
        endv_ty2 * ( toSciFiExtParams[3] + endv_tx * ( toSciFiExtParams[4] * factor + toSciFiExtParams[5] * endv_tx ) );

    // prediction of tracks x position with ut momentum estimate
    auto const xExt = qOpGeV * ( term1 + qOpGeV * ( toSciFiExtParams[6] * endv_tx + toSciFiExtParams[7] * qOpGeV ) );

    // 1/p, given a minPT cut and the tracks slopes
    auto const minInvPGeV = factor / minPTCutValue * sqrt( endv_txy2 ) / sqrt( 1.f + endv_txy2 );

    // given the above value of 1/p this should be the tracks x-position
    // we can use this to make our windows smaller given a minPT cut
    auto const min_pt_border =
        minInvPGeV * ( term1 + minInvPGeV * ( toSciFiExtParams[6] * endv_tx + toSciFiExtParams[7] * minInvPGeV ) );

    auto const a_qp                  = to_std_array( endv_qp );
    auto const a_factor              = to_std_array( factor );
    auto const a_upper_error         = to_std_array( upper_error );
    auto const a_lower_error         = to_std_array( lower_error );
    auto const a_xExt                = to_std_array( xExt );
    auto const a_min_pt_border       = to_std_array( min_pt_border );
    auto const a_track_in_upper_half = to_std_array( track_in_upper_half );

    auto const offset  = track.offset();
    auto const uttrack = utzipped.with<SIMDWrapper::InstructionSet::Scalar>();
    for ( size_t tr{0}; tr < track.width(); ++tr ) {
      auto const track_scalar = uttrack[offset + tr];
      if ( !testbit( loop_mask, tr ) ) { break; }
      // extrapolate velo track y-position into all 12 Layers
      // IIFE (Immediately invoked function expression) to keep yAtZ const
      // compiler explorer tested that this isn't producing overhead over simple loop
      auto const& scalar_endv_pos = track_scalar.StatePos();
      auto const& scalar_endv_dir = track_scalar.StateDir();

      auto const yAtZ = [&]() {
        std::array<float, 12> tmp;
        std::transform( cache.LayerZPos.begin() + up_down_offset[tr], cache.LayerZPos.begin() + 12 + up_down_offset[tr],
                        tmp.begin(), [=]( float const z ) {
                          return scalar_endv_pos.y.cast() + ( z - scalar_endv_pos.z.cast() ) * scalar_endv_dir.y.cast();
                        } );
        return tmp;
      }();

      mydebug( "yatz: ", yAtZ );

      auto bestcandidate{find_track<config_loop_1>( cache, hithandler, a_track_in_upper_half[tr], scalar_endv_pos,
                                                    scalar_endv_dir, a_qp[tr], a_factor[tr], a_upper_error[tr],
                                                    a_lower_error[tr], a_xExt[tr], a_min_pt_border[tr], yAtZ )};

      mydebug( "Search 1 done!" );
      if ( bestcandidate.ids.empty() && p_SecondLoop.value() ) {
        mydebug( "Search 2nd loop!" );
        ++buffer_2_loop;
        bestcandidate =
            find_track<config_loop_2>( cache, hithandler, a_track_in_upper_half[tr], scalar_endv_pos, scalar_endv_dir,
                                       a_qp[tr], a_factor[tr], 0.87f * a_upper_error[tr], 0.87f * a_lower_error[tr],
                                       0.87f * a_xExt[tr], 0.87f * a_min_pt_border[tr], yAtZ );
      }
      mydebug( "Search done!" );

      if ( !bestcandidate.ids.empty() ) {
        mydebug( "saving track: ", Output.size() );
        auto out = Output.emplace_back<SIMDWrapper::InstructionSet::Scalar>();

        out.field<TracksTag::trackVP>().set( track_scalar.trackVP() );
        out.field<TracksTag::trackUT>().set( offset + tr );
        out.field<TracksTag::trackSeed>().set( -1 );

        sF const qop = 1.f / bestcandidate.PQ;
        out.field<TracksTag::States>( 0 ).setQOverP( qop );
        out.field<TracksTag::States>( 1 ).setQOverP( qop );

        // EndV State
        out.field<TracksTag::States>( 0 ).setPosition( scalar_endv_pos );
        out.field<TracksTag::States>( 0 ).setDirection( scalar_endv_dir );

        // AtT State
        float const endT_z   = cache.LayerZPos[8];
        float const endT_x   = bestcandidate.newx0;
        auto        endT_y   = scalar_endv_pos.y + scalar_endv_dir.y * ( endT_z - scalar_endv_pos.z );
        auto        endT_pos = Vec3<sF>( endT_x, endT_y, endT_z );
        auto        endT_dir = Vec3<sF>( bestcandidate.finaltx, scalar_endv_dir.y, 1.f );
        out.field<TracksTag::States>( 1 ).setPosition( endT_pos );
        out.field<TracksTag::States>( 1 ).setDirection( endT_dir );

        // store Velo hit indices
        const int n_vphits = track_scalar.nVPHits().cast();
        const int n_uthits = track_scalar.nUTHits().cast();
        out.field<TracksTag::VPHits>().resize( n_vphits );
        out.field<TracksTag::UTHits>().resize( n_uthits );

        for ( auto idx{0}; idx < n_vphits; ++idx ) {
          out.field<TracksTag::VPHits>()[idx].template field<TracksTag::Index>().set( track_scalar.vp_index( idx ) );
          out.field<TracksTag::VPHits>()[idx].template field<TracksTag::LHCbID>().set( track_scalar.vp_lhcbID( idx ) );
        }

        for ( auto idx{0}; idx < n_uthits; ++idx ) {
          out.field<TracksTag::UTHits>()[idx].template field<TracksTag::Index>().set( track_scalar.ut_index( idx ) );
          out.field<TracksTag::UTHits>()[idx].template field<TracksTag::LHCbID>().set( track_scalar.ut_lhcbID( idx ) );
        }

        const sI n_fthits = bestcandidate.ids.size();
        out.field<TracksTag::FTHits>().resize( n_fthits );
        std::size_t n_hits = 0;
        /// FT hit indices & lhcbID
        for ( auto idx{bestcandidate.ids.begin()}; idx != bestcandidate.ids.end(); ++idx, ++n_hits ) {
          out.field<TracksTag::FTHits>()[n_hits].template field<TracksTag::Index>().set( *idx );
          out.field<TracksTag::FTHits>()[n_hits].template field<TracksTag::LHCbID>().set(
              LHCb::LHCbID( hithandler.IDs[*idx] ) );
        }
      } // save the bestcandidate if we built one
    }
  } // loop over the UT tracks

  m_counter_created_tracks += Output.size();
  return Output;
}

template <typename config>
[[using gnu: const, hot]] SciFiTrackForwardingCand const SciFiTrackForwarding::find_track(
    GeomCache const& cache, SciFiTrackForwardingHits const& hithandler, int const track_in_upper_half,
    Vec3<SIMDWrapper::scalar::float_v> const& endv_pos, Vec3<SIMDWrapper::scalar::float_v> const& endv_dir,
    float const endv_qp, float const factor, float const upper_error, float const lower_error, float const xExt,
    float const min_pt_border, std::array<float, 12> const& yAtZ ) const {

  using c = config;

  float const endv_x    = endv_pos.x.cast();
  float const endv_tx   = endv_dir.x.cast();
  float const endv_ty   = endv_dir.y.cast();
  float const endv_z    = endv_pos.z.cast();
  float const endv_tx2  = endv_tx * endv_tx;
  float const endv_ty2  = endv_ty * endv_ty;
  float const endv_txy2 = endv_tx2 + endv_ty2;
  float const endv_txy  = std::sqrt( endv_txy2 );
  float const endv_pq   = 1.f / endv_qp;
  float const endv_p    = std::abs( endv_pq );
  float const endv_pz   = endv_p / std::sqrt( 1.f + endv_txy2 );
  float const endv_pt   = endv_pz * endv_txy;

  SciFiTrackForwardingCand bestcandidate{};

  int const   up_down_offset        = track_in_upper_half * ( cache.LayerZPos.size() / 2 );
  int const   up_down_dxdy_offset   = track_in_upper_half * ( cache.dXdYFiber.size() / 2 );
  float const seed_station_z_pos_l0 = cache.LayerZPos[up_down_offset + c::seed_first_layer_z_idx];

  // extrapolate velo track x-position into magnet,
  float const x_at_magnet = endv_x + ( c::zPosKinkMagnet - endv_z ) * endv_tx;

  // extrapolate the velo straight into the SciFi
  float const straighExt = endv_x + ( seed_station_z_pos_l0 - endv_z ) * endv_tx;

  // depending on the sign of q set the lower and upper error in the right direction
  // on top of that make sure that the upper error isn't larger than the min_pt_border
  // and the lower error isn't bigger than xExt which means we don't look further than the straight line prediction
  float xMin = factor > 0 ? straighExt + ( ( xExt < lower_error ) ? 0 : ( xExt - lower_error ) )
                          : straighExt + std::max( xExt - upper_error, min_pt_border );

  float xMax = factor > 0 ? straighExt + std::min( xExt + upper_error, min_pt_border )
                          : straighExt + ( ( xExt > -lower_error ) ? 0 : ( xExt + lower_error ) );

  if ( endv_pt > p_wrongSignPT.value() ) {
    xMin = straighExt - std::abs( xExt ) - upper_error;
    xMax = straighExt + std::abs( xExt ) + upper_error;
  }

  mydebug( "==================" );
  mydebug( "straight: ", straighExt, "xExt ", xExt, " dxref ", min_pt_border, "lowError", lower_error,
           "upError: ", upper_error, "xmin", xMin, "xmax: ", xMax );

  // set start and end points in 1D hit array depending on if a track is in upper or lower detector half
  mydebug( "Selecting upper or lower  hits", track_in_upper_half );
  auto const Xzones  = track_in_upper_half ? LHCb::Detector::FT::xZonesUpper : LHCb::Detector::FT::xZonesLower;
  auto const UVzones = track_in_upper_half ? LHCb::Detector::FT::uvZonesUpper : LHCb::Detector::FT::uvZonesLower;

  // a range is a pair of start and size
  auto const seed_station_range_l0 = hithandler.zonerange[Xzones[c::seed_station_idx_first]];

  int const start_s3_l0 =
      seed_station_range_l0.first +
      get_closest_hit_idx_binary( {hithandler.hits.data() + seed_station_range_l0.first, seed_station_range_l0.second},
                                  xMin );

  int const end_s3_l0 = start_s3_l0 + get_closest_hit_idx_binary(
                                          {hithandler.hits.data() + start_s3_l0,
                                           seed_station_range_l0.second - start_s3_l0 + seed_station_range_l0.first},
                                          xMax );

  int start_seed_station_l1 = hithandler.zonerange[UVzones[c::seed_station_idx_first]].first;
  int start_seed_station_l2 = hithandler.zonerange[UVzones[c::seed_station_idx_last]].first;
  int start_seed_station_l3 = hithandler.zonerange[Xzones[c::seed_station_idx_last]].first;

  int start_next_station_l0 = hithandler.zonerange[Xzones[c::next_station_idx_first]].first;
  int start_next_station_l1 = hithandler.zonerange[UVzones[c::next_station_idx_first]].first;
  int start_next_station_l2 = hithandler.zonerange[UVzones[c::next_station_idx_last]].first;
  int start_next_station_l3 = hithandler.zonerange[Xzones[c::next_station_idx_last]].first;

  // each variable is a pair [startIdx, size]
  auto const last_station_range_l0 = hithandler.zonerange[Xzones[c::last_station_idx_first]];
  auto const last_station_range_l1 = hithandler.zonerange[UVzones[c::last_station_idx_first]];
  auto const last_station_range_l2 = hithandler.zonerange[UVzones[c::last_station_idx_last]];
  auto const last_station_range_l3 = hithandler.zonerange[Xzones[c::last_station_idx_last]];

  // z delta between the two x-layers of the seeding station
  float const zdelta_seed_doublet =
      cache.LayerZPos[up_down_offset + 3 + c::seed_first_layer_z_idx] - seed_station_z_pos_l0;

  // inverse of the above quantity
  float const inv_zdelta_seed_doublet = 1.f / zdelta_seed_doublet;

  // delta in Z between kink position in magnet and z position of seed station
  float const inv_zdelta_kink_station = 1.f / ( seed_station_z_pos_l0 - c::zPosKinkMagnet );

  // Loop over the hits in the first x-layer within the search window
  for ( auto idxL0{start_s3_l0}; idxL0 < end_s3_l0; ++idxL0 ) {

    mydebug( "#####################Start new hit################################" );

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //                  Start of search in the first station               //
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // get the x value of hit
    float const firstXhit = hithandler.hits[idxL0];

    // calculate a rough estimate of the expected slope,
    // depending on the current iteration
    float const tx = ( firstXhit - x_at_magnet ) * inv_zdelta_kink_station;

    // this is a small correction mainly for lower momentum tracks for which the zkink isn't perfect
    // and they'll have some bending. This is solely empirical.
    float const doubletCurve =
        ( tx - endv_tx ) *
        ( c::doublet_curve_c0 + ( endv_tx * c::doublet_curve_c1 + tx * c::doublet_curve_c2 ) * ( tx - endv_tx ) );

    // piece together the x-prediction in the last x layer
    float const secondXhitPred = firstXhit + tx * zdelta_seed_doublet + doubletCurve;

    mydebug( "firstXhit", firstXhit, "lhcbid", hithandler.IDs[idxL0], "secondXhitPred", secondXhitPred, "tx", tx );
    mydebug( doubletCurve, inv_zdelta_kink_station, seed_station_z_pos_l0, c::zPosKinkMagnet, zdelta_seed_doublet );

    start_seed_station_l3 = get_closest_hit_idx( hithandler.hits, start_seed_station_l3, secondXhitPred );
    auto const delta      = std::abs( hithandler.hits[start_seed_station_l3] - secondXhitPred );

    // x position of the second hit of the doublet
    float const secondXhit = hithandler.hits[start_seed_station_l3];

    // calculate x-slope of the found doublet
    float const newtx = ( secondXhit - firstXhit ) * inv_zdelta_seed_doublet;

    // deltaSlope is proportional to the momentum and used in many places to tune momentum dependant search windows
    float const deltaSlope = newtx - endv_tx;

    // the minimum value of 0.011 avoids too small search windows for tracks with P > 100 GeV
    float const absDSlope = std::max( 0.011f, std::abs( deltaSlope ) );
    mydebug( "absdslope", absDSlope );

    // check if the x-hit is within our search window acceptance, depending on the current iteration
    if ( delta > std::min( c::doublet_win_off + c::doublet_win_slope * absDSlope, c::doublet_win_max ) ) continue;

    // we have a doublet and are going to try and build a track
    // for that we need a couple tmp objects to keep track of stuff
    // tmpidvec will store the indices of the hit to refind it later
    boost::container::static_vector<int, 12> tmpidvec{start_seed_station_l3, idxL0};

    // this is a vector that will remember some values for the linear fit performed at the very end.
    // for each hit this will store 4 values, [xvalue, zdelta, curvature, residual]
    // where zdelta is the distance between the hit and the hit in layer 0 of the seeding station
    // curvature is the non-linear part of the prediction. See x3curve variable for example
    std::array<float, 24> fitvec{firstXhit, 0, 0, 0, secondXhit, zdelta_seed_doublet, doubletCurve, 0};

    // fitcnt simply keeps track of what's in the array
    int fitcnt = 8;

    float const direction = -1.f * std::copysign( 1.f, deltaSlope );

    // correct the straight line estimate of the y-position at station 3 layer 0
    float const ycorr =
        absDSlope * ( direction * deltaYParams[0] +
                      endv_ty * ( deltaYParams[1] + deltaYParams[5] * endv_ty2 +
                                  endv_tx * ( direction * ( deltaYParams[2] + deltaYParams[6] * endv_ty2 ) +
                                              endv_tx * ( deltaYParams[3] + deltaYParams[7] * endv_ty2 +
                                                          deltaYParams[4] * direction * endv_tx ) ) ) );

    // keep track of how many hits we find
    int numuvhits{0};
    int numxhits{2};

    // calculate the x-predictions for the U V layers in station 3
    float const seed_l1_pred =
        firstXhit +
        newtx * ( cache.LayerZPos[up_down_offset + 1 + c::seed_first_layer_z_idx] - seed_station_z_pos_l0 ) +
        ( yAtZ[1 + c::seed_first_layer_z_idx] - c::seed_scale_ycorr * ycorr ) *
            cache.dXdYFiber[up_down_dxdy_offset + c::seed_first_layer_z_idx / 2];

    mydebug( newtx, ( cache.LayerZPos[up_down_offset + 1 + c::seed_first_layer_z_idx] - seed_station_z_pos_l0 ),
             ( yAtZ[1 + c::seed_first_layer_z_idx] - c::seed_scale_ycorr * ycorr ),
             cache.dXdYFiber[up_down_dxdy_offset + c::seed_first_layer_z_idx / 2] );

    float const seed_l2_pred =
        firstXhit +
        newtx * ( cache.LayerZPos[up_down_offset + 2 + c::seed_first_layer_z_idx] - seed_station_z_pos_l0 ) +
        ( yAtZ[2 + c::seed_first_layer_z_idx] - ( c::seed_scale_ycorr + 0.02f ) * ycorr ) *
            cache.dXdYFiber[up_down_dxdy_offset + 1 + c::seed_first_layer_z_idx / 2];

    mydebug( newtx, ( cache.LayerZPos[up_down_offset + 2 + c::seed_first_layer_z_idx] - seed_station_z_pos_l0 ),
             ( yAtZ[2 + c::seed_first_layer_z_idx] - ( c::seed_scale_ycorr + 0.02f ) * ycorr ),
             cache.dXdYFiber[up_down_dxdy_offset + 1 + c::seed_first_layer_z_idx / 2] );

    // search for hits in layer
    start_seed_station_l1    = get_closest_hit_idx( hithandler.hits, start_seed_station_l1, seed_l1_pred );
    auto const seed_delta_l1 = hithandler.hits[start_seed_station_l1] - seed_l1_pred;

    start_seed_station_l2    = get_closest_hit_idx( hithandler.hits, start_seed_station_l2, seed_l2_pred );
    auto const seed_delta_l2 = hithandler.hits[start_seed_station_l2] - seed_l2_pred;

    mydebug( "deltas x1, x2: ", seed_delta_l1, seed_delta_l2, hithandler.IDs[start_seed_station_l1],
             hithandler.IDs[start_seed_station_l2] );

    // calculate acceptance window for uv hits based on momentum
    float const seed_station_uvwindow = ( c::seed_uv_win_off + c::seed_uv_win_slope * absDSlope );

    if ( std::abs( seed_delta_l2 ) < seed_station_uvwindow ) {
      tmpidvec.push_back( start_seed_station_l2 );
      ++numuvhits;
    }

    if ( std::abs( seed_delta_l1 ) < seed_station_uvwindow ) {
      tmpidvec.push_back( start_seed_station_l1 );
      ++numuvhits;
    } else if ( numuvhits < 1 ) {
      // didn't find any uv hits let's try our luck in the next loop iteration
      continue;
    }

    // this variable is saved for later as it is quite good for ghost rejection
    float const seed_uv_delta = numuvhits < 2 ? 1.0f : 1.0f + std::abs( seed_delta_l1 + seed_delta_l2 );

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //                  Start of search in the second station              //
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    mydebug( "starting next station" );

    // these variables will stay true if x-hits are found in both layers
    bool next_station_doublet = true;

    float const next_l0_curve = cache.ext_coeff[c::ext_coeff_start + up_down_offset + 6] * deltaSlope +
                                cache.ext_coeff[c::ext_coeff_start + up_down_offset + 7] * deltaSlope +
                                cache.ext_coeff[c::ext_coeff_start + up_down_offset + 8] * deltaSlope * endv_ty2;

    float const next_station_z_pos_l0 = cache.LayerZPos[up_down_offset + c::next_first_layer_z_idx];
    float const next_l0_pred = firstXhit + newtx * ( next_station_z_pos_l0 - seed_station_z_pos_l0 ) + next_l0_curve;

    float const next_l3_curve = deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 9] +
                                deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 10] +
                                deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 11] * endv_ty2;

    float const next_station_z_pos_l3 = cache.LayerZPos[up_down_offset + c::next_first_layer_z_idx + 3];
    float       next_l3_pred = firstXhit + newtx * ( next_station_z_pos_l3 - seed_station_z_pos_l0 ) + next_l3_curve;

    // get the best hit candidate
    start_next_station_l3    = get_closest_hit_idx( hithandler.hits, start_next_station_l3, next_l3_pred );
    auto const next_delta_l3 = hithandler.hits[start_next_station_l3] - next_l3_pred;

    mydebug( "next l3 pred", next_l3_pred, "curve", next_l3_curve, "delta", next_delta_l3, "idx",
             start_next_station_l3 );

    // did we find a hit?
    if ( std::abs( next_delta_l3 ) < ( c::next_l3_win_off + c::next_l3_win_slope * absDSlope ) ) {
      mydebug( "hit accepted" );

      tmpidvec.push_back( start_next_station_l3 );
      ++numxhits;

      // fill fit vec with the x-hit information
      fillFitVec( fitvec, fitcnt, hithandler.hits[start_next_station_l3], next_station_z_pos_l3 - seed_station_z_pos_l0,
                  next_l3_curve, next_delta_l3 );
    } else {
      // we didn't find a hit and thus also not a doublet
      next_station_doublet = false;

      // do we even really want to continue with this track or is it crap?
      if ( numxhits + numuvhits < 4 ) continue;
    }

    start_next_station_l0    = get_closest_hit_idx( hithandler.hits, start_next_station_l0, next_l0_pred );
    auto const next_delta_l0 = hithandler.hits[start_next_station_l0] - next_l0_pred;

    mydebug( "next l0 pred", next_l0_pred, "curve", next_l0_curve, "delta", next_delta_l0, "idx",
             start_next_station_l0 );

    // did we find a hit?
    if ( std::abs( next_delta_l0 ) < ( c::next_l0_win_off + c::next_l0_win_slope * absDSlope ) ) {
      mydebug( "hit accepted" );
      tmpidvec.push_back( start_next_station_l0 );
      ++numxhits;

      fillFitVec( fitvec, fitcnt, hithandler.hits[start_next_station_l0], next_station_z_pos_l0 - seed_station_z_pos_l0,
                  next_l0_curve, next_delta_l0 );
    } else {
      next_station_doublet = false;
      if ( numxhits + numuvhits < 5 ) continue;
    }

    // we only check for uv hits if we found a doublet
    // this choice was made since a doublet and it's slope allow for small search windows
    // FIXME we could potentially allow this search in other/all scenarios but this need more tweaking
    // Note that this search will only help with ghost rejection, not help with efficiencies
    // Thus as long as ghosts aren't a problem there isn't really a big reason to loosen this requirement
    float next_uv_delta{1.0f};

    if ( next_station_doublet ) {
      mydebug( "next uv search" );

      float const next_tx = ( hithandler.hits[start_next_station_l3] - hithandler.hits[start_next_station_l0] ) /
                            ( next_station_z_pos_l3 - next_station_z_pos_l0 );

      float const next_l1_pred =
          hithandler.hits[start_next_station_l0] +
          next_tx * ( cache.LayerZPos[up_down_offset + c::next_first_layer_z_idx + 1] - next_station_z_pos_l0 ) +
          ( yAtZ[1 + c::next_first_layer_z_idx] - c::next_scale_ycorr * ycorr ) *
              cache.dXdYFiber[up_down_dxdy_offset + c::next_first_layer_z_idx / 2];

      float const next_l2_pred =
          hithandler.hits[start_next_station_l0] +
          next_tx * ( cache.LayerZPos[up_down_offset + c::next_first_layer_z_idx + 2] - next_station_z_pos_l0 ) +
          ( yAtZ[2 + c::next_first_layer_z_idx] - ( c::next_scale_ycorr + 0.02f ) * ycorr ) *
              cache.dXdYFiber[up_down_dxdy_offset + 1 + c::next_first_layer_z_idx / 2];

      start_next_station_l1    = get_closest_hit_idx( hithandler.hits, start_next_station_l1, next_l1_pred );
      auto const next_delta_l1 = hithandler.hits[start_next_station_l1] - next_l1_pred;

      start_next_station_l2    = get_closest_hit_idx( hithandler.hits, start_next_station_l2, next_l2_pred );
      auto const next_delta_l2 = hithandler.hits[start_next_station_l2] - next_l2_pred;

      mydebug( "next l1 pred", next_l1_pred, "delta", next_delta_l1, "idx", start_next_station_l1 );
      mydebug( "next l2 pred", next_l2_pred, "delta", next_delta_l2, "idx", start_next_station_l2 );

      float const S2uvwindow = ( c::next_uv_win_off + c::next_uv_win_slope * absDSlope );

      if ( std::abs( next_delta_l1 ) > S2uvwindow and std::abs( next_delta_l2 ) > S2uvwindow ) { continue; }

      if ( std::abs( next_delta_l1 ) < S2uvwindow and std::abs( next_delta_l2 ) < S2uvwindow ) {
        next_uv_delta += std::abs( next_delta_l1 + next_delta_l2 );
        tmpidvec.push_back( start_next_station_l2 );
        ++numuvhits;
        tmpidvec.push_back( start_next_station_l1 );
        ++numuvhits;
      } else {
        if ( std::abs( next_delta_l2 ) < S2uvwindow ) {
          tmpidvec.push_back( start_next_station_l2 );
          ++numuvhits;
        }

        if ( std::abs( next_delta_l1 ) < S2uvwindow ) {
          tmpidvec.push_back( start_next_station_l1 );
          ++numuvhits;
        }
      }

    } // end of the UV work in the first iteration

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //                      Start of search in Station 1                   //
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    mydebug( "Start Station 1 Work" );

    float const last_l3_curve = deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 3] +
                                deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 4] +
                                deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 5] * endv_ty2;

    // the layer of the seeding hit, thus its z position, depends on the current iteration
    float const last_station_z_pos_l3 = cache.LayerZPos[up_down_offset + c::last_first_layer_z_idx + 3];
    float const last_l3_pred = firstXhit + newtx * ( last_station_z_pos_l3 - seed_station_z_pos_l0 ) + last_l3_curve;

    auto const last_l3_hit_idx =
        last_station_range_l3.first +
        get_closest_hit_idx_binary(
            {hithandler.hits.data() + last_station_range_l3.first, last_station_range_l3.second}, last_l3_pred );
    auto const last_delta_l3 = hithandler.hits[last_l3_hit_idx] - last_l3_pred;

    mydebug( "last_straight", firstXhit + newtx * ( last_station_z_pos_l3 - seed_station_z_pos_l0 ), "last_l3_curve",
             last_l3_curve, "pred", last_l3_pred, "delta", last_delta_l3, "idx", last_l3_hit_idx );

    bool last_station_doublet = true;

    if ( std::abs( last_delta_l3 ) < c::last_l3_win_off + c::last_l3_win_slope * absDSlope ) {
      tmpidvec.push_back( last_l3_hit_idx );
      ++numxhits;

      // the delta z wrt the seeding layer depends on the current iteration
      fillFitVec( fitvec, fitcnt, hithandler.hits[last_l3_hit_idx], last_station_z_pos_l3 - seed_station_z_pos_l0,
                  last_l3_curve, last_delta_l3 );
      mydebug( "Hit Selected" );
    } else {
      last_station_doublet = false;
      if ( numxhits < 4 or numuvhits < 3 ) continue;
    }

    float const last_l0_curve = deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 0] +
                                deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 1] +
                                deltaSlope * cache.ext_coeff[c::ext_coeff_start + up_down_offset + 2] * endv_ty2;

    // the layer of the seeding hit, thus its z position, depends on the current iteration
    float const last_station_z_pos_l0 = cache.LayerZPos[up_down_offset + c::last_first_layer_z_idx];
    float const last_l0_pred = firstXhit + newtx * ( last_station_z_pos_l0 - seed_station_z_pos_l0 ) + last_l0_curve;

    auto const last_l0_hit_idx =
        last_station_range_l0.first +
        get_closest_hit_idx_binary(
            {hithandler.hits.data() + last_station_range_l0.first, last_station_range_l0.second}, last_l0_pred );
    auto const last_delta_l0 = hithandler.hits[last_l0_hit_idx] - last_l0_pred;

    mydebug( "last_straight", firstXhit + newtx * ( last_station_z_pos_l0 - seed_station_z_pos_l0 ), "last_l0_curve",
             last_l0_curve, "pred", last_l0_pred, "delta", last_delta_l0, "idx", last_l0_hit_idx );

    if ( std::abs( last_delta_l0 ) < c::last_l0_win_off + c::last_l0_win_slope * absDSlope ) {
      tmpidvec.push_back( last_l0_hit_idx );
      ++numxhits;

      // the delta z wrt the seeding layer
      fillFitVec( fitvec, fitcnt, hithandler.hits[last_l0_hit_idx], last_station_z_pos_l0 - seed_station_z_pos_l0,
                  last_l0_curve, last_delta_l0 );

      mydebug( "Hit Selected" );
    } else {
      last_station_doublet = false;
      if ( numxhits < ( 5 ) or numuvhits < 3 ) continue;
    }

    float last_uv_delta = 1.0f;

    if ( last_station_doublet ) {
      mydebug( "last uv search" );
      float const last_tx = ( hithandler.hits[last_l3_hit_idx] - hithandler.hits[last_l0_hit_idx] ) /
                            ( last_station_z_pos_l3 - last_station_z_pos_l0 );

      float const last_l1_pred =
          hithandler.hits[last_l0_hit_idx] +
          last_tx * ( cache.LayerZPos[up_down_offset + c::last_first_layer_z_idx + 1] - last_station_z_pos_l0 ) +
          ( yAtZ[1 + c::last_first_layer_z_idx] - c::last_scale_ycorr * ycorr ) *
              cache.dXdYFiber[up_down_dxdy_offset + c::last_first_layer_z_idx / 2];

      float const last_l2_pred =
          hithandler.hits[last_l0_hit_idx] +
          last_tx * ( cache.LayerZPos[up_down_offset + c::last_first_layer_z_idx + 2] - last_station_z_pos_l0 ) +
          ( yAtZ[2 + c::last_first_layer_z_idx] - ( c::last_scale_ycorr + 0.02f ) * ycorr ) *
              cache.dXdYFiber[up_down_dxdy_offset + 1 + c::last_first_layer_z_idx / 2];

      auto const last_l1_hit_idx =
          last_station_range_l1.first +
          get_closest_hit_idx_binary(
              {hithandler.hits.data() + last_station_range_l1.first, last_station_range_l1.second}, last_l1_pred );
      auto const last_delta_l1 = hithandler.hits[last_l1_hit_idx] - last_l1_pred;

      auto const last_l2_hit_idx =
          last_station_range_l2.first +
          get_closest_hit_idx_binary(
              {hithandler.hits.data() + last_station_range_l2.first, last_station_range_l2.second}, last_l2_pred );
      auto const last_delta_l2 = hithandler.hits[last_l2_hit_idx] - last_l2_pred;

      // the UV window value depends on the current iteration
      float const last_uvwindow = c::last_uv_win_off + c::last_uv_win_slope * absDSlope;

      mydebug( "last l1 pred", last_l1_pred, "delta", last_delta_l1, "idx", last_l1_hit_idx );
      mydebug( "last l2 pred", last_l2_pred, "delta", last_delta_l2, "idx", last_l2_hit_idx );

      if ( std::abs( last_delta_l1 ) > last_uvwindow and std::abs( last_delta_l2 ) > last_uvwindow ) {
        mydebug( "none found;skipping" );
        continue;
      }
      if ( std::abs( last_delta_l1 ) < last_uvwindow and std::abs( last_delta_l2 ) < last_uvwindow ) {
        last_uv_delta += std::abs( last_delta_l1 + last_delta_l2 );
        tmpidvec.push_back( last_l2_hit_idx );
        ++numuvhits;
        tmpidvec.push_back( last_l1_hit_idx );
        ++numuvhits;
      } else {

        if ( std::abs( last_delta_l2 ) < last_uvwindow ) {
          tmpidvec.push_back( last_l2_hit_idx );
          ++numuvhits;
        }

        if ( std::abs( last_delta_l1 ) < last_uvwindow ) {
          tmpidvec.push_back( last_l1_hit_idx );
          ++numuvhits;
        }
      }
    }

    mydebug( "End Station 1 Work" );

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    //                      Start of finalization procedure                //
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    mydebug( "+++++++++Trackcandidate++++++++++++++++" );

    // mydebug( fitvec, fitcnt, newtx, firstXhit );
    auto [finaltx, newx0, chi2] = simple_fit( fitvec, fitcnt, newtx, firstXhit );
    float const candPQ          = calc_mom<std::is_same_v<c, config_loop_1>>( finaltx, endv_tx, endv_ty2,
                                                                     std::copysign( 1.f, endv_qp ) * factor, factor );
    auto const  quality =
        calc_quality( endv_pt, factor, firstXhit, straighExt, endv_pq, endv_qp, candPQ, p_wrongSignPT.value(), chi2,
                      endv_txy, seed_uv_delta, next_uv_delta, last_uv_delta, numuvhits, numxhits );

    if ( quality < bestcandidate.quality )
      bestcandidate = {tmpidvec, candPQ, finaltx, newx0, quality, numuvhits + numxhits};

    mydebug( candPQ, finaltx, endv_pq );
    mydebug( seed_uv_delta, next_uv_delta, last_uv_delta );
    mydebug( "chi2", chi2, "pt", std::abs( candPQ ) * endv_txy / std::sqrt( 1 + endv_txy2 ) );
    mydebug( "final quality", quality, "worst candidate", bestcandidate.quality );
  } // loop over possible station three doublets
  return bestcandidate;
}
