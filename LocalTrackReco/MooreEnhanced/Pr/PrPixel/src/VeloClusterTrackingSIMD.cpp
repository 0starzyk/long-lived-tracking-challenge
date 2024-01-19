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

#include <algorithm>
#include <array>
#include <tuple>
#include <vector>

// Gaudi
#include "LHCbAlgs/Transformer.h"

// LHCb
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "DetDesc/IConditionDerivationMgr.h"
#include "Event/PrHits.h"
#include "Event/PrVeloTracks.h"
#include "Event/RawEvent.h"
#include "Event/StateParameters.h"
#include "Kernel/LHCbID.h"
#include "Kernel/STLExtensions.h"
#include "Kernel/VPConstants.h"
#include "VPDAQ/VPRetinaClusterConstants.h"
#include "VPDet/DeVP.h"
#include "VPDet/VPDetPaths.h"
#include "VPKernel/PixelUtils.h"

// Local
#include "VPClusCache.h"
#include "VeloKalmanHelpers.h"

#include "GaudiKernel/StdArrayAsProperty.h"

/**
 * This code is used for the tracking, and clustering
 * of the (pixel) VELO detector. The output tracks
 * *have* been fitted using the simple, VELO-only
 * (pT = 400MeV) fit set in VeloKalmanHelpers
 *
 * The code effectively:
 * - prepares the clusters in a way that is accessible
 *   to the user, but also prepares the clusters
 *   in the format used by the tracking.
 *   This functionality can be found in the
 *   getHits functions in the anon namespaces
 * - uses these to form tracks, as done in the
 *   code in the ClusterTrackingSIMD class
 * - immediately also fits these tracks, as part of
 *   the same function (in ClusterTrackingSIMD)
 *   but relying on calls to VeloKalmanHelpers
 *
 * Information on the overall algorithm can be found
 * at https://arxiv.org/abs/1912.09901
 *
 * This algorithm is intended to be the "default"
 * in the reconstruction of VELO tracks for
 * HLT2 in Run 3.
 *
 * Further details regarding the implementation:
 *  - Sensor unbiasing:
 *    for the 'sensor unbiasing', it was tried to
 *    make the boolean-checks optional through
 *    a bool template argument, but this ran
 *    into problems with the specialization
 *
 **/
using TracksTag = LHCb::Pr::Velo::Tag;

namespace LHCb::Pr::Velo {

  namespace VPInfos {
    constexpr int NPlanes           = 26;
    constexpr int NModulesPerPlane  = 2;
    constexpr int NModules          = NPlanes * NModulesPerPlane;
    constexpr int NSensorsPerModule = 4;
    constexpr int NSensors          = NModules * NSensorsPerModule;
    constexpr int NSensorsPerPlane  = NModulesPerPlane * NSensorsPerModule;
    constexpr int NChipsPerSensor   = 3;
    constexpr int NRows             = 256;
    constexpr int NColumns          = 256;
  } // namespace VPInfos

  enum class SearchMode { Default, Fast, Full };

  template <SearchMode mode>
  struct SearchModeParam {
    static constexpr int planes               = VPInfos::NPlanes;
    static constexpr int planes_stride        = VPInfos::NSensorsPerPlane;
    static constexpr int skip_forward_default = 1;
    static constexpr int plane_window         = 3; // how many planes to keep in the sliding window
  };

  template <>
  struct SearchModeParam<SearchMode::Full> {
    static constexpr int planes               = VPInfos::NModules;
    static constexpr int planes_stride        = VPInfos::NSensorsPerModule;
    static constexpr int skip_forward_default = 4;
    static constexpr int plane_window         = 6;
  };

  namespace {
    //=============================================================================
    // Internal data structures:
    //=============================================================================

    constexpr static int reserved_tracks_size = 2048;
    constexpr static int max_hits             = LHCb::Pr::TracksInfo::MaxVPHits;

    namespace LightTrackTag {
      struct p0 : Event::vec3_field {}; // vector
      struct p1 : Event::vec3_field {};
      struct p2 : Event::vec3_field {};
      struct nHits : Event::int_field {};
      struct skipped : Event::int_field {};
      struct sumScatter : Event::float_field {};
      struct hit : Event::ints_field<max_hits> {};

      template <typename T>
      using light_t = Event::SOACollection<T, p0, p1, p2, nHits, skipped, sumScatter, hit>;
    } // namespace LightTrackTag

    struct LightTracksSoA : LightTrackTag::light_t<LightTracksSoA> {
      using base_t = typename LightTrackTag::light_t<LightTracksSoA>;
      using base_t::base_t;
    };

    constexpr static int reserved_hits_inPlane = 1024;
    namespace HitsTag {
      struct pos : Event::vec3_field {}; // vector
      struct phi : Event::float_field {};
      struct idx : Event::int_field {};

      template <typename T>
      using hits_t = Event::SOACollection<T, pos, phi, idx>;
    } // namespace HitsTag

    struct HitsPlane : HitsTag::hits_t<HitsPlane> {
      using base_t = typename HitsTag::hits_t<HitsPlane>;
      using base_t::base_t;
    };
    template <typename T, typename M>
    [[gnu::always_inline]] inline void removeHits( HitsPlane* hits, M mask, T indices ) {
      auto hit      = hits->scalar();
      int  new_size = hits->size();
      for ( int i = 0; i < new_size; i++ ) {
        if ( !none( mask && ( T( hit[i].get<HitsTag::idx>().cast() ) == indices ) ) ) {
          new_size--;
          hit[i].field<HitsTag::pos>().set( hit[new_size].get<HitsTag::pos>().vec3() );
          hit[i].field<HitsTag::phi>().set( hit[new_size].get<HitsTag::phi>() );
          hit[i].field<HitsTag::idx>().set( hit[new_size].get<HitsTag::idx>() );
          i--;
        }
      }
      hits->resize( new_size );
    }
    //=============================================================================
    // Union-Find functions used in clustering:
    //=============================================================================
    [[gnu::always_inline]] inline int find( int* L, int i ) {
      int ai = L[i];
      while ( ai != L[ai] ) ai = L[ai];
      return ai;
    }

    [[gnu::always_inline]] inline int merge( int* L, int ai, int j ) { // Union
      int aj = find( L, j );
      if ( ai < aj )
        L[aj] = ai;
      else {
        L[ai] = aj;
        ai    = aj;
      }
      return ai;
    }

    //=============================================================================
    // Check version
    //=============================================================================

    template <RawBank::BankType VP_banktype>
    bool check_version( const unsigned int version );

    template <>
    inline __attribute__( ( always_inline ) ) bool check_version<RawBank::VP>( const unsigned int version ) {
      if ( version != 2 && version != 4 ) { return false; }
      return true;
    }

    template <>
    inline __attribute__( ( always_inline ) ) bool
    check_version<RawBank::VPRetinaCluster>( const unsigned int version ) {
      if ( version != 2 && version != 3 && version != 4 ) { return false; }
      return true;
    }

    //=============================================================================
    // Sort super pixels by column (major) and row (minor)
    //=============================================================================
    auto SPLowerThan = []( unsigned int lhs, unsigned int rhs ) { return ( lhs & 0x7FFF00 ) < ( rhs & 0x7FFF00 ); };
    auto SPEqual     = []( unsigned int lhs, unsigned int rhs ) { return ( lhs & 0x7FFF00 ) == ( rhs & 0x7FFF00 ); };

    //=============================================================================
    // Sort Banks
    //=============================================================================

    template <RawBank::BankType VP_banktype>
    void sortBanks( const RawBank::View& rawBanks_in, RawEvent& rawEvent_out );

    template <>
    [[gnu::always_inline]] inline void sortBanks<RawBank::VP>( const RawBank::View& rawBanks_in,
                                                               RawEvent&            rawEvent_out ) {
      for ( auto iterBank : rawBanks_in ) {
        if ( iterBank->version() > 3 ) {
          const uint32_t        sensor0 = ( ( iterBank->sourceID() ) & 0x1FFU ) << 1;
          const uint32_t        sensor1 = sensor0 + 1;
          std::vector<uint32_t> data0;
          data0.reserve( ( iterBank->range<uint32_t>() ).size() );
          std::vector<uint32_t> data1;
          data1.reserve( ( iterBank->range<uint32_t>() ).size() );
          for ( auto word : iterBank->range<uint32_t>() ) {
            if ( ( ( word >> VPRetinaCluster::SPsensorID_shift ) & 0x1U ) ) { // check if SP belongs to sensor1
              data1.push_back( word );
            } else {
              data0.push_back( word );
            }
          }

          // sort super pixels column major on each sensor
          std::sort( data0.begin(), data0.end(), SPLowerThan );
          std::sort( data1.begin(), data1.end(), SPLowerThan );

          // Remove duplicate super-pixels
          data0.erase( std::unique( data0.begin(), data0.end(), SPEqual ), data0.end() );
          data1.erase( std::unique( data1.begin(), data1.end(), SPEqual ), data1.end() );

          rawEvent_out.addBank( sensor0, LHCb::RawBank::VP, VPRetinaCluster::c_SPBankVersion, data0 );
          rawEvent_out.addBank( sensor1, LHCb::RawBank::VP, VPRetinaCluster::c_SPBankVersion, data1 );
        } else {
          const uint32_t        sensor = iterBank->sourceID();
          std::vector<uint32_t> data;
          data.reserve( ( iterBank->range<uint32_t>() ).size() );
          for ( auto word : iterBank->range<uint32_t>() ) { data.push_back( word ); }

          // sort super pixels column major on each sensor
          if ( data.size() > 1 ) { std::sort( ++data.begin(), data.end(), SPLowerThan ); }

          rawEvent_out.addBank( sensor, LHCb::RawBank::VP, 2, data );
        }
      }
    }

    template <>
    [[gnu::always_inline]] inline void sortBanks<RawBank::VPRetinaCluster>( const RawBank::View& rawBanks_in,
                                                                            RawEvent&            rawEvent_out ) {
      for ( auto iterBank : rawBanks_in ) {
        if ( iterBank->version() > 3 ) {
          const uint32_t        sensor0 = ( ( iterBank->sourceID() ) & 0x1FFU ) << 1;
          const uint32_t        sensor1 = sensor0 + 1;
          std::vector<uint32_t> data0;
          data0.reserve( ( iterBank->range<uint32_t>() ).size() );
          std::vector<uint32_t> data1;
          data1.reserve( ( iterBank->range<uint32_t>() ).size() );
          for ( auto word : iterBank->range<uint32_t>() ) {
            if ( ( ( word >> VPRetinaCluster::sensorID_shift ) & 0x1U ) ) { // check if cluster belongs to sensor1
              data1.push_back( word );
            } else {
              data0.push_back( word );
            }
          }
          rawEvent_out.addBank( sensor0, LHCb::RawBank::VPRetinaCluster, VPRetinaCluster::c_bankVersion, data0 );
          rawEvent_out.addBank( sensor1, LHCb::RawBank::VPRetinaCluster, VPRetinaCluster::c_bankVersion, data1 );
        } else {
          const uint32_t        sensor = iterBank->sourceID();
          std::vector<uint32_t> data;
          data.reserve( ( iterBank->range<uint32_t>() ).size() );
          for ( auto word : iterBank->range<uint32_t>() ) { data.push_back( word ); }

          rawEvent_out.addBank( sensor, LHCb::RawBank::VPRetinaCluster, iterBank->version(), data );
        }
      }
    }

    //=============================================================================
    // Clustering Algorithm
    //=============================================================================

    template <RawBank::BankType VP_banktype>
    void getHits( const LHCb::RawBank** VPRawBanks, const int sensor0, const int sensor1, const DeVP& vp,
                  HitsPlane& Pout, VP::Hits& hits, const std::array<bool, VPInfos::NSensors>& masked_sensors );

    template <>
    [[gnu::always_inline]] inline void
    getHits<RawBank::VP>( const LHCb::RawBank** VPRawBanks, // List of input Super-Pixels
                          const int sensor0, const int sensor1, const DeVP& vp, HitsPlane& Pout, VP::Hits& hits,
                          const std::array<bool, VPInfos::NSensors>& masked_sensors ) {

      // Clustering buffers
      constexpr int MAX_CLUSTERS_PER_SENSOR = 1024;    // Max number of SP per bank
      uint32_t      pixel_SP[MAX_CLUSTERS_PER_SENSOR]; // SP to process
      int           pixel_L[MAX_CLUSTERS_PER_SENSOR];  // Label equivalences
      uint32_t      pixel_SX[MAX_CLUSTERS_PER_SENSOR]; // x sum of clusters' pixels
      uint32_t      pixel_SY[MAX_CLUSTERS_PER_SENSOR]; // y sum of clusters' pixels
      uint32_t      pixel_S[MAX_CLUSTERS_PER_SENSOR];  // number of clusters' pixels
      Pout.resize( 0 );
      int n_hits = 0; // Number of clusters after filtering
      int offset = hits.size();
      for ( int s = sensor0; s < sensor1; s++ ) {
        if ( VPRawBanks[s] == nullptr ) { continue; }
        int bank_version = VPRawBanks[s]->version();

        const uint32_t* bank = VPRawBanks[s]->data();
        uint32_t        nsp;
        if ( bank_version > 3 ) {
          nsp = ( VPRawBanks[s]->range<uint32_t>() ).size();
        } else {
          nsp = *bank++;
        }

        int pixel_N = 0; // Pixels to cluster count
        int labels  = 0; // Total number of generated clusters

        for ( unsigned int i = 0; i < nsp; ++i ) {
          const uint32_t sp_word = *bank++;
          uint8_t        sp      = SP_getPixels( sp_word );

          if ( 0 == sp ) continue; // protect against zero super pixels.

          const auto sp_row = SP_getRow( sp_word );
          const auto sp_col = SP_getCol( sp_word );

          if ( sp_row > ( VPInfos::NRows / 4 - 1 ) )
            continue; // protect against super pixels outside sensor coordinates.
          if ( sp_col > ( VPInfos::NColumns * VPInfos::NChipsPerSensor / 2 - 1 ) )
            continue; // protect against super pixels outside sensor coordinates.

          // This SP is isolated, skip clustering :
          if ( SP_isIsolated( sp_word ) ) {
            uint8_t mask = s_SPMasks[sp];

            auto n_kx_ky = s_SPn_kx_ky[sp & mask];
            auto n       = SPn_kx_ky_getN( n_kx_ky );  // number of pixel in this sp
            auto kx      = SPn_kx_ky_getKx( n_kx_ky ); // sum of x in this sp
            auto ky      = SPn_kx_ky_getKy( n_kx_ky ); // sum of y in this sp

            pixel_SX[labels] = sp_col * n * 2 + kx;
            pixel_SY[labels] = sp_row * n * 4 + ky;
            pixel_S[labels]  = n;
            labels++;

            if ( mask != 0xFF ) { // Add 2nd cluster
              n_kx_ky = s_SPn_kx_ky[sp & ( ~mask )];
              n       = SPn_kx_ky_getN( n_kx_ky );  // number of pixel in this sp
              kx      = SPn_kx_ky_getKx( n_kx_ky ); // sum of x in this sp
              ky      = SPn_kx_ky_getKy( n_kx_ky ); // sum of y in this sp

              pixel_SX[labels] = sp_col * n * 2 + kx;
              pixel_SY[labels] = sp_row * n * 4 + ky;
              pixel_S[labels]  = n;
              labels++;
            }

            continue;
          }

          // This one is not isolated, add it to clustering :
          uint32_t mask     = 0x7FFFFF00 | s_SPMasks[sp];
          pixel_SP[pixel_N] = sp_word & mask;
          pixel_L[pixel_N]  = pixel_N;
          pixel_N++;

          if ( mask != 0x7FFFFFFF ) {                      // Add 2nd cluster
            pixel_SP[pixel_N] = sp_word & ( mask ^ 0xFF ); // ~ of low 8 bits
            pixel_L[pixel_N]  = pixel_N;
            pixel_N++;
          }
        } // loop over super pixels in raw bank

        // SparseCCL: Connected Components Labeling and Analysis for sparse images
        // https://hal.archives-ouvertes.fr/hal-02343598/document
        // (This version assume SP are ordered by col then row)
        int start_j = 0;
        for ( int i = 0; i < pixel_N; i++ ) { // Pixel Scan
          uint32_t sp_word_i, sp_word_j;
          uint32_t ai;
          uint8_t  spi, spj;
          int      x_i, x_j, y_i, y_j;

          sp_word_i = pixel_SP[i];
          spi       = SP_getPixels( sp_word_i );
          y_i       = SP_getRow( sp_word_i );
          x_i       = SP_getCol( sp_word_i );

          pixel_L[i] = i;
          ai         = i;

          for ( int j = start_j; j < i; j++ ) {
            sp_word_j = pixel_SP[j];
            spj       = SP_getPixels( sp_word_j );
            y_j       = SP_getRow( sp_word_j );
            x_j       = SP_getCol( sp_word_j );
            if ( is_adjacent_8C_SP( x_i, y_i, x_j, y_j, spi, spj ) ) {
              ai = merge( pixel_L, ai, j );
            } else if ( x_j + 1 < x_i )
              start_j++;
          }
        }

        for ( int i = 0; i < pixel_N; i++ ) { // Transitive Closure + CCA
          uint32_t sp_word = pixel_SP[i];
          uint8_t  sp      = SP_getPixels( sp_word );
          uint32_t y_i     = SP_getRow( sp_word );
          uint32_t x_i     = SP_getCol( sp_word );

          auto n_kx_ky = s_SPn_kx_ky[sp];
          auto n       = SPn_kx_ky_getN( n_kx_ky );                // number of pixel in this sp
          auto kx      = x_i * n * 2 + SPn_kx_ky_getKx( n_kx_ky ); // sum of x in this sp
          auto ky      = y_i * n * 4 + SPn_kx_ky_getKy( n_kx_ky ); // sum of y in this sp

          uint32_t l;
          if ( pixel_L[i] == i ) {
            l           = labels++; // new label
            pixel_SX[l] = kx;
            pixel_SY[l] = ky;
            pixel_S[l]  = n;
          } else {
            l = pixel_L[pixel_L[i]]; // transitive closure
            pixel_SX[l] += kx;
            pixel_SY[l] += ky;
            pixel_S[l] += n;
          }
          pixel_L[i] = l;
        }

        const auto sensorID = VPRawBanks[s]->sourceID();
        const auto sensor   = LHCb::Detector::VPChannelID::SensorID( sensorID );
        auto       ltg      = vp.ltg( sensor );
        for ( int i = 0; i < labels; i++ ) {
          uint32_t n = pixel_S[i];

          uint32_t x = pixel_SX[i];
          uint32_t y = pixel_SY[i];

          const uint32_t cx   = x / n;
          const auto     cy   = LHCb::Detector::VPChannelID::RowID{y / n};
          const auto     scol = LHCb::Detector::VPChannelID::ScolID{cx};

          // store target (3D point for tracking)
          const float fx      = x / static_cast<float>( n ) - cx;
          const float fy      = y / static_cast<float>( n );
          const float local_x = vp.local_x( cx ) + fx * vp.x_pitch( cx );
          const float local_y = ( 0.5f + fy ) * vp.pixel_size();

          const float gx = ( ltg[0] * local_x + ltg[1] * local_y + ltg[9] );
          const float gy = ( ltg[3] * local_x + ltg[4] * local_y + ltg[10] );
          const float gz = ( ltg[6] * local_x + ltg[7] * local_y + ltg[11] );

          auto hit = hits.emplace_back<SIMDWrapper::InstructionSet::Scalar>();
          hit.field<VP::VPHitsTag::pos>().set( gx, gy, gz );
          hit.field<VP::VPHitsTag::ChannelId>().set(
              SIMDWrapper::scalar::int_v( LHCb::Detector::VPChannelID{sensor, scol, cy}.channelID() ) );

          if ( masked_sensors[sensorID] ) { continue; }

          auto pout = Pout.emplace_back<SIMDWrapper::InstructionSet::Scalar>();
          pout.field<HitsTag::pos>().set( gx, gy, gz );
          pout.field<HitsTag::idx>().set( offset + n_hits );
          n_hits++;
        }
      } // Loop over sensors

      // Pre-compute phi
      for ( auto pout : Pout.simd() ) {
        auto pos = pout.get<HitsTag::pos>().vec3();
        pout.field<HitsTag::phi>().set( pos.phi() );
      }
    }

    template <>
    [[gnu::always_inline]] inline void
    getHits<RawBank::VPRetinaCluster>( const LHCb::RawBank** VPRawBanks, // List of input Super-Pixels
                                       const int sensor0, const int sensor1, const DeVP& vp, HitsPlane& Pout,
                                       VP::Hits& hits, const std::array<bool, VPInfos::NSensors>& masked_sensors ) {

      Pout.resize( 0 );
      int n_hits = 0; // Number of clusters after filtering
      int offset = hits.size();
      for ( int s = sensor0; s < sensor1; s++ ) {

        if ( VPRawBanks[s] == nullptr ) { continue; }
        int bank_version = VPRawBanks[s]->version();

        const uint32_t* bank = VPRawBanks[s]->data();
        uint32_t        nrc;
        if ( bank_version > 3 ) {
          nrc = ( VPRawBanks[s]->range<uint32_t>() ).size();
        } else {
          nrc = *bank++;
        }
        const auto sensorID = VPRawBanks[s]->sourceID();
        const auto sensor   = LHCb::Detector::VPChannelID::SensorID( sensorID );

        auto ltg = vp.ltg( sensor );
        for ( unsigned int i = 0; i < nrc; ++i ) {
          const uint32_t rc_word = *bank++;

          if ( 0 == rc_word ) continue; // protect against zero clusters.

          uint32_t                            cx, cx_frac_half, cx_frac_quarter, cy_frac_half, cy_frac_quarter;
          float                               fx, fy;
          LHCb::Detector::VPChannelID::RowID  cy;
          LHCb::Detector::VPChannelID::OrfxID or_fx;
          LHCb::Detector::VPChannelID::OrfyID or_fy;
          LHCb::Detector::VPChannelID::ScolID scol;

          if ( bank_version == 2 ) {
            cx   = ( rc_word >> 14 ) & 0x3FF;
            fx   = ( ( rc_word >> 11 ) & 0x7 ) / 8.f;
            cy   = LHCb::Detector::VPChannelID::RowID{( rc_word >> 3 ) & 0xFF};
            fy   = ( rc_word & 0x7FF ) / 8.f;
            scol = LHCb::Detector::VPChannelID::ScolID{cx};

            or_fx = LHCb::Detector::VPChannelID::OrfxID{0};
            or_fy = LHCb::Detector::VPChannelID::OrfyID{0};
          } else {
            cx   = ( rc_word >> 12 ) & 0x3FF;
            fx   = ( ( rc_word >> 10 ) & 0x3 ) / 4.f;
            cy   = LHCb::Detector::VPChannelID::RowID{( rc_word >> 2 ) & 0xFF};
            fy   = ( rc_word & 0x3FF ) / 4.f;
            scol = LHCb::Detector::VPChannelID::ScolID{cx};

            cx_frac_half    = ( rc_word >> 11 ) & 0x1;
            cx_frac_quarter = ( rc_word >> 10 ) & 0x1;
            or_fx           = LHCb::Detector::VPChannelID::OrfxID{( cx_frac_half | cx_frac_quarter )};

            cy_frac_half    = ( rc_word >> 1 ) & 0x1;
            cy_frac_quarter = (rc_word)&0x1;
            or_fy           = LHCb::Detector::VPChannelID::OrfyID{( cy_frac_half | cy_frac_quarter )};
          }

          if ( to_unsigned( cy ) > ( VPInfos::NRows - 1 ) )
            continue; // protect against super pixels outside sensor coordinates.
          if ( cx > ( VPInfos::NColumns * VPInfos::NChipsPerSensor - 1 ) )
            continue; // protect against super pixels outside sensor coordinates.

          const float local_x = vp.local_x( cx ) + fx * vp.x_pitch( cx );
          const float local_y = ( 0.5f + fy ) * vp.pixel_size();

          const float gx = ( ltg[0] * local_x + ltg[1] * local_y + ltg[9] );
          const float gy = ( ltg[3] * local_x + ltg[4] * local_y + ltg[10] );
          const float gz = ( ltg[6] * local_x + ltg[7] * local_y + ltg[11] );

          auto hit = hits.emplace_back<SIMDWrapper::InstructionSet::Scalar>();
          hit.field<VP::VPHitsTag::pos>().set( gx, gy, gz );
          hit.field<VP::VPHitsTag::ChannelId>().set(
              SIMDWrapper::scalar::int_v( LHCb::Detector::VPChannelID{sensor, scol, cy, or_fx, or_fy}.channelID() ) );

          if ( masked_sensors[sensorID] ) { continue; }

          auto pout = Pout.emplace_back<SIMDWrapper::InstructionSet::Scalar>();
          pout.field<HitsTag::pos>().set( gx, gy, gz );
          pout.field<HitsTag::idx>().set( offset + n_hits );

          n_hits++;
        }
      } // Loop over sensors

      // Pre-compute phi
      for ( auto pout : Pout.simd() ) {
        auto pos = pout.get<HitsTag::pos>().vec3();
        pout.field<HitsTag::phi>().set( pos.phi() );
      }
    }

    // ===========================================================================
    // Tracking algorithm
    // ===========================================================================

    template <int N, typename F, typename I>
    [[gnu::always_inline]] inline void closestsHitsInPhi( HitsPlane* P0, F phi1, Vec3<F> p0_pos[N], I p0_idx[N] ) {

      F distances[N];

      // fill the first hits
      const auto hitP0 = P0->scalar();
      for ( int i = 0; i < N; i++ ) {
        auto sPos = hitP0[i].get<HitsTag::pos>().vec3();
        auto sIdx = hitP0[i].get<HitsTag::idx>();
        auto sPhi = hitP0[i].get<HitsTag::phi>();

        p0_pos[i]    = Vec3<F>( sPos.x, sPos.y, sPos.z );
        p0_idx[i]    = sIdx.cast();
        distances[i] = abs( phi1 - F( sPhi ) );
      }

      for ( int i = N; i < static_cast<int>( P0->size() ); i++ ) {
        auto sPos = hitP0[i].get<HitsTag::pos>().vec3();
        auto sIdx = hitP0[i].get<HitsTag::idx>();
        auto sPhi = hitP0[i].get<HitsTag::phi>();

        Vec3<F> pos = Vec3<F>( sPos.x, sPos.y, sPos.z );
        I       idx = sIdx.cast();
        F       dis = abs( phi1 - F( sPhi ) );

        for ( int j = 0; j < N; j++ ) {
          auto mask = dis < distances[j];

          swap( mask, distances[j], dis );
          swap( mask, p0_pos[j].x, pos.x );
          swap( mask, p0_pos[j].y, pos.y );
          swap( mask, p0_pos[j].z, pos.z );
          swap( mask, p0_idx[j], idx );
        }
      }
    }
  } // namespace

  /**
   * @class ClusterTrackingSIMD
   * Clustering and Velo Tracking algorithm using SIMD
   * https://arxiv.org/abs/1912.09901
   *
   * @author Arthur Hennequin (CERN, LIP6)
   */
  template <RawBank::BankType VP_banktype, SearchMode searchMode>
  class ClusterTrackingSIMD
      : public LHCb::Algorithm::MultiTransformer<std::tuple<VP::Hits, Tracks, Tracks>(
                                                     const EventContext&, const RawBank::View&, const DeVP& ),
                                                 LHCb::DetDesc::usesConditions<DeVP>> {

  public:
    //=============================================================================
    // Standard constructor, initializes variables
    //=============================================================================
    ClusterTrackingSIMD( const std::string& name, ISvcLocator* pSvcLocator )
        : MultiTransformer(
              name, pSvcLocator, {KeyValue{"RawBanks", "DAQ/RawBanks/VP"}, KeyValue{"DEVP", LHCb::Det::VP::det_path}},
              {KeyValue{"HitsLocation", "Raw/VP/Hits"}, KeyValue{"TracksBackwardLocation", "Rec/Track/VeloBackward"},
               KeyValue{"TracksLocation", "Rec/Track/Velo"}} ) {}

    [[gnu::flatten]] void TrackSeeding( HitsPlane* P0, HitsPlane* P1, HitsPlane* P2, LightTracksSoA* tracks ) const {
      using simd = SIMDWrapper::best::types;
      using I    = simd::int_v;
      using F    = simd::float_v;
      using sF   = SIMDWrapper::scalar::float_v;

      const int N_CANDIDATES = 3;
      Vec3<F>   p0_candidates_pos[N_CANDIDATES];
      I         p0_candidates_idx[N_CANDIDATES];

      if ( P0->size() == 0 || P1->size() == 0 || P2->size() == 0 ) return;
      int        P1size = P1->size();
      const auto hitP1  = P1->simd();
      P1->resize( 0 );
      for ( int h1 = 0; h1 < P1size; h1 += simd::size ) {
        auto loop_mask = simd::loop_mask( h1, P1size );

        const Vec3<F> p1   = hitP1[h1].get<HitsTag::pos>().vec3();
        const F       phi1 = hitP1[h1].get<HitsTag::phi>();
        const I       vh1  = hitP1[h1].get<HitsTag::idx>();

        F bestFit = m_max_scatter_seeding.value();
        I vh0     = 0;
        I vh2     = 0;

        Vec3<F> bestP0 = Vec3<F>( 0, 0, 0 );
        Vec3<F> bestP2 = Vec3<F>( 0, 0, 0 );

        closestsHitsInPhi<N_CANDIDATES>( P0, phi1, p0_candidates_pos, p0_candidates_idx );

        for ( const auto& hitP2 : P2->scalar() ) {
          const I        cid2       = hitP2.get<HitsTag::idx>().cast();
          const Vec3<sF> sp2        = hitP2.get<HitsTag::pos>().vec3();
          const Vec3<F>  p2         = Vec3<F>( sp2.x, sp2.y, sp2.z );
          auto           m_found_h0 = simd::mask_false();

          for ( int i = 0; i < N_CANDIDATES; i++ ) {
            if ( i >= static_cast<int>( P0->size() ) ) break;

            const Vec3<F> p0 = p0_candidates_pos[i];
            const I       h0 = p0_candidates_idx[i];

            const F vr = ( p2.z - p0.z ) / ( p1.z - p0.z );
            const F x  = ( p1.x - p0.x ) * vr + p0.x;
            const F y  = ( p1.y - p0.y ) * vr + p0.y;

            const F dx      = x - p2.x;
            const F dy      = y - p2.y;
            const F scatter = dx * dx + dy * dy;

            auto m_bestfit = loop_mask && ( scatter < bestFit );

            if ( none( m_bestfit ) ) continue;

            m_found_h0 = m_found_h0 || m_bestfit;

            bestFit  = select( m_bestfit, scatter, bestFit );
            bestP0.x = select( m_bestfit, p0.x, bestP0.x );
            bestP0.y = select( m_bestfit, p0.y, bestP0.y );
            bestP0.z = select( m_bestfit, p0.z, bestP0.z );
            vh0      = select( m_bestfit, h0, vh0 );
          } // n_0_candidates
          bestP2.x = select( m_found_h0, p2.x, bestP2.x );
          bestP2.y = select( m_found_h0, p2.y, bestP2.y );
          bestP2.z = select( m_found_h0, p2.z, bestP2.z );
          vh2      = select( m_found_h0, cid2, vh2 );
        } // h2

        auto bestH0 = loop_mask && ( bestFit < m_max_scatter_seeding.value() );

        // Remove used hits
        auto m_remove_p1 = loop_mask && !bestH0;
        auto newP1       = P1->compress_back<SIMDWrapper::InstructionSet::Best>( m_remove_p1 );
        newP1.field<HitsTag::pos>().set( p1 );
        newP1.field<HitsTag::phi>().set( phi1 );
        newP1.field<HitsTag::idx>().set( vh1 );
        if ( none( bestH0 ) ) continue;

        removeHits( P0, bestH0, vh0 );
        removeHits( P2, bestH0, vh2 );

        auto const track = tracks->compress_back( bestH0 );
        track.field<LightTrackTag::sumScatter>().set( bestFit );
        track.field<LightTrackTag::p0>().set( bestP0 );
        track.field<LightTrackTag::p1>().set( p1 );
        track.field<LightTrackTag::p2>().set( bestP2 );
        track.field<LightTrackTag::nHits>().set( 3 );
        track.field<LightTrackTag::skipped>().set( 0 );
        track.field<LightTrackTag::hit>( 0 ).set( vh0 );
        track.field<LightTrackTag::hit>( 1 ).set( vh1 );
        track.field<LightTrackTag::hit>( 2 ).set( vh2 );
      } // h1
    }

    template <size_t N>
    [[gnu::always_inline]] inline void TrackForwarding( const LightTracksSoA* tracks, std::array<HitsPlane*, N> P2s,
                                                        LightTracksSoA* tracks_forwarded,
                                                        LightTracksSoA* tracks_finalized ) const {
      using simd = SIMDWrapper::best::types;
      using I    = simd::int_v;
      using F    = simd::float_v;
      using sF   = SIMDWrapper::scalar::float_v;

      tracks_forwarded->resize( 0 );

      for ( const auto& track : tracks->simd() ) {
        const auto mask = track.loop_mask();

        Vec3<F> p0 = track.get<LightTrackTag::p1>().vec3();
        Vec3<F> p1 = track.get<LightTrackTag::p2>().vec3();

        F       bestFit = m_max_scatter_forwarding.value();
        Vec3<F> bestP2  = p1;
        I       bestH2  = 0;

        const F tx = ( p1.x - p0.x ) / ( p1.z - p0.z );
        const F ty = ( p1.y - p0.y ) / ( p1.z - p0.z );

        for ( HitsPlane* P2 : P2s ) {
          for ( const auto hitP2 : P2->scalar() ) {
            const Vec3<sF> sp2 = hitP2.get<HitsTag::pos>().vec3();
            const Vec3<F>  p2  = Vec3<F>( sp2.x, sp2.y, sp2.z );

            const F dz = p2.z - p0.z;
            const F x  = tx * dz + p0.x;
            const F y  = ty * dz + p0.y;

            const F dx      = x - p2.x;
            const F dy      = y - p2.y;
            const F scatter = dx * dx + dy * dy;

            const auto m_bestfit = mask && ( scatter < bestFit );

            if ( none( m_bestfit ) ) continue;

            bestFit  = select( m_bestfit, scatter, bestFit );
            bestP2.x = select( m_bestfit, p2.x, bestP2.x );
            bestP2.y = select( m_bestfit, p2.y, bestP2.y );
            bestP2.z = select( m_bestfit, p2.z, bestP2.z );
            bestH2   = select( m_bestfit, I( hitP2.get<HitsTag::idx>().cast() ), bestH2 );
          }
        }

        auto bestH0 = mask && ( bestFit < m_max_scatter_forwarding.value() );

        // Finish loading tracks
        const Vec3<F> tp0         = track.get<LightTrackTag::p0>().vec3();
        I             n_hits      = track.get<LightTrackTag::nHits>();
        I             skipped     = track.get<LightTrackTag::skipped>();
        F             sum_scatter = track.get<LightTrackTag::sumScatter>();

        // increment or reset to 0
        skipped = select( bestH0, 0, skipped + 1 );

        // increment only if we found a hit
        n_hits      = select( bestH0, n_hits + 1, n_hits );
        sum_scatter = select( bestH0, sum_scatter + bestFit, sum_scatter );

        p1.x = select( bestH0, p1.x, p0.x );
        p1.y = select( bestH0, p1.y, p0.y );
        p1.z = select( bestH0, p1.z, p0.z );

        const auto max_skips = m_skip_forward + 1;
        auto       m_forward = mask && ( skipped < I( max_skips ) );

        // Remove used hits
        for ( auto P2 : P2s ) removeHits( P2, bestH0, bestH2 );

        // Forward tracks
        {
          auto track_fwd = tracks_forwarded->compress_back<SIMDWrapper::InstructionSet::Best>( m_forward );

          track_fwd.template field<LightTrackTag::sumScatter>().set( sum_scatter );
          track_fwd.template field<LightTrackTag::p0>().set( tp0 );
          track_fwd.template field<LightTrackTag::p1>().set( p1 );
          track_fwd.template field<LightTrackTag::p2>().set( bestP2 );
          track_fwd.template field<LightTrackTag::nHits>().set( n_hits );
          track_fwd.template field<LightTrackTag::skipped>().set( skipped );

          int max_n_hits = n_hits.hmax( m_forward );
          for ( int j = 0; j < max_n_hits; j++ ) {
            auto push_hit = bestH0 && ( n_hits == I( j + 1 ) );
            I    id       = select( push_hit, bestH2, track.get<LightTrackTag::hit>( j ) );
            track_fwd.template field<LightTrackTag::hit>( j ).set( id );
          }
        }

        // Finalize track
        auto m_final = mask && !m_forward && ( n_hits > 3 || sum_scatter < m_max_scatter_3hits.value() );

        if ( none( m_final ) ) continue; // Nothing to finalize

        {
          auto track_final = tracks_finalized->compress_back<SIMDWrapper::InstructionSet::Best>( m_final );

          track_final.template field<LightTrackTag::sumScatter>().set( sum_scatter );
          track_final.template field<LightTrackTag::p0>().set( tp0 );
          track_final.template field<LightTrackTag::p1>().set( p1 );
          track_final.template field<LightTrackTag::p2>().set( bestP2 );
          track_final.template field<LightTrackTag::nHits>().set( n_hits );
          track_final.template field<LightTrackTag::skipped>().set( skipped );

          int max_n_hits = n_hits.hmax( m_final );
          for ( int j = 0; j < max_n_hits; j++ ) {
            I id = track.get<LightTrackTag::hit>( j );
            track_final.template field<LightTrackTag::hit>( j ).set( id );
          }
        }
      }
    }

    [[gnu::always_inline]] inline void copy_remaining( const LightTracksSoA* tracks_candidates,
                                                       LightTracksSoA*       tracks ) const {
      using simd = SIMDWrapper::best::types;
      using I    = simd::int_v;
      using F    = simd::float_v;
      for ( const auto& trackcand : tracks_candidates->simd() ) {
        auto loop_mask = trackcand.loop_mask();

        I n_hits      = trackcand.get<LightTrackTag::nHits>();
        F sum_scatter = trackcand.get<LightTrackTag::sumScatter>();

        auto m_final = loop_mask && ( n_hits > 3 || sum_scatter < m_max_scatter_3hits.value() );

        auto track = tracks->compress_back( m_final );

        track.template field<LightTrackTag::sumScatter>().set( sum_scatter );
        track.template field<LightTrackTag::p0>().set( trackcand.get<LightTrackTag::p0>().vec3() );
        track.template field<LightTrackTag::p1>().set( trackcand.get<LightTrackTag::p1>().vec3() );
        track.template field<LightTrackTag::p2>().set( trackcand.get<LightTrackTag::p2>().vec3() );
        track.template field<LightTrackTag::nHits>().set( n_hits );
        track.template field<LightTrackTag::skipped>().set( trackcand.get<LightTrackTag::skipped>() );

        int max_n_hits = n_hits.hmax( m_final );
        for ( int j = 0; j < max_n_hits; j++ ) {
          const auto id = trackcand.get<LightTrackTag::hit>( j );
          track.template field<LightTrackTag::hit>( j ).set( id );
        }
      }
    }

    //=============================================================================
    // Main execution
    //=============================================================================
    std::tuple<VP::Hits, Velo::Tracks, Velo::Tracks>
    operator()( const EventContext& evtCtx, const RawBank::View& tBanks, const DeVP& devp ) const override {
      using Tracks = Velo::Tracks;
      using Hits   = VP::Hits;

      std::tuple<Hits, Tracks, Tracks> result{
          Hits( Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx ) ),
          Tracks( true, Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx ) ),
          Tracks( Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx ) )};
      auto& [hits, tracksBackward, tracksForward] = result;
      hits.reserve( 10000 );

      if ( tBanks.empty() ) return result;

      const unsigned int version           = tBanks[0]->version();
      bool               supported_version = check_version<VP_banktype>( version );

      if ( !supported_version ) {
        warning() << "Unsupported raw bank version (" << version << ")" << endmsg;
        return result;
      }

      RawEvent rawEvent_sorted;
      rawEvent_sorted.reserve( VPInfos::NSensors );
      sortBanks<VP_banktype>( tBanks, rawEvent_sorted );

      const LHCb::RawBank* VPRawBanks[VPInfos::NSensors] = {};

      // Copy rawbanks pointers to protect against unordered data
      for ( auto iterBank = rawEvent_sorted.banks( VP_banktype ).begin();
            iterBank != rawEvent_sorted.banks( VP_banktype ).end(); iterBank++ ) {
        const uint32_t sensor = ( *iterBank )->sourceID();
        VPRawBanks[sensor]    = *iterBank;
      }

      if constexpr ( searchMode == SearchMode::Full ) {
        // Swap sensor ids so they are grouped by z
        for ( uint32_t s = 0; s < VPInfos::NSensors; s += VPInfos::NSensorsPerModule ) {
          std::swap( VPRawBanks[s + 1], VPRawBanks[s + 3] );
        }
      }

      // Tracking buffers
      constexpr auto num_P = SearchModeParam<searchMode>::plane_window;
      auto           raw_P =
          LHCb::make_object_array<HitsPlane, num_P>( Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx ) );
      std::array<HitsPlane*, num_P> P{nullptr};
      for ( std::size_t i = 0; i < num_P; i++ ) {
#if defined( __GNUC__ ) && __GNUC__ < 13
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Warray-bounds"
#endif
        raw_P[i].reserve( reserved_hits_inPlane );
        P[i] = &raw_P[i];
#if defined( __GNUC__ ) && __GNUC__ < 13
#  pragma GCC diagnostic pop
#endif
      }

      LightTracksSoA t_candidates{Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx )};
      LightTracksSoA t_forwarded{Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx )};
      LightTracksSoA tracks{Zipping::generateZipIdentifier(), LHCb::getMemResource( evtCtx )};
      t_candidates.reserve( reserved_tracks_size );
      t_forwarded.reserve( reserved_tracks_size );
      tracks.reserve( reserved_tracks_size );
      LightTracksSoA *tracks_candidates = &t_candidates, *tracks_forwarded = &t_forwarded;

      // Do tracking backward
      if constexpr ( searchMode == SearchMode::Fast ) {
        auto raw_P2 = LHCb::make_object_array<HitsPlane, num_P>( Zipping::generateZipIdentifier(),
                                                                 LHCb::getMemResource( evtCtx ) );
        std::array<HitsPlane*, num_P> P2{nullptr};
        for ( std::size_t i = 0; i < num_P; i++ ) {
#if defined( __GNUC__ ) && __GNUC__ < 13
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Warray-bounds"
#endif
          raw_P2[i].reserve( reserved_hits_inPlane );
          P2[i] = &raw_P2[i];
#if defined( __GNUC__ ) && __GNUC__ < 13
#  pragma GCC diagnostic pop
#endif
        }

        const int i0 = VPInfos::NPlanes - 1, i1 = VPInfos::NPlanes - 2, i2 = VPInfos::NPlanes - 3;

        // Side 1
        int sensor0 = i0 * VPInfos::NSensorsPerPlane;
        int sensor1 = i1 * VPInfos::NSensorsPerPlane;
        int sensor2 = i2 * VPInfos::NSensorsPerPlane;
        getHits<VP_banktype>( VPRawBanks, sensor0, sensor0 + VPInfos::NSensorsPerModule, devp, *P[0], hits,
                              m_sensorMasks );
        getHits<VP_banktype>( VPRawBanks, sensor1, sensor1 + VPInfos::NSensorsPerModule, devp, *P[1], hits,
                              m_sensorMasks );
        getHits<VP_banktype>( VPRawBanks, sensor2, sensor2 + VPInfos::NSensorsPerModule, devp, *P[2], hits,
                              m_sensorMasks );
        TrackSeeding( P[0], P[1], P[2], tracks_candidates );

        // Side 2
        sensor0 += VPInfos::NSensorsPerModule;
        sensor1 += VPInfos::NSensorsPerModule;
        sensor2 += VPInfos::NSensorsPerModule;

        getHits<VP_banktype>( VPRawBanks, sensor0, sensor0 + VPInfos::NSensorsPerModule, devp, *P2[0], hits,
                              m_sensorMasks );
        getHits<VP_banktype>( VPRawBanks, sensor1, sensor1 + VPInfos::NSensorsPerModule, devp, *P2[1], hits,
                              m_sensorMasks );
        getHits<VP_banktype>( VPRawBanks, sensor2, sensor2 + VPInfos::NSensorsPerModule, devp, *P2[2], hits,
                              m_sensorMasks );
        TrackSeeding( P2[0], P2[1], P2[2], tracks_candidates );

        std::rotate( P.begin(), P.begin() + 1, P.end() );
        std::rotate( P2.begin(), P2.begin() + 1, P2.end() );

        for ( int p = VPInfos::NPlanes - 4; p >= 0; p-- ) {
          const int sensor = p * VPInfos::NSensorsPerPlane;
          getHits<VP_banktype>( VPRawBanks, sensor, sensor + VPInfos::NSensorsPerModule, devp, *P[2], hits,
                                m_sensorMasks );
          getHits<VP_banktype>( VPRawBanks, sensor + VPInfos::NSensorsPerModule,
                                sensor + 2 * VPInfos::NSensorsPerModule, devp, *P2[2], hits, m_sensorMasks );

          TrackForwarding( tracks_candidates, std::array<HitsPlane*, 2>( {P[2], P2[2]} ), tracks_forwarded, &tracks );

          std::swap( tracks_candidates, tracks_forwarded );

          TrackSeeding( P[0], P[1], P[2], tracks_candidates );
          TrackSeeding( P2[0], P2[1], P2[2], tracks_candidates );

          std::rotate( P.begin(), P.begin() + 1, P.end() );
          std::rotate( P2.begin(), P2.begin() + 1, P2.end() );
        }

        copy_remaining( tracks_candidates, &tracks );

      } else if constexpr ( searchMode == SearchMode::Full ) {

        // Prologue
        constexpr auto planes = SearchModeParam<searchMode>::planes;
        constexpr auto stride = SearchModeParam<searchMode>::planes_stride;

        int sensor0 = ( planes - 1 ) * stride;
        int sensor1 = ( planes - 2 ) * stride;
        int sensor2 = ( planes - 3 ) * stride;
        getHits<VP_banktype>( VPRawBanks, sensor0, sensor0 + stride, devp, *P[0], hits, m_sensorMasks );
        getHits<VP_banktype>( VPRawBanks, sensor1, sensor1 + stride, devp, *P[1], hits, m_sensorMasks );
        getHits<VP_banktype>( VPRawBanks, sensor2, sensor2 + stride, devp, *P[2], hits, m_sensorMasks );

        TrackSeeding( P[0], P[1], P[2], tracks_candidates );
        int last_sensor = 3;

        for ( int p = planes - 4; p >= 0; p-- ) {
          const int sensor = p * stride;
          getHits<VP_banktype>( VPRawBanks, sensor, sensor + stride, devp, *P[last_sensor], hits, m_sensorMasks );
          TrackForwarding( tracks_candidates, std::array<HitsPlane*, 1>( {P[last_sensor]} ), tracks_forwarded,
                           &tracks );
          std::swap( tracks_candidates, tracks_forwarded );

          for ( int j = 1; j < last_sensor; j++ ) {
            for ( int k = 0; k < j; k++ ) { TrackSeeding( P[k], P[j], P[last_sensor], tracks_candidates ); }
          }

          if ( last_sensor < num_P - 1 )
            last_sensor++;
          else
            std::rotate( P.begin(), P.begin() + 1, P.end() );
        }

        copy_remaining( tracks_candidates, &tracks );
      } else {

        const int i0 = VPInfos::NPlanes - 1, i1 = VPInfos::NPlanes - 2, i2 = VPInfos::NPlanes - 3;

        const int sensor0 = i0 * VPInfos::NSensorsPerPlane;
        const int sensor1 = i1 * VPInfos::NSensorsPerPlane;
        const int sensor2 = i2 * VPInfos::NSensorsPerPlane;

        getHits<VP_banktype>( VPRawBanks, sensor0, sensor0 + VPInfos::NSensorsPerPlane, devp, *P[0], hits,
                              m_sensorMasks );
        getHits<VP_banktype>( VPRawBanks, sensor1, sensor1 + VPInfos::NSensorsPerPlane, devp, *P[1], hits,
                              m_sensorMasks );
        getHits<VP_banktype>( VPRawBanks, sensor2, sensor2 + VPInfos::NSensorsPerPlane, devp, *P[2], hits,
                              m_sensorMasks );
        TrackSeeding( P[0], P[1], P[2], tracks_candidates );
        std::rotate( P.begin(), P.begin() + 1, P.end() );

        for ( int p = VPInfos::NPlanes - 4; p >= 0; p-- ) {
          const int sensor = p * VPInfos::NSensorsPerPlane;
          getHits<VP_banktype>( VPRawBanks, sensor, sensor + VPInfos::NSensorsPerPlane, devp, *P[2], hits,
                                m_sensorMasks );
          TrackForwarding( tracks_candidates, std::array<HitsPlane*, 1>( {P[2]} ), tracks_forwarded, &tracks );
          std::swap( tracks_candidates, tracks_forwarded );
          TrackSeeding( P[0], P[1], P[2], tracks_candidates );
          std::rotate( P.begin(), P.begin() + 1, P.end() );
        }

        copy_remaining( tracks_candidates, &tracks );
      } // backward tracking

      // Special recovery step for hits in same module but different sensors:
      if constexpr ( searchMode == SearchMode::Full ) {
        // Find start-end of each module
        std::vector<int> module_start;
        int              prev_module = -1;
        for ( const auto hit : hits.scalar() ) {
          auto id  = Detector::VPChannelID( hit.get<VP::VPHitsTag::ChannelId>().cast() );
          int  mod = id.module();
          if ( prev_module != mod ) {
            module_start.push_back( hit.offset() );
            prev_module = mod;
          }
        }
        module_start.push_back( hits.size() );

        std::vector<int> newhits;
        newhits.reserve( VPInfos::NModules );

        // for every track check extra hits in same modules
        for ( auto track : tracks.scalar() ) {
          auto p0 = track.get<LightTrackTag::p0>().vec3();
          auto p2 = track.get<LightTrackTag::p2>().vec3();
          auto d  = p0 - p2;

          auto tx = d.x / d.z;
          auto ty = d.y / d.z;
          auto x0 = p2.x - p2.z * tx;
          auto y0 = p2.y - p2.z * ty;

          auto nHits = track.get<LightTrackTag::nHits>().cast();

          newhits.clear();
          for ( int i = 0; i < nHits; i++ ) {
            auto hit_idx = track.get<LightTrackTag::hit>( i ).cast();
            newhits.push_back( hit_idx );

            auto mod_start = *( std::upper_bound( module_start.begin(), module_start.end(), hit_idx ) - 1 );
            auto mod_end   = *std::upper_bound( module_start.begin(), module_start.end(), hit_idx );

            auto sensor =
                Detector::VPChannelID( hits.scalar()[hit_idx].template get<VP::VPHitsTag::ChannelId>().cast() )
                    .sensor();
            for ( int hit_idx2 = mod_start; hit_idx2 < mod_end; hit_idx2++ ) {
              // check if sensor !=
              auto sensor2 =
                  Detector::VPChannelID( hits.scalar()[hit_idx2].template get<VP::VPHitsTag::ChannelId>().cast() )
                      .sensor();
              if ( sensor == sensor2 ) continue;

              // check if point on track
              auto p        = hits.scalar()[hit_idx2].template get<VP::VPHitsTag::pos>().vec3();
              auto dx       = x0 + tx * p.z - p.x;
              auto dy       = y0 + ty * p.z - p.y;
              auto distance = dx * dx + dy * dy;
              if ( distance.cast() > m_max_scatter_forwarding.value() ) continue;

              newhits.push_back( hit_idx2 );
              break;
            }
          }

          if ( (int)newhits.size() == nHits ) continue;

          std::sort( newhits.begin(), newhits.end() );
          for ( int i = 0; i < (int)newhits.size(); i++ ) { track.field<LightTrackTag::hit>( i ).set( newhits[i] ); }
          track.field<LightTrackTag::nHits>().set( newhits.size() );
        }
      }

      // Make LHCb tracks
      tracksBackward.reserve( reserved_tracks_size ); // reserve capacity at a reasonable large size
      tracksForward.reserve( reserved_tracks_size );
      using simd = SIMDWrapper::best::types;
      using I    = simd::int_v;
      using F    = simd::float_v;
      for ( const auto track : tracks.simd() ) {
        auto loop_mask = track.loop_mask();

        // Simple fit
        auto p1 = track.get<LightTrackTag::p0>().vec3();
        auto p2 = track.get<LightTrackTag::p2>().vec3();
        auto d  = p1 - p2;

        auto tx = d.x / d.z;
        auto ty = d.y / d.z;
        auto x0 = p2.x - p2.z * tx;
        auto y0 = p2.y - p2.z * ty;

        F z_beam = p1.z;
        F denom  = tx * tx + ty * ty;
        z_beam   = select( denom < 0.001f * 0.001f, z_beam, -( x0 * tx + y0 * ty ) / denom );

        const auto hitproxy = hits.simd();
        auto       n_hits   = track.get<LightTrackTag::nHits>();

        // Store backward tracks
        auto        backwards = ( z_beam > p2.z ) && loop_mask;
        auto const& bwd       = tracksBackward.compress_back( backwards );
        auto        bwd_hits  = bwd.template field<TracksTag::Hits>();
        bwd_hits.resize( n_hits );
        for ( int h = 0; h < n_hits.hmax( backwards ); h++ ) {
          auto hit_index = track.get<LightTrackTag::hit>( h );
          auto hit       = hitproxy.gather( hit_index, backwards && h < n_hits );
          auto id        = LHCb::Event::lhcbid_v<simd>(
              LHCbID::make( LHCbID::channelIDtype::VP, hit.template get<VP::VPHitsTag::ChannelId>() ) );
          bwd_hits[h].template field<TracksTag::Index>().set( hit_index );
          bwd_hits[h].template field<TracksTag::LHCbID>().set( id );
        }
        bwd.template field<TracksTag::States>( 0 ).set( 1.f, 1.f, 1.f, tx, ty );

        // Store forward tracks
        auto        forwards = ( !backwards ) && loop_mask;
        auto const& fwd      = tracksForward.compress_back( forwards );
        auto        fwd_hits = fwd.template field<TracksTag::Hits>();
        fwd_hits.resize( n_hits );
        for ( int h = 0; h < n_hits.hmax( forwards ); h++ ) {
          auto hit_index = track.get<LightTrackTag::hit>( h );
          auto hit       = hitproxy.gather( hit_index, forwards && h < n_hits );
          auto id        = LHCb::Event::lhcbid_v<simd>(
              LHCbID::make( LHCbID::channelIDtype::VP, hit.template get<VP::VPHitsTag::ChannelId>() ) );
          fwd_hits[h].template field<TracksTag::Index>().set( hit_index );
          fwd_hits[h].template field<TracksTag::LHCbID>().set( id );
        }
        fwd.template field<TracksTag::States>( 1 ).set( x0 + StateParameters::ZEndVelo * tx,
                                                        y0 + StateParameters::ZEndVelo * ty, StateParameters::ZEndVelo,
                                                        tx, ty );
      } // loop all tracks

      // Fit forwards
      for ( auto const& track : tracksForward.simd() ) {
        auto    loop_mask     = track.loop_mask();
        I       nhits         = track.nHits();
        Vec3<F> dir           = track.StateDir( 1 );
        auto    vp_index      = track.vp_indices();
        auto    closestToBeam = fitBackward<F, I>( loop_mask, nhits, hits, dir, vp_index );
        closestToBeam.transportTo( closestToBeam.zBeam() );
        track.StatePosDir( 0 ).setPosition( closestToBeam.pos() );
        track.StatePosDir( 0 ).setDirection( closestToBeam.dir() );
        track.setStateCovXY( 0, closestToBeam.covX(), closestToBeam.covY() );
        track.setStateCovXY( 1, Vec3<F>( 100.f, 0.f, 1.f ), Vec3<F>( 100.f, 0.f, 1.f ) ); //  avoid NAN in chechers
      }
      // Fit backwards
      for ( auto const& track : tracksBackward.simd() ) {
        auto    loop_mask     = track.loop_mask();
        I       nhits         = track.nHits();
        Vec3<F> dir           = track.StateDir( 0 );
        auto    vp_index      = track.vp_indices();
        auto    closestToBeam = fitForward<F, I>( loop_mask, nhits, hits, dir, vp_index );
        closestToBeam.transportTo( closestToBeam.zBeam() );
        track.StatePosDir( 0 ).setPosition( closestToBeam.pos() );
        track.StatePosDir( 0 ).setDirection( closestToBeam.dir() );
        track.setStateCovXY( 0, closestToBeam.covX(), closestToBeam.covY() );
        track.setStateCovXY( 1, Vec3<F>( 100.f, 0.f, 1.f ), Vec3<F>( 100.f, 0.f, 1.f ) ); //  avoid NAN in chechers
      }
      m_nbClustersCounter += hits.size();
      m_nbTracksCounter += tracks.size();
      return result;
    }

  private:
    mutable Gaudi::Accumulators::SummingCounter<>        m_nbClustersCounter{this, "Nb of Produced Clusters"};
    mutable Gaudi::Accumulators::SummingCounter<>        m_nbTracksCounter{this, "Nb of Produced Tracks"};
    Gaudi::Property<std::array<bool, VPInfos::NSensors>> m_sensorMasks{this, "SensorMasks", {}};
    Gaudi::Property<float>                               m_max_scatter_seeding{this, "MaxScatterSeeding", 0.1f};
    Gaudi::Property<float>                               m_max_scatter_forwarding{this, "MaxScatterForwarding", 0.1f};
    Gaudi::Property<float>                               m_max_scatter_3hits{this, "MaxScatter3hits", 0.02f};
    Gaudi::Property<int> m_skip_forward{this, "SkipForward", SearchModeParam<searchMode>::skip_forward_default};
  };

  using VeloClusterTrackingSIMD             = ClusterTrackingSIMD<RawBank::VP, SearchMode::Default>;
  using VeloClusterTrackingSIMDFaster       = ClusterTrackingSIMD<RawBank::VP, SearchMode::Fast>;
  using VeloClusterTrackingSIMDFull         = ClusterTrackingSIMD<RawBank::VP, SearchMode::Full>;
  using VeloRetinaClusterTrackingSIMD       = ClusterTrackingSIMD<RawBank::VPRetinaCluster, SearchMode::Default>;
  using VeloRetinaClusterTrackingSIMDFaster = ClusterTrackingSIMD<RawBank::VPRetinaCluster, SearchMode::Fast>;
  using VeloRetinaClusterTrackingSIMDFull   = ClusterTrackingSIMD<RawBank::VPRetinaCluster, SearchMode::Full>;
  DECLARE_COMPONENT_WITH_ID( VeloClusterTrackingSIMD, "VeloClusterTrackingSIMD" )
  DECLARE_COMPONENT_WITH_ID( VeloClusterTrackingSIMDFaster, "VeloClusterTrackingSIMDFaster" )
  DECLARE_COMPONENT_WITH_ID( VeloClusterTrackingSIMDFull, "VeloClusterTrackingSIMDFull" )
  DECLARE_COMPONENT_WITH_ID( VeloRetinaClusterTrackingSIMD, "VeloRetinaClusterTrackingSIMD" )
  DECLARE_COMPONENT_WITH_ID( VeloRetinaClusterTrackingSIMDFaster, "VeloRetinaClusterTrackingSIMDFaster" )
  DECLARE_COMPONENT_WITH_ID( VeloRetinaClusterTrackingSIMDFull, "VeloRetinaClusterTrackingSIMDFull" )
} // namespace LHCb::Pr::Velo
