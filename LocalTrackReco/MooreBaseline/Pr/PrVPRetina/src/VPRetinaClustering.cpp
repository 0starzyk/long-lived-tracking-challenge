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
// LHCb
#include "DetDesc/IConditionDerivationMgr.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"

// Rec
#include "PrKernel/PrPixelUtils.h"

#include "GaudiKernel/Transform3DTypes.h"

// Local
#include "VPDet/VPDetPaths.h"
#include "VPRetinaClustering.h"
#include <iomanip>

DECLARE_COMPONENT( VPRetinaClustering )

namespace {
  using namespace PrPixel;
}

namespace {
  struct SPCache {
    std::array<float, 4> fxy;
    unsigned char        pattern;
    unsigned char        nx1;
    unsigned char        nx2;
    unsigned char        ny1;
    unsigned char        ny2;
  };
  //=========================================================================
  // Cache Super Pixel cluster patterns.
  //=========================================================================
  auto create_SPPatterns() {
    std::array<SPCache, 256> SPCaches;
    // create a cache for all super pixel cluster patterns.
    // this is an unoptimized 8-way flood fill on the 8 pixels
    // in the super pixel.
    // no point in optimizing as this is called once in
    // initialize() and only takes about 20 us.

    // define deltas to 8-connectivity neighbours
    const int dx[] = {-1, 0, 1, -1, 0, 1, -1, 1};
    const int dy[] = {-1, -1, -1, 1, 1, 1, 0, 0};

    // clustering buffer for isolated superpixels.
    unsigned char sp_buffer[8];

    // SP index buffer and its size for single SP clustering
    unsigned char sp_idx[8];
    unsigned char sp_idx_size = 0;

    // stack and stack pointer for single SP clustering
    unsigned char sp_stack[8];
    unsigned char sp_stack_ptr = 0;

    // loop over all possible SP patterns
    for ( unsigned int sp = 0; sp < 256; ++sp ) {
      sp_idx_size = 0;
      for ( unsigned int shift = 0; shift < 8; ++shift ) {
        const unsigned char p = sp & ( 1 << shift );
        sp_buffer[shift]      = p;
        if ( p ) { sp_idx[sp_idx_size++] = shift; }
      }

      // loop over pixels in this SP and use them as
      // cluster seeds.
      // note that there are at most two clusters
      // in a single super pixel!
      unsigned char clu_idx = 0;
      for ( unsigned int ip = 0; ip < sp_idx_size; ++ip ) {
        unsigned char idx = sp_idx[ip];

        if ( 0 == sp_buffer[idx] ) { continue; } // pixel is used

        sp_stack_ptr             = 0;
        sp_stack[sp_stack_ptr++] = idx;
        sp_buffer[idx]           = 0;
        unsigned char x          = 0;
        unsigned char y          = 0;
        unsigned char n          = 0;

        while ( sp_stack_ptr ) {
          idx                     = sp_stack[--sp_stack_ptr];
          const unsigned char row = idx % 4;
          const unsigned char col = idx / 4;
          x += col;
          y += row;
          ++n;

          for ( unsigned int ni = 0; ni < 8; ++ni ) {
            const char ncol = col + dx[ni];
            if ( ncol < 0 || ncol > 1 ) continue;
            const char nrow = row + dy[ni];
            if ( nrow < 0 || nrow > 3 ) continue;
            const unsigned char nidx = ( ncol << 2 ) | nrow;
            if ( 0 == sp_buffer[nidx] ) continue;
            sp_stack[sp_stack_ptr++] = nidx;
            sp_buffer[nidx]          = 0;
          }
        }

        const uint32_t cx = x / n;
        const uint32_t cy = y / n;
        const float    fx = x / static_cast<float>( n ) - cx;
        const float    fy = y / static_cast<float>( n ) - cy;

        // store the centroid pixel
        SPCaches[sp].pattern |= ( ( cx << 2 ) | cy ) << 4 * clu_idx;

        // set the two cluster flag if this is the second cluster
        SPCaches[sp].pattern |= clu_idx << 3;

        // set the pixel fractions
        SPCaches[sp].fxy[2 * clu_idx]     = fx;
        SPCaches[sp].fxy[2 * clu_idx + 1] = fy;

        // increment cluster count. note that this can only become 0 or 1!
        ++clu_idx;
      }
    }
    return SPCaches;
  }
  // SP pattern buffers for clustering, cached once.
  // There are 256 patterns and there can be at most two
  // distinct clusters in an SP.
  static const std::array<SPCache, 256> s_SPCaches = create_SPPatterns();

  //=============================================================================
  // Format a Cluster word from pixel weighed sum (shift)
  //=============================================================================
  uint32_t FormatCluster( uint32_t col, uint32_t row, uint32_t shift_values ) {
    return ( col << 20 | row << 12 | shift_values );
  }

  uint32_t FormatCluster( uint32_t col, uint32_t row, uint32_t shift_col, uint32_t shift_row, uint32_t n ) {
    uint32_t shift_values = shift_col << 8 | shift_row << 4 | n;
    return FormatCluster( col, row, shift_values );
  }
} // namespace

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
VPRetinaClustering::VPRetinaClustering( const std::string& name, ISvcLocator* pSvcLocator )
    : MultiTransformer(
          name, pSvcLocator,
          {KeyValue{"RawEventLocation", LHCb::RawEventLocation::Default}, KeyValue{"DEVP", LHCb::Det::VP::det_path}},
          {KeyValue{"ClusterLocation", LHCb::VPClusterLocation::Light},
           KeyValue{"ClusterOffsets", LHCb::VPClusterLocation::Offsets}} ) {}

//=============================================================================
// Main execution
//=============================================================================
std::tuple<std::vector<LHCb::VPLightCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>> VPRetinaClustering::
                                                                                                 operator()( const LHCb::RawEvent& rawEvent, const DeVP& devp ) const {
  auto result = std::tuple<std::vector<LHCb::VPLightCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>>{};

  const auto& tBanks = rawEvent.banks( LHCb::RawBank::VP );
  if ( tBanks.empty() ) return result;

  const unsigned int version = tBanks[0]->version();
  if ( version != 2 ) {
    warning() << "Unsupported raw bank version (" << version << ")" << endmsg;
    return result;
  }

  // Since the pool is local, to first preallocate the pool, then count hits per module,
  // and then preallocate per module and move hits might not be faster than adding
  // directly to the PrPixelModuleHits (which would require more allocations, but
  // not many if we start with a sensible default)
  auto& [pool, offsets] = result;
  const unsigned int startSize = 10000U;
  pool.reserve( startSize );

  // Loop over VP RawBanks
  for ( auto iterBank : tBanks ) {

    const auto         sensor = LHCb::Detector::VPChannelID::SensorID( iterBank->sourceID() );
    const unsigned int module = 1 + ( to_unsigned( sensor ) / VP::NSensorsPerModule );
    if ( m_modulesToSkipMask[module - 1] ) { continue; }

    const uint32_t* bank = iterBank->data();

    auto retinaClusters = makeRetinaClusters( bank );

    std::transform( begin( retinaClusters ), end( retinaClusters ), std::back_inserter( pool ),
                    [&]( uint32_t iter ) -> LHCb::VPLightCluster {
                      const uint32_t col = iter >> 20;
                      const auto     row = LHCb::Detector::VPChannelID::RowID{( iter >> 12 ) & 0xFF};
                      const uint32_t shift_col = ( iter >> 8 ) & 0xF;
                      const uint32_t shift_row = ( iter >> 4 ) & 0xF;
                      const uint32_t n = iter & 0xF;

                      const uint32_t cx = col + shift_col / n;
                      const auto     cy = LHCb::Detector::VPChannelID::RowID{to_unsigned( row ) + shift_row / n};
                      const float    fx = shift_col / static_cast<float>( n ) - shift_col / n;
                      const float    fy = shift_row / static_cast<float>( n ) - shift_row / n;

                      const auto chip = LHCb::Detector::VPChannelID::ChipID{cx / CHIP_COLUMNS};
                      const auto ccol = LHCb::Detector::VPChannelID::ColumnID{cx % CHIP_COLUMNS};

                      LHCb::Detector::VPChannelID cid( sensor, chip, ccol, cy );
                      const float                 local_x = devp.local_x( cx ) + fx * devp.x_pitch( cx );
                      const float                 local_y = ( to_unsigned( cy ) + 0.5 + fy ) * devp.pixel_size();

                      const auto& ltg = devp.ltg( sensor );
                      const float gx = ( ltg[0] * local_x + ltg[1] * local_y + ltg[2] );
                      const float gy = ( ltg[3] * local_x + ltg[4] * local_y + ltg[5] );
                      const float gz = ( ltg[6] * local_x + ltg[7] * local_y + ltg[8] );

                      return {1, 1, gx, gy, gz, cid};
                    } ); // transform over all retinaClusters
    offsets[module] += retinaClusters.size();
  } // loop over all banks

  std::partial_sum( offsets.begin(), offsets.end(), offsets.begin() );

  bool odd = false;
  for ( size_t moduleID = 0; moduleID < VeloInfo::Numbers::NModules; ++moduleID ) {
    // if( msgLevel(MSG::DEBUG)){
    //   debug()<<"Sorting hits in moduleID by phi, usign x,y information "<<moduleID<<endmsg;
    // }
    odd = moduleID % 2 == 1;
    // In even modules you fall in the branching at -180, 180 degrees, you want to do that continuos
    std::stable_sort( pool.begin() + offsets[moduleID], pool.begin() + offsets[moduleID + 1],
                      [&odd]( const LHCb::VPLightCluster& a, const LHCb::VPLightCluster& b ) {
                        // sorting in phi for even modules
                        return (
                            // odd modules, change in y value
                            ( odd && ( a.y() < 0.f && b.y() > 0.f ) ) ||
                            // or even modules, change in y value, but swap logic
                            ( !odd && ( a.y() > 0.f && b.y() < 0.f ) ) ||
                            // same y side even and odd modules, check y1/x1 < y2/x2
                            ( ( a.y() * b.y() ) > 0.f && ( a.y() * b.x() < b.y() * a.x() ) ) );
                      } );
  }
  return result;
}

//=============================================================================
// make RetinaClusters from bank
//=============================================================================
std::vector<uint32_t> VPRetinaClustering::makeRetinaClusters( const uint32_t* bank ) const {
  const uint32_t nsp = *bank++;

  std::vector<VPRetinaMatrix> RetinaMatrixVector;
  RetinaMatrixVector.reserve( 20 );
  std::vector<uint32_t> RetinaCluster;
  RetinaCluster.reserve( 40 );

  // Read super pixel
  for ( unsigned int i = 0; i < nsp; ++i ) {
    const uint32_t sp_word = *bank++;

    uint8_t sp = sp_word & 0xFFU;

    if ( 0 == sp ) continue; // protect against zero super pixels.

    const uint32_t sp_addr          = ( sp_word & 0x007FFF00U ) >> 8;
    const uint32_t sp_row           = sp_addr & 0x3FU;
    const uint32_t sp_col           = ( sp_addr >> 6 );
    const uint32_t no_sp_neighbours = sp_word & 0x80000000U;

    // if a super pixel is isolated the clustering boils
    // down to a simple pattern look up.
    if ( no_sp_neighbours ) {

      // there is always at least one cluster in the super
      // pixel. look up the pattern and add it.

      /* // lines commented because caches needs to be rewrited
      RetinaCluster.push_back(FormatCluster(sp_col*2, sp_row*4, m_SPCaches[sp] & 0x00000FFF));

      // if there is a second cluster for this pattern
      // add it as well.
      if (m_SPCaches[sp] & 0x00FFF000) {
        line commented because caches needs to be rewrited //RetinaCluster.push_back(FormatCluster(sp_col*2, sp_row*4,
      m_SPCaches[sp] & 0x00FFF000));
      } */

      // remove after caches rewrite
      const auto&    spcache = s_SPCaches[sp];
      const uint32_t idx     = spcache.pattern;

      const uint32_t row = idx & 0x03U;
      const uint32_t col = ( idx >> 2 ) & 1;

      const uint32_t n  = 10;
      const uint32_t fx = ( uint32_t )( spcache.fxy[0] * n );
      const uint32_t fy = ( uint32_t )( spcache.fxy[1] * n );

      RetinaCluster.push_back( FormatCluster( sp_col * 2 + col, sp_row * 4 + row, fx, fy, n ) );
      if ( idx & 8 ) {
        const uint32_t row = ( idx >> 4 ) & 3;
        const uint32_t col = ( idx >> 6 ) & 1;

        const uint32_t fx = ( uint32_t )( spcache.fxy[2] * n );
        const uint32_t fy = ( uint32_t )( spcache.fxy[3] * n );

        RetinaCluster.push_back( FormatCluster( sp_col * 2 + col, sp_row * 4 + row, fx, fy, n ) );
      }

      continue;
    }

    // this one is not isolated, we fill a Retina

    // we look for already created Retina
    auto iterRetina = std::find_if( RetinaMatrixVector.begin(), RetinaMatrixVector.end(),
                                    [&]( const VPRetinaMatrix& m ) { return m.IsInRetina( sp_row, sp_col ); } );
    if ( iterRetina != RetinaMatrixVector.end() ) {
      ( *iterRetina ).AddSP( sp_row, sp_col, sp );
    } else {
      RetinaMatrixVector.emplace_back( sp_row, sp_col, sp );
    }
  } // loop over super pixels in raw bank

  // searchRetinaCluster
  for ( auto& m : RetinaMatrixVector ) {
    const auto& clusters = m.SearchCluster();
    RetinaCluster.insert( end( RetinaCluster ), begin( clusters ), end( clusters ) );
  };

  return RetinaCluster;
}
