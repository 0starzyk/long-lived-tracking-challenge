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
#include "Event/RawEvent.h"
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Kernel/VPConstants.h"
#include "PrKernel/PrPixelUtils.h"
#include "VPDAQ/VPRetinaClusterConstants.h"

// Local
#include "VPClus.h"
#include <VPDet/VPDetPaths.h>
#include <algorithm>
#include <iomanip>

DECLARE_COMPONENT_WITH_ID( LHCb::Pr::Velo::VPClus, "VPClus" )

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
  //=============================================================================
  // Sort super pixels by column (major) and row (minor)
  //=============================================================================
  auto SPLowerThan = []( unsigned int lhs, unsigned int rhs ) { return ( lhs & 0x7FFF00 ) < ( rhs & 0x7FFF00 ); };
  auto SPEqual     = []( unsigned int lhs, unsigned int rhs ) { return ( lhs & 0x7FFF00 ) == ( rhs & 0x7FFF00 ); };

  // SP pattern buffers for clustering, cached once.
  // There are 256 patterns and there can be at most two
  // distinct clusters in an SP.
  static const std::array<SPCache, 256> s_SPCaches = create_SPPatterns();

  using namespace PrPixel;
} // namespace

namespace LHCb::Pr::Velo {
  //=============================================================================
  // Standard constructor, initializes variables
  //=============================================================================
  VPClus::VPClus( const std::string& name, ISvcLocator* pSvcLocator )
      : MultiTransformer(
            name, pSvcLocator,
            {KeyValue{"RawEventLocation", RawEventLocation::Default}, KeyValue{"DEVP", LHCb::Det::VP::det_path}},
            {KeyValue{"ClusterLocation", VPClusterLocation::Light},
             KeyValue{"ClusterOffsets", VPClusterLocation::Offsets}} ) {}

  //=============================================================================
  // Main execution
  //=============================================================================
  std::tuple<std::vector<VPLightCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>> VPClus::
                                                                                             operator()( const EventContext&, const RawEvent& rawEvent, const DeVP& devp ) const {
    auto result = std::tuple<std::vector<VPLightCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>>{};

    const auto& tBanks = rawEvent.banks( RawBank::VP );
    if ( tBanks.empty() ) return result;

    const unsigned int version = tBanks[0]->version();
    if ( version != 2 && version != 4 ) {
      warning() << "Unsupported raw bank version (" << version << ")" << endmsg;
      return result;
    }
    // WARNING:
    // This is a rather long function. Please refrain from breaking this
    // up into smaller functions as this will severely impact the
    // timing performance. And yes, this has been measured. Just don't.

    // Clustering buffers
    std::vector<uint8_t>  buffer( VP::NPixelsPerSensor );
    std::vector<uint32_t> pixel_idx;
    std::vector<uint32_t> stack;

    // reserve a minimal stack
    stack.reserve( 64 );

    // Since the pool is local, to first preallocate the pool, then count hits per module,
    // and then preallocate per module and move hits might not be faster than adding
    // directly to the PrPixelModuleHits (which would require more allocations, but
    // not many if we start with a sensible default)
    auto& [pool, offsets] = result;
    pool.reserve( 10000U );

    const LHCb::RawBank* VPRawBanks[VP::NSensors] = {};

    RawEvent rawEvent_sorted;
    rawEvent_sorted.reserve( VP::NSensors );

    if ( version > 3 ) {
      for ( auto iterBank : tBanks ) {
        const uint32_t        sensor0 = ( ( iterBank->sourceID() ) & 0x1FFU ) << 1;
        const uint32_t        sensor1 = sensor0 + 1;
        std::vector<uint32_t> data0;
        data0.reserve( ( iterBank->range<uint32_t>() ).size() );
        std::vector<uint32_t> data1;
        data1.reserve( ( iterBank->range<uint32_t>() ).size() );
        for ( auto word : iterBank->range<uint32_t>() ) {
          if ( ( ( word >> 23 ) & 0x1U ) ) { // check if SP belongs to sensor1
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

        rawEvent_sorted.addBank( sensor0, LHCb::RawBank::VP, VPRetinaCluster::c_SPBankVersion, data0 );
        rawEvent_sorted.addBank( sensor1, LHCb::RawBank::VP, VPRetinaCluster::c_SPBankVersion, data1 );
      }
    }

    if ( version > 3 ) {
      for ( auto iterBank = rawEvent_sorted.banks( RawBank::VP ).begin();
            iterBank != rawEvent_sorted.banks( RawBank::VP ).end(); iterBank++ ) {
        const uint32_t sensor = ( *iterBank )->sourceID();
        VPRawBanks[sensor] = *iterBank;
      }
    } else {
      for ( auto iterBank = tBanks.begin(); iterBank != tBanks.end(); iterBank++ ) {
        const uint32_t sensor = ( *iterBank )->sourceID();
        VPRawBanks[sensor] = *iterBank;
      }
    }

    // Loop over VP RawBanks
    for ( uint32_t s = 0; s < VP::NSensors; s++ ) {

      if ( VPRawBanks[s] == nullptr ) { continue; }

      const auto         sensor = Detector::VPChannelID::SensorID( s );
      const unsigned int module = 1 + ( s / VP::NSensorsPerModule );
      if ( m_modulesToSkipMask[module - 1] ) { continue; }
      // reset and then fill the super pixel buffer for a sensor
      // memset(m_sp_buffer,0,256*256*3*sizeof(unsigned char));
      // the memset is too slow here. the occupancy is so low that
      // resetting a few elements is *much* faster.
      const unsigned int nrc = pixel_idx.size();
      for ( unsigned int irc = 0; irc < nrc; ++irc ) { buffer[pixel_idx[irc]] = false; }
      pixel_idx.clear();

      const auto& ltg = devp.ltg( sensor );

      const uint32_t* bank = ( VPRawBanks[s] )->data();
      uint32_t        nsp;
      if ( version > 3 ) {
        nsp = ( ( VPRawBanks[s] )->range<uint32_t>() ).size();
      } else {
        nsp = *bank++;
      }

      for ( unsigned int i = 0; i < nsp; ++i ) {
        const uint32_t sp_word = *bank++;
        uint8_t        sp = sp_word & 0xFFU;

        if ( 0 == sp ) continue; // protect against zero super pixels.

        const uint32_t sp_addr = ( sp_word & 0x007FFF00U ) >> 8;
        const uint32_t sp_row = sp_addr & 0x3FU;
        const uint32_t sp_col = ( sp_addr >> 6 );
        const uint32_t no_sp_neighbours = sp_word & 0x80000000U;

        if ( sp_row > ( VP::NRows / 4 - 1 ) ) continue; // protect against super pixels outside sensor coordinates.
        if ( sp_col > ( VP::NColumns * VP::NChipsPerSensor / 2 - 1 ) )
          continue; // protect against super pixels outside sensor coordinates.

        // if a super pixel is isolated, the clustering boils
        // down to a simple pattern look up.
        // don't do this if we run in offline mode where we want to record all
        // contributing channels; in that scenario a few more us are negligible
        // compared to the complication of keeping track of all contributing
        // channel IDs.
        if ( no_sp_neighbours ) {
          const auto&    spcache = s_SPCaches[sp];
          const uint32_t idx = spcache.pattern;

          // there is always at least one cluster in the super
          // pixel. look up the pattern and add it.
          const uint32_t row = idx & 0x03U;
          const uint32_t col = ( idx >> 2 ) & 1;
          const uint32_t cx = sp_col * 2 + col;
          const auto     cy = Detector::VPChannelID::RowID{sp_row * 4 + row};
          const auto     chip = Detector::VPChannelID::ChipID{cx / CHIP_COLUMNS};
          const auto     ccol = Detector::VPChannelID::ColumnID{cx % CHIP_COLUMNS};

          Detector::VPChannelID cid( sensor, chip, ccol, cy );

          const float fx = spcache.fxy[0];
          const float fy = spcache.fxy[1];
          const float local_x = devp.local_x( cx ) + fx * devp.x_pitch( cx );
          const float local_y = ( to_unsigned( cy ) + 0.5f + fy ) * devp.pixel_size();

          // gx
          const float gx = ( ltg[0] * local_x + ltg[1] * local_y + ltg[9] );
          // gy
          const float gy = ( ltg[3] * local_x + ltg[4] * local_y + ltg[10] );
          // gz
          const float gz = ( ltg[6] * local_x + ltg[7] * local_y + ltg[11] );

          pool.emplace_back( 1, 1, gx, gy, gz, cid );
          ++offsets[module];

          // if there is a second cluster for this pattern
          // add it as well.
          if ( idx & 8 ) {
            const uint32_t row = ( idx >> 4 ) & 3;
            const uint32_t col2 = ( idx >> 6 ) & 1;
            const auto     cy = Detector::VPChannelID::RowID{sp_row * 4 + row};

            // to get rid of the div use old ccol
            const auto ccol2 = Detector::VPChannelID::ColumnID{to_unsigned( ccol ) + ( col2 - col )};

            Detector::VPChannelID cid( sensor, chip, ccol2, cy );

            const float fx = spcache.fxy[2];
            const float fy = spcache.fxy[3];
            const float local_x = devp.local_x( cx ) + fx * devp.x_pitch( cx );
            const float local_y = ( to_unsigned( cy ) + 0.5f + fy ) * devp.pixel_size();

            // gx
            const float gx = ( ltg[0] * local_x + ltg[1] * local_y + ltg[9] );
            // gy
            const float gy = ( ltg[3] * local_x + ltg[4] * local_y + ltg[10] );
            // gz
            const float gz = ( ltg[6] * local_x + ltg[7] * local_y + ltg[11] );
            pool.emplace_back( 1, 1, gx, gy, gz, cid );
            ++offsets[module];
          }
          continue; // move on to next super pixel
        }

        // this one is not isolated or we are targeting clusters; record all
        // pixels.
        for ( uint32_t shift = 0; shift < 8; ++shift ) {
          const uint8_t pixel = sp & 1;
          if ( pixel ) {
            const uint32_t row = sp_row * 4 + shift % 4;
            const uint32_t col = sp_col * 2 + shift / 4;
            const uint32_t idx = ( col << 8 ) | row;
            buffer[idx] = pixel;
            pixel_idx.push_back( idx );
          }
          sp = sp >> 1;
          if ( 0 == sp ) break;
        }
      } // loop over super pixels in raw bank

      // the sensor buffer is filled, perform the clustering on
      // clusters that span several super pixels.
      const unsigned int nidx = pixel_idx.size();
      for ( unsigned int irc = 0; irc < nidx; ++irc ) {

        const uint32_t idx = pixel_idx[irc];

        if ( !buffer[idx] ) continue; // pixel is used in another cluster

        // 8-way row scan optimized seeded flood fill from here.
        stack.clear();

        // mark seed as used
        buffer[idx] = false;

        // initialize sums
        unsigned int x = 0;
        unsigned int y = 0;
        unsigned int n = 0;

        // push seed on stack
        stack.push_back( idx );

        // as long as the stack is not exhausted:
        // - pop the stack and add popped pixel to cluster
        // - scan the row to left and right, adding set pixels
        //   to the cluster and push set pixels above and below
        //   on the stack (and delete both from the pixel buffer).
        while ( !stack.empty() ) {

          // pop pixel from stack and add it to cluster
          const uint32_t idx = stack.back();
          stack.pop_back();
          const uint32_t row = idx & 0xFFU;
          const uint32_t col = ( idx >> 8 ) & 0x3FFU;
          x += col;
          y += row;
          ++n;

          // check up and down
          uint32_t u_idx = idx + 1;
          if ( row < VP::NRows - 1 && buffer[u_idx] ) {
            buffer[u_idx] = false;
            stack.push_back( u_idx );
          }
          uint32_t d_idx = idx - 1;
          if ( row > 0 && buffer[d_idx] ) {
            buffer[d_idx] = false;
            stack.push_back( d_idx );
          }

          // scan row to the right
          for ( unsigned int c = col + 1; c < VP::NSensorColumns; ++c ) {
            const uint32_t nidx = ( c << 8 ) | row;
            // check up and down
            u_idx = nidx + 1;
            if ( row < VP::NRows - 1 && buffer[u_idx] ) {
              buffer[u_idx] = false;
              stack.push_back( u_idx );
            }
            d_idx = nidx - 1;
            if ( row > 0 && buffer[d_idx] ) {
              buffer[d_idx] = false;
              stack.push_back( d_idx );
            }
            // add set pixel to cluster or stop scanning
            if ( buffer[nidx] ) {
              buffer[nidx] = false;
              x += c;
              y += row;
              ++n;
            } else {
              break;
            }
          }

          // scan row to the left
          for ( int c = col - 1; c >= 0; --c ) {
            const uint32_t nidx = ( c << 8 ) | row;
            // check up and down
            u_idx = nidx + 1;
            if ( row < VP::NRows - 1 && buffer[u_idx] ) {
              buffer[u_idx] = false;
              stack.push_back( u_idx );
            }
            d_idx = nidx - 1;
            if ( row > 0 && buffer[d_idx] ) {
              buffer[d_idx] = false;
              stack.push_back( d_idx );
            }
            // add set pixel to cluster or stop scanning
            if ( buffer[nidx] ) {
              buffer[nidx] = false;
              x += c;
              y += row;
              ++n;
            } else {
              break;
            }
          }
        } // while the stack is not empty

        // we are done with this cluster, calculate
        // centroid pixel coordinate and fractions.
        if ( n <= m_maxClusterSize ) {
          // if the pixel is smaller than the max cluster size, store it for the tracking
          const unsigned int cx = x / n;
          const auto         cy = Detector::VPChannelID::RowID{y / n};

          const auto chip = Detector::VPChannelID::ChipID{cx / CHIP_COLUMNS};
          const auto ccol = Detector::VPChannelID::ColumnID{cx % CHIP_COLUMNS};

          // store target (3D point for tracking)
          Detector::VPChannelID cid( sensor, chip, ccol, cy );

          const float fx = x / static_cast<float>( n ) - cx;
          const float fy = y / static_cast<float>( n ) - to_unsigned( cy );
          const float local_x = devp.local_x( cx ) + fx * devp.x_pitch( cx );
          const float local_y = ( to_unsigned( cy ) + 0.5f + fy ) * devp.pixel_size();

          // gx
          const float gx = ( ltg[0] * local_x + ltg[1] * local_y + ltg[9] );
          // gy
          const float gy = ( ltg[3] * local_x + ltg[4] * local_y + ltg[10] );
          // gz
          const float gz = ( ltg[6] * local_x + ltg[7] * local_y + ltg[11] );

          pool.emplace_back( 1, 1, gx, gy, gz, cid );
          ++offsets[module];
        }
      } // loop over all potential seed pixels
    } // loop over all banks

    std::partial_sum( offsets.begin(), offsets.end(), offsets.begin() );

    // sorting in phi for even modules
    auto cmp_phi_for_odd_modules = []( const VPLightCluster& a, const VPLightCluster& b ) {
      return ( a.y() < 0.f && b.y() > 0.f ) ||
             // same y side even and odd modules, check y1/x1 < y2/x2
             ( ( a.y() * b.y() ) > 0.f && ( a.y() * b.x() < b.y() * a.x() ) );
    };

    // sorting in phi for odd modules
    auto cmp_phi_for_even_modules = []( const VPLightCluster& a, const VPLightCluster& b ) {
      return ( a.y() > 0.f && b.y() < 0.f ) ||
             // same y side even and odd modules, check y1/x1 < y2/x2
             ( ( a.y() * b.y() ) > 0.f && ( a.y() * b.x() < b.y() * a.x() ) );
    };

    auto sort_module = [pool = std::ref( pool ), offsets = std::ref( offsets )]( auto id, auto cmp ) {
      std::sort( pool.get().begin() + offsets.get()[id], pool.get().begin() + offsets.get()[id + 1], cmp );
    };

    for ( size_t moduleID = 0; moduleID < VeloInfo::Numbers::NModules; ++moduleID ) {
      // if( msgLevel(MSG::DEBUG)){
      //   debug()<<"Sorting hits in moduleID by phi, usign x,y information "<<moduleID<<endmsg;
      // }
      // In even modules you fall in the branching at -180, 180 degrees, you want to do that continuos
      if ( moduleID % 2 == 1 ) {
        sort_module( moduleID, cmp_phi_for_odd_modules );
      } else {
        sort_module( moduleID, cmp_phi_for_even_modules );
      }
    }
    // if( msgLevel(MSG::DEBUG)){
    //   info()<<"Priting clusters information, checking sorting by atan2_approximation1, which is used in tracking
    //   !"<<endmsg; float phi; std::vector<float> phivec; for( size_t moduleID = 0 ; moduleID<
    //   VeloInfo::Numbers::NModules; ++moduleID){
    //     bool odd = moduleID%2 == 1;
    //     info()<<"Printing next module with :" << offsets[moduleID+1] - offsets[moduleID] <<" hits" << std::endl;
    //     for( size_t hit = offsets[moduleID]; hit< offsets[moduleID+1]; ++hit){
    //       phi = atan2_approximation1( pool[hit].y(), pool[hit].x() );
    //       phi = phi + ( phi<0. && !odd ? 360.f : 0.f);
    //       phivec.push_back( phi);
    //     }
    //     bool sortedbyphi = std::is_sorted( phivec.begin() + offsets[moduleID], phivec.begin() +offsets[moduleID+1]);
    //     if( !sortedbyphi){
    //       info()<<"Phi sorting OK"<<moduleID<<endmsg;
    //     }else{
    //       warning()<<"Phi sorting BUG "<<moduleID<<endmsg;
    //     }
    //   }
    //   //printout for comparing reference during dev ( old sorting by x )
    //   // debug()<<"printing next module with: " << offsets[i + 1] - offsets[i] << " hits" << endmsg;
    //   //   for ( size_t hit = offsets[i]; hit < offsets[i + 1]; ++hit ) {
    //   //     const auto& cl = pool[hit];
    //   //     std::cout << cl.x() << " " << cl.y() << " " << cl.z() << " " << cl.channelID().channelID() << std::endl;
    //   //   }
    //   //   //std::cout.precision(10);
    //   //   //for ( const auto& cl : pool )
    //   //   //std::cout << cl.x() << " " << cl.y() << " " << cl.z() << " " << cl.channelID().channelID() << std::endl;
    // }
    m_nbClustersCounter += pool.size();
    return result;
  }
} // namespace LHCb::Pr::Velo
