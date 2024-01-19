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
#include "Event/StateParameters.h"
#include "Event/Track.h"
#include "Kernel/VPConstants.h"
#include "VPDAQ/VPRetinaClusterConstants.h"
#include "VPKernel/PixelUtils.h"

// Local
#include "VPClusFull.h"
#include <VPDet/VPDetPaths.h>
#include <iomanip>

DECLARE_COMPONENT_WITH_ID( LHCb::Pr::Velo::VPClusFull, "VPClusFull" )

namespace LHCb::Pr::Velo {
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

    // SP pattern buffers for clustering, cached once.
    // There are 256 patterns and there can be at most two
    // distinct clusters in an SP.
    static const std::array<SPCache, 256> s_SPCaches = create_SPPatterns();

    using namespace Pixel;
  } // namespace
  //=============================================================================
  // Standard constructor, initializes variables
  //=============================================================================
  VPClusFull::VPClusFull( const std::string& name, ISvcLocator* pSvcLocator )
      : MultiTransformer( name, pSvcLocator,
                          {KeyValue{"RawEventLocation",
                                    Gaudi::Functional::concat_alternatives(
                                        RawEventLocation::Velo, RawEventLocation::Default, RawEventLocation::Other )},
                           KeyValue{"DEVP", LHCb::Det::VP::det_path}},
                          {KeyValue{"ClusterLocation", VPFullClusterLocation::Default},
                           KeyValue{"ClusterOffsets", VPFullClusterLocation::Offsets}} ) {}

  //=============================================================================
  // Main execution
  //=============================================================================
  std::tuple<std::vector<VPFullCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>> VPClusFull::
                                                                                            operator()( const EventContext&, const RawEvent& rawEvent, const DeVP& devp ) const {
    // WARNING:
    // This clustering algorithm is designed to perform the Offline Clustering without checking any isolation flag
    auto result = std::tuple<std::vector<VPFullCluster>, std::array<unsigned, VeloInfo::Numbers::NOffsets>>{};

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
    std::array<bool, VP::NPixelsPerSensor> buffer{};
    std::vector<uint32_t>                  pixel_idx;
    std::vector<uint32_t>                  stack;

    // reserve a minimal stack
    stack.reserve( 64 );

    // Since the pool is local, to first preallocate the pool, then count hits per module,
    // and then preallocate per module and move hits might not be faster than adding
    // directly to the PixelModuleHits (which would require more allocations, but
    // not many if we start with a sensible default)
    auto& [pool, offsets] = result;
    const unsigned int startSize = 10000U;
    pool.reserve( startSize );

    // 1 module = N Cluster -> N x M channelIDs
    std::vector<std::vector<Detector::VPChannelID>> channelIDs;
    channelIDs.reserve( 400 );

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
      channelIDs.clear();
      // reset and then fill the super pixel buffer for a sensor
      // memset(m_sp_buffer,0,256*256*3*sizeof(unsigned char));
      // the memset is too slow here. the occupancy is so low that
      // resetting a few elements is *much* faster.
      const unsigned int nrc = pixel_idx.size();
      for ( unsigned int irc = 0; irc < nrc; ++irc ) { buffer[pixel_idx[irc]] = 0; }
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

        if ( sp_row > ( VP::NRows / 4 - 1 ) ) continue; // protect against super pixels outside sensor coordinates.
        if ( sp_col > ( VP::NColumns * VP::NChipsPerSensor / 2 - 1 ) )
          continue; // protect against super pixels outside sensor coordinates.

        // const uint32_t no_sp_neighbours = sp_word & 0x80000000U;
        // Skip the check of no_sp_neighbours here, we want to record all pixels here for MCLinking purposes [offline
        // Mode, full pixels recording]

        // Record all pixels
        for ( uint32_t shift = 0; shift < 8; ++shift ) {
          const uint8_t pixel = sp & 1;
          if ( pixel ) {
            const uint32_t row = sp_row * 4 + shift % 4;
            const uint32_t col = sp_col * 2 + shift / 4;
            const uint32_t idx = ( col << 8 ) | row;
            // std::cout << (int)pixel << std::endl;
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
        const uint8_t  pixel = buffer[idx];
        if ( 0 == pixel ) continue; // pixel is used in another cluster

        // 8-way row scan optimized seeded flood fill from here.
        stack.clear();

        // mark seed as used
        buffer[idx] = 0;

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
          const auto     row = Detector::VPChannelID::RowID{idx & 0xFFu};
          const uint32_t col = ( idx >> 8 ) & 0x3FFu;
          x += col;
          y += to_unsigned( row );
          ++n;

          const auto            chip = Detector::VPChannelID::ChipID( col / CHIP_COLUMNS );
          Detector::VPChannelID cid( sensor, chip, Detector::VPChannelID::ColumnID( col % CHIP_COLUMNS ), row );
          if ( n == 1 ) {
            channelIDs.push_back( {{cid}} );
          } else {
            channelIDs.back().push_back( cid );
          }
          // check up and down
          uint32_t u_idx = idx + 1;
          if ( to_unsigned( row ) < VP::NRows - 1 && buffer[u_idx] ) {
            buffer[u_idx] = 0;
            stack.push_back( u_idx );
          }
          uint32_t d_idx = idx - 1;
          if ( to_unsigned( row ) > 0 && buffer[d_idx] ) {
            buffer[d_idx] = 0;
            stack.push_back( d_idx );
          }

          // scan row to the right
          for ( unsigned int c = col + 1; c < VP::NSensorColumns; ++c ) {
            const uint32_t nidx = ( c << 8 ) | to_unsigned( row );
            // check up and down
            u_idx = nidx + 1;
            if ( to_unsigned( row ) < VP::NRows - 1 && buffer[u_idx] ) {
              buffer[u_idx] = 0;
              stack.push_back( u_idx );
            }
            d_idx = nidx - 1;
            if ( to_unsigned( row ) > 0 && buffer[d_idx] ) {
              buffer[d_idx] = 0;
              stack.push_back( d_idx );
            }
            // add set pixel to cluster or stop scanning
            if ( buffer[nidx] ) {
              buffer[nidx] = 0;
              x += c;
              y += to_unsigned( row );
              ++n;
              const auto chip = Detector::VPChannelID::ChipID( c / CHIP_COLUMNS );
              channelIDs.back().emplace_back( sensor, chip, Detector::VPChannelID::ColumnID( c % CHIP_COLUMNS ), row );
            } else {
              break;
            }
          }

          // scan row to the left
          for ( int c = col - 1; c >= 0; --c ) {
            const uint32_t nidx = ( c << 8 ) | to_unsigned( row );
            // check up and down
            u_idx = nidx + 1;
            if ( to_unsigned( row ) < VP::NRows - 1 && buffer[u_idx] ) {
              buffer[u_idx] = 0;
              stack.push_back( u_idx );
            }
            d_idx = nidx - 1;
            if ( to_unsigned( row ) > 0 && buffer[d_idx] ) {
              buffer[d_idx] = 0;
              stack.push_back( d_idx );
            }
            // add set pixel to cluster or stop scanning
            if ( buffer[nidx] ) {
              buffer[nidx] = 0;
              x += c;
              y += to_unsigned( row );
              ++n;
              const auto chip = Detector::VPChannelID::ChipID( c / CHIP_COLUMNS );
              channelIDs.back().emplace_back( sensor, chip, Detector::VPChannelID::ColumnID( c % CHIP_COLUMNS ), row );
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
          pool.emplace_back( fx, fy, gx, gy, gz, cid, std::move( channelIDs.back() ) );
          ++offsets[module];
        }
      } // loop over all potential seed pixels
    } // loop over all banks

    std::partial_sum( offsets.begin(), offsets.end(), offsets.begin() );
    // Do we actually need to sort the hits ? [ depends, if the offline clustering will be re-run and tracking will use
    // those clusters , yes, otherwise no ]
    for ( size_t moduleID = 0; moduleID < VeloInfo::Numbers::NModules; ++moduleID ) {
      // In even modules you fall in the branching at -180, 180 degrees, you want to do that continuos
      std::stable_sort(
          pool.begin() + offsets[moduleID], pool.begin() + offsets[moduleID + 1],
          []( const VPFullCluster& a, const VPFullCluster& b ) { return a.channelID() < b.channelID(); } );
    }
    if ( msgLevel( MSG::DEBUG ) ) {
      for ( auto& cl : pool ) {
        info() << "----" << endmsg;
        info() << cl << endmsg;
        info() << " [fx,fy] " << cl.xfraction() << "," << cl.yfraction() << endmsg;
        info() << " [x,y,z] " << cl.x() << "," << cl.y() << "," << cl.z() << endmsg;
        info() << "pixels" << endmsg;
        for ( auto& pixel : cl.pixels() ) { info() << "\t" << pixel << endmsg; }
      }
    }
    if ( msgLevel( MSG::DEBUG ) ) { info() << "N VPFullCluster produced :  " << pool.size() << endmsg; }
    return result;
  }
} // namespace LHCb::Pr::Velo
