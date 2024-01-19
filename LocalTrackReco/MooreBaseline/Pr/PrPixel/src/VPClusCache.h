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

#pragma once

#include <bitset>

namespace LHCb::Pr::Velo {
  namespace {
    //=========================================================================
    // Cache Super-Pixel cluster mask
    //=========================================================================
    // SP cluster mask for clustering, cached once.
    // There are 256 patterns and there can be at most two
    // distinct clusters in an SP.
    auto create_SPMasks() {
      std::array<uint8_t, 256> SPMasks;
      for ( int data = 0; data < 256; data++ ) {
        std::bitset<8> sp( data );
        if ( ( ( sp[1] || sp[5] ) && ( sp[2] || sp[6] ) ) || // 2 halfs are linked
             ( data & 0x33 ) == 0 || ( data & 0xCC ) == 0    // or one of the half is empty
        ) {
          SPMasks[data] = 0xFF; // Default mask for 1 cluster in SP
        } else {
          SPMasks[data] = 0x33; // Mask for half SP (2 clusters in SP)
        }
      }
      return SPMasks;
    }
    static const std::array<uint8_t, 256> s_SPMasks = create_SPMasks();

    //=========================================================================
    // Cache test for North-South connection test
    //=========================================================================
    auto create_linkNS() {
      std::array<uint8_t, 256> linkNS;
      for ( int data = 0; data < 256; data++ ) {
        std::bitset<8> sp( data );
        linkNS[data] = ( sp[4] && ( sp[0] || sp[1] ) ) || ( sp[5] && ( sp[0] || sp[1] || sp[2] ) ) ||
                       ( sp[6] && ( sp[1] || sp[2] || sp[3] ) ) || ( sp[7] && ( sp[2] || sp[3] ) );
      }
      return linkNS;
    }
    static const std::array<uint8_t, 256> s_linkNS = create_linkNS();

    //=========================================================================
    // Cache n, kx, ky for a given SP pattern
    //=========================================================================
    // n is the number of pixels
    // kx is the sum of pixel x
    // ky is the sum of pixel y
    auto create_SPn_kx_ky() {
      std::array<uint16_t, 256> SPn_kx_ky;
      for ( int data = 0; data < 256; data++ ) {
        std::bitset<8> sp( data );

        int cx0 = sp[0] + sp[1] + sp[2] + sp[3];
        int cx1 = sp[4] + sp[5] + sp[6] + sp[7];

        int cy1 = sp[1] + sp[5];
        int cy2 = sp[2] + sp[6];
        int cy3 = sp[3] + sp[7];

        uint8_t kx = cx1;
        uint8_t ky = cy1 + 2 * cy2 + 3 * cy3;
        uint8_t n  = cx0 + cx1;

        SPn_kx_ky[data] = ( ky << 8 ) | ( kx << 4 ) | n;
      }
      return SPn_kx_ky;
    }
    static const std::array<uint16_t, 256> s_SPn_kx_ky = create_SPn_kx_ky();

    constexpr auto SPn_kx_ky_getN( uint16_t n_kx_ky ) { return n_kx_ky & 0xF; }

    constexpr auto SPn_kx_ky_getKx( uint16_t n_kx_ky ) { return ( n_kx_ky >> 4 ) & 0xF; }

    constexpr auto SPn_kx_ky_getKy( uint16_t n_kx_ky ) { return n_kx_ky >> 8; }

    //=============================================================================
    // Link Horizontal : test if there is a link between two horizontal adjacent SP
    // 0 1 2 3 | 0 1 2 3
    // 4 5 6 7 | 4 5 6 7
    //=============================================================================
    constexpr auto linkH( uint8_t l, uint8_t r ) { return ( l & 0x88 ) && ( r & 0x11 ); }

    //=============================================================================
    // Link Diagonal Forward : test if there is a link between two SW-NE adjacent
    // SP
    //         | 0 1 2 3
    //         | 4 5 6 7
    // --------|--------
    // 0 1 2 3 |
    // 4 5 6 7 |
    //=============================================================================
    constexpr auto linkDF( uint8_t sw, uint8_t ne ) { return ( sw & 0x08 ) && ( ne & 0x10 ); }

    //=============================================================================
    // Link Diagonal Backward : test if there is a link between two NW-SE adjacent
    // SP
    // 0 1 2 3 |
    // 4 5 6 7 |
    // --------|--------
    //         | 0 1 2 3
    //         | 4 5 6 7
    //=============================================================================
    constexpr auto linkDB( uint8_t nw, uint8_t se ) { return ( nw & 0x80 ) && ( se & 0x01 ); }

    //=============================================================================
    // Link Vertical : test if there is a link between two vertical adjacent SP
    // 0 1 2 3
    // 4 5 6 7
    // -------
    // 0 1 2 3
    // 4 5 6 7
    //=============================================================================
    inline auto linkV( uint8_t b, uint8_t t ) { return s_linkNS[( b & 0xF0 ) | ( t & 0x0F )]; }

    //=============================================================================
    // Test for adjacency in 8-Connectivity
    // Assume that (xi * w + yi) is always > (xj * w + yj)
    //=============================================================================
    inline bool is_adjacent_8C_SP( int xi, int yi, int xj, int yj, int spi, int spj ) {
      if ( xi == xj ) {
        // if (yi == yj-1) return linkH(spi, spj); else
        if ( yi == yj + 1 ) return linkH( spj, spi );
      } else
          /*if (xi == xj-1) {
            if (yi == yj)   return linkV(spi, spj); else
            if (yi == yj-1) return linkDF(spj, spi); else
            if (yi == yj+1) return linkDB(spi, spj);
          } else*/
          if ( xi == xj + 1 ) {
        if ( yi == yj )
          return linkV( spj, spi );
        else if ( yi == yj - 1 )
          return linkDF( spi, spj );
        else if ( yi == yj + 1 )
          return linkDB( spj, spi );
      }
      return false;
    }

    //=============================================================================
    // Helper functions to manipulate SP :
    //=============================================================================

    constexpr auto SP_getPixels( uint32_t sp ) { return sp & 0xFFU; }

    constexpr auto SP_getRow( uint32_t sp ) { return ( sp >> 8 ) & 0x3FU; }

    constexpr auto SP_getCol( uint32_t sp ) { return ( sp >> 14 ) & 0x1FFU; }

    constexpr auto SP_isIsolated( uint32_t sp ) { return sp & 0x80000000U; }

  } // namespace
} // namespace LHCb::Pr::Velo
