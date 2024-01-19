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
#include "Event/ZipUtils.h"
#include "Kernel/EventLocalAllocator.h"
#include "Kernel/STLExtensions.h"
#include "PrKernel/UTHit.h"
#include "UTDAQ/UTInfo.h"
#include <cassert>
#include <cstdint>
#include <vector>

/**
 *  UTHitHandler contains the hits in the UT detector and the accessor to them
 */
namespace UT {
  class HitHandler {
  public:
    /// Allocator type
    using allocator_type = LHCb::Allocators::EventLocal<UT::Hit>;

  private:
    /// Type for hit storage
    using HitVector = std::vector<UT::Hit, allocator_type>;

    /// Internal indices storage for ranges
    using HitIndices = std::pair<std::size_t, std::size_t>;

  public:
    class HitsInUT {
      std::array<HitIndices, static_cast<int>( UTInfo::DetectorNumbers::Stations ) *
                                 static_cast<int>( UTInfo::DetectorNumbers::Layers ) *
                                 static_cast<int>( UTInfo::DetectorNumbers::Regions ) *
                                 static_cast<int>( UTInfo::DetectorNumbers::Sectors )>
          m_data = {};

    public:
      constexpr static int idx( unsigned int station, unsigned int layer, unsigned int region, unsigned int sector ) {
        assert( station != 0 && station <= static_cast<int>( UTInfo::DetectorNumbers::Stations ) );
        assert( layer != 0 && layer <= static_cast<int>( UTInfo::DetectorNumbers::Layers ) );
        assert( region != 0 && region <= static_cast<int>( UTInfo::DetectorNumbers::Regions ) );
        assert( sector != 0 && sector <= static_cast<int>( UTInfo::DetectorNumbers::Sectors ) );
        return ( ( ( station - 1 ) * static_cast<int>( UTInfo::DetectorNumbers::Layers ) + ( layer - 1 ) ) *
                     static_cast<int>( UTInfo::DetectorNumbers::Regions ) +
                 ( region - 1 ) ) *
                   static_cast<int>( UTInfo::DetectorNumbers::Sectors ) +
               ( sector - 1 );
      }
      const HitIndices& operator()( unsigned int station, unsigned int layer, unsigned int region,
                                    unsigned int sector ) const {
        return m_data[idx( station, layer, region, sector )];
      }
      HitIndices& operator()( unsigned int station, unsigned int layer, unsigned int region, unsigned int sector ) {
        return m_data[idx( station, layer, region, sector )];
      }
      const HitIndices& operator()( unsigned int idx ) const { return m_data[idx]; }
      HitIndices&       operator()( unsigned int idx ) { return m_data[idx]; }
    };

    using HitRange = LHCb::span<const UT::Hit>;

    HitHandler( Zipping::ZipFamilyNumber, allocator_type alloc = {} ) : m_allhits{alloc} {}
    HitHandler( allocator_type alloc = {} ) : m_allhits{alloc} {}

    UT::Hit& emplace_back( const DeUTSector& aSector, unsigned int fullChanIdx, unsigned int strip, double fracStrip,
                           LHCb::Detector::UT::ChannelID chanID, unsigned int size, bool highThreshold ) {
      double dxDy{0};
      double dzDy{0};
      double xAtYEq0{0};
      double zAtYEq0{0};
      double yBegin{0};
      double yEnd{0};
      //--- this method allow to set the values
      const auto fracStripOvfour = fracStrip / 4;
      aSector.trajectory( strip, fracStripOvfour, dxDy, dzDy, xAtYEq0, zAtYEq0, yBegin, yEnd );
      const auto cos   = aSector.cosAngle();
      const auto error = aSector.pitch() / std::sqrt( 12.0 );

      // NB : The way the indices are setup here assumes all hits for a given
      //      station, layer, region and sector come in order, which appears
      //      to be the case. This must remain so...
      //
      //      Currently, what I am seeing from the MC has this sorting.
      //      But this would also need to be the case for real data.

      // get the indices for this region
      auto& indices = m_indices( fullChanIdx );

      // if first for this range, set the begin and end indices
      if ( &indices != m_last_indices ) {
        // check to see if this range has been filled previously.
        // If it has, assumed ordering is broken
        assert( indices.first == indices.second );
        // reset indices to current end of container
        indices = {m_allhits.size(), m_allhits.size()};
        // update used last index cache
        m_last_indices = &indices;
      }

      // add a new hit
      auto& hit = m_allhits.emplace_back( chanID, size, highThreshold, dxDy, xAtYEq0, zAtYEq0, yBegin, yEnd, cos, error,
                                          strip, fracStripOvfour );
      // increment the end index for current range
      ++indices.second;
      return hit;
    }

    /// Reserve size in the overall hit container
    void reserve( const std::size_t nHits ) { m_allhits.reserve( nHits ); }

    /// Access the range for a given set of hits
    [[nodiscard]] HitRange hits( unsigned int station, unsigned int layer, unsigned int region,
                                 unsigned int sector ) const noexcept {
      const auto& indices = m_indices( station, layer, region, sector );
      return LHCb::make_span( m_allhits.begin() + indices.first, m_allhits.begin() + indices.second );
    }

    [[nodiscard]] const auto& indices( const int fullChanIdx ) const { return m_indices( fullChanIdx ); }
    [[nodiscard]] const auto& hits() const { return m_allhits; }

    [[nodiscard]] const UT::Hit* hit( LHCb::Detector::UT::ChannelID id ) const noexcept {
      auto hs  = hits( id.station(), id.layer(), id.detRegion(), id.sector() );
      auto hit = std::find_if( hs.begin(), hs.end(), [id]( UT::Hit const& hit ) { return hit.chanID() == id; } );
      return hit != hs.end() ? &*hit : nullptr;
    }

    /// get the total number of hits
    [[nodiscard]] auto nbHits() const noexcept { return m_allhits.size(); }

    HitHandler& clear() {
      m_indices = {};
      m_allhits.clear();
      m_last_indices = nullptr;
      return *this;
    }

  private:
    // Indices for each range
    HitsInUT m_indices;

    // single vector of all hits
    HitVector m_allhits;

    // cache pointer to last indices used
    HitIndices* m_last_indices = nullptr;
  };
} // namespace UT
