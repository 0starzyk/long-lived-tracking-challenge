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
/** @file TrackCloneData.h
 *
 * @author Wouter Hulsbergen
 * @date 2012-12-11
 */

#ifndef TRACKKERNEL_TRACKCLONEDATA_H
#define TRACKKERNEL_TRACKCLONEDATA_H

#include <array>
#include <cassert>
#include <memory>
#include <type_traits>

#include "Event/Track.h"
#include "Kernel/HitPattern.h"
#include "Kernel/LHCbID.h"
#include "LHCbMath/BloomFilter.h"
#include "Relations/Range.h"

namespace LHCb {
  namespace TrackCloneDataUtils {
    /// find number of common entries in sorted containers
    template <class Container>
    size_t nCommonEntries( const Container& lhs, const Container& rhs ) {
      size_t rc( 0 );
      auto   first1 = lhs.begin();
      auto   last1  = lhs.end();
      auto   first2 = rhs.begin();
      auto   last2  = rhs.end();
      while ( first1 != last1 && first2 != last2 ) {
        if ( *first1 < *first2 )
          ++first1;
        else if ( *first2 < *first1 )
          ++first2;
        else {
          ++first1;
          ++first2;
          ++rc;
        }
      }
      return rc;
    }

    /** @brief base class for LHCb::TrackCloneData (BloomFilter based)
     *
     * @author Wouter Hulsbergen
     * @date 2012-12-11
     * - Initial release
     *
     * @author Manuel Schiller
     * @data 2014-11-18
     * - factorise for better modularity
     *
     * This version of the class does all overlap calculations using only
     * BloomFilter estimates.
     *
     * Tuned for a false positive rate below 2.2 percent for 24 hits per
     * detector type or less
     */
    template <class BLOOMFILTER = BloomFilter<LHCb::LHCbID, 24, 22638, 1 << 20>>
    class TrackCloneDataBaseBloomSlice {
    public:
      /// hit types to distinguish
      enum HitType { // take care to put frequent comparisons first (cache!)
        T,           ///< T station hits
        VP,          ///< Velo hit
        UT,          ///< UT hits
        MaxType      ///< must be last
      };
      /// range of LHCbIDs of one type
      typedef Relations::Range_<LHCb::Track::LHCbIDContainer> LHCbIDs;
      /// per-type ranges of LHCbIDs
      typedef std::array<LHCbIDs, MaxType> PerTypeLHCbIDs;

      /// number of common hits
      size_t nCommon( const TrackCloneDataBaseBloomSlice& rhs, HitType type ) const {
        const auto overlap = ( m_bloomfilters[type] & rhs.m_bloomfilters[type] );
        if ( overlap.empty() ) return 0;
        return ( m_bloomfilters[type] & rhs.m_bloomfilters[type] ).size();
      }

      /// fraction of common hits
      double overlapFraction( const TrackCloneDataBaseBloomSlice& rhs, HitType type ) const {
        const auto overlap = ( m_bloomfilters[type] & rhs.m_bloomfilters[type] );
        if ( overlap.empty() ) return 0;
        size_t minsize = std::min( m_bloomfilters[type].size(), rhs.m_bloomfilters[type].size() );
        return ( minsize > 0 ) ? overlap.size() / double( minsize ) : 0;
      }

    protected:
      /** @brief Bloom filters per detector type for fast overlap checks
       */
      std::array<BLOOMFILTER, MaxType> m_bloomfilters;
    };

    /** @brief base class for LHCb::TrackCloneData (BloomFilter based)
     *
     * @author Wouter Hulsbergen
     * @date 2012-12-11
     * - Initial release
     *
     * @author Manuel Schiller
     * @data 2014-11-18
     * - factorise for better modularity
     *
     * This class uses BloomFilters to detect the very common no-overlap case,
     * but reverts back to determining the intersection of sorted lists of
     * LHCbIDs in case of (potential) overlaps
     *
     * Tuned for a false positive rate below 7.9 percent for 24 hits per
     * detector type or less
     */
    template <class BLOOMFILTER = BloomFilter<LHCb::LHCbID, 24, 81682, 1 << 20>>
    class TrackCloneDataBaseBloomSliceWithLHCbIDs
        : public LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSlice<BLOOMFILTER> {
    public:
      using typename LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSlice<BLOOMFILTER>::HitType;
      using typename LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSlice<BLOOMFILTER>::LHCbIDs;
      using typename LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSlice<BLOOMFILTER>::PerTypeLHCbIDs;
      /// LHCb IDs of a certain type
      const LHCbIDs& lhcbIDs( HitType type ) const { return m_ids[type]; }

      /// number of common hits
      size_t nCommon( const TrackCloneDataBaseBloomSliceWithLHCbIDs& rhs, HitType type ) const {
        const auto overlap( this->m_bloomfilters[type] & rhs.m_bloomfilters[type] );
        if ( overlap.empty() ) return 0;
        return nCommonEntries( m_ids[type], rhs.m_ids[type] );
      }

      /// fraction of common hits
      double overlapFraction( const TrackCloneDataBaseBloomSliceWithLHCbIDs& rhs, HitType type ) const {
        const auto overlap( this->m_bloomfilters[type] & rhs.m_bloomfilters[type] );
        if ( overlap.empty() ) return 0;
        size_t minsize = std::min( m_ids[type].size(), rhs.m_ids[type].size() );
        return ( minsize > 0 ) ? nCommonEntries( m_ids[type], rhs.m_ids[type] ) / double( minsize ) : 0;
      }

    protected:
      /// per hit type ranges of LHCbIDs
      PerTypeLHCbIDs m_ids;
    };

    /// help switch between TrackCloneDataBase and TrackCloneDataBaseWithLHCbIDs
    template <bool BloomFilterBasedOverlapCalc>
    class TrackCloneDataBaseBloomSliceSwitch;
    /// specialisation using TrackCloneDataBase
    template <>
    class TrackCloneDataBaseBloomSliceSwitch<true> : public TrackCloneDataBaseBloomSlice<> {};
    /// specialisation using TrackCloneDataBaseWithLHCbIDs
    template <>
    class TrackCloneDataBaseBloomSliceSwitch<false> : public TrackCloneDataBaseBloomSliceWithLHCbIDs<> {};

    /** @brief base class for LHCb::TrackCloneData ("track pointer slice")
     *
     * @author Wouter Hulsbergen
     * @date 2012-12-11
     * - Initial release
     *
     * @author Manuel Schiller
     * @data 2014-11-18
     * - factorise for better modularity
     *
     * this provides the "track pointer" based functionality of
     * TrackCloneData, and allows switching between unique_ptr (owning the
     * contained track), and plain LHCb::Track* (where any resource management
     * needs to be taken care of by the caller)
     */
    template <bool OwnTrackObjects>
    class TrackCloneDataBaseTrackPtrSlice {
    public:
      /// owning/non-owning track pointer type
      typedef typename std::conditional<OwnTrackObjects, std::unique_ptr<LHCb::Track>, LHCb::Track*>::type TrackPtr;

      /// return reference to track
      LHCb::Track& track() { return *m_track; }

      /// return reference to track
      const LHCb::Track& track() const { return *m_track; }

      /** @brief return TrackPtr to track (transfers ownership to caller if
       * OwnTrackObjects == true)
       */
      TrackPtr trackptr() && { return std::move( m_track ); }

    protected:
      /// pointer to track object
      TrackPtr m_track;

      /// constructor for use in derived classes
      TrackCloneDataBaseTrackPtrSlice( TrackPtr track ) : m_track( std::move( track ) ) {}
    };

    /** @brief base class for LHCb::TrackCloneData ("flags slice")
     *
     * @author Wouter Hulsbergen
     * @date 2012-12-11
     * - Initial release
     *
     * @author Manuel Schiller
     * @data 2014-11-18
     * - factorise for better modularity
     */
    class TrackCloneDataBaseFlagsSlice {
    private:
      enum { WeightBits = 28, UserFlagBits = 4 };

    public:
      /// type to record track weight
      typedef unsigned WeightType;

      /// sort by number of LHCbIDs
      bool operator<( const TrackCloneDataBaseFlagsSlice& rhs ) const { return m_weight > rhs.m_weight; }

      /// set the weight for ordering
      void setWeight( WeightType w ) {
        assert( w < WeightType( 1 << WeightBits ) );
        m_weight = w;
      }

      /// return the weight
      WeightType weight() const { return m_weight; }

    protected:
      /// weight of the track
      unsigned m_weight : WeightBits;
      /// user flags (derived classes are free to use these)
      unsigned m_userflags : UserFlagBits;

      /// default constructor for use in derived classes
      TrackCloneDataBaseFlagsSlice() : m_weight( 0 ), m_userflags( 0 ) {}

      /// return user flags
      unsigned userFlags() const { return m_userflags; }
      /// set user flags
      void setUserFlags( unsigned userflags ) {
        assert( userflags < ( 1u << UserFlagBits ) );
        m_userflags = userflags;
      }
    };
  } // namespace TrackCloneDataUtils

  /** @brief class to save some data for each track for use in clone killing
   *
   * @author Wouter Hulsbergen
   * @date 2012-12-11
   *  - Initial release
   *
   * @author Manuel Schiller <Manuel.Schiller@cern.ch>
   * @date 2014-10-17
   *  - Added BloomFilter to speed up the non-overlapping case
   *  - Move to ranges of LHCbIDs to reduce object size
   *  - Added option to do probabilistic overlap calculation based on
   *    BloomFilter alone
   * @author Manuel Schiller <Manuel.Schiller@cern.ch>
   * @date 2014-11-14
   *  - switch between BloomFilter-based and sorted-list based overlap
   *    calculations using template argument instead of preprocessor
   * @date 2014-11-18
   *    - modularise further, default to using unique_ptr for track
   *
   * @tparam OwnTrackObjects
   *    if true, TrackCloneData will take ownership of track pointer
   *    on construction using unique_ptr
   * @tparam BloomFilterBasedOverlapCalc
   *    if true, overlap calculations will be done based on the
   *    BloomFilters alone, the result is an estimate of the size of
   *    the overlap; if false, the class keeps the sorted ranges of
   *    the various LHCbIDs, and bases overlap calculations on that
   */
  template <bool OwnTrackObjects = true, bool BloomFilterBasedOverlapCalc = false>
  class TrackCloneData :
      // weights and userflags are the hottest objects in cache
      public LHCb::TrackCloneDataUtils::TrackCloneDataBaseFlagsSlice,
      public LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSliceSwitch<BloomFilterBasedOverlapCalc>,
      public LHCb::TrackCloneDataUtils::TrackCloneDataBaseTrackPtrSlice<OwnTrackObjects> {
  public:
    typedef typename LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSliceSwitch<BloomFilterBasedOverlapCalc>::HitType
        HitType;
    typedef typename LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSliceSwitch<BloomFilterBasedOverlapCalc>::LHCbIDs
        LHCbIDs;
    typedef typename LHCb::TrackCloneDataUtils::TrackCloneDataBaseBloomSliceSwitch<
        BloomFilterBasedOverlapCalc>::PerTypeLHCbIDs PerTypeLHCbIDs;
    typedef typename LHCb::TrackCloneDataUtils::TrackCloneDataBaseTrackPtrSlice<OwnTrackObjects>::TrackPtr TrackPtr;
    typedef LHCb::TrackCloneDataUtils::TrackCloneDataBaseFlagsSlice::WeightType                            WeightType;

    /// constructor
    TrackCloneData( TrackPtr track )
        : LHCb::TrackCloneDataUtils::TrackCloneDataBaseTrackPtrSlice<OwnTrackObjects>( std::move( track ) ) {
      doInit();
    }

  private:
    /// little helper for doInit() below which does the hard work
    void doInit( PerTypeLHCbIDs& pertypeids );
    /// little helper for the constructor (specialised on BloomFilterBasedOverlapCalc)
    template <class VOID = void>
    typename std::enable_if<BloomFilterBasedOverlapCalc, VOID>::type doInit() {
      PerTypeLHCbIDs ids;
      doInit( ids );
    }
    template <class VOID = void>
    typename std::enable_if<!BloomFilterBasedOverlapCalc, VOID>::type doInit() {
      doInit( this->m_ids );
    }
  };

  template <bool OwnTrackObjects, bool BloomFilterBasedOverlapCalc>
  void TrackCloneData<OwnTrackObjects, BloomFilterBasedOverlapCalc>::doInit( PerTypeLHCbIDs& pertypeids ) {
    const LHCb::Track&                  track( this->track() );
    const LHCb::Track::LHCbIDContainer& trackids( track.lhcbIDs() );
    for ( auto it = trackids.cbegin(); trackids.cend() != it; ++it ) {
      switch ( it->detectorType() ) {
      case LHCb::LHCbID::channelIDtype::VP:
        if ( it->isVP() ) {
          auto jt = it + 1;
          for ( ; trackids.cend() != jt && ( LHCb::LHCbID::channelIDtype::VP == jt->detectorType() ) && jt->isVP();
                ++jt )
            ;
          pertypeids[HitType::VP] = LHCbIDs( it, jt );
          this->m_bloomfilters[HitType::VP].insert( pertypeids[HitType::VP] );
          it = --jt;
        }
        break;
      case LHCb::LHCbID::channelIDtype::UT: {
        auto jt = it + 1;
        for ( ; trackids.cend() != jt && ( LHCb::LHCbID::channelIDtype::UT == jt->detectorType() ); ++jt )
          ;
        pertypeids[HitType::UT] = LHCbIDs( it, jt );
        this->m_bloomfilters[HitType::UT].insert( pertypeids[HitType::UT] );
        it = --jt;
      } break;
      case LHCb::LHCbID::channelIDtype::FT: {
        auto jt = it + 1;
        for ( ; trackids.cend() != jt && ( LHCb::LHCbID::channelIDtype::FT == jt->detectorType() ); ++jt )
          ;
        pertypeids[HitType::T] = LHCbIDs( it, jt );
        this->m_bloomfilters[HitType::T].insert( pertypeids[HitType::T] );
        it = --jt;
      } break;
      default:
        break;
      }
    }

    // compute a weight for sorting.
    LHCb::HitPattern hp( trackids );
    // first sort by type
    const int               fieldwidth         = 7; // wide enough to avoid overflow
    static const WeightType weightFromType[10] = {0, 11, 10, 20, 15, 16, 12, 2, 3};
    m_weight                                   = weightFromType[static_cast<int>( track.type() )] << ( 3 * fieldwidth );
    // then sort by the number of UT planes
    m_weight += hp.numUT() << ( 2 * fieldwidth );
    // then sort by number of T layers / velo _clusters_
    const WeightType numVeloLayers = hp.numVeloStations();
    m_weight += ( numVeloLayers + hp.numTLayers() ) << fieldwidth;
    // only finally sort by total number of hits
    m_weight += trackids.size();
  }
} // namespace LHCb

#endif // TRACKKERNEL_TRACKCLONEDATA_H
