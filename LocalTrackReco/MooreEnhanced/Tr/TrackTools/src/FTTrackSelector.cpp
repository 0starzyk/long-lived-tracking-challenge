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
/** @file FTTrackSelector.cpp
 *
 *  Implementation file for reconstruction tool FTTrackSelector
 *
 *  @author Sophie Hollitt
 */

#include "Detector/FT/FTConstants.h"
#include "TrackSelector.h"

class FTTrackSelector : public TrackSelector {

public:
  /// Constructor
  using TrackSelector::TrackSelector;

  /** Returns if the given track is selected or not
   *
   *  @param track Reference to the track to test
   *
   *  @return boolean indicating if the track is selected or not
   *  @retval true  Track is selected
   *  @retval false Track is rejected
   */
  bool accept( const LHCb::Track& track ) const override;

private:
  Gaudi::Property<size_t> m_minHitsCSide{this, "MinHitsCSide", 0};
  Gaudi::Property<size_t> m_minHitsASide{this, "MinHitsASide", 0};
  Gaudi::Property<size_t> m_minHitsTopHalf{this, "MinHitsTopHalf", 0};
  Gaudi::Property<size_t> m_minHitsBottomHalf{this, "MinHitsBottomHalf", 0};
  Gaudi::Property<size_t> m_minHits{this, "MinHits", 0};
  Gaudi::Property<size_t> m_nExcludedChannelsFromMatEdge{this, "NExcludedChannelsFromMatEdge", 0};
  Gaudi::Property<size_t> m_nExcludedChannelsFromSiPMEdge{this, "NExcludedChannelsFromSiPMEdge", 0};
  Gaudi::Property<size_t> m_nExcludedChannelsFromSiPMCenter{this, "NExcludedChannelsFromSiPMCenter", 0};
};

DECLARE_COMPONENT( FTTrackSelector )

//-----------------------------------------------------------------------------

bool FTTrackSelector::accept( const LHCb::Track& track ) const {
  std::array<size_t, 5> numHits{};
  const unsigned int    maxChannels = LHCb::Detector::FT::nChannels;
  const unsigned int    maxInDie    = maxChannels / 2; // FIXME: collect value from Detector
  for ( const LHCb::LHCbID lhcbid : track.lhcbIDs() ) {
    if ( !lhcbid.isFT() ) continue;
    const LHCb::Detector::FTChannelID ftid = lhcbid.ftID();
    const unsigned int                side = ftid.getSide();

    // Check position of hits for excluded channels
    if ( m_nExcludedChannelsFromSiPMEdge > 0 || m_nExcludedChannelsFromMatEdge > 0 ||
         m_nExcludedChannelsFromSiPMCenter ) {
      const unsigned int localChannelIdx = ftid.channel();
      const unsigned int localSiPMIdx    = ftid.sipm();

      // remove all of outer SiPM if required
      if ( m_nExcludedChannelsFromMatEdge >= maxChannels ) {
        if ( localSiPMIdx == 0 || localSiPMIdx == 3 ) {
          return false;
        }
        // remove remaining channels from inner SiPMs
        else if ( m_nExcludedChannelsFromMatEdge > maxChannels ) {
          const unsigned int localExcludedChannels = m_nExcludedChannelsFromMatEdge - maxChannels;
          if ( localSiPMIdx == 1 && localChannelIdx < localExcludedChannels ) {
            return false;
          } else if ( localSiPMIdx == 2 && localChannelIdx >= maxChannels - localExcludedChannels ) {
            return false;
          }
        }
      }
      // remove channels on outer SiPMs
      else if ( m_nExcludedChannelsFromMatEdge > 0 ) {
        if ( localSiPMIdx == 0 && localChannelIdx < m_nExcludedChannelsFromMatEdge ) {
          return false;
        } else if ( localSiPMIdx == 3 && localChannelIdx >= maxChannels - m_nExcludedChannelsFromMatEdge ) {
          return false;
        }
      }

      // remove outer channels on all SiPMs
      if ( m_nExcludedChannelsFromSiPMEdge > 0 ) {
        if ( localChannelIdx < m_nExcludedChannelsFromSiPMEdge ||
             localChannelIdx >= maxChannels - m_nExcludedChannelsFromSiPMEdge ) {
          return false;
        };
      };

      if ( m_nExcludedChannelsFromSiPMCenter > 0 ) {
        if ( localChannelIdx < maxInDie && localChannelIdx >= maxInDie - m_nExcludedChannelsFromSiPMCenter ) {
          return false;
        } else if ( localChannelIdx >= maxInDie && localChannelIdx - maxInDie < m_nExcludedChannelsFromSiPMCenter ) {
          return false;
        };
      };
    };

    // Count hits for SciFi region selection
    ++numHits[side];
    if ( ftid.isTop() ) ++numHits[2];
    if ( ftid.isBottom() ) ++numHits[3];
    ++numHits[4];
  }
  return numHits[0] >= m_minHitsCSide && numHits[1] >= m_minHitsASide && numHits[2] >= m_minHitsTopHalf &&
         numHits[3] >= m_minHitsBottomHalf && numHits[4] >= m_minHits && TrackSelector::accept( track );
}
