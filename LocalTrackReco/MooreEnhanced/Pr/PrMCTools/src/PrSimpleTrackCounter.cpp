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

/** @class PrSimpleTrackCounter PrSimpleTrackCounter.cpp
 *
 *  This algorithm demonstrates how the LinksByKey created by PrTrackAssociator is used
 *  together with the corresponding track and MC particles containers to obtain very basic
 *  efficiency information and to select the MC particle which matches with the highest weight.
 *
 *  The class is templated on the input track type. It is mandatory that the same identifier is used
 *  when reading the LinksByKey object created by PrTrackAssociatior. For KeyedContainers the key of an object
 *  is used as identifier. For plain vectors the index is used.
 *
 *  Original author Sascha Stahl, see revision history for other contributors:
 *  @author Sascha Stahl
 *  @date   2019-07-15
 */

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/Track.h"

#include "LHCbAlgs/Consumer.h"

namespace {
  using KeyValue = std::pair<std::string, std::string>;
}

template <typename Tracks>
class PrSimpleTrackCounter
    : public LHCb::Algorithm::Consumer<void( Tracks const&, LHCb::MCParticles const&, LHCb::LinksByKey const& )> {
public:
  PrSimpleTrackCounter( const std::string& name, ISvcLocator* pSvcLocator )
      : LHCb::Algorithm::Consumer<void( Tracks const&, LHCb::MCParticles const&, LHCb::LinksByKey const& )>(
            name, pSvcLocator,
            {KeyValue{"Tracks", ""}, KeyValue{"MCParticles", LHCb::MCParticleLocation::Default},
             KeyValue{"Links", ""}} ){};

  void operator()( Tracks const& tracks, LHCb::MCParticles const& mc_particles,
                   LHCb::LinksByKey const& links ) const override {
    unsigned int index{0};
    for ( auto const& track : tracks ) {
      unsigned int key{0};
      if constexpr ( std::is_same_v<Tracks, LHCb::Track::Range> ) {
        key = track->index();
      } else {
        key = index;
      }

      LHCb::MCParticle const* mcparticle{nullptr};
      float                   max_weight{0};
      unsigned int            n_mcparticles{0};
      links.applyToLinks( key, [&n_mcparticles, &max_weight, &mcparticle,
                                &mc_particles]( unsigned int /* srcIndex */, unsigned int mcPartKey, float weight ) {
        n_mcparticles++;
        if ( weight > max_weight ) {
          max_weight = weight;
          mcparticle = static_cast<LHCb::MCParticle const*>( mc_particles.containedObject( mcPartKey ) );
        }
      } );
      bool is_associated = ( mcparticle != nullptr );
      m_efficiency += is_associated;
      if ( is_associated ) {
        m_mcparticles_per_track += n_mcparticles;
        m_weight_per_track += max_weight;
      }

      index++;
    }
  };

private:
  mutable Gaudi::Accumulators::BinomialCounter<>  m_efficiency{this, "Efficiency"};
  mutable Gaudi::Accumulators::AveragingCounter<> m_mcparticles_per_track{this, "MC particles per track"};
  mutable Gaudi::Accumulators::AveragingCounter<> m_weight_per_track{this, "Max weight per track"};
};

DECLARE_COMPONENT_WITH_ID( PrSimpleTrackCounter<LHCb::Track::Range>, "PrSimpleTrackCounter" )
DECLARE_COMPONENT_WITH_ID( PrSimpleTrackCounter<std::vector<LHCb::Event::v2::Track>>, "PrSimpleV2TrackCounter" )
