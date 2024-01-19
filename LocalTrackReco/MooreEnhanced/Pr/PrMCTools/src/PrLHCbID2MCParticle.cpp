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

#include "Event/FTLiteCluster.h"
#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/PrHits.h"
#include "Event/UTCluster.h"
#include "Event/VPFullCluster.h"
#include "Event/VPLightCluster.h"
#include "PrKernel/UTHitHandler.h"

#include "LHCbAlgs/Transformer.h"

#include <boost/numeric/conversion/cast.hpp>

namespace {

  template <typename ContainerType>
  std::string getBaseName();

  template <>
  std::string getBaseName<std::vector<LHCb::VPFullCluster>>() {
    return "VPFullClusters";
  }

  template <>
  std::string getBaseName<std::vector<LHCb::VPLightCluster>>() {
    return "VPLightCluster";
  }

  template <>
  std::string getBaseName<UT::HitHandler>() {
    return "UTHits";
  }

  template <>
  std::string getBaseName<LHCb::FTLiteCluster::FTLiteClusters>() {
    return "FTLiteClusters";
  }

  template <>
  std::string getBaseName<MuonHitContainer>() {
    return "MuonHits";
  }

  /// build a Location string for a given type
  template <typename ContainerType>
  std::string getLocation() {
    return getBaseName<ContainerType>() + "Location";
  }

  /// build a Link Location string for a given type
  template <typename ContainerType>
  std::string getLinkLocation() {
    return getBaseName<ContainerType>() + "LinkLocation";
  }

  /**
   * type used to replace LinksByKey in template parameters pack expansion
   * as it needs to depend on a template argument ot be expanded
   */
  template <typename Container>
  using LinksByKeyT = LHCb::LinksByKey;

  using CounterType  = Gaudi::Accumulators::SummingCounter<>;
  using BufferedType = std::decay_t<decltype( std::declval<CounterType>().buffer() )>;

  /// link all particles to the specified id
  void linkAll( const LHCb::LinksByKey& ilink, LHCb::LinksByKey& olink, const LHCb::MCParticles& mcParts,
                LHCb::LHCbID id, const std::vector<unsigned int>& ids, BufferedType& nullMCParticles_counter ) {
    std::vector<const LHCb::MCParticle*> partList;
    partList.reserve( ids.size() );
    for ( auto subID : ids ) {
      ilink.applyToLinks( subID, [&mcParts, &partList]( unsigned int, unsigned int tgtIndex, float ) {
        const LHCb::MCParticle* mcPart = static_cast<const LHCb::MCParticle*>( mcParts.containedObject( tgtIndex ) );
        partList.push_back( mcPart );
      } );
    }

    // Remove any null entries -- the muon linking seems to produce a few of these...
    {
      auto new_end = std::remove( partList.begin(), partList.end(), nullptr );
      nullMCParticles_counter += std::distance( new_end, partList.end() );
      partList.erase( new_end, partList.end() );
    }

    // SORTING:
    // THe access pattern in the ilink is determined by the sorting of the underlying Clusters
    // For the Velo we do offline Clustering storing vector< VPFullCluster > which can be in a different order w.r.t the
    // tracking clusters The underlying linker of lhcbID 2 MCParticle gets filled in a different order using weights.
    // For equal weights the access patters is first in -> first out (from the linker ) We need to sort this parList
    // again to have a stable PrTrackAssociator algorithm using the linker produced in this algorithm
    std::sort( partList.begin(), partList.end(), []( const LHCb::MCParticle* a, const LHCb::MCParticle* b ) {
      if ( a->key() != b->key() ) {
        return a->key() < b->key();
      } else {
        return a->pt() > b->pt();
      }
    } );
    // remove any possible duplicates
    partList.erase( std::unique( partList.begin(), partList.end() ), partList.end() );
    for ( const auto& part : partList ) { olink.link( id.lhcbID(), mcParts( part->index() ), 1.0 ); }
  }

  /**
   * generic functor handling set of clusters and for each set, looping over clusters and calling linkAll
   * for each individual cluster
   * This instantiation is the one looping over the set of clusters and corresponding set of links via
   * template recursion. When a single cluster and link are used, specializations are available
   */
  template <typename ClusterContainer, typename... OtherContainers>
  struct LoopOverClusters {
    void operator()( const ClusterContainer& clusters, const OtherContainers&... otherClusters,
                     const LHCb::LinksByKey& clustersLink, const LinksByKeyT<OtherContainers>&... otherLinks,
                     LHCb::LinksByKey& idLink, const LHCb::MCParticles& mcParts,
                     BufferedType& nullMCParticles_counter ) {
      LoopOverClusters<ClusterContainer>{}( clusters, clustersLink, idLink, mcParts, nullMCParticles_counter );
      LoopOverClusters<OtherContainers...>{}( otherClusters..., otherLinks..., idLink, mcParts,
                                              nullMCParticles_counter );
    }
  };

  /// generic method looping over a set of clusters and calling linkAll for each of them
  template <typename ClusterContainer>
  struct LoopOverClusters<ClusterContainer> {
    void operator()( const ClusterContainer& clusters, const LHCb::LinksByKey& clusterLink, LHCb::LinksByKey& idLink,
                     const LHCb::MCParticles& mcParts, BufferedType& nullMCParticles_counter ) {
      for ( const auto& clus : clusters ) {
        linkAll( clusterLink, idLink, mcParts, clus.channelID(), {clus.channelID().channelID()},
                 nullMCParticles_counter );
      }
    }
  };

  template <>
  struct LoopOverClusters<UT::HitHandler> {
    void operator()( const UT::HitHandler& hitHandler, const LHCb::LinksByKey& clusterLink, LHCb::LinksByKey& idLink,
                     const LHCb::MCParticles& mcParts, BufferedType& nullMCParticles_counter ) {
      for ( unsigned int station = 1; station < 3; station++ ) {
        for ( unsigned int layer = 1; layer < 3; layer++ ) {
          for ( unsigned int region = 1; region < 4; region++ ) {
            for ( unsigned int sector = 1; sector < 99; sector++ ) {
              for ( const auto& hit : hitHandler.hits( station, layer, region, sector ) ) {
                linkAll( clusterLink, idLink, mcParts, hit.chanID(),
                         {boost::numeric_cast<unsigned int>( hit.chanID().channelID() )}, nullMCParticles_counter );
              }
            }
          }
        }
      }
    }
  };

  template <>
  struct LoopOverClusters<LHCb::FTLiteCluster::FTLiteClusters> {
    void operator()( const LHCb::FTLiteCluster::FTLiteClusters& clusters, const LHCb::LinksByKey& clusterLink,
                     LHCb::LinksByKey& idLink, const LHCb::MCParticles& mcParts,
                     BufferedType& nullMCParticles_counter ) {
      for ( const auto& clus : clusters.range() ) {
        linkAll( clusterLink, idLink, mcParts, clus.channelID(), {clus.channelID().channelID()},
                 nullMCParticles_counter );
      }
    }
  };

  template <>
  struct LoopOverClusters<MuonHitContainer> {
    void operator()( MuonHitContainer const& muonHits, LHCb::LinksByKey const& muonDigitsLink, LHCb::LinksByKey& idLink,
                     LHCb::MCParticles const& mcParts, BufferedType& nullMCParticles_counter ) {
      std::vector<unsigned int> sub_ids;
      for ( auto n_station = 0; n_station < 4; ++n_station ) {
        auto hits = muonHits.hits( n_station );
        for ( auto const& hit : hits ) {
          sub_ids.clear();
          for ( auto const& subtile : hit.subtiles() ) { sub_ids.push_back( static_cast<unsigned long>( subtile ) ); }
          linkAll( muonDigitsLink, idLink, mcParts, hit.tile(), sub_ids, nullMCParticles_counter );
        }
      }
    }
  };
} // namespace

/**
 * This algorithm takes a set of cluster/hit containers and corresponding LinksByKey relations
 * linking them to MCParticles to build another LinksByKey relation linking all LHcbIds of all
 * input containers to their MCParticle
 * Note that it is templated on the input containers so that it can work with any number of
 * containers. Typical usages are with one container (e.g. only velo clusters) and up to 3
 * (velo clusters, UT hits and FT light clusters).
 * The way it handles the unknown number of containers is by calling the LoopOverClusters
 * functor that is templated to take any number of inputs and uses template recursion to
 * handle its inputs one through specializations for the one input case.
 * In order to be able to use this class in python, it needs to be instantiated using
 * DECLARE_COMPONENT_WITH_ID (or DECLARE_COMPONENT but then its name is ugly). This thus
 * needs to be done for each possible instantiation of it.
 */
template <typename... ContainerType>
class PrLHCbID2MCParticle
    : public LHCb::Algorithm::Transformer<LHCb::LinksByKey( const LHCb::MCParticles&, const ContainerType&...,
                                                            const LinksByKeyT<ContainerType>&... )> {

public:
  using Transformer = LHCb::Algorithm::Transformer<LHCb::LinksByKey( const LHCb::MCParticles&, const ContainerType&...,
                                                                     const LinksByKeyT<ContainerType>&... )>;
  using KeyValue    = typename Transformer::KeyValue;
  PrLHCbID2MCParticle( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator,
                     {
                         KeyValue( "MCParticlesLocation", LHCb::MCParticleLocation::Default ),
                         KeyValue( getLocation<ContainerType>(), "" )...,
                         KeyValue( getLinkLocation<ContainerType>(), "" )...,
                     },
                     KeyValue( "TargetName", "Link/Pr/LHCbID" ) ) {}

  LHCb::LinksByKey operator()( const LHCb::MCParticles& mcParts, const ContainerType&... clusters,
                               const LinksByKeyT<ContainerType>&... links ) const override {
    LHCb::LinksByKey lhcbLink{std::in_place_type<ContainedObject>, std::in_place_type<LHCb::MCParticle>,
                              LHCb::LinksByKey::Order::decreasingWeight};
    auto             buffered_counter = m_nullMCParticles.buffer();
    LoopOverClusters<ContainerType...>{}( clusters..., links..., lhcbLink, mcParts, buffered_counter );
    return lhcbLink;
  }

  mutable CounterType m_nullMCParticles{this, "#removed null MCParticles"};
};

// PrLHCbID2MCParticle for a single subdetector
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticle<std::vector<LHCb::VPFullCluster>>, "PrLHCbID2MCParticleVP" )
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticle<std::vector<LHCb::VPLightCluster>>, "PrLHCbID2MCParticleVPL" )
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticle<UT::HitHandler>, "PrLHCbID2MCParticleUT" )
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticle<LHCb::FTLiteCluster::FTLiteClusters>, "PrLHCbID2MCParticleFT" )

// PrLHCbID2MCParticle for the HLT1 chain
typedef PrLHCbID2MCParticle<std::vector<LHCb::VPLightCluster>, UT::HitHandler> PrLHCbID2MCParticleVPLUT;
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticleVPLUT, "PrLHCbID2MCParticleVPLUT" )
typedef PrLHCbID2MCParticle<std::vector<LHCb::VPFullCluster>, UT::HitHandler> PrLHCbID2MCParticleVPUT;
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticleVPUT, "PrLHCbID2MCParticleVPUT" )
typedef PrLHCbID2MCParticle<std::vector<LHCb::VPLightCluster>, UT::HitHandler, LHCb::FTLiteCluster::FTLiteClusters>
    PrLHCbID2MCParticleVPLUTFT;
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticleVPLUTFT, "PrLHCbID2MCParticleVPLUTFT" )
typedef PrLHCbID2MCParticle<std::vector<LHCb::VPFullCluster>, UT::HitHandler, LHCb::FTLiteCluster::FTLiteClusters>
    PrLHCbID2MCParticleVPUTFT;
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticleVPUTFT, "PrLHCbID2MCParticle" )
typedef PrLHCbID2MCParticle<std::vector<LHCb::VPFullCluster>, UT::HitHandler, LHCb::FTLiteCluster::FTLiteClusters,
                            MuonHitContainer>
    PrLHCbID2MCParticleVPUTFTMU;
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticleVPUTFTMU, "PrLHCbID2MCParticleVPUTFTMU" )
typedef PrLHCbID2MCParticle<std::vector<LHCb::VPFullCluster>, LHCb::FTLiteCluster::FTLiteClusters>
    PrLHCbID2MCParticleVPFT;
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticleVPFT, "PrLHCbID2MCParticleVPFT" )
typedef PrLHCbID2MCParticle<std::vector<LHCb::VPFullCluster>, LHCb::FTLiteCluster::FTLiteClusters, MuonHitContainer>
    PrLHCbID2MCParticleVPFTMU;
DECLARE_COMPONENT_WITH_ID( PrLHCbID2MCParticleVPFTMU, "PrLHCbID2MCParticleVPFTMU" )
