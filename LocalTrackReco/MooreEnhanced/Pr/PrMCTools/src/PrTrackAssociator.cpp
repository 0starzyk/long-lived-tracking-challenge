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
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/MCVertex.h"
#include "Event/Track.h"
#include "GaudiAlg/FunctionalDetails.h"
#include "Kernel/LHCbID.h"
#include "Kernel/STLExtensions.h"
#include "LHCbAlgs/Transformer.h"

namespace {
  using KeyValue = std::pair<std::string, std::string>;

  struct TruthCounter {
    TruthCounter() = default;
    TruthCounter( const LHCb::MCParticle* part ) : particle( part ) {}
    const LHCb::MCParticle* particle{nullptr};
    unsigned int            nVelo{0};
    unsigned int            nTT{0};
    unsigned int            nT{0};
    unsigned int            nMuon{0};
    unsigned int            nMeasurements() const { return nVelo + nTT + nT; };
  };

  TruthCounter& getCounter( const LHCb::MCParticle* part, std::vector<TruthCounter>& truthCounters ) {
    auto it = std::find_if( begin( truthCounters ), end( truthCounters ),
                            [part]( auto& item ) { return item.particle == part; } );
    return it == truthCounters.end() ? truthCounters.emplace_back( part ) : *it;
  }
  void incrementVelo( const LHCb::MCParticle* part, std::vector<TruthCounter>& truthCounters ) {
    ( getCounter( part, truthCounters ).nVelo )++;
  }
  void incrementTT( const LHCb::MCParticle* part, std::vector<TruthCounter>& truthCounters ) {
    ( getCounter( part, truthCounters ).nTT )++;
  }
  void incrementT( const LHCb::MCParticle* part, std::vector<TruthCounter>& truthCounters ) {
    ( getCounter( part, truthCounters ).nT )++;
  }
  void incrementMuon( const LHCb::MCParticle* part, std::vector<TruthCounter>& truthCounters ) {
    ( getCounter( part, truthCounters ).nMuon )++;
  }

  std::tuple<TruthCounter, std::vector<TruthCounter>> match_track( LHCb::span<LHCb::LHCbID const> lhcbIDs,
                                                                   const LHCb::MCParticles&       mcParts,
                                                                   const LHCb::LinksByKey&        idlinks ) {
    /// total number of measurements
    TruthCounter total;
    /// vector of counters of associated MC particles
    std::vector<TruthCounter> truthCounters;

    // Loop over collection of LHCbIDs
    for ( const auto& ids : lhcbIDs ) {
      if ( ids.isVP() ) {
        ++total.nVelo; // Count number of Velo hits
        idlinks.applyToLinks( ids.lhcbID(), [&mcParts, &truthCounters]( unsigned int, unsigned int mcPartKey, float ) {
          auto* part = static_cast<const LHCb::MCParticle*>( mcParts.containedObject( mcPartKey ) );
          if ( part && &mcParts == part->parent() ) incrementVelo( part, truthCounters );
        } );
      } else if ( ids.isUT() ) {
        ++total.nTT; // Count number of TT hits
        idlinks.applyToLinks( ids.lhcbID(), [&mcParts, &truthCounters]( unsigned int, unsigned int mcPartKey, float ) {
          auto* part = static_cast<const LHCb::MCParticle*>( mcParts.containedObject( mcPartKey ) );
          if ( part && &mcParts == part->parent() ) incrementTT( part, truthCounters );
        } );
      } else if ( ids.isFT() ) {
        ++total.nT; // Count number of T hits
        idlinks.applyToLinks( ids.lhcbID(), [&mcParts, &truthCounters]( unsigned int, unsigned int mcPartKey, float ) {
          auto* part = static_cast<const LHCb::MCParticle*>( mcParts.containedObject( mcPartKey ) );
          if ( part && &mcParts == part->parent() ) incrementT( part, truthCounters );
        } );
      } else if ( ids.isMuon() ) {
        ++total.nMuon; // Count number of Muon hits
        idlinks.applyToLinks( ids.lhcbID(), [&mcParts, &truthCounters]( unsigned int, unsigned int mcPartKey, float ) {
          auto* part = static_cast<const LHCb::MCParticle*>( mcParts.containedObject( mcPartKey ) );
          if ( part && &mcParts == part->parent() ) incrementMuon( part, truthCounters );
        } );
      }
    }

    // If the Track has total # Velo hits > 2 AND total # T hits > 2, cumul mother and daughter
    if ( ( 2 < total.nVelo ) && ( 2 < total.nT ) ) {
      for ( auto& counter1 : truthCounters ) {
        if ( counter1.nT == 0 ) continue;
        const LHCb::MCVertex* vOrigin = counter1.particle->originVertex();
        if ( 0 != vOrigin ) {
          const LHCb::MCParticle* mother = vOrigin->mother();
          if ( 0 == mother ) continue; // no ancestor;
          for ( auto& counter2 : truthCounters ) {
            if ( mother == counter2.particle ) {
              if ( counter2.nVelo == 0 ) continue;

              // if ( msgLevel( MSG::DEBUG ) )
              //   debug() << "  *** Particle " << counter1.particle->key() << "[" <<
              //   counter1.particle->particleID().pid()
              //           << "] (" << counter1.nVelo << "," << counter1.nTT << "," << counter1.nT << ")"
              //           << " is daughter of " << counter2.particle->key() << "["
              //           << counter2.particle->particleID().pid() << "] (" << counter2.nVelo << "," << counter2.nTT
              //           << "," << counter2.nT << ")"
              //           << " type " << vOrigin->type() << ". Merge hits to tag both." << endmsg;

              //== Daughter hits are added to mother.
              counter2.nVelo += counter1.nVelo;
              counter2.nTT += counter1.nTT;
              counter2.nT += counter1.nT;
              if ( counter2.nVelo > total.nVelo ) counter2.nVelo = total.nVelo;
              if ( counter2.nTT > total.nTT ) counter2.nTT = total.nTT;
              if ( counter2.nT > total.nT ) counter2.nT = total.nT;

              //== Mother hits overwrite Daughter hits
              counter1.nVelo = counter2.nVelo;
              counter1.nTT   = counter2.nTT;
              counter1.nT    = counter2.nT;
            }
          }
        }
      }
    }
    return {total, truthCounters};
  }

  std::optional<double> matchingFraction( double fractionOK, TruthCounter const& counter, TruthCounter const& total ) {
    //===============================================================
    // Association definition
    // Velo matching:
    //      * either less than 2 hits in total (no Velo on track)
    //      * or at least 'fractionOK' of the Velo hits are associated to this MCParticle.
    // TT matching:
    //      * this MC particle has at least nTT-2 hits.
    //      * or this track has both Velo and T stations: TT not used for long tracks.
    // T matching:
    //      * either less than 2 hits in total (no T station hits on track)
    //      * or at least 'fractionOK' of the T station hits are associated to this MCParticle.
    //
    // Muon matching needed in tracking efficiency study
    //     * this MC particle has at least nMuon 3 hits
    //     * Muon hits only counted for VeloMuon / MuonUT / standalone Muon tracks
    // A MCParticle matches if all 4 criteria are OK.
    //
    // The weight is given by the number of hits associated to this MCParticle divided by the total number of hits.
    //
    //===============================================================

    bool veloOK = true;
    if ( 2 < total.nVelo ) {
      veloOK       = false;
      double ratio = (double)counter.nVelo / total.nVelo;
      if ( fractionOK <= ratio ) { veloOK = true; }
    }

    bool tOK = true;
    if ( 2 < total.nT ) {
      tOK          = false;
      double ratio = (double)counter.nT / total.nT;
      if ( fractionOK <= ratio ) { tOK = true; }
    }

    bool ttOK = ( counter.nTT + 2 > total.nTT );
    if ( 2 < total.nVelo && 2 < total.nT ) { ttOK = true; }

    bool MuonOK = ( total.nMuon > 2 && total.nMuon < counter.nMuon + 2 );

    // for standalone Muon/ VeloMuon /MuonUT match for tracking efficiency study
    if ( total.nMeasurements() == 0 && MuonOK ) { return (double)( counter.nMuon ) / total.nMuon; }
    if ( total.nTT == 0 && total.nT == 0 && veloOK && MuonOK ) {
      return (double)( counter.nMuon + counter.nVelo ) / ( total.nVelo + total.nMuon );
    }
    if ( total.nVelo == 0 && total.nT == 0 && ttOK && MuonOK ) {
      return (double)( counter.nMuon + counter.nTT ) / ( total.nTT + total.nMuon );
    }

    if ( veloOK && tOK && ttOK && ( 0 < total.nMeasurements() ) ) {
      return (double)( counter.nVelo + counter.nTT + counter.nT ) / total.nMeasurements();
    }
    return {};
  }

} // namespace

/** @class PrTrackAssociator PrTrackAssociator.cpp
 *
 *  This algorithm computes the link between a Track and a MCParticle.
 *  The requirement is a match of both the Velo/VP and the T part of the
 *  Track. If there are not enough coordinates, the match is assumed so that
 *  a Velo only or a T only are matched properly.
 *  The required fraction of hits is a jobOption 'FractionOK', default 0.70.
 *
 *  Rewritten for the upgrade, handles all containers in one instance
 *
 *  The class is templated on the track type. The current requirement is that the underlying container provides
 *  a way to iterate over single objects, a unique integer per object and the single objects provides
 *  a function which returns a view into the list of lhcb IDs of a track. For KeyedContainers the key of an object
 *  is used as identifier. For plain vectors the index is used. It is mandatory that the same identifier is used
 *  when reading the LinksByKey object created by this algorithm.
 *
 *  Note: The dependency on MCVertices is quite implicit. They have to be unpacked before this algorithms runs.
 *  Otherwise mcparticle->originVertex() always returns a null pointer, even if it should not.
 *
 *  Original author Olivier Callot, see revision history for other contributors:
 *  @author Olivier Callot
 *  @date   2012-04-04
 */
template <typename Tracks>
class PrTrackAssociator
    : public LHCb::Algorithm::Transformer<LHCb::LinksByKey( const LHCb::MCParticles&, const LHCb::MCVertices&,
                                                            const Tracks&, const LHCb::LinksByKey& )> {
public:
  PrTrackAssociator( const std::string& name, ISvcLocator* pSvcLocator )
      : LHCb::Algorithm::Transformer<LHCb::LinksByKey( const LHCb::MCParticles&, const LHCb::MCVertices&, const Tracks&,
                                                       const LHCb::LinksByKey& )>(
            name, pSvcLocator,
            {KeyValue{"MCParticleLocation", LHCb::MCParticleLocation::Default},
             KeyValue{"MCVerticesInput", LHCb::MCVertexLocation::Default}, KeyValue{"SingleContainer", ""},
             KeyValue{"LinkerLocationID", "Link/Pr/LHCbID"}},
            KeyValue{"OutputLocation", ""} ){};
  LHCb::LinksByKey operator()( const LHCb::MCParticles& mcParts, const LHCb::MCVertices& /* */, const Tracks& tracks,
                               const LHCb::LinksByKey& idlinks ) const override {
    // Create the Linker table from Track to MCParticle
    // Sorted by decreasing weight, so first retrieved has highest weight
    // This has to be done, even if there are no tracks in the event, to satisfy the DST writer
    LHCb::LinksByKey result{
        std::in_place_type<std::conditional_t<std::is_same_v<Tracks, LHCb::Track::Range>, LHCb::Track, void>>,
        std::in_place_type<LHCb::MCParticle>, LHCb::LinksByKey::Order::decreasingWeight};

    // Loop over the Tracks
    for ( auto const& [index, tr] : LHCb::range::enumerate( tracks ) ) {
      /// total number of measurements and
      /// vector of counters of associated MC particles
      auto [total, truthCounters] = match_track( Gaudi::Functional::details::deref( tr ).lhcbIDs(), mcParts, idlinks );

      bool         is_associated = false;
      unsigned int n_mcparticles = 0;
      for ( auto const& counter : truthCounters ) {
        auto matched = matchingFraction( m_fractionOK, counter, total );
        //=== Decision. Fill Linker
        if ( !matched ) continue;
        is_associated = true;
        n_mcparticles++;
        if constexpr ( std::is_same_v<Tracks, LHCb::Track::Range> ) {
          result.link( tr, counter.particle, *matched );
        } else {
          result.link( index, counter.particle, *matched );
        }
      }
      m_efficiency += is_associated;
      if ( is_associated ) { m_mcparticles_per_track += n_mcparticles; }
    } // End loop over Tracks
    return result;
  };

private:
  Gaudi::Property<double> m_fractionOK{this, "FractionOK", 0.70, "minimal good matching fraction"};
  mutable Gaudi::Accumulators::BinomialCounter<>  m_efficiency{this, "Efficiency"};
  mutable Gaudi::Accumulators::AveragingCounter<> m_mcparticles_per_track{this, "MC particles per track"};
};

DECLARE_COMPONENT_WITH_ID( PrTrackAssociator<LHCb::Track::Range>, "PrTrackAssociator" )
DECLARE_COMPONENT_WITH_ID( PrTrackAssociator<std::vector<LHCb::Event::v2::Track>>, "PrV2TrackAssociator" )
