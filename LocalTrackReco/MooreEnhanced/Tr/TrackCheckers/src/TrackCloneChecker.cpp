/*****************************************************************************\
* (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "LHCbAlgs/Consumer.h"
#include "TrackCheckerBase.h"

#include "boost/format.hpp"

#include <map>

//-----------------------------------------------------------------------------
// Implementation file for class : TrackCloneChecker
//
// 2007-09-13 : Chris Jones
//-----------------------------------------------------------------------------

/** @class TrackCloneChecker TrackCloneChecker.h
 *
 *  Produce some simple plots for the Clone linker information
 *
 *  @author Chris Jones
 *  @date   2007-09-13
 */

class TrackCloneChecker : public LHCb::Algorithm::Consumer<void( const LHCb::Tracks&, const LHCb::MCParticles&,
                                                                 const LHCb::LinksByKey&, const LHCb::LinksByKey& ),
                                                           Gaudi::Functional::Traits::BaseClass_t<TrackCheckerBase>> {

public:
  /// Standard constructor
  TrackCloneChecker( const std::string& name, ISvcLocator* pSvcLocator );

  /** Algorithm execute */
  void operator()( const LHCb::Tracks&, const LHCb::MCParticles&, const LHCb::LinksByKey&,
                   const LHCb::LinksByKey& ) const override;

  StatusCode finalize() override; ///< Algorithm finalize

private:
  /** @class TrackTally TrackCloneChecker.h
   *
   *  Counts track information for clones
   *
   *  @author Chris Jones
   *  @date   2007-09-13
   */
  struct TrackTally {
    unsigned long totalClones{0};
    unsigned long totalNonClones{0};
    unsigned long totalGhosts{0};
    unsigned long rejectedClones{0};
    unsigned long rejectedNonClones{0};
    unsigned long rejectedGhosts{0};
    /// Map for one tally object per track history type
    typedef std::map<LHCb::Track::History, TrackTally> Map;
  };

private:
  /// Get efficiency
  inline std::pair<double, double> getEff1( const double top, const double bot ) const {
    return {( bot > 0 ? 100.0 * top / bot : 0 ),
            ( bot > 0 ? 100.0 * sqrt( ( top / bot ) * ( 1. - top / bot ) / bot ) : 0 )};
  }

  /// Get efficiency
  inline std::pair<double, double> getEff2( const double top, const double bot ) const {
    return {( bot > 0 ? top / bot : 0 ), ( bot > 0 ? sqrt( top ) / bot : 0 )};
  }

private:
  /// Summary map XXXX This is not thread safe. Should be rewritten using Gaudi counters
  mutable TrackTally::Map m_trackMap;

  /// KL distance cut
  Gaudi::Property<double> m_klCut{this, "CloneCut", 5000};

  /// Event count
  mutable Gaudi::Accumulators::Counter<> m_nEvtsCounter{this, "Nb events"};
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( TrackCloneChecker )

TrackCloneChecker::TrackCloneChecker( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator,
                {KeyValue{"TracksInContainer", LHCb::TrackLocation::Default},
                 KeyValue{"MCParticleInContainer", LHCb::MCParticleLocation::Default},
                 KeyValue{"LinkerInTable", "Link/" + LHCb::TrackLocation::Default},
                 KeyValue{"CloneLinkerLocation", "Link/" + LHCb::TrackLocation::Default + "Clones"}} ) {}

void TrackCloneChecker::operator()( const LHCb::Tracks& tracks, const LHCb::MCParticles& mcParts,
                                    const LHCb::LinksByKey& links, const LHCb::LinksByKey& cloneLinks ) const {
  ++m_nEvtsCounter;

  // Builds a lookup table to know which MC particle has clones
  // It holds how many reconstructed tracks we have for each MCParticle, identified by key
  std::map<int, unsigned int> nbReconstructedTracks;
  links.applyToAllLinks(
      [&nbReconstructedTracks]( unsigned int, unsigned int mcPartKey, float ) { ++nbReconstructedTracks[mcPartKey]; } );

  // loop over tracks
  for ( const auto* track : tracks ) {
    // MCP for main track
    const LHCb::MCParticle* mcP = mcTruth( *track, mcParts, links );
    if ( !selected( mcP ) ) continue;

    // pick up the clone info for this track
    cloneLinks.applyToLinks(
        track->key(), [&tracks, &mcParts, &links, this, &mcP]( unsigned int, unsigned int cloneKey, float weight ) {
          const LHCb::Track* cloneTrack = static_cast<const LHCb::Track*>( tracks.containedObject( cloneKey ) );
          if ( mcP ) {
            // MCP for clone track
            const LHCb::MCParticle* mcP_clone = this->mcTruth( *cloneTrack, mcParts, links );
            const bool              mcSel     = ( mcP_clone ? selected( mcP_clone ) : false );

            // log10(KLdistance)
            const double logFLdist = log10( weight );
            // const bool cloneID = ( weight < m_klCut );

            if ( mcP_clone && mcSel ) {
              if ( mcP == mcP_clone ) {
                plot1D( logFLdist, "KLDtrueClones", "Log10(KLDistance) | True Clones", -5, 13, 100 );
              } else {
                plot1D( logFLdist, "KLDnotClones", "Log10(KLDistance) | Not Clones", -5, 13, 100 );
              }
            } else if ( mcP_clone && !mcSel ) {
              plot1D( logFLdist, "KLDrejMCPs", "Log10(KLDistance) | Rejected MCParticles", -5, 13, 100 );
            } else {
              plot1D( logFLdist, "KLDghosts", "Log10(KLDistance) | Ghosts", -5, 13, 100 );
            }
          }
        } );

    // clone ID
    const bool cloneID = ( track->info( LHCb::Track::AdditionalInfo::CloneDist, 9e99 ) < m_klCut );

    // real clone ?
    const bool hasMCclones = ( nbReconstructedTracks.at( mcP->key() ) > 1 );

    // tally object
    TrackTally& tally = m_trackMap[track->history()];

    // tallies
    if ( mcP ) {
      if ( hasMCclones ) ++tally.totalClones;
      if ( !hasMCclones ) ++tally.totalNonClones;
      if ( hasMCclones && cloneID ) ++tally.rejectedClones;
      if ( !hasMCclones && cloneID ) ++tally.rejectedNonClones;
    } else {
      ++tally.totalGhosts;
      if ( cloneID ) ++tally.rejectedGhosts;
    }

  } // track loop
}

StatusCode TrackCloneChecker::finalize() {
  const std::string& lines =
      "============================================================================================";
  always() << lines << endmsg;
  always() << "      Clone summary for '" << inputLocation<3>() << "' IDed clones with KLdist<" << m_klCut << endmsg;
  always() << lines << endmsg;

  std::pair<double, double> r1, r2, r3, r4, r5;

  always() << "   Track type     NonClones/Evt  Clones/Evt     ClonesID/%     NonClonesID/%  GhostsID/%" << endmsg;
  for ( auto iM = m_trackMap.begin(); iM != m_trackMap.end(); ++iM ) {
    const TrackTally& tally = iM->second;
    r1                      = getEff1( tally.rejectedClones, tally.totalClones / 2.0 );
    r2                      = getEff1( tally.rejectedNonClones, tally.totalNonClones );
    r3                      = getEff1( tally.rejectedGhosts, tally.totalGhosts );
    r4                      = getEff2( tally.totalNonClones, m_nEvtsCounter.nEntries() );
    r5                      = getEff2( tally.totalClones / 2.0, m_nEvtsCounter.nEntries() );
    always() << boost::format( "%15s %6.2f +-%5.2f %6.2f +-%5.2f %6.2f +-%5.2f %6.2f +-%5.2f %6.2f +-%5.2f" ) %
                    iM->first % r4.first % r4.second % r5.first % r5.second % r1.first % r1.second % r2.first %
                    r2.second % r3.first % r3.second
             << endmsg;
  }

  always() << lines << endmsg;

  return Consumer::finalize();
}
