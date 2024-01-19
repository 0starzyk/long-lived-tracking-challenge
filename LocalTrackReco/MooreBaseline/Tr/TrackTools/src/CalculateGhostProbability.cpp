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
#include "Event/GhostProbability.h"
#include "Event/PartialChiSquareds.h"
#include "Event/PrHits.h"
#include "Event/SOACollectionMerger.h"
#include "Event/Track.h"
#include "Event/TrackFitResult.h"
#include "Event/Track_v3.h"
#include "Event/UTCluster.h"
#include "Event/VPLightCluster.h"
#include "GaudiKernel/DataHandle.h"
#include "LHCbAlgs/Transformer.h"
#include "PrKernel/UTHitHandler.h"
#include "TMVA/TMVA_Downstream_25nsLL_2017_MLP.h"
#include "TMVA/TMVA_Long_25nsLL_2017_MLP.h"
#include "TMVA/TMVA_T_25nsLL_2017_MLP.h"
#include "TMVA/TMVA_Upstream_25nsLL_2017_MLP.h"
#include "TMVA/TMVA_Velo_25nsLL_2017_MLP.h"
#include "TMVA/UpdateFlattenDownstream.h"
#include "TMVA/UpdateFlattenLong.h"
#include "TMVA/UpdateFlattenTtrack.h"
#include "TMVA/UpdateFlattenUpstream.h"
#include "TMVA/UpdateFlattenVelo.h"
#include "TMath.h"
#include "TrackKernel/TrackFunctors.h"
//-----------------------------------------------------------------------------
// Implementation file for class : UpgradeGhostId
//
// 2014-12-30 : Paul Seyfert following an earlier version by Angelo Di Canto
// 2019-06-18 : Menglin Xu
//-----------------------------------------------------------------------------

/** @class UpgradeGhostId UpgradeGhostId.h
 *
 *  @author Paul Seyfert
 *  @date   30-12-2014
 */

namespace {

  static const std::array<std::string_view, 8> veloVars = {"UpgradeGhostInfo_obsVP",
                                                           "UpgradeGhostInfo_FitVeloChi2",
                                                           "UpgradeGhostInfo_FitVeloNDoF",
                                                           "UpgradeGhostInfo_veloHits",
                                                           "UpgradeGhostInfo_utHits",
                                                           "TRACK_CHI2",
                                                           "TRACK_NDOF",
                                                           "TRACK_ETA"};

  static const std::array<std::string_view, 11> upstreamVars = {"UpgradeGhostInfo_obsVP",
                                                                "UpgradeGhostInfo_FitVeloChi2",
                                                                "UpgradeGhostInfo_FitVeloNDoF",
                                                                "UpgradeGhostInfo_obsUT",
                                                                "UpgradeGhostInfo_UToutlier",
                                                                "UpgradeGhostInfo_veloHits",
                                                                "UpgradeGhostInfo_utHits",
                                                                "TRACK_CHI2",
                                                                "TRACK_NDOF",
                                                                "TRACK_PT",
                                                                "TRACK_ETA"};

  static const std::array<std::string_view, 10> downstreamVars = {"UpgradeGhostInfo_obsFT",
                                                                  "UpgradeGhostInfo_FitTChi2",
                                                                  "UpgradeGhostInfo_FitTNDoF",
                                                                  "UpgradeGhostInfo_obsUT",
                                                                  "UpgradeGhostInfo_UToutlier",
                                                                  "UpgradeGhostInfo_veloHits",
                                                                  "UpgradeGhostInfo_utHits",
                                                                  "TRACK_CHI2",
                                                                  "TRACK_PT",
                                                                  "TRACK_ETA"};

  static const std::array<std::string_view, 14> longVars = {"UpgradeGhostInfo_obsVP",
                                                            "UpgradeGhostInfo_FitVeloChi2",
                                                            "UpgradeGhostInfo_FitVeloNDoF",
                                                            "UpgradeGhostInfo_obsFT",
                                                            "UpgradeGhostInfo_FitTChi2",
                                                            "UpgradeGhostInfo_FitTNDoF",
                                                            "UpgradeGhostInfo_obsUT",
                                                            "UpgradeGhostInfo_FitMatchChi2",
                                                            "UpgradeGhostInfo_UToutlier",
                                                            "UpgradeGhostInfo_veloHits",
                                                            "UpgradeGhostInfo_utHits",
                                                            "TRACK_CHI2",
                                                            "TRACK_PT",
                                                            "TRACK_ETA"};

  static const std::array<std::string_view, 9> ttrackVars = {"UpgradeGhostInfo_obsFT",
                                                             "UpgradeGhostInfo_FitTChi2",
                                                             "UpgradeGhostInfo_FitTNDoF",
                                                             "UpgradeGhostInfo_veloHits",
                                                             "UpgradeGhostInfo_utHits",
                                                             "TRACK_CHI2",
                                                             "TRACK_NDOF",
                                                             "TRACK_PT",
                                                             "TRACK_ETA"};
} // namespace

namespace LHCb::Event::v3 {

  class CalculateGhostProbability
      : public LHCb::Algorithm::Transformer<LHCb::Event::v3::GhostProbabilities(
            const LHCb::Event::v3::Tracks&, const LHCb::Event::v3::Track::PartialChiSquareds& info )> {

    using simd = SIMDWrapper::best::types;

  public:
    /// Standard constructor
    CalculateGhostProbability( const std::string& name, ISvcLocator* pSvcLocator );

    StatusCode initialize() override;

    LHCb::Event::v3::GhostProbabilities
    operator()( const LHCb::Event::v3::Tracks&,
                const LHCb::Event::v3::Track::PartialChiSquareds& info ) const override final;

  private:
    std::vector<std::unique_ptr<::Track::TabulatedFunction1D>> m_flatters;

    DataObjectReadHandle<LHCb::VPLightClusters> m_vpClusters{this, "VPClusterLocation", LHCb::VPClusterLocation::Light,
                                                             "Location of VP light clusters"};
    DataObjectReadHandle<LHCb::Pr::Hits<LHCb::Pr::HitType::UT>> m_utClusters{this, "UTClusterLocation", "",
                                                                             "Location of UT clusters"};

    ReadGhostProbabilityLong       m_longClassifier{longVars};
    ReadGhostProbabilityDownstream m_downstreamClassifier{downstreamVars};
    ReadGhostProbabilityVelo       m_veloClassifier{veloVars};
    ReadGhostProbabilityUpstream   m_upstreamClassifier{upstreamVars};
    ReadGhostProbabilityT          m_tClassifier{ttrackVars};

    mutable Gaudi::Accumulators::MsgCounter<MSG::ERROR> m_unknownTrackType{
        this, "Unknown track type - setting ghost probability to +inf"};
  };
  DECLARE_COMPONENT( CalculateGhostProbability )
} // namespace LHCb::Event::v3

namespace LHCb::Event::v3 {
  CalculateGhostProbability::CalculateGhostProbability( const std::string& name, ISvcLocator* pSvcLocator )
      : Transformer( name, pSvcLocator, {KeyValue{"InputTracksName", "Rec/Track/Velo"}, KeyValue{"InputInfoName", ""}},
                     KeyValue{"OutputGhostProbName", "Rec/Track/GhostProb"} ) {}

  /// Initialization
  StatusCode CalculateGhostProbability::initialize() {

    return Transformer::initialize().andThen( [&] {
      m_flatters.clear();
      m_flatters.resize( static_cast<int>( LHCb::Track::Types::Last ) );
      m_flatters[static_cast<int>( LHCb::Track::Types::Velo )]       = Update_VeloTable();
      m_flatters[static_cast<int>( LHCb::Track::Types::Upstream )]   = Update_UpstreamTable();
      m_flatters[static_cast<int>( LHCb::Track::Types::Downstream )] = Update_DownstreamTable();
      m_flatters[static_cast<int>( LHCb::Track::Types::Long )]       = Update_LongTable();
      m_flatters[static_cast<int>( LHCb::Track::Types::Ttrack )]     = Update_TtrackTable();
    } );
  }

  GhostProbabilities CalculateGhostProbability::
                     operator()( const LHCb::Event::v3::Tracks&                    inputTracks,
              const LHCb::Event::v3::Track::PartialChiSquareds& info ) const {

    using SL = LHCb::Event::v3::Tracks::StateLocation;

    GhostProbabilities ghostProbs{inputTracks.zipIdentifier()};
    if ( inputTracks.empty() ) { return ghostProbs; }

    // -- we'll make it scalar for now
    auto trackInfoZip = LHCb::Event::make_zip<SIMDWrapper::Scalar>( inputTracks, info );

    for ( auto const& track : trackInfoZip ) {

      float netresponse = std::numeric_limits<float>::max();

      switch ( track.type() ) {
      case LHCb::Track::Types::Long: {
        std::array<float, longVars.size()> variables{};
        variables[0]  = track.nVPHits().cast();
        variables[1]  = track.FitVeloChi2().cast();
        variables[2]  = track.FitVeloNDoF().cast();
        variables[3]  = track.nFTHits().cast();
        variables[4]  = track.FitTChi2().cast();
        variables[5]  = track.FitTNDoF().cast();
        variables[6]  = track.nUTHits().cast();
        variables[7]  = track.FitMatchChi2().cast();
        variables[8]  = track.NUTOutliers().cast();
        variables[9]  = m_vpClusters.get()->size();
        variables[10] = m_utClusters.get()->nHits();
        variables[11] = track.chi2().cast();
        variables[12] = track.pt( SL::ClosestToBeam ).cast();
        variables[13] = track.pseudoRapidity( SL::ClosestToBeam ).cast();
        netresponse   = m_longClassifier.GetMvaValue( variables );
      } break;
      case LHCb::Track::Types::Downstream: {
        std::array<float, downstreamVars.size()> variables{};
        variables[0] = track.nFTHits().cast();
        variables[1] = track.FitTChi2().cast();
        variables[2] = track.FitTNDoF().cast();
        variables[3] = track.nUTHits().cast();
        variables[4] = track.NUTOutliers().cast();
        variables[5] = m_vpClusters.get()->size();
        variables[6] = m_utClusters.get()->nHits();
        variables[7] = track.chi2().cast();
        variables[8] = track.pt( SL::ClosestToBeam ).cast();
        variables[9] = track.pseudoRapidity( SL::ClosestToBeam ).cast();
        netresponse  = m_downstreamClassifier.GetMvaValue( variables );
      } break;
      case LHCb::Track::Types::Velo: {
        std::array<float, veloVars.size()> variables{};
        variables[0] = track.nVPHits().cast();
        variables[1] = track.FitVeloChi2().cast();
        variables[2] = track.FitVeloNDoF().cast();
        variables[3] = m_vpClusters.get()->size();
        variables[4] = m_utClusters.get()->nHits();
        variables[5] = track.chi2().cast();
        variables[6] = track.nDoF().cast();
        variables[7] = track.pseudoRapidity( SL::ClosestToBeam ).cast();
        netresponse  = m_veloClassifier.GetMvaValue( variables );
      } break;
      case LHCb::Track::Types::Upstream: {
        std::array<float, upstreamVars.size()> variables{};
        variables[0]  = track.nVPHits().cast();
        variables[1]  = track.FitVeloChi2().cast();
        variables[2]  = track.FitVeloNDoF().cast();
        variables[3]  = track.nUTHits().cast();
        variables[4]  = track.NUTOutliers().cast();
        variables[5]  = m_vpClusters.get()->size();
        variables[6]  = m_utClusters.get()->nHits();
        variables[7]  = track.chi2().cast();
        variables[8]  = track.nDoF().cast();
        variables[9]  = track.pt( SL::ClosestToBeam ).cast();
        variables[10] = track.pseudoRapidity( SL::ClosestToBeam ).cast();
        netresponse   = m_upstreamClassifier.GetMvaValue( variables );
      } break;
      case LHCb::Track::Types::Ttrack: {
        std::array<float, ttrackVars.size()> variables{};
        variables[0] = track.nFTHits().cast();
        variables[1] = track.FitTChi2().cast();
        variables[2] = track.FitTNDoF().cast();
        variables[3] = m_vpClusters.get()->size();
        variables[4] = m_utClusters.get()->nHits();
        variables[5] = track.chi2().cast();
        variables[6] = track.nDoF().cast();
        variables[7] = track.pt( SL::ClosestToBeam ).cast();
        variables[8] = track.pseudoRapidity( SL::ClosestToBeam ).cast();
        netresponse  = m_tClassifier.GetMvaValue( variables );
      } break;
      default:
        ++m_unknownTrackType;
        break;
      }

      netresponse   = m_flatters[static_cast<int>( track.type() )]->value( netresponse );
      auto const gp = ghostProbs.emplace_back<SIMDWrapper::InstructionSet::Scalar>();
      gp.field<GhostProbabilityTag::GhostProbability>().set( 1.0 - netresponse );
    }

    return ghostProbs;
  }
  DECLARE_COMPONENT_WITH_ID( SOACollectionMerger<LHCb::Event::v3::GhostProbabilities>, "GhostProbabilityMerger" )

} // namespace LHCb::Event::v3
