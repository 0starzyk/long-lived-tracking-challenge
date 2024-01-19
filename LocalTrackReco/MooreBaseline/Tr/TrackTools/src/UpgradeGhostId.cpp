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
#include "Event/FTCluster.h"
#include "Event/GhostTrackInfo.h"
#include "Event/Particle.h"
#include "Event/PrHits.h"
#include "Event/Track.h"
#include "Event/TrackFitResult.h"
#include "Event/UTCluster.h"
#include "Event/VPLightCluster.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/DataHandle.h"
#include "GaudiKernel/IIncidentListener.h"
#include "Kernel/HitPattern.h"
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
#include "TrackInterfaces/IGhostProbability.h"
#include "TrackKernel/TrackFunctors.h"
#include <boost/container/static_vector.hpp>

//-----------------------------------------------------------------------------
// Implementation file for class : UpgradeGhostId
//
// 2014-12-30 : Paul Seyfert following an earlier version by Angelo Di Canto
// 2019-06-18 : Menglin Xu
//-----------------------------------------------------------------------------
//=============================================================================

namespace {

  inline constexpr std::array<std::string_view, 8> veloVars = {"UpgradeGhostInfo_obsVP",
                                                               "UpgradeGhostInfo_FitVeloChi2",
                                                               "UpgradeGhostInfo_FitVeloNDoF",
                                                               "UpgradeGhostInfo_veloHits",
                                                               "UpgradeGhostInfo_utHits",
                                                               "TRACK_CHI2",
                                                               "TRACK_NDOF",
                                                               "TRACK_ETA"};

  inline constexpr std::array<std::string_view, 11> upstreamVars = {"UpgradeGhostInfo_obsVP",
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

  inline constexpr std::array<std::string_view, 10> downstreamVars = {"UpgradeGhostInfo_obsFT",
                                                                      "UpgradeGhostInfo_FitTChi2",
                                                                      "UpgradeGhostInfo_FitTNDoF",
                                                                      "UpgradeGhostInfo_obsUT",
                                                                      "UpgradeGhostInfo_UToutlier",
                                                                      "UpgradeGhostInfo_veloHits",
                                                                      "UpgradeGhostInfo_utHits",
                                                                      "TRACK_CHI2",
                                                                      "TRACK_PT",
                                                                      "TRACK_ETA"};

  inline constexpr std::array<std::string_view, 14> longVars = {"UpgradeGhostInfo_obsVP",
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

  inline constexpr std::array<std::string_view, 9> ttrackVars = {"UpgradeGhostInfo_obsFT",
                                                                 "UpgradeGhostInfo_FitTChi2",
                                                                 "UpgradeGhostInfo_FitTNDoF",
                                                                 "UpgradeGhostInfo_veloHits",
                                                                 "UpgradeGhostInfo_utHits",
                                                                 "TRACK_CHI2",
                                                                 "TRACK_NDOF",
                                                                 "TRACK_PT",
                                                                 "TRACK_ETA"};
} // namespace

/** @class UpgradeGhostId UpgradeGhostId.h
 *
 *  @author Paul Seyfert
 *  @date   30-12-2014
 */
class UpgradeGhostId : public extends<GaudiTool, IGhostProbability> {

public:
  /// Standard constructor
  using extends::extends;

  // using extends::extends;

  StatusCode finalize() override;   ///< Algorithm initialization
  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute( LHCb::Track& aTrack ) const override;
  /** reveal the variable names for a track type */
  std::vector<std::string_view> variableNames( LHCb::Track::Types type ) const override;

  /** reveal the variable values for a track */
  std::vector<float> netInputs( LHCb::Track& ) const override {
    fatal() << "UpgradeGhostId::netInputs is NOT IMPLEMENTED" << endmsg;
    return {};
  }

private:
  std::vector<std::unique_ptr<Track::TabulatedFunction1D>> m_flatters;

  std::vector<std::string> m_inNames;

  DataObjectReadHandle<LHCb::VPLightClusters> m_vpClusters{this, "VPClusterLocation", LHCb::VPClusterLocation::Light,
                                                           "Location of VP light clusters"};
  DataObjectReadHandle<LHCb::Pr::Hits<LHCb::Pr::HitType::UT>> m_utClusters{this, "UTClusterLocation",
                                                                           "/Event/UT/PrUTHits", "Location of UT hits"};

  ReadGhostProbabilityLong       m_longClassifier{longVars};
  ReadGhostProbabilityDownstream m_downstreamClassifier{downstreamVars};
  ReadGhostProbabilityVelo       m_veloClassifier{veloVars};
  ReadGhostProbabilityUpstream   m_upstreamClassifier{upstreamVars};
  ReadGhostProbabilityT          m_tClassifier{ttrackVars};
};

DECLARE_COMPONENT( UpgradeGhostId )

StatusCode UpgradeGhostId::finalize() {

  std::for_each( m_flatters.begin(), m_flatters.end(), []( auto& up ) { up.reset(); } );
  return GaudiTool::finalize();
}

StatusCode UpgradeGhostId::initialize() {
  return GaudiTool::initialize().andThen( [&] {
    m_flatters.clear();
    m_flatters.resize( static_cast<int>( LHCb::Track::Types::Last ) );
    m_flatters[static_cast<int>( LHCb::Track::Types::Velo )]         = Update_VeloTable();
    m_flatters[static_cast<int>( LHCb::Track::Types::VeloBackward )] = Update_VeloTable();
    m_flatters[static_cast<int>( LHCb::Track::Types::Upstream )]     = Update_UpstreamTable();
    m_flatters[static_cast<int>( LHCb::Track::Types::Downstream )]   = Update_DownstreamTable();
    m_flatters[static_cast<int>( LHCb::Track::Types::Long )]         = Update_LongTable();
    m_flatters[static_cast<int>( LHCb::Track::Types::Ttrack )]       = Update_TtrackTable();
  } );
}

namespace {
  class SubDetHits {
    std::array<int, 3>   m_hits = {};
    constexpr static int idx( LHCb::LHCbID::channelIDtype t ) {
      switch ( t ) {
      case LHCb::LHCbID::channelIDtype::VP:
        return 0;
      case LHCb::LHCbID::channelIDtype::UT:
        return 1;
      case LHCb::LHCbID::channelIDtype::FT:
        return 2;
      default:
        throw;
      }
    }
    static constexpr bool validType( LHCb::LHCbID::channelIDtype t ) {
      return t == LHCb::LHCbID::channelIDtype::VP || t == LHCb::LHCbID::channelIDtype::UT ||
             t == LHCb::LHCbID::channelIDtype::FT;
    }

  public:
    SubDetHits( const LHCb::Track& aTrack ) {
      for ( auto lhcbid : aTrack.lhcbIDs() ) {
        if ( !validType( lhcbid.detectorType() ) ) continue; // may be a hit in a non-tracking detector
        ++m_hits[idx( lhcbid.detectorType() )];
      }
    }

    double operator[]( LHCb::LHCbID::channelIDtype t ) const { return m_hits[idx( t )]; }
  };
} // namespace

//=============================================================================
StatusCode UpgradeGhostId::execute( LHCb::Track& aTrack ) const {
  auto obsarray = SubDetHits{aTrack};

  boost::container::static_vector<float, 15> variables;
  if ( aTrack.hasVelo() ) {
    variables.push_back( obsarray[LHCb::LHCbID::channelIDtype::VP] );
    variables.push_back( aTrack.info( LHCb::Track::AdditionalInfo::FitVeloChi2, -999 ) );
    variables.push_back( aTrack.info( LHCb::Track::AdditionalInfo::FitVeloNDoF, -999 ) );
  }
  if ( aTrack.hasT() ) {
    variables.push_back( obsarray[LHCb::LHCbID::channelIDtype::FT] );
    variables.push_back( aTrack.info( LHCb::Track::AdditionalInfo::FitTChi2, -999 ) );
    variables.push_back( aTrack.info( LHCb::Track::AdditionalInfo::FitTNDoF, -999 ) );
  }
  if ( aTrack.hasUT() ) { // includes longtracks w/o ut hits
    variables.push_back( obsarray[LHCb::LHCbID::channelIDtype::UT] );
  }
  if ( LHCb::Track::Types::Long == aTrack.type() ) {
    variables.push_back( aTrack.info( LHCb::Track::AdditionalInfo::FitMatchChi2, -999 ) );
  }
  if ( aTrack.hasUT() ) {
    const LHCb::TrackFitResult* fit = fitResult( aTrack );
    variables.push_back(
        fit ? ( fit->nMeasurements<LHCb::Measurement::UT>() - fit->nActiveMeasurements<LHCb::Measurement::UT>() )
            : 0 ); // "UpgradeGhostInfo_UToutlier",'F'
  }
  variables.push_back( m_vpClusters.get()->size() );
  variables.push_back( m_utClusters.get()->nHits() );
  variables.push_back( aTrack.chi2() );
  if ( LHCb::Track::Types::Long != aTrack.type() && LHCb::Track::Types::Downstream != aTrack.type() ) {
    variables.push_back( aTrack.nDoF() );
  }
  if ( ( LHCb::Track::Types::Velo != aTrack.type() ) && !aTrack.isVeloBackward() ) {
    variables.push_back( aTrack.pt() );
  }
  variables.push_back( aTrack.pseudoRapidity() );

  float netresponse = std::numeric_limits<float>::max();

  switch ( aTrack.type() ) {
  case LHCb::Track::Types::Long:
    netresponse = m_longClassifier.GetMvaValue( LHCb::span{variables}.first<longVars.size()>() );
    break;
  case LHCb::Track::Types::Downstream:
    netresponse = m_downstreamClassifier.GetMvaValue( LHCb::span{variables}.first<downstreamVars.size()>() );
    break;
  case LHCb::Track::Types::Velo:
  case LHCb::Track::Types::VeloBackward:
    netresponse = m_veloClassifier.GetMvaValue( LHCb::span{variables}.first<veloVars.size()>() );
    break;
  case LHCb::Track::Types::Upstream:
    netresponse = m_upstreamClassifier.GetMvaValue( LHCb::span{variables}.first<upstreamVars.size()>() );
    break;
  case LHCb::Track::Types::Ttrack:
    netresponse = m_tClassifier.GetMvaValue( LHCb::span{variables}.first<ttrackVars.size()>() );
    break;
  default:
    Error( "Track type: " + std::to_string( static_cast<int>( aTrack.type() ) ) +
               " not known to the ghost probability tool. Will set ghost prob value to inf",
           StatusCode::SUCCESS, 10 )
        .ignore();
  }

  // float netresponse = m_readers[aTrack.type()]->GetRarity(variables); // TODO rarity would be nice, see
  // https://sft.its.cern.ch/jira/browse/ROOT-7050
  netresponse = m_flatters[static_cast<int>( aTrack.type() )]->value( netresponse );
  aTrack.setGhostProbability( 1. - netresponse );

  return StatusCode::SUCCESS;
}

std::vector<std::string_view> UpgradeGhostId::variableNames( LHCb::Track::Types type ) const {
  switch ( type ) {
  case LHCb::Track::Types::Velo:
  case LHCb::Track::Types::VeloBackward:
    return std::vector<std::string_view>( veloVars.begin(), veloVars.end() );
  case LHCb::Track::Types::Long:
    return std::vector<std::string_view>( longVars.begin(), longVars.end() );
  case LHCb::Track::Types::Upstream:
    return std::vector<std::string_view>( upstreamVars.begin(), upstreamVars.end() );
  case LHCb::Track::Types::Downstream:
    return std::vector<std::string_view>( downstreamVars.begin(), downstreamVars.end() );
  case LHCb::Track::Types::Ttrack:
    return std::vector<std::string_view>( ttrackVars.begin(), ttrackVars.end() );
  default:
    return {};
  }
}
