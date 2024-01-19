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
// Include files

// from Gaudi
#include "GaudiKernel/AnyDataHandle.h"
#include "LHCbAlgs/Consumer.h"

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"

#include "TrackInterfaces/ITrackExtrapolator.h"

#include "MCInterfaces/IIdealStateCreator.h"

#include "Associators/Associators.h"

#include <TTree.h>

#include "Event/MCParticle.h"
#include "Event/MCTrackInfo.h"
#include "Event/ODIN.h"
#include "Event/State.h"
#include "Event/Track.h"

#include <TFile.h>

/** @class CompareTracks CompareTracks.h
 *
 *  Algorithm that compares two fitted tracks to MC truth
 *
 *
 *  Parameters:
 *  - OutputFile:       Output location and file name for the
 *
 *  @author Simon Stemmle
 *  @date   2017-11-08
 */

// struct that contains variables for the tupling
struct tupleVars {
  int m_MC_status;

  double m_true_qop_vertex;

  std::array<std::array<double, 5>, 3>  m_One_x;
  std::array<std::array<double, 15>, 3> m_One_P;
  std::array<std::array<double, 5>, 3>  m_One_true_x;
  std::array<double, 3>                 m_One_z;
  double                                m_One_chi2;
  double                                m_One_ndof;

  std::array<std::array<double, 5>, 3>  m_Two_x;
  std::array<std::array<double, 15>, 3> m_Two_P;
  std::array<std::array<double, 5>, 3>  m_Two_true_x;
  std::array<double, 3>                 m_Two_z;
  double                                m_Two_chi2;
  double                                m_Two_ndof;
};

template <typename TrackListType1, typename TrackListType2>
using base_class =
    LHCb::Algorithm::Consumer<void( TrackListType1 const&, TrackListType2 const&, LHCb::ODIN const&,
                                    LHCb::LinksByKey const&, LHCb::MCProperty const&, DetectorElement const& ),
                              LHCb::DetDesc::usesConditions<DetectorElement>>;
template <typename TrackListType1, typename TrackListType2>
class CompareTracks : public base_class<TrackListType1, TrackListType2> {
public:
  using base_class<TrackListType1, TrackListType2>::msgLevel;
  using base_class<TrackListType1, TrackListType2>::debug;
  using base_class<TrackListType1, TrackListType2>::warning;
  using base_class<TrackListType1, TrackListType2>::error;
  using base_class<TrackListType1, TrackListType2>::info;
  using base_class<TrackListType1, TrackListType2>::msgSvc;
  using typename base_class<TrackListType1, TrackListType2>::KeyValue;
  /// Standard constructor
  CompareTracks( const std::string& name, ISvcLocator* pSvcLocator )
      : base_class<TrackListType1, TrackListType2>{
            name,
            pSvcLocator,
            {std::pair<std::string, std::string>{"InputTracks1", "Rec/Track/ForwardFastFitted_TMP"},
             std::pair<std::string, std::string>{"InputTracks2", "Rec/Track/ForwardFastFitted"},
             std::pair<std::string, std::string>{"ODINLocation", LHCb::ODINLocation::Default},
             std::pair<std::string, std::string>{"LinkerLocation", Links::location( "Rec/Track/ForwardFastFitted" )},
             std::pair<std::string, std::string>{"MCProperty", LHCb::MCPropertyLocation::TrackInfo},
             std::pair<std::string, std::string>{"StandardGeometryTop", LHCb::standard_geometry_top}}} {};

  /// Algorithm execution
  void operator()( TrackListType1 const& tracks1, TrackListType2 const& tracks2, LHCb::ODIN const& odin,
                   LHCb::LinksByKey const& links, LHCb::MCProperty const&, DetectorElement const& ) const override;

private:
  Gaudi::Property<std::string> m_FileName{this, "OutputFile", "CompareTracks"};

  // ideal state creator for tuning and performance checks
  ToolHandle<IIdealStateCreator> m_idealStateCreator = {"IdealStateCreator", this};

  // extrapolator
  ToolHandle<ITrackExtrapolator> m_extrapolator = {"TrackMasterExtrapolator/extr", this};

  //#################
  // 1. Level methods
  //#################

  /// Create trees that should be filled for tuning and perfomance checks
  void addBranches( TTree& trees, tupleVars* vars ) const;

  /// Fill information for the comparison of two tracks
  void FillNtuple( LHCb::Track const& track1, LHCb::Track const& track2, LHCb::LinksByKey const& links,
                   MCTrackInfo const&, tupleVars* vars, double z, int nPos, IGeometryInfo const& geometry,
                   bool closeToVertex = false ) const;

  //#######################################
  // Further methods for the Kalman filter
  //#######################################

  /// Check if a MC particle is linked to this track
  int MatchesMC( const LHCb::Track& track, const LHCb::LinksByKey& links, const MCTrackInfo& trackInfo ) const;

  /// Get true state at a given z position
  bool TrueState( double zpos, double& trueX, double& trueY, double& truetX, double& truetY, double& trueqop,
                  LHCb::Track const& track, LHCb::LinksByKey const& links, IGeometryInfo const& geometry,
                  bool initialQop = true, bool closeToVertex = false ) const;
};

// Declaration of the Algorithm Factory

typedef CompareTracks<std::vector<LHCb::Track>, std::vector<LHCb::Track>> CompareTracksVecVec;
typedef CompareTracks<std::vector<LHCb::Track>, LHCb::Tracks>             CompareTracksVecTr;
typedef CompareTracks<LHCb::Tracks, std::vector<LHCb::Track>>             CompareTracksTrVec;
typedef CompareTracks<LHCb::Tracks, LHCb::Tracks>                         CompareTracksTrTr;

DECLARE_COMPONENT_WITH_ID( CompareTracksVecVec, "CompareTracksVecVec" )
DECLARE_COMPONENT_WITH_ID( CompareTracksVecTr, "CompareTracksVecTr" )
DECLARE_COMPONENT_WITH_ID( CompareTracksTrVec, "CompareTracksTrVec" )
DECLARE_COMPONENT_WITH_ID( CompareTracksTrTr, "CompareTracksTrTr" )

//==================================================================================================
// Main execution
//==================================================================================================
template <typename TrackListType1, typename TrackListType2>
void CompareTracks<TrackListType1, TrackListType2>::
     operator()( TrackListType1 const& tracks1, TrackListType2 const& tracks2, LHCb::ODIN const& odin,
            LHCb::LinksByKey const& links, LHCb::MCProperty const& mcprop, DetectorElement const& lhcb ) const {
  const auto trackInfo = MCTrackInfo{mcprop};

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  if ( tracks1.empty() || tracks2.empty() ) { return; }

  std::vector<LHCb::Track> tracksV1{};
  std::vector<LHCb::Track> tracksV2{};
  tracksV1.reserve( tracks1.size() );
  tracksV2.reserve( tracks2.size() );

  std::for_each( std::begin( tracks1 ), std::end( tracks1 ), [&]( auto& track ) {
    if ( Gaudi::Functional::details::deref( track ).hasKey() ) {
      tracksV1.emplace_back( Gaudi::Functional::details::deref( track ),
                             Gaudi::Functional::details::deref( track ).key() );
    } else {
      tracksV1.emplace_back( Gaudi::Functional::details::deref( track ) );
    }
  } );

  std::for_each( std::begin( tracks2 ), std::end( tracks2 ), [&]( auto& track ) {
    if ( Gaudi::Functional::details::deref( track ).hasKey() ) {
      tracksV2.emplace_back( Gaudi::Functional::details::deref( track ),
                             Gaudi::Functional::details::deref( track ).key() );
    } else {
      tracksV2.emplace_back( Gaudi::Functional::details::deref( track ) );
    }
  } );

  // Create output tuple
  // Create a new file for every event in order to be threadsafe
  TFile outputFile( ( m_FileName.toString() + "_" + std::to_string( odin.runNumber() ) + "_" +
                      std::to_string( odin.eventNumber() ) + ".root" )
                        .c_str(),
                    "RECREATE" );

  // create the varaibles to be filled
  tupleVars vars;

  // create the tree
  TTree tree( "compare", "compare" );
  addBranches( tree, &vars );

  // Loop over tracks1
  for ( auto const& track1 : tracksV1 ) {
    // Search for the respective track in tracks2
    auto i = std::find_if( std::begin( tracksV2 ), std::end( tracksV2 ),
                           [&]( auto t ) { return track1.nCommonLhcbIDs( t ) == track1.nLHCbIDs(); } );
    if ( i == std::end( tracksV2 ) ) {
      warning() << "No matching track in second container found" << endmsg;
      return;
    }

    // compare the tracks and fill the tuple

    // 1.
    // Near the beam pipe

    // Get the z position of the state near the beam pipe
    double             zBeam  = 0;
    const LHCb::State* state1 = track1.stateAt( LHCb::State::Location::ClosestToBeam );
    const LHCb::State* state2 = ( *i ).stateAt( LHCb::State::Location::ClosestToBeam );
    if ( state1 != nullptr && state2 != nullptr ) zBeam = 0.5 * ( state1->z() + state2->z() );

    FillNtuple( track1, *( i ), links, trackInfo, &vars, zBeam, 0, *lhcb.geometry(), true );

    // 2.
    // Get the z position of the first measurement
    double             zFirstMeas = 0;
    const LHCb::State* state1FM   = track1.stateAt( LHCb::State::Location::FirstMeasurement );
    const LHCb::State* state2FM   = ( *i ).stateAt( LHCb::State::Location::FirstMeasurement );
    if ( state1FM != nullptr && state2FM != nullptr )
      zFirstMeas = 0.5 * ( state1FM->z() + state2FM->z() );
    else if ( state1FM != nullptr )
      zFirstMeas = state1FM->z();
    else if ( state2FM != nullptr )
      zFirstMeas = state2FM->z();

    FillNtuple( track1, *( i ), links, trackInfo, &vars, zFirstMeas, 1, *lhcb.geometry(), false );

    // 3.
    // Get the z position of the last measurement
    double             zLastMeas = 0;
    const LHCb::State* state1LM  = track1.stateAt( LHCb::State::Location::LastMeasurement );
    const LHCb::State* state2LM  = ( *i ).stateAt( LHCb::State::Location::LastMeasurement );
    if ( state1LM != nullptr && state2LM != nullptr )
      zLastMeas = 0.5 * ( state1LM->z() + state2LM->z() );
    else if ( state1LM != nullptr )
      zLastMeas = state1LM->z();
    else if ( state2LM != nullptr )
      zLastMeas = state2LM->z();

    FillNtuple( track1, *( i ), links, trackInfo, &vars, zLastMeas, 2, *lhcb.geometry(), false );

    tree.Fill();
    if ( msgLevel( MSG::DEBUG ) ) {
      debug() << "Compare qop: " << vars.m_One_x[0][4] << " " << vars.m_One_x[1][4] << " " << vars.m_One_x[2][4]
              << endmsg;
      debug() << "         vs: " << vars.m_Two_x[0][4] << " " << vars.m_Two_x[1][4] << " " << vars.m_Two_x[2][4]
              << endmsg;
      debug() << "         vs: " << vars.m_Two_true_x[0][4] << " " << vars.m_Two_true_x[1][4] << " "
              << vars.m_Two_true_x[2][4] << endmsg;
    }
  }

  outputFile.Write();
  outputFile.Close();

  return;
}

//==================================================================================================
//  Add branches to the tuple that will be filled with the comparison information
//==================================================================================================
template <typename TrackListType1, typename TrackListType2>
void CompareTracks<TrackListType1, TrackListType2>::addBranches( TTree& tree, tupleVars* treeVars ) const {

  // Set the branches
  tree.Branch( "One_sigmaxx_T", &( treeVars->m_One_P[2][0] ), "One_sigmaxx_T/D" );
  tree.Branch( "One_sigmayy_T", &( treeVars->m_One_P[2][2] ), "One_sigmayy_T/D" );
  tree.Branch( "One_sigmatxtx_T", &( treeVars->m_One_P[2][5] ), "One_sigmatxtx_T/D" );
  tree.Branch( "One_sigmatyty_T", &( treeVars->m_One_P[2][9] ), "One_sigmatyty_T/D" );
  tree.Branch( "One_sigmaqopqop_T", &( treeVars->m_One_P[2][14] ), "One_sigmaqopqop_T/D" );
  tree.Branch( "One_sigmaxy_T", &( treeVars->m_One_P[2][1] ), "One_sigmaxy_T/D" );
  tree.Branch( "One_sigmaxtx_T", &( treeVars->m_One_P[2][3] ), "One_sigmaxtx_T/D" );
  tree.Branch( "One_sigmaxty_T", &( treeVars->m_One_P[2][6] ), "One_sigmaxty_T/D" );
  tree.Branch( "One_sigmaxqop_T", &( treeVars->m_One_P[2][10] ), "One_sigmaxqop_T/D" );
  tree.Branch( "One_sigmaytx_T", &( treeVars->m_One_P[2][4] ), "One_sigmaytx_T/D" );
  tree.Branch( "One_sigmayty_T", &( treeVars->m_One_P[2][7] ), "One_sigmayty_T/D" );
  tree.Branch( "One_sigmayqop_T", &( treeVars->m_One_P[2][11] ), "One_sigmayqop_T/D" );
  tree.Branch( "One_sigmatxty_T", &( treeVars->m_One_P[2][8] ), "One_sigmatxty_T/D" );
  tree.Branch( "One_sigmatxqop_T", &( treeVars->m_One_P[2][12] ), "One_sigmatxqop_T/D" );
  tree.Branch( "One_sigmatyqop_T", &( treeVars->m_One_P[2][13] ), "One_sigmatyqop_T/D" );

  tree.Branch( "One_x_T", &( treeVars->m_One_x[2][0] ), "One_x_T/D" );
  tree.Branch( "One_y_T", &( treeVars->m_One_x[2][1] ), "One_y_T/D" );
  tree.Branch( "One_tx_T", &( treeVars->m_One_x[2][2] ), "One_tx_T/D" );
  tree.Branch( "One_ty_T", &( treeVars->m_One_x[2][3] ), "One_ty_T/D" );
  tree.Branch( "One_qop_T", &( treeVars->m_One_x[2][4] ), "One_qop_T/D" );

  tree.Branch( "One_z_T", &( treeVars->m_One_z[2] ), "One_z_T/D" );

  tree.Branch( "One_true_x_T", &( treeVars->m_One_true_x[2][0] ), "One_true_x_T/D" );
  tree.Branch( "One_true_y_T", &( treeVars->m_One_true_x[2][1] ), "One_true_y_T/D" );
  tree.Branch( "One_true_tx_T", &( treeVars->m_One_true_x[2][2] ), "One_true_tx_T/D" );
  tree.Branch( "One_true_ty_T", &( treeVars->m_One_true_x[2][3] ), "One_true_ty_T/D" );
  tree.Branch( "One_true_qop_T", &( treeVars->m_One_true_x[2][4] ), "One_true_qop_T/D" );

  tree.Branch( "One_sigmaxx_V", &( treeVars->m_One_P[0][0] ), "One_sigmaxx_V/D" );
  tree.Branch( "One_sigmayy_V", &( treeVars->m_One_P[0][2] ), "One_sigmayy_V/D" );
  tree.Branch( "One_sigmatxtx_V", &( treeVars->m_One_P[0][5] ), "One_sigmatxtx_V/D" );
  tree.Branch( "One_sigmatyty_V", &( treeVars->m_One_P[0][9] ), "One_sigmatyty_V/D" );
  tree.Branch( "One_sigmaqopqop_V", &( treeVars->m_One_P[0][14] ), "One_sigmaqopqop_V/D" );
  tree.Branch( "One_sigmaxy_V", &( treeVars->m_One_P[0][1] ), "One_sigmaxy_V/D" );
  tree.Branch( "One_sigmaxtx_V", &( treeVars->m_One_P[0][3] ), "One_sigmaxtx_V/D" );
  tree.Branch( "One_sigmaxty_V", &( treeVars->m_One_P[0][6] ), "One_sigmaxty_V/D" );
  tree.Branch( "One_sigmaxqop_V", &( treeVars->m_One_P[0][10] ), "One_sigmaxqop_V/D" );
  tree.Branch( "One_sigmaytx_V", &( treeVars->m_One_P[0][4] ), "One_sigmaytx_V/D" );
  tree.Branch( "One_sigmayty_V", &( treeVars->m_One_P[0][7] ), "One_sigmayty_V/D" );
  tree.Branch( "One_sigmayqop_V", &( treeVars->m_One_P[0][11] ), "One_sigmayqop_V/D" );
  tree.Branch( "One_sigmatxty_V", &( treeVars->m_One_P[0][8] ), "One_sigmatxty_V/D" );
  tree.Branch( "One_sigmatxqop_V", &( treeVars->m_One_P[0][12] ), "One_sigmatxqop_V/D" );
  tree.Branch( "One_sigmatyqop_V", &( treeVars->m_One_P[0][13] ), "One_sigmatyqop_V/D" );

  tree.Branch( "One_x_V", &( treeVars->m_One_x[0][0] ), "One_x_V/D" );
  tree.Branch( "One_y_V", &( treeVars->m_One_x[0][1] ), "One_y_V/D" );
  tree.Branch( "One_tx_V", &( treeVars->m_One_x[0][2] ), "One_tx_V/D" );
  tree.Branch( "One_ty_V", &( treeVars->m_One_x[0][3] ), "One_ty_V/D" );
  tree.Branch( "One_qop_V", &( treeVars->m_One_x[0][4] ), "One_qop_V/D" );

  tree.Branch( "One_z_V", &( treeVars->m_One_z[0] ), "One_z_V/D" );

  tree.Branch( "One_true_x_V", &( treeVars->m_One_true_x[0][0] ), "One_true_x_V/D" );
  tree.Branch( "One_true_y_V", &( treeVars->m_One_true_x[0][1] ), "One_true_y_V/D" );
  tree.Branch( "One_true_tx_V", &( treeVars->m_One_true_x[0][2] ), "One_true_tx_V/D" );
  tree.Branch( "One_true_ty_V", &( treeVars->m_One_true_x[0][3] ), "One_true_ty_V/D" );
  tree.Branch( "One_true_qop_V", &( treeVars->m_One_true_x[0][4] ), "One_true_qop_V/D" );

  tree.Branch( "One_sigmaxx", &( treeVars->m_One_P[1][0] ), "One_sigmaxx/D" );
  tree.Branch( "One_sigmayy", &( treeVars->m_One_P[1][2] ), "One_sigmayy/D" );
  tree.Branch( "One_sigmatxtx", &( treeVars->m_One_P[1][5] ), "One_sigmatxtx/D" );
  tree.Branch( "One_sigmatyty", &( treeVars->m_One_P[1][9] ), "One_sigmatyty/D" );
  tree.Branch( "One_sigmaqopqop", &( treeVars->m_One_P[1][14] ), "One_sigmaqopqop/D" );
  tree.Branch( "One_sigmaxy", &( treeVars->m_One_P[1][1] ), "One_sigmaxy/D" );
  tree.Branch( "One_sigmaxtx", &( treeVars->m_One_P[1][3] ), "One_sigmaxtx/D" );
  tree.Branch( "One_sigmaxty", &( treeVars->m_One_P[1][6] ), "One_sigmaxty/D" );
  tree.Branch( "One_sigmaxqop", &( treeVars->m_One_P[1][10] ), "One_sigmaxqop/D" );
  tree.Branch( "One_sigmaytx", &( treeVars->m_One_P[1][4] ), "One_sigmaytx/D" );
  tree.Branch( "One_sigmayty", &( treeVars->m_One_P[1][7] ), "One_sigmayty/D" );
  tree.Branch( "One_sigmayqop", &( treeVars->m_One_P[1][11] ), "One_sigmayqop/D" );
  tree.Branch( "One_sigmatxty", &( treeVars->m_One_P[1][8] ), "One_sigmatxty/D" );
  tree.Branch( "One_sigmatxqop", &( treeVars->m_One_P[1][12] ), "One_sigmatxqop/D" );
  tree.Branch( "One_sigmatyqop", &( treeVars->m_One_P[1][13] ), "One_sigmatyqop/D" );

  tree.Branch( "One_x", &( treeVars->m_One_x[1][0] ), "One_x/D" );
  tree.Branch( "One_y", &( treeVars->m_One_x[1][1] ), "One_y/D" );
  tree.Branch( "One_tx", &( treeVars->m_One_x[1][2] ), "One_tx/D" );
  tree.Branch( "One_ty", &( treeVars->m_One_x[1][3] ), "One_ty/D" );
  tree.Branch( "One_qop", &( treeVars->m_One_x[1][4] ), "One_qop/D" );

  tree.Branch( "One_z", &( treeVars->m_One_z[1] ), "One_z/D" );
  tree.Branch( "One_chi2", &( treeVars->m_One_chi2 ), "One_chi2/D" );
  tree.Branch( "One_ndof", &( treeVars->m_One_ndof ), "One_ndof/D" );

  tree.Branch( "One_true_x", &( treeVars->m_One_true_x[1][0] ), "One_true_x/D" );
  tree.Branch( "One_true_y", &( treeVars->m_One_true_x[1][1] ), "One_true_y/D" );
  tree.Branch( "One_true_tx", &( treeVars->m_One_true_x[1][2] ), "One_true_tx/D" );
  tree.Branch( "One_true_ty", &( treeVars->m_One_true_x[1][3] ), "One_true_ty/D" );
  tree.Branch( "One_true_qop", &( treeVars->m_One_true_x[1][4] ), "One_true_qop/D" );

  tree.Branch( "Two_sigmaxx_T", &( treeVars->m_Two_P[2][0] ), "Two_sigmaxx_T/D" );
  tree.Branch( "Two_sigmayy_T", &( treeVars->m_Two_P[2][2] ), "Two_sigmayy_T/D" );
  tree.Branch( "Two_sigmatxtx_T", &( treeVars->m_Two_P[2][5] ), "Two_sigmatxtx_T/D" );
  tree.Branch( "Two_sigmatyty_T", &( treeVars->m_Two_P[2][9] ), "Two_sigmatyty_T/D" );
  tree.Branch( "Two_sigmaqopqop_T", &( treeVars->m_Two_P[2][14] ), "Two_sigmaqopqop_T/D" );
  tree.Branch( "Two_sigmaxy_T", &( treeVars->m_Two_P[2][1] ), "Two_sigmaxy_T/D" );
  tree.Branch( "Two_sigmaxtx_T", &( treeVars->m_Two_P[2][3] ), "Two_sigmaxtx_T/D" );
  tree.Branch( "Two_sigmaxty_T", &( treeVars->m_Two_P[2][6] ), "Two_sigmaxty_T/D" );
  tree.Branch( "Two_sigmaxqop_T", &( treeVars->m_Two_P[2][10] ), "Two_sigmaxqop_T/D" );
  tree.Branch( "Two_sigmaytx_T", &( treeVars->m_Two_P[2][4] ), "Two_sigmaytx_T/D" );
  tree.Branch( "Two_sigmayty_T", &( treeVars->m_Two_P[2][7] ), "Two_sigmayty_T/D" );
  tree.Branch( "Two_sigmayqop_T", &( treeVars->m_Two_P[2][11] ), "Two_sigmayqop_T/D" );
  tree.Branch( "Two_sigmatxty_T", &( treeVars->m_Two_P[2][8] ), "Two_sigmatxty_T/D" );
  tree.Branch( "Two_sigmatxqop_T", &( treeVars->m_Two_P[2][12] ), "Two_sigmatxqop_T/D" );
  tree.Branch( "Two_sigmatyqop_T", &( treeVars->m_Two_P[2][13] ), "Two_sigmatyqop_T/D" );

  tree.Branch( "Two_x_T", &( treeVars->m_Two_x[2][0] ), "Two_x_T/D" );
  tree.Branch( "Two_y_T", &( treeVars->m_Two_x[2][1] ), "Two_y_T/D" );
  tree.Branch( "Two_tx_T", &( treeVars->m_Two_x[2][2] ), "Two_tx_T/D" );
  tree.Branch( "Two_ty_T", &( treeVars->m_Two_x[2][3] ), "Two_ty_T/D" );
  tree.Branch( "Two_qop_T", &( treeVars->m_Two_x[2][4] ), "Two_qop_T/D" );

  tree.Branch( "Two_z_T", &( treeVars->m_Two_z[2] ), "Two_z_T/D" );

  tree.Branch( "Two_true_x_T", &( treeVars->m_Two_true_x[2][0] ), "Two_true_x_T/D" );
  tree.Branch( "Two_true_y_T", &( treeVars->m_Two_true_x[2][1] ), "Two_true_y_T/D" );
  tree.Branch( "Two_true_tx_T", &( treeVars->m_Two_true_x[2][2] ), "Two_true_tx_T/D" );
  tree.Branch( "Two_true_ty_T", &( treeVars->m_Two_true_x[2][3] ), "Two_true_ty_T/D" );
  tree.Branch( "Two_true_qop_T", &( treeVars->m_Two_true_x[2][4] ), "Two_true_qop_T/D" );

  tree.Branch( "Two_sigmaxx_V", &( treeVars->m_Two_P[0][0] ), "Two_sigmaxx_V/D" );
  tree.Branch( "Two_sigmayy_V", &( treeVars->m_Two_P[0][2] ), "Two_sigmayy_V/D" );
  tree.Branch( "Two_sigmatxtx_V", &( treeVars->m_Two_P[0][5] ), "Two_sigmatxtx_V/D" );
  tree.Branch( "Two_sigmatyty_V", &( treeVars->m_Two_P[0][9] ), "Two_sigmatyty_V/D" );
  tree.Branch( "Two_sigmaqopqop_V", &( treeVars->m_Two_P[0][14] ), "Two_sigmaqopqop_V/D" );
  tree.Branch( "Two_sigmaxy_V", &( treeVars->m_Two_P[0][1] ), "Two_sigmaxy_V/D" );
  tree.Branch( "Two_sigmaxtx_V", &( treeVars->m_Two_P[0][3] ), "Two_sigmaxtx_V/D" );
  tree.Branch( "Two_sigmaxty_V", &( treeVars->m_Two_P[0][6] ), "Two_sigmaxty_V/D" );
  tree.Branch( "Two_sigmaxqop_V", &( treeVars->m_Two_P[0][10] ), "Two_sigmaxqop_V/D" );
  tree.Branch( "Two_sigmaytx_V", &( treeVars->m_Two_P[0][4] ), "Two_sigmaytx_V/D" );
  tree.Branch( "Two_sigmayty_V", &( treeVars->m_Two_P[0][7] ), "Two_sigmayty_V/D" );
  tree.Branch( "Two_sigmayqop_V", &( treeVars->m_Two_P[0][11] ), "Two_sigmayqop_V/D" );
  tree.Branch( "Two_sigmatxty_V", &( treeVars->m_Two_P[0][8] ), "Two_sigmatxty_V/D" );
  tree.Branch( "Two_sigmatxqop_V", &( treeVars->m_Two_P[0][12] ), "Two_sigmatxqop_V/D" );
  tree.Branch( "Two_sigmatyqop_V", &( treeVars->m_Two_P[0][13] ), "Two_sigmatyqop_V/D" );

  tree.Branch( "Two_x_V", &( treeVars->m_Two_x[0][0] ), "Two_x_V/D" );
  tree.Branch( "Two_y_V", &( treeVars->m_Two_x[0][1] ), "Two_y_V/D" );
  tree.Branch( "Two_tx_V", &( treeVars->m_Two_x[0][2] ), "Two_tx_V/D" );
  tree.Branch( "Two_ty_V", &( treeVars->m_Two_x[0][3] ), "Two_ty_V/D" );
  tree.Branch( "Two_qop_V", &( treeVars->m_Two_x[0][4] ), "Two_qop_V/D" );

  tree.Branch( "Two_z_V", &( treeVars->m_Two_z[0] ), "Two_z_V/D" );

  tree.Branch( "Two_true_x_V", &( treeVars->m_Two_true_x[0][0] ), "Two_true_x_V/D" );
  tree.Branch( "Two_true_y_V", &( treeVars->m_Two_true_x[0][1] ), "Two_true_y_V/D" );
  tree.Branch( "Two_true_tx_V", &( treeVars->m_Two_true_x[0][2] ), "Two_true_tx_V/D" );
  tree.Branch( "Two_true_ty_V", &( treeVars->m_Two_true_x[0][3] ), "Two_true_ty_V/D" );
  tree.Branch( "Two_true_qop_V", &( treeVars->m_Two_true_x[0][4] ), "Two_true_qop_V/D" );

  tree.Branch( "Two_sigmaxx", &( treeVars->m_Two_P[1][0] ), "Two_sigmaxx/D" );
  tree.Branch( "Two_sigmayy", &( treeVars->m_Two_P[1][2] ), "Two_sigmayy/D" );
  tree.Branch( "Two_sigmatxtx", &( treeVars->m_Two_P[1][5] ), "Two_sigmatxtx/D" );
  tree.Branch( "Two_sigmatyty", &( treeVars->m_Two_P[1][9] ), "Two_sigmatyty/D" );
  tree.Branch( "Two_sigmaqopqop", &( treeVars->m_Two_P[1][14] ), "Two_sigmaqopqop/D" );
  tree.Branch( "Two_sigmaxy", &( treeVars->m_Two_P[1][1] ), "Two_sigmaxy/D" );
  tree.Branch( "Two_sigmaxtx", &( treeVars->m_Two_P[1][3] ), "Two_sigmaxtx/D" );
  tree.Branch( "Two_sigmaxty", &( treeVars->m_Two_P[1][6] ), "Two_sigmaxty/D" );
  tree.Branch( "Two_sigmaxqop", &( treeVars->m_Two_P[1][10] ), "Two_sigmaxqop/D" );
  tree.Branch( "Two_sigmaytx", &( treeVars->m_Two_P[1][4] ), "Two_sigmaytx/D" );
  tree.Branch( "Two_sigmayty", &( treeVars->m_Two_P[1][7] ), "Two_sigmayty/D" );
  tree.Branch( "Two_sigmayqop", &( treeVars->m_Two_P[1][11] ), "Two_sigmayqop/D" );
  tree.Branch( "Two_sigmatxty", &( treeVars->m_Two_P[1][8] ), "Two_sigmatxty/D" );
  tree.Branch( "Two_sigmatxqop", &( treeVars->m_Two_P[1][12] ), "Two_sigmatxqop/D" );
  tree.Branch( "Two_sigmatyqop", &( treeVars->m_Two_P[1][13] ), "Two_sigmatyqop/D" );

  tree.Branch( "Two_x", &( treeVars->m_Two_x[1][0] ), "Two_x/D" );
  tree.Branch( "Two_y", &( treeVars->m_Two_x[1][1] ), "Two_y/D" );
  tree.Branch( "Two_tx", &( treeVars->m_Two_x[1][2] ), "Two_tx/D" );
  tree.Branch( "Two_ty", &( treeVars->m_Two_x[1][3] ), "Two_ty/D" );
  tree.Branch( "Two_qop", &( treeVars->m_Two_x[1][4] ), "Two_qop/D" );

  tree.Branch( "Two_z", &( treeVars->m_Two_z[1] ), "Two_z/D" );
  tree.Branch( "Two_chi2", &( treeVars->m_Two_chi2 ), "Two_chi2/D" );
  tree.Branch( "Two_ndof", &( treeVars->m_Two_ndof ), "Two_ndof/D" );

  tree.Branch( "Two_true_x", &( treeVars->m_Two_true_x[1][0] ), "Two_true_x/D" );
  tree.Branch( "Two_true_y", &( treeVars->m_Two_true_x[1][1] ), "Two_true_y/D" );
  tree.Branch( "Two_true_tx", &( treeVars->m_Two_true_x[1][2] ), "Two_true_tx/D" );
  tree.Branch( "Two_true_ty", &( treeVars->m_Two_true_x[1][3] ), "Two_true_ty/D" );
  tree.Branch( "Two_true_qop", &( treeVars->m_Two_true_x[1][4] ), "Two_true_qop/D" );

  tree.Branch( "true_qop_vertex", &( treeVars->m_true_qop_vertex ), "true_qop_vertex/D" );

  tree.Branch( "MCstatus", &( treeVars->m_MC_status ), "MCstatus/I" );
}

//==================================================================================================
//  Get the states informations at position z and fill the tree variables
//==================================================================================================
template <typename TrackListType1, typename TrackListType2>
void CompareTracks<TrackListType1, TrackListType2>::FillNtuple( LHCb::Track const& track1, LHCb::Track const& track2,
                                                                LHCb::LinksByKey const& links,
                                                                MCTrackInfo const& trackInfo, tupleVars* vars, double z,
                                                                int nPos, IGeometryInfo const& geometry,
                                                                bool closeToVertex ) const {
  vars->m_MC_status = MatchesMC( track1, links, trackInfo );

  // Get the states closest to the desired z position
  LHCb::State state1 = track1.closestState( z );
  LHCb::State state2 = track2.closestState( z );

  // extrapolate to the exact z position
  StatusCode sc1 = m_extrapolator->propagate( state1, z, geometry );
  StatusCode sc2 = m_extrapolator->propagate( state2, z, geometry );

  // covariance matrices
  Gaudi::TrackSymMatrix covMat1;
  covMat1 = state1.covariance();

  Gaudi::TrackSymMatrix covMat2;
  covMat2 = state2.covariance();

  if ( !sc1.isSuccess() || !sc2.isSuccess() ) return;

  // Get the true position at this z position
  double trueX( 0 ), trueY( 0 ), trueTX( 0 ), trueTY( 0 ), trueQoP( 0 ), trueQoP_Vertex( 0 );
  // Get the QoP at the vertex
  TrueState( z, trueX, trueY, trueTX, trueTY, trueQoP_Vertex, track1, links, geometry, true, closeToVertex );
  // Get the actual state at the z position
  TrueState( z, trueX, trueY, trueTX, trueTY, trueQoP, track1, links, geometry, false, closeToVertex );

  // Set the state variables
  vars->m_One_x[nPos][0] = state1.x();
  vars->m_One_x[nPos][1] = state1.y();
  vars->m_One_x[nPos][2] = state1.tx();
  vars->m_One_x[nPos][3] = state1.ty();
  vars->m_One_x[nPos][4] = state1.qOverP();

  vars->m_One_true_x[nPos][0] = trueX;
  vars->m_One_true_x[nPos][1] = trueY;
  vars->m_One_true_x[nPos][2] = trueTX;
  vars->m_One_true_x[nPos][3] = trueTY;
  vars->m_One_true_x[nPos][4] = trueQoP;

  vars->m_Two_x[nPos][0] = state2.x();
  vars->m_Two_x[nPos][1] = state2.y();
  vars->m_Two_x[nPos][2] = state2.tx();
  vars->m_Two_x[nPos][3] = state2.ty();
  vars->m_Two_x[nPos][4] = state2.qOverP();

  vars->m_Two_true_x[nPos][0] = trueX;
  vars->m_Two_true_x[nPos][1] = trueY;
  vars->m_Two_true_x[nPos][2] = trueTX;
  vars->m_Two_true_x[nPos][3] = trueTY;
  vars->m_Two_true_x[nPos][4] = trueQoP;

  vars->m_true_qop_vertex = trueQoP_Vertex;

  vars->m_One_z[nPos] = z;

  vars->m_Two_z[nPos] = z;

  int k = 0;
  for ( int i = 0; i < 5; i++ ) {
    for ( int j = 0; j <= i; j++ ) {
      vars->m_One_P[nPos][k] = covMat1( i, j );
      vars->m_Two_P[nPos][k] = covMat2( i, j );
      k++;
    }
  }

  // Set other track varaibles
  vars->m_One_ndof = track1.nDoF();
  vars->m_One_chi2 = track1.chi2();

  vars->m_Two_ndof = track2.nDoF();
  vars->m_Two_chi2 = track2.chi2();
}

//==================================================================================================
// Check if a MC particle is linked to this track
//==================================================================================================
template <typename TrackListType1, typename TrackListType2>
int CompareTracks<TrackListType1, TrackListType2>::MatchesMC( const LHCb::Track& track, const LHCb::LinksByKey& links,
                                                              const MCTrackInfo& trackInfo ) const {
  InputLinks<ContainedObject, LHCb::MCParticle> TrackParticleLinks( links );
  // Look for an associated MC particle
  auto trackLinks = TrackParticleLinks.from( track.key() );
  if ( trackLinks.empty() ) {
    debug() << "No links for track key " << track.key() << endmsg;
    return 0;
  }
  auto mcpart = std::max_element( trackLinks.begin(), trackLinks.end(), [&]( const auto& a, const auto& b ) {
                  return a.to()->momentum().P() < b.to()->momentum().P();
                } )->to();
  if ( !mcpart ) return 0;

  // check quality of matching
  if ( 0 == trackInfo.fullInfo( mcpart ) ) return 2;
  bool isLong = trackInfo.hasVeloAndT( mcpart );
  isLong      = isLong && ( abs( mcpart->particleID().pid() ) != 11 ); // and not electron
  if ( !isLong ) return 2;
  bool eta25 = ( mcpart->momentum().Eta() > 1.8 && mcpart->momentum().Eta() < 5.3 );
  if ( !eta25 ) return 2;

  if ( std::fabs( track.pseudoRapidity() - mcpart->momentum().Eta() ) > 0.05 ) return 2;
  return 1;
}

//==================================================================================================
// Get true state at a given z position
//==================================================================================================
template <typename TrackListType1, typename TrackListType2>
bool CompareTracks<TrackListType1, TrackListType2>::TrueState( double zpos, double& trueX, double& trueY,
                                                               double& truetX, double& truetY, double& trueqop,
                                                               LHCb::Track const& track, LHCb::LinksByKey const& links,
                                                               IGeometryInfo const& geometry, bool initialQop,
                                                               bool closeToVertex ) const {
  InputLinks<ContainedObject, LHCb::MCParticle> TrackParticleLinks( links );
  // Look for an associated MC particle
  auto trackLinks = TrackParticleLinks.from( track.key() );
  if ( trackLinks.empty() ) {
    debug() << "No links for track key " << track.key() << endmsg;
    return 0;
  }
  auto mcpart = std::max_element( trackLinks.begin(), trackLinks.end(), [&]( const auto& a, const auto& b ) {
                  return a.to()->momentum().P() < b.to()->momentum().P();
                } )->to();
  if ( !mcpart ) return false;

  LHCb::State state;
  // create the true state from the MC hits using the ideal state creator
  if ( !closeToVertex ) {
    StatusCode sc = m_idealStateCreator->createState( mcpart, zpos, state, geometry );
    if ( !sc.isSuccess() ) error() << "No ideal state could be created" << endmsg;
  }
  // use the MCParticle information to get the true state near the vertex
  else {
    StatusCode sc = m_idealStateCreator->createStateVertex( mcpart, state );
    if ( !sc.isSuccess() ) error() << "No ideal state could be created" << endmsg;
    m_extrapolator->propagate( state, zpos, geometry ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
  }
  trueX  = state.x();
  trueY  = state.y();
  truetX = state.tx();
  truetY = state.ty();

  if ( !initialQop )
    trueqop = state.qOverP();
  else
    trueqop = mcpart->particleID().threeCharge() * 1. / 3 * 1. / mcpart->momentum().P();

  return true;
}
