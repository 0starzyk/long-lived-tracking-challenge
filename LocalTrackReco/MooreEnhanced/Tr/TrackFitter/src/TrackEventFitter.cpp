/*****************************************************************************\
* (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include "TrackInterfaces/ITrackFitter.h"
#include "TrackKernel/TrackFunctors.h"

#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/Track.h"
#include "Kernel/STLExtensions.h"

#include "Gaudi/Algorithm.h"
#include "GaudiAlg/FunctionalDetails.h"
#include "GaudiKernel/ToolHandle.h"
#include "LHCbAlgs/Transformer.h"

//-----------------------------------------------------------------------------
// Implementation file for class : TrackEventFitter
//
// 2005-05-30 : Eduardo Rodrigues
//-----------------------------------------------------------------------------

namespace {
  template <typename TrackListType>
  using MyBase =
      LHCb::Algorithm::Transformer<LHCb::Tracks( const TrackListType&, const DetectorElement& ),
                                   LHCb::DetDesc::usesBaseAndConditions<FixTESPath<Gaudi::Algorithm>, DetectorElement>>;
}

template <typename TrackListType>
class TrackEventFitter : public MyBase<TrackListType> {
public:
  using MyBase<TrackListType>::msgLevel;
  using MyBase<TrackListType>::info;
  using MyBase<TrackListType>::debug;
  using MyBase<TrackListType>::always;
  using MyBase<TrackListType>::verbose;
  using MyBase<TrackListType>::inputLocation;
  using MyBase<TrackListType>::outputLocation;

  /// Standard constructor
  TrackEventFitter( const std::string& name, ISvcLocator* pSvcLocator )
      : MyBase<TrackListType>(
            name, pSvcLocator,
            {typename MyBase<TrackListType>::KeyValue{"TracksInContainer", LHCb::TrackLocation::Default},
             typename MyBase<TrackListType>::KeyValue{"StandardGeometryTop", LHCb::standard_geometry_top}},
            {"TracksOutContainer", LHCb::TrackLocation::Default} ) {}

  LHCb::Tracks operator()( const TrackListType&, const DetectorElement& ) const override;
  /// resets the internal fitter
  StatusCode initialize() override { return MyBase<TrackListType>::initialize(); }
  StatusCode stop() override {
    m_tracksFitter->reset();
    return StatusCode::SUCCESS;
  }

private:
  /// interface to tracks fitter tool
  ToolHandle<ITrackFitter> m_tracksFitter{this, "Fitter", "TrackMasterFitter/Fitter"};

  Gaudi::Property<bool>     m_skipFailedFitAtInput{this, "SkipFailedFitAtInput", true};
  Gaudi::Property<double>   m_maxChi2DoF{this, "MaxChi2DoF", 9999999,
                                       "Max chi2 per track when output is a new container"};
  Gaudi::Property<unsigned> m_batchSize{this, "BatchSize", 50, "Size of batches"};

  mutable Gaudi::Accumulators::StatCounter<unsigned int> m_nTracksCount{this, "nTracks"};
  mutable Gaudi::Accumulators::StatCounter<unsigned int> m_nFittedCount{this, "nFitted"};
  mutable Gaudi::Accumulators::StatCounter<unsigned int> m_nBadInputCount{this, "nBadInput"};
  mutable Gaudi::Accumulators::StatCounter<double>       m_chisqprobSum{this, "chisqprobSum"};
  mutable Gaudi::Accumulators::StatCounter<unsigned int> m_badChisq{this, "badChisq"};
  mutable Gaudi::Accumulators::StatCounter<unsigned int> m_flipCharge{this, "flipCharge"};
  mutable Gaudi::Accumulators::StatCounter<unsigned int> m_numOutliers{this, "numOutliers"};
  mutable Gaudi::Accumulators::StatCounter<unsigned int> m_rejectedChisqCut{this, "rejectedChisqCut"};
  mutable Gaudi::Accumulators::BinomialCounter<>         m_fitPerf{this, "Fit Failure Rate"};
};

DECLARE_COMPONENT_WITH_ID( TrackEventFitter<std::vector<LHCb::Track>>, "VectorOfTracksFitter" )
DECLARE_COMPONENT_WITH_ID( TrackEventFitter<LHCb::Tracks>, "TrackEventFitter" )
DECLARE_COMPONENT_WITH_ID( TrackEventFitter<LHCb::Track::Selection>, "SharedTrackEventFitter" )

//=============================================================================
// Main execution
//=============================================================================

template <typename TrackListType>
LHCb::Tracks TrackEventFitter<TrackListType>::operator()( const TrackListType&   tracksCont,
                                                          const DetectorElement& lhcb ) const {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  LHCb::Tracks retTracks;
  retTracks.reserve( tracksCont.size() );

  // Loop over the tracks and fit them
  // ---------------------------------

  unsigned int nFitFail{0}, nBadInput{0};

  std::vector<LHCb::Track> tracks{};
  tracks.reserve( tracksCont.size() );
  std::vector<double> qopBeforeVector;
  qopBeforeVector.reserve( tracksCont.size() );

  std::for_each( std::begin( tracksCont ), std::end( tracksCont ), [&]( auto const& tr ) {
    auto const& track = Gaudi::Functional::details::deref( tr );
    if ( track.nStates() == 0 || track.checkFlag( LHCb::Track::Flags::Invalid ) ||
         ( m_skipFailedFitAtInput.value() && track.fitStatus() == LHCb::Track::FitStatus::FitFailed ) ) {
      // don't put failures on the output container. this is how they want it in HLT.
      ++nBadInput;
    } else {
      if ( track.hasKey() ) {
        tracks.emplace_back( track, track.key() );
      } else {
        tracks.emplace_back( track );
      }
      qopBeforeVector.push_back( track.firstState().qOverP() );
    }
  } );

  ( *m_tracksFitter )( tracks, *lhcb.geometry(), LHCb::Tr::PID::Pion() ).ignore();

  // counter buffers
  auto chisqprobSum     = m_chisqprobSum.buffer();
  auto badChisq         = m_badChisq.buffer();
  auto flipCharge       = m_flipCharge.buffer();
  auto numOutliers      = m_numOutliers.buffer();
  auto rejectedChisqCut = m_rejectedChisqCut.buffer();

  for ( auto&& [wrapper, qopBefore] : ranges::views::zip( tracks, qopBeforeVector ) ) {
    LHCb::Track& track = wrapper;

    if ( track.fitStatus() == LHCb::Track::FitStatus::Fitted ) {
      // Update counters
      if ( track.nDoF() > 0 ) {
        const auto chisqprob = track.probChi2();
        chisqprobSum += chisqprob;
        badChisq += ( chisqprob < 0.01 );
      }
      flipCharge += ( qopBefore * track.firstState().qOverP() < 0 );
      numOutliers += nMeasurementsRemoved( track );
      // Add the track to the new Tracks container if it passes the chi2-cut.
      // -----------------------------------------
      if ( m_maxChi2DoF > track.chi2PerDoF() ) {
        if ( track.hasKey() ) {
          retTracks.add( new LHCb::Track( std::move( track ), track.key() ) );
        } else {
          retTracks.add( new LHCb::Track( std::move( track ) ) );
        }
      } else {
        rejectedChisqCut += 1;
      }
    } else {
      ++nFitFail;
    }
  }

  // Update counters
  // ---------------
  const auto nTracks = tracksCont.size();
  const auto totFail = nFitFail + nBadInput;
  m_nTracksCount += nTracks;
  m_nFittedCount += ( nTracks - totFail );
  m_nBadInputCount += nBadInput;
  m_fitPerf += {totFail, nTracks};

  return retTracks;
}
