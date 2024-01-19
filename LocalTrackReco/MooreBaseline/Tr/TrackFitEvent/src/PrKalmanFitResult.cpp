/*****************************************************************************\
* (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "Event/PrKalmanFitResult.h"

namespace LHCb {
  std::unique_ptr<LHCb::ITrackFitResult> LHCb::PrKalmanFitResult::clone() const {
    return std::unique_ptr<PrKalmanFitResult>( new PrKalmanFitResult{*this} );
  }

  unsigned int PrKalmanFitResult::nActiveMeasurements() const {
    return std::count_if( fitnodes.begin(), fitnodes.end(),
                          [&]( const Pr::Tracks::Fit::Node& node ) { return node.isHitOnTrack(); } );
  }
  ChiSquare PrKalmanFitResult::chi2() const {
    auto chi2_bkwd = ChiSquare{0, -5};
    auto chi2_fwd  = ChiSquare{0, -5};
    for ( auto const& node : fitnodes ) {
      chi2_bkwd += node.delta_chi2[Pr::Tracks::Fit::Node::backward];
      chi2_fwd += node.delta_chi2[Pr::Tracks::Fit::Node::forward];
    }
    return chi2_fwd.chi2() < chi2_bkwd.chi2() ? chi2_fwd : chi2_bkwd;
  }

  ChiSquare PrKalmanFitResult::chi2Velo() const {
    auto velo_chi2 = ChiSquare{0, -5};
    for ( auto const& node : fitnodes ) {
      if ( node.type() == Pr::Tracks::Fit::Node::Type::VPHit ) {
        velo_chi2 += node.delta_chi2[Pr::Tracks::Fit::Node::backward];
      }
    }
    return velo_chi2;
  }

  ChiSquare PrKalmanFitResult::chi2Upstream() const {
    auto upstream_chi2 = ChiSquare{};
    for ( auto const& node : fitnodes ) {
      if ( node.type() == Pr::Tracks::Fit::Node::Type::UTHit ) {
        upstream_chi2 += node.delta_chi2[Pr::Tracks::Fit::Node::backward];
      }
    }
    upstream_chi2 += chi2Velo();
    return upstream_chi2;
  }

  ChiSquare PrKalmanFitResult::chi2Downstream() const {
    auto down_chi2 = ChiSquare{0, -5};
    for ( auto const& node : fitnodes ) {
      if ( node.type() == Pr::Tracks::Fit::Node::Type::FTHit ) {
        down_chi2 += node.delta_chi2[Pr::Tracks::Fit::Node::forward];
      }
    }
    return down_chi2;
  }

  ChiSquare PrKalmanFitResult::chi2Match() const { return chi2() - chi2Upstream() - chi2Downstream(); }

  /// TODO PrKalmanFilter does not support muon tracks at the moment
  ChiSquare PrKalmanFitResult::chi2Muon() const { return ChiSquare{}; }
} // namespace LHCb
