/*****************************************************************************\
* (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once

#include "Event/ITrackFitResult.h"
#include "Event/PrFitNode.h"
#include "Event/TrackTypes.h"
#include "Kernel/LHCbID.h"
#include "Kernel/STLExtensions.h"
#include <vector>

namespace LHCb {
  struct PrKalmanFitResult final : public ITrackFitResult {
    using NodeType = Pr::Tracks::Fit::Node;
    std::unique_ptr<ITrackFitResult> clone() const override;

    std::vector<Pr::Tracks::Fit::Node> fitnodes{};
    std::vector<Gaudi::TrackMatrix>    gain_matrices{};
    double                             scattering_momentum =
        std::numeric_limits<double>::signaling_NaN(); ///< Momentum used for computing material corrections
    int number_of_iter = -1;                          ///< Number of iterations in track fit

    /// For common interface with TrackFitResult
    /// Momentum used for computing material corrections
    double pScatter() const { return scattering_momentum; }
    /// For common interface with TrackFitResult
    /// Number of iterations in track fit
    auto nIter() const { return number_of_iter; }

    /// Get the number of active (non-outlier) measurements on the track
    unsigned int nActiveMeasurements() const;

    /// Track chi square obtained by summing delta chi2 from all fit nodes.
    /// Gives real chi2 only if the contributions are uncorrelated.
    LHCb::ChiSquare chi2() const;

    /// VELO segment chi square obtained by summing delta chi2 from all VPHit fit nodes.
    /// Gives real chi2 only if the contributions are uncorrelated.
    ChiSquare chi2Velo() const;

    /// Upstream segment chi square obtained by summing delta chi2 from all VPHit and UTHit fit nodes.
    /// Gives real chi2 only if the contributions are uncorrelated. Note that calculating chi2 of UT segment
    /// alone is not possible unless there are no Velo hits or no SciFi and Muon hits.
    ChiSquare chi2Upstream() const;

    /// Downstream segment chi square obtained by summing delta chi2 from all FTHit fit nodes.
    /// Gives real chi2 only if the contributions are uncorrelated.
    /// Note that the chi2 of SciFi segment is not correct if there are Muon hits.
    ChiSquare chi2Downstream() const;

    /// Upstream versus downstream segment chi square calculated as chi2Track - chi2Upstream - chi2Downstream.
    ChiSquare chi2Match() const;

    /// Muon segment chi square
    /// TODO PrKalmanFilter does not support muon tracks at the moment
    ChiSquare chi2Muon() const;
  }; // class TrackFitResult

  inline auto nodes( PrKalmanFitResult const& fr ) { return LHCb::span<const Pr::Tracks::Fit::Node>{fr.fitnodes}; }

} // namespace LHCb
