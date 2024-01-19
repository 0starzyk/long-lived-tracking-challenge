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
#pragma once

#include "Event/PrHits.h"
#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

#include "DetDesc/DetectorElement.h"

#include "Event/PartialChiSquareds.h"
#include "Event/PrKalmanFitResult.h"
#include "Event/Track_v3.h"

namespace LHCb::Pr {
  namespace Long {
    struct Tracks;
  }
  namespace Downstream {
    struct Tracks;
  }
  namespace Seeding {
    struct Tracks;
  }
  namespace Velo {
    struct Tracks;
  }
  namespace Upstream {
    struct Tracks;
  }
} // namespace LHCb::Pr

using V3FullOutput = std::tuple<LHCb::Event::v3::Tracks, LHCb::Event::v3::Track::PartialChiSquareds,
                                std::vector<LHCb::PrKalmanFitResult>>;

struct IPrFitterTool : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IPrFitterTool, 1, 1 );

  // Interface for fitting different supported track types (Velo, Downstream, Long, Upstream)
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Pr::Long::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // Fit PrDownstreamTracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Pr::Downstream::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // Fit PrSeedTracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Pr::Seeding::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // Fit PrVeloTracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Pr::Velo::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>& ) const = 0;

  // Fit PrUpstreamTracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Pr::Upstream::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>& ) const = 0;
  // fit v1 long tracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Event::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // fit v1 long tracks w/o UT hits
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Event::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // fit v1 longmuon tracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Event::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>& ) const = 0;
  // fit v1 longmuon tracks w/o UT hits
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Event::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>& ) const = 0;
  // fit v1 Ttrack tracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Event::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // fit v1 velo tracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Event::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>& ) const = 0;
  // fit v1 downstream tracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Event::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // fit v1 upstream tracks
  virtual std::tuple<LHCb::Event::Tracks> operator()( const LHCb::Event::Tracks&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                                      const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>& ) const = 0;

  // Same fit methods as above but now returning v3::Tracks
  virtual V3FullOutput fitted_v3_tracks( const LHCb::Pr::Long::Tracks&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;

  virtual V3FullOutput fitted_v3_tracks( const LHCb::Pr::Downstream::Tracks&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;

  virtual V3FullOutput fitted_v3_tracks( const LHCb::Pr::Seeding::Tracks&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;

  virtual V3FullOutput fitted_v3_tracks( const LHCb::Pr::Velo::Tracks&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>& ) const = 0;

  virtual V3FullOutput fitted_v3_tracks( const LHCb::Pr::Upstream::Tracks&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>& ) const = 0;

  // fit v1 long track
  virtual LHCb::Event::Track operator()( const LHCb::Event::Track&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // fit v1 longmuon track
  virtual LHCb::Event::Track operator()( const LHCb::Event::Track&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::Muon>& ) const = 0;
  // fit v1 Ttrack track
  virtual LHCb::Event::Track operator()( const LHCb::Event::Track&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // fit v1 velo track
  virtual LHCb::Event::Track operator()( const LHCb::Event::Track&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>& ) const = 0;
  // fit v1 downstream track
  virtual LHCb::Event::Track operator()( const LHCb::Event::Track&, const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::FT>& ) const = 0;
  // fit v1 upstream track
  virtual LHCb::Event::Track operator()( const LHCb::Event::Track&, const LHCb::Pr::Hits<LHCb::Pr::HitType::VP>&,
                                         const LHCb::Pr::Hits<LHCb::Pr::HitType::UT>& ) const = 0;
};