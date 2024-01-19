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

#include <vector>

// Gaudi
#include "LHCbAlgs/Transformer.h"

// LHCb
#include "Event/PrVeloTracks.h"
#include "Event/PrimaryVertices.h"
#include "Event/RecVertex_v2.h"
#include "Kernel/AllocatorUtils.h"

namespace {
  using Tracks = LHCb::Pr::Velo::Tracks;
  // using Vertices = LHCb::Event::v2::RecVertices;
} // namespace

/**
 * Vectorized IPFilter for TracksVP
 *
 * @author Arthur Hennequin (CERN, LIP6)
 */
template <typename Vertices>
struct FilterIP : public LHCb::Algorithm::Transformer<Tracks( Tracks const&, Vertices const& )> {
  using base_class = LHCb::Algorithm::Transformer<Tracks( Tracks const&, Vertices const& )>;
  FilterIP( const std::string& name, ISvcLocator* pSvcLocator )
      : base_class( name, pSvcLocator,
                    {{"Input", "Rec/Track/Velo"}, {"InputVertices", LHCb::Event::v2::RecVertexLocation::Primary}},
                    {"Output", ""} ) {}

  Tracks operator()( const Tracks& tracks, const Vertices& vertices ) const override {
    using dType = SIMDWrapper::best::types;
    using F     = dType::float_v;

    const F ip_cut_value = m_ipcut.value();

    auto tracks_out =
        LHCb::make_obj_propagating_allocator<LHCb::Pr::Velo::Tracks>( tracks, Zipping::generateZipIdentifier() );
    tracks_out.reserve( tracks.size() );

    for ( auto const& track : tracks.simd() ) {
      auto loop_mask = track.loop_mask(); // true if track.offset()<tracks.size() else false

      Vec3<F> B = track.StatePos( 0 ); // get the origin and direction
      Vec3<F> u = track.StateDir( 0 ); // of the state

      F min_d = 10e3;
      for ( int j = 0; j < (int)vertices.size(); j++ ) {
        auto    PV = vertices[j].position();
        Vec3<F> A  = Vec3<F>( PV.x(), PV.y(), PV.z() );
        auto    d  = ( B - A ).cross( u ).mag2();
        min_d      = min( min_d, d );
      }

      auto d    = sqrt( min_d ) / u.mag(); // distance from closest PV to line
      auto mask = ip_cut_value < d;        // ip cut

      auto const i = track.offset();
      tracks_out.copy_back<dType>( tracks, i, mask && loop_mask ); // conditional push_back to the output
    }

    m_nbTracksCounter += tracks_out.size();
    return tracks_out;
  };

private:
  mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
  Gaudi::Property<float>                        m_ipcut{this, "IPcut", 0.0};
};

DECLARE_COMPONENT_WITH_ID( FilterIP<LHCb::Event::v2::RecVertices>, "PrFilterIPSoAV2" )
DECLARE_COMPONENT_WITH_ID( FilterIP<LHCb::Event::PV::PrimaryVertexContainer>, "PrFilterIPSoA" )
