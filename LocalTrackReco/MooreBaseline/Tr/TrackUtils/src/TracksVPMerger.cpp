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

#include "Event/PrVeloTracks.h"
#include "LHCbAlgs/Transformer.h"
#include <vector>

/**
 * Merge two TracksVP containers into one.
 *
 * @author Arthur Hennequin (CERN, LIP6)
 */

namespace LHCb {

  class TracksVPMerger
      : public Algorithm::Transformer<Pr::Velo::Tracks( const Pr::Velo::Tracks&, const Pr::Velo::Tracks& )> {
    using Tracks = Pr::Velo::Tracks;

  public:
    TracksVPMerger( const std::string& name, ISvcLocator* pSvcLocator )
        : Transformer(
              name, pSvcLocator,
              {KeyValue{"TracksLocation1", "Rec/Track/VeloBackward"}, KeyValue{"TracksLocation2", "Rec/Track/Velo"}},
              KeyValue{"OutputTracksLocation", "Rec/Track/VeloMerged"} ) {}

    Tracks operator()( const Tracks& tracks1, const Tracks& tracks2 ) const override {
      Tracks     out;
      const auto N = tracks1.size() + tracks2.size();
      out.reserve( N );
      m_nbTracksCounter += N;

      using dType = SIMDWrapper::best::types;

      for ( auto const& track : tracks1.simd() ) {
        auto       loop_mask = track.loop_mask();
        auto const t         = track.offset();
        out.copy_back<dType>( tracks1, t, loop_mask );
      }

      for ( auto const& track : tracks2.simd() ) {
        auto       loop_mask = track.loop_mask();
        auto const t         = track.offset();
        out.copy_back<dType>( tracks2, t, loop_mask );
      }

      return out;
    };

  private:
    mutable Gaudi::Accumulators::SummingCounter<> m_nbTracksCounter{this, "Nb of Produced Tracks"};
  };

  DECLARE_COMPONENT_WITH_ID( TracksVPMerger, "TracksVPMerger" )

} // namespace LHCb
