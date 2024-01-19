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
#pragma once

#include "Event/SOACollection.h"
#include "Kernel/EventLocalAllocator.h"
#include "LHCbMath/SIMDWrapper.h"
#include "UTDAQ/UTInfo.h"

/** Mutable UT hit class for internal use in pattern recognition algorithms
 *
 *  @author Michel De Cian
 *  @date   2020-04-06
 */

namespace LHCb::Pr::UT {

  namespace Mut {

    namespace HitTag {
      struct xs : Event::float_field {};
      struct zs : Event::float_field {};
      struct coss : Event::float_field {};
      struct sins : Event::float_field {};
      struct weights : Event::float_field {};
      struct projections : Event::float_field {};
      struct channelIDs : Event::int_field {};
      struct indexs : Event::int_field {};

      template <typename T>
      using muthit_t = Event::SOACollection<T, xs, zs, coss, sins, weights, projections, channelIDs, indexs>;
    } // namespace HitTag

    struct Hits : HitTag::muthit_t<Hits> {
      using base_t = typename HitTag::muthit_t<Hits>;
      using base_t::base_t;

      std::array<int, static_cast<int>( UTInfo::DetectorNumbers::TotalLayers )> layerIndices;

      template <SIMDWrapper::InstructionSet simd, LHCb::Pr::ProxyBehaviour behaviour, typename ContainerType>
      struct MutHitProxy : Event::Proxy<simd, behaviour, ContainerType> {
        using Event::Proxy<simd, behaviour, ContainerType>::Proxy;
        [[nodiscard]] auto x() const { return this->template get<HitTag::xs>(); }
        [[nodiscard]] auto z() const { return this->template get<HitTag::zs>(); }
        [[nodiscard]] auto cos() const { return this->template get<HitTag::coss>(); }
        [[nodiscard]] auto sin() const { return this->template get<HitTag::sins>(); }
        [[nodiscard]] auto weight() const { return this->template get<HitTag::weights>(); }
        [[nodiscard]] auto projection() const { return this->template get<HitTag::projections>(); }
        [[nodiscard]] auto channelID() const { return this->template get<HitTag::channelIDs>(); }
        [[nodiscard]] auto index() const { return this->template get<HitTag::indexs>(); }

        /// Retrieve the plane code
        auto planeCode() const {
          auto id      = channelID();
          auto station = ( id & static_cast<int>( UTInfo::MasksBits::StationMask ) ) >>
                         static_cast<int>( UTInfo::MasksBits::StationBits );
          auto layer = ( id & static_cast<int>( UTInfo::MasksBits::LayerMask ) ) >>
                       static_cast<int>( UTInfo::MasksBits::LayerBits );

          return 2 * ( station - 1 ) + ( layer - 1 );
        }
      };
      template <SIMDWrapper::InstructionSet simd, LHCb::Pr::ProxyBehaviour behaviour, typename ContainerType>
      using proxy_type = MutHitProxy<simd, behaviour, ContainerType>;
    };

  } // namespace Mut
} // namespace LHCb::Pr::UT
