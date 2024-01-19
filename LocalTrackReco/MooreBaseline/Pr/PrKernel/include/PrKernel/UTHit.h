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
#include "Event/UTHit.h"
#include "Kernel/EventLocalAllocator.h"
#include "LHCbMath/bit_cast.h"

namespace UT {

  namespace Mut {
    // Mutable UTHit to allow a modifiable object, needed in some algorithms.
    // Maybe it is also usefull to use a template and merge with the ModPrHit from the seeding
    struct Hit {

      const UT::Hit* HitPtr;
      float          x, z;
      float          projection;

      Hit( const UT::Hit* ptr, float x, float z ) : HitPtr( ptr ), x( x ), z( z ) {}

      // Accessors to HitPtr
      [[nodiscard]] int  lhcbID() const { return bit_cast<int>( HitPtr->lhcbID().lhcbID() ); }
      [[nodiscard]] auto planeCode() const { return HitPtr->planeCode(); }
      [[nodiscard]] auto weight() const { return HitPtr->weight(); }
      [[nodiscard]] auto sin() const { return HitPtr->tanT() * HitPtr->cos(); }
      [[nodiscard]] auto zAtYEq0() const { return HitPtr->zAtYEq0(); }
      [[nodiscard]] bool isYCompatible( const float y, const float tol ) const {
        return HitPtr->isYCompatible( y, tol );
      }
      [[nodiscard]] auto xAt( const float y ) const { return HitPtr->xAt( y ); }
    };

    using Hits = std::vector<Hit, LHCb::Allocators::EventLocal<Hit>>;

    inline constexpr auto IncreaseByProj = []( const Hit& lhs, const Hit& rhs ) {
      if ( lhs.projection < rhs.projection ) return true;
      if ( rhs.projection < lhs.projection ) return false;
      return lhs.lhcbID() < rhs.lhcbID();
    };
  } // namespace Mut
} // namespace UT
