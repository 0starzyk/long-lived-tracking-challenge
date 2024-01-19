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
#include "GaudiKernel/SystemOfUnits.h"
#include "TrackProjector.h"

/** @class TrajProjector TrajProjector.h TrajProjector.h
 *
 *  Projects into the measurement space
 *
 *  @author Jeroen van Tilburg
 *  @date   2006-02-21
 */

template <typename T>
struct TrajProjector : TrackProjector {

  /// FIXME/TODO: project given a traj instead of a state

  /// Standard constructor
  TrajProjector( const std::string& type, const std::string& name, const IInterface* parent );
};

//-----------------------------------------------------------------------------
/// Standard constructor, initializes variables
//-----------------------------------------------------------------------------
template <typename T>
TrajProjector<T>::TrajProjector( const std::string& type, const std::string& name, const IInterface* parent )
    : TrackProjector( type, name, parent ) {
  setProperty( "Tolerance", T::defaultTolerance() ).ignore();
}

namespace TrajProj {
  /// declare and instantiate ST and Velo projectors...
  struct ST {
    static constexpr double defaultTolerance() { return 0.002 * Gaudi::Units::mm; }
  };
  struct UT {
    static constexpr double defaultTolerance() { return 0.002 * Gaudi::Units::mm; }
  };
  struct Velo {
    static constexpr double defaultTolerance() { return 0.0005 * Gaudi::Units::mm; }
  };
  struct VP {
    static constexpr double defaultTolerance() { return 0.0005 * Gaudi::Units::mm; }
  };
  struct FT {
    static constexpr double defaultTolerance() { return 0.002 * Gaudi::Units::mm; }
  };
  struct Muon {
    static constexpr double defaultTolerance() { return 0.002 * Gaudi::Units::mm; }
  };
} // namespace TrajProj

typedef TrajProjector<TrajProj::ST> TrajSTProjector;
DECLARE_COMPONENT( TrajSTProjector )

typedef TrajProjector<TrajProj::UT> TrajUTProjector;
DECLARE_COMPONENT( TrajUTProjector )

typedef TrajProjector<TrajProj::Velo> TrajVeloProjector;
DECLARE_COMPONENT( TrajVeloProjector )

typedef TrajProjector<TrajProj::VP> TrajVPProjector;
DECLARE_COMPONENT( TrajVPProjector )

typedef TrajProjector<TrajProj::FT> TrajFTProjector;
DECLARE_COMPONENT( TrajFTProjector )

typedef TrajProjector<TrajProj::Muon> TrajMuonProjector;
DECLARE_COMPONENT( TrajMuonProjector )
