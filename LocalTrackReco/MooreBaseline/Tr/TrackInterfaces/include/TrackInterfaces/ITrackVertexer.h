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

#include "Event/Track.h"
#include "Kernel/STLExtensions.h"

#include "GaudiKernel/IAlgTool.h"

#include <memory>
#include <vector>

#include "DetDesc/IGeometryInfo.h"
namespace LHCb {
  class TwoProngVertex;
  class State;
  class RecVertex;
} // namespace LHCb

/**
 *  @author Wouter HULSBERGEN
 *  @date   2007-11-07
 *
 */
struct ITrackVertexer : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackVertexer, 4, 0 );

  virtual std::unique_ptr<LHCb::TwoProngVertex> fit( const LHCb::State& stateA, const LHCb::State& stateB,
                                                     IGeometryInfo const& geometry ) const = 0;
  virtual std::unique_ptr<LHCb::RecVertex>      fit( LHCb::span<const LHCb::State* const> states,
                                                     IGeometryInfo const&                 geometry ) const = 0;
  virtual std::unique_ptr<LHCb::RecVertex>      fit( LHCb::span<LHCb::Track const* const> tracks,
                                                     IGeometryInfo const&                 geometry ) const = 0;
  virtual bool computeDecayLength( const LHCb::TwoProngVertex& vertex, const LHCb::RecVertex& pv, double& chi2,
                                   double& decaylength, double& decaylengtherr ) const     = 0;

#if defined( __GNUC__ ) && ( __GNUC__ < 10 )
  // cppgsl 3.x gcc 9.x workaround. To be removed when gcc 9 or older no longer supported
  auto fit( std::vector<LHCb::State*> const& states, IGeometryInfo const& geometry ) const {
    using StateSpan = LHCb::span<const LHCb::State* const>;
    return fit( StateSpan{states.data(), static_cast<StateSpan::size_type>( states.size() )}, geometry );
  }
#endif

  /// Return the ip chi2 for a track (uses stateprovider, not good for
  /// HLT: better call routine below with track->firstState())
  virtual double ipchi2( const LHCb::Track& track, const LHCb::RecVertex& pv, IGeometryInfo const& geometry ) const = 0;

  /// Return the ip chi2 for a track state
  virtual double ipchi2( const LHCb::State& state, const LHCb::RecVertex& pv, IGeometryInfo const& geometry ) const = 0;
};
