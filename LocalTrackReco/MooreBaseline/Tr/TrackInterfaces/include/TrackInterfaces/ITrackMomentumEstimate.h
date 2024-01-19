/*****************************************************************************\
* (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once

// Include files
// -------------
// from Gaudi
#include "Event/SOACollection.h"
#include "GaudiKernel/IAlgTool.h"
#include "LHCbMath/SIMDWrapper.h"
#include <Magnet/DeMagnet.h>
// forward declarations
namespace LHCb {
  class State;
}

/** @class ITrackMomentumEstimate ITrackMomentumEstimate.h TrackInterfaces/ITrackMomentumEstimate.h
 *
 *
 *  @author Stephanie Hansmann-Menzemer
 *  @date   2007-10-30
 */
struct ITrackMomentumEstimate : extend_interfaces<IAlgTool> {
  using simd = SIMDWrapper::best::types;

  // Return the interface ID
  DeclareInterfaceID( ITrackMomentumEstimate, 2, 0 );

  // Estimate the momentum P of a State in T at ZAtMidT
  virtual StatusCode calculate( const DeMagnet& magnet, const LHCb::State* TState, double& qOverP, double& sigmaQOverP,
                                bool cubical = 0 ) const = 0;

  // Estimate the momentum P of a velo State and a State in T at ZAtMidT
  virtual StatusCode calculate( const DeMagnet& magnet, const LHCb::State* veloState, const LHCb::State* tState,
                                double& qOverP, double& sigmaQOverP, bool cubical = 0 ) const = 0;

  // Estimate the momentum P of a velo State and a State in T at ZAtMidT
  virtual StatusCode calculate( const DeMagnet& magnet, const simd::float_v txT, const simd::float_v txV,
                                const simd::float_v tyV, simd::float_v& qOverP, simd::float_v& sigmaQOverP,
                                bool cubical = 0 ) const = 0;

  // Estimate the momentum P of a velo State and a State in T at ZAtMidT
  virtual StatusCode calculate( const DeMagnet& magnet, const simd::float_v txT, const simd::float_v txY,
                                const simd::float_v xT, const simd::float_v zT, simd::float_v& qOverP,
                                simd::float_v& sigmaQOverP, bool cubical = 0 ) const = 0;
};
