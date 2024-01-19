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

#include "DetDesc/ITransportSvc.h"
#include "GaudiKernel/IAlgTool.h"
#include <any>

// Forward declarations
class Material;
namespace LHCb {
  class State;
}

/**
 *  Interface for state correction tools
 *
 *  @author Eduardo Rodrigues
 *  @date   2006-08-18
 */
struct IStateCorrectionTool : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IStateCorrectionTool, 2, 0 );
  /// Correct a State
  virtual void correctState( LHCb::State& state, const MaterialPtr material, std::any& cache, double wallThickness = 0,
                             bool upstream = true, double mass = 0 ) const = 0;

  virtual std::any createBuffer() const { return std::any(); }
};
