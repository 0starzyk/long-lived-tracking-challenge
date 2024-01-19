/*****************************************************************************\
* (c) Copyright 2022 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#pragma once

#include "PrKernel/IPrDebugTrackingTool.h"

#include "GaudiAlg/GaudiTupleTool.h"
#include <string_view>

struct PrDebugTrackingToolBase : public extends<GaudiTupleTool, IPrDebugTrackingTool> {

  // inherit standard constructors
  using extends::extends;

  /**
   * @brief This is the default implementation for storing data in the tool interface.
   *
   * @param vars_and_values
   * @param tuple_name
   *
   * @note Supports tupling of int, float, std::vector<int>, std::vector<float>. The vectors
   * have an arbitrary size limit of 1024 entries.
   */

  virtual void storeData( LHCb::span<const VariableDef> vars_and_values, std::string_view tuple_name ) const override {

    Tuple tuple = nTuple( std::string{tuple_name} );
    for ( auto [name, value] : vars_and_values ) {
      if ( std::holds_alternative<int>( value ) ) {
        tuple->column( name, std::get<int>( value ) ).ignore();
      } else if ( std::holds_alternative<float>( value ) ) {
        tuple->column( name, std::get<float>( value ) ).ignore();
      } else if ( std::holds_alternative<double>( value ) ) {
        tuple->column( name, std::get<double>( value ) ).ignore();
      } else if ( std::holds_alternative<std::vector<int>>( value ) ) {
        tuple->farray( name, std::get<std::vector<int>>( value ), std::string{name} + "_length", 1024 ).ignore();
      } else if ( std::holds_alternative<std::vector<float>>( value ) ) {
        tuple->farray( name, std::get<std::vector<float>>( value ), std::string{name} + "_length", 1024 ).ignore();
      } else {
        tuple->farray( name, std::get<std::vector<double>>( value ), std::string{name} + "_length", 1024 ).ignore();
      }
    }
    tuple->write().ignore();
  }
};
