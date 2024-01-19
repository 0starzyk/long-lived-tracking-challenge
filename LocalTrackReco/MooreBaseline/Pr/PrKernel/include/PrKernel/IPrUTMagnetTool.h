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
#ifndef IPRUTMAGNETTOOL_H
#define IPRUTMAGNETTOOL_H 1

// Include files
#include "Kernel/LUTForFunction.h"
#include <GaudiKernel/IAlgTool.h>
#include <GaudiKernel/SystemOfUnits.h>

/** @class IPrDebugTool IPrUTMagnetTool.h PrKernel/IPrUTMagnetTool.h
 *  Interface to the tool used by the upstream pattern recognition to
 *  obtain integrals of the B-field as a function of parameters of a
 *  Velo track
 *
 *  @author Roel Aaij
 *  @date   2019-05-29
 */
class IPrUTMagnetTool : public extend_interfaces<IAlgTool> {
public:
  struct Cache {
    Cache()
        : lutBdl{std::make_unique<LUTForFunction<30, 10, 10>>(
              LUT::Range{-0.3, +0.3}, LUT::Range{-250. * Gaudi::Units::mm, 250. * Gaudi::Units::mm},
              LUT::Range{0. * Gaudi::Units::mm, 800. * Gaudi::Units::mm} )}
        , lutZHalfBdl{std::make_unique<LUTForFunction<30, 10, 10>>(
              LUT::Range{-0.3, +0.3}, LUT::Range{-250. * Gaudi::Units::mm, 250. * Gaudi::Units::mm},
              LUT::Range{0. * Gaudi::Units::mm, 800. * Gaudi::Units::mm} )}
        , lutDxLay{std::make_unique<LUTForFunction<3, 30>>( LUT::Range{0., 3.}, LUT::Range{0.0, 0.3} )}
        , lutDxToMom{std::make_unique<LUTForFunction<30>>( LUT::Range{0.0, 0.3} )} {}

    std::array<float, 4> zLayers = {0};

    float zCenterUT = std::numeric_limits<float>::signaling_NaN();
    float zMidField = std::numeric_limits<float>::signaling_NaN();
    float dist2mom  = std::numeric_limits<float>::signaling_NaN();

    std::unique_ptr<LUTForFunction<30, 10, 10>> lutBdl;
    std::unique_ptr<LUTForFunction<30, 10, 10>> lutZHalfBdl;
    std::unique_ptr<LUTForFunction<3, 30>>      lutDxLay;
    std::unique_ptr<LUTForFunction<30>>         lutDxToMom;

    bool noField = false;
  };

  DeclareInterfaceID( IPrUTMagnetTool, 1, 0 );

  virtual float bdlIntegral( float ySlopeVelo, float zOrigin, float zVelo ) const = 0;
  virtual float zBdlMiddle( float ySlopeVelo, float zOrigin, float zVelo ) const  = 0;
  virtual float dist2mom( float ySlope ) const                                    = 0;

  virtual float zMidUT() const    = 0;
  virtual float zMidField() const = 0;

  virtual float averageDist2mom() const = 0;

  virtual LUTForFunction<3, 30> const&      DxLayTable() const = 0;
  virtual LUTForFunction<30, 10, 10> const& BdlTable() const   = 0;

  virtual const std::string& cacheLocation() const = 0;
};

#endif
