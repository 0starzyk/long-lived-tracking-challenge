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
#ifndef IPRFITTOOL_H
#define IPRFITTOOL_H 1

// Include files
#include <optional>

// from Gaudi
#include "GaudiKernel/IAlgTool.h"
#include "GaudiKernel/Point3DTypes.h"
#include "PrFitParams/LinParFit.h"

struct IPrFitTool : extend_interfaces<IAlgTool> {
  DeclareInterfaceID( IPrFitTool, 2, 0 );

  enum class XY { X, Y };

  virtual std::optional<std::tuple<double, double>> fitLine( const std::vector<Gaudi::XYZPoint>& hit, XY mode,
                                                             double z0 ) const = 0;

  virtual std::optional<std::tuple<double, double, double>> fitParabola( const std::vector<Gaudi::XYZPoint>& hit,
                                                                         XY mode, double z0 ) const = 0;

  virtual std::optional<std::tuple<double, double, double, double>> fitCubic( const std::vector<Gaudi::XYZPoint>& hit,
                                                                              XY mode, double z0 ) const = 0;
};
#endif // IPRFITTOOL_H
