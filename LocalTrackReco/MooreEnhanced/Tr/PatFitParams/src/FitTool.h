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
#ifndef FITTOOL_H
#define FITTOOL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"

#include "GaudiKernel/Point3DTypes.h"

#include "LinParFit.h"

static const InterfaceID IID_FitTool( "FitTool", 1, 0 );

/** @class FitTool FitTool.h
 *
 *
 *  @author Olivier Callot
 *  @date   2006-12-08
 */
class FitTool : public GaudiTool {
public:
  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_FitTool; }

  /// Standard constructor
  FitTool( const std::string& type, const std::string& name, const IInterface* parent );

  virtual ~FitTool(); ///< Destructor

  bool fitLine( const std::vector<Gaudi::XYZPoint>& hit, int mode, double z0, double& a, double& b );

  bool fitParabola( const std::vector<Gaudi::XYZPoint>& hit, int mode, double z0, double& a, double& b, double& c );

  bool fitCubic( const std::vector<Gaudi::XYZPoint>& hit, int mode, double x0, double& a, double& b, double& c,
                 double& d );

protected:
private:
  // avoid frequent reallocation of vectors etc.
  LinParFit<double> m_fit2;
  LinParFit<double> m_fit3;
  LinParFit<double> m_fit4;
};
#endif // FITTOOL_H
