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
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/Point3DTypes.h"
#include "PrFitParams/IPrFitTool.h"
#include "PrFitParams/LinParFit.h"

/** @class PrFitTool PrFitTool.h
 *
 *
 *  @author Olivier Callot
 *  @date   2006-12-08
 */
class PrFitTool final : public extends<GaudiTool, IPrFitTool> {
public:
  /// Standard constructor
  using extends::extends;

  std::optional<std::tuple<double, double>> fitLine( const std::vector<Gaudi::XYZPoint>& hit, XY mode,
                                                     double z0 ) const override;

  std::optional<std::tuple<double, double, double>> fitParabola( const std::vector<Gaudi::XYZPoint>& hit, XY mode,
                                                                 double z0 ) const override;

  std::optional<std::tuple<double, double, double, double>> fitCubic( const std::vector<Gaudi::XYZPoint>& hit, XY mode,
                                                                      double z0 ) const override;

private:
  mutable LinParFit<double> m_fit2{2};
  mutable LinParFit<double> m_fit3{3};
  mutable LinParFit<double> m_fit4{4};
};

DECLARE_COMPONENT( PrFitTool )

//=========================================================================
//  Fit a simple line in the specified projection
//=========================================================================
std::optional<std::tuple<double, double>> PrFitTool::fitLine( const std::vector<Gaudi::XYZPoint>& hits, XY mode,
                                                              double z0 ) const {
  m_fit2.clear();
  for ( const auto& p : hits ) m_fit2.accumulate( ( mode == XY::Y ? p.y() : p.x() ), 1., p.z() - z0 );
  if ( !m_fit2.solve() ) return {};
  return std::tuple{m_fit2[0], m_fit2[1]};
}
//=========================================================================
//  Fit a parabola in the specified projection
//=========================================================================
std::optional<std::tuple<double, double, double>> PrFitTool::fitParabola( const std::vector<Gaudi::XYZPoint>& hits,
                                                                          XY mode, double z0 ) const {
  m_fit3.clear();
  for ( const auto& p : hits ) m_fit3.accumulate( ( mode == XY::Y ? p.y() : p.x() ), 1., 1e-3 * ( p.z() - z0 ) );
  if ( !m_fit3.solve() ) return {};
  return std::tuple{m_fit3[0], 1e-3 * m_fit3[1], 1e-6 * m_fit3[2]};
}
//=========================================================================
//  Fit a Cubic in the specified projection
//=========================================================================
std::optional<std::tuple<double, double, double, double>> PrFitTool::fitCubic( const std::vector<Gaudi::XYZPoint>& hits,
                                                                               XY mode, double z0 ) const {
  m_fit4.clear();
  for ( const auto& p : hits ) m_fit4.accumulate( ( mode == XY::Y ? p.y() : p.x() ), 1., 1e-3 * ( p.z() - z0 ) );
  if ( !m_fit4.solve() ) return {};
  return std::tuple{m_fit4[0], 1e-3 * m_fit4[1], 1e-6 * m_fit4[2], 1e-9 * m_fit4[3]};
}

//=============================================================================
