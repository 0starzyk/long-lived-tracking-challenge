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
#include "Event/State.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "TrackInterfaces/IPtTransporter.h" // Interface
#include <cmath>
#include <string>

//-----------------------------------------------------------------------------
// Implementation file for class : PtTransporter
//
// 2008-05-08 : Johannes Albrecht
//-----------------------------------------------------------------------------

/** @class PtTransporter PtTransporter.h
 *
 * provide a fast way to calculate a pt estimate in the velo, given a
 * state after the magnet
 *
 * @author Manuel Tobias Schiller <schiller@phys.uni-heidelberg.de>
 * @date   2008-04-16
 */
class PtTransporter : public extends<GaudiTool, IPtTransporter> {
public:
  /// Standard constructor
  using extends::extends;

  double ptAtOrigin( double zref, double xref, double yref, double tx, double ty, double p ) const override;
  double ptAtOrigin( const LHCb::State& state ) const override;

private:
  Gaudi::Property<double> m_zMagnet{this, "zMagnet", 5300. * Gaudi::Units::mm};
};
// Declaration of the Tool Factory
DECLARE_COMPONENT( PtTransporter )

//=============================================================================

double PtTransporter::ptAtOrigin( double zref, double xref, double /* yref */, double tx, double ty, double p ) const {
  // assume B field conserves magnitude of p, and assume that py is
  // not altered; model effect of magnetic field using kick in center
  // of magnet plane
  const double xmag = xref + tx * ( m_zMagnet - zref );
  const double r    = std::sqrt( xmag * xmag + m_zMagnet * m_zMagnet );
  double       py   = p * ty / std::sqrt( 1. + tx * tx + ty * ty );
  p *= xmag / r;
  p *= p;
  py *= m_zMagnet / r;
  py *= py;
  return std::sqrt( p + py );
}

double PtTransporter::ptAtOrigin( const LHCb::State& state ) const {
  // protect against division by zero
  if ( std::abs( state.qOverP() ) < 1e-42 ) return HUGE_VAL;
  return ptAtOrigin( state.z(), state.x(), state.y(), state.tx(), state.ty(), 1. / std::abs( state.qOverP() ) );
}
