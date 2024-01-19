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
#ifndef TRACKINTERFACES_IPTTRANSPORTER_H
#define TRACKINTERFACES_IPTTRANSPORTER_H 1

// Include files

// from Gaudi
#include "GaudiKernel/IAlgTool.h"

namespace LHCb {
  class State;
}

/** @class IPtTransporter IPtTransporter.h TrackInterfaces/IPtTransporter.h
 *
 *  calculate pt at origin from a given state at T
 *
 *  @author Johannes Albrecht
 *  @date   2008-05-08
 */
struct IPtTransporter : extend_interfaces<IAlgTool> {

  // Return the interface ID
  DeclareInterfaceID( IPtTransporter, 2, 0 );

  virtual double ptAtOrigin( double zref, double xref, double yref, double tx, double ty, double p ) const = 0;

  virtual double ptAtOrigin( const LHCb::State& stateAtT ) const = 0;
};
#endif // TRACKINTERFACES_IPTTRANSPORTER_H
