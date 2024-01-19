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
#ifndef PRKERNEL_IPRDEBUGTOOL_H
#define PRKERNEL_IPRDEBUGTOOL_H 1

// Include files
// from STL
#include <string>

// from Gaudi
#include "GaudiKernel/IAlgTool.h"
#include "Kernel/LHCbID.h"

static const InterfaceID IID_IPrDebugTool( "IPrDebugTool", 1, 0 );

/** @class IPrDebugTool IPrDebugTool.h PrKernel/IPrDebugTool.h
 *  Interface to the Pattern debug tool for the upgrade
 *
 *  @author Olivier Callot
 *  @date   2012-03-22
 */
class IPrDebugTool : virtual public IAlgTool {
public:
  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_IPrDebugTool; }

  virtual bool matchKey( LHCb::LHCbID id, int key ) const = 0;

  virtual void printKey( MsgStream& msg, LHCb::LHCbID id ) const = 0;
};
#endif // PRKERNEL_IPRDEBUGTOOL_H
