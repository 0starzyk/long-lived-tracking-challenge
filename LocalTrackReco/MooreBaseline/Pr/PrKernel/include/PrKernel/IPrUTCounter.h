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
#ifndef PATKERNEL_IPAUTTCOUNTER_H
#define PATKERNEL_IPAUTTCOUNTER_H 1

// Include files
// from STL
#include <string>

// from Gaudi
#include "GaudiKernel/IAlgTool.h"

#include "Kernel/LHCbID.h"

namespace LHCb {
  class MCParticle;
} // namespace LHCb

static const InterfaceID IID_IPrUTCounter( "IPrUTCounter", 1, 1 );

/** @class IPrUTCounter IPrUTCounter.h PrKernel/IPrUTCounter.h
 *
 *
 *  @author Wenbin Qian
 *  @date   2011-03-21
 */
class IPrUTCounter : virtual public IAlgTool {
public:
  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_IPrUTCounter; }

  virtual void initEvent() = 0;

  virtual void count( const LHCb::MCParticle* part, std::vector<bool> flags, std::vector<LHCb::LHCbID>& ids ) = 0;

  virtual void setContainer( std::string name ) = 0;

  virtual void addSelection( std::string name ) = 0;

  virtual void printStatistics() = 0;
};
#endif // PATKERNEL_IPAUTTCOUNTER_H
