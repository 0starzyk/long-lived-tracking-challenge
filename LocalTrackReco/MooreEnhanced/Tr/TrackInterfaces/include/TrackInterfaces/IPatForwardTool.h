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
#ifndef TRACKINTERFACES_IPATFORWARDTOOL_H
#define TRACKINTERFACES_IPATFORWARDTOOL_H 1

// Include files
// from Gaudi
#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

/** @class IPatForwardTool IPatForwardTool.h TrackInterfaces/IPatForwardTool.h
 *  Interface to the forward pattern tool
 *
 *  @author Olivier Callot
 *  @date   2005-10-04
 */
class IPatForwardTool : public extend_interfaces<IAlgTool> {
public:
  DeclareInterfaceID( IPatForwardTool, 2, 0 );

  virtual void forwardTrack( const LHCb::Track& track, LHCb::Tracks& output ) const = 0;
  virtual void setNNSwitch( bool nnSwitch )                                         = 0;
};
#endif // TRACKINTERFACES_IPATFORWARDTOOL_H
