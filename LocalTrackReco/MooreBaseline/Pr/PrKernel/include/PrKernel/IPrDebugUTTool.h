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
#pragma once

#include "Event/Track.h"
#include "PrKernel/UTHit.h"

#include "GaudiKernel/IAlgTool.h"

class MsgStream;
#include "DetDesc/IGeometryInfo.h"

static const InterfaceID IID_IPrDebugUTTool( "IPrDebugUTTool", 1, 0 );

/** @class IPrDebugUTTool IPrDebugUTTool.h PrKernel/IPrDebugUTTool.h
 *
 *
 *  @author Olivier Callot
 *  @date   2007-10-22
 *
 *  @2017-03-01: Christoph Hasse (adapt to future framework)
 */
class IPrDebugUTTool : virtual public IAlgTool {
public:
  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_IPrDebugUTTool; }

  virtual void debugUTClusterOnTrack( const LHCb::Track* track, const UT::Mut::Hits::const_iterator beginCoord,
                                      const UT::Mut::Hits::const_iterator endCoord ) = 0;

  virtual void debugUTCluster( MsgStream& msg, const UT::Mut::Hit& hit ) = 0;

  virtual bool isTrueHit( const LHCb::Track* track, const UT::Mut::Hit& hit ) = 0;

  virtual double fracGoodHits( const LHCb::Track* track, const UT::Mut::Hits& hits ) = 0;

  virtual bool isTrueTrack( const LHCb::Track* track, const UT::Mut::Hits& hits ) = 0;

  virtual void chi2Tuple( const double p, const double chi2, const unsigned int nHits ) = 0;

  // added by AD 2/1/16 for efficiency vs step

  virtual void initializeSteps( std::vector<std::string> steps ) = 0; // initialize all steps in the process

  virtual void recordStepInProcess( std::string step, bool result ) = 0; // record the result of a step in the process

  virtual void resetflags() = 0; // reset all flags

  virtual void forceMCHits( UT::Mut::Hits& hits,
                            LHCb::Track*   track ) = 0; // Special. Force only MC matched hits in the track.

  virtual void tuneFisher( const LHCb::Track* seedTrack ) = 0;

  virtual void tuneDeltaP( const LHCb::Track* seedTrack, const double deltaP, const double momentum ) = 0;

  virtual void tuneFinalMVA( const LHCb::Track* seedTrack, const bool goodTrack, std::vector<double> vals ) = 0;

  virtual void getMagnetError( const LHCb::Track* seedTrack, IGeometryInfo const& geometry ) = 0;
};
