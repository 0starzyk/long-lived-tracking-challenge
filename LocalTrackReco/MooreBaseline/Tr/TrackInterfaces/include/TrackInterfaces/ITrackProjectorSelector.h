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
#ifndef TRACKINTERFACES_ITRACKPROJECTORSELECTOR_H
#define TRACKINTERFACES_ITRACKPROJECTORSELECTOR_H 1

// Include files
// -------------
// from Gaudi
#include "GaudiKernel/IAlgTool.h"

// Forward Declaration - from TrackInterfaces
struct ITrackProjector;
namespace LHCb {
  class Measurement;
}

/** @class ITrackProjectorSelector ITrackProjectorSelector.h TrackInterfaces/ITrackProjectorSelector.h
 *
 *  Interface class to select which Projector to use.
 *
 */

class ITrackProjectorSelector : public extend_interfaces<IAlgTool> {
public:
  DeclareInterfaceID( ITrackProjectorSelector, 2, 0 );
  virtual ITrackProjector* projector( const LHCb::Measurement& m ) const = 0;
};

#endif // TRACKINTERFACES_ITRACKPROJECTORSELECTOR_H
