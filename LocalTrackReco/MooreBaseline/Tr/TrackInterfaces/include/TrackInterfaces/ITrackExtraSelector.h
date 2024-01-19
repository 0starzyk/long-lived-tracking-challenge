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
#ifndef TRACKINTERFACES_ITRACKEXTRASELECTOR_H
#define TRACKINTERFACES_ITRACKEXTRASELECTOR_H 1

// Include files
// -------------
// from Gaudi
#include "GaudiKernel/IAlgTool.h"

// Forward Declaration - from TrackInterfaces
struct ITrackExtrapolator;

/** @class ITrackExtraSelector ITrackExtraSelector.h TrackInterfaces/ITrackExtraSelector
 *
 *  Interface class to select which extrapolator to use.
 *
 */

struct ITrackExtraSelector : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( ITrackExtraSelector, 2, 0 );

  virtual const ITrackExtrapolator* select( double zStart, double zEnd ) const = 0;
};

#endif // TRACKINTERFACES_ITRACKEXTRASELECTOR_H
