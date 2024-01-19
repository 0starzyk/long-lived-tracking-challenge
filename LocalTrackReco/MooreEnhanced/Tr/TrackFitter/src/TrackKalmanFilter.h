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
#ifndef TRACKFITTER_TRACKKALMANFILTER_H
#define TRACKFITTER_TRACKKALMANFILTER_H 1

// Include files
// -------------
// from Gaudi
#include "GaudiAlg/GaudiTool.h"

// from TrackInterfaces
#include "TrackInterfaces/ITrackKalmanFilter.h"

// from TrackEvent
#include "Event/Track.h"

// From LoKi
#include "GaudiKernel/Range.h"

/** @class TrackKalmanFilter TrackKalmanFilter.h
 *
 *
 *  @author Jose Angel Hernando Morata, Eduardo Rodrigues
 *  @date   2005-04-15
 *  reusing the previous code
 *  @author Rutger van der Eijk  07-04-1999
 *  @author Mattiew Needham
 */

class TrackKalmanFilter : public extends<GaudiTool, ITrackKalmanFilter> {
public:
  /// Standard constructor
  using extends::extends;

  //! fit a track
  StatusCode fit( LHCb::Track& track ) const override;

private:
  void printErrMeasures( LHCb::Track& track ) const;

  void printStates( LHCb::Track& track ) const;

  // job options
  Gaudi::Property<bool> m_forceBiDirectionalFit{this, "ForceBiDirectionalFit",
                                                true};             ///< Flag for forcing bidirectional fit
  Gaudi::Property<bool> m_forceSmooth{this, "ForceSmooth", false}; ///< Flag for force the smoothing (for debug reason)
  Gaudi::Property<unsigned int> m_DoF{this, "DoF", 5u};

  //! helper to print a failure comment
  StatusCode failure( const std::string& comment ) const;
};
#endif // TRACKFITTER_TRACKKALMANFILTER_H
