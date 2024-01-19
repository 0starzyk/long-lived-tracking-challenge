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
#include "GaudiAlg/GaudiHistoAlg.h"
#include "Kernel/ITrajPoca.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include <map>
#include <string>

/** @class TrackMonitorBase TrackMonitorBase.h "TrackCheckers/TrackMonitorBase"
 *
 *  Base class for track monitoring: essentially a 'box' of common tools

 *  @author M. Needham.
 *  @date   6-5-2007
 */

class TrackMonitorBase : public GaudiHistoAlg {

public:
  /** Standard construtor */
  using GaudiHistoAlg::GaudiHistoAlg;

  /** Algorithm initialization */
  StatusCode initialize() override;

protected:
  /** Get a pointer to the track extrapolator
   *  @return extrapolator
   */
  const ITrackExtrapolator* extrapolator() const { return m_extrapolator.get(); }

private:
  ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator",
                                                "TrackMasterExtrapolator"}; ///< Pointer to extrapolator
};
