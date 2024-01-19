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
#include "GaudiAlg/GaudiTupleAlg.h"
#include "Kernel/ITrajPoca.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include <map>
#include <string>

/** @class TrackMonitorTupleBase TrackMonitorTupleBase.h "TrackCheckers/TrackMonitorTupleBase"
 *
 *  Base class for track monitoring: essentially a 'box' of common tools

 *  @author Ch. Elsasser
 *  @date   16-7-2013
 */

class TrackMonitorTupleBase : public GaudiTupleAlg {

public:
  /** Standard construtor */
  using GaudiTupleAlg::GaudiTupleAlg;

  /** Algorithm initialization */
  StatusCode initialize() override;

protected:
  /** Get a pointer to the poca tool
   *  @return poca tool
   */
  const ITrajPoca* pocaTool() const { return m_poca.get(); }

  /** Get a pointer to the track extrapolator
   *  @return extrapolator
   */
  const ITrackExtrapolator* extrapolator() const { return m_extrapolator.get(); }

  /** Whether to split by algorithm
   *  @return splitByAlgorithm true or false
   */
  bool splitByAlgorithm() const { return m_splitByAlgorithm; }
  bool splitByType() const { return m_splitByType; }

protected:
  void setSplitByType( bool b ) { m_splitByType = b; }

  std::string histoDirName( const LHCb::Track& track ) const;

private:
  Gaudi::Property<bool>          m_splitByAlgorithm{this, "SplitByAlgorithm", false};
  Gaudi::Property<bool>          m_splitByType{this, "SplitByType", true};
  ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator",
                                                "TrackMasterExtrapolator"}; ///< Pointer to extrapolator
  PublicToolHandle<ITrajPoca>    m_poca{this, "TrajPoca", "TrajPoca"};      ///< Pointer to the ITrajPoca interface
};
