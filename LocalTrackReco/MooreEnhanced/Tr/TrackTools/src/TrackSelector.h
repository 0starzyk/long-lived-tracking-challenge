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

//-----------------------------------------------------------------------------
/** @file TrackSelector.h
 *
 *  Header file for reconstruction tool : TrackSelector
 *
 *  @author M.Needham Matt.Needham@cern.ch
 *  @author Chris Jones   Christopher.Rob.Jones@cern.ch
 *  @date   30/12/2005
 */
//-----------------------------------------------------------------------------

#ifndef TRACKTOOLS_TrackSelector_H
#define TRACKTOOLS_TrackSelector_H

//-----------------------------------------------------------------------------
/** @class TrackSelector TrackSelector.h
 *
 *  General track Selection tool
 *
 *  Cuts can be applied on various quantities like p, hits, chi^2, pt, and track type.
 *
 *  @author M.Needham Matt.Needham@cern.ch
 *  @author C. Jones  Christopher.Rob.Jones@cern.ch
 *
 *  @date   30/12/2005
 */
//-----------------------------------------------------------------------------

// STL
#include <sstream>
#include <string>

// base class
#include "TrackSelectorBase.h"

// boost
#include <limits>

class TrackSelector : public TrackSelectorBase {

public:
  /// constructer
  using TrackSelectorBase::TrackSelectorBase;

  /** Returns if the given track is selected or not
   *
   *  @param aTrack Reference to the Track to test
   *
   *  @return boolean indicating if the track is selected or not
   *  @retval true  Track is selected
   *  @retval false Track is rejected
   */
  bool accept( const LHCb::Track& aTrack ) const override;

private:
  Gaudi::Property<double> m_minChi2Cut{this, "MinChi2Cut", -1, "Min chi2 cut"};
  Gaudi::Property<double> m_maxChi2Cut{this, "MaxChi2Cut", -1, "Max chi2 cut"};

  Gaudi::Property<double> m_minPCut{this, "MinPCut", 0.0, "Min p cut in GeV"};
  Gaudi::Property<double> m_maxPCut{this, "MaxPCut", -1, "Max p cut in GeV"};

  Gaudi::Property<double> m_minPtCut{this, "MinPtCut", 0.0, "Min p cut in GeV"};
  Gaudi::Property<double> m_maxPtCut{this, "MaxPtCut", -1, "Max p cut in GeV"};

  Gaudi::Property<int> m_minNDoF{this, "MinNDoF", 0, "Minimum number of dofs"};
  Gaudi::Property<int> m_maxNDoF{this, "MaxNDoF", std::numeric_limits<int>::max(), "Maximum number of dofs"};

  Gaudi::Property<double> m_minEtaCut{this, "MinEtaCut", std::numeric_limits<double>::lowest(),
                                      "Minimum track eta cut"};
  Gaudi::Property<double> m_maxEtaCut{this, "MaxEtaCut", std::numeric_limits<double>::max(), "Maximum track eta cut"};

  Gaudi::Property<double> m_minPhiCut{this, "MinPhiCut", std::numeric_limits<double>::lowest(),
                                      "Minimum track phi cut"};
  Gaudi::Property<double> m_maxPhiCut{this, "MaxPhiCut", std::numeric_limits<double>::max(), "Maximum track phi cut"};

  Gaudi::Property<double> m_minLikCut{this, "MinLikelihoodCut", std::numeric_limits<double>::lowest()};
  Gaudi::Property<double> m_maxLikCut{this, "MaxLikelihoodCut", std::numeric_limits<double>::max()};

  Gaudi::Property<bool>   m_acceptClones{this, "AcceptClones", true, "Flag to turn on/off reject of clones"};
  Gaudi::Property<double> m_minCloneCut{this, "MinCloneDistCut", -1e10, "Minimum Clone distance cut"};
  Gaudi::Property<double> m_maxCloneCut{this, "MaxCloneDistCut", std::numeric_limits<double>::max(),
                                        "Maximum Clone distance cut"};

  Gaudi::Property<double> m_minGhostProb{this, "MinGhostProbCut", std::numeric_limits<double>::lowest(),
                                         "minimum ghost probability cut"};
  Gaudi::Property<double> m_maxGhostProb{this, "MaxGhostProbCut", std::numeric_limits<double>::max(),
                                         "maximum ghost probability cut"};

  Gaudi::Property<double> m_maxChi2Velo{this, "MaxChi2PerDoFVelo", -1};
  Gaudi::Property<double> m_maxChi2Upstream{this, "MaxChi2PerDoFUpstream", -1};
  Gaudi::Property<double> m_maxChi2Downstream{this, "MaxChi2PerDoFDownstream", -1};
  Gaudi::Property<double> m_maxChi2Match{this, "MaxChi2PerDoFMatch", -1};

  Gaudi::Property<int> m_minNVeloHits{this, "MinNVeloHits", 0};
  Gaudi::Property<int> m_minNTHits{this, "MinNTHits", 0};
  Gaudi::Property<int> m_minNVeloLayers{this, "MinNVeloLayers", 0};
  Gaudi::Property<int> m_minNVeloALayers{this, "MinNVeloALayers", 0};
  Gaudi::Property<int> m_minNVeloCLayers{this, "MinNVeloCLayers", 0};
  Gaudi::Property<int> m_minNVeloOverlap{this, "MinNVeloOverlap", 0};
  Gaudi::Property<int> m_minNTLayers{this, "MinNTLayers", 0};
  Gaudi::Property<int> m_minNUTLayers{this, "MinNUTLayers", 0};
  Gaudi::Property<int> m_maxNVeloHoles{this, "MaxNVeloHoles", 999, "Maximum number of holes in Velo segment"};
  Gaudi::Property<int> m_maxNTHoles{this, "MaxNTHoles", 999, "Maximum number of holes in T segment"};
};

#endif // TRACKTOOLS_TrackSelector_H
