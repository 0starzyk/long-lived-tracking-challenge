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
#ifndef __TRACKKERNEL_TRACKTRAJVERTEX_H__
#define __TRACKKERNEL_TRACKTRAJVERTEX_H__

#include "Event/ZTrajectory.h"
#include "TrackKernel/TrackStateVertex.h"

namespace LHCb {

  class TrajVertex : public TrackStateVertex {
  public:
    TrajVertex( const std::vector<const ZTrajectory<double>*>& trajectories, double zseed,
                double ztolerance = 10, // mm
                double maxdchisq = 0.01, size_t maxiterations = 10 );
    ~TrajVertex();

    typedef std::vector<const ZTrajectory<double>*> Trajectories;
    const Trajectories&                             trajectories() const { return m_trajectories; }

    /// fit until converged
    FitStatus fit( double ztolerance = 10, double maxdchisq = 0.01, size_t maxiterations = 10 );

    /// adapative fit. downweight tracks with chi2 contribution larger than maxtrkchi2
    FitStatus fitAdaptive( double maxtrkchi2 = 4, double ztolerance = 10, double maxdchisq = 0.01,
                           size_t maxiterations = 10 );

  private:
    void updateStates( double z );

  private:
    Trajectories m_trajectories;
  };

} // namespace LHCb
#endif
