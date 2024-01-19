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
//============================================================================
/** @file Track1DTabFunc.h
 *
 *  Header file for utility class : Track::TabulatedFunction1D
 *
 *  @author Paul Seyfert      Paul.Seyfert@cern.ch
 *  @date   2016-12-20
 */
//============================================================================

#ifndef TRACKKERNEL_TRACK1DTABFUNC_H
#define TRACKKERNEL_TRACK1DTABFUNC_H 1

#include <vector>

namespace Track {

  //============================================================================
  /** @class Track::TabulatedFunction1D Track1DTabFunc.h
   *
   *  A class describing a tabulated function with equidistant y-binning from 0 to 1.
   *
   *  @author Paul Seyfert      Paul.Seyfert@cern.ch
   *  @date   2016-12-20
   */
  //============================================================================

  class TabulatedFunction1D final {

  public:
    /** Default Constructor with optional interpolator type argument
     *
     *  @param x         braced initializer list for the x bins (N bins = N+1 x values)
     */
    TabulatedFunction1D( std::initializer_list<float> x );

    /** Computes the function value (y) for the given parameter (x) value
     *  with linear interpolation between bin edges
     *  does out of range check
     *
     *  @param x The parameter value
     *
     *  @return The value of the function at the given parameter value
     */
    float value( float x ) const;

  private: // data
    // edges of x bins
    const std::vector<float> m_xedges;

    // width of y bins
    float m_width;
  };

} // namespace Track

#endif // TRACKKERNEL_TRACK1DTABFUNC_H
