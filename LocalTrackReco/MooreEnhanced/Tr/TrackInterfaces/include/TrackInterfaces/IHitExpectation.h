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
#include "GaudiKernel/IAlgTool.h"
#include "Kernel/LHCbID.h"

#include <vector>

#include "DetDesc/IGeometryInfo.h"

/**
 *  How many hits (of a given type) do we expect on a track ?
 *
 *  @author M.Needham
 *  @date   07/09/2007
 */
struct IHitExpectation : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IHitExpectation, 2, 0 );

  /** small struct returning hit info....
   * nExpected --> number of expected hits
   * number of the expected hits that are found
   * likelihood -> what these hits would contribute to likelihood
   * (In fact this only matters for OT where eff = function of r)
   */
  typedef struct {
    unsigned int nExpected;
    unsigned int nFound;
    double       likelihood;
  } Info;

  /** Returns number of hits expected, from zFirst to inf
   *
   *  @param aTrack Reference to the Track to test
   *
   *  @return unsigned int number of hits expected
   */
  virtual unsigned int nExpected( const LHCb::Track& aTrack, IGeometryInfo const& geometry ) const = 0;

  /** Returns number of hits expected
   *
   *  @param aTrack Reference to the Track to test
   *
   *  @return Info info including likelihood
   */
  virtual Info expectation( const LHCb::Track& aTrack, IGeometryInfo const& geometry ) const = 0;

  /** Collect all the expected hits
   *
   * @param aTrack Reference to the Track to test
   * @param hits collected lhcbIDs
   *
   **/
  virtual void collect( const LHCb::Track& aTrack, std::vector<LHCb::LHCbID>& ids,
                        IGeometryInfo const& geometry ) const = 0;
};
