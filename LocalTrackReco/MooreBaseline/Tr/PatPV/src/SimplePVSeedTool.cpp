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
#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "IPVSeeding.h" // Interface

/**
 *  @author Mariusz Witek
 *  @date   2005-11-19
 */
class SimplePVSeedTool : public extends<GaudiTool, IPVSeeding> {
public:
  /// Standard constructor
  using extends::extends;

  std::vector<Gaudi::XYZPoint> getSeeds( LHCb::span<const LHCb::Track* const> inputTracks,
                                         const Gaudi::XYZPoint&               beamspot,
                                         IGeometryInfo const&                 geometry ) const override;
};

DECLARE_COMPONENT( SimplePVSeedTool )

std::vector<Gaudi::XYZPoint> SimplePVSeedTool::getSeeds( LHCb::span<const LHCb::Track* const> inputTracks,
                                                         const Gaudi::XYZPoint& beamspot, IGeometryInfo const& ) const {

  std::vector<Gaudi::XYZPoint> seeds;
  if ( inputTracks.size() < 3 ) return seeds;

  if ( msgLevel( MSG::DEBUG ) ) { debug() << " Beam spot is ignored. BS: " << beamspot << endmsg; }

  seeds.emplace_back( 0., 0., 0. );
  return seeds;
}
