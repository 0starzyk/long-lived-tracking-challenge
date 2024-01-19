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
#include "Event/Measurement.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/VectorMap.h"
#include "TrackInterfaces/ITrackProjector.h"
#include "TrackInterfaces/ITrackProjectorSelector.h"
#include <map>
#include <string>

using namespace Gaudi;
using namespace LHCb;

/** @class TrackProjectorSelector TrackProjectorSelector.h TrackProjectors/TrackProjectorSelector.h
 *
 *  TrackProjectorSelector decides which projection to use for
 *  a given (type of) measurement
 *
 *  @author Gerhard Raven
 *  @date   2006-06-22
 */
class TrackProjectorSelector : public extends<GaudiTool, ITrackProjectorSelector> {

public:
  /// Standard constructor
  TrackProjectorSelector( const std::string& type, const std::string& name, const IInterface* parent );

  StatusCode initialize() override;

  ITrackProjector* projector( const LHCb::Measurement& ) const override;

private:
  std::map<LHCb::Measurement::Type, std::string>                   m_projNames;
  GaudiUtils::VectorMap<LHCb::Measurement::Type, ITrackProjector*> m_projectors;
};

// Declaration of the Tool Factory
DECLARE_COMPONENT( TrackProjectorSelector )
//-----------------------------------------------------------------------------
/// Standard constructor, initializes variables
//-----------------------------------------------------------------------------
TrackProjectorSelector::TrackProjectorSelector( const std::string& type, const std::string& name,
                                                const IInterface* parent )
    : base_class( type, name, parent ) {
  // FIXME: as soon as the warnings in GaudiAlg on multiple tools are gone, we
  //       can remove the different names for the
  declareProperty( "VP", m_projNames[Measurement::Type::VP] = "TrajProjector<TrajProj::VP>/TrajVPProjector" );
  declareProperty( "UT", m_projNames[Measurement::Type::UT] = "TrajProjector<TrajProj::UT>/TrajUTProjector" );
  declareProperty( "FT", m_projNames[Measurement::Type::FT] = "TrajProjector<TrajProj::FT>/TrajFTProjector" );
  declareProperty( "Muon", m_projNames[Measurement::Type::Muon] = "TrajProjector<TrajProj::Muon>/TrajMuonProjector" );
}

//-----------------------------------------------------------------------------
/// Initialize
//-----------------------------------------------------------------------------
StatusCode TrackProjectorSelector::initialize() {
  return base_class::initialize().andThen( [&] {
    m_projectors.clear();
    std::for_each( m_projNames.begin(), m_projNames.end(),
                   [&]( const auto& i ) { m_projectors.insert( i.first, this->tool<ITrackProjector>( i.second ) ); } );
    if ( msgLevel( MSG::DEBUG ) ) {
      std::for_each( m_projNames.begin(), m_projNames.end(),
                     [&]( const auto& i ) { debug() << " projector for " << i.first << " : " << i.second << endmsg; } );
    }
  } );
}

//-----------------------------------------------------------------------------
/// select the projector;
/// TODO: return an object which represents the binding of the measurement
///       and projector (taking care of any downcasting here, when creating
///       such an object)
//-----------------------------------------------------------------------------
ITrackProjector* TrackProjectorSelector::projector( const LHCb::Measurement& m ) const {
  auto i = m_projectors.find( m.type() );
  return i != m_projectors.end() ? i->second : nullptr;
}
