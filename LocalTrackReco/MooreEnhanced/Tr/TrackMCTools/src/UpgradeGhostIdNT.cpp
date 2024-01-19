/*****************************************************************************\
* (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "Event/GhostTrackInfo.h"
#include "Event/MCParticle.h"
#include "Event/Track.h"
#include "GaudiAlg/GaudiTupleTool.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/ToolHandle.h"
#include "Linker/LinkerTable.h"
#include "MCInterfaces/ITrackGhostClassification.h"
#include "TMath.h"
#include "TrackInterfaces/IGhostProbability.h"

//-----------------------------------------------------------------------------
// Implementation file for class : UpgradeGhostIdNT
//
//  To Generate ntuples for the ghost probability study
//
//-----------------------------------------------------------------------------

/*@class UpgradeGhostIdNT UpgradeGhostIdNT.h

  To Generate ntuples for the ghost probability study

*/
class UpgradeGhostIdNT : public extends<GaudiTupleTool, IGhostProbability> {
public:
  /// Standard constructor
  using extends::extends;

  StatusCode                    execute( LHCb::Track& aTrack ) const override;
  std::vector<std::string_view> variableNames( LHCb::Track::Types type ) const override {
    return m_ghostTool->variableNames( type );
  };
  std::vector<float> netInputs( LHCb::Track& aTrack ) const override { return m_ghostTool->netInputs( aTrack ); };

private:
  ToolHandle<IGhostProbability> m_ghostTool{this, "Tool", "UpgradeGhostId"};
};

DECLARE_COMPONENT( UpgradeGhostIdNT )

//=============================================================================
StatusCode UpgradeGhostIdNT::execute( LHCb::Track& aTrack ) const {
  if ( m_ghostTool->execute( aTrack ).isFailure() ) { return StatusCode::SUCCESS; }

  std::vector<float>            variables = m_ghostTool->netInputs( aTrack );
  std::vector<std::string_view> varnames  = m_ghostTool->variableNames( aTrack.type() );
  std::vector<std::string_view> varnames_new;

  SmartDataPtr<LHCb::LinksByKey> links( evtSvc(), LHCb::LinksByKey::linkerName( "Rec/Track/Best" ) );
  if ( links ) links->resolveLinks( evtSvc() );
  const LinkerTable<LHCb::Track, LHCb::MCParticle> table{links};
  auto                                             range = table.relations( &aTrack );

  Tuples::Tuple tup = GaudiTupleTool::nTuple( "tracks", CLID_ColumnWiseTuple );

  if ( LHCb::Track::Types::Long == aTrack.type() ) {
    for ( unsigned ivar = 0; ivar < varnames.size(); ivar++ ) {
      std::string_view tmp_var = varnames[ivar];
      if ( ivar != ( varnames.size() - 3 ) ) {
        varnames_new.push_back( tmp_var );
      } else {
        varnames_new.push_back( tmp_var );
        varnames_new.push_back( "TRACK_NDOF" );
      }
    }
  } else {
    for ( unsigned ivar = 0; ivar < varnames.size(); ivar++ ) {
      std::string_view tmp_var = varnames[ivar];
      varnames_new.push_back( tmp_var );
    }
  }
  if ( varnames_new.size() != variables.size() )
    fatal() << "ALARM  " << varnames_new.size() << " != " << variables.size() << "  " << endmsg;
  for ( unsigned i = 0; i < varnames_new.size(); ++i ) {
    tup->column( std::string( varnames_new[i] ).c_str(), variables[i] ).ignore();
  }

  tup->column( "ghostprob", (float)aTrack.ghostProbability() ).ignore();
  tup->column( "tracks_PP_TrackHistory", (float)aTrack.history() ).ignore();
  tup->column( "tracks_TRACK_Type", (float)aTrack.type() ).ignore();
  tup->column( "tracks_assoc", (float)(int)( !( range.empty() ) ) ).ignore();
  tup->column( "mctruepid", (float)(int)( ( range.empty() ) ? 0 : ( range.begin()->to()->particleID().pid() ) ) )
      .ignore();

  return tup->write();
}
