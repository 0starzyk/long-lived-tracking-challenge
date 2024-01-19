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
#include "GaudiAlg/GaudiAlgorithm.h"
#include "gsl/gsl_sys.h"
#include <cmath>

//-----------------------------------------------------------------------------
// Implementation file for class : HltInsertTrackErrParam
//
// 2007-07-16 : Patrick Koppenburg
//-----------------------------------------------------------------------------

/** @class HltInsertTrackErrParam HltInsertTrackErrParam.h
 *
 *  Reimplementation of Hugo's TrgInsertTrackErrParam.
 *  Hacks error matrix of tracks according to 1/Pt.
 *  See note <A NAME="http://cdsweb.cern.ch/record-restricted/832444/">LHCb-2005-012</a> (internal).
 *
 *  @tod: Should this be an HltAlgorithm?
 *
 *  @author Patrick Koppenburg, Hugo Ruiz
 *  @date   2007-07-16
 */
class HltInsertTrackErrParam : public GaudiAlgorithm {
public:
  /// Standard constructor
  using GaudiAlgorithm::GaudiAlgorithm;

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution

private:
  StatusCode insertParamInTrack( LHCb::Track* tr );
  StatusCode insertParamInState( LHCb::State* s );
  void       printParams( LHCb::span<const double>, std::string );

  Gaudi::Property<std::vector<double>> m_xParams{this,
                                                 "XParams",
                                                 {0.0214 * Gaudi::Units::mm, -0.000801 * Gaudi::Units::mm,
                                                  0.0108 * Gaudi::Units::mm, -0.00122 * Gaudi::Units::mm,
                                                  0.0000547 * Gaudi::Units::mm}};
  Gaudi::Property<std::vector<double>> m_yParams{this, "YParams",
                                                 m_xParams.value()}; // y is equal to x by default ///< y parameters
  DataObjectReadHandle<LHCb::Tracks>   m_inputLocation{this, "InputLocation", "Undefined"}; ///< Location for input
  /// Location for output. If empty, tracks will be modified in situ
  DataObjectWriteHandle<LHCb::Tracks> m_outputLocation{this, "OutputLocation", ""};
  Gaudi::Property<bool>               m_keepMomentum{this, "KeepMomentumCovariance", true,
                                       "Keep Momentum terms of the covariance matrix"}; ///< keep momentum covariance
                                                                                                      ///< terms
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( HltInsertTrackErrParam )

//=============================================================================
// Initialization
//=============================================================================
StatusCode HltInsertTrackErrParam::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;              // error printed already by GaudiAlgorithm

  printParams( m_xParams.value(), "X" );
  printParams( m_yParams.value(), "Y" );

  info() << "Tracks will be taken from " << m_inputLocation << endmsg;
  info() << " and stored in ";
  if ( !m_outputLocation.objKey().empty() )
    info() << m_outputLocation << endmsg;
  else
    info() << "the same location" << endmsg;

  return StatusCode::SUCCESS;
}

//=========================================================================
//
//=========================================================================
void HltInsertTrackErrParam::printParams( LHCb::span<const double> params, std::string coord ) {
  info() << "Your Pt-dependent " << coord << " error parameterisation will be " << endmsg;
  int p = 0;
  info() << "Error on " << coord << " = ";
  for ( auto i : params ) {
    if ( ( p > 0 ) && ( i > 0. ) ) info() << " +";
    info() << i;
    if ( p > 0 ) info() << "*pt^(-" << p << ")";
    ++p;
  }
  info() << endmsg;
}
//=============================================================================
// Main execution
//=============================================================================
StatusCode HltInsertTrackErrParam::execute() {

  if ( msgLevel( MSG::DEBUG ) ) debug() << "==> Execute" << endmsg;

  /// @todo Should use some HltAlgorithm method to get access to tracks anywhere

  LHCb::Tracks* tracks = m_inputLocation.getIfExists();
  if ( !tracks ) {
    setFilterPassed( false );
    Warning( "No tracks at " + m_inputLocation.objKey() ).ignore();
    return StatusCode::SUCCESS;
  }
  auto newTracks = ( !m_outputLocation.objKey().empty() ? std::make_unique<LHCb::Tracks>()
                                                        : std::unique_ptr<LHCb::Tracks>{nullptr} );

  if ( msgLevel( MSG::DEBUG ) ) debug() << "Loaded " << tracks->size() << " tracks" << endmsg;

  for ( auto& track : *tracks ) {
    if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Looping on track " << track->key() << endmsg;

    LHCb::Track* tk = ( newTracks ? new LHCb::Track( *track ) : track );

    StatusCode sc = insertParamInTrack( tk );
    if ( !sc ) {
      Warning( " call to insertParamInTrack failed -- abandoning all new tracks" ).ignore();
      if ( newTracks ) {
        std::for_each( newTracks->begin(), newTracks->end(), std::default_delete<LHCb::Track>() );
        newTracks->clear();
      }
      return sc;
    }
    if ( newTracks ) newTracks->insert( tk );
  }
  if ( newTracks ) m_outputLocation.put( std::move( newTracks ) );

  return StatusCode::SUCCESS;
}
//=========================================================================
//  Loop over states
//=========================================================================
StatusCode HltInsertTrackErrParam::insertParamInTrack( LHCb::Track* tr ) {
  for ( const auto& s : tr->states() ) {
    LHCb::State::Location loc = s->location();
    if ( ( loc == LHCb::State::Location::ClosestToBeam ) || ( loc == LHCb::State::Location::FirstMeasurement ) ||
         ( loc == LHCb::State::Location::EndVelo ) ) {
      if ( msgLevel( MSG::DEBUG ) ) debug() << "State at " << loc << " " << s->position().z() << endmsg;
      /// @todo add counter later
      StatusCode sc = insertParamInState( tr->stateAt( loc ) );
      if ( !sc ) return sc;
    }
  }
  return StatusCode::SUCCESS;
}
//=========================================================================
//  The actual thing
//=========================================================================
StatusCode HltInsertTrackErrParam::insertParamInState( LHCb::State* state ) {

  if ( 0 == gsl_fcmp( state->qOverP(), 0, 1.e-9 ) ) { // that needs a 1000 TeV track
    err() << "Track state at " << state->position() << " (" << state->location() << ") has q/P = " << state->qOverP()
          << endmsg;
    return StatusCode::FAILURE;
  }

  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << state->position().x() << " " << state->position().y() << " " << state->position().z() << " "
            << 1. / state->qOverP() << endmsg;
    debug() << "Pt: " << state->pt() << endmsg;
    debug() << "Existent matrix: \n" << state->covariance() << endmsg;
  }

  double invPt = 1. / ( fabs( state->pt() / Gaudi::Units::GeV ) );

  if ( msgLevel( MSG::DEBUG ) ) debug() << "Inverse of Pt: " << invPt << endmsg;

  double sigmaX = 0;
  for ( unsigned int p = 0; p < m_xParams.size(); ++p ) {
    sigmaX += m_xParams[p] * std::pow( invPt, static_cast<int>( p ) );
  }
  double sigmaY = 0;
  for ( unsigned int p = 0; p < m_yParams.size(); ++p ) {
    sigmaY += m_yParams[p] * std::pow( invPt, static_cast<int>( p ) );
  }

  if ( msgLevel( MSG::DEBUG ) ) debug() << "sigmaX, sigmaY: " << sigmaX << " " << sigmaY << endmsg;

  Gaudi::TrackSymMatrix newMatrix;
  if ( m_keepMomentum.value() ) newMatrix = state->covariance();

  newMatrix( 0, 0 ) = sigmaX * sigmaX;
  newMatrix( 1, 1 ) = sigmaY * sigmaY;
  newMatrix( 1, 0 ) = 0; // xy term
  newMatrix( 0, 1 ) = 0; // xy term

  /// @todo could be done with insert...
  if ( m_keepMomentum.value() ) { // set all cros terms to 0
    newMatrix( 0, 2 ) = 0;
    newMatrix( 0, 3 ) = 0;
    newMatrix( 0, 4 ); // x - tx, ty, q/P
    newMatrix( 1, 2 ) = 0;
    newMatrix( 1, 3 ) = 0;
    newMatrix( 1, 4 ); // y - tx, ty, q/P
    newMatrix( 2, 0 ) = 0;
    newMatrix( 3, 0 ) = 0;
    newMatrix( 4, 0 ); // x - tx, ty, q/P
    newMatrix( 2, 1 ) = 0;
    newMatrix( 3, 1 ) = 0;
    newMatrix( 4, 1 ); // y - tx, ty, q/P
  } else {             // keep diagonals
    newMatrix( 2, 2 ) = state->covariance()( 2, 2 );
    newMatrix( 3, 3 ) = state->covariance()( 3, 3 );
    newMatrix( 4, 4 ) = state->covariance()( 4, 4 );
  }

  if ( msgLevel( MSG::DEBUG ) ) debug() << "New matrix: \n" << newMatrix << endmsg;

  state->setCovariance( newMatrix );

  return StatusCode::SUCCESS;
}
