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
#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "Kernel/LHCbID.h"
#include "TrackInterfaces/ITrackCloneFinder.h"

//-----------------------------------------------------------------------------
// Implementation file for class : TrackCloneFinder
//
// 2005-12-08 : Eduardo Rodrigues
// 2005-05-05 : Adrian Perieanu
// 2009-05-08: Georg Krocker
// 2009-09-10: Kostyantyn Holubyev
//-----------------------------------------------------------------------------

/** @class TrackCloneFinder TrackCloneFinder.h
 *
 *
 *  @author Eduardo Rodrigues
 *  @date   2005-12-08
 *  Modified for speed & clone rate
 *  @author Adrian Perieanu
 *  @date   2008-05-05
 */
class TrackCloneFinder : public extends<GaudiTool, ITrackCloneFinder> {
public:
  /// Standard constructor
  using extends::extends;

  /** Compare two input Tracks and find whether one is a clone
   *  of the other based on some "overlap criteria".
   *  Note: the method ignores whether the Tracks themselves have been
   *        previously flagged as clones! It merely does a comparison.
   *  @param  track1 input 1st track
   *  @param  track2 input 2nd track
   *  @param  setFlag input parameter indicates whether the clone track
   *          is to be set as such (default = false)
   */
  bool areClones( const LHCb::Track& track1, const LHCb::Track& track2 ) const override;

  /** Compare two input Tracks and find whether one is a clone
   *  of the other based on some "overlap criteria".
   *  The corresponding flag may be set accordingly (NOT DONE BY DEFAULT)
   *  depending on the value of the "setFlag" argument.
   *  Note: the method ignores whether the Tracks themselves have been
   *        previously flagged as clones! It merely does a comparison.
   *  @param  track1 input 1st track
   *  @param  track2 input 2nd track
   *  @param  setFlag input parameter indicates whether the clone track
   *          is to be set as such (default = false)
   */
  bool flagClones( LHCb::Track& track1, LHCb::Track& track2 ) const override;

private:
  /** Calculate the number of common hits of a given LHCb type
   *  between two input Tracks.
   *  Note: hits can here mean either Measurements or LHCbIDs,
   *        depending on the settings for the "CompareAtLHCbIDsLevel"
   *        tool property:
   *        1) CompareAtLHCbIDsLevel = false (default):
   *           the comparison is based on the Measurements and it is
   *           therefore assumed that the Tracks have been fitted
   *        2) CompareAtLHCbIDsLevel = true:
   *           the comparison is based on the LHCbIDs and is therefore
   *           done on the Tracks as output by the pattern recognition
   *  @param[in]  track1      input 1st track
   *  @param[in]  track2      input 2nd track
   *  @param[out] nCommonHits number of shared hits
   */
  bool areTracksClose( const LHCb::Track& track1, const LHCb::Track& track2 ) const;

private:
  // tool properties
  // ---------------
  /* fraction of hits in common that defines the overlap between the
   * two tracks and their being clones or not */
  Gaudi::Property<double> m_matchingFraction{this, "MatchingFraction", 0.7};
  Gaudi::Property<double> m_matchingFractionT{this, "MatchingFractionT", 0.50};
  Gaudi::Property<bool>   m_compareLDT{this, "CompareLDT", false};
  // In some cases, such as HLT we want to search for clones only if the tracks
  // are close. To do so, set RestrictedSearch to true
  Gaudi::Property<bool> m_restrictedSearch{this, "RestrictedSearch", false};

  // The search windows for which two tracks in the velo are supposed to be
  // close.
  // TODO: The Value of this parameters is currently determined iteratively
  // and may bot be optimal. This has to be investigated further.
  Gaudi::Property<double> m_xseparationV{this, "VeloXSeparation", 100};
  Gaudi::Property<double> m_yseparationV{this, "VeloYSeparation", 100};
  Gaudi::Property<double> m_txseparationV{this, "VeloTXSeparation", 3e-2};
  Gaudi::Property<double> m_tyseparationV{this, "VeloTYSeparation", 2e-2};
};

DECLARE_COMPONENT( TrackCloneFinder )

//=============================================================================
// Compare two input Tracks and find whether one is a clone
// of the other based on some "overlap criteria".
// The corresponding flags are set accordingly.
//=============================================================================
bool TrackCloneFinder::flagClones( LHCb::Track& track1, LHCb::Track& track2 ) const {
  bool theyAreClones = TrackCloneFinder::areClones( track1, track2 );
  if ( theyAreClones ) {
    size_t n1 = track1.nLHCbIDs();
    size_t n2 = track2.nLHCbIDs();
    if ( n1 > n2 ) {
      track2.setFlag( LHCb::Track::Flags::Clone, true );
    } else if ( n2 > n1 ) {
      track1.setFlag( LHCb::Track::Flags::Clone, true );
      // In some cases the tracks might not be fitted so we can not use the chi2 as
      // criterion. So, check if the track is fitted..
    } else if ( track1.fitStatus() == LHCb::Track::FitStatus::Fitted &&
                track2.fitStatus() == LHCb::Track::FitStatus::Fitted ) {
      const double chi1 = track1.chi2PerDoF();
      const double chi2 = track2.chi2PerDoF();
      chi1 < chi2 ? track2.setFlag( LHCb::Track::Flags::Clone, true )
                  : track1.setFlag( LHCb::Track::Flags::Clone, true );
    } else {
      if ( msgLevel( MSG::DEBUG ) )
        debug() << "At least one of the tracks is not fitted, selecting a clone randomly" << endmsg;
      // for fast reco sequence the clone killer is run before the fit,
      // so chi2 is not available. However, in most cases the decision is
      // made using the number n of LHCb IDs, so for the rare case when n1==n2
      // we simply select the track to mark as a clone randomly
      track2.setFlag( LHCb::Track::Flags::Clone, true );
    }
  }
  return theyAreClones;
}

//=============================================================================
// Compare two input Tracks and find whether one is a clone
// of the other based on some "overlap criteria".
// The corresponding flags are set accordingly.
//=============================================================================
bool TrackCloneFinder::areClones( const LHCb::Track& track1, const LHCb::Track& track2 ) const {
  if ( msgLevel( MSG::DEBUG ) ) {
    debug() << "Looking at tracks " << track1.key() << " in " << track1.parent()->name() << " and " << track2.key()
            << " in " << track2.parent()->name() << endmsg;
  }

  // If we want to speed up the time for clonesearch we can look only on clones
  // which are physcially close
  bool theyAreClones( false );
  if ( !m_restrictedSearch || areTracksClose( track1, track2 ) ) {

    size_t nHitsCommon = track1.nCommonLhcbIDs( track2 );
    size_t n1          = track1.nLHCbIDs();
    size_t n2          = track2.nLHCbIDs();
    size_t nTrackMin   = std::min( n1, n2 );
    theyAreClones      = nHitsCommon > m_matchingFraction * nTrackMin;

    if ( !theyAreClones && m_compareLDT && ( track1.type() == LHCb::Track::Types::Long ) &&
         ( track2.type() == LHCb::Track::Types::Downstream || track2.type() == LHCb::Track::Types::Ttrack ) ) {
      theyAreClones = nHitsCommon > m_matchingFractionT * nTrackMin;
    }

    if ( msgLevel( MSG::DEBUG ) ) debug() << "-> areClones = " << theyAreClones << endmsg;
  }
  return theyAreClones;
}

//=============================================================================
// Look if two tracks are really close to each other
//=============================================================================
bool TrackCloneFinder::areTracksClose( const LHCb::Track& tr1, const LHCb::Track& tr2 ) const {
  // We only check wheter the tracks are close in the velo or not so far
  // before we check if they are ghosts, just to save some time

  // first we get the states in the velo
  const LHCb::State* vstate1 = tr1.stateAt( LHCb::State::Location::ClosestToBeam );
  const LHCb::State* vstate2 = tr2.stateAt( LHCb::State::Location::ClosestToBeam );

  // and check if they exist
  if ( tr1.hasStateAt( LHCb::State::Location::ClosestToBeam ) &&
       tr2.hasStateAt( LHCb::State::Location::ClosestToBeam ) ) {
    // We check if one of the states lies outside the velo, then our model
    // of a linear track aproximation is wrong and we want to check the
    // tracks anyway for clones
    if ( vstate1->z() > 990. || vstate2->z() > 990. ) return true;

    // As LHCb::State::Location::ClosestToBeam is not always at the same z position,
    // we extrapolate the tracks linear to a z position in the middle of
    // the two states
    double extrapolateTo = vstate1->z() + ( vstate2->z() - vstate1->z() / 2 );

    // if they are not close together in x,y or if the slopes are very
    // different they should not be clones. We first check for the slopes
    if ( fabs( vstate1->tx() - vstate2->tx() ) > m_txseparationV ) return false;
    if ( fabs( vstate1->ty() - vstate2->ty() ) > m_tyseparationV ) return false;
    // Afterwards we extrapolate the tracks
    if ( fabs( ( vstate1->x() - ( vstate1->tx() * ( vstate1->z() - extrapolateTo ) ) ) -
               ( vstate2->x() - ( vstate2->tx() * ( vstate2->z() - extrapolateTo ) ) ) ) > m_xseparationV )
      return false;
    if ( fabs( ( vstate1->y() - ( vstate1->ty() * ( vstate1->z() - extrapolateTo ) ) ) -
               ( vstate2->y() - ( vstate2->ty() * ( vstate2->z() - extrapolateTo ) ) ) ) > m_yseparationV )
      return false;
  }

  // ok, they are really close so they can be clones
  return true;
}

//=============================================================================
// vim:sw=4:tw=78:ft=cpp
