/*****************************************************************************\
* (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#ifndef TRACKCHECKERBASE_H
#define TRACKCHECKERBASE_H 1

// from Gaudi
#include "GaudiAlg/GaudiHistoAlg.h"

#include <string>

// interfaces
#include "GaudiKernel/IMagneticFieldSvc.h"
#include "MCInterfaces/IIdealStateCreator.h"
#include "MCInterfaces/IMCReconstructible.h"
#include "MCInterfaces/ITrackGhostClassification.h"
#include "MCInterfaces/IVisPrimVertTool.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/Track.h"

/** @class TrackCheckerBase TrackCheckerBase.h "TrackCheckers/TrackCheckerBase"
 *
 *  Base class for track monitoring: essentially a 'box' of common tools

 *  @author M. Needham.
 *  @date   7-5-2007
 */

class TrackCheckerBase : public GaudiHistoAlg {

public:
  /** Standard construtor */
  using GaudiHistoAlg::GaudiHistoAlg;

  /** Algorithm initialization */
  StatusCode initialize() override;

public:
  /** Get a pointer to Magnetic field service
   *  @return field service
   */
  const IMagneticFieldSvc* fieldSvc() const;

  /** Get a pointer to the track selection tool
   *  @return field service
   */
  const IMCReconstructible* selector() const;

  /** Get a pointer to the idealStateCreator
   *  @return IdealStateCreator
   */
  const IIdealStateCreator* idealStateCreator() const;

  /** Get a pointer to the track extrapolator
   *  @return extrapolator
   */
  const ITrackExtrapolator* extrapolator() const;

  /** small struct for link info */
  struct LinkInfo {
    LinkInfo( const LHCb::Track* t, unsigned int c, double p ) : track( t ), clone( c ), purity( p ) {}
    LinkInfo() = default;
    const LHCb::Track* track{nullptr};
    unsigned int       clone{0};
    double             purity{-1};
  };

  /** link to truth
   * @param  aTrack track
   * @return linked particle
   */
  const LHCb::MCParticle* mcTruth( const LHCb::Track&, const LHCb::MCParticles&, const LHCb::LinksByKey& ) const;

  /** Selected as type
   *
   * @return bool
   */
  bool selected( const LHCb::MCParticle* particle ) const;

  /** Whether to split by algorithm
   *  @return splitByAlgorithm true or false
   */
  bool splitByAlgorithm() const;

  /** Whether to split by algorithm
   *  @return splitByType true or false
   */
  bool splitByType() const { return m_splitByType.value(); }

  /** Pointer to the visible primary vertex tool
   *  @return IVisPrimVertTool
   */
  const IVisPrimVertTool* visPrimVertTool() const;

  /** Pointer to ghost classification tool
   *  @return ITrackGhostClassification
   */
  const ITrackGhostClassification* ghostClassification() const;

  /** Is a b child ? ie has b quark somewhere in history
   * @param  mcPart MC particle
   * @return bool true/false
   */
  bool bAncestor( const LHCb::MCParticle* mcPart ) const;

  /** Is a lambda/ks
   * @param  mcPart MC particle
   * @return bool true/false
   */
  bool ksLambdaAncestor( const LHCb::MCParticle* mcPart ) const;

  /** are all stable daughters of this particle reconstructible?
   * @param  mcPart MC particle
   * @return bool true/false
   */

  bool allDaughtersReconstructible( const LHCb::MCParticle* mcPart ) const;

  bool bAncestorWithReconstructibleDaughters( const LHCb::MCParticle* mcPart ) const;

private:
  Gaudi::Property<std::string>    m_selectionCriteria{this, "SelectionCriteria", "ChargedLong"};
  IMCReconstructible::RecCategory m_recCat; ///<  Pointer to selector

  ServiceHandle<IMagneticFieldSvc>      m_pIMF{this, "MagneticFieldService",
                                          "MagneticFieldSvc"}; ///<  Pointer to the magn. field service
  ToolHandle<ITrackGhostClassification> m_ghostClassification{
      this, "GhostTool", "LongGhostClassification/GhostTool"};                               ///< Pointer to ghost tool
  ToolHandle<IMCReconstructible> m_selector{this, "Selector", "MCReconstructible/Selector"}; ///<  Pointer to selector
  PublicToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator",
                                                      "TrackMasterExtrapolator"}; ///<  Pointer to extrapolator
  PublicToolHandle<IIdealStateCreator> m_stateCreator{this, "StateCreator",
                                                      "IdealStateCreator"}; ///<  IdealStateCreator
  PublicToolHandle<IVisPrimVertTool>   m_visPrimVertTool{this, "VisPrimVertTool",
                                                       "VisPrimVertTool"}; ///< Visible primary vertices..

  Gaudi::Property<bool> m_splitByAlgorithm{this, "SplitByAlgorithm", false};
  Gaudi::Property<bool> m_splitByType{this, "SplitByType", false};
};

inline const IMagneticFieldSvc* TrackCheckerBase::fieldSvc() const { return m_pIMF.get(); }

inline const IMCReconstructible* TrackCheckerBase::selector() const { return m_selector.get(); }

inline bool TrackCheckerBase::selected( const LHCb::MCParticle* particle ) const {
  return selector()->isReconstructibleAs( m_recCat, particle );
}

inline const ITrackExtrapolator* TrackCheckerBase::extrapolator() const { return m_extrapolator.get(); }

inline const IIdealStateCreator* TrackCheckerBase::idealStateCreator() const { return m_stateCreator.get(); }

inline const IVisPrimVertTool* TrackCheckerBase::visPrimVertTool() const { return m_visPrimVertTool.get(); }

inline bool TrackCheckerBase::splitByAlgorithm() const { return m_splitByAlgorithm.value(); }

inline const ITrackGhostClassification* TrackCheckerBase::ghostClassification() const {
  return m_ghostClassification.get();
}

#endif // TRACKCHECKERBASE_H
