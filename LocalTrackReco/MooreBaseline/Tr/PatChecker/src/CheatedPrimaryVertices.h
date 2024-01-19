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
#ifndef CHEATEDPRIMARYVERTICES_H
#define CHEATEDPRIMARYVERTICES_H 1
// Include files:
// from Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
// From PhysEvent
#include "Event/MCParticle.h"
#include "Event/Particle.h"
#include "Event/ProtoParticle.h"
#include "Event/RecVertex.h"
#include "MCInterfaces/IForcedBDecayTool.h"

/** @class PVEff CheatedPrimaryVertices.h
 *
 *
 *  @author Marcin Kucharczyk
 *  @date   2008-12-07
 */

class CheatedPrimaryVertices : public GaudiAlgorithm {
public:
  // Standard constructor
  CheatedPrimaryVertices( const std::string& name, ISvcLocator* pSvcLocator );
  StatusCode execute() override;

private:
  LHCb::Tracks*      m_inputTracks;
  LHCb::RecVertices* m_outputVertices;
  std::string        m_inputTracksName;
  std::string        m_outputVerticesName;
  int                fromMCVertex( const LHCb::MCParticle* mcParticle, const LHCb::MCVertex* mcVertex );
  double             chi2Tr( const LHCb::Track* pvTr, const LHCb::RecVertex* invt );
};
#endif // CHEATEDPRIMARYVERTICES_H
