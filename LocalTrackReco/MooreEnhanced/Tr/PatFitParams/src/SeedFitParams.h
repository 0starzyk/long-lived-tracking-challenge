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
#ifndef SEEDFITPARAMS_H
#define SEEDFITPARAMS_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTupleAlg.h"

#include "FitTool.h"
#include "FwdParameters.h"

/** @class SeedFitParams SeedFitParams.h
 *
 *
 *  @author Olivier Callot
 *  @date   2006-12-08
 */
class SeedFitParams : public GaudiTupleAlg {
public:
  /// Standard constructor
  SeedFitParams( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution
  StatusCode finalize() override;   ///< Algorithm finalization

protected:
private:
  FitTool*    m_fitTool;
  std::string m_tupleName;
  double      m_zRef;
  double      m_zSeed;
  double      m_zTT;
  int         m_nEvent;
  int         m_nTrack;

  std::vector<double> m_momentumScaleParams;
  std::vector<double> m_initialArrowParams;
  std::vector<double> m_zMagParams;

  FwdParameters m_momentumScalePar;
  FwdParameters m_initialArrowPar;
  FwdParameters m_zMagPar;

  std::vector<double> m_dRatio;
  FwdParameters       m_dRatioPar;

  std::vector<double> m_yCorrection;
  FwdParameters       m_yCorrectionPar;
};
#endif // SEEDFITPARAMS_H
