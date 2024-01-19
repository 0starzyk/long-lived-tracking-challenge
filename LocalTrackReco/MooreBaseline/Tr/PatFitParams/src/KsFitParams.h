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
#ifndef KSFITPARAMS_H
#define KSFITPARAMS_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/Point3DTypes.h"

#include "FitTool.h"
#include "FwdParameters.h"

/** @class KsFitParams KsFitParams.h
 *  PArameterize the KShort tracks
 *
 *  @author Olivier Callot
 *  @date   2002-11-02
 */
class KsFitParams : public GaudiTupleAlg {
public:
  /// Standard constructor
  KsFitParams( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~KsFitParams(); ///< Destructor

  StatusCode initialize() override; ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution
  StatusCode finalize() override;   ///< Algorithm finalization

protected:
private:
  FitTool*    m_fitTool;
  std::string m_tupleName;
  double      m_zTT1;
  double      m_zRef;

  std::vector<double> m_zMagParams;
  std::vector<double> m_momParams;

  FwdParameters m_zMagPar;
  FwdParameters m_momPar;

  int m_nEvent;
  int m_nTrack;
};
#endif // KSFITPARAMS_H
