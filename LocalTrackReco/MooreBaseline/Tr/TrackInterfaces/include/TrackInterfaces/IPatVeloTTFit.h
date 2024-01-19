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
#ifndef TRACKINTERFACES_IPATVELOTTFIT_H
#define TRACKINTERFACES_IPATVELOTTFIT_H 1

#include "Event/Track.h"
#include "GaudiKernel/IAlgTool.h"

class PatTTHit;

/** @class IPatVeloTTFit IPatVeloTTFit.h
 *
 * provide a convenient interface to the internal fit used in the PatVeloTTFit
 * algorithm in the pattern reco
 *
 * @author Pavel Krokovny <krokovny@physi.uni-heidelberg.de>
 * @date   2009-03-09
 */
struct IPatVeloTTFit : extend_interfaces<IAlgTool> {

  DeclareInterfaceID( IPatVeloTTFit, 2, 0 );

  virtual StatusCode fitVTT( LHCb::Track& track ) const = 0;

  virtual void finalFit( const std::vector<PatTTHit*>& theHits, const std::array<float, 8>& vars,
                         std::array<float, 3>& params ) const = 0;
};
#endif // INCLUDE_IPATVELOTTFIT_H
