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
#ifndef PRGECFILTER_H
#define PRGECFILTER_H 1

// Include files
// from Gaudi
#include "Event/RawEvent.h"
#include "GaudiAlg/FilterPredicate.h"

/** @class PrGECFilter PrGECFilter.h
 *  \brief give decision concerning GEC
 */
class PrGECFilter
    : public Gaudi::Functional::FilterPredicate<bool( const LHCb::RawBank::View&, const LHCb::RawBank::View& )> {

public:
  /// Standard constructor
  PrGECFilter( const std::string& name, ISvcLocator* pSvcLocator );

  /// Algorithm execution
  bool operator()( const LHCb::RawBank::View&, const LHCb::RawBank::View& ) const override;

private:
  Gaudi::Property<unsigned int> m_nFTUTClusters{this, "NumberFTUTClusters", 0};
  Gaudi::Property<unsigned int> m_clusterMaxWidth{this, "ClusterMaxWidth", 4, "Maximal cluster width"};

  mutable Gaudi::Accumulators::Counter<> m_eventsProcessedCounter{this, "Nb Events Processed"};
  mutable Gaudi::Accumulators::Counter<> m_eventsRejectedCounter{this, "Nb events removed"};
};
#endif // PRGECFILTER_H
