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
// Include files

// local
#include "PrGECFilter.h"
#include "FTDAQ/FTDAQHelper.h"
#include "UTDAQ/UTDAQHelper.h"

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( PrGECFilter )

/// Standard constructor, initializes variables
PrGECFilter::PrGECFilter( const std::string& name, ISvcLocator* pSvcLocator )
    : FilterPredicate( name, pSvcLocator,
                       {KeyValue{"FTRawBanks", "DAQ/RawBanks/FTCluster"}, KeyValue{"UTRawBanks", "DAQ/RawBanks/UT"}} ) {
}

bool PrGECFilter::operator()( const LHCb::RawBank::View& ftBanks, const LHCb::RawBank::View& utBanks ) const {
  ++m_eventsProcessedCounter;

  // do not work for nothing !
  if ( m_nFTUTClusters <= 0 ) { return true; }

  // check UT clusters
  auto nbUTClusters = LHCb::UTDAQ::nbUTClusters( utBanks, m_nFTUTClusters );
  if ( !nbUTClusters ) {
    ++m_eventsRejectedCounter;
    return false;
  }

  // check FT clusters
  auto nbFTClusters = LHCb::FTDAQ::nbFTClusters( ftBanks );
  if ( nbFTClusters + nbUTClusters.value() < m_nFTUTClusters ) { return true; }
  ++m_eventsRejectedCounter;
  return false;
}
