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
/** @class TrackCompetition TrackCompetition.h
 *
 *  Copy a container of tracks. By default do not copy tracks that failed the fit
 *
 *  @author M.Needham
 *  @date   30/05/2006
 */

#include "Event/Track.h"
#include "LHCbAlgs/Transformer.h"
#include <numeric>
#include <string>

using namespace LHCb;

class TrackCompetition : public LHCb::Algorithm::Transformer<Tracks( const Tracks& )> {

public:
  TrackCompetition( const std::string& name, ISvcLocator* pSvcLocator );
  Tracks operator()( const Tracks& ) const override;

private:
  Gaudi::Property<double> m_fracUsed{this, "fracUsed", 0.25};
};

DECLARE_COMPONENT( TrackCompetition )

TrackCompetition::TrackCompetition( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer{
          name, pSvcLocator, {"inputLocation", TrackLocation::Default}, {"outputLocation", TrackLocation::Default}} {}

Tracks TrackCompetition::operator()( const Tracks& in ) const {

  Tracks outCont;
  outCont.reserve( in.size() );

  // sort
  auto inCont = std::vector<const LHCb::Track*>{in.begin(), in.end()};
  std::sort( inCont.begin(), inCont.end(), []( const LHCb::Track* lhs, const LHCb::Track* rhs ) {
    return std::tuple{lhs->nLHCbIDs(), rhs->chi2()} > std::tuple{rhs->nLHCbIDs(), lhs->chi2()};
  } );

  // TODO: maybe a set is better? This is O(N^2)..
  // or keep used sorted (N log N), and use the fact that ids
  // is also sorted so the lookup is done O(N) instead of O(N^2)
  // or use a bloom filter (which would need to be updated...)
  std::vector<LHCb::LHCbID> used;
  used.reserve( 4000 );
  for ( const auto& t : inCont ) {
    const auto& ids   = t->lhcbIDs();
    auto        nUsed = std::count_if( ids.begin(), ids.end(), [&]( const LHCbID& id ) {
      return std::find( used.begin(), used.end(), id ) != used.end();
    } );

    if ( nUsed < m_fracUsed * ids.size() ) {
      used.insert( used.end(), ids.begin(), ids.end() );
      outCont.insert( new LHCb::Track( *t ) );
    }
  } // for each

  return outCont;
}
