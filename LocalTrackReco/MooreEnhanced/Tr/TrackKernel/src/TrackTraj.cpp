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
#include "TrackKernel/TrackTraj.h"
#include "Event/FitNode.h"
#include "Event/TrackFitResult.h"
#include "Event/TrackTags.h"
#include "GaudiKernel/reverse.h"
#include "Kernel/STLExtensions.h"
#include "TrackKernel/TrackFunctors.h"
#include "gsl/gsl_integration.h"
#include "range/v3/view/reverse.hpp"
#include "range/v3/view/transform.hpp"
#include <algorithm>

namespace {
  template <typename T>
  decltype( auto ) front( LHCb::span<T, gsl::dynamic_extent> s ) {
    return s.first( 1 )[0];
  }

  template <typename T>
  decltype( auto ) back( LHCb::span<T, gsl::dynamic_extent> s ) {
    return s.last( 1 )[0];
  }

  const auto append_range = []( auto& container, auto&& range ) {
    return container.insert( container.end(), range.begin(), range.end() );
  };
  const auto take_address = ranges::views::transform( []( const auto& i ) { return &i; } );
} // namespace

namespace LHCb {

  inline bool compareStateZ( const State* lhs, const State* rhs ) { return lhs->z() < rhs->z(); }

  inline bool equalStateZ( const State* lhs, const State* rhs ) {
    return std::abs( lhs->z() - rhs->z() ) < TrackParameters::propagationTolerance;
  }

  TrackTraj::TrackTraj( const Track& track, const DeMagnet* magfieldsvc )
      : ZTrajectory<double>{}, m_bfield{0, 0, 0}, m_cachedindex{InvalidCacheIndex} {
    // we rely on the fact that nodes and states in a track are already sorted.

    const auto& ts = track.states();
    // first add the states from the track.nodes(). make sure these
    // are ordered in increasing z.
    auto* fit = fitResult( track );
    if ( fit ) {
      const auto& nodes = fit->nodes();
      if ( !nodes.empty() ) {
        m_states.reserve( ts.size() + nodes.size() );
        if ( nodes.front()->z() < nodes.back()->z() ) {
          // nodes in right order
          for ( const auto& n : nodes ) m_states.push_back( &( n->state() ) );
        } else {
          // nodes in wrong order
          for ( const auto& n : reverse( nodes ) ) m_states.push_back( &( n->state() ) );
        }
      }
    } else {
      m_states.reserve( ts.size() );
    }
    if ( !ts.empty() ) {
      // states on backward tracks are in reverse order
      auto pivot = ( ts.front()->z() < ts.back()->z() ? m_states.insert( m_states.end(), ts.begin(), ts.end() )
                                                      : m_states.insert( m_states.end(), ts.rbegin(), ts.rend() ) );
      std::inplace_merge( m_states.begin(), pivot, m_states.end(), compareStateZ );
    }

    // check states and initialize cache
    init( magfieldsvc );
  }

  TrackTraj::TrackTraj( span<const FitNode* const> nodes, const DeMagnet* magfieldsvc )
      : ZTrajectory<double>{}, m_bfield{0, 0, 0}, m_cachedindex{InvalidCacheIndex} {
    // first add the states from the track.nodes(). make sure these
    // are ordered in increasing z.
    if ( !nodes.empty() ) {
      m_states.reserve( nodes.size() );
      if ( front( nodes )->z() < back( nodes )->z() ) {
        // nodes in right order
        for ( const auto& n : nodes ) m_states.push_back( &( n->state() ) );
      } else {
        // nodes in wrong order
        for ( const auto& n : reverse( nodes ) ) m_states.push_back( &( n->state() ) );
      }
    }

    // check states and initialize cache
    init( magfieldsvc );
  }

  TrackTraj::TrackTraj( span<const LHCb::State* const> states, const DeMagnet* magfieldsvc )
      : ZTrajectory<double>{}
      , m_states{states.begin(), states.end()}
      , m_bfield{0, 0, 0}
      , m_cachedindex{InvalidCacheIndex} {
    // sort
    std::sort( m_states.begin(), m_states.end(), compareStateZ );
    // check states and initialize cache
    init( magfieldsvc );
  }

  TrackTraj::TrackTraj( span<const LHCb::State* const> states, LHCb::Tag::State::AssumeSorted_tag,
                        const DeMagnet*                magfieldsvc )
      : ZTrajectory<double>{}, m_bfield{0, 0, 0}, m_cachedindex{InvalidCacheIndex} {
    if ( !states.empty() ) {
      if ( front( states )->z() < back( states )->z() ) {
        m_states.insert( m_states.begin(), states.begin(), states.end() );
      } else {
        m_states.insert( m_states.begin(), states.rbegin(), states.rend() );
      }
    }
    // check states and initialize cache
    init( magfieldsvc );
  }

  TrackTraj::TrackTraj( span<const LHCb::State> states, LHCb::Tag::State::AssumeSorted_tag,
                        const DeMagnet*         magfieldsvc )
      : ZTrajectory<double>{}, m_bfield{0, 0, 0}, m_cachedindex{InvalidCacheIndex} {
    if ( !states.empty() ) {
      if ( front( states ).z() < back( states ).z() ) {
        append_range( m_states, states | take_address );
      } else {
        append_range( m_states, states | ranges::views::reverse | take_address );
      }
    }
    // check states and initialize cache
    init( magfieldsvc );
  }

  void TrackTraj::init( const DeMagnet* magfieldsvc ) {
    // add this points the vector of states must be sorted!  remove
    // any states with equal z
    m_states.erase( std::unique( m_states.begin(), m_states.end(), equalStateZ ), m_states.end() );

    // test that there are sufficient states left
    if ( m_states.empty() )
      throw GaudiException( "TrackTraj: not enough states for interpolation!", "TrackTraj::TrackTraj",
                            StatusCode::FAILURE );

    // set the range of the trajectory
    this->setRange( m_states.front()->z(), m_states.back()->z() );

    // set the field at the first state
    if ( magfieldsvc ) m_bfield = magfieldsvc->fieldVector( m_states.front()->position() );

    // invalidate the cache
    invalidateCache();
  }

  void TrackTraj::updatecache( double z ) const {
    // m_cachedindex==0: before first state
    // m_cachedindex==[1,...,numstates-1] --> between states
    // m_cachedindex==numstates --> after last state
    // m_cachedindex==INVALIDCACHEINDEX --> cache is not valid
    bool cacheisvalid =
        ( m_cachedindex != InvalidCacheIndex ) &&
        ( ( m_cachedindex == 0 && z <= m_states.front()->z() ) ||
          ( m_cachedindex == m_states.size() && z >= m_states.back()->z() ) ||
          ( m_cachedindex != 0 && m_states[m_cachedindex - 1]->z() <= z && z < m_states[m_cachedindex]->z() ) );

    if ( !cacheisvalid ) {
      if ( z <= m_states.front()->z() ) {
        m_cachedindex = 0;
        m_cachedinterpolation.init( *m_states.front(), m_bfield );
      } else if ( z >= m_states.back()->z() ) {
        m_cachedindex = m_states.size();
        m_cachedinterpolation.init( *m_states.back(), Gaudi::XYZVector( 0, 0, 0 ) );
      } else {
        m_cachedindex = 1;
        while ( m_cachedindex < m_states.size() - 1 && m_states[m_cachedindex]->z() <= z ) ++m_cachedindex;
        m_cachedinterpolation.init( *m_states[m_cachedindex - 1], *m_states[m_cachedindex] );
      }
    }
  }

  // Copied from Gerhard
  class ArcLengthComputer {
  private:
    typedef LHCb::TrackTraj    param_t;
    gsl_integration_workspace* m_workspace;
    size_t                     m_limit;
    static double GSLgluefun( double z, void* x ) { return static_cast<param_t*>( x )->dArclengthDMu( z ); }

  public:
    ArcLengthComputer() : m_limit( 1000 ) { m_workspace = gsl_integration_workspace_alloc( m_limit ); }
    ~ArcLengthComputer() { gsl_integration_workspace_free( m_workspace ); }
    double compute( const TrackTraj& traj, double z1, double z2 ) const {
      const double epsAbs = 1 * Gaudi::Units::cm;
      const double epsRel = 1e-3;
      gsl_function f;
      f.function = &ArcLengthComputer::GSLgluefun;
      f.params   = const_cast<param_t*>( &traj );
      double result, error;
      // size_t neval;
      // gsl_integration_qng(&f, z1, z2, epsAbs, epsRel,&result, &error,&neval);
      // std::cout << "ArcLengthComputer::compute: " << error << " " << neval << std::endl ;
      const int key = 2;
      gsl_integration_qag( &f, z1, z2, epsAbs, epsRel, m_limit, key, m_workspace, &result, &error );
      // std::cout << "ArcLengthComputer::compute: " << error << std::endl ;
      return result;
    }
  };

  double TrackTraj::arclength( double z1, double z2 ) const {
    static const ArcLengthComputer computer{};
    return computer.compute( *this, z1, z2 );
  }

  std::vector<StateVector> TrackTraj::refStateVectors() const {
    std::vector<StateVector> states;
    states.reserve( m_states.size() );
    for ( const auto& s : m_states ) states.emplace_back( s->stateVector(), s->z() );
    return states;
  }

  double TrackTraj::distTo1stError( double z, double tolerance, int pathDirection ) const {
    // WH: if timing is ever an issue, we should probably just return 'a' and not care about boundaries.
    updatecache( z );
    const double a           = m_cachedinterpolation.distTo1stError( z, tolerance, pathDirection );
    const bool   extrapolate = z <= m_states.front()->z() || z >= m_states.back()->z();
    // add fudge factor to make sure we step across boundaries
    const double fudgefactor = 1.01;
    return fudgefactor * ( extrapolate ? a
                                       : std::min( a, pathDirection > 0 ? m_states[m_cachedindex]->z() - z
                                                                        : z - m_states[m_cachedindex - 1]->z() ) );
  }

  double TrackTraj::distTo2ndError( double z, double tolerance, int pathDirection ) const {
    // WH: if timing is ever an issue, we should probably just return 'a' and not care about boundaries.
    updatecache( z );
    const double a           = m_cachedinterpolation.distTo2ndError( z, tolerance, pathDirection );
    const bool   extrapolate = z <= m_states.front()->z() || z >= m_states.back()->z();
    // add fudge factor to make sure we step across boundaries
    const double fudgefactor = 1.01;
    return fudgefactor * ( extrapolate ? a
                                       : std::min( a, pathDirection > 0 ? m_states[m_cachedindex]->z() - z
                                                                        : z - m_states[m_cachedindex - 1]->z() ) );
  }

} // namespace LHCb
