/*****************************************************************************\
# (c) Copyright 2000-2021 CERN for the benefit of the LHCb Collaboration      #
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once
#include <chrono>
namespace LHCb::TrackKernel {
  class Timer {
    typedef std::chrono::high_resolution_clock high_resolution_clock;
    typedef std::chrono::milliseconds          nanoseconds;

  public:
    explicit Timer() : m_total( 0 ), m_max( 0 ), m_numcalls( 0 ), m_start( high_resolution_clock::now() ) {}
    void   start() { m_start = high_resolution_clock::now(); }
    double stop() {
      auto   diff    = high_resolution_clock::now() - m_start;
      double elapsed = std::chrono::duration<double, std::nano>( diff ).count();
      m_total += elapsed;
      if ( elapsed > m_max ) m_max = elapsed;
      ++m_numcalls;
      return m_total;
    }
    double total() const { return m_total; }
    double average() const { return m_numcalls > 0 ? m_total / m_numcalls : 0; }
    double maximum() const { return m_max; }
    size_t numcalls() const { return m_numcalls; }
    void   reset() {
      m_total = m_max = 0;
      m_numcalls      = 0;
    }

  private:
    double                            m_total;
    double                            m_max;
    size_t                            m_numcalls;
    high_resolution_clock::time_point m_start;
  };
} // namespace LHCb::TrackKernel
