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
#ifndef TRACKFITEVENT_KALMANFITRESULT_H
#define TRACKFITEVENT_KALMANFITRESULT_H

#include "Event/ChiSquare.h"
#include "Event/Measurement.h"
#include "Event/TrackFitResult.h"
#include "Event/TrackMaterialIntersection.h"
#include "Event/TrackTypes.h"

namespace LHCb {
  class FitNode;

  class KalmanFitResult final : public TrackFitResult {
  public:
    using Intersections = std::vector<LHCb::TrackMaterialIntersection>;

    enum algoType { Predict, Filter, Smooth, ComputeResidual, WeightedAverage };
    enum errorType { Initialization, MatrixInversion, AlgError, Other };
    enum directionType { Forward, Backward, BiDirection };

    // default constructor. do nothing.
    KalmanFitResult() = default;

    // copy constructor
    KalmanFitResult( const KalmanFitResult& rhs );

    // copy from TrackFitResult
    KalmanFitResult( const TrackFitResult& rhs );

    // clone
    std::unique_ptr<ITrackFitResult> clone() const override;

    // set links between nodes and to parent
    void establishNodeLinks();

    // reset the kalman filter in all the nodes to 'initialized'
    void resetFilterStatus();

    // reset the cache
    void resetCache() { m_chi2CacheValid = false; }

    // get the seed covariance
    const Gaudi::TrackSymMatrix& seedCovariance() const { return m_seedCovariance; }

    // set the seed covariance
    void setSeedCovariance( const Gaudi::TrackSymMatrix& rhs ) { m_seedCovariance = rhs; }

    // return the number of track parameters
    int nTrackParameters() const { return m_nTrackParameters; }

    // set the number of track parameters and reset the cache
    void setNTrackParameters( int n ) {
      m_nTrackParameters = n;
      m_chi2CacheValid   = false;
    }

    KalmanFitResult& unConst() const { return const_cast<KalmanFitResult&>( *this ); }

    // return (chisq,dof) for this track
    const ChiSquare& chi2() const override {
      if ( !m_chi2CacheValid ) unConst().computeChiSquares();
      return m_chi2;
    }

    // return (chisq,dof) for the velo part of this track
    const ChiSquare& chi2Velo() const override {
      if ( !m_chi2CacheValid ) unConst().computeChiSquares();
      return m_chi2Velo;
    }

    // return (chisq,dof) for the segment downstream of the magnet (T + Muon)
    const ChiSquare& chi2Downstream() const override {
      if ( !m_chi2CacheValid ) unConst().computeChiSquares();
      return m_chi2Downstream;
    }

    // return (chisq,dof) for the segment upstream of the magnet
    const ChiSquare& chi2Upstream() const override {
      if ( !m_chi2CacheValid ) unConst().computeChiSquares();
      return m_chi2Upstream;
    }

    // return (chisq,dof) for upstream versus downstream segment
    ChiSquare chi2Match() const override {
      if ( !m_chi2CacheValid ) unConst().computeChiSquares();
      return chi2() - m_chi2Upstream - m_chi2Downstream;
    }

    // return (chisq,dof) for the velo-TT-T segment, so everything excluding muon
    const ChiSquare& chi2Long() const override {
      if ( !m_chi2CacheValid ) unConst().computeChiSquares();
      return m_chi2Long;
    }

    // return (chisq,dof) for the muon segment
    const ChiSquare& chi2Muon() const override {
      if ( !m_chi2CacheValid ) unConst().computeChiSquares();
      return m_chi2Muon;
    }

    // return (chisq,dof) for the muon - T match
    ChiSquare chi2MuonTMatch() const {
      if ( !m_chi2CacheValid ) unConst().computeChiSquares();
      return chi2() - m_chi2Long - m_chi2Muon;
    }

    // return (chisq,dof) of the velo-T match without TT hits. this is _not_ fast.
    ChiSquare chi2VeloTMatch() const;

    // return (chisq,dof) for the contribution of TT (or UT) hits
    ChiSquare chi2TTHits() const;

    // set the error flag
    void setErrorFlag( unsigned short extrainfo, unsigned short algnum, unsigned short errnum );

    // check if there is an error
    bool inError() const;

    // check the type of error
    std::string getError() const;

    // check if the fit is bidirectionnal or classical(false)
    bool biDirectionnalSmoother() const { return m_bidirectionalSmoother; };

    // set the type of fit: bidirectionnal or classical(false)
    void setBiDirectionnalSmoother( bool bidir ) { m_bidirectionalSmoother = bidir; };

    // the number of hits for which OT drifttime is used
    unsigned int nActiveOTTimes() const { return 0; };

    // accessors to intersections
    const auto& intersections() const { return m_intersections; }
    auto&       intersections() { return m_intersections; }

  public:
    void computeChiSquares();

  private:
    enum errorBits {
      typeBits   = 0, // error type bit position
      algBits    = 2, // function bit position
      dirBits    = 5, // direction bit position
      globalBits = 7  // global error bit position
    };
    enum errorMasks {
      typeMask   = 0x03, // error type mask
      algMask    = 0x1C, // function mask
      dirMask    = 0x60, // direction mask
      globalMask = 0x80  // direction mask
    };

  private:
    Gaudi::TrackSymMatrix                  m_seedCovariance;
    int                                    m_nTrackParameters      = 5;
    mutable bool                           m_chi2CacheValid        = false;
    unsigned short                         m_errorFlag             = 0x00;
    bool                                   m_bidirectionalSmoother = true;
    std::vector<TrackMaterialIntersection> m_intersections;
  };

  inline auto nodes( KalmanFitResult const& fr ) {
    struct Range {
      std::vector<LHCb::FitNode*> const* parent;
      struct Iterator {
        std::vector<LHCb::FitNode*>::const_iterator iter;

        using difference_type [[maybe_unused]]   = std::ptrdiff_t;
        using value_type [[maybe_unused]]        = const LHCb::FitNode;
        using pointer [[maybe_unused]]           = const LHCb::FitNode*;
        using reference [[maybe_unused]]         = const LHCb::FitNode&;
        using iterator_category [[maybe_unused]] = std::random_access_iterator_tag;

        Iterator& operator++() {
          ++iter;
          return *this;
        }
        Iterator& operator--() {
          --iter;
          return *this;
        }
        // Dereference both the iterator AND the pointer to FitNode.
        // Key part of making loops over LHCb::FitNodes and PrFitNodes
        // look the same. Used when templating algorithms that use fit nodes.
        const LHCb::FitNode& operator*() const { return **iter; }
        bool                 operator!=( Iterator const& rhs ) { return iter != rhs.iter; }
      };
      Iterator begin() const { return {parent->begin()}; }
      Iterator end() const { return {parent->end()}; }
      auto     front() const { return *( parent->front() ); }
      auto     back() const { return *( parent->back() ); }
      auto     size() const { return parent->size(); }
    };
    return Range{&fr.nodes()};
  }

} // namespace LHCb

#endif
