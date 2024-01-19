/*****************************************************************************\
* (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#ifndef TrackEvent_TrackFitResult_H
#define TrackEvent_TrackFitResult_H 1

// Include files
#include "Event/ChiSquare.h"
#include "Event/FitNode.h"
#include "Event/ITrackFitResult.h"
#include "Event/Measurement.h"
#include "Kernel/LHCbID.h"
#include "Kernel/STLExtensions.h"

namespace LHCb {

  /** @class TrackFitResult TrackFitResult.h
   *
   * TrackFitResult holds transient information from the track fit
   *
   * @author Jose Hernando, Eduardo Rodrigues, Wouter Hulsbergen, Daniel Campora
   * created Sun Jan 20 00:21:21 2019
   *
   */

  class TrackFitResult : public ITrackFitResult {
  public:
    /// Container for LHCb::Measurements on track
    using MeasurementContainer = std::vector<LHCb::Measurement>;
    /// define the node type (for use in templates)
    using NodeType = LHCb::FitNode;
    /// Container for LHCb::FitNodes on track
    using NodeContainer = std::vector<NodeType*>;

    /// Constructor
    TrackFitResult() = default;

    /// Destructor
    ~TrackFitResult() override;

    /// Track chi square
    virtual const LHCb::ChiSquare& chi2() const { return m_chi2; }

    /// VELO segment chi square
    virtual const LHCb::ChiSquare& chi2Velo() const { return m_chi2Velo; }

    /// Upstream segment chi square
    virtual const LHCb::ChiSquare& chi2Upstream() const { return m_chi2Upstream; }

    /// Velo-TT-T segment chi square
    virtual const LHCb::ChiSquare& chi2Long() const { return m_chi2Long; }

    /// Muon segment chi square
    virtual const LHCb::ChiSquare& chi2Muon() const { return m_chi2Muon; }

    /// Downstream segment chi square
    virtual const LHCb::ChiSquare& chi2Downstream() const { return m_chi2Downstream; }

    /// Upstream versus downstream segment chi square
    virtual ChiSquare chi2Match() const;

    /// Clone the track (you take ownership of the pointer)
    std::unique_ptr<LHCb::ITrackFitResult> clone() const override;

    /// Get the number of outliers on the track
    unsigned int nOutliers() const;

    /// Get the number of active (non-outlier) measurements on the track
    unsigned int nActiveMeasurements() const;

    /// set the measurements of this track. Track takes ownership.
    void setMeasurements( std::vector<LHCb::Measurement>&& measurements );

    /// Clear all measurements on track
    void clearMeasurements();

    /// Clear all nodes on track
    void clearNodes();

    /// Retrieve the number of Measurements on the track
    unsigned int nMeasurements() const;

    /// Retrieve the number of measurements of a certain type
    template <typename Type>
    unsigned int nMeasurements() const {
      return std::count_if( m_measurements.begin(), m_measurements.end(),
                            []( const Measurement& m ) { return m.is<Type>(); } );
    }

    /// Retrieve the number of active measurements of a certain type
    template <typename Type>
    unsigned int nActiveMeasurements() const {
      return std::count_if( nodes().begin(), nodes().end(), [&]( const FitNode* node ) {
        return node->type() == LHCb::FitNode::Type::HitOnTrack && node->measurement().is<Type>();
      } );
    }

    /// Return pointer to measurement on the track corresponding to the input LHCbID if present. If not, return nullptr
    /// (ONLY for tracking code, not for general use.)
    const LHCb::Measurement* measurement( const LHCb::LHCbID& value ) const;

    /// Clear the track before re-use
    virtual void reset();

    /// Retrieve const  Number of iterations in track fit
    int nIter() const { return m_nIter; }

    /// Update  Number of iterations in track fit
    void setNIter( int value );

    /// Retrieve const  Momentum used for computing material corrections
    double pScatter() const { return m_pScatter; }

    /// Update  Momentum used for computing material corrections
    void setPScatter( double value );

    /// Update  Track chi square
    void setChi2( const LHCb::ChiSquare& value );

    /// Update  VELO segment chi square
    void setChi2Velo( const LHCb::ChiSquare& value );

    /// Update  Upstream segment chi square
    void setChi2Upstream( const LHCb::ChiSquare& value );

    /// Update  Velo-TT-T segment chi square
    void setChi2Long( const LHCb::ChiSquare& value );

    /// Update  Muon segment chi square
    void setChi2Muon( const LHCb::ChiSquare& value );

    /// Update  Downstream segment chi square
    void setChi2Downstream( const LHCb::ChiSquare& value );

    /// Retrieve const  Container of Measurements
    const std::vector<LHCb::Measurement>& measurements() const { return m_measurements; }

    /// Retrieve const  Container of Nodes
    const std::vector<LHCb::FitNode*>& nodes() const { return m_nodes; }

    /// Retrieve  Container of Nodes
    std::vector<LHCb::FitNode*>& nodes() { return m_nodes; }

  protected:
    /// Copy constructor (hidden). Use clone method instead.
    TrackFitResult( const LHCb::TrackFitResult& rhs );

    LHCb::ChiSquare m_chi2;           ///< Track chi square
    LHCb::ChiSquare m_chi2Velo;       ///< VELO segment chi square
    LHCb::ChiSquare m_chi2Upstream;   ///< Upstream segment chi square
    LHCb::ChiSquare m_chi2Long;       ///< Velo-TT-T segment chi square
    LHCb::ChiSquare m_chi2Muon;       ///< Muon segment chi square
    LHCb::ChiSquare m_chi2Downstream; ///< Downstream segment chi square

  private:
    int                            m_nIter    = -1; ///< Number of iterations in track fit
    double                         m_pScatter = 0;  ///< Momentum used for computing material corrections
    std::vector<LHCb::Measurement> m_measurements;  ///< Container of Measurements
    std::vector<LHCb::FitNode*>    m_nodes;         ///< Container of Nodes

  }; // class TrackFitResult

  inline auto nodes( TrackFitResult const& fr ) {
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

// -----------------------------------------------------------------------------
// end of class
// -----------------------------------------------------------------------------

inline void LHCb::TrackFitResult::setNIter( int value ) { m_nIter = value; }

inline void LHCb::TrackFitResult::setPScatter( double value ) { m_pScatter = value; }

inline void LHCb::TrackFitResult::setChi2( const LHCb::ChiSquare& value ) { m_chi2 = value; }

inline void LHCb::TrackFitResult::setChi2Velo( const LHCb::ChiSquare& value ) { m_chi2Velo = value; }

inline void LHCb::TrackFitResult::setChi2Upstream( const LHCb::ChiSquare& value ) { m_chi2Upstream = value; }

inline void LHCb::TrackFitResult::setChi2Long( const LHCb::ChiSquare& value ) { m_chi2Long = value; }

inline void LHCb::TrackFitResult::setChi2Muon( const LHCb::ChiSquare& value ) { m_chi2Muon = value; }

inline void LHCb::TrackFitResult::setChi2Downstream( const LHCb::ChiSquare& value ) { m_chi2Downstream = value; }

inline LHCb::ChiSquare LHCb::TrackFitResult::chi2Match() const { return m_chi2 - m_chi2Upstream - m_chi2Downstream; }

inline unsigned int LHCb::TrackFitResult::nMeasurements() const { return m_measurements.size(); }

#endif /// TrackEvent_TrackFitResult_H
