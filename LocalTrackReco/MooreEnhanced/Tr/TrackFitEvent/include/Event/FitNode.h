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
#ifndef TRACKFITEVENT_FITNODE_H
#define TRACKFITEVENT_FITNODE_H 1

#ifdef USE_DD4HEP
#  include "Core/TaggedBool.h"
#endif

// from TrackEvent
#include "Event/ChiSquare.h"
#include "Event/Measurement.h"
#include "Event/State.h"
#include "GaudiKernel/compose.h"
#include "LHCbMath/ValueWithError.h"

#include "GaudiKernel/SerializeSTL.h"
// From LHCbMath
#include "LHCbMath/MatrixManip.h"
#include "LHCbMath/Similarity.h"

namespace LHCb {

  /** @class FitNode FitNode.h Event/FitNode.h
   *
   *  This File contains the declaration of the FitNode.
   *
   *  A FitNode is a basket of objects at a given z position.
   *  The information inside the FitNode has to be sufficient
   *  to allow smoothing and refitting.
   *
   *  At the moment a FitNode contains or allows you to access
   *  info on the the (kth) measurement,
   *  transport from k --> k + 1 , predicted state at k+1
   *  (predicted from filter step)  and the best state at k
   *  (notation note filter procedes from k -> k + 1 -> k + 2 ......)
   *
   *  @author Eduardo Rodrigues
   *  @date   2005-04-15 (adaptations to the new track event model)
   *
   *  @author Matthew Needham
   *  @date   19-08-1999
   */

  class KalmanFitResult;

  namespace Enum::FitNode {
    meta_enum_class( Type, unsigned char, Unknown, HitOnTrack, Outlier, Reference )
  }

  namespace Enum::nDim {
    //  meta_enum_class( Type, unsigned char, Unknown, one=1, two=2 )
    enum Type { one = 1, two = 2 };
  } // namespace Enum::nDim

  class FitNode final {
  public:
    /// enumerator for the type of Node
    using Type    = Enum::FitNode::Type;
    using TypeDim = Enum::nDim::Type;

#ifdef USE_DD4HEP
    using AssumeHitIsUsed = LHCb::Detector::tagged_bool_ns::tagged_bool<struct AssumeHitIsUsed_tag>;
#else
    using AssumeHitIsUsed = tagged_bool<struct AssumeHitIsUsed_tag>;
#endif

    // important note: for the Forward fit, smoothed means
    // 'classical'. for the backward fit, it means 'bidirectional'.
    enum FilterStatus { Uninitialized, Initialized, Predicted, Filtered, Smoothed };
    enum CachedBool { False = 0, True = 1, Unknown = 2 };
    enum Direction { Forward = 0, Backward = 1 };

    /// Retrieve the local chi^2
    LHCb::ChiSquare chi2() const;

    /// retrieve the state, overloading the inline function in Node
    const State& state() const;

    /// Update the state vector
    void setState( const Gaudi::TrackVector& stateVector, const Gaudi::TrackSymMatrix& stateCovariance,
                   const double& z );

    /// Update the reference vector
    void setRefVector( const Gaudi::TrackVector& refVector );

    /// Update the reference vector
    void setRefVector( const LHCb::StateVector& refVector );

    /// Retrieve const  the reference to the measurement
    const LHCb::Measurement& measurement() const { return *m_measurement; }

    /// Set the measurement
    void setMeasurement( const LHCb::Measurement& meas ) { m_measurement = &meas; }

    /// Return true if this Node has a valid pointer to measurement
    bool hasMeasurement() const { return m_measurement != nullptr; }

    /// z position of Node
    double z() const { return m_refVector.z(); }

    /// Position of the State on the Node
    Gaudi::XYZPoint position() const { return m_state.position(); }

    /// Set the location of the state in this node.
    void setLocation( LHCb::State::Location& location ) { m_state.setLocation( location ); }

    /// Location of the state in this node
    auto location() const { return m_state.location(); }

    /// Retrieve const  type of node
    const Type& type() const { return m_type; }

    /// Update  type of node
    void setType( const Type& value ) { m_type = value; }

    /// Check type to be hit on track - for common interface with PrFitNode
    bool isHitOnTrack() const { return m_type == Type::HitOnTrack; }

    /// Check type to be outlier - for common interface with PrFitNode
    bool isOutlier() const { return m_type == Type::Outlier; }

    /// Check type of measurement - for common interface with PrFitNode
    bool isVP() const {
      return m_measurement &&
             ( m_measurement->is<LHCb::Measurement::VP>() || m_measurement->is<LHCb::Measurement::VP2D>() );
    }
    bool isUT() const { return m_measurement && m_measurement->is<LHCb::Measurement::UT>(); }
    bool isFT() const { return m_measurement && m_measurement->is<LHCb::Measurement::FT>(); }
    bool isMuon() const { return m_measurement && m_measurement->is<LHCb::Measurement::Muon>(); }

    /// Update  state
    void setState( const LHCb::State& value ) { m_state = value; }

    /// Retrieve const  the reference vector
    const LHCb::StateVector& refVector() const { return m_refVector; }

    /// Retrieve const  flag for the reference vector
    bool refIsSet() const { return m_refIsSet; }

    /// Retrieve const  Signed doca (of ref-traj). For ST/velo this is equal to minus (ref)residual
    double doca() const { return m_doca; }

    /// Update  Signed doca (of ref-traj). For ST/velo this is equal to minus (ref)residual
    void setDoca( double value ) { m_doca = value; }

    /// Retrieve const  Unit vector perpendicular to state and measurement
    const Gaudi::XYZVector& pocaVector() const { return m_pocaVector; }

    /// Update  Unit vector perpendicular to state and measurement
    void setPocaVector( const Gaudi::XYZVector& value ) { m_pocaVector = value; }

    /// retrieve transport matrix
    const Gaudi::TrackMatrix& transportMatrix() const { return m_transportMatrix; }

    /// retrieve invert transport matrix
    const Gaudi::TrackMatrix& invertTransportMatrix() const { return m_invertTransportMatrix; }

    /// retrieve transport vector
    const Gaudi::TrackVector& transportVector() const { return m_transportVector; }

    /// retrieve noise matrix
    const Gaudi::TrackSymMatrix& noiseMatrix() const { return m_noiseMatrix; }

    /// retrieve noise matrix
    Gaudi::TrackSymMatrix& noiseMatrix() { return m_noiseMatrix; }

    /// set transport matrix
    void setTransportMatrix( const Gaudi::TrackMatrix& transportMatrix );

    /// set transport vector
    void setTransportVector( const Gaudi::TrackVector& transportVector ) { m_transportVector = transportVector; }

    /// set noise matrix
    void setNoiseMatrix( const Gaudi::TrackSymMatrix& noiseMatrix ) { m_noiseMatrix = noiseMatrix; }

    /// set the seed matrix
    void setSeedCovariance( const Gaudi::TrackSymMatrix& cov ) {
      m_predictedState[Forward].covariance()  = cov;
      m_predictedState[Backward].covariance() = cov;
    }

    /// retrieve state predicted by the kalman filter step
    const State& predictedStateForward() const { return predictedState( Forward ); }

    /// retrieve predicted state from backward filter
    const State& predictedStateBackward() const { return predictedState( Backward ); }

    /// retrieve state filtered by the kalman filter step
    const State& filteredStateForward() const { return filteredState( Forward ); }

    /// retrieve state filtered by the kalman filter step
    const State& filteredStateBackward() const { return filteredState( Backward ); }

    /// retrieve chisq contribution in upstream filter
    const LHCb::ChiSquare& deltaChi2Forward() const {
      filteredStateForward();
      return m_deltaChi2[Forward];
    }

    /// retrieve chisq contribution in downstream filter
    const LHCb::ChiSquare& deltaChi2Backward() const {
      filteredStateBackward();
      return m_deltaChi2[Backward];
    }

    /// retrieve the total chi2 of the filter including this node
    const LHCb::ChiSquare& totalChi2( int direction ) const {
      filteredState( direction );
      return m_totalChi2[direction % 2];
    }

    /// set the delta-energy
    void setDeltaEnergy( double e ) { m_deltaEnergy = e; }

    /// get the delta-energy
    double deltaEnergy() const { return m_deltaEnergy; }

    // get the smoother gain matrix
    const Gaudi::TrackMatrix& smootherGainMatrix() const { return m_smootherGainMatrix; }

    /// get the filter status (only useful for debugging)
    FilterStatus filterStatus( int direction ) const { return m_filterStatus[direction]; }

    /// return whether or not this node has active nodes upstream
    bool hasInfoUpstream( int direction ) const;

    /// Deactivate this node (outlier)
    void deactivateMeasurement( bool deactivate = true );

    /// Get the index of this node. For debugging only.
    int index() const;

    /// set previous node
    void setPreviousNode( FitNode* previousNode ) {
      m_prevNode = previousNode;
      if ( m_prevNode ) m_prevNode->m_nextNode = this;
    }

    /// Unlink this node
    void unLink() {
      m_prevNode = m_nextNode = nullptr;
      m_parent                = nullptr;
    }

    /// set the parent
    void setParent( KalmanFitResult* p ) { m_parent = p; }

    /// get the parent
    KalmanFitResult* getParent() { return m_parent; }

    void updateResidual( const LHCb::State& smoothedState );

  private:
    Type                     m_type;                  ///< type of node
    LHCb::State              m_state;                 ///< state
    LHCb::StateVector        m_refVector;             ///< the reference vector
    bool                     m_refIsSet    = false;   ///< flag for the reference vector
    const LHCb::Measurement* m_measurement = nullptr; ///< pointer to the measurement (not owner)
    double           m_doca = 0;   ///< Signed doca (of ref-traj). For ST/velo this is equal to minus (ref)residual
    Gaudi::XYZVector m_pocaVector; ///< Unit vector perpendicular to state and measurement

    // ! check that the contents of the cov matrix are fine
    bool isPositiveDiagonal( const Gaudi::TrackSymMatrix& mat ) const;

  public:
    const FitNode* prevNode( int direction ) const { return direction == Forward ? m_prevNode : m_nextNode; }
    const FitNode* nextNode( int direction ) const { return direction == Forward ? m_nextNode : m_prevNode; }

    /// retrieve the predicted state
    const LHCb::State& predictedState( int direction ) const {
      if ( m_filterStatus[direction] < Predicted ) unConst().computePredictedState( direction );
      return m_predictedState[direction];
    }

    /// retrieve the filtered state
    const LHCb::State& filteredState( int direction ) const {
      if ( m_filterStatus[direction] < Filtered ) unConst().computeFilteredState( direction );
      return m_filteredState[direction];
    }

    /// retrieve the bismoothed state
    const LHCb::State& biSmoothedState() const {
      if ( m_filterStatus[Backward] < Smoothed ) unConst().computeBiSmoothedState();
      return m_state;
    }

    /// retrieve the classically smoothed state
    const LHCb::State& classicalSmoothedState() const {
      if ( m_filterStatus[Forward] < Smoothed ) unConst().computeClassicalSmoothedState();
      return m_classicalSmoothedState;
    }

    /// This is used from the projectors (or from any set method?)
    void resetFilterStatus( FilterStatus s = Initialized ) {
      resetFilterStatus( Forward, s );
      resetFilterStatus( Backward, s );
    }

  private:
    void computePredictedState( int direction );
    void computeFilteredState( int direction );
    void computeBiSmoothedState();
    void computeClassicalSmoothedState();
    void smooth() const;

  private:
    /// unconst this node
    FitNode& unConst() const { return const_cast<FitNode&>( *this ); }

    /// reset the cache for the previous function
    void resetHasInfoUpstream( int direction );

    /// reset the filter status
    void resetFilterStatus( int direction, FilterStatus s = Initialized );

  private:
    Gaudi::TrackMatrix m_transportMatrix;       ///< transport matrix for propagation from previous node to this one
    Gaudi::TrackMatrix m_invertTransportMatrix; ///< transport matrix for propagation from this node to the previous one
    Gaudi::TrackVector m_transportVector;       ///< transport vector for propagation from previous node to this one
    Gaudi::TrackSymMatrix m_noiseMatrix;        ///< noise in propagation from previous node to this one
    double                m_deltaEnergy = 0;    ///< change in energy in propagation from previous node to this one
    FilterStatus          m_filterStatus[2];    ///< Status of the Node in the fit process
    CachedBool            m_hasInfoUpstream[2]; ///< Are the nodes with active measurement upstream of this node?
    State                 m_predictedState[2];  ///< predicted state of forward/backward filter
    State                 m_filteredState[2];   ///< filtered state of forward/backward filter
    LHCb::State           m_classicalSmoothedState;
    LHCb::ChiSquare       m_deltaChi2[2];       ///< chisq contribution in forward filter
    LHCb::ChiSquare       m_totalChi2[2];       ///< total chi2 after this filterstep
    Gaudi::TrackMatrix    m_smootherGainMatrix; ///< smoother gain matrix (smoothedfit only)
    FitNode*              m_prevNode = nullptr; ///< Previous Node
    FitNode*              m_nextNode = nullptr; ///< Next Node
    KalmanFitResult*      m_parent   = nullptr; ///< Owner

  public:
    // Measurement dimension specific info and functionality
    template <TypeDim dim, typename Float = double>
    struct DimInfos final {
      // DimInfos( FitNode* _owner ) : owner{_owner} {}

      using SVectorT          = ROOT::Math::SVector<Float, dim>;
      using SMatrixT          = ROOT::Math::SMatrix<Float, dim, dim>;
      using ProjectionMatrixT = ROOT::Math::SMatrix<Float, dim, 5>;

      static constexpr inline TypeDim typedim = dim;
      // FitNode*                        owner   = nullptr;
      // void setOwner( FitNode* fnode ) { owner = fnode; }

      SVectorT          m_residual;         ///< the residual value
      SMatrixT          m_errResidual;      ///< the residual error
      SVectorT          m_errMeasure;       ///< the measure error
      SVectorT          m_refResidual;      ///< residual of the reference
      ProjectionMatrixT m_projectionMatrix; ///< the projection matrix

      /// retrieve the residual, overloading the function in Node
      SVectorT residual( const FitNode& owner ) const {
        owner.smooth();
        return m_residual;
      };

      /// retrieve the residual, overloading the function in Node
      SMatrixT errResidual( const FitNode& owner ) const {
        owner.smooth();
        return m_errResidual;
      };

      SMatrixT errResidual2( const FitNode& owner ) const {
        auto     errRes = errResidual( owner );
        SMatrixT returnvar; // initialized with 0
        std::transform( errRes.begin(), errRes.end(), returnvar.begin(), []( Float x ) { return x * x; } );
        return returnvar;
      }

      SMatrixT errResidualinv( const FitNode& owner ) const {
        auto     errRes = errResidual( owner );
        int      ifail  = 1;
        SMatrixT Rinv   = errRes.Inverse( ifail );
        if ( ifail == 1 ) {
          // inversion fails if elements are 0. This happens when type is not HitOnTrack.
          // previously caught when dividing by 1D value: check if not zero (otherwise do nothing).
          // Residuals are not interesting anyway for a non-measurement FitNode -
          //  responsibility of the caller to realize this in current implementation.
          // throw std::logic_error( "Inversion of errResidual failed" );
          SMatrixT zero; // initialized as 0.
          return zero;
        }
        return Rinv;
      }

      SMatrixT errResidual2inv( const FitNode& owner ) const {
        auto     errRes2 = errResidual2( owner );
        int      ifail   = 1;
        SMatrixT Rinv    = errRes2.Inverse( ifail );
        if ( ifail == 1 ) {
          // throw std::logic_error( "Inversion of errResidual2 failed" );
          SMatrixT zero; // initialized as 0.
          return zero;
        }
        return Rinv;
      }

      /// Retrieve const  the measure error
      SVectorT errMeasure() const { return m_errMeasure; }

      /// Return the measure error squared
      SVectorT errMeasure2() const {
        SVectorT returnval;
        std::transform( m_errMeasure.begin(), m_errMeasure.end(), returnval.begin(), []( Float x ) { return x * x; } );
        return returnval;
      }

      /// Retrieve const  the projection matrix
      const ProjectionMatrixT& projectionMatrix() const { return m_projectionMatrix; }

      /// retrieve the residual of the reference
      SVectorT refResidual() const { return m_refResidual; }

      void updateProjection( FitNode& owner, ROOT::Math::SMatrix<Float, dim, 5>& H,
                             ROOT::Math::SVector<Float, dim> refresidual, ROOT::Math::SVector<Float, dim> errmeasure ) {
        m_projectionMatrix = H;
        m_refResidual      = refresidual;
        m_errMeasure       = errmeasure;
        for ( unsigned int i = 0; i < dim; i++ ) {
          m_residual[i] = 0;
          for ( unsigned int j = 0; j < dim; j++ ) { m_errResidual[i][j] = 0; }
        }
        owner.resetFilterStatus( Predicted );
      }

      SVectorT unbiasedResidual( const FitNode& owner ) const {
        if ( owner.type() == Type::HitOnTrack ) {
          SVectorT returnvar; // initialized with 0
          auto     errResidual2 = this->errResidual2( owner );
          // 1D: residual * sqrt(errMeasure2 / errResidual2)
          for ( int i = 0; i < dim; i++ ) {
            double f = 0;
            // sum up all (off-diagonal) errors due to HCH
            for ( int j = 0; j < dim; j++ ) { f += errResidual2[i][j]; }
            if ( f == 0 ) {
              returnvar[i] = 0;
            } else {
              returnvar[i] = m_residual[i] * std::sqrt( m_errMeasure[i] * m_errMeasure[i] / f );
            }
          }
          return returnvar;
        } else {
          return m_residual;
        }
      }

      /// Retrieve error on unbiased residual
      SMatrixT errUnbiasedResidual( const FitNode& owner ) const {
        if ( owner.type() == Type::HitOnTrack ) {
          SMatrixT returnvar;
          auto     errMeasure2 = this->errMeasure2();
          // 1D: errMeasure2 / errResidual
          // 2D: as we grouped together effects from off-diagonal HCH
          //     for the unbiasedResidual, we do the same here.
          for ( int i = 0; i < dim; i++ ) {
            double f = 0;
            for ( int j = 0; j < dim; j++ ) { f += m_errResidual[i][j]; }
            if ( f == 0 ) {
              returnvar[i][i] = 0;
            } else {
              returnvar[i][i] = errMeasure2[i] / f;
            }
          }
          return returnvar;
        } else {
          return m_errResidual;
        }
      }

      /// Retrieve the unbiased smoothed state at this position
      LHCb::State unbiasedState( const FitNode& owner ) const {
        // we can redo this routine by smoothing the unfiltered states

        // This performs an inverse kalman filter step.
        // First calculate the gain matrix
        const auto& H       = m_projectionMatrix;
        const auto& biasedC = owner.state().covariance();

        const Gaudi::TrackSymMatrix unit = ROOT::Math::SMatrixIdentity();
        Gaudi::TrackSymMatrix       unbiasedC;
        LHCb::State                 returnstate;
        auto                        r    = m_residual;
        auto                        Rinv = this->errResidual2inv( owner );

        ROOT::Math::SMatrix<double, 5, dim> K = ( biasedC * ROOT::Math::Transpose( H ) ) * Rinv;
        // Gaudi::TrackVector unbiasedX = m_state.stateVector() - rmat * K;
        Gaudi::TrackVector unbiasedX = owner.state().stateVector();
        for ( int i = 0; i < dim; i++ ) { unbiasedX -= K.Col( i ) * r[i]; }
        ROOT::Math::AssignSym::Evaluate( unbiasedC, ( unit + K * H ) * biasedC );
        return LHCb::State{unbiasedX, unbiasedC, owner.z(), owner.state().location()};
      }

      struct computeResidualResult {
        SVectorT residual;
        SMatrixT errResidual;
      };
      void set( computeResidualResult& result ) {
        m_residual    = result.residual;
        m_errResidual = sqrt( result.errResidual );
      }
      bool isValid( computeResidualResult& result ) const {
        bool covcheck = true;
        for ( int i = 0; i < dim; i++ ) {
          if ( result.errResidual[i][i] <= 0 ) covcheck = false;
        }
        return covcheck;
      }

      /// update node residual using a smoothed state
      computeResidualResult computeResidual( const FitNode& owner, const LHCb::State& state,
                                             AssumeHitIsUsed biased ) const {
        SVectorT r;
        SMatrixT R;
        // We should put this inside the This->
        assert( owner.hasMeasurement() );
        const Gaudi::TrackVector& refX = owner.refVector().parameters();
        const auto                H    = m_projectionMatrix;
        auto                      HCH  = LHCb::Math::Similarity( H, state.covariance() ); // dim x dim
        double                    sign = biased ? -1 : 1;
        for ( int i = 0; i < dim; i++ ) {
          r[i] = m_refResidual[i] + ROOT::Math::Dot( H.Row( i ), ( refX - state.stateVector() ) );
          // As measurements only have errors in their measurement directions [dim],
          //  all off-diagonal terms in errRes come from the track covariance projection.
          R[i][i] = m_errMeasure[i] * m_errMeasure[i];
          for ( int j = 0; j < dim; j++ ) { R[i][j] += sign * HCH[i][j]; }
        }
        return computeResidualResult{r, R};
      }

    }; // end of dimension-specific struct

  public:
    // Use a variant for the various FitNode/Measurement dimensions
    using DimInfo = std::variant<DimInfos<TypeDim::one, double>, DimInfos<TypeDim::two, double>>;

  private:
    DimInfo m_dim{FitNode::DimInfos<Enum::nDim::Type::one>{}};

  public:
    template <typename R, typename... Fs>
    R visit_r( Fs&&... fs ) {
      R r = static_cast<R>( std::visit( Gaudi::overload( std::forward<Fs>( fs )... ), m_dim ) );
      return r;
    }
    template <typename R, typename... Fs>
    R visit_r( Fs&&... fs ) const {
      R r = static_cast<R>( std::visit( Gaudi::overload( std::forward<Fs>( fs )... ), m_dim ) );
      return r;
    }
    template <typename... Fs>
    decltype( auto ) visit( Fs&&... fs ) {
      return std::visit( Gaudi::overload( std::forward<Fs>( fs )... ), m_dim );
    }
    template <typename... Fs>
    decltype( auto ) visit( Fs&&... fs ) const {
      return std::visit( Gaudi::overload( std::forward<Fs>( fs )... ), m_dim );
    }
    // Same as above, but allow for an additional variant as first argument of visit
    template <typename Variant, typename... Fs>
    decltype( auto ) visit2( Variant&& v, Fs&&... fs ) {
      return std::visit( Gaudi::overload( std::forward<Fs>( fs )... ), m_dim, std::forward<Variant>( v ) );
    }
    template <typename Variant, typename... Fs>
    decltype( auto ) visit2( Variant&& v, Fs&&... fs ) const {
      return std::visit( Gaudi::overload( std::forward<Fs>( fs )... ), m_dim, std::forward<Variant>( v ) );
    }

    /// Get the dimension-specific information
    template <typename DimI>
    const DimI* getIf() const {
      return std::get_if<DimI>( &m_dim );
    }

    /// Check whether this FitNode 'is-a' specific dimension
    template <typename DimI>
    bool is() const {
      return std::holds_alternative<DimI>( m_dim );
    }

    /// Templated interfaces hiding DimInfos visits
    template <typename T = double>
    T residual() const;
    template <typename T = double>
    T errResidual() const;
    template <typename T = double>
    T errResidual2() const;
    template <typename T = double>
    T errMeasure2() const;
    template <typename T = double>
    T unbiasedResidual() const;

    /// Default constructor
    FitNode();

    /// Constructor from a z-position
    FitNode( double z );

    /// Constructor from a z position and a location
    FitNode( double zPos, LHCb::State::Location location = LHCb::State::Location::LocationUnknown );

    /// Constructor from a Measurement
    struct FitNode1DType {};
    struct FitNode2DType {};
    FitNode( const Measurement& aMeas, FitNode1DType );
    FitNode( const Measurement& aMeas, FitNode2DType );
    FitNode( const Measurement& aMeas, DimInfo dimi );

    /// Destructor
    ~FitNode();

  public:
    /**
     * @brief      Helper functions to print FitNodes in a compact form.
     */
    std::ostream& print( std::ostream& s, const Gaudi::TrackVector& v ) const {
      return GaudiUtils::details::ostream_joiner( s, v, " " );
    }

    template <auto N>
    std::ostream& print( std::ostream& s, const ROOT::Math::SMatrix<double, N, 5>& v ) const {
      return GaudiUtils::details::ostream_joiner( s, v, " " );
    }

    std::ostream& print( std::ostream& s, const Gaudi::TrackSymMatrix& v ) const {
      return GaudiUtils::details::ostream_joiner( s, v, " " );
    }

    std::ostream& print( std::ostream& s, const Gaudi::TrackMatrix& v ) const {
      return GaudiUtils::details::ostream_joiner( s, v, " " );
    }

  public:
    /**
     * @brief      Prints a FitNode in a compact form.
     */
    friend std::ostream& operator<<( std::ostream& s, const LHCb::FitNode& node );
  };

  // Functions bellow are used to ease templating of various algorithms on fit node type.
  // Each function has a twin in namespace LHCb::Pr::Tracks::Fit in PrFitNode.h
  // which returns the apropriate PrFitNode alternative. This allows making loops over fit nodes
  // look the same and thus helps with templating, while keeping the PrFitNode object "lightweight".
  inline const LHCb::LHCbID      id( const FitNode& node ) { return node.measurement().lhcbID(); }
  inline const LHCb::State&      state( const FitNode& node ) { return node.state(); }
  inline const Gaudi::XYZVector& pocaVector( const FitNode& node ) { return node.pocaVector(); }
  inline const LHCb::State&      predictedStateForward( const FitNode& node ) { return node.predictedStateForward(); }
  inline const LHCb::State&      predictedStateBackward( const FitNode& node ) { return node.predictedStateBackward(); }
  inline const LHCb::State&      filteredStateForward( const FitNode& node ) { return node.filteredStateForward(); }
  inline const LHCb::State&      filteredStateBackward( const FitNode& node ) { return node.filteredStateBackward(); }
  inline const Gaudi::TrackSymMatrix& smoothedStateCovariance( const FitNode& node, bool classical = true ) {
    if ( classical )
      return node.classicalSmoothedState().covariance();
    else
      return node.biSmoothedState().covariance();
  };
  inline const LHCb::State& smoothedState( const FitNode& node, bool classical = true ) {
    if ( classical )
      return node.classicalSmoothedState();
    else
      return node.biSmoothedState();
  }
} // namespace LHCb

inline LHCb::FitNode::FitNode( double z ) : m_type( Type::Reference ) { m_refVector.setZ( z ); }

inline void LHCb::FitNode::setState( const Gaudi::TrackVector&    stateVector,
                                     const Gaudi::TrackSymMatrix& stateCovariance, const double& z ) {
  m_state.setState( stateVector );
  m_state.setCovariance( stateCovariance );
  m_state.setZ( z );
}

inline void LHCb::FitNode::setRefVector( const Gaudi::TrackVector& refVector ) {
  m_refIsSet               = true;
  m_refVector.parameters() = refVector;
}

inline void LHCb::FitNode::setRefVector( const LHCb::StateVector& refVector ) {
  m_refIsSet  = true;
  m_refVector = refVector;
}

/// Templated interfaces hiding DimInfos visits
template <typename T>
T LHCb::FitNode::residual() const {
  return std::visit( [&]( auto& f ) -> T { return f.residual( *this ); }, m_dim );
}
template <typename T>
T LHCb::FitNode::errResidual() const {
  return std::visit( [&]( auto& f ) -> T { return f.errResidual( *this ); }, m_dim );
}
template <typename T>
T LHCb::FitNode::errResidual2() const {
  return std::visit( [&]( auto& f ) -> T { return f.errResidual2( *this ); }, m_dim );
}
template <typename T>
T LHCb::FitNode::errMeasure2() const {
  return std::visit( [&]( auto& f ) -> T { return f.errMeasure2(); }, m_dim );
}
template <typename T>
T LHCb::FitNode::unbiasedResidual() const {
  return std::visit( [&]( auto& f ) -> T { return f.unbiasedResidual( *this ); }, m_dim );
}

/// Interfaces for 'backward compatible' 1D float returns
template <>
inline double LHCb::FitNode::residual<double>() const {
  return visit( [&]( auto& f ) { return f.residual( *this )[0]; } );
}
template <>
inline double LHCb::FitNode::errResidual<double>() const {
  return visit( [&]( auto& f ) { return f.errResidual( *this )[0][0]; } );
}
template <>
inline double LHCb::FitNode::errResidual2<double>() const {
  return visit( [&]( auto& f ) { return f.errResidual2( *this )[0][0]; } );
}
template <>
inline double LHCb::FitNode::errMeasure2<double>() const {
  return visit( [&]( auto& f ) { return f.errMeasure2()[0]; } );
}
template <>
inline double LHCb::FitNode::unbiasedResidual<double>() const {
  return visit( [&]( auto& f ) { return f.unbiasedResidual( *this )[0]; } );
}

#endif // TRACKFITEVENT_FITNODE_H
