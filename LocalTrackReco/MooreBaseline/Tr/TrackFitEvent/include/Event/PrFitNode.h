/*****************************************************************************\
* (c) Copyright 2020 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#pragma once

// LHCb
#include "Event/ChiSquare.h"
#include "Event/State.h"
#include "Event/StateParameters.h"
#include "Event/StateVector.h"
#include "Event/Track.h"
#include "Event/TrackTypes.h"
#include "Kernel/LHCbID.h"
#include "Kernel/STLExtensions.h"
#include "LHCbMath/Similarity.h"
// Gaudi
#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/Kernel.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "GaudiKernel/reverse.h"
// std
#include <array>
#include <limits>
#include <memory>

namespace LHCb::Pr::Tracks::Fit {
  struct Node final {
    enum class Type {
      Unknown = 0,
      VPHit,
      UTHit,
      FTHit,
      MuonHit,
      BegRich1,
      EndRich1
    }; // the real measurements are first in the enum

    int constexpr static forward  = 0;
    int constexpr static backward = 1;

    Node( double pos_x, double pos_y, double pos_z, double dir_x, double dir_y, double dir_z, double err, Type type,
          LHCb::LHCbID lhcbid )
        : measurement_error{err}
        , m_z{pos_z}
        , m_type{type}
        , lhcbID{lhcbid}
        , measurement_pos{pos_x, pos_y, pos_z}
        , measurement_dir{dir_x, dir_y, dir_z} {}

    double                               ref_residual{0.0};
    double                               measurement_error{std::numeric_limits<double>::signaling_NaN()};
    double                               delta_energy{0.0};
    double                               m_z{std::numeric_limits<double>::signaling_NaN()};
    Type                                 m_type{Type::Unknown};
    LHCb::LHCbID                         lhcbID{};
    bool                                 m_is_outlier{false};
    LHCb::StateVector                    lhcb_ref_vector{{}, m_z};
    std::array<LHCb::ChiSquare, 2>       delta_chi2{};
    Gaudi::TrackProjectionMatrix1D       projection{}; // H in EKF formalism
    std::array<Gaudi::TrackVector, 2>    predicted_state_vec{};
    std::array<Gaudi::TrackSymMatrix, 2> predicted_state_cov{};
    Gaudi::TrackSymMatrix                noise_matrix{};     // Q in EKF
    Gaudi::TrackMatrix                   transport_matrix{}; // F in EKF formalism
    Gaudi::TrackMatrix                   transport_matrix_inverse{};
    Gaudi::TrackVector                   transport_vector{};
    std::array<Gaudi::TrackVector, 2>    filtered_state_vec{};
    std::array<Gaudi::TrackSymMatrix, 2> filtered_state_cov{};
    Gaudi::Vector3                       measurement_pos{};
    Gaudi::Vector3                       measurement_dir{};
    Gaudi::TrackVector                   final_state_vec{};
    Gaudi::TrackSymMatrix                final_state_cov{};

    [[nodiscard]] Gaudi::TrackVector& ref_vector() { return lhcb_ref_vector.parameters(); }

    // needed for common interface with LHCb::FitNode in param scatter
    [[nodiscard]] double                   z() const { return m_z; }
    [[nodiscard]] Type                     type() const { return m_type; }
    [[nodiscard]] LHCb::StateVector const& refVector() const { return lhcb_ref_vector; }

    // needed for common interface with LHCb::FitNode in alignment
    [[nodiscard]] double                  errMeasure2() const { return measurement_error * measurement_error; }
    const Gaudi::TrackProjectionMatrix1D& projectionMatrix() const { return projection; }

    void set_delta_energy( double dE ) { delta_energy = dE; }
    void set_noise_matrix( LHCb::State const& s ) { noise_matrix = s.covariance(); }
    void set_transport( Gaudi::TrackMatrix const& F, Gaudi::TrackVector const& transportvec ) {
      transport_matrix         = F;
      transport_matrix_inverse = F;
      // no B-field means -> straight line -> easy to invert
      if ( F( 0, 4 ) == 0 ) {
        transport_matrix_inverse( 0, 2 ) = -F( 0, 2 );
        transport_matrix_inverse( 1, 3 ) = -F( 1, 3 );
      } else {
        // transport_matrix_inverse(0,0) = transport_matrix_inverse(1,1) = transport_matrix_inverse(4,4) = 1 ;
        // write
        //      ( 1  0 |  S00 S01 | U0 )
        //      ( 0  1 |  S10 S01 | U1 )
        // F =  ( 0  0 |  T00 T01 | V0 )
        //      ( 0  0 |  T10 T11 | V1 )
        //      ( 0  0 |   0   0  | 1  )
        // then we have
        // Tinv = T^{-1}
        double det                       = F( 2, 2 ) * F( 3, 3 ) - F( 2, 3 ) * F( 3, 2 );
        transport_matrix_inverse( 2, 2 ) = F( 3, 3 ) / det;
        transport_matrix_inverse( 3, 3 ) = F( 2, 2 ) / det;
        transport_matrix_inverse( 2, 3 ) = -F( 2, 3 ) / det;
        transport_matrix_inverse( 3, 2 ) = -F( 3, 2 ) / det;
        // Vinv = - T^-1 * V
        transport_matrix_inverse( 2, 4 ) =
            -transport_matrix_inverse( 2, 2 ) * F( 2, 4 ) - transport_matrix_inverse( 2, 3 ) * F( 3, 4 );
        transport_matrix_inverse( 3, 4 ) =
            -transport_matrix_inverse( 3, 2 ) * F( 2, 4 ) - transport_matrix_inverse( 3, 3 ) * F( 3, 4 );
        // Uinv = S * T^-1 * V - U = - S * Vinv - U
        transport_matrix_inverse( 0, 4 ) =
            -F( 0, 4 ) - F( 0, 2 ) * transport_matrix_inverse( 2, 4 ) - F( 0, 3 ) * transport_matrix_inverse( 3, 4 );
        transport_matrix_inverse( 1, 4 ) =
            -F( 1, 4 ) - F( 1, 2 ) * transport_matrix_inverse( 2, 4 ) - F( 1, 3 ) * transport_matrix_inverse( 3, 4 );
        // Sinv  = - S * T^{-1}
        transport_matrix_inverse( 0, 2 ) =
            -F( 0, 2 ) * transport_matrix_inverse( 2, 2 ) - F( 0, 3 ) * transport_matrix_inverse( 3, 2 );
        transport_matrix_inverse( 0, 3 ) =
            -F( 0, 2 ) * transport_matrix_inverse( 2, 3 ) - F( 0, 3 ) * transport_matrix_inverse( 3, 3 );
        transport_matrix_inverse( 1, 2 ) =
            -F( 1, 2 ) * transport_matrix_inverse( 2, 2 ) - F( 1, 3 ) * transport_matrix_inverse( 3, 2 );
        transport_matrix_inverse( 1, 3 ) =
            -F( 1, 2 ) * transport_matrix_inverse( 2, 3 ) - F( 1, 3 ) * transport_matrix_inverse( 3, 3 );
      }
      transport_vector = transportvec;
    }

    void project_reference() {
      Gaudi::Vector3 const r_pos{lhcb_ref_vector.x(), lhcb_ref_vector.y(), z()};
      Gaudi::Vector3 const r_dir{lhcb_ref_vector.tx(), lhcb_ref_vector.ty(), 1.0};

      Gaudi::Vector3 const delta_pos = r_pos - measurement_pos;
      Gaudi::Vector3 const v         = Cross( measurement_dir, r_dir );
      // I need the norms so calculate them here
      double const inv_norm2 = 1. / Dot( v, v );
      double const inv_norm  = std::sqrt( inv_norm2 );
      // and use them to normalize the vector instead of calling the member funtion
      // which would not reuse the inverted norms above
      Gaudi::Vector3 const n             = v * inv_norm;
      Gaudi::Vector3 const meas_dir_norm = measurement_dir * inv_norm;

      // calculate d(n_{x,y,z})/d_tx
      double const dv_dtx   = meas_dir_norm[2] * v[1] - meas_dir_norm[1] * v[2];
      double const dv_dtx_n = inv_norm2 * dv_dtx;
      double const dnx_dtx  = -v[0] * dv_dtx_n;
      double const dny_dtx  = -v[1] * dv_dtx_n + meas_dir_norm[2];
      double const dnz_dtx  = -v[2] * dv_dtx_n - meas_dir_norm[1];

      // calculate d(n_{x,y,z})/d_ty
      double const dv_dty   = meas_dir_norm[0] * v[2] - meas_dir_norm[2] * v[0];
      double const dv_dty_n = inv_norm2 * dv_dty;
      double const dnx_dty  = -v[0] * dv_dty_n - meas_dir_norm[2];
      double const dny_dty  = -v[1] * dv_dty_n;
      double const dnz_dty  = -v[2] * dv_dty_n + meas_dir_norm[0];

      // projection[0,1] = nx, ny
      projection( 0, 0 ) = n[0];
      projection( 0, 1 ) = n[1];
      // projection[2] = delta_x * dnx/dtx + delta_y * dny/dtx + delta_z * dnz/dtx
      projection( 0, 2 ) = delta_pos[0] * dnx_dtx + delta_pos[1] * dny_dtx + delta_pos[2] * dnz_dtx;
      // projection[3] = delta_x * dnx/dty + delta_y * dny/dty + delta_z * dnz/dty
      projection( 0, 3 ) = delta_pos[0] * dnx_dty + delta_pos[1] * dny_dty + delta_pos[2] * dnz_dty;
      // projection marix is independent of q/p
      projection( 0, 4 ) = 0;

      ref_residual = -Dot( delta_pos, n );
    }

    bool isHitOnTrack() const { return ( !isOutlier() && hasMeasurement() ); }

    bool isOutlier() const { return m_is_outlier; }

    bool isVP() const { return m_type == Type::VPHit; }

    bool isUT() const { return m_type == Type::UTHit; }

    bool isFT() const { return m_type == Type::FTHit; }

    bool isMuon() const { return m_type == Type::MuonHit; }

    bool hasMeasurement() const { return ( isVP() || isUT() || isFT() || isMuon() ); }

    [[nodiscard]] double residual() const {
      const Gaudi::TrackVector& refX = refVector().parameters();
      return ref_residual + ( projection * ( refX - final_state_vec ) )( 0 );
    }

    [[nodiscard]] double errResidual() const { return std::sqrt( errResidual2() ); }

    [[nodiscard]] double errResidual2() const {
      double V   = errMeasure2();
      double HCH = LHCb::Math::Similarity( projection, final_state_cov )[0][0];
      if ( isHitOnTrack() )
        return V - HCH;
      else
        return V + HCH;
    }

    [[nodiscard]] double unbiasedResidual() const {
      auto r = residual();
      if ( isHitOnTrack() )
        return r * errMeasure2() / errResidual2();
      else
        return r;
    }

    [[nodiscard]] double unbiasedErrResidual() const {
      auto R = errResidual();
      if ( isHitOnTrack() ) return errMeasure2() / R;
      return R;
    }

    inline LHCb::State unbiasedState() const {
      const auto& H       = projection;
      const auto& biasedC = final_state_cov;
      double      r       = residual();
      double      R       = errResidual2();

      ROOT::Math::SMatrix<double, 5, 1> K = ( biasedC * Transpose( H ) ) / R;

      Gaudi::TrackVector unbiasedX = final_state_vec - K.Col( 0 ) * r;

      const Gaudi::TrackSymMatrix unit = ROOT::Math::SMatrixIdentity();
      Gaudi::TrackSymMatrix       unbiasedC;
      ROOT::Math::AssignSym::Evaluate( unbiasedC, ( unit + K * H ) * biasedC );
      return LHCb::State{unbiasedX, unbiasedC, m_z, LHCb::State::Location::LocationUnknown};
    }
  };

  // Functions bellow are used to ease templating of various algorithms on fit node type.
  // Each function has a twin in namespace LHCb in FitNode.h
  // which returns the apropriate LHCb::FitNode alternative. This allows making loops over fit nodes
  // look the same and thus helps with templating, while keeping the PrFitNode object "lightweight".
  inline const LHCb::LHCbID     id( const Node& node ) { return node.lhcbID; }
  inline const Gaudi::XYZVector pocaVector( const Node& node ) {
    Gaudi::Vector3 d0 = node.measurement_dir;
    Gaudi::Vector3 d1{node.lhcb_ref_vector.tx(), node.lhcb_ref_vector.ty(), 1.0};
    auto           unit_poca = Cross( d0, d1 ).Unit();
    return Gaudi::XYZVector{unit_poca[0], unit_poca[1], unit_poca[2]};
  }
  inline const LHCb::State state( const Node& node ) {
    return LHCb::State( node.final_state_vec, node.final_state_cov, node.m_z, LHCb::State::Location::LocationUnknown );
  }
  inline const LHCb::State predictedStateForward( const Node& node ) {
    return LHCb::State( node.predicted_state_vec[0], node.predicted_state_cov[0], node.m_z,
                        LHCb::State::Location::LocationUnknown );
  }
  inline const LHCb::State predictedStateBackward( const Node& node ) {
    return LHCb::State( node.predicted_state_vec[1], node.predicted_state_cov[1], node.m_z,
                        LHCb::State::Location::LocationUnknown );
  }
  inline const LHCb::State filteredStateForward( const Node& node ) {
    return LHCb::State( node.filtered_state_vec[0], node.filtered_state_cov[0], node.m_z,
                        LHCb::State::Location::LocationUnknown );
  }
  inline const LHCb::State filteredStateBackward( const Node& node ) {
    return LHCb::State( node.filtered_state_vec[1], node.filtered_state_cov[1], node.m_z,
                        LHCb::State::Location::LocationUnknown );
  }

  inline const Gaudi::TrackSymMatrix& smoothedStateCovariance( const Node& node ) { return node.final_state_cov; };
  inline const LHCb::State            smoothedState( const Node& node ) {
    return LHCb::State{node.final_state_vec, node.final_state_cov, node.m_z, LHCb::State::Location::LocationUnknown};
  };

} // namespace LHCb::Pr::Tracks::Fit
