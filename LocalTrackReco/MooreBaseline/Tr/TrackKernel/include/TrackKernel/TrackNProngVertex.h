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
#ifndef TRACKKERNEL_TRACKNPRONGVERTEX_H
#define TRACKKERNEL_TRACKNPRONGVERTEX_H

// Include files
#include "Event/State.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "GaudiKernel/Vector4DTypes.h"
#include "LHCbMath/MatrixTransforms.h"
#include "TrackKernel/TrackCompactVertex.h"
#include <limits>

namespace LHCb {
  namespace TrackKernel {
    /** @class TrackNProngVertex TrackNProngVertex.h
     *
     * Vertex(ing) class for multitrack vertices
     *
     * @author Wouter Hulsbergen
     * created Wed Dec 12 14:58:51 2017
     *
     */

    /// class TrackNProngVertexDaughter. Helper class that represents a track in a vertex
    template <class STATETYPE = LHCb::State, typename FTYPE = double>
    class TrackNProngVertexDaughter final {
    public:
      using State              = STATETYPE;
      using Vector2            = ROOT::Math::SVector<FTYPE, 2>;
      using Vector3            = ROOT::Math::SVector<FTYPE, 3>;
      using SymMatrix2x2       = ROOT::Math::SMatrix<FTYPE, 2, 2, ROOT::Math::MatRepSym<FTYPE, 2>>;
      using SymMatrix3x3       = ROOT::Math::SMatrix<FTYPE, 3, 3, ROOT::Math::MatRepSym<FTYPE, 3>>;
      using PositionCovariance = SymMatrix3x3;
      using PositionParameters = ROOT::Math::SVector<FTYPE, 3>;
      using MomentumParameters = ROOT::Math::SVector<FTYPE, 3>;
      using TrackVector        = ROOT::Math::SVector<FTYPE, 5>;
      using XYZVector =
          ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<FTYPE>, ROOT::Math::DefaultCoordinateSystemTag>;

    private:
      const State m_state;
      FTYPE       m_weight{1.0};

      // In "fruhwirth's" notation
      SymMatrix2x2                     m_V; // cov matrix of (x,y) of state at vertex pos
      SymMatrix2x2                     m_G; // weight matrix of (x,y) of state
      ROOT::Math::SMatrix<FTYPE, 2, 3> m_A; // projection matrix for vertex position
      ROOT::Math::SMatrix<FTYPE, 3, 2> m_K; // momentum gain matrix (WBG in Fruhwirth)
      ROOT::Math::SVector<FTYPE, 3>    m_q; // predicted/fitted slope (tx,ty,q/p)

    public:
      // Constructor
      TrackNProngVertexDaughter( const State& state, FTYPE weight = 1.0F )
          : m_state{state}, m_weight{weight}, m_q{state.tx(), state.ty(), state.qOverP()} {}

      /// set the weight
      void setweight( FTYPE w ) { m_weight = w; }

      /// retrieve the weight
      FTYPE weight() const { return m_weight; }

      /// retrieve the momentum parameters (tx,ty,q/p)
      const auto& mom() const { return m_q; }

      /// input state
      State const& state() const { return m_state; }

      /// Position projection matrix (needed by alignment)
      const ROOT::Math::SMatrix<FTYPE, 5, 3>& A() const { return m_A; }

      /// compute chisq derivatives. returns 1 if fine.
      int project( const ROOT::Math::SVector<FTYPE, 3>& vertexpos, ROOT::Math::SVector<FTYPE, 3>& halfDChisqDX,
                   SymMatrix3x3& halfD2ChisqDX2, FTYPE& chi2 );

      /// update momentum for given change in position
      void updateSlopes( const ROOT::Math::SVector<FTYPE, 3>& vertexpos );

      /// add contribution of this track to vertexed four-vector
      // Important note: the part added to momcov is incomplete! This
      // is on purpose: the caller must take care of the
      // correlations using the gainmatrix. See my notes.
      void addToFourVector( const FTYPE mass, Gaudi::LorentzVector& sump4, Gaudi::SymMatrix4x4& sump4cov,
                            ROOT::Math::SMatrix<FTYPE, 4, 3>& sumgainmatrix ) const;

      /// add contribution of this track to vertexed
      /// three-vector. includes the 'momentum' part of the daughter
      /// that can later be used to compute the invariant mass once
      /// daughter masses are assigned.
      void addToThreeVector( XYZVector& sump3, SymMatrix3x3& sump3cov,
                             ROOT::Math::SMatrix<FTYPE, 3, 3>& summomgainmatrix, FTYPE& P, FTYPE& Pcov,
                             ROOT::Math::SMatrix<FTYPE, 1, 3>& Pmomcov,
                             ROOT::Math::SMatrix<FTYPE, 1, 3>& Pgainmatrix ) const;

      /// chisquare (ignores the weight)
      FTYPE chisq( const ROOT::Math::SVector<FTYPE, 3>& vertexpos ) const;
    };

    /// class TrackNProngVertex.
    template <std::size_t SIZE = Dynamic, typename FTYPE = double, typename STATETYPE = LHCb::State>
    class TrackNProngVertex final {
    public:
      inline constexpr static size_t StorageSize = SIZE;
      // size_type Size() const noexcept { return SIZE ; }
      using Ftype                = FTYPE;
      using State                = STATETYPE;
      using VertexTrack          = TrackNProngVertexDaughter<State, FTYPE>;
      using VertexTrackContainer = typename detail::ArrayTrait<VertexTrack, SIZE>::Storage;
      using MomentumParameters   = typename VertexTrack::MomentumParameters;
      using PositionParameters   = typename VertexTrack::PositionParameters;
      using PositionCovariance   = typename VertexTrack::PositionCovariance;
      using MomentumCovariance   = typename VertexTrack::PositionCovariance;
      enum FitStatus { FitSuccess, FitFailure, UnFitted };
      enum ErrCode {
        NoError              = 0,
        NotConverged         = 1,
        BadInputTrack        = 2,
        InsufficientTracks   = 4,
        ProjectionFailure    = 8,
        InversionFailure     = 16,
        NotConvergedAdaptive = 32
      };

    private:
      VertexTrackContainer m_daughters;
      PositionParameters   m_pos;
      PositionCovariance   m_poscov;
      FitStatus            m_fitStatus{UnFitted};
      ErrCode              m_error{NoError};
      // bool m_posInitialized = false ;
      FTYPE              m_chi2{-1};
      unsigned short     m_niter{0};
      PositionParameters m_refpos;    // position of reference position
      PositionCovariance m_refweight; // weight (inverse cov) of reference position
      PositionCovariance m_posweight; // inverse of the position covariance matrix

    public:
      /// Construct vertex from set of states.
      template <typename StateContainer>
      TrackNProngVertex( const StateContainer& states ) : m_daughters{states} {
        initPos();
      }

      /// Construct vertex from two states.
      TrackNProngVertex( const State& state1, const State& state2 ) : m_daughters{state1, state2} { initPos(); }

      /// Construct from a reference vertex. Then add track states
      /// later. If you use inverse of cov matrix, specify 'isweight=true'.
      // TrackNProngVertex(const Gaudi::XYZPoint& refposition, const Gaudi::SymMatrix3x3& refcovariance, bool
      // isweight=false) ;

      /// fit a single iteration. returns the delta-chisquare.
      FTYPE fitOneStep();

      /// fit until converged
      FitStatus fit( FTYPE maxdchisq = 0.01, size_t maxiterations = 10 );

      /// return the fit status
      FitStatus fitStatus() const { return m_fitStatus; }

      /// errcode gives some information on the type of error for failed fits
      ErrCode error() const { return m_error; }

      /// returns the chisquare of the vertex fit
      FTYPE chi2() const { return m_chi2; }

      /// number of dofs in vertex fit
      int nDoF() const { return nDaughters() * 2 - 3; }

      /// chi2/dof, accessor named to match Track and the new functors
      FTYPE chi2PerDoF() const { return m_chi2 / nDoF(); }

      /// Number of iterations that the vertex fit used
      auto nIter() const { return m_niter; }
      /// number of tracks in vertex
      auto nDaughters() const { return m_daughters.size(); }
      /// access the tracks in the vertex
      const auto& tracks() const { return m_daughters; }
      /// Position of the vertex
      auto position() const { return Gaudi::XYZPoint( m_pos( 0 ), m_pos( 1 ), m_pos( 2 ) ); }
      /// Position of the vertex as SVector
      const auto& pos() const { return m_pos; }
      /// Position covariance of the vertex
      const auto& covMatrix() const { return m_poscov; }
      /// Inverse of position covariance of the vertex
      const auto& weightMatrix() const { return m_posweight; }
      /// For   return j<=i ? momMomCovMatrixFast(i,j) : ROOT::Math::Transpose( momMomCovMatrixFast(j,i) ) ; }
      template <typename MassHypos>
      void computeParticleParams( const MassHypos& daughtermasses, Gaudi::XYZPoint& pos, Gaudi::LorentzVector& p4,
                                  Gaudi::SymMatrix3x3& poscov, Gaudi::SymMatrix4x4& p4cov,
                                  Gaudi::Matrix4x3& momposcov ) const;
      /// compute just the mother p4
      template <typename MassHypos>
      Gaudi::LorentzVector p4( const MassHypos& daughtermasses ) const;

      /// get the weight of a track after the adaptive fit
      auto weight( size_t i ) const { return m_daughters[i].weight(); }

      /// set the weight of a track
      void setWeight( size_t i, FTYPE w ) { m_daughters[i].setWeight( w ); }

    private:
      void initPos();
    };

    // FIXME: Here I didn't manage yet. I would like to write a (set of)
    // functions that returns the proper object depending on the
    // container type:
    // * for a vector, create one with dynamic storage
    // * for an array, an initializer list { } , or explicit arguments
    //   create one with fixed size storage

    /*
      template<typename FTYPE, typename STATETYPE, STATETYPE ... Arguments>
      auto makeVertex( Arguments... ) {
      return TrackNProngVertex<sizeof...(Arguments),STATETYPE,FTYPE>(STATETYPE(Argyments){...}) ;
      }
    */
    template <typename STATETYPE, typename FTYPE = double>
    auto makeVertex( const std::vector<STATETYPE>& states ) {
      return TrackNProngVertex<Dynamic, STATETYPE, FTYPE>( states );
    }
    template <int NDAUGHTERS, typename STATETYPE, typename FTYPE = double>
    auto makeVertex( const std::array<STATETYPE, NDAUGHTERS>& states ) {
      return TrackNProngVertex<NDAUGHTERS, STATETYPE, FTYPE>( states );
    }
    // template<typename STATETYPE, typename FTYPE=double>
    // auto makeVertex(std::initializer_list<STATETYPE> states ) {
    // return TrackNProngVertex<states.size(),STATETYPE,FTYPE>(states) ;
    //}
    template <int NDAUGHTERS, typename STATETYPE, typename FTYPE = double>
    auto makeVertex( const STATETYPE ( &states )[NDAUGHTERS] ) {
      return TrackNProngVertex<NDAUGHTERS, STATETYPE, FTYPE>( states );
    }
  } // namespace TrackKernel
} // namespace LHCb

// The rest could go into a icpp file
#include "Event/TrackVertexUtils.h"

namespace LHCb {
  namespace TrackKernel {
    /* Implementation for TrackNProngVertexDaughter */
    template <class STATETYPE, typename FTYPE>
    int TrackNProngVertexDaughter<STATETYPE, FTYPE>::project( const ROOT::Math::SVector<FTYPE, 3>& vertexpos,
                                                              ROOT::Math::SVector<FTYPE, 3>&       halfDChisqDX,
                                                              SymMatrix3x3& halfD2ChisqDX2, FTYPE& chi2 ) {
      // compute the weight matrix
      const auto  dz  = vertexpos( 2 ) - m_state.z();
      const auto  dz2 = dz * dz;
      const auto& V   = m_state.covariance();
      m_V( 0, 0 )     = V( 0, 0 ) + 2 * dz * V( 2, 0 ) + dz2 * V( 2, 2 );
      m_V( 1, 0 )     = V( 1, 0 ) + dz * ( V( 3, 0 ) + V( 2, 1 ) ) + dz2 * V( 3, 2 );
      m_V( 1, 1 )     = V( 1, 1 ) + 2 * dz * V( 3, 1 ) + dz2 * V( 3, 3 );
      m_G             = m_V;
      const auto rc   = m_G.Invert();
      // compute residual
      Vector2 res{vertexpos( 0 ) - ( m_state.x() + m_state.tx() * dz ),
                  vertexpos( 1 ) - ( m_state.y() + m_state.ty() * dz )};

      // fill the projection matrix: use the fitted momentum!
      m_A( 0, 0 ) = m_A( 1, 1 ) = 1;
      m_A( 0, 2 )               = -m_q( 0 );
      m_A( 1, 2 )               = -m_q( 1 );

      // I tried to make this faster by writing it out, but H does
      // not contain sufficiently manby zeroes. Better to
      // parallelize.
      halfD2ChisqDX2 += m_weight * ROOT::Math::Similarity( ROOT::Math::Transpose( m_A ), m_G );
      halfDChisqDX += m_weight * ( ROOT::Math::Transpose( m_A ) * ( m_G * res ) );
      chi2 += m_weight * ROOT::Math::Similarity( res, m_G );
      return rc;
    }

    template <class STATETYPE, typename FTYPE>
    void TrackNProngVertexDaughter<STATETYPE, FTYPE>::updateSlopes( const ROOT::Math::SVector<FTYPE, 3>& vertexpos ) {
      // first update the residual. (note the subtle difference with that in project!)
      const FTYPE dz = vertexpos( 2 ) - m_state.z();
      // compute residual
      Vector2 res{vertexpos( 0 ) - ( m_state.x() + m_state.tx() * dz ),
                  vertexpos( 1 ) - ( m_state.y() + m_state.ty() * dz )};

      // get the matrix that is the correlation of (x,y) and (tx,ty,qop)
      // This is the only place we need the 'extrapolated Vtx,x. Make sure to do it right!
      const auto&                      cov = m_state.covariance();
      ROOT::Math::SMatrix<FTYPE, 3, 2> Vba;
      Vba( 0, 0 ) = cov( 2, 0 ) + dz * cov( 2, 2 );
      Vba( 1, 0 ) = cov( 3, 0 ) + dz * cov( 3, 2 );
      Vba( 2, 0 ) = cov( 4, 0 ) + dz * cov( 4, 2 );
      Vba( 0, 1 ) = cov( 2, 1 ) + dz * cov( 3, 2 );
      Vba( 1, 1 ) = cov( 3, 1 ) + dz * cov( 3, 3 );
      Vba( 2, 1 ) = cov( 4, 1 ) + dz * cov( 4, 3 );
      // compute the corresponding gain matrix (this is WBG in the BFR fit, but here it is Vba * Vaa^-1)
      m_K = Vba * m_G;
      // compute the momentum vector
      m_q( 0 ) = m_state.tx();
      m_q( 1 ) = m_state.ty();
      m_q( 2 ) = m_state.qOverP();
      m_q += m_K * res;
    }

    template <class STATETYPE, typename FTYPE>
    void TrackNProngVertexDaughter<STATETYPE, FTYPE>::addToFourVector(
        const FTYPE mass, Gaudi::LorentzVector& sump4, Gaudi::SymMatrix4x4& sump4cov,
        ROOT::Math::SMatrix<FTYPE, 4, 3>& sumgainmatrix ) const {
      Gaudi::LorentzVector p4tmp;
      Gaudi::Math::geo2LA( m_q, mass, p4tmp );
      sump4 += p4tmp;
      ROOT::Math::SMatrix<FTYPE, 4, 3> dP4dMom;
      Gaudi::Math::JacobdP4dMom( m_q, mass, dP4dMom );
      ROOT::Math::SMatrix<FTYPE, 4, 2> FK = dP4dMom * m_K;
      sump4cov += ROOT::Math::Similarity( dP4dMom, m_state.covariance().template Sub<SymMatrix3x3>( 2, 2 ) );
      sump4cov -= ROOT::Math::Similarity( FK, m_V );
      sumgainmatrix += FK * m_A;
    }

    template <class STATETYPE, typename FTYPE>
    void TrackNProngVertexDaughter<STATETYPE, FTYPE>::addToThreeVector(
        XYZVector& sump3, SymMatrix3x3& sump3cov, ROOT::Math::SMatrix<FTYPE, 3, 3>& summomgainmatrix, FTYPE& P,
        FTYPE& Pcov, ROOT::Math::SMatrix<FTYPE, 1, 3>& Pmomcov, ROOT::Math::SMatrix<FTYPE, 1, 3>& Pgainmatrix ) const {
      // This could be moved to LHCb math (see also the case for p4 above)
      const auto& tx  = m_q( 0 );
      const auto& ty  = m_q( 1 );
      const auto& qop = m_q( 2 );
      const auto  p   = 1 / std::abs( qop );
      const auto  n2  = 1 + tx * tx + ty * ty;
      const auto  n   = std::sqrt( n2 );
      const auto  pz  = p / n;
      const auto  px  = pz * tx;
      const auto  py  = pz * ty;
      sump3 += XYZVector{px, py, pz};
      ROOT::Math::SMatrix<FTYPE, 3, 3> dP3dMom;
      const auto                       n3 = n2 * n;
      dP3dMom( 0, 0 )                     = p * ( 1 + ty * ty ) / n3; // dpx/dtx
      dP3dMom( 0, 1 )                     = p * tx * -ty / n3;        // dpx/dty
      dP3dMom( 0, 2 )                     = -px / qop;                // dpx/dqop
      dP3dMom( 1, 0 )                     = p * ty * -tx / n3;        // dpy/dtx
      dP3dMom( 1, 1 )                     = p * ( 1 + tx * tx ) / n3; // dpy/dty
      dP3dMom( 1, 2 )                     = -py / qop;                // dpy/dqop
      dP3dMom( 2, 0 )                     = pz * -tx / n2;            // dpz/dtx
      dP3dMom( 2, 1 )                     = pz * -ty / n2;            // dpz/dtx
      dP3dMom( 2, 2 )                     = -pz / qop;                // dpz/dqop
      SymMatrix3x3 covT                   = m_state.covariance().template Sub<SymMatrix3x3>( 2, 2 );
      sump3cov += ROOT::Math::Similarity( dP3dMom, covT );
      ROOT::Math::SMatrix<FTYPE, 3, 2> FK = dP3dMom * m_K;
      sump3cov -= ROOT::Math::Similarity( FK, m_V );
      summomgainmatrix += FK * m_A;
      // compute P and the gain matrix for p2
      P = p;
      // This can be optimized quite a bit: We can compute m_K V m_K^T upfront and then reuse it.
      auto        KA      = m_K * m_A;
      const FTYPE dPdqop  = -P / m_q( 2 );
      Pgainmatrix( 0, 0 ) = dPdqop * KA( 2, 0 );
      Pgainmatrix( 0, 1 ) = dPdqop * KA( 2, 1 );
      Pgainmatrix( 0, 2 ) = dPdqop * KA( 2, 2 );
      // as for the other one, this still misses the contribution from cov(x): the caller should add that.
      SymMatrix3x3 KVK = ROOT::Math::Similarity( m_K, m_V );
      Pcov             = dPdqop * ( covT( 2, 2 ) - KVK( 2, 2 ) ) * dPdqop;
      // auto covTdP3dMomT = covT * ROOT::Math::Transpose( dP3dMom ) ;
      ROOT::Math::SMatrix<FTYPE, 3, 3> covTMinusKVKdP3dMomT = ( covT - KVK ) * ROOT::Math::Transpose( dP3dMom );
      Pmomcov( 0, 0 )                                       = dPdqop * covTMinusKVKdP3dMomT( 2, 0 );
      Pmomcov( 0, 1 )                                       = dPdqop * covTMinusKVKdP3dMomT( 2, 1 );
      Pmomcov( 0, 2 )                                       = dPdqop * covTMinusKVKdP3dMomT( 2, 2 );
    }

    template <class STATETYPE, typename FTYPE>
    FTYPE TrackNProngVertexDaughter<STATETYPE, FTYPE>::chisq( const ROOT::Math::SVector<FTYPE, 3>& vertexpos ) const {
      const FTYPE dz = vertexpos( 2 ) - m_state.z();
      Vector2     res{vertexpos( 0 ) - ( m_state.x() + m_state.tx() * dz ),
                  vertexpos( 1 ) - ( m_state.y() + m_state.ty() * dz )};
      return ROOT::Math::Similarity( res, m_G );
    }

    /* Implementation for TrackNProngVertex */

    template <std::size_t NDAUGHTERS, typename FTYPE, typename STATETYPE>
    void TrackNProngVertex<NDAUGHTERS, FTYPE, STATETYPE>::initPos() {
      // we gain about one iteration by initializing the position with the doca
      Gaudi::XYZPoint docapoint;
      int success = LHCb::TrackVertexUtils::poca( m_daughters[0].state(), m_daughters[1].state(), docapoint );
      if ( !success ) m_fitStatus = FitFailure;
      m_pos( 0 ) = docapoint.x();
      m_pos( 1 ) = docapoint.y();
      m_pos( 2 ) = docapoint.z();
    }

    template <std::size_t NDAUGHTERS, typename FTYPE, typename STATETYPE>
    FTYPE TrackNProngVertex<NDAUGHTERS, FTYPE, STATETYPE>::fitOneStep() {
      // add one to the counter
      m_niter += 1;
      // invalidate some caches
      // m_validMomCovCache = false ;

      // This implements my optimized algorithm
      // adds the reference position
      m_posweight                                = m_refweight;
      ROOT::Math::SVector<FTYPE, 3> halfDChisqDX = m_refweight * ( m_pos - m_refpos );
      m_chi2                                     = 0;
      // add all the tracks
      for ( auto& itrack : m_daughters ) itrack.project( m_pos, halfDChisqDX, m_posweight, m_chi2 );

      // calculate the covariance and the change in the position
      FTYPE dchisq  = -1;
      m_poscov      = m_posweight;
      const bool ok = m_poscov.InvertChol();
      if ( !ok ) {
        m_fitStatus = FitFailure;
        m_error     = InversionFailure;
      } else {
        // compute the delta
        PositionParameters dpos = -m_poscov * halfDChisqDX;
        // update the position
        m_pos += dpos;
        // update the momenta
        for ( auto& itrack : m_daughters ) itrack.updateSlopes( m_pos );
        // return the delta-chisquare
        dchisq = ROOT::Math::Dot( dpos, halfDChisqDX );
        m_chi2 += dchisq;
        m_fitStatus = FitSuccess;
      }
      return dchisq;
    }

    template <std::size_t NDAUGHTERS, typename FTYPE, typename STATETYPE>
    typename TrackNProngVertex<NDAUGHTERS, FTYPE, STATETYPE>::FitStatus
    TrackNProngVertex<NDAUGHTERS, FTYPE, STATETYPE>::fit( FTYPE maxdchisq, size_t maxnumiter ) {
      bool   converged( false );
      size_t iter( 0 );
      m_fitStatus = FitSuccess;
      for ( ; iter < maxnumiter && !converged && m_fitStatus == FitSuccess; ++iter ) {
        auto dchisq = fitOneStep();
        converged   = -dchisq < maxdchisq;
      }
      if ( !converged ) {
        m_fitStatus = FitFailure;
        m_error     = NotConverged;
        // std::cout << "TrackNProngVertex fit failure: not converging" << std::endl ;
      }
      return m_fitStatus;
    }

    template <std::size_t NDAUGHTERS, typename FTYPE, typename STATETYPE>
    template <typename MASSHYPOS>
    void TrackNProngVertex<NDAUGHTERS, FTYPE, STATETYPE>::computeParticleParams(
        const MASSHYPOS& masshypos, Gaudi::XYZPoint& pos, Gaudi::LorentzVector& p4, Gaudi::SymMatrix3x3& poscov,
        Gaudi::SymMatrix4x4& p4cov, Gaudi::Matrix4x3& momposcov ) const {
      pos = position();
      // make sure they are empty?
      p4              = Gaudi::LorentzVector{};
      p4cov           = Gaudi::SymMatrix4x4{};
      using Matrix4x3 = ROOT::Math::SMatrix<FTYPE, 4, 3>;
      Matrix4x3    gainmatrix;
      const size_t dim = m_daughters.size();
      for ( size_t index = 0; index < dim; ++index )
        m_daughters[index].addToFourVector( masshypos[index], p4, p4cov, gainmatrix );
      p4cov += ROOT::Math::Similarity( gainmatrix, m_poscov );
      momposcov = gainmatrix * m_poscov;
      poscov    = m_poscov;
    }

    template <std::size_t NDAUGHTERS, typename FTYPE, typename STATETYPE>
    template <typename MASSHYPOS>
    Gaudi::LorentzVector TrackNProngVertex<NDAUGHTERS, FTYPE, STATETYPE>::p4( const MASSHYPOS& daughtermasses ) const {
      Gaudi::LorentzVector p4tot;
      for ( size_t index = 0; index < m_daughters.size(); ++index ) {
        Gaudi::LorentzVector p4tmp;
        Gaudi::Math::geo2LA( m_daughters[index].mom(), daughtermasses[index], p4tmp );
        p4tot += p4tmp;
      }
      return p4tot;
    }

    template <typename TrackNProngVertex>
    auto convertToCompactVertex( const TrackNProngVertex& vertex ) {
      using Ftype                  = typename TrackNProngVertex::Ftype;
      using TrackCompactVertexType = TrackCompactVertex<TrackNProngVertex::StorageSize, Ftype>;
      TrackCompactVertexType cvertex{vertex.nDaughters()};
      cvertex.pos = vertex.position();
      cvertex.setChi2( vertex.chi2() );
      cvertex.poscov = vertex.covMatrix();
      // std::cout << "TCV: " << cvertex.daughters.NumDaughters << " "
      //<< cvertex.daughters.NumCols << " " << cvertex.daughters.data.size()
      //	  << " " << sizeof(cvertex.daughters) << std::endl ;
      // cvertex.momcov = ThisTrackCompactVertex::SymMatrix3x3{} ;
      // cvertex.mom    = ThisTrackCompactVertex::XYZVector{} ;
      typename TrackCompactVertexType::Matrix3x3 momposgainmatrix;
      // temporary structure to holds parts we need for every daughter
      struct DaughterPParts {
        Ftype                            cov;
        ROOT::Math::SMatrix<Ftype, 1, 3> momcov;
        ROOT::Math::SMatrix<Ftype, 1, 3> posgainmatrix;
      };
      auto dauPparts = detail::ArrayTrait<DaughterPParts, TrackNProngVertex::StorageSize>::make( vertex.nDaughters() );
      for ( size_t idau{0}; idau < vertex.nDaughters(); ++idau )
        vertex.tracks()[idau].addToThreeVector( cvertex.mom, cvertex.momcov, momposgainmatrix,
                                                cvertex.daughters.P( idau ), dauPparts[idau].cov,
                                                dauPparts[idau].momcov, dauPparts[idau].posgainmatrix );
      cvertex.momcov += ROOT::Math::Similarity( momposgainmatrix, cvertex.poscov );
      cvertex.momposcov = momposgainmatrix * cvertex.poscov;
      // fill the various parts of the daughter covariance matrices
      for ( size_t idau{0}; idau < vertex.nDaughters(); ++idau ) {
        // position
        auto Pposcov = dauPparts[idau].posgainmatrix * cvertex.poscov;
        for ( size_t ipos = 0; ipos < 3; ++ipos ) cvertex.daughters.Pposcov( idau, ipos ) = Pposcov( 0, ipos );
        // const auto Pmomcov = dauPparts[idau].momcov + dauPparts[idau].posgainmatrix * cvertex.poscov *
        // ROOT::Math::Transpose( momposgainmatrix ) ;
        ROOT::Math::SMatrix<Ftype, 1, 3> Pmomcov =
            dauPparts[idau].momcov +
            dauPparts[idau].posgainmatrix * cvertex.poscov * ROOT::Math::Transpose( momposgainmatrix );
        for ( size_t imom = 0; imom < 3; ++imom ) cvertex.daughters.Pmomcov( idau, imom ) = Pmomcov( 0, imom );
        // P,P cov
        cvertex.daughters.Pcov( idau, idau ) =
            dauPparts[idau].cov + ( ROOT::Math::Similarity( dauPparts[idau].posgainmatrix, cvertex.poscov ) )( 0, 0 );
        // note that we fill them both, even though it is symmetric. that's just to simplify things
        for ( size_t jdau{0}; jdau < idau; ++jdau ) {
          cvertex.daughters.Pcov( idau, jdau ) = cvertex.daughters.Pcov( jdau, idau ) =
              ( dauPparts[idau].posgainmatrix * cvertex.poscov *
                ROOT::Math::Transpose( dauPparts[jdau].posgainmatrix ) )( 0, 0 );
        }
      }
      return cvertex;
    }
  } // namespace TrackKernel
} // namespace LHCb

#endif
