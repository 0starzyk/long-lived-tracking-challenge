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
#ifndef TRACKKERNEL_TRACKCOMPACTVERTEX_H
#define TRACKKERNEL_TRACKCOMPACTVERTEX_H

#include "Event/ZipUtils.h"
#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "GaudiKernel/Vector4DTypes.h"
#include "Kernel/Traits.h"
#include "LHCbMath/MatVec.h"
#include "LHCbMath/Vec3.h"
#include "Math/Point3D.h"
#include "Math/SMatrix.h"
#include "Math/Vector4D.h"
#include "TrackKernel/ChildRelation.h"
#include <array>
#include <vector>

namespace LHCb {
  namespace TrackKernel {
    // this can be replaced by std::dynamic_extent (defined from #include <span>) in c++20
    // OL 190704: used to be -1, but cling doesn't like converting that to std::size_t
    constexpr size_t Dynamic = std::numeric_limits<std::size_t>::max();

    namespace detail {
      // A helper class to create a container that can be either array or vector
      // FIXME: we can also extend this to have a dynamic size fixed
      // allocation. we then need another template field 'MaxSize',
      // just like Eigen uses.
      template <typename T, std::size_t Size>
      struct ArrayTrait {
        using Storage = std::array<T, Size>;
        static auto make( size_t /*size*/ ) { return Storage{}; }
      };

      template <typename T>
      struct ArrayTrait<T, Dynamic> {
        using Storage = std::vector<T>;
        static auto make( size_t size ) { return Storage( size ); }
      };

      // This container keeps the data for the daughters in the
      // vertex. The data for one daughter depends on the number of
      // daughters (since there is the NxN covariance matrix.) Since
      // we don't want more than one memory allocation, this makes the
      // navigation a bit more cumbersome, especially after I added
      // the track pointer to the data as well.
#define ITSEEMSFASTERWITHOUTTHECASTS
#ifdef ITSEEMSFASTERWITHOUTTHECASTS
      template <std::size_t N, typename FTYPE, typename TRACKPTR>
      struct DaughterDataStorage {
        DaughterDataStorage( size_t /*numdaughters*/ ) {}
        static constexpr size_t                   StorageSize  = N;
        static constexpr size_t                   NumDaughters = N;
        static constexpr size_t                   NumCols      = 1 + 3 + 3 + N;
        std::array<FTYPE, NumDaughters * NumCols> data;
      };
      template <typename FTYPE, typename TRACKPTR>
      struct DaughterDataStorage<LHCb::TrackKernel::Dynamic, FTYPE, TRACKPTR> {
        DaughterDataStorage( size_t numdaughters )
            : NumDaughters{numdaughters}, NumCols{1 + 3 + 3 + numdaughters}, data( NumDaughters * NumCols ) {}
        static constexpr size_t StorageSize = LHCb::TrackKernel::Dynamic;
        const size_t            NumDaughters;
        const size_t            NumCols;
        std::vector<FTYPE>      data;
      };

      // FIXME: It is a mystery to me why I need to specify the base class to access its members ...
      template <std::size_t N, typename FTYPE, typename TRACKPTR>
      struct DaughterData : public DaughterDataStorage<N, FTYPE, TRACKPTR> {
        using Base = DaughterDataStorage<N, FTYPE, TRACKPTR>;
        DaughterData( size_t numdaughters ) : Base{numdaughters} {}
        constexpr std::size_t size() const { return Base::NumDaughters; }
        /// scalar momentum squared (p2) of daughter
        FTYPE&       P( size_t dauindex ) { return Base::data[dauindex * Base::NumCols]; }
        const FTYPE& P( size_t dauindex ) const { return Base::data[dauindex * Base::NumCols]; }
        /// cov(x, daughter p2)
        FTYPE& Pposcov( size_t dauindex, size_t posindex ) {
          return Base::data[dauindex * Base::NumCols + 1 + posindex];
        }
        const FTYPE& Pposcov( size_t dauindex, size_t posindex ) const {
          return Base::data[dauindex * Base::NumCols + 1 + posindex];
        }
        /// cov(mother px,py,pz, daughter p2)
        FTYPE& Pmomcov( size_t dauindex, size_t momindex ) {
          return Base::data[dauindex * Base::NumCols + 4 + momindex];
        }
        const FTYPE& Pmomcov( size_t dauindex, size_t momindex ) const {
          return Base::data[dauindex * Base::NumCols + 4 + momindex];
        }
        /// cov( daughter p2 )
        FTYPE& Pcov( size_t dauindex, size_t daujndex ) { return Base::data[dauindex * Base::NumCols + 7 + daujndex]; }
        const FTYPE& Pcov( size_t dauindex, size_t daujndex ) const {
          return Base::data[dauindex * Base::NumCols + 7 + daujndex];
        }
      };

#else
      template <std::size_t N, typename FTYPE, typename TRACKPTR>
      struct DaughterDataStorage {
        DaughterDataStorage( size_t /*numdaughters*/ ) {}
        static constexpr size_t                    StorageSize  = N;
        static constexpr size_t                    NumDaughters = N;
        static constexpr size_t                    RowLength = ( 1 + 3 + 3 + N ) * sizeof( FTYPE ) + sizeof( TRACKPTR );
        std::array<char, NumDaughters * RowLength> data;
      };
      template <typename FTYPE, typename TRACKPTR>
      struct DaughterDataStorage<LHCb::TrackKernel::Dynamic, FTYPE, TRACKPTR> {
        DaughterDataStorage( size_t numdaughters )
            : NumDaughters{numdaughters}
            , RowLength{( 1 + 3 + 3 + numdaughters ) * sizeof( FTYPE ) + sizeof( TRACKPTR )}
            , data( NumDaughters * RowLength ) {}
        static constexpr size_t StorageSize = LHCb::TrackKernel::Dynamic;
        const size_t            NumDaughters;
        const size_t            RowLength;
        std::vector<char>       data;
      };

      // FIXME: It is a mystery to me why I need to specify the base class to access its members ...
      template <std::size_t N, typename FTYPE, typename TRACKPTR>
      struct DaughterData : public DaughterDataStorage<N, FTYPE, TRACKPTR> {
        using Base = DaughterDataStorage<N, FTYPE, TRACKPTR>;
        DaughterData( size_t numdaughters ) : Base{numdaughters} {}
        ///
        FTYPE* fdata( size_t dauindex ) {
          return reinterpret_cast<FTYPE*>( Base::data.data() + dauindex * Base::RowLength );
        }
        const FTYPE* fdata( size_t dauindex ) const {
          return reinterpret_cast<const FTYPE*>( Base::data.data() + dauindex * Base::RowLength );
        }
        /// scalar momentum squared (p2) of daughter
        FTYPE&       P( size_t dauindex ) { return fdata( dauindex )[0]; }
        const FTYPE& P( size_t dauindex ) const { return fdata( dauindex )[0]; }
        /// cov(x, daughter p2)
        FTYPE&       Pposcov( size_t dauindex, size_t posindex ) { return fdata( dauindex )[1 + posindex]; }
        const FTYPE& Pposcov( size_t dauindex, size_t posindex ) const { return fdata( dauindex )[1 + posindex]; }
        /// cov(mother px,py,pz, daughter p2)
        FTYPE&       Pmomcov( size_t dauindex, size_t momindex ) { return fdata( dauindex )[4 + momindex]; }
        const FTYPE& Pmomcov( size_t dauindex, size_t momindex ) const { return fdata( dauindex )[4 + momindex]; }
        /// cov( daughter p2 )
        FTYPE&       Pcov( size_t dauindex, size_t daujndex ) { return fdata( dauindex )[7 + daujndex]; }
        const FTYPE& Pcov( size_t dauindex, size_t daujndex ) const { return fdata( dauindex )[7 + daujndex]; }
        TRACKPTR&    track( size_t dauindex ) {
          return *reinterpret_cast<TRACKPTR*>( Base::data.data() + dauindex * Base::RowLength + Base::RowLength -
                                               sizeof( TRACKPTR ) );
        }
        const TRACKPTR& track( size_t dauindex ) const {
          return *reinterpret_cast<const TRACKPTR*>( Base::data.data() + dauindex * Base::RowLength + Base::RowLength -
                                                     sizeof( TRACKPTR ) );
        }
      };
#endif
    } // namespace detail

    // This would be the event model class for an N-prong vertex
    // FIXME: should we make this a class with private data?
    // FIXME: make this look a bit more like the interface of VertexBase (without actually using that as base class)
    using namespace ROOT::Math;
    class Track;
    template <std::size_t NDAUGHTERS, typename FTYPE = double, typename TRACKPTR = size_t>
    struct TrackCompactVertex {
      using XYZVector                          = DisplacementVector3D<Cartesian3D<FTYPE>, DefaultCoordinateSystemTag>;
      using XYZPoint                           = PositionVector3D<Cartesian3D<FTYPE>, DefaultCoordinateSystemTag>;
      using SymMatrix3x3                       = ROOT::Math::SMatrix<FTYPE, 3, 3, ROOT::Math::MatRepSym<FTYPE, 3>>;
      using Matrix3x3                          = ROOT::Math::SMatrix<FTYPE, 3, 3>;
      static constexpr std::size_t NumChildren = NDAUGHTERS;
      XYZVector                    mom; // momentum (px,py,pz)
      XYZPoint                     pos; // vertex   (x,y,z)
      SymMatrix3x3                 momcov;
      SymMatrix3x3                 poscov;
      Matrix3x3                    momposcov;
      FTYPE                        m_chi2;
      // so how can we do this such that we only have a single dynamic
      // allocation? just invent a confusing way of packing?
      detail::DaughterData<NDAUGHTERS, FTYPE, TRACKPTR> daughters;
      // TODO this needs reworking to support dynamic sizes
      std::array<ChildRelation, NumChildren> m_child_relations;

      TrackCompactVertex( size_t ndaughters ) : daughters{ndaughters} {}
      int                 nDoF() const { return daughters.size() * 2 - 3; }
      FTYPE               chi2() const { return m_chi2; }
      void                setChi2( FTYPE chi2 ) { m_chi2 = chi2; }
      FTYPE               chi2PerDoF() const { return m_chi2 / nDoF(); }
      XYZVector const&    threeMomentum() const { return mom; }
      SymMatrix3x3 const& threeMomCovMatrix() const { return momcov; }
      Matrix3x3 const&    threeMomPosCovMatrix() const { return momposcov; }
      // Interface compatibility with VertexBase
      XYZPoint const& position() const { return pos; }
      // Interface compatibility with Particle
      SymMatrix3x3 const& posCovMatrix() const { return poscov; }
      XYZPoint const&     referencePoint() const { return pos; }
      XYZPoint const&     endVertex() const { return pos; }

      // Interface compatibility with ThOr
      //
      friend LinAlg::Vec3<SIMDWrapper::scalar::float_v> referencePoint( const TrackCompactVertex& vtx ) {
        return {vtx.x(), vtx.y(), vtx.z()};
      }

      friend LinAlg::Vec3<SIMDWrapper::scalar::float_v> slopes( const TrackCompactVertex& vtx ) {
        if ( vtx.mom.Z() == 0 ) { return {0, 0, 1}; }
        return {vtx.mom.X() / vtx.mom.Z(), vtx.mom.Y() / vtx.mom.Z(), 1};
      }

      friend LinAlg::Vec3<SIMDWrapper::scalar::float_v> endVertexPos( const TrackCompactVertex& vtx ) {
        return {vtx.x(), vtx.y(), vtx.z()};
      }
      friend LinAlg::Vec3<SIMDWrapper::scalar::float_v> threeMomentum( const TrackCompactVertex& vtx ) {
        return {vtx.px(), vtx.py(), vtx.pz()};
      }
      friend auto posCovMatrix( const TrackCompactVertex& vtx ) {
        return LHCb::LinAlg::convert<SIMDWrapper::scalar::float_v>( vtx.poscov );
      }
      friend auto threeMomPosCovMatrix( const TrackCompactVertex& vtx ) {
        return LHCb::LinAlg::convert<SIMDWrapper::scalar::float_v>( vtx.momposcov );
      }
      friend auto threeMomCovMatrix( const TrackCompactVertex& vtx ) {
        return LHCb::LinAlg::convert<SIMDWrapper::scalar::float_v>( vtx.momcov );
      }

      friend auto decayProducts( const TrackCompactVertex& vtx ) {
        struct children_t {
          TrackCompactVertex const* vtx;
          struct Sentinel {};
          struct Proxy {
            bool                      operator!=( Sentinel ) const { return i != vtx->numChildren(); }
            TrackCompactVertex const* vtx;
            unsigned                  i;
            Proxy&                    operator++() {
              ++i;
              return *this;
            }
            const auto& operator*() const { return *this; }
            // hack - hack -- client code sofar only needs mag2 of a vector....
            // so instead we need another customization point which default to returning threeMomentum(x).mag2()...
            auto threeMomentum() const {
              return LHCb::LinAlg::Vec3<SIMDWrapper::scalar::float_v>{0, 0, vtx->daughters.P( i )};
            }
          };
          Proxy    begin() const { return {vtx, 0}; }
          Sentinel end() const { return {}; }
          Proxy    operator[]( unsigned i ) { return {vtx, i}; }
          auto     size() const { return vtx->numChildren(); }
        };
        return children_t{&vtx};
      }

      static constexpr auto canBeExtrapolatedDownstream = Event::CanBeExtrapolatedDownstream::no;

      XYZVector slopes() const {
        if ( mom.Z() ) {
          return XYZVector( mom.X() / mom.Z(), mom.Y() / mom.Z(), 1 );
        } else {
          return XYZVector( 0, 0, 1 );
        }
      }
      FTYPE                     x() const { return pos.X(); }
      FTYPE                     y() const { return pos.Y(); }
      FTYPE                     z() const { return pos.Z(); }
      FTYPE                     p() const { return mom.R(); }
      FTYPE                     pt() const { return mom.Rho(); }
      FTYPE                     px() const { return mom.X(); }
      FTYPE                     py() const { return mom.Y(); }
      FTYPE                     pz() const { return mom.Z(); }
      FTYPE                     phi() const { return mom.Phi(); }
      FTYPE                     pseudoRapidity() const { return mom.Eta(); }
      auto const&               child_relations() const { return m_child_relations; }
      [[nodiscard]] static auto numChildren() { return NumChildren; }
      [[nodiscard]] auto childRelationIndex( std::size_t n_child ) const { return m_child_relations[n_child].index(); }
      [[nodiscard]] auto childRelationFamily( std::size_t n_child ) const {
        return m_child_relations[n_child].zipIdentifier();
      }

      void setChildIndex( std::size_t child, std::size_t index ) { m_child_relations[child].setIndex( index ); }
      void setChildZipIdentifier( std::size_t child, Zipping::ZipFamilyNumber family ) {
        m_child_relations[child].setZipIdentifier( family );
      }
      constexpr std::size_t size() const { return daughters.size(); }
    };

    // global function to get everything we would like to fill in a particle
    template <typename TRACKCOMPACTVERTEX, typename MASSHYPOS>
    void computeParticleParams( const TRACKCOMPACTVERTEX& vertex, const MASSHYPOS& daughtermasses, Gaudi::XYZPoint& pos,
                                Gaudi::LorentzVector& p4, Gaudi::SymMatrix3x3& poscov, Gaudi::SymMatrix4x4& p4cov,
                                Gaudi::Matrix4x3& momposcov ) {
      // first fill in everything that can just be copied
      pos    = vertex.pos;
      poscov = vertex.poscov;
      for ( size_t ipos = 0; ipos < 3; ++ipos ) {
        for ( size_t imom = 0; imom < 3; ++imom ) momposcov( imom, ipos ) = vertex.momposcov( imom, ipos );
        momposcov( 3, ipos ) = 0;
      }
      for ( size_t imom = 0; imom < 3; ++imom ) {
        for ( size_t jmom = 0; jmom <= imom; ++jmom ) p4cov( imom, jmom ) = vertex.momcov( imom, jmom );
        p4cov( 3, imom ) = 0;
      }
      p4cov( 3, 3 ) = 0;
      p4.SetPx( vertex.mom.X() );
      p4.SetPy( vertex.mom.Y() );
      p4.SetPz( vertex.mom.Z() );
      // now add the mass dependent part. compute E and the missing pieces of the covariance matrix
      auto const N = vertex.size();
      double     E{0};
      // we need to cache daughter derivatives to P
      auto dEdPi = detail::ArrayTrait<double, TRACKCOMPACTVERTEX::NumChildren>::make( N );
      for ( size_t idau = 0; idau < N; ++idau ) {
        const double Pi = vertex.daughters.P( idau );
        const double Ei = std::sqrt( Pi * Pi + daughtermasses[idau] * daughtermasses[idau] );
        E += Ei;
        dEdPi[idau] = Pi / Ei;
      }
      p4.SetE( E );
      for ( size_t idau = 0; idau < N; ++idau ) {
        // position
        for ( size_t ipos = 0; ipos < 3; ++ipos )
          momposcov( 3, ipos ) += dEdPi[idau] * vertex.daughters.Pposcov( idau, ipos );
        // momentum
        for ( size_t imom = 0; imom < 3; ++imom )
          p4cov( 3, imom ) += dEdPi[idau] * vertex.daughters.Pmomcov( idau, imom );
        // energy (make sure not to double count)
        for ( size_t jdau = 0; jdau <= idau; ++jdau )
          p4cov( 3, 3 ) += dEdPi[idau] * vertex.daughters.Pcov( idau, jdau ) * dEdPi[jdau];
      }
    }
  } // namespace TrackKernel
} // namespace LHCb

#endif
