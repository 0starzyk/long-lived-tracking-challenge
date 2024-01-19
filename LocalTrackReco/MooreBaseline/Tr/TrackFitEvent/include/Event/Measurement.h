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
#pragma once
#include "Event/StateZTraj.h"
#include "Event/TrackTypes.h"
#include "GaudiKernel/compose.h"
#include "Kernel/LHCbID.h"
#include "Kernel/LineTraj.h"
#include "Kernel/STLExtensions.h"
#include "Kernel/Trajectory.h"
#include "Kernel/meta_enum.h"
#include <type_traits>
#include <variant>

// VP
#include "VPDet/DeVPSensor.h"

// UT
#include "UTDet/DeUTSector.h"

// FT
#include "FTDet/DeFTMat.h"

// Muon
#include "MuonDet/DeMuonChamber.h"

namespace LHCb {
  /** @class Measurement Measurement.h
   *
   * Measurement defines the functionality required for the
   * Kalman Filter, Alignment and Monitoring
   *
   * @author Jose Hernando, Eduardo Rodrigues
   * created Fri Jan  4 20:26:46 2019
   *
   */
  namespace Enum::Measurement {
    /// enumerator for the type of Measurement
    meta_enum_class( Type, unsigned char, Unknown, Unknown2, Unknown3, Unknown4, Unknown5, Unknown6, Unknown7, Unknown8,
                     Muon, Unknown9, Unknown10, VP, Calo, Origin, FT, UT = 19, VP2D )
  } // namespace Enum::Measurement

  namespace details::Measurement {

    template <typename Dir>
    bool isX( const Dir& dir ) {
      return dir.x() > dir.y();
    }

    template <typename T, typename Variant, std::size_t... I>
    constexpr int index_helper( std::index_sequence<I...> ) {
      auto b = std::array{std::is_same_v<T, std::variant_alternative_t<I, Variant>>...};
      for ( int i = 0; i != static_cast<int>( size( b ) ); ++i ) {
        if ( b[i] ) return i;
      }
      return -1;
    }

    template <typename T, typename Variant>
    constexpr int index() {
      constexpr auto idx = index_helper<T, Variant>( std::make_index_sequence<std::variant_size_v<Variant>>() );
      static_assert( idx != -1, "T does not appear in Variant" );
      return idx;
    }

  } // namespace details::Measurement

  using Point =
      typename ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>, ROOT::Math::DefaultCoordinateSystemTag>;

  // Scaffolding for visit overloading with another variant
  template <typename T>
  struct is_variant_ : std::false_type {};
  template <typename... Ts>
  struct is_variant_<std::variant<Ts...>> : std::true_type {};
  template <typename T>
  constexpr bool is_variant_v = is_variant_<std::decay_t<T>>::value;

  class Measurement final {
  public:
    struct FT {
      LHCb::LineTraj<double> trajectory;
#ifdef USE_DD4HEP
      DeFTMat detElem{};
      auto&   mat() const { return detElem; }
#else
      const DeFTMat*       detElem = nullptr;
      auto&                mat() const { return *detElem; }
#endif
      double errMeasure;
    };
    struct Muon {
      LHCb::LineTraj<double> trajectory;
#ifdef USE_DD4HEP
      DeMuonChamber detElem{};
      auto&         chamber() const { return detElem; }
#else
      const DeMuonChamber* detElem = nullptr;
      auto&                chamber() const { return *detElem; }
#endif
      double errMeasure;
    };
    struct UT {
      LHCb::LineTraj<double> trajectory;
#ifdef USE_DD4HEP
      DeUTSector detElem{};
      auto&      sec() const { return detElem; }
#else
      const DeUTSector*    detElem = nullptr;
      auto&                sec() const { return *detElem; }
#endif
      double errMeasure;
    };
    struct VP {
      /// x or y measurement
      enum class Projection { X = 1, Y = 2, XY = 3 };
      LHCb::LineTraj<double> trajectory;
#ifdef USE_DD4HEP
      DeVPSensor detElem{};
      auto&      sensor() const { return detElem; }
#else
      const DeVPSensor*    detElem = nullptr;
      auto&                sensor() const { return *detElem; }
#endif
      double     errMeasure;
      bool       isX;
      Projection projection() const {
        return isX ? VP::Projection::Y // measurement is perpendicular to direction
                   : VP::Projection::X;
      }
    };
    struct VP2D {
      Point trajectory;
#ifdef USE_DD4HEP
      DeVPSensor detElem{};
#else
      const DeVPSensor*    detElem = nullptr;
#endif
      auto& sensor() const {
#ifdef USE_DD4HEP
        return detElem;
#else
        return *detElem;
#endif
      }
      ROOT::Math::SVector<double, 2> errMeasure; ///< the measurement error
    };

  private:
    // Note: for now, use a variant. Alternatively, store this information 'elsewhere,
    //       and put a "pointer" (which could be "container + index" on top of SOA storage)
    //       to 'elsewhere' into the variant instead.
    using SubInfo = std::variant<VP, UT, FT, Muon, VP2D>;
    double  m_z;      ///< the z-position of the measurement
    LHCbID  m_lhcbID; ///< the corresponding LHCbID
    SubInfo m_sub;    ///< subdetector specific information

    template <typename SubI, typename = std::enable_if_t<std::is_convertible_v<SubI, SubInfo>>>
    Measurement( LHCbID id, double z, SubI&& subi ) : m_z{z}, m_lhcbID{id}, m_sub{std::forward<SubI>( subi )} {}

  public:
    /// VP constructor
#ifdef USE_DD4HEP
    Measurement( LHCbID id, double z, LHCb::LineTraj<double> traj, double errMeas, const DeVPSensor elem )
#else
    Measurement( LHCbID id, double z, LHCb::LineTraj<double> traj, double errMeas, const DeVPSensor* elem )
#endif
        : Measurement{id, z, VP{std::move( traj ), elem, errMeas, details::Measurement::isX( traj.direction( 0 ) )}} {
    }

    /// UT constructor
#ifdef USE_DD4HEP
    Measurement( LHCbID id, double z, LHCb::LineTraj<double> traj, double errMeas, const DeUTSector elem )
#else
    Measurement( LHCbID id, double z, LHCb::LineTraj<double> traj, double errMeas, const DeUTSector* elem )
#endif
        : Measurement{id, z, UT{std::move( traj ), elem, errMeas}} {
    }

    /// FT constructor
#ifdef USE_DD4HEP
    Measurement( LHCbID id, double z, LHCb::LineTraj<double> traj, double errMeas, const DeFTMat elem )
#else
    Measurement( LHCbID id, double z, LHCb::LineTraj<double> traj, double errMeas, const DeFTMat* elem )
#endif
        : Measurement{id, z, FT{std::move( traj ), elem, errMeas}} {
    }

    /// Muon constructor
#ifdef USE_DD4HEP
    Measurement( LHCbID id, double z, LHCb::LineTraj<double> traj, double errMeas, const DeMuonChamber elem )
#else
    Measurement( LHCbID id, double z, LHCb::LineTraj<double> traj, double errMeas, const DeMuonChamber* elem )
#endif
        : Measurement{id, z, Muon{std::move( traj ), elem, errMeas}} {
    }

    /// VP2D constructor
#ifdef USE_DD4HEP
    Measurement( LHCbID id, double z, Point traj, ROOT::Math::SVector<double, 2> errMeas, const DeVPSensor elem )
#else
    Measurement( LHCbID id, double z, Point traj, ROOT::Math::SVector<double, 2> errMeas, const DeVPSensor* elem )
#endif
        : Measurement{id, z, VP2D{std::move( traj ), elem, errMeas}} {
    }

    // Observers

    /// invoke the matching callable which can be invoked for the detector-specific content in this measurement
    // Same as above, but allow for an additional variant as first argument of visit
    template <typename F, typename... Fs>
    decltype( auto ) visit( F&& f, Fs&&... fs ) const {
      if constexpr ( is_variant_v<F> ) {
        return std::visit( Gaudi::overload( std::forward<Fs>( fs )... ), m_sub, std::forward<F>( f ) );
      } else {
        return std::visit( Gaudi::overload( std::forward<F>( f ), std::forward<Fs>( fs )... ), m_sub );
      }
    }

    /// Get the sub-detector specific information
    template <typename SubI>
    const SubI* getIf() const {
      return std::get_if<SubI>( &m_sub );
    }

    /// Check whether this Measurements 'is-a' specific 'sub' measurement
    template <typename SubI>
    bool is() const {
      return std::holds_alternative<SubI>( m_sub );
    }

    /// Retrieve const  the corresponding LHCbID
    LHCbID lhcbID() const { return m_lhcbID; }

    /// Retrieve const  the z-position of the measurement
    double z() const { return m_z; }

    /// This set of methods is wrapping DetectorElement ones so that user of Measurements
    /// do not have to see the DetectorElement class and thus can use the DetDesc or the
    /// Detector version the same way. Not really nice as an interface but the long term
    /// goal begin to drop Measurements, this is an easy solution
    Gaudi::XYZPoint  toLocal( const Gaudi::XYZPoint& globalPoint ) const;
    Gaudi::XYZVector toLocal( const Gaudi::XYZVector& globalPoint ) const;
    Gaudi::XYZPoint  toGlobal( const Gaudi::XYZPoint& globalPoint ) const;
    Gaudi::XYZVector toGlobal( const Gaudi::XYZVector& globalPoint ) const;
    std::string      name() const;
    bool             isSameDetectorElement( const Measurement& other ) const;

    using Type = Enum::Measurement::Type;

    template <typename SubI>
    static constexpr auto index = details::Measurement::index<SubI, SubInfo>();

    /// Retrieve the measurement type
    Type type() const {
      auto idx = m_sub.index();
      ASSUME( idx < std::variant_size_v<SubInfo> );
      switch ( idx ) {
      case index<VP>:
        return Type::VP;
      case index<UT>:
        return Type::UT;
      case index<FT>:
        return Type::FT;
      case index<Muon>:
        return Type::Muon;
      case index<VP2D>:
        return Type::VP2D;
      default:
        return Type::Unknown;
      }
    };

    /// Modifiers

    template <typename SubI>
    SubI* getIf() {
      return std::get_if<SubI>( &m_sub );
    }

  }; //< class Measurement

  struct Minimize1DResult {
    double                         zState;
    double                         sMeas;
    double                         doca;
    Gaudi::TrackProjectionMatrix1D H;
    Gaudi::XYZVector               unitPocaVector;
  };
  struct Minimize2DResult {
    double                         zState;
    double                         sMeas;
    double                         doca;
    Gaudi::TrackProjectionMatrix2D H;
    Gaudi::XYZVector               unitPocaVector;
  };

  using MinimizeResult = std::variant<Minimize1DResult, Minimize2DResult>;
  MinimizeResult minimize( const LHCb::Measurement& m, const LHCb::StateZTraj<double>& refTraj, double zState );

} // namespace LHCb
