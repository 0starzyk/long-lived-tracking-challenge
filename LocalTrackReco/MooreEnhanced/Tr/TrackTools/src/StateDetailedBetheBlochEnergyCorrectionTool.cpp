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
#include "Core/Material.h"
#include "Event/State.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "TrackInterfaces/IStateCorrectionTool.h"
#include "vdt/exp.h"
#include "vdt/log.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>

using namespace Gaudi::Units;

namespace {
  // tell gcc/clang/icc that vdt_pow has no side-effects whatsoever
  double vdt_pow( const double x, const double y ) noexcept __attribute__( ( const ) );
  double vdt_pow( const double x, const double y ) noexcept { return vdt::fast_exp( y * vdt::fast_log( x ) ); }

  const auto log_2_me = std::log( 2 * electron_mass_c2 );

#ifdef USE_DD4HEP
  using InternalPtr = decltype( std::declval<MaterialPtr>().access() );
  using MaterialRef = const LHCb::Detector::Material;
#else
  using InternalPtr = const Material*;
  using MaterialRef = const Material&;
#endif

  class MaterialCache {
    struct CacheData {
      double X0, C, X1, a, m, DensityFactor, LogI;
      CacheData( MaterialRef mat )
          : X0{mat.X0()}
          , C{mat.C()}
          , X1{mat.X1()}
          , a{mat.a()}
          , m{mat.m()}
          , DensityFactor{( 30.71 * MeV * mm2 / mole ) * mat.Z() * mat.density() / mat.A()}
          , LogI{std::log( mat.I() )} {}
    };
    using Material2FactorMap = std::unordered_map<InternalPtr, CacheData>;

    Material2FactorMap           mat2factors{};
    Material2FactorMap::iterator last = end( mat2factors );

  public:
    const CacheData& operator()( MaterialPtr material ) {
#ifdef USE_DD4HEP
      assert( material.isValid() );
      InternalPtr ptr = material.access();
#else
      assert( material != nullptr );
      InternalPtr ptr = material;
#endif
      // material lookups are slightly different: the "hot spot" occupied by the
      // last material used is likely to be missed (about 85% of the time), but if
      // we can save the hash by looking at the last used element first, that's a
      // good thing
      if ( end( mat2factors ) == last || last->first != ptr ) {
        last = mat2factors.find( ptr );
        // cache material quantities:
        // - X0, C, X1, a, m
        // - const. * Z / A * density ("DensityFactor")
        // - log(I)
        // we trade nine virtual function calls plus associated memory accesses, four
        // floating point operations and a logarithm for a lookup
        if ( end( mat2factors ) == last ) {
          // coming here is the very rare exception - after the first couple of
          // events, the material cache miss rate is essentially zero
          last = mat2factors
                     .try_emplace( ptr,
#ifdef USE_DD4HEP
                                   material
#else
                                   *material
#endif
                                   )
                     .first;
        }
      }
      return last->second;
    }
  };
} // namespace

//-----------------------------------------------------------------------------
// Implementation file for class : StateDetailedBetheBlochEnergyCorrectionTool
//
// 2006-08-18 : Eduardo Rodrigues
//-----------------------------------------------------------------------------

/** @class StateDetailedBetheBlochEnergyCorrectionTool
 *
 *  This state correction tool applies a dE/dx energy loss correction
 *  with the full version of the Bethe-Bloch equation.
 *
 *  @author Stephanie Hansmann-Menzemer
 *  @date   2008-05-02
 *
 */
class StateDetailedBetheBlochEnergyCorrectionTool : public extends<GaudiTool, IStateCorrectionTool> {
public:
  /// Standard constructor
  using extends::extends;

  /// Correct a State for dE/dx energy losses with a simplified Bethe-Bloch equiaton
  void correctState( LHCb::State& state, const MaterialPtr material, std::any& cache, double wallThickness,
                     bool upstream, double mass ) const override;

  std::any createBuffer() const override { return MaterialCache(); }

private:
  // Job options
  Gaudi::Property<double> m_energyLossCorr{this, "EnergyLossFactor", 1.0}; ///< tunable energy loss correction
  Gaudi::Property<double> m_maxEnergyLoss{this, "MaximumEnergyLoss",
                                          100. * Gaudi::Units::MeV}; ///< maximum energy loss in dE/dx correction
};

// Declaration of the Tool Factory
DECLARE_COMPONENT( StateDetailedBetheBlochEnergyCorrectionTool )

//=============================================================================
// Correct a State for dE/dx energy losses
//=============================================================================
void StateDetailedBetheBlochEnergyCorrectionTool::correctState( LHCb::State& state, const MaterialPtr material,
                                                                std::any& cache, double wallThickness, bool upstream,
                                                                const double mass ) const {
  // mass zero is quite unlikely here, and we'll save a lot of computation in
  // that case
  if ( 0 == mass ) return;

  const auto& mat = std::any_cast<MaterialCache&>( cache )( material );

  // apply correction - note: for now only correct the state vector
  auto& qOverP = state.stateVector()[4];

  const auto eta2_inv = ( qOverP * mass ) * ( qOverP * mass );
  const auto x4       = -vdt::fast_log( eta2_inv );
  const auto x        = x4 / 4.606;

  double rho_2 = 0;
  if ( x > mat.X0 ) { // branch taken in about 3/4 of cases
    rho_2 = x4 - mat.C;
    if ( x < mat.X1 ) { // branch taken in about 9/10 cases
      rho_2 += mat.a * vdt_pow( mat.X1 - x, mat.m );
    }
    rho_2 *= 0.5;
  }
  // correct wall thickness for angle of incidence
  const auto tx = state.tx(), ty = state.ty();
  wallThickness *= std::sqrt( 1 + tx * tx + ty * ty );
  const auto beta2_inv = 1 + eta2_inv;
  auto       eLoss =
      m_energyLossCorr * wallThickness * mat.DensityFactor * ( beta2_inv * ( log_2_me + x4 - mat.LogI - rho_2 ) - 1 );

  eLoss = std::copysign( std::min( m_maxEnergyLoss.value(), eLoss ), ( 2 * int( upstream ) - 1 ) * qOverP );

  qOverP /= 1 + eLoss * qOverP;
}
