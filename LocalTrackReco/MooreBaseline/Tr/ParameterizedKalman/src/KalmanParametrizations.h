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

#pragma once

// Include files
#include "GaudiKernel/GenericMatrixTypes.h"
#include "GaudiKernel/GenericVectorTypes.h"
#include "GaudiKernel/SymmetricMatrixTypes.h"
#include "Kernel/STLExtensions.h"
#include "Kernel/Trajectory.h"
#include "Math/Vector3D.h"
#include <stdio.h>

#include "KalmanParametrizationsCoef.h"

enum class Polarity { Up, Down };

// Dimensions of the parameters {n, m} : n sets with m parameters
constexpr std::pair<unsigned int, unsigned int> ParDim_predictV     = {2, 10};
constexpr std::pair<unsigned int, unsigned int> ParDim_predictVUT   = {2, 30};
constexpr std::pair<unsigned int, unsigned int> ParDim_predictUT    = {7, 20};
constexpr std::pair<unsigned int, unsigned int> ParDim_predictUTFUT = {1, 1};
constexpr std::pair<unsigned int, unsigned int> ParDim_predictUTTF  = {2, 20};
constexpr std::pair<unsigned int, unsigned int> ParDim_predictTFT   = {2, 20};
constexpr std::pair<unsigned int, unsigned int> ParDim_predictT     = {46, 20};

constexpr std::pair<unsigned int, unsigned int> ParDim_TLayer  = {2, 12};
constexpr std::pair<unsigned int, unsigned int> ParDim_UTLayer = {1, 4};

/** @class KalmanParametrizations KalmanParametrizations.h
 *  Contains a set of extrapolation methods that are used in the ParameterizedKalmanFit
 *  It uses parametrizations for the extrapolation through material and the magnetic field
 *
 *
 *  @author Simon Stemmle
 *  @date   2017-10-26
 */

template <std::size_t BATCHNUMBER, std::size_t BATCHSIZE>
struct KalmanParameters {
  static const unsigned int batchN = BATCHNUMBER;
  static const unsigned int batchS = BATCHSIZE;

  std::array<double, BATCHSIZE * BATCHNUMBER> parameters;

  double&                  operator()( int i, int j ) { return parameters[i * BATCHSIZE + j]; }
  LHCb::span<const double> operator[]( int i ) const { return {&parameters[i * BATCHSIZE], BATCHSIZE}; }
};

class KalmanParametrizations {
public:
  using FTYPE = double;

  /// Constructor from parameters
  KalmanParametrizations( const std::string& ParamFileLocation, const Polarity polarity,
                          bool useOneParameterSet = false ) {
    SetParameters( ParamFileLocation, polarity, useOneParameterSet );
  }

  /// Extrapolate inside the VELO
  void ExtrapolateInV( double zFrom, double zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// Extrapolate VELO <-> UT
  bool ExtrapolateVUT( double zFrom, double zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// Extrapolate VELO <-> UT (traj)
  bool ExtrapolateVUT( double zFrom, ROOT::Math::XYZPoint point, ROOT::Math::XYZVector dir, double& zTo,
                       Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// Get noise for VELO <- UT
  void GetNoiseVUTBackw( double zFrom, double zTo, const Gaudi::Vector5& x, Gaudi::SymMatrix5x5& Q ) const;

  /// Predict UT <-> UT (traj)
  void ExtrapolateInUT( double zFrom, int nLayer, ROOT::Math::XYZPoint point, ROOT::Math::XYZVector dir, double& zTo,
                        Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// Predict UT <-> UT
  void ExtrapolateInUT( double zFrom, int nLayer, double zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                        Gaudi::SymMatrix5x5& Q ) const;

  /// Predict around the last UT layer to the last UT hit
  void ExtrapolateUTFUT( double& zFrom, ROOT::Math::XYZPoint point, ROOT::Math::XYZVector dir, Gaudi::Vector5& x,
                         Gaudi::Matrix5x5& F ) const;

  /// Predict to the start position of the UTTF extrapolation
  void ExtrapolateUTFUTDef( double& zFrom, Gaudi::Vector5& x, Gaudi::Matrix5x5& F ) const;

  /// Predict around the last UT layer
  void ExtrapolateUTFUT( double zFrom, double zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F ) const;

  /// Extrapolate UT (fixed z) -> T (fixed z)
  void ExtrapolateUTT( Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// actual parametrization for the step UT (fixed z) -> T (fixed z)
  int extrapUTT( double zi, double zf, int quad_interp, double& x, double& y, double& tx, double& ty, double qOp,
                 double* der_tx, double* der_ty, double* der_qop ) const;

  /// Get noise for UT (fixed z) <- T (fixed z)
  void GetNoiseUTTBackw( const Gaudi::Vector5& x, Gaudi::SymMatrix5x5& Q ) const;

  /// Extrapolate T <-> T (traj)
  void ExtrapolateInT( double zFrom, int nLayer, ROOT::Math::XYZPoint point, ROOT::Math::XYZVector dir, double& zTo,
                       Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// Extrapolate T <-> T (no hit)
  void ExtrapolateInT( double zFrom, int nLayer, double& zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                       Gaudi::SymMatrix5x5& Q ) const;

  /// Extrapolate T <-> T
  void ExtrapolateInT( double zFrom, int nLayer, double zTo, double DzDy, double DzDty, Gaudi::Vector5& x,
                       Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// Extrapolate T(fixed z) <-> first T layer (traj)
  void ExtrapolateTFT( double zFrom, ROOT::Math::XYZPoint point, ROOT::Math::XYZVector dir, double& zTo,
                       Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// Extrapolate T(fixed z) <-> first T layer (no hit)
  void ExtrapolateTFTDef( double zFrom, double& zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F,
                          Gaudi::SymMatrix5x5& Q ) const;

  /// Extrapolate T(fixed z) <-> first T layer
  void ExtrapolateTFT( double zFrom, double zTo, Gaudi::Vector5& x, Gaudi::Matrix5x5& F, Gaudi::SymMatrix5x5& Q ) const;

  /// Set the parameters (if needed) for the given magnet polarity
  void SetParameters( const std::string& ParamFileLocation, const Polarity polarity, bool useOneParameterSet = false );

  /// Returns the z end-position of the extrapolation through the magnet
  double UTTExtrEndZ() const;

  /// Returns the z begin-position of the extrapolation through the magnet
  double UTTExtrBeginZ() const;

  /// Returns the default z position of the first UT layer
  double VUTExtrEndZ() const;

private:
  ///  Read extrapolation parameters from file
  template <std::size_t BN, std::size_t BS>
  void read_params( std::string_view file, KalmanParameters<BN, BS>& params );

  ///  Read extrapolation parameters from file - here for the UT -> T extrapolation
  void read_params_UTT( std::string_view file );

  /// This switches all parameters that linearly/cubicly/... depend on q/p
  template <std::size_t SIZE, std::size_t BN, std::size_t BS>
  void SwitchParamsForPolarity( KalmanParameters<BN, BS>& params, const std::array<unsigned int, SIZE> list );

  ///  set hard coded parameters fro UTT extrapolation
  void setUTTParameters( Polarity polarity );

  // Tracks the magnet polarity, the current parameters correspond to.
  Polarity m_Polarity   = Polarity::Up;
  bool     paramsLoaded = false;

  //###########################
  // parameter vectors
  //###########################

  // Parameters that change sign under polarity flip
  const std::array<unsigned int, 1> flip_Par_predictV     = {4};
  const std::array<unsigned int, 4> flip_Par_predictVUT   = {0, 8, 9, 10};
  const std::array<unsigned int, 4> flip_Par_predictUT    = {5, 6, 7, 10};
  const std::array<unsigned int, 2> flip_Par_predictUTFUT = {0};
  const std::array<unsigned int, 8> flip_Par_predictUTTF  = {0, 2, 3, 5, 6, 8, 9, 11};
  const std::array<unsigned int, 2> flip_Par_predictTFT   = {5, 6};
  const std::array<unsigned int, 3> flip_Par_predictT     = {5, 6, 7};

  bool m_qop_flip = false;

  // predict params
  KalmanParameters<ParDim_predictV.first, ParDim_predictV.second>         Par_predictV;
  KalmanParameters<ParDim_predictVUT.first, ParDim_predictVUT.second>     Par_predictVUT;
  KalmanParameters<ParDim_predictUT.first, ParDim_predictUT.second>       Par_predictUT;
  KalmanParameters<ParDim_predictUTFUT.first, ParDim_predictUTFUT.second> Par_predictUTFUT;
  KalmanParameters<ParDim_predictUTTF.first, ParDim_predictUTTF.second>   Par_predictUTTF;
  KalmanParameters<ParDim_predictTFT.first, ParDim_predictTFT.second>     Par_predictTFT;
  KalmanParameters<ParDim_predictT.first, ParDim_predictT.second>         Par_predictT;
  KalmanParameters<ParDim_TLayer.first, ParDim_TLayer.second>             Par_TLayer;
  KalmanParameters<ParDim_UTLayer.first, ParDim_UTLayer.second>           Par_UTLayer;

  // parameters for pierres method
  // they are initialized in the constructor depending on the magnet polarity

#define NBINXMAX 60
#define NBINYMAX 50
#define QUADRATICINTERPOLATION 1
  std::array<std::array<KalmanParametrizationsCoef, NBINYMAX>, NBINXMAX> C;
  double                                                                 ZINI{0}, ZFIN{0};
  double                                                                 PMIN{0};
  double BENDX{0}, BENDX_X2{0}, BENDX_Y2{0}, BENDY_XY{0};
  double Txmax{0}, Tymax{0}, XFmax{0}, Xmax{0}, Ymax{0};
  double Dtxy{0};
  double step{0};

  int Nbinx{0}, Nbiny{0};
  int XGridOption{0}, YGridOption{0};
  int DEGX1{0}, DEGX2{0}, DEGY1{0}, DEGY2{0};
};
