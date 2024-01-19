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
// Rec
#include "Event/PrFitNode.h"
// LHCb
#include "Event/FitNode.h"
#include "Event/TrackTypes.h"

/**
 * This file contains the implementation for a ``simple'' per-node
 * scattering parametrisation, intended to be used in the
 * TrackMasterFitter.
 *
 * Information about its implementation can be found here:
 * https://surfdrive.surf.nl/files/index.php/s/B4pOFnL1rDRFfvO
 *
 * @author Wouter Hulsbergen
 **/
namespace TrackFit {
  namespace param_scatter_impl {
    enum NodeType {
      ClosestToBeam,
      VPHit,
      EndVelo,
      BegRich1,
      EndRich1,
      UTHit,
      AtT,
      FTHit,
      BegRich2,
      EndRich2,
      MuonHit,
      RefNode,
      HitNode,
      NTypes
    };

    struct ParametrisedScatter {
      float Q{0};        // effective scattering thickness in tx (usually in range 1-5)
      float etaxx{0.5};  // always in [0,1]
      float etaxtx{0.5}; // always in [0,1]
      float eloss{1.0};  // energy loss in MeV
    };

    using ParametrisedScatters = std::map<int, std::map<int, ParametrisedScatter>>;

    inline NodeType nodetype( LHCb::Pr::Tracks::Fit::Node const& node ) {
      switch ( node.type() ) {
      case LHCb::Pr::Tracks::Fit::Node::Type::VPHit:
        return NodeType::VPHit;
      case LHCb::Pr::Tracks::Fit::Node::Type::UTHit:
        return NodeType::UTHit;
      case LHCb::Pr::Tracks::Fit::Node::Type::FTHit:
        return NodeType::FTHit;
      case LHCb::Pr::Tracks::Fit::Node::Type::MuonHit:
        return NodeType::MuonHit;
      case LHCb::Pr::Tracks::Fit::Node::Type::BegRich1:
        return NodeType::BegRich1;
      case LHCb::Pr::Tracks::Fit::Node::Type::EndRich1:
        return NodeType::EndRich1;
      default:
        std::cout << "Unknown node type : " << node.z() << std::endl;
        return NodeType::HitNode;
      }
    }

    inline NodeType nodetype( const LHCb::FitNode& node ) {
      NodeType rc = NTypes;
      if ( node.hasMeasurement() ) {
        rc = node.measurement().visit(
            []( LHCb::Measurement::FT const& ) { return FTHit; }, []( LHCb::Measurement::VP const& ) { return VPHit; },
            []( LHCb::Measurement::VP2D const& ) { return VPHit; },
            []( LHCb::Measurement::UT const& ) { return UTHit; },
            []( LHCb::Measurement::Muon const& ) { return MuonHit; },
            [&]( ... ) {
              std::cout << "Unknown node type for measurement: " << node.measurement().type() << std::endl;
              return HitNode;
            } );

      } else {
        switch ( node.location() ) {
        case LHCb::State::ClosestToBeam:
          rc = ClosestToBeam;
          break;
        case LHCb::State::EndVelo:
          rc = EndVelo;
          break;
        case LHCb::State::BegRich1:
          rc = BegRich1;
          break;
        case LHCb::State::EndRich1:
          rc = EndRich1;
          break;
        case LHCb::State::BegRich2:
          rc = BegRich2;
          break;
        case LHCb::State::EndRich2:
          rc = EndRich2;
          break;
        case LHCb::State::AtT:
          rc = AtT;
          break;
        default:
          rc = RefNode;
          std::cout << "Unknown reference node type : " << node.z() << std::endl;
        }
      }
      return rc;
    }

    inline ParametrisedScatters fillParametrisedScatters() {
      ParametrisedScatters rc;

      // tuned on 10k events of upgrade-magdown-sim10-up08-30000000-digi from TestFileDB
      rc[VPHit][ClosestToBeam] = {2.91, 0.808, 0.793, 1.29};
      rc[VPHit][VPHit]         = {1.48, 0.643, 0.526, 0.592};
      rc[EndVelo][VPHit]       = {4.2, 0.294, 0.233, 1.6};
      rc[BegRich1][EndVelo]    = {7.52, 0.495, 0.469, 2.53};
      rc[EndRich1][BegRich1]   = {5.71, 0.576, 0.496, 3.46};
      rc[UTHit][EndRich1]      = {2.62, 0.722, 0.689, 1.81};
      rc[UTHit][UTHit]         = {1.2, 0.634, 0.503, 0.735};
      rc[AtT][EndRich1]        = {74.6, 0.249, 0.175, 26.5};
      rc[AtT][UTHit]           = {5.99, 0.419, 0.283, 3.49};
      rc[FTHit][AtT]           = {0.676, 0.909, 0.884, 0.594};
      rc[FTHit][FTHit]         = {1.17, 0.62, 0.499, 1.05};
      rc[BegRich2][FTHit]      = {0.501, 0.313, 0.253, 0.504};
      rc[BegRich2][AtT]        = {14.7, 0.636, 0.568, 13.8};
      rc[EndRich2][BegRich2]   = {22.4, 0.696, 0.611, 10.8};
      rc[MuonHit][EndRich2]    = {1.58e+04, 0.571, 0.539, 1.76e+03};
      rc[MuonHit][MuonHit]     = {1.04e+04, 0.543, 0.506, 120};
      rc[MuonHit][FTHit]       = {1.58e+04, 0.743, 0.735, 1.77e+03};

      // needed for the state creation after fit in PrKalmanFilter
      rc[BegRich1][VPHit] = {13.1, 0.535, 0.478, 4.49};
      rc[FTHit][EndRich1] = {75, 0.298, 0.201, 27.1};

      // simple fitter doesn't use reference states so we need these transitions
      rc[UTHit][VPHit] = {23.3, 0.482, 0.357, 10.5};
      rc[FTHit][VPHit] = {99.3, 0.306, 0.244, 35.1};
      rc[FTHit][UTHit] = {6.97, 0.53, 0.383, 4.18};

      // fill the inverse table as well. I hope that I have the formulas right
      for ( int i = ClosestToBeam; i < NTypes; ++i )
        for ( int j = ClosestToBeam; j < i; ++j ) {
          const ParametrisedScatter s1 = rc[i][j];
          ParametrisedScatter       s2 = s1;
          s2.etaxtx                    = s1.etaxtx - 1; // take into account sign change in dz
          s2.etaxx                     = std::sqrt( 1 + s2.etaxx * s2.etaxx - 2 * s1.etaxtx );
          rc[j][i]                     = s2;
        }

      return rc;
    }

    inline std::pair<Gaudi::TrackSymMatrix, double>
    computeNoiseAndDeltaE( NodeType thisnodetype, double const x, double const y, double const z, double const tx,
                           double const ty, NodeType prevnodetype, double const prevnode_z, double const pscatter ) {
      static const ParametrisedScatters scatters = fillParametrisedScatters();
      Gaudi::TrackSymMatrix             Q;
      double                            deltaE{0};
      auto const                        dz = z - prevnode_z;

      // this cna be optimized later
      if ( std::abs( dz ) <= 0.5 ) return {Q, deltaE};

      // let's complete the correction for a thin scatterer, then use
      // that to take the known effects out of the noise.
      const auto tx2  = tx * tx;
      const auto ty2  = ty * ty;
      const auto n2   = 1 + tx2 + ty2;
      const auto n    = std::sqrt( n2 );
      const auto invp = 1 / pscatter;
      const auto norm = n2 * invp * invp * n; // I believe that LHCb tools are missing the last factor n

      const ParametrisedScatter scatter = scatters.at( prevnodetype ).at( thisnodetype );

      auto normCms = norm * scatter.Q;
      // add a bit extra for tracks inside the rf foil. need to find a more efficient way to do this.
      if ( prevnodetype == VPHit && thisnodetype == VPHit ) {
        const auto xprime = y + x;
        const auto yprime = y - x;
        bool infoil = ( yprime > -15 && xprime >= 0 && xprime < 15 ) || ( yprime < 15 && xprime > -15 && xprime <= 0 );
        const float rffoilscatter = 0.6;
        if ( infoil ) normCms += norm * rffoilscatter;
      }
      // else if( prevnodetype == EndVelo && thisnodetype == VPHit ) {
      // 	const auto x = state.x() ;
      // 	const auto y = state.y() ;
      // 	const auto R2 = x*x+y*y ;
      // 	if( R2 >28*28 )  normCms += norm * 9. ;
      // } else if( prevnodetype == BegRich1 && thisnodetype == EndVelo ) {
      // 	//const auto x = state.x() ;
      // 	//const auto y = state.y() ;
      // 	//const auto R2 = x*x+y*y ;
      // 	//if( R2 <40*40 )  normCms += norm * 4. ;
      // 	// I don;t think
      // 	if( t2 < 0.06*0.06 ) normCms += norm * 10 ;
      // }
      else if ( prevnodetype == AtT && thisnodetype == UTHit ) {
        if ( tx2 + ty2 < 0.02 * 0.02 ) normCms += norm * 10;
      }

      Q( 2, 2 ) = ( 1 + tx2 ) * normCms;
      Q( 3, 3 ) = ( 1 + ty2 ) * normCms;
      Q( 3, 2 ) = tx * ty * normCms;

      // x,tx part
      Q( 0, 0 ) = Q( 2, 2 ) * dz * dz * scatter.etaxx * scatter.etaxx;
      Q( 2, 0 ) = Q( 2, 2 ) * dz * scatter.etaxtx;
      // y,ty part
      Q( 1, 1 ) = Q( 3, 3 ) * dz * dz * scatter.etaxx * scatter.etaxx;
      Q( 3, 1 ) = Q( 3, 3 ) * dz * scatter.etaxtx;
      // x,y part
      Q( 1, 0 ) = Q( 3, 2 ) * dz * dz * scatter.etaxx * scatter.etaxx;
      Q( 3, 0 ) = Q( 2, 1 ) = Q( 3, 2 ) * dz * scatter.etaxtx;

      // energyloss part
      deltaE = ( dz < 0 ? 1 : -1 ) * scatter.eloss * n;
      // the landau distribution is wide: assign full Eloss as error in cov matrix
      // (since we add many small contributions, it will be small effect in the end)
      // const auto qop = state.qOverP() ;
      // const auto sigmaQOverP = qop*qop*deltaE ;
      // Q(4,4) += sigmaQOverP*sigmaQOverP ;
      return std::make_pair( Q, deltaE );
    }

    template <typename FitNodeType>
    std::pair<Gaudi::TrackSymMatrix, double> computeNoiseAndDeltaE( FitNodeType& prevnode, FitNodeType& node,
                                                                    double pscatter ) {
      const auto& state = node.refVector();
      return computeNoiseAndDeltaE( nodetype( node ), state.x(), state.y(), node.z(), state.tx(), state.ty(),
                                    nodetype( prevnode ), prevnode.z(), pscatter );
    }
  } // namespace param_scatter_impl
} // namespace TrackFit
