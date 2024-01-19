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
#include "Event/KalmanFitResult.h"
#include "Event/FitNode.h"
#include "Event/Measurement.h"
#include <algorithm>

namespace LHCb {
  // copy constructor
  KalmanFitResult::KalmanFitResult( const KalmanFitResult& rhs )
      : TrackFitResult( rhs )
      , m_seedCovariance( rhs.m_seedCovariance )
      , m_nTrackParameters( rhs.m_nTrackParameters )
      , m_chi2CacheValid( rhs.m_chi2CacheValid )
      , m_errorFlag( rhs.m_errorFlag )
      , m_bidirectionalSmoother( rhs.m_bidirectionalSmoother ) {
    establishNodeLinks();
  }

  // copy from TrackFitResult
  KalmanFitResult::KalmanFitResult( const TrackFitResult& rhs )
      : TrackFitResult( rhs )
      , m_nTrackParameters( 5 )
      , m_chi2CacheValid( false )
      , m_errorFlag( 0x00 )
      , m_bidirectionalSmoother( true ) {}

  // clone
  std::unique_ptr<ITrackFitResult> KalmanFitResult::clone() const { return std::make_unique<KalmanFitResult>( *this ); }

  // set the error flag out of direction, algorithm and error type identifiers
  void KalmanFitResult::setErrorFlag( unsigned short direction, unsigned short algnum, unsigned short errnum ) {
    m_errorFlag = ( ( (unsigned short)1 ) << globalBits ) + ( ( (unsigned short)direction ) << dirBits ) +
                  ( ( (unsigned short)algnum ) << algBits ) + ( ( (unsigned short)errnum ) << typeBits );
  }

  // check the global error status of the node
  bool KalmanFitResult::inError() const { return m_errorFlag != 0; }

  // get the error description
  std::string KalmanFitResult::getError() const {
    if ( m_errorFlag == 0 ) return std::string( "No error" );
    unsigned short     direction = ( m_errorFlag & dirMask ) >> dirBits;
    unsigned short     algnum    = ( m_errorFlag & algMask ) >> algBits;
    unsigned short     errnum    = ( m_errorFlag & typeMask );
    std::ostringstream errMsg;
    std::ostringstream dir;
    // Set the direction
    switch ( direction ) {
    case Forward:
      dir << "forward ";
      break;
    case Backward:
      dir << "backward ";
      break;
    default:
      dir << "";
      break;
    }
    // Set the algorithm
    switch ( algnum ) {
    case Predict:
      errMsg << "Error in predict " << dir.str() << "function: ";
      if ( errnum == Initialization )
        errMsg << "seed covariance is not set!";
      else if ( errnum == AlgError )
        errMsg << "something goes wrong in the prediction";
      else
        errMsg << "unknown error";
      break;
    case Filter:
      errMsg << "Error in filter " << dir.str() << "function: ";
      if ( errnum == Initialization )
        errMsg << "projection matrix is not set!";
      else if ( errnum == AlgError )
        errMsg << "something goes wrong in the filtering";
      else
        errMsg << "unknown error";
      break;
    case Smooth:
      errMsg << "Error in smooth function: ";
      if ( errnum == MatrixInversion )
        errMsg << "error in matrix inversion";
      else if ( errnum == AlgError )
        errMsg << "non positive diagonal element in coveriance matrix";
      else if ( errnum == Other )
        errMsg << "problem with HCH.";
      else
        errMsg << "unknown error";
      break;
    case ComputeResidual:
      errMsg << "Error in compute residual: ";
      if ( errnum == Other )
        errMsg << " non positive variance.";
      else
        errMsg << "unknown error";
      break;
    case WeightedAverage:
      errMsg << "Error in weighted average: ";
      if ( errnum == Other )
        errMsg << " non positive variance.";
      else
        errMsg << "unknown error";
      break;
    default:
      errMsg << "Unknown error...";
      break;
    }
    return errMsg.str();
  }

  // return (chisq,dof) out of the differnet contribution
  void KalmanFitResult::computeChiSquares() {
    // This routine calculates the chisquare contributions from
    // different segments of the track. It uses the 'delta-chisquare'
    // contributions from the bi-directional kalman fit. Summing these
    // leads to a real chisquare only if the contributions are
    // uncorrelated. For a Velo-TT-T track you can then calculate:
    //
    // - the chisuare of the T segment and the T-TT segment by using the
    // 'upstream' contributions
    //
    // - the chisquare of the Velo segment and the Velo-TT segment by
    // using the 'downstream' contributions
    //
    // Note that you cannot calculate the contribution of the TT segment
    // seperately (unless there are no T or no Velo hits). Also, if
    // there are Muon hits, you cannot calculate the T station part, so
    // for now this only works for tracks without muon hits.

    // first reset everything

    if ( nodes().empty() ) {
      m_chi2 = m_chi2Muon = m_chi2Downstream = m_chi2Velo = m_chi2Upstream = m_chi2Long = ChiSquare();
    } else {
      auto fitnodes  = nodes();
      auto firstVelo = fitnodes.end();
      auto firstTT   = fitnodes.end();
      auto firstT    = fitnodes.end();
      auto firstMuon = fitnodes.end();
      auto lastVelo  = fitnodes.end();
      auto lastTT    = fitnodes.end();
      auto lastT     = fitnodes.end();
      auto lastMuon  = fitnodes.end();

      for ( auto it = fitnodes.begin(); it != fitnodes.end(); ++it ) {
        if ( ( *it )->type() != LHCb::FitNode::Type::HitOnTrack ) continue;

        ( *it )->measurement().visit( [&]( const auto& arg ) {
          using arg_t = std::decay_t<decltype( arg )>;
          if constexpr ( std::is_base_of_v<LHCb::Measurement::VP, arg_t> ) {
            if ( firstVelo == fitnodes.end() ) firstVelo = it;
            lastVelo = it;
          } else if constexpr ( std::is_base_of_v<LHCb::Measurement::VP2D, arg_t> ) {
            if ( firstVelo == fitnodes.end() ) firstVelo = it;
            lastVelo = it;
          } else if constexpr ( std::is_base_of_v<LHCb::Measurement::UT, arg_t> ) {
            if ( firstTT == fitnodes.end() ) firstTT = it;
            lastTT = it;
          } else if constexpr ( std::is_base_of_v<LHCb::Measurement::FT, arg_t> ) {
            if ( firstT == fitnodes.end() ) firstT = it;
            lastT = it;
          } else if constexpr ( std::is_base_of_v<LHCb::Measurement::Muon, arg_t> ) {
            if ( firstMuon == fitnodes.end() ) firstMuon = it;
            lastMuon = it;
          }
        } );
      }

      bool upstream = nodes().front()->z() > nodes().back()->z();

      m_chi2Muon = ( firstMuon == fitnodes.end() ? LHCb::ChiSquare{}
                                                 : upstream ? ( **lastMuon ).totalChi2( FitNode::Forward )
                                                            : ( **firstMuon ).totalChi2( FitNode::Backward ) );

      m_chi2Downstream = ( firstT == fitnodes.end() ? m_chi2Muon
                                                    : upstream ? ( **lastT ).totalChi2( FitNode::Forward )
                                                               : ( **firstT ).totalChi2( FitNode::Backward ) );

      m_chi2Velo = ( firstVelo == fitnodes.end() ? LHCb::ChiSquare{}
                                                 : upstream ? ( **firstVelo ).totalChi2( FitNode::Backward )
                                                            : ( **lastVelo ).totalChi2( FitNode::Forward ) );

      m_chi2Upstream = ( firstTT == fitnodes.end() ? m_chi2Velo
                                                   : upstream ? ( **firstTT ).totalChi2( FitNode::Backward )
                                                              : ( **lastTT ).totalChi2( FitNode::Forward ) );

      m_chi2Long = ( firstT == fitnodes.end() ? m_chi2Upstream
                                              : upstream ? ( **firstT ).totalChi2( FitNode::Backward )
                                                         : ( **lastT ).totalChi2( FitNode::Forward ) );

      const LHCb::ChiSquare& chi2A = nodes().front()->totalChi2( LHCb::FitNode::Backward );
      const LHCb::ChiSquare& chi2B = nodes().back()->totalChi2( LHCb::FitNode::Forward );
      m_chi2                       = ( chi2A.chi2() > chi2B.chi2() ) ? chi2A : chi2B;
    }
    m_chi2CacheValid = true;
  }

  /// setup the link to previous/next fitnode
  void KalmanFitResult::establishNodeLinks() {
    LHCb::FitNode* prev = nullptr;
    for ( auto& fitnode : nodes() ) {
      fitnode->setPreviousNode( prev );
      fitnode->setParent( this );
      prev = fitnode;
    }
  }

  LHCb::ChiSquare KalmanFitResult::chi2VeloTMatch() const {
    LHCb::ChiSquare velo = chi2Velo();
    LHCb::ChiSquare T    = chi2Downstream();
    LHCb::ChiSquare tot  = chi2();
    return tot - T - velo;
  }

  LHCb::ChiSquare KalmanFitResult::chi2TTHits() const {
    LHCb::ChiSquare rc;
    const auto&     fitnodes = nodes();
    auto            isUT     = []( const auto* n ) {
      return n->type() == LHCb::FitNode::Type::HitOnTrack && ( n->measurement().template is<LHCb::Measurement::UT>() );
    };
    auto                 ifnode    = std::find_if( fitnodes.begin(), fitnodes.end(), isUT );
    const LHCb::FitNode* firstnode = ( ifnode != fitnodes.end() ? *ifnode : nullptr );
    if ( firstnode ) {
      // this is guaranteed find a result prior to hitting `rend()`
      const LHCb::FitNode* lastnode = *std::find_if( fitnodes.rbegin(), fitnodes.rend(), isUT );
      LHCb::ChiSquare      chi2WithOutTTHits;
      // first add the contributions untill the TT hits
      if ( lastnode->hasInfoUpstream( LHCb::FitNode::Backward ) )
        chi2WithOutTTHits += lastnode->prevNode( LHCb::FitNode::Backward )->totalChi2( LHCb::FitNode::Backward );
      if ( firstnode->hasInfoUpstream( LHCb::FitNode::Forward ) )
        chi2WithOutTTHits += firstnode->prevNode( LHCb::FitNode::Forward )->totalChi2( LHCb::FitNode::Forward );

      // now compute the contribution of the velo-T match, if there is any
      if ( lastnode->hasInfoUpstream( LHCb::FitNode::Backward ) &&
           firstnode->hasInfoUpstream( LHCb::FitNode::Forward ) ) {
        // that means we need to take off the hits.

        // for now, let's do it the lazy way: make a clone, then flag
        // all hits as outliers. the cloning is real slow .
        auto copy = std::unique_ptr<KalmanFitResult>( static_cast<KalmanFitResult*>( this->clone().release() ) );

        // set all TT hits as ourliers
        LHCb::FitNode *copylastnode( nullptr ), *copyfirstnode( nullptr );
        for ( auto& node : copy->nodes() ) {
          if ( isUT( node ) ) {
            node->deactivateMeasurement();
            copylastnode = node;
            if ( !copyfirstnode ) copyfirstnode = node;
          }
        }
        LHCb::State         stateA = copylastnode->predictedState( LHCb::FitNode::Forward );
        LHCb::State         stateB = copylastnode->predictedState( LHCb::FitNode::Backward );
        Gaudi::SymMatrix5x5 weight = stateA.covariance() + stateB.covariance();
        weight.InvertChol();
        Gaudi::Vector5 res = stateA.stateVector() - stateB.stateVector();
        chi2WithOutTTHits += LHCb::ChiSquare( ROOT::Math::Similarity( res, weight ), 5 );
      }
      rc = chi2() - chi2WithOutTTHits;
    }
    return rc;
  }

} // namespace LHCb
