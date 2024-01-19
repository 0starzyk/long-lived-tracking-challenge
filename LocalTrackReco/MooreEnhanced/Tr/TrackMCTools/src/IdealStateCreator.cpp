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
#include "Event/LinksByKey.h"
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Event/MCVertex.h"
#include "Event/State.h"
#include "Event/StateVector.h"
#include "Event/Track.h"
#include "Event/TrackParameters.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Kernel/ILHCbMagnetSvc.h"
#include "MCInterfaces/IIdealStateCreator.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include <cmath>
#include <string>
#include <vector>

namespace {
  double qOverP( const LHCb::MCParticle* mcPart ) {
    const double charge = ( mcPart->particleID().threeCharge() ) / 3.;
    const double p      = mcPart->p();
    return ( p < TrackParameters::lowTolerance ? 0.0 : charge / p );
  }

  double qOverP( const LHCb::MCHit* mcHit ) {
    const double charge = ( mcHit->mcParticle()->particleID().threeCharge() ) / 3.;
    const double p      = mcHit->p();
    return ( p < TrackParameters::lowTolerance ? qOverP( mcHit->mcParticle() ) : charge / p );
  }

} // namespace

/** @class IdealStateCreator IdealStateCreator.h "TrackMCTools/IdealStateCreator.h"
 *
 *  An IdealStateCreator is an IIdealStateCreator tool that creates a
 *  State. There are two methods: one creates a state at a certain
 *  z-position using the closest two extry/exit points from a MCParticle,
 *  and a track extrapolator. The other creates a state at the vertex, using
 *  the properties of the MCParticles and MCVertex.
 *
 *  WARNING: This does NOT return an "ideal" state for arbitrary z-positions.
 *           As said above, the tool finds the closest MCHit and extrapolates
 *           from there to the z-position you wish to create a state at.
 *           Depending on the traversed material and the distance that was extrapolated,
 *           this can lead to substantial differences between the created state
 *           and the actual true position of the MCParticle at said z-position in the Simulation.
 *           But unfortunatley we don't have access to this true information for arbitraty Z-positions.
 *           The returned covariance of the state can be used as indicator
 *           to determine how likely the created state is significantly wrong.
 *
 *
 *  The diagonal elements of the covariance matrix are set with
 *  the job options. Note that "eP" is dp/p.
 *
 *  Moved to LHCb v20r0. Adapted code to use updated Det packages.
 *  @author Edwin Bos
 *  @date   2006-02-02
 *
 *  @author Eduardo Rodrigues (adaptations to new track event model)
 *  @date   2005-04-06
 *
 *  @author Rutger van der Eijk, Jeroen van Tilburg
 *  @date   3-7-2002
 */
class IdealStateCreator : public extends<GaudiTool, IIdealStateCreator> {
public:
  /// Standard constructor
  using extends::extends;

  Gaudi::Property<float> m_eX2{this, "ErrorX2", 0};
  Gaudi::Property<float> m_eY2{this, "ErrorY2", 0};
  Gaudi::Property<float> m_eTx2{this, "ErrorTx2", 0};
  Gaudi::Property<float> m_eTy2{this, "ErrorTy2", 0};
  Gaudi::Property<float> m_eP{this, "ErrorP", 0};
  Gaudi::Property<bool>  m_correctSlopes{this, "CorrectSlopes", false};

private:
  ServiceHandle<ILHCbMagnetSvc>  m_magSvc{this, "MagneticFieldService", "MagneticFieldSvc"};
  ToolHandle<ITrackExtrapolator> m_extrapolator{this, "Extrapolator", "TrackRungeKuttaExtrapolator"};

  DataObjectReadHandle<LHCb::MCVertices> m_mc_vertices{this, "MCVertices", LHCb::MCVertexLocation::Default};
  DataObjectReadHandle<LHCb::MCHits>     m_mc_vp_hits{this, "VPMCHits", LHCb::MCHitLocation::VP};
  DataObjectReadHandle<LHCb::LinksByKey> m_mc_vp_hits_links{this, "VPMCHitLinks", ""};
  DataObjectReadHandle<LHCb::MCHits>     m_mc_ut_hits{this, "UTMCHits", LHCb::MCHitLocation::UT};
  DataObjectReadHandle<LHCb::LinksByKey> m_mc_ut_hits_links{this, "UTMCHitLinks", ""};
  DataObjectReadHandle<LHCb::MCHits>     m_mc_ft_hits{this, "FTMCHits", LHCb::MCHitLocation::FT};
  DataObjectReadHandle<LHCb::LinksByKey> m_mc_ft_hits_links{this, "FTMCHitLinks", ""};

  //=============================================================================
  // Correct slopes for magnetic field given an MCHit and a MCParticle
  //=============================================================================
  void correctSlopes( LHCb::MCHit const* mcHit, double& tx, double& ty ) const {
    // TODO: I hope this method can be removed as soon as the displacement vector
    // in the MCHit is calculated in Gauss using the momentum direction of the
    // *entry point*. (JvT: 27/10/2006).
    using namespace Gaudi::Units;

    // Get magnetic field vector
    Gaudi::XYZVector B;
    m_magSvc->fieldVector( mcHit->midPoint(), B ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

    // Calculate new displacement vector and tx,ty slopes
    Gaudi::XYZVector d    = mcHit->displacement();
    Gaudi::XYZVector dNew = d - ( 0.5 * d.R() * qOverP( mcHit ) * d.Cross( B ) * eplus * c_light );
    tx                    = dNew.x() / dNew.z();
    ty                    = dNew.y() / dNew.z();
  }

public:
  /** This method creates a state at z position from a MCParticle
   *  using the entry/exit points of the MCHits.
   *  @return StatusCode
   *  @param  mcPart The MCParticle from which the state will be created
   *  @param  zRec   The z-position at which the state will be created
   *  @param  state The ref to the State which is created.
   */
  StatusCode createState( LHCb::MCParticle const* mcPart, double zRec, LHCb::State& state,
                          IGeometryInfo const& geometry ) const override {

    // Check if MCParticle exists
    if ( mcPart == nullptr ) return StatusCode::FAILURE;

    // Get the closest MCHit
    LHCb::MCHit const* closestHit = findClosestHit( mcPart, zRec );
    if ( !closestHit ) return Error( "No closest MCHit found!!" );

    return createState( closestHit, zRec, state, geometry );
  }

  /** This method creates a state at z position from a MCParticle
   *  using the entry/exit points of the MCHits.
   *  @return StatusCode
   *  @param  mcPart The MCParticle from which the state will be created
   *  @param  zRec   The z-position at which the state will be created
   *  @param  state  The StateVector which is created.
   */
  StatusCode createStateVector( LHCb::MCParticle const* mcPart, double zRec, LHCb::StateVector& state,
                                IGeometryInfo const& geometry ) const override {

    // Check if MCParticle exists
    if ( mcPart == 0 ) return StatusCode::FAILURE;

    // Get the closest MCHit
    LHCb::MCHit const* closestHit = findClosestHit( mcPart, zRec );
    if ( !closestHit ) return Error( "No closest MCHit found!!" );

    return createStateVector( closestHit, zRec, state, geometry );
  }

  /** This method creates a state at z position from a MCHit
   *  using the entry/exit points of the MCHit.
   *  @return StatusCode
   *  @param  aHit The MCHit from which the state will be created
   *  @param  zRec   The z-position at which the state will be created
   *  @param  state The ref to the State which is created.
   */
  StatusCode createState( LHCb::MCHit const* aHit, double zRec, LHCb::State& state,
                          IGeometryInfo const& geometry ) const override {

    LHCb::StateVector pVec;
    createStateVector( aHit, pVec );

    state.setState( pVec.parameters() );
    state.setZ( pVec.z() );

    // set covariance matrix
    auto cov    = Gaudi::TrackSymMatrix();
    cov( 0, 0 ) = m_eX2;
    cov( 1, 1 ) = m_eY2;
    cov( 2, 2 ) = m_eTx2;
    cov( 3, 3 ) = m_eTy2;
    cov( 4, 4 ) = std::pow( m_eP * state.qOverP(), 2 );
    state.setCovariance( cov );

    // transport to the z we want
    StatusCode sc = m_extrapolator->propagate( state, zRec, geometry );
    if ( sc.isFailure() ) {
      warning() << "Extrapolation of True State from z = " << state.z() << " to z = " << zRec << " failed!" << endmsg;
    }

    return sc;
  }

  /** This method creates a state at z position from a MCHit
   *  using the entry/exit points of the MCHit.
   *  @return StatusCode
   *  @param  aHit The MCHit from which the state will be created
   *  @param  zRec  The z-position at which the state will be created
   *  @param  pVec  The StateVector which is created.
   */
  StatusCode createStateVector( LHCb::MCHit const* aHit, double zRec, LHCb::StateVector& pVec,
                                IGeometryInfo const& geometry ) const override {

    // createState vector
    createStateVector( aHit, pVec );

    // extrapolate state to exact z position
    //  TrackVector& trackVec = pVec.parameters();
    // StatusCode sc = m_extrapolator->propagate(trackVec,pVec.z() ,zRec );
    StatusCode sc = m_extrapolator->propagate( pVec, zRec, geometry );
    if ( sc.isFailure() ) {
      warning() << "Extrapolation of True State from z = " << pVec.z() << " to z = " << zRec << " failed!" << endmsg;
    }
    //  pVec.setZ(zRec);

    return StatusCode::SUCCESS;
  }

  /** This method creates a state at the z position of the MCHit
   *  using the entry point of the MCHit.
   *  @return StatusCode
   *  @param  aHit The MCHit from which the state will be created
   *  @param  pVec The ref to the StateVector which is created.
   */

  void createStateVector( LHCb::MCHit const* aHit, LHCb::StateVector& pVec ) const {
    // Correct tx and ty from the MCHit for the magnetic field
    double tx = aHit->dxdz();
    double ty = aHit->dydz();
    if ( m_correctSlopes ) correctSlopes( aHit, tx, ty );
    Gaudi::XYZVector direction( tx, ty, 1.0 );

    // determine Q/P
    const double trueQOverP = qOverP( aHit );

    // construct true State
    pVec = LHCb::StateVector( aHit->entry(), direction, trueQOverP );
  }

  /** This method creates a state at the origin vertex from a MCParticle
   *  using the entry/exit points of the MCHits.
   *  @return StatusCode
   *  @param  mcParticle The MCParticle from which the state will be created
   *  @param  state The ref to the State which is created.
   */
  StatusCode createStateVertex( LHCb::MCParticle const* mcParticle, LHCb::State& state ) const override {

    // Check if MCParticle exists
    if ( mcParticle == 0 ) return StatusCode::FAILURE;

    LHCb::StateVector pVec;
    StatusCode        sc = createStateVectorVertex( mcParticle, pVec );
    if ( sc.isFailure() ) { return Warning( "Failed to create state vector", StatusCode::SUCCESS ); }

    state.setZ( pVec.z() );
    state.setState( pVec.parameters() );

    // set covariance matrix
    auto cov    = Gaudi::TrackSymMatrix();
    cov( 0, 0 ) = m_eX2;
    cov( 1, 1 ) = m_eY2;
    cov( 2, 2 ) = m_eTx2;
    cov( 3, 3 ) = m_eTy2;
    cov( 4, 4 ) = std::pow( m_eP * state.qOverP(), 2 );
    state.setCovariance( cov );

    return StatusCode::SUCCESS;
  }

  /** This method creates a state at the origin vertex from a MCParticle
   *  using the entry/exit points of the MCHits.
   *  @return StatusCode
   *  @param  mcParticle The MCParticle from which the state will be created
   *  @param  pVec the StateVector which is created.
   */
  StatusCode createStateVectorVertex( LHCb::MCParticle const* mcParticle, LHCb::StateVector& pVec ) const override {

    // Check if MCParticle exists
    if ( mcParticle == 0 ) return StatusCode::FAILURE;

    // retrieve true MC particle info
    auto const* mcVertex = mcParticle->originVertex();
    auto const& mcPos    = mcVertex->position();
    auto const& mc4Mom   = mcParticle->momentum();

    // determine Q/P
    const double trueQOverP = qOverP( mcParticle );

    // construct true State
    pVec = LHCb::StateVector( mcPos, Gaudi::XYZVector( mc4Mom ), trueQOverP );

    return StatusCode::SUCCESS;
  }

  //=============================================================================
  // Find the z-closest MCHit associated to an MCParticle
  // looping over the hits in all the tracking detectors
  //=============================================================================
  LHCb::MCHit const* findClosestHit( LHCb::MCParticle const* mcPart, const double zRec ) const {
    LHCb::MCHit const* closestHit = nullptr;
    double             closestZ   = 1000000.0;

    auto make_closest_hit_lambda = [&closestZ, &closestHit, zRec]( auto const* mc_hit_container ) {
      return
          [&closestZ, &closestHit, zRec, mc_hit_container]( unsigned /*srcIdx*/, unsigned tgtIdx, float /*weight*/ ) {
            auto const* mchit  = static_cast<LHCb::MCHit const*>( mc_hit_container->containedObject( tgtIdx ) );
            auto const  deltaZ = std::abs( mchit->midPoint().z() - zRec );
            if ( deltaZ < closestZ || closestHit == nullptr ) {
              closestHit = mchit;
              closestZ   = deltaZ;
            }
          };
    };

    m_mc_vp_hits_links.get()->applyToLinks( mcPart->index(), make_closest_hit_lambda( m_mc_vp_hits.get() ) );
    // m_mc_vp_hits_links.get()->applyToLinks( mcPart->index(), get_closest_hit_lambda );
    // If we are already pretty close we can't get closer in a different detector
    // TODO what is the biggest save value here?
    if ( closestZ < 100 ) return closestHit;
    m_mc_ut_hits_links.get()->applyToLinks( mcPart->index(), make_closest_hit_lambda( m_mc_ut_hits.get() ) );
    if ( closestZ < 100 ) return closestHit;
    m_mc_ft_hits_links.get()->applyToLinks( mcPart->index(), make_closest_hit_lambda( m_mc_ft_hits.get() ) );
    return closestHit;
  }

  /** Return a vector containing a LHCb::State for each MCHits related to MCParticle
   *  @return StatusCode
   *  @param  mcPart The MCParticle for which the MCHits will be added
   *  @param  states vector used as return paramter of the states.
   */
  StatusCode getMCHitStates( const LHCb::MCParticle& mcPart, std::vector<LHCb::State>& states ) const override {
    // adds states at all MC hit positions

    auto make_add_states_lambda = [&states, this]( auto const* mc_hit_container ) {
      return [&states, mc_hit_container, this]( unsigned /*srcIdx*/, unsigned tgtIdx, float /*weight*/ ) {
        auto const* mchit = static_cast<LHCb::MCHit const*>( mc_hit_container->containedObject( tgtIdx ) );
        if ( !mchit ) {
          error() << "MCHit is nullptr, this should not happen! Please make sure you are using the correct LinksByKey!!"
                  << endmsg;
          error() << "A MCParticle -> MCHits LinksByKey Object is expected" << endmsg;
          return;
        }
        LHCb::StateVector statevec;
        createStateVector( mchit, statevec );
        states.emplace_back( statevec );
      };
    };

    m_mc_vp_hits_links.get()->applyToLinks( mcPart.index(), make_add_states_lambda( m_mc_vp_hits.get() ) );
    m_mc_ut_hits_links.get()->applyToLinks( mcPart.index(), make_add_states_lambda( m_mc_ut_hits.get() ) );
    m_mc_ft_hits_links.get()->applyToLinks( mcPart.index(), make_add_states_lambda( m_mc_ft_hits.get() ) );

    return StatusCode::SUCCESS;
  }
};

DECLARE_COMPONENT( IdealStateCreator )
