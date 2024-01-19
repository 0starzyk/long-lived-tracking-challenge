/*****************************************************************************\
* (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
#include "DetDesc/DetectorElement.h"
#include "DetDesc/GenericConditionAccessorHolder.h"
#include "Event/MCHit.h"
#include "Event/MCParticle.h"
#include "Event/Track.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "LHCbAlgs/Consumer.h"
#include "Linker/LinkedFrom.h"
#include "TrackCheckerBase.h"

using namespace Gaudi;
using namespace Gaudi::Units;
using namespace LHCb;

/** @class ExtrapolatorChecker ExtrapolatorChecker.h
 *
 *  This algorithm checks the performance of the TrackExtrapolators. This
 *  algorithm can be run in Boole or Brunel, since it relies only on the
 *  presence of MCHits. Algorithm description:
 *  @li create a State from a MCHit,
 *  @li extrapolate this State to the next MCHit on the track,
 *  @li compare extrapolated state with MCHit parameters,
 *  @li and make the corresponding resolution and pull plots.
 *
 *  The plots are split up in:
 *  @li first MCHit (extrapolated from true vertex)
 *  @li from VELO MCHit to VELO MCHit
 *  @li from TT MCHit to TT MCHit
 *  @li from IT MCHit to IT MCHit
 *  @li from OT MCHit to OT MCHit
 *  @li through RICH1: from VELO MCHit to TT MCHit
 *  @li through magnet: from TT MCHit to IT/OT MCHit
 *
 *  Note that the extrapolations are done starting from the @b entry @b point
 *  of the MCHit and extrapolating to the entry point of the next MCHit. A
 *  correction for the magnetic field is applied on the slopes (@e tx and @e ty)
 *  of the MCHit. This is required since the displacement vector of the MCHit
 *  is calculated from the difference between the exit and entry point, and
 *  therefore this displacement vector already accounts for half the effect of
 *  the magnetic field. This effect is subtracted in the correction method.
 *
 *  @author Jeroen van Tilburg
 *  @date   2006-07-06
 */

class ExtrapolatorChecker
    : public LHCb::Algorithm::Consumer<void( MCParticles const&, DetectorElement const& ),
                                       LHCb::DetDesc::usesBaseAndConditions<TrackCheckerBase, DetectorElement>> {

public:
  /// Standard constructor
  ExtrapolatorChecker( const std::string& name, ISvcLocator* pSvcLocator );

  /// Algorithm execution
  void operator()( const MCParticles&, DetectorElement const& ) const override;

private:
  /// Find the next MCHit belonging to the same MCParticle starting from a z-pos
  int findNextHit( const LHCb::MCParticle* mcPart, const double zRec, LHCb::MCHit const*& closestHit,
                   std::string& detectorName ) const;

  /// Helper function to find the next MCHit
  LHCb::MCHit const* findNextXxxHit( const LHCb::MCParticle* mcPart, const double zRec,
                                     LHCb::LinksByKey const* links ) const;

  /// Get the q/p for a given MCHit
  double qOverP( const LHCb::MCParticle* mcPart, const LHCb::MCHit* mcHit ) const;

  /// Correct slopes for magnetic field given an MCHit and a MCParticle
  void correctSlopes( const LHCb::MCParticle* mcPart, const LHCb::MCHit* mcHit, double& tx, double& ty ) const;

  /// String of the available detectors.
  Gaudi::Property<std::vector<std::string>> m_dets{this, "Detectors", {"IT", "OT", "TT", "Velo"}};
};

DECLARE_COMPONENT( ExtrapolatorChecker )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
ExtrapolatorChecker::ExtrapolatorChecker( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer{name,
               pSvcLocator,
               {KeyValue{"MCParticleLocation", MCParticleLocation::Default},
                KeyValue{"MCParticleInContainer", LHCb::MCParticleLocation::Default}}} {}

//=============================================================================
// Main execution
//=============================================================================
void ExtrapolatorChecker::operator()( const MCParticles& particles, DetectorElement const& geometry ) const {

  // Loop over the MCParticles
  for ( MCParticle const* mcPart : particles ) {
    // Decide whether the MCParticle will be checked
    if ( selected( mcPart ) ) {

      // Get the state at the vertex
      State stateVtx;
      idealStateCreator()
          ->createStateVertex( mcPart, stateVtx )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

      // Find first MCHit
      std::string  detName;
      MCHit const* mcHit = nullptr;
      int          detID = findNextHit( mcPart, stateVtx.z(), mcHit, detName );

      // Get the entry point of the MCHit
      XYZPoint entryP = mcHit->entry();

      // Extrapolate through RF foil
      extrapolator()
          ->propagate( stateVtx, entryP.z(), geometry )
          .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
      TrackVector    vec = stateVtx.stateVector();
      TrackSymMatrix cov = stateVtx.covariance();

      // Correct tx and ty from the MCHit for the magnetic field
      double tx, ty;
      correctSlopes( mcPart, mcHit, tx, ty );

      // Determine the resolutions
      double dx  = vec( 0 ) - entryP.x();
      double dy  = vec( 1 ) - entryP.y();
      double dtx = vec( 2 ) - tx;
      double dty = vec( 3 ) - ty;
      double dqp = vec( 4 ) - qOverP( mcPart, mcHit );

      // fill the histograms
      plot1D( dx, 1, "X resolution at 1st meas", -0.5, 0.5, 100 );
      plot1D( dy, 2, "Y resolution at 1st meas", -0.5, 0.5, 100 );
      plot1D( dtx, 3, "Tx resolution at 1st meas", -0.01, 0.01, 100 );
      plot1D( dty, 4, "Ty resolution at 1st meas", -0.01, 0.01, 100 );
      plot1D( stateVtx.p() - mcHit->p(), 5, "dp at 1st meas", -0.3, 0.3, 100 );
      if ( cov( 0, 0 ) != 0 && cov( 1, 1 ) != 0 && cov( 2, 2 ) != 0 && cov( 3, 3 ) != 0 ) {
        plot1D( dx / sqrt( cov( 0, 0 ) ), 11, "X pull at 1st meas", -5., 5., 100 );
        plot1D( dy / sqrt( cov( 1, 1 ) ), 12, "Y pull at 1st meas", -5., 5., 100 );
        plot1D( dtx / sqrt( cov( 2, 2 ) ), 13, "Tx pull at 1st meas", -5., 5., 100 );
        plot1D( dty / sqrt( cov( 3, 3 ) ), 14, "Ty pull at 1st meas", -5., 5., 100 );
      }
      if ( cov( 4, 4 ) != 0 ) plot1D( dqp / sqrt( cov( 4, 4 ) ), 15, "q/p pull at 1st meas", -5., 5., 100 );

      const Gaudi::TrackSymMatrix zeroCov;
      State                       state;
      state.setState( entryP.x(), entryP.y(), entryP.z(), tx, ty, qOverP( mcPart, mcHit ) );
      state.setCovariance( zeroCov );

      detID = findNextHit( mcPart, state.z(), mcHit, detName );

      while ( mcHit ) {

        entryP    = mcHit->entry();
        double dz = entryP.z() - state.z();

        // Extrapolate to next MCHit
        extrapolator()
            ->propagate( state, entryP.z(), geometry )
            .ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );
        TrackVector    vec = state.stateVector();
        TrackSymMatrix cov = state.covariance();

        // Correct tx and ty from the MCHit for the magnetic field
        correctSlopes( mcPart, mcHit, tx, ty );

        // Determine the resolutions
        double dx  = vec( 0 ) - entryP.x();
        double dy  = vec( 1 ) - entryP.y();
        double dtx = vec( 2 ) - tx;
        double dty = vec( 3 ) - ty;
        double dqp = vec( 4 ) - qOverP( mcPart, mcHit );

        // Define the ranges for the resolution plots
        double xr = 0.04;
        double tr = 0.001;
        double pr = 0.3;

        // Determine the histogram title
        int         ID    = 100 * detID;
        std::string title = " at " + detName + " hits";
        if ( dz > 4. * m && ( detName == "OT" || detName == "IT" ) ) {
          ID    = 1100;
          title = " after magnet extrapolation";
          xr    = 5.;
          tr    = 0.005;
          pr    = 5.;
        } else if ( dz > 1. * m && detName == "TT" ) {
          ID    = 1000;
          title = " after RICH1 extrapolation";
          xr    = 5.;
          tr    = 0.005;
          pr    = 10.;
        }

        // fill the histograms
        plot1D( dx, ID + 1, "X resolution" + title, -xr, xr, 100 );
        plot1D( dy, ID + 2, "Y resolution" + title, -xr, xr, 100 );
        plot1D( dtx, ID + 3, "Tx resolution" + title, -tr, tr, 100 );
        plot1D( dty, ID + 4, "Ty resolution" + title, -tr, tr, 100 );
        plot1D( state.p() - mcHit->p(), ID + 5, "dp" + title, -pr, pr, 100 );
        if ( cov( 0, 0 ) != 0 && cov( 1, 1 ) != 0 && cov( 2, 2 ) != 0 && cov( 3, 3 ) != 0 ) {
          plot1D( dx / sqrt( cov( 0, 0 ) ), ID + 11, "X pull" + title, -5., 5., 100 );
          plot1D( dy / sqrt( cov( 1, 1 ) ), ID + 12, "Y pull" + title, -5., 5., 100 );
          plot1D( dtx / sqrt( cov( 2, 2 ) ), ID + 13, "Tx pull" + title, -5., 5., 100 );
          plot1D( dty / sqrt( cov( 3, 3 ) ), ID + 14, "Ty pull" + title, -5., 5., 100 );
        }
        if ( cov( 4, 4 ) != 0 ) plot1D( dqp / sqrt( cov( 4, 4 ) ), ID + 15, "q/p pull" + title, -5., 5., 100 );

        state.setState( entryP.x(), entryP.y(), entryP.z(), tx, ty, qOverP( mcPart, mcHit ) );
        state.setCovariance( zeroCov );

        detID = findNextHit( mcPart, entryP.z(), mcHit, detName );
      }
    }
  } // End loop over MCParticles
}

//=============================================================================
// Find the next MCHit starting from a given z
// looping over the hits in all the tracking detectors
//=============================================================================
int ExtrapolatorChecker::findNextHit( const MCParticle* mcPart, const double zRec, MCHit const*& closestHit,
                                      std::string& detectorName ) const {
  detectorName = "Not found!";
  int detID    = 0;
  closestHit   = nullptr;
  ;
  double closestZ = 100000.;

  for ( auto itDets = m_dets.begin(); itDets != m_dets.end(); ++itDets ) {
    auto links = SmartDataPtr<LHCb::LinksByKey>{
        evtSvc(), LHCb::LinksByKey::linkerName( MCParticleLocation::Default + "2MC" + *itDets + "Hits" )};

    MCHit const* tmpClosestHit = findNextXxxHit( mcPart, zRec, links );
    if ( tmpClosestHit && tmpClosestHit->entry().z() > zRec && tmpClosestHit->entry().z() < closestZ ) {
      closestHit   = tmpClosestHit;
      closestZ     = tmpClosestHit->entry().z();
      detectorName = *itDets;
      detID        = 1 + itDets - m_dets.begin();
    }
  }
  return detID;
}

//=============================================================================
// Find the next MCHit of type Xxx starting from a given z
//=============================================================================
const MCHit* ExtrapolatorChecker::findNextXxxHit( const MCParticle* mcPart, const double zRec,
                                                  const LHCb::LinksByKey* links ) const {
  // Retrieve MCParticle to MCHit linker tables
  double       closestZ   = 100000.;
  const MCHit* closestHit = nullptr;
  for ( const auto& mcHit : LinkedFrom<MCHit>{links}.range( mcPart ) ) {
    double zOfHit = mcHit.entry().z();
    // get the closest hit
    if ( zOfHit > zRec + 0.1 && zOfHit < closestZ ) {
      closestHit = &mcHit;
      closestZ   = zOfHit;
    }
  }

  return closestHit;
}

//=============================================================================
// Determine q/p given an MCHit and a MCParticle
//=============================================================================
double ExtrapolatorChecker::qOverP( const MCParticle* mcPart, const MCHit* mcHit ) const {
  double charge   = ( mcPart->particleID().threeCharge() ) / 3.;
  double momentum = mcPart->p();
  if ( mcHit != NULL && mcHit->p() != 0. ) momentum = mcHit->p();
  return ( momentum > TrackParameters::lowTolerance ) ? charge / momentum : 0.;
}

//=============================================================================
// Correct slopes for magnetic field given an MCHit and a MCParticle
//=============================================================================
void ExtrapolatorChecker::correctSlopes( const MCParticle* mcPart, const MCHit* mcHit, double& tx, double& ty ) const {
  // TODO: I hope this method can be removed as soon as the displacement vector
  // in the MCHit is calculated in Gauss using the momentum direction of the
  // *entry point*. (JvT: 27/10/2006).

  // Get magnetic field vector
  Gaudi::XYZVector B;
  fieldSvc()->fieldVector( mcHit->midPoint(), B ).ignore( /* AUTOMATICALLY ADDED FOR gaudi/Gaudi!763 */ );

  // Calculate new displacement vector and tx,ty slopes
  Gaudi::XYZVector d    = mcHit->displacement();
  Gaudi::XYZVector dNew = d - ( 0.5 * d.R() * qOverP( mcPart, mcHit ) * d.Cross( B ) * eplus * c_light );
  tx                    = dNew.x() / dNew.z();
  ty                    = dNew.y() / dNew.z();
}
