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
#include "Event/Track.h"
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "IPVSeeding.h" // Interface
#include <cmath>

// auxiliary class for searching of clusters of tracks
struct vtxCluster final {

  double z          = 0; // z of the cluster
  double sigsq      = 0; // sigma**2 of the cluster
  double sigsqmin   = 0; // minimum sigma**2 of the tracks forming cluster
  int    ntracks    = 0; // number of tracks in the cluster
  int    not_merged = 0; // flag for iterative merging

  vtxCluster() = default;
};

/** @class PVSeedTool PVSeedTool.h tmp/PVSeedTool.h
 *
 *
 *  @author Mariusz Witek
 *  @date   2005-11-19
 */
class PVSeedTool : public extends<GaudiTool, IPVSeeding> {
public:
  /// Standard constructor
  using extends::extends;

  std::vector<Gaudi::XYZPoint> getSeeds( LHCb::span<const LHCb::Track* const> inputTracks,
                                         const Gaudi::XYZPoint&               beamspot,
                                         IGeometryInfo const&                 geometry ) const override;

private:
  std::vector<double> findClusters( std::vector<vtxCluster>& vclus ) const;
  void                errorForPVSeedFinding( double tx, double ty, double& sigzaq ) const;

  double zCloseBeam( const LHCb::Track* track, const Gaudi::XYZPoint& beamspot ) const;

  // steering parameters for merging procedure
  Gaudi::Property<double> m_maxChi2Merge{this, "maxChi2Merge", 25.};
  Gaudi::Property<double> m_factorToIncreaseErrors{this, "factorToIncreaseErrors", 15.};

  // steering parameters for final cluster selection
  Gaudi::Property<int>    m_minClusterMult{this, "minClusterMult", 3};
  Gaudi::Property<double> m_dzCloseTracksInCluster{this, "dzCloseTracksInCluster", 5. * Gaudi::Units::mm};
  Gaudi::Property<int>    m_minCloseTracksInCluster{this, "minCloseTracksInCluster", 3};
  Gaudi::Property<int>    m_highMult{this, "highMult", 10};
  Gaudi::Property<double> m_ratioSig2HighMult{this, "ratioSig2HighMult", 1.0};
  Gaudi::Property<double> m_ratioSig2LowMult{this, "ratioSig2LowMult", 0.9};

  double                  m_scatCons = 0; // calculated from m_x0MS
  Gaudi::Property<double> m_x0MS{this, "x0MS", 0.01,
                                 [=]( auto& ) {
                                   double X0        = this->m_x0MS;
                                   this->m_scatCons = ( 13.6 * sqrt( X0 ) * ( 1. + 0.038 * log( X0 ) ) );
                                 },
                                 Gaudi::Details::Property::ImmediatelyInvokeHandler{true}}; // X0 (tunable) of MS to add
                                                                                            // for extrapolation of
                                                                                            // track parameters to PV
};

namespace {
  bool vtxcomp( vtxCluster* first, vtxCluster* second ) { return first->z < second->z; }
  bool multcomp( vtxCluster* first, vtxCluster* second ) { return first->ntracks > second->ntracks; }

  // auxiliary class for merging procedure of tracks/clusters
  struct pair_to_merge final {

    vtxCluster* first    = nullptr; // pointer to first cluster to be merged
    vtxCluster* second   = nullptr; // pointer to second cluster to be merged
    double      chi2dist = 10.e+10; // a chi2dist = zdistance**2/(sigma1**2+sigma2**2)

    pair_to_merge( vtxCluster* f, vtxCluster* s, double chi2 ) : first( f ), second( s ), chi2dist( chi2 ) {}
  };

  bool paircomp( const pair_to_merge& first, const pair_to_merge& second ) { return first.chi2dist < second.chi2dist; }

  // void print_clusters(const std::vector<vtxCluster*>& pvclus) {
  //  std::cout << pvclus.size() << " clusters at this iteration" << std::endl;
  //  for(const auto& vc : pvclus) {
  //     std::cout << vc->ntracks << " " << vc->z << " "
  //            <<  vc->sigsq << " " <<  vc->sigsqmin << std::endl;
  //  }
  //}

  constexpr static const int s_p2mstatic = 5000;
} // namespace

//-----------------------------------------------------------------------------
// Implementation file for class : PVSeedTool
//
// 2005-11-19 : Mariusz Witek
//-----------------------------------------------------------------------------

DECLARE_COMPONENT( PVSeedTool )

//=============================================================================
// getSeeds
//=============================================================================
std::vector<Gaudi::XYZPoint> PVSeedTool::getSeeds( LHCb::span<const LHCb::Track* const> inputTracks,
                                                   const Gaudi::XYZPoint& beamspot, IGeometryInfo const& ) const {

  std::vector<Gaudi::XYZPoint> seeds;
  if ( inputTracks.size() < 3 ) return seeds;

  std::vector<vtxCluster> vclusters;

  for ( const auto& trk : inputTracks ) {

    double sigsq;
    double zclu;

    zclu = zCloseBeam( trk, beamspot );
    errorForPVSeedFinding( trk->firstState().tx(), trk->firstState().ty(), sigsq );

    if ( fabs( zclu ) > 2000. ) continue;
    vtxCluster clu;
    clu.z        = zclu;
    clu.sigsq    = sigsq;
    clu.sigsqmin = clu.sigsq;
    clu.ntracks  = 1;
    vclusters.push_back( clu );
  }

  auto zseeds = findClusters( vclusters );

  seeds.reserve( zseeds.size() );
  std::transform( zseeds.begin(), zseeds.end(), std::back_inserter( seeds ), [&]( double z ) {
    return Gaudi::XYZPoint{beamspot.X(), beamspot.Y(), z};
  } );

  return seeds;
}

std::vector<double> PVSeedTool::findClusters( std::vector<vtxCluster>& vclus ) const {

  std::vector<double> zclusters;
  if ( vclus.empty() ) return zclusters;

  std::vector<vtxCluster*> pvclus;
  pvclus.reserve( vclus.size() );

  for ( auto& itvtx : vclus ) {
    itvtx.sigsq *= m_factorToIncreaseErrors * m_factorToIncreaseErrors; // blow up errors
    itvtx.sigsqmin = itvtx.sigsq;
    pvclus.push_back( &itvtx );
  }

  std::sort( pvclus.begin(), pvclus.end(), vtxcomp );
  //  print_clusters(pvclus);

  bool more_merging = true;
  while ( more_merging ) {
    // find pair of clusters for merging

    // refresh flag for this iteration
    for ( auto ivc = pvclus.begin(); ivc != pvclus.end(); ivc++ ) { ( *ivc )->not_merged = 1; }

    // for a given cluster look only up to a few consequtive ones to merge
    // "a few" might become a property?
    auto n_consequtive = std::min( 5, static_cast<int>( pvclus.size() ) );
    auto ivc2up        = std::next( pvclus.begin(), n_consequtive );

    std::vector<pair_to_merge> vecp2m;
    vecp2m.reserve( std::min( static_cast<int>( pvclus.size() ) * n_consequtive, s_p2mstatic ) );
    for ( auto ivc1 = pvclus.begin(); ivc1 != pvclus.end() - 1; ivc1++ ) {
      if ( ivc2up != pvclus.end() ) ++ivc2up;
      for ( auto ivc2 = std::next( ivc1 ); ivc2 != ivc2up; ivc2++ ) {
        double zdist    = ( *ivc1 )->z - ( *ivc2 )->z;
        double chi2dist = zdist * zdist / ( ( *ivc1 )->sigsq + ( *ivc2 )->sigsq );
        if ( chi2dist < m_maxChi2Merge && vecp2m.size() < s_p2mstatic ) {
          vecp2m.emplace_back( *ivc1, *ivc2, chi2dist );
        }
      }
    }

    // done if no more pairs to merge
    if ( vecp2m.empty() ) {
      more_merging = false;
    } else {
      // sort if number of pairs reasonable. Sorting increases efficency.
      if ( vecp2m.size() < 100 ) std::sort( vecp2m.begin(), vecp2m.end(), paircomp );

      // merge pairs
      for ( auto itp2m = vecp2m.begin(); itp2m != vecp2m.end(); itp2m++ ) {
        vtxCluster* pvtx1 = itp2m->first;
        vtxCluster* pvtx2 = itp2m->second;
        if ( pvtx1->not_merged == 1 && pvtx2->not_merged == 1 ) {

          double z1       = pvtx1->z;
          double z2       = pvtx2->z;
          double s1       = pvtx1->sigsq;
          double s2       = pvtx2->sigsq;
          double s1min    = pvtx1->sigsqmin;
          double s2min    = pvtx2->sigsqmin;
          double sigsqmin = s1min;
          if ( s2min < s1min ) sigsqmin = s2min;

          double w_inv  = ( s1 * s2 / ( s1 + s2 ) );
          double zmerge = w_inv * ( z1 / s1 + z2 / s2 );

          pvtx1->z        = zmerge;
          pvtx1->sigsq    = w_inv;
          pvtx1->sigsqmin = sigsqmin;
          pvtx1->ntracks += pvtx2->ntracks;
          pvtx2->ntracks    = 0; // mark second cluster as used
          pvtx1->not_merged = 0;
          pvtx2->not_merged = 0;
        }
      }

      // remove clusters which where used
      pvclus.erase(
          std::remove_if( pvclus.begin(), pvclus.end(), []( const vtxCluster* cl ) { return cl->ntracks < 1; } ),
          pvclus.end() );
    }
  }

  // End of clustering.

  // Sort according to multiplicity

  std::sort( pvclus.begin(), pvclus.end(), multcomp );

  // Select good clusters.

  for ( auto ivc = pvclus.begin(); ivc != pvclus.end(); ivc++ ) {
    int n_tracks_close = 0;
    for ( auto itvtx = vclus.begin(); itvtx != vclus.end(); itvtx++ ) {
      if ( fabs( itvtx->z - ( *ivc )->z ) < m_dzCloseTracksInCluster ) n_tracks_close++;
    }

    double dist_to_closest = 1000000.;
    if ( pvclus.size() > 1 ) {
      for ( auto ivc1 = pvclus.begin(); ivc1 != pvclus.end(); ivc1++ ) {
        if ( ivc != ivc1 && ( fabs( ( *ivc1 )->z - ( *ivc )->z ) < dist_to_closest ) ) {
          dist_to_closest = fabs( ( *ivc1 )->z - ( *ivc )->z );
        }
      }
    }

    // ratio to remove clusters made of one low error track and many large error ones
    double rat     = ( *ivc )->sigsq / ( *ivc )->sigsqmin;
    bool   igood   = false;
    int    ntracks = ( *ivc )->ntracks;
    if ( ntracks >= m_minClusterMult ) {
      if ( dist_to_closest > 10. && rat < 0.95 ) igood = true;
      if ( ntracks >= m_highMult && rat < m_ratioSig2HighMult ) igood = true;
      if ( ntracks < m_highMult && rat < m_ratioSig2LowMult ) igood = true;
    }
    // veto
    if ( n_tracks_close < m_minCloseTracksInCluster ) igood = false;
    if ( igood ) zclusters.push_back( ( *ivc )->z );
  }

  //  print_clusters(pvclus);
  return zclusters;
}

void PVSeedTool::errorForPVSeedFinding( double tx, double ty, double& sigz2 ) const {

  // the seeding results depend weakly on this eror parametrization

  double invPMean2 = 1. / ( 3000. * Gaudi::Units::MeV * 3000. * Gaudi::Units::MeV );

  double tanTheta2    = tx * tx + ty * ty;
  double invSinTheta2 = 1. / tanTheta2 + 1.;

  // assume that first hit in VD at 8 mm
  double distr        = 8. * Gaudi::Units::mm;
  double dist2        = distr * distr * invSinTheta2;
  double sigma_ms2    = m_scatCons * m_scatCons * dist2 * invPMean2;
  double fslope2      = 0.0005 * 0.0005;
  double sigma_slope2 = fslope2 * dist2;

  sigz2 = ( sigma_ms2 + sigma_slope2 ) * invSinTheta2;
}

double PVSeedTool::zCloseBeam( const LHCb::Track* track, const Gaudi::XYZPoint& beamspot ) const {

  Gaudi::XYZPoint  tpoint = track->position();
  Gaudi::XYZVector tdir   = track->slopes();

  double wx = ( 1. + tdir.x() * tdir.x() ) / track->firstState().errX2();
  double wy = ( 1. + tdir.y() * tdir.y() ) / track->firstState().errY2();

  double x0      = tpoint.x() - tpoint.z() * tdir.x() - beamspot.X();
  double y0      = tpoint.y() - tpoint.z() * tdir.y() - beamspot.Y();
  double den     = wx * tdir.x() * tdir.x() + wy * tdir.y() * tdir.y();
  double zAtBeam = -( wx * x0 * tdir.x() + wy * y0 * tdir.y() ) / den;

  double xb       = tpoint.x() + tdir.x() * ( zAtBeam - tpoint.z() ) - beamspot.X();
  double yb       = tpoint.y() + tdir.y() * ( zAtBeam - tpoint.z() ) - beamspot.Y();
  double r2AtBeam = xb * xb + yb * yb;

  return r2AtBeam < 0.5 * 0.5 ? zAtBeam : 10e8;
}

//=============================================================================