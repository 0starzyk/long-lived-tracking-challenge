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
#include "DetDesc/DetectorElement.h"
#include "Event/Measurement.h"
#include "Event/State.h"
#include "Event/Track.h"
#include "Event/UTCluster.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Kernel/LHCbID.h"
#include "Kernel/SiChargeFun.h"
#include "Kernel/UTNames.h"
#include "Map.h"
#include "TrackKernel/TrackFunctors.h"
#include "TrackMonitorBase.h"
#include "UTDet/DeUTSector.h"

#include "LHCbAlgs/Consumer.h"

#include "AIDA/IHistogram1D.h"

#include "range/v3/view/transform.hpp"

#include <map>

using namespace LHCb;
using namespace Gaudi;

/**
 * Class for UT track monitoring
 *  @author M. Needham.
 *  @date   16-1-2009
 */
class UTTrackMonitor
    : public LHCb::Algorithm::Consumer<void( LHCb::Track::Range const&, LHCb::UTClusters const&,
                                             DetectorElement const& ),
                                       LHCb::DetDesc::usesBaseAndConditions<TrackMonitorBase, DetectorElement>> {

public:
  UTTrackMonitor( const std::string& name, ISvcLocator* pSvcLocator );
  void operator()( LHCb::Track::Range const&, LHCb::UTClusters const&, DetectorElement const& ) const override;

private:
  void fillHistograms( LHCb::Track const& track, LHCb::UTClusters const& fullclusters,
                       std::vector<LHCb::LHCbID> const& utIDs, IGeometryInfo const& geometry ) const;

  double m_xMax;
  double m_yMax;

  Gaudi::Property<double> m_refZ{this,
                                 "ReferenceZ",
                                 2500.0,
                                 [=]( auto& ) {
                                   // guess that the size of the histograms at ref is ~ 0.35 by 0.3
                                   this->m_xMax = 0.35 * this->m_refZ / Gaudi::Units::cm;
                                   this->m_yMax = 0.3 * this->m_refZ / Gaudi::Units::cm;
                                 },
                                 Gaudi::Details::Property::ImmediatelyInvokeHandler{true},
                                 "midpoint of UT"};

  Gaudi::Property<unsigned int> m_minNumUTHits{this, "minNumUTHits", 2u};
};

namespace {
  unsigned int histoBin( const LHCb::Detector::UT::ChannelID& chan ) {

    // convert layer and station to a flat number
    unsigned int layer = ( chan.station() == 1u ? chan.layer() : ( chan.layer() + 2 ) );
    return layer * 400 + ( chan.detRegion() - 1 ) * 120 + chan.sector();
  }
} // namespace

DECLARE_COMPONENT( UTTrackMonitor )

UTTrackMonitor::UTTrackMonitor( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator,
                {{"TracksInContainer", LHCb::TrackLocation::Default},
                 {"InputData", LHCb::UTClusterLocation::UTClusters},
                 {"StandardGeometryTop", LHCb::standard_geometry_top}} ) {}

void UTTrackMonitor::operator()( LHCb::Track::Range const& tracks, LHCb::UTClusters const& clusters,
                                 const DetectorElement& lhcb ) const {
  auto& geometry = *lhcb.geometry();

  // tmp container for ids
  std::vector<unsigned int> usedIDs;
  usedIDs.reserve( clusters.size() );

  // histograms per track
  for ( const LHCb::Track* track : tracks )
    if ( track->hasUT() ) {

      // find the IT hits on the track
      const auto&               ids = track->lhcbIDs();
      std::vector<LHCb::LHCbID> ttHits;
      ttHits.reserve( ids.size() );
      std::copy_if( ids.begin(), ids.end(), std::back_inserter( ttHits ),
                    []( const LHCb::LHCbID& id ) { return id.isUT(); } );

      if ( ids.size() < m_minNumUTHits ) continue;

      // insert into tmp container
      std::transform( ids.begin(), ids.end(), std::back_inserter( usedIDs ),
                      []( const LHCb::LHCbID& id ) { return id.utID(); } );

      fillHistograms( *track, clusters, ids, geometry );
    }

  // see how many hits in the IT were actually used
  if ( !clusters.empty() ) {
    std::sort( usedIDs.begin(), usedIDs.end() );
    std::unique( usedIDs.begin(), usedIDs.end() );
    plot( usedIDs.size() / (double)clusters.size(), "fraction of used UT hits", 0., 1., 50 );
  }
}

namespace {
  size_t ttUniqueSectorID( const LHCb::Detector::UT::ChannelID& id ) {
    // std::cout << "tt: " << id.station() << " " << id.layer() << " " << id.detRegion() << " " << id.sector() <<
    // std::endl ;
    return ( id.detRegion() - 1 ) * 24 + id.sector() - 1;
  }
} // namespace

void UTTrackMonitor::fillHistograms( LHCb::Track const& track, UTClusters const& fullclusters,
                                     std::vector<LHCb::LHCbID> const& utIDs, IGeometryInfo const& geometry ) const {
  // track parameters at some reference z
  LHCb::StateVector aState;
  extrapolator()->propagate( track, m_refZ, aState, geometry ).ignore();
  if ( fullDetail() ) {
    plot2D( aState.x() / Gaudi::Units::cm, aState.y() / Gaudi::Units::cm, "xy", "x vs y", -m_xMax, m_xMax, -m_yMax,
            m_yMax, 50, 50 );
  }
  plot( aState.tx(), "tx", "tx", -0.2, 0.2, 200 );
  plot( aState.ty(), "ty", "ty", -0.2, 0.2, 200 );
  plot( track.p() / Gaudi::Units::GeV, "momentum", "momentum", -5., 205., 21 );
  plot( aState.x() / Gaudi::Units::cm, "xdist", "x dist", -m_xMax, m_xMax );
  plot( aState.y() / Gaudi::Units::cm, "ydist", "y dist", -m_yMax, m_yMax );

  // Loop over the nodes to get the hits variables
  std::vector<const UTCluster*> clusters;
  clusters.reserve( 24 );
  std::vector<const LHCb::FitNode*> nodesByUTLayer[4];

  unsigned int nHigh = 0u;
  for ( const LHCb::FitNode* node : nodes( track ) ) {
    // Only loop on hits with measurement that is IT type
    if ( !node->hasMeasurement() ) continue;

    // monitor overlapping tracks
    const LHCb::Measurement&     measurement = node->measurement();
    const LHCb::Measurement::UT* hit         = node->measurement().getIf<LHCb::Measurement::UT>();
    if ( !hit ) continue;
    LHCb::LHCbID lhcbID = measurement.lhcbID();
    assert( lhcbID.isUT() );

    // unbiased residuals and biased residuals
    Detector::UT::ChannelID chan        = lhcbID.utID();
    unsigned int            uniquelayer = ( chan.station() - 1 ) * 2 + chan.layer() - 1;
    nodesByUTLayer[uniquelayer].push_back( node );

    const std::string layerName = UTNames().UniqueLayerToString( chan );
    plot( node->unbiasedResidual(), "unbiasedResidual", "unbiasedResidual", -2., 2., 200 );
    plot( node->residual(), "biasedResidual", "biasedResidual", -2., 2., 200 );
    // rms unbiased residual. that's the one you want to look at.
    double residual = node->residual() * std::sqrt( node->errMeasure2() / node->errResidual2() );
    plot( residual, layerName + "/Residual", "Residual (rms unbiased)", -0.5, 0.5, 100 );

    // 2D plots in full detail mode
    if ( fullDetail() ) {
      const unsigned int bin = histoBin( chan );
      plot2D( bin, node->unbiasedResidual(), "unbiasedResSector" + layerName, "unbiasedResSector" + layerName, 99.5,
              500.5, -2., 2., 401, 200 );
      plot2D( bin, node->residual(), "biasedResSector" + layerName, "/biasedResSector" + layerName, 99.5, 500.5, -2.,
              2., 401, 200 );

      double noise = hit->sec().noise( chan );
      if ( noise > 0 ) {
        const UTCluster* cluster       = fullclusters.object( chan );
        const double     signalToNoise = cluster->totalCharge() / hit->sec().noise( chan );
        plot2D( bin, signalToNoise, "SNSector" + layerName, "SNSector" + layerName, 99.5, 500.5, -0.25, 100.25, 401,
                201 );
        plot2D( bin, cluster->totalCharge(), "CSector" + layerName, "CSector" + layerName, 99.5, 500.5, -0.5, 200.5,
                401, 201 );
      }
    }

    const UTCluster* cluster = fullclusters.object( chan );

    // get the measurement and plot UT related quantities
    plot( cluster->totalCharge(), "charge", "clusters charge", 0., 200., 100 );
    plot( cluster->size(), "size", "cluster size", -0.5, 10.5, 11 );

    // for this one we actually need the component perpendicular to B field.
    Gaudi::XYZVector dir = node->state().slopes();
    profile1D( dir.x(), cluster->size(), "cluster size vs tx", "cluster size vs tx", -0.3, 0.3 );
    profile1D( std::sqrt( dir.x() * dir.x() + dir.y() * dir.y() ) / dir.z(), cluster->totalCharge(),
               "cluster charge vs slope", "cluster charge vs local slope", 0, 0.4 );

    if ( cluster->highThreshold() ) ++nHigh;
    clusters.push_back( cluster );

  } // nodes

  // plot high fraction and also shorth
  plot( nHigh / double( utIDs.size() ), "highFrac", "high fraction", 0., 1., 100 );
  const double shorth = SiChargeFun::shorth( clusters.begin(), clusters.end() );
  plot( shorth, "shorth", "shorth charge", 0., 100., 200 );
  const double gm = SiChargeFun::generalizedMean( clusters.begin(), clusters.end() );
  plot( gm, "generalized mean", "generalized mean charge", 0., 100, 200 );

  // make overlap histograms
  for ( size_t ilay = 0; ilay < 4; ++ilay ) {
    const char* layname[4] = {"UTaX", "UTaU", "UTbV", "UTbX"};
    std::string prefix     = std::string( "/" ) + layname[ilay] + "/";
    plot( nodesByUTLayer[ilay].size(), prefix + "Number of hits", "Number of hits", -0.5, 5.5 );

    if ( nodesByUTLayer[ilay].size() == 2 ) {
      const LHCb::FitNode* firstnode  = nodesByUTLayer[ilay].front();
      const LHCb::FitNode* secondnode = nodesByUTLayer[ilay].back();

      if ( firstnode->measurement().isSameDetectorElement( secondnode->measurement() ) ) {

        if ( firstnode->measurement().toGlobal( Gaudi::XYZPoint() ).x() >
             secondnode->measurement().toGlobal( Gaudi::XYZPoint() ).x() )
          std::swap( firstnode, secondnode );

        auto   firstTrajectory  = firstnode->measurement().getIf<LHCb::Measurement::UT>()->trajectory;
        auto   secondTrajectory = secondnode->measurement().getIf<LHCb::Measurement::UT>()->trajectory;
        int    firstsign        = firstTrajectory.direction( 0 ).y() > 0 ? 1 : -1;
        int    secondsign       = secondTrajectory.direction( 0 ).y() > 0 ? 1 : -1;
        double firstresidual    = firstsign * firstnode->residual();
        double secondresidual   = secondsign * secondnode->residual();

        double diff = firstresidual - secondresidual;

        // let's still make the correction for the track angle, such that we really look in the wafer plane.
        Gaudi::XYZVector localdir = firstnode->measurement().toLocal( firstnode->state().slopes() );
        double           localTx  = localdir.x() / localdir.z();
        diff *= std::sqrt( 1 + localTx * localTx );

        const std::string layerName = UTNames().UniqueLayerToString( firstnode->measurement().lhcbID().utID() );
        plot( diff, prefix + "Overlap residual", std::string( "Overlap residuals in " ) + layerName, -1.0, 1.0, 100 );
        plot( firstresidual * std::sqrt( firstnode->errMeasure2() / firstnode->errResidual2() ),
              prefix + "Residuals in overlaps (left)", std::string( "Residuals in overlaps (left) in " ) + layerName,
              -0.5, 0.5, 100 );
        plot( secondresidual * std::sqrt( secondnode->errMeasure2() / secondnode->errResidual2() ),
              prefix + "Residuals in overlaps (right)", std::string( "Residuals in overlaps (right) in " ) + layerName,
              -0.5, 0.5, 100 );

        // this needs to be fixed: can we somehow can a consecutive ladder ID?
        size_t sectorID = ttUniqueSectorID( firstnode->measurement().lhcbID().utID() );
        profile1D( double( sectorID ), diff, prefix + "Overlap residual versus sector ID",
                   "Average overlap residual versus sector ID", -0.5, 47.5, 48 );

        if ( fullDetail() ) {
          std::string firstname  = firstnode->measurement().name().substr( 55, 10 );
          std::string secondname = secondnode->measurement().name().substr( 55, 10 );
          plot( diff, prefix + "Overlap residual for " + firstname + " and " + secondname,
                "Overlap residual for " + firstname + " and " + secondname, -1.0, 1.0, 50 );
        }
      }
    }
  }
}
