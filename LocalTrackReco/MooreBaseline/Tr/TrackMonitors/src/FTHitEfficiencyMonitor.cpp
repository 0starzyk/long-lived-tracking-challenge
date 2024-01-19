/*****************************************************************************\
* (c) Copyright 2022 CERN for the benefit of the LHCb Collaboration      *
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
#include "Detector/FT/FTChannelID.h"
#include "Detector/FT/FTConstants.h"
#include "Event/PrHits.h"
#include "FTDet/DeFTDetector.h"

#include "Event/PrHits.h"
#include "Event/Track.h"

#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "LHCbAlgs/Consumer.h"
#include <Gaudi/Accumulators/Histogram.h>

#include "GaudiAlg/GaudiHistoAlg.h"

#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackInterfaces/ITrackInterpolator.h"

/**
 * Algorithm based on the track-based efficiencies, in particular the
 * OTHitEfficiencyMonitor
 * conceptualised by Wouter Hulsbergen, modified by Francesco Dettori,
 * and others in the next generations.
 *
 * Written a SciFi-version by Laurent Dufour & Zehua Xu.
 *
 * To use this, please take proper care of the hit unbiasing beforehand.
 * This is most adequately done by unbiasing the pattern recognition.
 * In addition, due to the complex 3D histograms, it is not intended
 * to be ran online as-is.
 *
 * A lightweight version, in which only the 1D histograms are
 * created, can be used.
 *
 * When setting the compiler define FTHITEFFMON_ADD_DEBUG_COUNTERS
 * additional counters will be added to this algorithm that
 * could help debugging why certain hits are not considered.
 **/
class FTHitEfficiencyMonitor
    : public LHCb::Algorithm::Consumer<void( const LHCb::Track::Range&, const LHCb::Pr::FT::Hits&,
                                             const DetectorElement&, const DeFT& ),
                                       LHCb::DetDesc::usesBaseAndConditions<GaudiHistoAlg, DetectorElement, DeFT>> {
private:
  Gaudi::Property<unsigned int> m_layerUnderStudy{this, "LayerUnderStudy"};

  Gaudi::Property<float> m_tolerance{this, "ResidualTolerance", 2.5 * Gaudi::Units::mm}; // tolerance in mm
  Gaudi::Property<float> m_maxTrackError{this, "MaxTrackCov", 0.2 * Gaudi::Units::mm};

  Gaudi::Property<float> m_track_min_momentum{this, "TrackMinP", 5000 * Gaudi::Units::MeV};
  Gaudi::Property<float> m_track_min_pt{this, "TrackMinPT", 400 * Gaudi::Units::MeV};

  ToolHandle<ITrackInterpolator> m_interpolator{this, "Interpolator", "TrackInterpolator"};
  ToolHandle<ITrackExtrapolator> m_linearextrapolator{this, "Extrapolator", "TrackLinearExtrapolator"};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_efficiencyPerModule{
      this,
      "efficiencyVsModuleIdx",
      "efficiency per module",
      {FTConstants::nModulesTotal, -0.5, FTConstants::nModulesTotal - 0.5}};

  mutable Gaudi::Accumulators::ProfileHistogram<2> m_efficiencyPerXY{
      this,
      "efficiencyVsXY",
      "efficiency per XY, layer under study",
      {660, -3200. * Gaudi::Units::mm, 3200. * Gaudi::Units::mm},
      {124, -2500 * Gaudi::Units::mm, 2500 * Gaudi::Units::mm}};

  mutable Gaudi::Accumulators::Histogram<3> m_residualPerXY{this,
                                                            "residualVsXY",
                                                            "residual per XY, layer under study",
                                                            {100, -3200. * Gaudi::Units::mm, 3200. * Gaudi::Units::mm},
                                                            {124, -2500 * Gaudi::Units::mm, 2500 * Gaudi::Units::mm},
                                                            {100, -8.0 * Gaudi::Units::mm, 8.0 * Gaudi::Units::mm}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_efficiencyPerX{
      this, "efficiencyVsX", "efficiency per X", {400, -3200. * Gaudi::Units::mm, 3200. * Gaudi::Units::mm}};

#ifdef FTHITEFFMON_ADD_DEBUG_COUNTERS
  mutable Gaudi::Accumulators::Counter<> m_wrongStationIdCounter{this, "Wrong station ID of hit"};
  mutable Gaudi::Accumulators::Counter<> m_wrongQuarterIdCounter{this, "Wrong quarter ID of hit"};
  mutable Gaudi::Accumulators::Counter<> m_wrongLayerIdCounter{this, "Wrong layer ID of hit"};
  mutable Gaudi::Accumulators::Counter<> m_wrongMatIdCounter{this, "Wrong mat ID of hit"};
  mutable Gaudi::Accumulators::Counter<> m_correctStationIdCounter{this, "Correct station ID of hit"};
  mutable Gaudi::Accumulators::Counter<> m_correctQuarterIdCounter{this, "Correct quarter ID of hit"};
  mutable Gaudi::Accumulators::Counter<> m_correctLayerIdCounter{this, "Correct layer ID of hit"};
  mutable Gaudi::Accumulators::Counter<> m_correctMatIdCounter{this, "Correct mat ID of hit"};
#endif

public:
  FTHitEfficiencyMonitor( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"TrackLocation", ""}, KeyValue{"PrFTHitsLocation", ""},
                   KeyValue{"LHCbGeoLocation", LHCb::standard_geometry_top},
                   KeyValue{"FTLocation", DeFTDetectorLocation::Default}} ) {}

  void operator()( const LHCb::Track::Range& tracks, const LHCb::Pr::FT::Hits& hits, const DetectorElement& lhcb,
                   const DeFT& ftDet ) const override {
    // in case you're thinking: wow, this is  absurd
    // you're correct.
    // we don't have a function FTDet::getLayer( ... );
    Gaudi::Plane3D my_layer_plane;
    float          z_position = -1;

#ifdef USE_DD4HEP
    auto func = [&]( const DeFTLayer& layer ) {
      const auto id = layer.layerIdx();
      if ( id == m_layerUnderStudy.value() ) {
        z_position     = layer.globalZ();
        my_layer_plane = layer.plane();
      }
    };

    ftDet.applyToAllLayers( func );
#else
    for ( auto station : ftDet.stations() ) {
      for ( auto layer : station->layers() ) {
        const auto id = 4 * ( station->stationID() - 1 ) + layer->layerID();
        if ( id == m_layerUnderStudy.value() ) {
          z_position     = layer->globalZ();
          my_layer_plane = layer->plane();
        }
      }
    }
#endif

    if ( z_position < 0 ) {
      error() << "Could not find this SciFi layer." << endmsg;
      return;
    }

    for ( const auto& track : tracks ) {
      if ( m_track_min_momentum.value() > 0 && track->p() < m_track_min_momentum.value() ) continue;

      if ( m_track_min_pt.value() > 0 && track->pt() < m_track_min_pt.value() ) continue;

      LHCb::State state;
      if ( !m_interpolator->interpolate( *track, z_position, state, lhcb ).isSuccess() ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Could not interpolate track." << endmsg;
        continue;
      }

      // only consider if state has reasonably small error. we should
      // still fine-tune this a little bit, make a propor projecion
      // etc. the best is to cut only when you have unbiased the
      // state. cutting too tight could introduce a bias in
      // the efficiency measurement?
      if ( std::sqrt( state.covariance()( 0, 0 ) ) > m_maxTrackError.value() ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Too large cov element" << endmsg;
        continue;
      }

      // to get a slightly more precise answer, intersect with the plane
      if ( !m_linearextrapolator->propagate( state, my_layer_plane, lhcb ).isSuccess() ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Could not extrapolate track" << endmsg;
        continue;
      }

      const auto& ft_mat = ftDet.findMat( state.position() );

      if ( !ft_mat ) {
        if ( msgLevel( MSG::VERBOSE ) ) verbose() << "Could not find a corresponding FTMat" << endmsg;
        continue;
      }

      // this is a float to help with the Profile-histogram call...
      // could well be a int/bool
      float foundHit = 0;

      /**
       * This is very suboptimal:
       * what one wants to do, is to compute the channelIDs which
       * fall under the tolerance. Then check if these are in the
       * hits container. Now the x position is approximately calcualted
       * and then a tolerance is applied... this is probably
       * quite inefficient computing wise.
       **/
      for ( unsigned int i = 0; i < hits.size(); i++ ) {
        if ( hits.id( i ) == LHCb::Detector::FTChannelID{} ) continue;
        auto ftChannelID = hits.id( i );

#ifdef FTHITEFFMON_ADD_DEBUG_COUNTERS
        if ( unsigned( ftChannelID.mat() ) != unsigned( ft_mat->matID() ) ) {
          ++m_wrongMatIdCounter;
        } else {
          ++m_correctMatIdCounter;
        }

        if ( unsigned( ftChannelID.layer() ) != unsigned( ft_mat->layerID() ) ) {
          ++m_wrongLayerIdCounter;
        } else {
          ++m_correctLayerIdCounter;
        }

        if ( unsigned( ftChannelID.station() ) != unsigned( ft_mat->stationID() ) ) {
          ++m_wrongStationIdCounter;
        } else {
          ++m_correctStationIdCounter;
        }

        if ( unsigned( ftChannelID.quarter() ) != unsigned( ft_mat->quarterID() ) ) {
          ++m_wrongQuarterIdCounter;
        } else {
          ++m_correctQuarterIdCounter;
        }
#endif

        if ( ftChannelID.globalMatID() != ft_mat->elementID().globalMatID() ) continue;

        const float residual = hits.distance( i, state.x(), state.y() ); // x_at_state - state.x();
        if ( fabs( residual ) < m_tolerance ) foundHit = 1;

        ++m_residualPerXY[{state.x(), state.y(), residual}];
      }

      m_efficiencyPerModule[ft_mat->elementID().globalModuleIdx()] += foundHit;
      m_efficiencyPerXY[{state.x(), state.y()}] += foundHit; // Fill( state.x(), state.y(), foundHit );
      m_efficiencyPerX[state.x()] += foundHit;               // Fill( state.x(), state.y(), foundHit );
    }
  }
};
DECLARE_COMPONENT( FTHitEfficiencyMonitor )
