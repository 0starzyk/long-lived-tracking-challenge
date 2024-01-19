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
#include "Event/LinksByKey.h"
#include "Event/MCParticle.h"
#include "Event/MCVertex.h"
#include "Event/RecVertex.h"
#include "Event/Track.h"
#include "Event/TrackVertexUtils.h"
#include "Gaudi/Accumulators/Histogram.h"
#include "GaudiAlg/GaudiHistoAlg.h"
#include "LHCbAlgs/Consumer.h"
#include "TrackCheckerBase.h"
#include "TrackKernel/CubicStateInterpolationTraj.h"
#include "TrackKernel/TrackFunctors.h"
#include <map>

class TriggerObjectsCompatibilityProfileChecker final
    : public LHCb::Algorithm::Consumer<void( const LHCb::Tracks&, const LHCb::Tracks&, const LHCb::MCParticles&,
                                             const LHCb::LinksByKey&, const LHCb::LinksByKey&, const LHCb::RecVertices&,
                                             const LHCb::RecVertices& ),
                                       Gaudi::Functional::Traits::BaseClass_t<TrackCheckerBase>> {
public:
  TriggerObjectsCompatibilityProfileChecker( const std::string& name, ISvcLocator* pSvcLocator )
      : Consumer( name, pSvcLocator,
                  {KeyValue{"TrackContainerHLT1", LHCb::TrackLocation::Default},
                   KeyValue{"TrackContainerHLT2", LHCb::TrackLocation::Default},
                   KeyValue{"MCParticleInput", LHCb::MCParticleLocation::Default},
                   KeyValue{"LinkerLocationHLT1", "Link/Pr/LHCbID"}, KeyValue{"LinkerLocationHLT2", "Link/Pr/LHCbID"},
                   KeyValue{"PVContainerHLT1", LHCb::RecVertexLocation::Primary},
                   KeyValue{"PVContainerHLT2", LHCb::RecVertexLocation::Primary}} ) {}
  void operator()( const LHCb::Tracks& hlt1_tracks, const LHCb::Tracks& hlt2_track,
                   const LHCb::MCParticles& mcparticles, const LHCb::LinksByKey& linker_hlt1,
                   const LHCb::LinksByKey& linker_hlt2, const LHCb::RecVertices& hlt1_pvs,
                   const LHCb::RecVertices& hlt2_pvs ) const override;

  mutable Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, int> m_hlt2_fitted_track{
      this, "hlt2_fitted_track", "hlt2 fitted track", {2, 0, 2}};
  mutable Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, int> m_hlt1_fitted_track{
      this, "hlt1_fitted_track", "hlt1 fitted track", {2, 0, 2}};

  mutable Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, int> m_numPVs_hlt1{
      this, "number_reconstructed_PVs_hlt1", "number_reconstructed_PVs_hlt1", {20, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, int> m_numPVs_hlt2{
      this, "number_reconstructed_PVs_hlt2", "number_reconstructed_PVs_hlt2", {20, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1> m_probChi2_hlt1{this, "probChi2_hlt1", "probChi2_hlt1", {50, 0, 1.}};
  mutable Gaudi::Accumulators::Histogram<1> m_probChi2_hlt2{this, "probChi2_hlt2", "probChi2_hlt2", {50, 0, 1.}};

  mutable Gaudi::Accumulators::Histogram<1> m_chi2_hlt1{this, "chi2_hlt1", "chi2_hlt1", {100, 0, 200}};
  mutable Gaudi::Accumulators::Histogram<1> m_chi2_hlt2{this, "chi2_hlt2", "chi2_hlt2", {100, 0, 200}};

  mutable Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, int> m_ndof_hlt1{
      this, "ndof_hlt1", "ndof_hlt1", {30, 0, 30}};
  mutable Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, int> m_ndof_hlt2{
      this, "ndof_hlt2", "ndof_hlt2", {30, 0, 30}};

  mutable Gaudi::Accumulators::Histogram<1> m_chi2_ndof_hlt1{this, "chi2_ndof_hlt1", "chi2_ndof_hlt1", {100, 0, 50}};
  mutable Gaudi::Accumulators::Histogram<1> m_chi2_ndof_hlt2{this, "chi2_ndof_hlt2", "chi2_ndof_hlt2", {100, 0, 50}};

  mutable Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, int> m_type_hlt1{
      this, "track_type_hlt1", "track_type_hlt1", {10, 0, 10}};
  mutable Gaudi::Accumulators::Histogram<1, Gaudi::Accumulators::atomicity::full, int> m_type_hlt2{
      this, "track_type_hlt2", "track_type_hlt2", {10, 0, 10}};

  mutable Gaudi::Accumulators::Histogram<1> m_deltaX{
      this, "deltaX_first_state", "hlt2 - hlt1 x firstate", {100, -0.2, 0.2}};
  mutable Gaudi::Accumulators::Histogram<1> m_deltaY{
      this, "deltaY_first_state", "hlt2 - hlt1 y firstate", {100, -0.2, 0.2}};
  mutable Gaudi::Accumulators::Histogram<1> m_deltaZ{
      this, "deltaZ_first_state", "hlt2 - hlt1 z firstate", {100, -2., 2.}};
  mutable Gaudi::Accumulators::Histogram<1> m_deltaTX{
      this, "deltaTX_first_state", "hlt2 - hlt1 tx firstate", {100, -0.005, 0.005}};
  mutable Gaudi::Accumulators::Histogram<1> m_deltaTY{
      this, "deltaTY_first_state", "hlt2 - hlt1 ty firstate", {100, -0.005, 0.005}};
  mutable Gaudi::Accumulators::Histogram<1> m_deltaQOP{
      this, "deltaQOP_first_state", "hlt2 - hlt1 qop firstate", {100, -0.00005, 0.00005}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_2d_X{
      this, "x_first_state_2d", "firstate_x;hlt1;hlt2", {100, -1., 1.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_2d_Y{
      this, "y_first_state_2d", "firstate_y;hlt1;hlt2", {100, -1., 1.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_2d_Z{
      this, "z_first_state_2d", "firstate_z;hlt1;hlt2", {100, -200., 600.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_2d_TX{
      this, "tx_first_state_2d", "firstate_tx;hlt1;hlt2", {100, -0.3, 0.3}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_2d_TY{
      this, "ty_first_state_2d", "firstate_ty;hlt1;hlt2", {100, -0.2, 0.2}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_2d_QOP{
      this, "qop_first_state_2d", "firstate_qop;hlt1;hlt2", {100, -0.001, 0.001}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_2d_PT{
      this, "pt_first_state_2d", "firstate_pt;hlt1;hlt2", {100, 0., 10000.}};

  mutable Gaudi::Accumulators::Histogram<1> m_deltaIP3D{this, "recdeltaIP3D", "delta_recIP3D", {100, -0.5, 0.5}};
  mutable Gaudi::Accumulators::Histogram<1> m_deltaIPx{this, "recdeltaIPx", "delta_recIPx", {100, -0.5, 0.5}};
  mutable Gaudi::Accumulators::Histogram<1> m_deltaIPy{this, "recdeltaIPy", "delta_recIPy", {100, -0.5, 0.5}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_IP3D{this, "recIP3D", "recIP3D;hlt1;hlt2", {100, 0., 1.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_IPx{this, "recIPx", "recIPx;hlt1;hlt2", {100, -1., 1.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_IPy{this, "recIPy", "recIPy;hlt1;hlt2", {100, -1., 1.}};

  mutable Gaudi::Accumulators::Histogram<1> m_deltaIPxErr{this, "delta_recIPxErr", "delta_recIPxErr", {100, -0.2, 0.2}};
  mutable Gaudi::Accumulators::Histogram<1> m_deltaIPyErr{this, "delta_recIPyErr", "delta_recIPyErr", {100, -0.2, 0.2}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_IPxErr{this, "recIPxErr", "recIPxErr;hlt1;hlt2", {100, 0., 0.5}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_IPyErr{this, "recIPyErr", "recIPyErr;hlt1;hlt2", {100, 0., 0.5}};

  mutable Gaudi::Accumulators::Histogram<1> m_deltaIPChi2{this, "delta_recIPChi2", "delta_recIPChi2", {100, -10., 10.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_IPChi2{this, "recIPChi2", "recIPChi2;hlt1;hlt2", {40, 0., 20.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_IPChi2_zoom{
      this, "recIPChi2_zoom", "recIPChi2_zoom;hlt1;hlt2", {20, 0., 2.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1, Gaudi::Accumulators::atomicity::full, int> m_nTracksPerPV{
      this, "nTracksPerPV", "nTracksPV;hlt1;hlt2", {50, 0, 50}};
  mutable Gaudi::Accumulators::ProfileHistogram<1, Gaudi::Accumulators::atomicity::full, int> m_nTracks{
      this, "nTracks", "nTracks;hlt1;hlt2", {50, 0, 100}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_X{this, "x_PV", "PVx;hlt1;hlt2", {100, -0.5, 0.5}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_Y{this, "y_PV", "PVy;hlt1;hlt2", {100, -0.5, 0.5}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_Z{this, "z_PV", "PVz;hlt1;hlt2", {100, -200., 200.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_cov00{this, "cov00_PV", "PVcov00;hlt1;hlt2", {50, 0., 0.0001}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_cov01{this, "cov01_PV", "PVcov01;hlt1;hlt2", {50, -0.05, 0.05}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_cov02{this, "cov02_PV", "PVcov02;hlt1;hlt2", {50, -0.05, 0.05}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_cov11{this, "cov11_PV", "PVcov11;hlt1;hlt2", {50, -0.1, 0.1}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_cov12{this, "cov12_PV", "PVcov12;hlt1;hlt2", {50, -0.1, 0.1}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_PV_cov22{this, "cov22_PV", "PVcov12;hlt1;hlt2", {50, 0., 1.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_X{
      this, "x_stAtVtx", "stAtVtx_x;hlt1;hlt2", {100, -1., 1.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_Y{
      this, "y_stAtVtx", "stAtVtx_y;hlt1;hlt2", {100, -1., 1.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_Z{
      this, "z_stAtVtx", "stAtVtx_z;hlt1;hlt2", {100, -200., 200.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_TX{
      this, "tx_stAtVtx", "stAtVtx_tx;hlt1;hlt2", {100, -0.3, 0.3}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_TY{
      this, "ty_stAtVtx", "stAtVtx_ty;hlt1;hlt2", {100, -0.2, 0.2}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_QOP{
      this, "qop_stAtVtx", "stAtVtx_qop;hlt1;hlt2", {100, -0.001, 0.001}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov00{
      this, "cov00_stAtVtx", "stAtVtxcov00;hlt1;hlt2", {50, 0., 0.01}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov01{
      this, "cov01_stAtVtx", "stAtVtxcov01;hlt1;hlt2", {50, -1, 10.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov02{
      this, "cov02_stAtVtx", "stAtVtxcov02;hlt1;hlt2", {50, -500., 500.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov03{
      this, "cov03_stAtVtx", "stAtVtxcov03;hlt1;hlt2", {50, -10., 10.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov04{
      this, "cov04_stAtVtx", "stAtVtxcov04;hlt1;hlt2", {50, -10., 10.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov11{
      this, "cov11_stAtVtx", "stAtVtxcov11;hlt1;hlt2", {50, 0., 10.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov12{
      this, "cov12_stAtVtx", "stAtVtxcov12;hlt1;hlt2", {50, -500., 500.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov13{
      this, "cov13_stAtVtx", "stAtVtxcov13;hlt1;hlt2", {50, -1., 1.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov14{
      this, "cov14_stAtVtx", "stAtVtxcov14;hlt1;hlt2", {50, -1., 0.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov22{
      this, "cov22_stAtVtx", "stAtVtxcov22;hlt1;hlt2", {50, -500., 500.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov23{
      this, "cov23_stAtVtx", "stAtVtxcov23;hlt1;hlt2", {50, -500., 500.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov24{
      this, "cov24_stAtVtx", "stAtVtxcov24;hlt1;hlt2", {50, -500., 500.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov33{
      this, "cov33_stAtVtx", "stAtVtxcov33;hlt1;hlt2", {50, 0., 500.}};
  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov34{
      this, "cov34_stAtVtx", "stAtVtxcov34;hlt1;hlt2", {50, -1000., 1000.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1> m_stAtVtx_cov44{
      this, "cov44_stAtVtx", "stAtVtxcov44;hlt1;hlt2", {50, 0., 800.}};

  mutable Gaudi::Accumulators::ProfileHistogram<1, Gaudi::Accumulators::atomicity::full, int> m_nVPHits{
      this, "nVPHits", "nVPHits;hlt1;hlt2", {21, 0, 21}};
  mutable Gaudi::Accumulators::ProfileHistogram<1, Gaudi::Accumulators::atomicity::full, int> m_nUTHits{
      this, "nUTHits", "nUTHits;hlt1;hlt2", {4, 0, 4}};
  mutable Gaudi::Accumulators::ProfileHistogram<1, Gaudi::Accumulators::atomicity::full, int> m_nFTHits{
      this, "nFTHits", "nFTHits;hlt1;hlt2", {16, 0, 16}};
};
// Difference between LHCb::Tracks and LHCb::Track::Range?
DECLARE_COMPONENT( TriggerObjectsCompatibilityProfileChecker )
void TriggerObjectsCompatibilityProfileChecker::
     operator()( const LHCb::Tracks& hlt1_tracks, const LHCb::Tracks& hlt2_tracks, const LHCb::MCParticles& mcparticles,
            const LHCb::LinksByKey& linker_hlt1, const LHCb::LinksByKey& linker_hlt2, const LHCb::RecVertices& hlt1_pvs,
            const LHCb::RecVertices& hlt2_pvs ) const {

  std::map<const LHCb::MCVertex*, int> truepvs;
  int                                  itrack_hlt2 = 0;
  for ( const auto& hlt2_track : hlt2_tracks ) {

    ++itrack_hlt2;
    const LHCb::MCParticle* hlt2_mcparticle = mcTruth( *hlt2_track, mcparticles, linker_hlt2 ); // CheckerBase
    auto                    fit             = fitResult( *hlt2_track );
    int                     trackWasFitted  = fit && !fit->nodes().empty();

    ++m_hlt2_fitted_track[trackWasFitted];

    if ( !selected( hlt2_mcparticle ) ) continue;
    if ( !( hlt2_track->hasVelo() ) ) continue;
    if ( int( hlt2_track->type() ) != 3 ) continue;

    const auto& hlt2_state = trackState( *hlt2_track );
    double      hlt2_x     = hlt2_state.x();
    double      hlt2_y     = hlt2_state.y();
    double      hlt2_z     = hlt2_state.z();
    double      hlt2_tx    = hlt2_state.tx();
    double      hlt2_ty    = hlt2_state.ty();
    double      hlt2_qop   = hlt2_state.qOverP();
    double      hlt2_pt    = hlt2_state.pt();

    ++m_numPVs_hlt2[int( hlt2_pvs.size() )];
    ++m_numPVs_hlt1[int( hlt1_pvs.size() )];
    m_nTracks[int( hlt1_tracks.size() )] += int( hlt2_tracks.size() );

    LHCb::CubicStateInterpolationTraj hlt2_tracktraj( hlt2_state, Gaudi::XYZVector() );
    Gaudi::XYZPoint                   hlt2_trkpos( hlt2_state.position() );
    Gaudi::XYZVector                  hlt2_trkdir( hlt2_state.slopes().Unit() );

    const LHCb::RecVertex* hlt2_recpv( nullptr );
    double                 hlt2_bestip2( 10000 ); // is this correct?

    for ( const auto& thispv : hlt2_pvs ) {
      Gaudi::XYZVector dx    = hlt2_trkpos - thispv->position();
      Gaudi::XYZVector delta = dx - hlt2_trkdir * dx.Dot( hlt2_trkdir );
      double           ip2   = delta.Mag2();
      if ( ( hlt2_recpv == 0 || ip2 < hlt2_bestip2 ) && thispv->isPrimary() ) {
        hlt2_bestip2 = ip2;
        hlt2_recpv   = thispv;
      }
    }

    for ( const auto& hlt1_track : hlt1_tracks ) {
      const LHCb::MCParticle* hlt1_mcparticle = mcTruth( *hlt1_track, mcparticles, linker_hlt1 );

      auto fit                 = fitResult( *hlt1_track );
      int  hlt1_trackWasFitted = fit && !fit->nodes().empty();

      ++m_hlt1_fitted_track[hlt1_trackWasFitted];

      if ( !selected( hlt1_mcparticle ) ) continue;

      if ( hlt2_mcparticle == hlt1_mcparticle ) {
        const auto& hlt1_state = trackState( *hlt1_track );
        double      hlt1_x     = hlt1_state.x();
        double      hlt1_y     = hlt1_state.y();
        double      hlt1_z     = hlt1_state.z();
        double      hlt1_tx    = hlt1_state.tx();
        double      hlt1_ty    = hlt1_state.ty();
        double      hlt1_qop   = hlt1_state.qOverP();
        double      hlt1_pt    = hlt1_state.pt();

        LHCb::CubicStateInterpolationTraj hlt1_tracktraj( hlt1_state, Gaudi::XYZVector() );
        Gaudi::XYZPoint                   hlt1_trkpos( hlt1_state.position() );
        Gaudi::XYZVector                  hlt1_trkdir( hlt1_state.slopes().Unit() );
        const LHCb::RecVertex*            hlt1_recpv( nullptr );
        double                            hlt1_bestip2( 10000 ); // is this correct?

        for ( const auto& thispv : hlt1_pvs ) {
          Gaudi::XYZVector dx    = hlt1_trkpos - thispv->position();
          Gaudi::XYZVector delta = dx - hlt1_trkdir * dx.Dot( hlt1_trkdir );
          double           ip2   = delta.Mag2();
          if ( hlt1_recpv == 0 || ip2 < hlt1_bestip2 ) {
            hlt1_bestip2 = ip2;
            hlt1_recpv   = thispv;
          }
        }

        ++m_probChi2_hlt1[hlt1_track->probChi2()];
        ++m_probChi2_hlt2[hlt2_track->probChi2()];

        ++m_chi2_hlt1[hlt1_track->chi2()];
        ++m_chi2_hlt2[hlt2_track->chi2()];

        ++m_ndof_hlt1[hlt1_track->nDoF()];
        ++m_ndof_hlt2[hlt2_track->nDoF()];

        ++m_chi2_ndof_hlt1[hlt1_track->chi2() / hlt1_track->nDoF()];
        ++m_chi2_ndof_hlt2[hlt2_track->chi2() / hlt2_track->nDoF()];

        ++m_type_hlt1[int( hlt1_track->type() )];
        ++m_type_hlt2[int( hlt2_track->type() )];

        ++m_deltaX[hlt2_x - hlt1_x];
        ++m_deltaY[hlt2_y - hlt1_y];
        ++m_deltaZ[hlt2_z - hlt1_z];
        ++m_deltaTX[hlt2_tx - hlt1_tx];
        ++m_deltaTY[hlt2_ty - hlt1_ty];
        ++m_deltaQOP[hlt2_qop - hlt1_qop];

        m_2d_X[hlt1_x] += hlt2_x;
        m_2d_Y[hlt1_y] += hlt2_y;
        m_2d_Z[hlt1_z] += hlt2_z;
        m_2d_TX[hlt1_tx] += hlt2_tx;
        m_2d_TY[hlt1_ty] += hlt2_ty;
        m_2d_QOP[hlt1_qop] += hlt2_qop;

        m_2d_PT[hlt1_pt] += hlt2_pt;
        m_nVPHits[hlt1_track->nVPHits()] += hlt2_track->nVPHits();
        m_nUTHits[hlt1_track->nUTHits()] += hlt2_track->nUTHits();
        m_nFTHits[hlt1_track->nFTHits()] += hlt2_track->nFTHits();

        if ( hlt2_recpv && hlt1_recpv ) {

          LHCb::State hlt2_stateAtVtx = hlt2_tracktraj.state( hlt2_recpv->position().z() );
          LHCb::State hlt1_stateAtVtx = hlt1_tracktraj.state( hlt1_recpv->position().z() );

          double hlt2_tx        = hlt2_stateAtVtx.tx();
          double hlt2_recipxerr = std::sqrt( hlt2_state.covariance()( 0, 0 ) + hlt2_recpv->covMatrix()( 0, 0 ) +
                                             2 * hlt2_tx * hlt2_recpv->covMatrix()( 0, 2 ) +
                                             hlt2_tx * hlt2_tx * hlt2_recpv->covMatrix()( 2, 2 ) );
          double hlt2_ty        = hlt2_stateAtVtx.ty();
          double hlt2_recipyerr = std::sqrt( hlt2_state.covariance()( 1, 1 ) + hlt2_recpv->covMatrix()( 1, 1 ) +
                                             2 * hlt2_ty * hlt2_recpv->covMatrix()( 1, 2 ) +
                                             hlt2_ty * hlt2_ty * hlt2_recpv->covMatrix()( 2, 2 ) );

          double hlt1_tx        = hlt1_stateAtVtx.tx();
          double hlt1_recipxerr = std::sqrt( hlt1_state.covariance()( 0, 0 ) + hlt1_recpv->covMatrix()( 0, 0 ) +
                                             2 * hlt1_tx * hlt1_recpv->covMatrix()( 0, 2 ) +
                                             hlt1_tx * hlt1_tx * hlt1_recpv->covMatrix()( 2, 2 ) );
          double hlt1_ty        = hlt1_stateAtVtx.ty();

          double hlt1_recipyerr = std::sqrt( hlt1_state.covariance()( 1, 1 ) + hlt1_recpv->covMatrix()( 1, 1 ) +
                                             2 * hlt1_ty * hlt1_recpv->covMatrix()( 1, 2 ) +
                                             hlt1_ty * hlt1_ty * hlt1_recpv->covMatrix()( 2, 2 ) );

          ++m_deltaIP3D[std::sqrt( hlt2_bestip2 ) - std::sqrt( hlt1_bestip2 )];
          ++m_deltaIPx[( hlt2_stateAtVtx.x() - hlt2_recpv->position().x() ) -
                       ( hlt1_stateAtVtx.x() - hlt1_recpv->position().x() )];
          ++m_deltaIPy[( hlt2_stateAtVtx.y() - hlt2_recpv->position().y() ) -
                       ( hlt1_stateAtVtx.y() - hlt1_recpv->position().y() )];

          ++m_deltaIPxErr[hlt2_recipxerr - hlt1_recipxerr];
          ++m_deltaIPyErr[hlt2_recipyerr - hlt1_recipyerr];

          m_IP3D[std::sqrt( hlt1_bestip2 )] += std::sqrt( hlt2_bestip2 );
          m_IPx[( hlt1_stateAtVtx.x() - hlt1_recpv->position().x() )] +=
              ( hlt2_stateAtVtx.x() - hlt2_recpv->position().x() );
          m_IPy[( hlt1_stateAtVtx.y() - hlt1_recpv->position().y() )] +=
              ( hlt2_stateAtVtx.y() - hlt2_recpv->position().y() );

          m_IPxErr[hlt1_recipxerr] += hlt2_recipxerr;
          m_IPyErr[hlt1_recipyerr] += hlt2_recipyerr;

          auto hlt2_bestipchi2 =
              LHCb::TrackVertexUtils::vertexChi2( hlt2_stateAtVtx, hlt2_recpv->position(), hlt2_recpv->covMatrix() );
          auto hlt1_bestipchi2 =
              LHCb::TrackVertexUtils::vertexChi2( hlt1_stateAtVtx, hlt1_recpv->position(), hlt1_recpv->covMatrix() );

          ++m_deltaIPChi2[hlt2_bestipchi2 - hlt1_bestipchi2];
          m_IPChi2[hlt1_bestipchi2] += hlt2_bestipchi2;
          m_IPChi2_zoom[hlt1_bestipchi2] += hlt2_bestipchi2;
          ++m_deltaIPChi2[( hlt1_bestipchi2 - hlt2_bestipchi2 )];

          m_nTracksPerPV[hlt1_recpv->tracks().size()] += hlt2_recpv->tracks().size();

          m_PV_X[hlt1_recpv->position().x()] += hlt2_recpv->position().x();
          m_PV_Y[hlt1_recpv->position().y()] += hlt2_recpv->position().y();
          m_PV_Z[hlt1_recpv->position().z()] += hlt2_recpv->position().z();

          m_PV_cov00[hlt1_recpv->covMatrix()( 0, 0 )] += hlt2_recpv->covMatrix()( 0, 0 );
          m_PV_cov01[hlt1_recpv->covMatrix()( 0, 1 )] += hlt2_recpv->covMatrix()( 0, 1 );
          m_PV_cov02[hlt1_recpv->covMatrix()( 0, 2 )] += hlt2_recpv->covMatrix()( 0, 2 );

          m_PV_cov11[hlt1_recpv->covMatrix()( 1, 1 )] += hlt2_recpv->covMatrix()( 1, 1 );
          m_PV_cov12[hlt1_recpv->covMatrix()( 1, 2 )] += hlt2_recpv->covMatrix()( 1, 2 );

          m_PV_cov22[hlt1_recpv->covMatrix()( 2, 2 )] += hlt2_recpv->covMatrix()( 2, 2 );

          m_stAtVtx_X[hlt1_stateAtVtx.position().x()] += hlt2_stateAtVtx.position().x();
          m_stAtVtx_Y[hlt1_stateAtVtx.position().y()] += hlt2_stateAtVtx.position().y();
          m_stAtVtx_Z[hlt1_stateAtVtx.position().z()] += hlt2_stateAtVtx.position().z();
          m_stAtVtx_TX[hlt1_stateAtVtx.tx()] += hlt2_stateAtVtx.tx();
          m_stAtVtx_TY[hlt1_stateAtVtx.ty()] += hlt2_stateAtVtx.ty();
          m_stAtVtx_QOP[hlt1_stateAtVtx.qOverP()] += hlt2_stateAtVtx.qOverP();

          m_stAtVtx_cov00[hlt1_stateAtVtx.posMomCovariance()( 0, 0 )] += hlt2_stateAtVtx.posMomCovariance()( 0, 0 );
          m_stAtVtx_cov01[hlt1_stateAtVtx.posMomCovariance()( 0, 1 )] += hlt2_stateAtVtx.posMomCovariance()( 0, 1 );
          m_stAtVtx_cov02[hlt1_stateAtVtx.posMomCovariance()( 0, 2 )] += hlt2_stateAtVtx.posMomCovariance()( 0, 2 );
          m_stAtVtx_cov03[hlt1_stateAtVtx.posMomCovariance()( 0, 3 )] += hlt2_stateAtVtx.posMomCovariance()( 0, 3 );
          m_stAtVtx_cov04[hlt1_stateAtVtx.posMomCovariance()( 0, 4 )] += hlt2_stateAtVtx.posMomCovariance()( 0, 4 );

          m_stAtVtx_cov11[hlt1_stateAtVtx.posMomCovariance()( 1, 1 )] += hlt2_stateAtVtx.posMomCovariance()( 1, 1 );
          m_stAtVtx_cov12[hlt1_stateAtVtx.posMomCovariance()( 1, 2 )] += hlt2_stateAtVtx.posMomCovariance()( 1, 2 );
          m_stAtVtx_cov13[hlt1_stateAtVtx.posMomCovariance()( 1, 3 )] += hlt2_stateAtVtx.posMomCovariance()( 1, 3 );
          m_stAtVtx_cov14[hlt1_stateAtVtx.posMomCovariance()( 1, 4 )] += hlt2_stateAtVtx.posMomCovariance()( 1, 4 );
          m_stAtVtx_cov22[hlt1_stateAtVtx.posMomCovariance()( 2, 2 )] += hlt2_stateAtVtx.posMomCovariance()( 2, 2 );
          m_stAtVtx_cov23[hlt1_stateAtVtx.posMomCovariance()( 2, 3 )] += hlt2_stateAtVtx.posMomCovariance()( 2, 3 );

          m_stAtVtx_cov24[hlt1_stateAtVtx.posMomCovariance()( 2, 4 )] += hlt2_stateAtVtx.posMomCovariance()( 2, 4 );

          m_stAtVtx_cov33[hlt1_stateAtVtx.posMomCovariance()( 3, 3 )] += hlt2_stateAtVtx.posMomCovariance()( 3, 3 );
          m_stAtVtx_cov34[hlt1_stateAtVtx.posMomCovariance()( 3, 4 )] += hlt2_stateAtVtx.posMomCovariance()( 3, 4 );

          m_stAtVtx_cov44[hlt1_stateAtVtx.posMomCovariance()( 4, 4 )] += hlt2_stateAtVtx.posMomCovariance()( 4, 4 );

        } // if reconstructd pv in both cases
        else {

          ++m_deltaIP3D[-999999.];
          ++m_deltaIPx[-999999.];
          ++m_deltaIPy[-999999.];

          m_IP3D[-999999.] += -999999.;
          m_IPx[-999999.] += -999999.;
          m_IPy[-999999.] += -999999.;

          m_IPxErr[-999999.] += -999999.;
          m_IPyErr[-999999.] += -999999.;

          ++m_deltaIPChi2[-999999.];
          m_IPChi2[-999999.] += -999999.;
          m_IPChi2_zoom[-999999.] += -999999.;
          ++m_deltaIPChi2[-99999.];

          m_nTracksPerPV[-9999.] += -9999.;
          m_PV_X[-9999.] += -9999.;
          m_PV_Y[-9999.] += -9999.;
          m_PV_Z[-9999.] += -9999.;

          m_PV_cov00[-9999.] += -9999.;
          m_PV_cov01[-9999.] += -9999.;
          m_PV_cov02[-9999.] += -9999.;

          m_PV_cov11[-9999.] += -9999.;
          m_PV_cov12[-9999.] += -9999.;

          m_PV_cov22[-9999.] += -9999.;

          m_stAtVtx_X[-9999.] += -9999.;
          m_stAtVtx_Y[-9999.] += -9999.;
          m_stAtVtx_Z[-9999.] += -9999.;
          m_stAtVtx_TX[-9999.] += -9999.;
          m_stAtVtx_TY[-9999.] += -9999.;
          m_stAtVtx_QOP[-9999.] += -9999.;

          m_stAtVtx_cov00[-9999.] += -9999.;
          m_stAtVtx_cov01[-9999.] += -9999.;
          m_stAtVtx_cov02[-9999.] += -9999.;
          m_stAtVtx_cov03[-9999.] += -9999.;
          m_stAtVtx_cov04[-9999.] += -9999.;

          m_stAtVtx_cov11[-9999.] += -9999.;
          m_stAtVtx_cov12[-9999.] += -9999.;
          m_stAtVtx_cov13[-9999.] += -9999.;
          m_stAtVtx_cov14[-9999.] += -9999.;

          m_stAtVtx_cov22[-9999.] += -9999.;
          m_stAtVtx_cov23[-9999.] += -9999.;
          m_stAtVtx_cov24[-9999.] += -9999.;

          m_stAtVtx_cov33[-9999.] += -9999.;
          m_stAtVtx_cov34[-9999.] += -9999.;

          m_stAtVtx_cov44[-9999.] += -9999.;

        } // if not resonctructed pvs

        // We also want to monitor the reconstructed IP, so the IP with respect to the reconstructed PVs.

      } // if mcparticles are the same for hlt1 and hlt2
    }   // loop over hlt1 tracks

  } // loop over hlt2 tracks
};
