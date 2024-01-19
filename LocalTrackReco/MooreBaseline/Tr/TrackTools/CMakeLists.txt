###############################################################################
# (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
#[=======================================================================[.rst:
Tr/TrackTools
-------------
#]=======================================================================]

gaudi_add_module(TrackTools
    SOURCES
	    src/FTTrackSelector.cpp
        src/CountVeloTracks.cpp
        src/DelegatingTrackSelector.cpp
        src/FastMomentumEstimate.cpp
        src/FTHitExpectation.cpp
        src/FTMeasurementProvider.cpp
        src/Line3D.cpp
        src/MeasurementProvider.cpp
        src/MeasurementProviderProjector.cpp
        src/MuonMeasurementProvider.cpp
        src/MuonSeeding.cpp
        src/MuonTrackSelector.cpp
        src/MuonUTTracking.cpp
        src/PtTransporter.cpp
        src/SeedMuonBuilder.cpp
        src/SelectTrackInVertex.cpp
        src/StandaloneMuonRec.cpp
        src/StateDetailedBetheBlochEnergyCorrectionTool.cpp
        src/StateElectronEnergyCorrectionTool.cpp
        src/StateSimpleBetheBlochEnergyCorrectionTool.cpp
        src/StateThickMSCorrectionTool.cpp
        src/StateThinMSCorrectionTool.cpp
        src/TrackCaloMatch.cpp
        src/TrackChi2Calculator.cpp
        src/TrackCloneFinder.cpp
        src/TrackInterpolator.cpp
        src/TrackMuonMatching.cpp
        src/TrackPtKick.cpp
        src/TrackSelector.cpp
        src/TrackUsedLHCbID.cpp
        src/TrackVertexer.cpp
        src/UpgradeGhostId.cpp
        src/UTHitExpectation.cpp
        src/UTMeasurementProvider.cpp
        src/VeloMuonBuilder.cpp
	    src/VPHitExpectation.cpp
        src/CalculateGhostProbability.cpp
        src/VPClusterPosition.cpp
        src/VPMeasurementProvider.cpp
        src/VPTrackSelector.cpp
    LINK
        Boost::container
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::CaloFutureInterfaces
        LHCb::CaloFutureUtils
        LHCb::DAQEventLib
        LHCb::DetDescLib
        LHCb::DigiEvent
        LHCb::FTDetLib
        LHCb::FTEvent
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::LHCbTrackInterfaces
        LHCb::LinkerEvent
        LHCb::MagnetLib
        LHCb::MCEvent
        LHCb::MCInterfaces
        LHCb::MuonDAQLib
        LHCb::MuonDetLib
        LHCb::PhysEvent
        LHCb::RecEvent
        LHCb::RelationsLib
        LHCb::TrackEvent
        LHCb::UTDetLib
        LHCb::UTKernelLib
        LHCb::VPDetLib
        Rec::MuonInterfacesLib
        Rec::PrKernel
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib
        Rec::TrackKernel
        ROOT::MathCore
        VDT::vdt
)
