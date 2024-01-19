###############################################################################
# (c) Copyright 2000-2021 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
#[=======================================================================[.rst:
Pr/PrMCTools
------------
#]=======================================================================]

gaudi_add_module(PrMCTools
    SOURCES
        src/PrCheatedLongTracking.cpp
        src/PrCheatedSciFiTracking.cpp
        src/PrCheatedVP.cpp
        src/PrChecker.cpp
        src/PrMCDebugMatchToolNN.cpp
        src/PrMCDebugReconstructibleLong.cpp
        src/PrMCDebugForwardTool.cpp
        src/PrDebugTrackingTool.cpp
        src/PrDebugTrackingLosses.cpp
        src/PrLHCbID2MCParticle.cpp
        src/PrMultiplicityChecker.cpp
        src/PrPixelDebugTool.cpp
        src/PrSimpleTrackCounter.cpp
        src/PrTStationDebugTool.cpp
        src/PrTrackAssociator.cpp
        src/PrTrackConverter.cpp
        src/PrTrackCounter.cpp
        src/PrTrackRecoDumper.cpp
        src/PrTrackerDumper2.cpp
        src/PrUTCounter.cpp
        src/PrVPHitsChecker.cpp
        src/PrUTHitsChecker.cpp
        src/PrFTHitsChecker.cpp
        src/PrVeloUTChecker.cpp
        src/VPClusterEfficiency.cpp
        src/PrGhostDumper.cpp
        src/PrDownstreamDumper.cpp
    LINK
        Boost::filesystem
        Boost::headers
        fmt::fmt
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::AssociatorsBase
        LHCb::DAQEventLib
        LHCb::DetDescLib
        LHCb::DigiEvent
        LHCb::FTDetLib
        LHCb::FTEvent
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::LinkerEvent
        LHCb::LoKiCoreLib
        LHCb::LoKiMCLib
        LHCb::MCEvent
        LHCb::MCInterfaces
        LHCb::MuonDAQLib
        LHCb::MuonDetLib
        LHCb::PartPropLib
        LHCb::RecEvent
        LHCb::TrackEvent
        LHCb::UTDetLib
        LHCb::VPDetLib
        Rec::PrFitParamsLib
        Rec::PrKernel
        Rec::TrackInterfacesLib
        ROOT::Core
        ROOT::RIO
        ROOT::Tree
        )

if(USE_DD4HEP)
    target_link_libraries(PrMCTools PUBLIC Detector::DetectorLib)
endif()
