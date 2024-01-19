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
Tr/TrackMonitors
----------------
#]=======================================================================]

gaudi_add_module(TrackMonitors
    SOURCES
	    src/FTTrackMonitor.cpp
	src/FTHitEfficiencyMonitor.cpp
        src/PrSciFiHitsMonitor.cpp
        src/Map.cpp
        src/MonitorDetectorCorrelations.cpp
        src/TrackCaloMatchMonitor.cpp
        src/TrackDiMuonMonitor.cpp
        src/TrackFitMatchMonitor.cpp
        src/TrackMonitor.cpp
        src/TrackMonitorBase.cpp
        src/TrackMonitorTupleBase.cpp
        src/TrackMuonMatchMonitor.cpp
        src/TrackPV2HalfMonitor.cpp
        src/TrackParticleMonitor.cpp
        src/TrackTune.cpp
        src/TrackV0Monitor.cpp
        src/TrackVertexMonitor.cpp
        src/TrackVPOverlapMonitor.cpp
        src/UTTrackMonitor.cpp
        src/VPTrackMonitor.cpp
        src/VertexCompare.cpp
    LINK
        AIDA::aida
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        Gaudi::GaudiUtilsLib
        GSL::gsl
        LHCb::CaloDetLib
        LHCb::DAQEventLib
        LHCb::DetDescLib
        LHCb::DigiEvent
        LHCb::GenEvent
        LHCb::HltEvent
        LHCb::LHCbAlgsLib
        LHCb::LHCbKernel
        LHCb::LinkerEvent
        LHCb::MagnetLib
        LHCb::MCEvent
        LHCb::MCInterfaces
        LHCb::MuonDetLib
        LHCb::PartPropLib
        LHCb::PhysEvent
        LHCb::RecEvent
        LHCb::TrackEvent
        LHCb::UTDetLib
        LHCb::VPDetLib
        Rangev3::rangev3
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib
        Rec::TrackKernel
        ROOT::Hist
        ROOT::MathCore
)

gaudi_install(PYTHON)
