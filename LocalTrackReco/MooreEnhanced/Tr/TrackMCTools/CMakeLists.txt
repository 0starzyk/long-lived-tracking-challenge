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
Tr/TrackMCTools
---------------
#]=======================================================================]

gaudi_add_module(TrackMCTools
    SOURCES
        src/DownstreamGhostClassification.cpp
        src/IdealStateCreator.cpp
        src/LHCbIDsToMCHits.cpp
        src/LHCbIDsToMCParticles.cpp
        src/LongGhostClassification.cpp
        src/MeasureIPResolution.cpp
        src/TrackGhostClassificationBase.cpp
        src/UpgradeGhostIdNT.cpp
        src/UpstreamGhostClassification.cpp
        src/VeloGhostClassification.cpp
    LINK
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::DigiEvent
        LHCb::FTEvent
        LHCb::LHCbKernel
        LHCb::LinkerEvent
        LHCb::MCEvent
        LHCb::MCInterfaces
        LHCb::PartPropLib
        LHCb::RecEvent
        LHCb::TrackEvent
        Rec::TrackInterfacesLib
        ROOT::MathCore
)
