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
Tr/TrackCheckers
----------------
#]=======================================================================]

gaudi_add_module(TrackCheckers
    SOURCES
        src/ExtrapolatorChecker.cpp
        src/Map.cpp
        src/TrackCheckerBase.cpp
        src/TrackCloneChecker.cpp
        src/TrackEffChecker.cpp
        src/TrackIPResolutionChecker.cpp
        src/TrackIPResolutionCheckerNT.cpp
        src/TriggerObjectsCompatibilityChecker.cpp
        src/TriggerObjectsCompatibilityProfileChecker.cpp 
        src/TrackResChecker.cpp
        src/VertexChecker.cpp
	src/TrackPVChecker.cpp
    LINK
        AIDA::aida
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        Gaudi::GaudiUtilsLib
        LHCb::DetDescLib
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::LinkerEvent
        LHCb::MCEvent
        LHCb::MCInterfaces
        LHCb::RecEvent
        LHCb::TrackEvent
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib
        Rec::TrackKernel
        ROOT::Hist
)
