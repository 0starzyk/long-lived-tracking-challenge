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
Tr/PatPV
--------
#]=======================================================================]

gaudi_add_module(PatPV
    SOURCES
        src/AdaptivePV3DFitter.cpp
        src/LSAdaptPV3DFitter.cpp
        src/LSAdaptPVFitter.cpp
        src/PVOfflineRecalculate.cpp
        src/PVOfflineTool.cpp
        src/PVSeed3DOfflineTool.cpp
        src/PVSeed3DTool.cpp
        src/PVSeedTool.cpp
        src/PVsEmptyProducer.cpp
        src/PatPV3D.cpp
        src/PatPV3DFuture.cpp
        src/PatPVOffline.cpp
        src/SimplePVFitter.cpp
        src/SimplePVSeedTool.cpp
        src/TrackBeamLineVertexFinder.cpp
        src/TrackBeamLineVertexFinderSoA.cpp
	src/TrackUnbiasedPVFinderSoA.cpp
	src/PatPVLeftRightFinderMerger.cpp
    LINK
        AIDA::aida
        Boost::container
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::DetDescLib
        LHCb::LHCbAlgsLib
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::RecEvent
        LHCb::TrackEvent
        Rec::PrKernel
        Rec::TrackInterfacesLib
	Rec::TrackKernel
)

gaudi_install(PYTHON)
