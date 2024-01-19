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
Pr/PrPixel
----------
#]=======================================================================]

gaudi_add_module(PrPixel
    SOURCES
        src/VPClus.cpp
        src/VPClusFull.cpp
        src/VeloClusterTrackingSIMD.cpp
        src/VeloKalman.cpp
    LINK
        Gaudi::GaudiAlgLib
        LHCb::DAQEventLib
        LHCb::DetDescLib
        LHCb::DigiEvent
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::TrackEvent
        LHCb::VPDetLib
        LHCb::VPKernel
        LHCb::VPDAQLib
        Rec::PrKernel
)
