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
Tr/PrKalmanFilter
-----------------
#]=======================================================================]

gaudi_add_header_only_library(PrKalmanFilterLib
    LINK
        Boost::container
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::LHCbMathLib
        LHCb::TrackEvent
        Rec::PrKernel
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib
)

gaudi_add_module(PrKalmanFilter
    SOURCES
        src/KalmanFilter.cpp
        src/KalmanFilterTool.cpp
        src/ToolExample.cpp
    LINK
        PrKalmanFilterLib
        Boost::container
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::LHCbAlgsLib
        LHCb::LHCbMathLib
        LHCb::TrackEvent
        Rec::PrKernel
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib
)

