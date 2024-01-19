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
Pr/PrFitParams
--------------
#]=======================================================================]

gaudi_add_header_only_library(PrFitParamsLib
    LINK
        Gaudi::GaudiKernel
        ROOT::MathCore
)

gaudi_add_module(PrFitParams
    SOURCES
        src/PrFitKsParams.cpp
        src/PrFitMatchParams.cpp
        src/PrFitParameters.cpp
        src/PrFitSeedParams.cpp
        src/PrParameterisationData.cpp
        src/PrFitTool.cpp
    LINK
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::LHCbAlgsLib
        LHCb::LHCbKernel
        LHCb::LinkerEvent
        LHCb::MCEvent
        LHCb::TrackEvent
        PrFitParamsLib
        PrKernel
        TrackInterfacesLib
)
