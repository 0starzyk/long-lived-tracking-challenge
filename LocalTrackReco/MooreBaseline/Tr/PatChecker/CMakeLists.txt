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
Tr/PatChecker
-------------
#]=======================================================================]

gaudi_add_module(PatChecker
    SOURCES
        src/CheatedPrimaryVertices.cpp
        src/PatKShortChecker.cpp
        src/PrimaryVertexChecker.cpp
    LINK
        AIDA::aida
        LHCb::VPDetLib
        fmt::fmt
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        Gaudi::GaudiUtilsLib
        GSL::gsl
        LHCb::DetDescLib
        LHCb::LHCbMathLib
        LHCb::LinkerEvent
        LHCb::MCEvent
        LHCb::MCInterfaces
        LHCb::PhysEvent
        LHCb::RecEvent
        LHCb::TrackEvent
)

 
        
        