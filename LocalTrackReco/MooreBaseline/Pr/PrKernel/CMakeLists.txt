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
Pr/PrKernel
-----------
#]=======================================================================]

gaudi_add_library(PrKernel
    SOURCES
        src/FTMatsCache.cpp
    LINK
        PUBLIC
            Boost::container
            Boost::headers
            Gaudi::GaudiKernel
            LHCb::DetDescLib
            LHCb::DigiEvent
            LHCb::EventBase
            LHCb::FTDAQLib
            LHCb::FTDetLib
            LHCb::LHCbKernel
            LHCb::LHCbMathLib
            LHCb::TrackEvent
            LHCb::UTDAQLib
            LHCb::UTDetLib
)

gaudi_add_dictionary(PrKernelDict
    HEADERFILES dict/PrKernelDict.h
    SELECTION dict/PrKernelDict.xml
    LINK PrKernel
    OPTIONS ${REC_DICT_GEN_DEFAULT_OPTS}
)

gaudi_add_executable(TestSelection
    SOURCES
        tests/src/TestSelection.cpp
    LINK
        Boost::unit_test_framework
        LHCb::TrackEvent
        PrKernel
    TEST
)
