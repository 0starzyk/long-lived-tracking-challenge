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
Pr/SciFiTrackForwarding
-----------------------
#]=======================================================================]

gaudi_add_module(SciFiTrackForwarding
    SOURCES
        src/SciFiTrackForwarding.cpp
        src/SciFiTrackForwardingStoreHit.cpp
    LINK
        Boost::container
        Boost::headers
        Gaudi::GaudiAlgLib
        LHCb::DetDescLib
        LHCb::FTDetLib
        LHCb::FTEvent
        LHCb::MagnetLib        
        LHCb::LHCbKernel
        LHCb::TrackEvent
        Rec::PrKernel
        VDT::vdt
)

# when building just Rec we can be more picky
# (within a superproject we get too many warnings from Gaudi)
if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    target_compile_options(SciFiTrackForwarding
        PRIVATE
            -Wdouble-promotion
            -Wredundant-decls
            -Wcast-align
            -Wmissing-declarations
            -Wswitch-enum
            -Winvalid-pch
            -Wformat=2
            -Wmissing-format-attribute
            -Wformat-nonliteral
            -Wnull-dereference
            -Wshadow
            -Wold-style-cast
            -Wcast-qual
    )

    # These flags only for gcc, they are not known by clang8
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(SciFiTrackForwarding
            PRIVATE
                -Wduplicated-cond
                -Wduplicated-branches
                -Wlogical-op
                -Wrestrict
                -Wuseless-cast
        )
    endif()

    # good for development, but too verbose for production
    # I'd like to have Winline, but the gaudi algorithm trigger to many warnings :/
    #target_compile_options(SciFiTrackForwarding PRIVATE -Wconversion -Wsign-conversion -Wfloat-equal -Winline)
endif()
