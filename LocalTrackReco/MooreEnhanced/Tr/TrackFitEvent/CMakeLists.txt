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
Tr/TrackFitEvent
----------------
#]=======================================================================]

gaudi_add_library(TrackFitEvent
    SOURCES
        src/FitNode.cpp
        src/KalmanFitResult.cpp
        src/Measurement.cpp
        src/PrKalmanFitResult.cpp
        src/TrackFitResult.cpp
        src/TriangularInverter.cpp
    LINK
        PUBLIC
            Gaudi::GaudiKernel
            LHCb::DetDescLib
            LHCb::FTDetLib
            LHCb::LHCbKernel
            LHCb::LHCbMathLib
            LHCb::TrackEvent
            LHCb::UTDetLib
            LHCb::VPDetLib
            LHCb::MuonDetLib
)

gaudi_add_executable(test_traj
    SOURCES
        tests/src/test_traj.cpp
    LINK
        Boost::unit_test_framework
        LHCb::LHCbMathLib
        TrackFitEvent
    TEST
)

gaudi_add_executable(test_measurement_minimize
    SOURCES
        tests/src/test_measurement_minimize.cpp
    LINK
        Boost::unit_test_framework
        TrackFitEvent
    TEST
)
