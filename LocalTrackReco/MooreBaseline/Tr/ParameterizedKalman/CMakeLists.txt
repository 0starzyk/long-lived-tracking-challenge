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
Tr/ParameterizedKalman
----------------------
#]=======================================================================]

gaudi_add_module(ParameterizedKalman
    SOURCES
        src/CompareTracks.cpp
        src/KalmanParametrizations.cpp
        src/KalmanParametrizationsCoef.cpp
        src/ParameterizedKalmanFit.cpp
        src/ParameterizedKalmanFit_Checker.cpp
        src/ParameterizedKalmanFit_Methods.cpp
        src/SerializeTrack.cpp
    LINK
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::AssociatorsBase
        LHCb::DAQEventLib
        LHCb::FTDetLib
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::MagnetLib
        LHCb::MCEvent
        LHCb::MCInterfaces
        LHCb::TrackEvent
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib
        Rec::TrackKernel
        ROOT::RIO
        ROOT::Tree
)
