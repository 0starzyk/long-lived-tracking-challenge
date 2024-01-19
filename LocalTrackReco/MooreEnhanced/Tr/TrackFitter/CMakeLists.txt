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
Tr/TrackFitter
--------------
#]=======================================================================]

gaudi_add_module(TrackFitter
    SOURCES
        src/TrackEventFitter.cpp
        src/TrackKalmanFilter.cpp
        src/TrackMasterFitter.cpp
    LINK
        Boost::container
        Boost::headers
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::LHCbMathLib
        LHCb::TrackEvent
	LHCb::LHCbKernel
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib
        Rec::TrackKernel
)

gaudi_install(PYTHON)
