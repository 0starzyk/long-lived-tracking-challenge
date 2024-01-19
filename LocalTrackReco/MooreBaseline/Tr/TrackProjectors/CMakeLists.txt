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
Tr/TrackProjectors
------------------
#]=======================================================================]

gaudi_add_module(TrackProjectors
    SOURCES
        src/TrackProjector.cpp
        src/TrackProjectorSelector.cpp
        src/TrajProjector.cpp
    LINK
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::LHCbKernel
        LHCb::TrackEvent
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib    
)
