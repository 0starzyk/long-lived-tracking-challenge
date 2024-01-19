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
Tr/TrackKernel
--------------
#]=======================================================================]

gaudi_add_library(TrackKernel
    SOURCES
        src/CubicStateInterpolationTraj.cpp
        src/CubicStateVectorInterpolationTraj.cpp
        src/LineDifTraj.cpp
        src/StateTraj.cpp
        src/Track1DTabFunc.cpp
        src/TrackStateVertex.cpp
        src/TrackTraj.cpp
        src/TrajVertex.cpp
    LINK
        PUBLIC
            Gaudi::GaudiKernel
            LHCb::LHCbKernel
            LHCb::LHCbMathLib
            LHCb::RelationsLib
            LHCb::TrackEvent
            ROOT::MathCore
            Rec::TrackFitEvent
	    LHCb::MagnetLib
        PRIVATE
            Boost::headers
            GSL::gsl
            Rangev3::rangev3
)

gaudi_add_dictionary(TrackKernelDict
    HEADERFILES dict/TrackKernelDict.h
    SELECTION dict/TrackKernelDict.xml
    LINK TrackKernel
    OPTIONS ${REC_DICT_GEN_DEFAULT_OPTS}
)
