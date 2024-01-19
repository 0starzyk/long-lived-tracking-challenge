###############################################################################
# (c) Copyright 2000-2022 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
#[=======================================================================[.rst:
Tr/TrackInterfaces
------------------
#]=======================================================================]

gaudi_add_header_only_library(TrackInterfacesLib
    LINK
        Gaudi::GaudiKernel
        LHCb::DetDescLib
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::MagnetLib
        LHCb::PartPropLib
        LHCb::RecEvent
        LHCb::TrackEvent
        LHCb::VPDetLib
        Rec::TrackFitEvent
        ROOT::MathCore
)

gaudi_add_dictionary(TrackInterfacesDict
    HEADERFILES dict/TrackInterfacesDict.h
    SELECTION dict/TrackInterfacesDict.xml
    LINK TrackInterfacesLib
    OPTIONS ${REC_DICT_GEN_DEFAULT_OPTS}
)

gaudi_add_dictionary(TrackExtrapolatorDict
    HEADERFILES dict/TrackExtrapolatorDict.h
    SELECTION dict/TrackExtrapolatorDict.xml
    LINK TrackInterfacesLib
    OPTIONS ${REC_DICT_GEN_DEFAULT_OPTS}
)
