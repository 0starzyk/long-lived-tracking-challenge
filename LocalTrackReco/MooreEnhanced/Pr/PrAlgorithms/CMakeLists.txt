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
Pr/PrAlgorithms
---------------
#]=======================================================================]

gaudi_add_header_only_library(PrAlgorithmsLib
    LINK
        Gaudi::GaudiAlgLib
        LHCb::LHCbAlgsLib
        LHCb::LHCbKernel
        Rec::FunctorCoreLib
        Rec::SelKernelLib
)

gaudi_add_module(PrAlgorithms
    SOURCES
        src/PrAddUTHitsTool.cpp
        src/PrIgnoreUTHitsTool.cpp
        src/PrCheckEmptyFilter.cpp
        src/PrFilter.cpp
        src/PrFilterIPSoA.cpp
        src/PrForwardTracking.cpp
        src/PrGECFilter.cpp
        src/PrHybridSeeding.cpp
        src/PrLongLivedTracking.cpp
        src/PrMakeSelection.cpp
        src/PrMatchNN.cpp
        src/PrResidualPrUTHits.cpp
        src/PrResidualSciFiHits.cpp
        src/PrResidualSeeding.cpp
        src/PrResidualUTHits.cpp
        src/PrResidualVeloTracks.cpp
        src/PrStoreFTHit.cpp
        src/PrStoreSciFiHits.cpp
        src/PrStoreUTHit.cpp
        src/PrUpstreamFromVelo.cpp
        src/HitsEmptyProducer.cpp
    LINK
        PrAlgorithmsLib
        Boost::container
        Boost::headers
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        LHCb::LHCbAlgsLib
        LHCb::DAQEventLib
        LHCb::DetDescLib
        LHCb::FTDAQLib
        LHCb::FTDetLib
        LHCb::FTEvent
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::LoKiCoreLib
        LHCb::MagnetLib
        LHCb::PhysEvent
        LHCb::RecEvent
        LHCb::TrackEvent
        LHCb::UTDAQLib
        LHCb::UTDetLib
        LHCb::UTKernelLib
        LHCb::VPDetLib
        Rangev3::rangev3
        Rec::PrKernel
        Rec::TrackInterfacesLib
        ROOT::MathCore
        Vc::Vc
        VDT::vdt
)
