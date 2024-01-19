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
Tr/TrackUtils
-------------
#]=======================================================================]

gaudi_add_module(TrackUtils
    SOURCES
        src/CreateFastTrackCollection.cpp
        src/FilterDownstreamTracks.cpp
        src/FilterMatchTracks.cpp
        src/FilterSeedTracks.cpp
        src/HltInsertTrackErrParam.cpp
        src/MakeView.cpp
        src/MakeZipContainer.cpp
        src/MoveRecVertexTracks.cpp
        src/NumberOfTracksFilter.cpp
        src/PrCloneKiller.cpp
	src/PrKalmanFitResultMerger.cpp
        src/TTrackFromLong.cpp
        src/TrackAddNNGhostId.cpp
        src/TrackBestTrackCreator.cpp
        src/TrackCloneCleaner.cpp
        src/TrackCloneKiller.cpp
        src/TrackCompetition.cpp
        src/TrackContainerCleaner.cpp
        src/TrackContainerCopy.cpp
        src/TrackContainerSplitter.cpp
        src/TrackContainersMerger.cpp
        src/TrackEraseExtraInfo.cpp
        src/TrackEventCloneKiller.cpp
        src/TrackFromDST.cpp
        src/TrackListFilter.cpp
        src/TrackListMerger.cpp
        src/TrackListPrinter.cpp
        src/TrackListRefiner.cpp
        src/TrackSelectionToContainer.cpp
        src/TrackToDST.cpp
        src/TrackV0Finder.cpp
        src/TracksDownstreamConverter.cpp
        src/TracksEmptyProducer.cpp
        src/TracksFTConverter.cpp
        src/TracksFitConverter.cpp
        src/TracksMatchConverter.cpp
        src/TracksToSelection.cpp
        src/TracksUTConverter.cpp
        src/TracksVPConverter.cpp
        src/TracksVPMerger.cpp
        src/VertexListRefiner.cpp
    LINK
        Boost::headers
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
        Gaudi::GaudiPluginService
        GSL::gsl
        LHCb::DigiEvent
        LHCb::LHCbAlgsLib
        LHCb::LHCbKernel
        LHCb::LHCbMathLib
        LHCb::LinkerEvent
        LHCb::PartPropLib
        LHCb::PhysEvent
        LHCb::RecEvent
        LHCb::TrackEvent
        Rangev3::rangev3
        Rec::TrackFitEvent
        Rec::TrackInterfacesLib
        Rec::TrackKernel
        Rec::FunctorCoreLib
)

gaudi_add_executable(TestZipInfrastructure
    SOURCES
        tests/src/TestZipInfrastructure.cpp
    LINK
        Boost::unit_test_framework
        Gaudi::GaudiKernel
        LHCb::EventBase
        LHCb::RecEvent
        LHCb::TrackEvent
        Rec::SelKernelLib
    TEST
)
