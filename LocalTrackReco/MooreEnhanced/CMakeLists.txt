cmake_minimum_required(VERSION 3.15)

project(MooreEnhanced
    VERSION 0
    LANGUAGES CXX)

# Enable tests support
include(CTest)

# Load "lhcb_find_package" function
include(/cvmfs/lhcb.cern.ch/lib/var/lib/LbEnv/2930/stable/linux-64/lib/python3.9/site-packages/LbDevTools/data/cmake/LHCbFindPackage.cmake)

# Get the list of projects we are working on
execute_process(
    COMMAND git remote
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE enabled_projects
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
# (this is a trick to split a atring into a list of words)
string(REGEX MATCHALL "[a-zA-Z0-9]+" enabled_projects "${enabled_projects}")
# Set flags to inherit private dependencies of the projects we work on
foreach(p IN LISTS enabled_projects)
    set(WITH_${p}_PRIVATE_DEPENDENCIES TRUE)
endforeach()

list(PREPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/MooreEnhancedDependencies.cmake)

set(GAUDI_PREFER_LOCAL_TARGETS TRUE)
lhcb_add_subdirectories(
    # -- begin: list of subdirectories --
    Pr/PrAlgorithms
    Pr/PrConverters
    Pr/PrFitParams
    Pr/PrKernel
    Pr/PrMCTools
    Pr/PrPixel
    Pr/PrVPRetina
    Pr/PrVeloUT
    Pr/SciFiTrackForwarding
    Tr/ParameterizedKalman
    Tr/PatChecker
    Tr/PatFitParams
    Tr/PatPV
    Tr/PrKalmanFilter
    Tr/TrackAssociators
    Tr/TrackCheckers
    Tr/TrackExtrapolators
    Tr/TrackFitEvent
    Tr/TrackFitter
    Tr/TrackInterfaces
    Tr/TrackKernel
    Tr/TrackMCTools
    Tr/TrackMonitors
    Tr/TrackProjectors
    Tr/TrackTools
    Tr/TrackUtils
    # -- end: list of subdirectories --
)

# Final configuration steps
set(GAUDI_LEGACY_CMAKE_SUPPORT TRUE)
lhcb_finalize_configuration(NO_EXPORT)

# Enable Ganga integration
include(/cvmfs/lhcb.cern.ch/lib/var/lib/LbEnv/2930/stable/linux-64/lib/python3.9/site-packages/LbDevTools/data/cmake/GangaTools.cmake)
enable_ganga_integration()
