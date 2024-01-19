###############################################################################
# (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
"""
Test the conversion of AoS-like tracks aka (v1::Track and v2::Track) among
themselves and to SoA track containers.
"""
from Configurables import ApplicationMgr, LHCbApp, UniqueIDGeneratorAlg, UnpackTrack
from Configurables import (
    LHCb__Converters__Track__v2__fromV1TrackV2Track as fromV1TrackV2Track,
    LHCb__Converters__Track__v1__fromV2TrackV1Track as fromV2TrackV1Track,
    LHCb__Converters__Track__SOA__fromV1Track as fromV1TrackFittedGenericTrack)
from Gaudi.Configuration import ERROR
from PRConfig import TestFileDB

# Unpack tracks
track_unpacker = UnpackTrack()

# Track converters (old versions)
track_converter_v1_to_v2 = fromV1TrackV2Track('TrackConverter_v1_to_v2')
track_converter_v1_to_v2.InputTracksName = 'Rec/Track/Best'
track_converter_v1_to_v2.OutputTracksName = 'Rec/Track/BestV2'

track_converter_v2_to_v1 = fromV2TrackV1Track('TrackConverter_v2_to_v1')
track_converter_v2_to_v1.InputTracksName = 'Rec/Track/BestV2'
track_converter_v2_to_v1.OutputTracksName = 'Rec/Track/BestV1'

# Track converters (SIMD-friendly version)
track_converter_v1_to_SOA = fromV1TrackFittedGenericTrack(
    'TrackConverter_v1_to_SOA')
track_converter_v1_to_SOA.InputTracks = track_converter_v2_to_v1.OutputTracksName
track_converter_v1_to_SOA.OutputTracks = 'Rec/Track/BestFromV1SOA'
track_converter_v1_to_SOA.Relations = 'Rec/Track/BestFromV1SOARelations'
track_converter_v1_to_SOA.RestrictToType = 'Long'
track_converter_v1_to_SOA.OutputLevel = ERROR  # avoid invalid states warnings

# Define the sequence
ApplicationMgr().TopAlg = [
    UniqueIDGeneratorAlg(),  # so UniqueIDGenerator is built
    track_unpacker,
    track_converter_v1_to_v2,
    track_converter_v2_to_v1,
    # see issue #410
    # track_converter_v1_to_SOA
]

app = LHCbApp()
app.EvtMax = 100

# Pick a file that has the reconstruction available
f = TestFileDB.test_file_db["upgrade_minbias_hlt1_filtered"]
f.setqualifiers(configurable=app)
f.run(configurable=app, withDB=True)
