###############################################################################
# (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
from Configurables import (GaudiSequencer, TrackMonitor, TrackVertexMonitor,
                           TrackFitMatchMonitor, TrackV0Monitor,
                           TrackDiMuonMonitor, TrackEcalMatchMonitor,
                           TrackMuonMatchMonitor, TrackPV2HalfMonitor)
from Configurables import (RecSysConf, RecMoniConf)


def ConfiguredTrackMonitorSequence(Name="TrackMonitorSequence",
                                   HistoPrint=False):

    # figure out detectors
    seq = GaudiSequencer(Name)
    subDets = None  #default, all monitors
    from Configurables import LHCbApp
    if hasattr(LHCbApp(), "Detectors"):
        if LHCbApp().isPropertySet("Detectors"):
            subDets = LHCbApp().getProp("Detectors")
    seq.Members.append(TrackMonitor())
    seq.Members.append(TrackDiMuonMonitor(HistoPrint=HistoPrint))
    seq.Members.append(TrackVertexMonitor(HistoPrint=HistoPrint))

    if RecMoniConf().getProp("Histograms") != "Online":
        seq.Members.append(
            TrackV0Monitor(
                HistoPrint=HistoPrint,
                RequireObjects=["/Event/Rec/Vertex/V0"]))
        seq.Members.append(TrackFitMatchMonitor())
        seq.Members.append(TrackPV2HalfMonitor())
        if "CALO" in RecSysConf().RecoSequence:
            if (subDets is None or "Ecal" in subDets):
                seq.Members.append(
                    TrackEcalMatchMonitor(HistoPrint=HistoPrint))
        if "MUON" in RecSysConf().RecoSequence:
            if (subDets is None or "Muon" in subDets):
                seq.Members.append(
                    TrackMuonMatchMonitor(
                        "TrackMuonMatchMonitor", HistoPrint=HistoPrint))

    return seq
