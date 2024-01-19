###############################################################################
# (c) Copyright 2020-2022 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
from Configurables import (LHCbApp, ExtrapolatorTester, ApplicationMgr,
                           TrackRungeKuttaExtrapolator, TrackKiselExtrapolator,
                           TrackHerabExtrapolator, TrackMasterExtrapolator,
                           TrackSimpleExtraSelector)
from PRConfig import TestFileDB
from DDDB.CheckDD4Hep import UseDD4Hep

app = LHCbApp()
app.EvtMax = 1
app.DataType = "Upgrade"
app.Simulation = True
all = []

if not UseDD4Hep:
    from Configurables import CondDB
    CondDB().Upgrade = True
else:
    from Configurables import LHCb__Tests__FakeRunNumberProducer as FET
    from Configurables import LHCb__Det__LbDD4hep__IOVProducer as IOVProducer
    odin_path = '/Event/DummyODIN'
    all = [
        FET('FakeRunNumber', ODIN=odin_path, Start=42, Step=20),
        IOVProducer("ReserveIOVDD4hep", ODIN=odin_path)
    ]

ex = ExtrapolatorTester()
ex.Extrapolators = []
add = lambda x: ex.Extrapolators.append(ex.addTool(x))
add(TrackRungeKuttaExtrapolator("Reference"))
add(
    TrackRungeKuttaExtrapolator(
        "BogackiShampine3", RKScheme="BogackiShampine3"))
add(TrackRungeKuttaExtrapolator("Verner7", RKScheme="Verner7"))
add(TrackRungeKuttaExtrapolator("Verner9", RKScheme="Verner9"))
add(
    TrackRungeKuttaExtrapolator(
        "Tsitouras5", RKScheme="Tsitouras5", OutputLevel=1))
add(TrackKiselExtrapolator("Kisel"))
add(TrackHerabExtrapolator("Herab"))
all.append(ex)

mgr = ApplicationMgr(TopAlg=all)
if UseDD4Hep:
    from Configurables import LHCb__Det__LbDD4hep__DD4hepSvc as DD4hepSvc
    dd4hep = DD4hepSvc()
    dd4hep.ConditionsVersion = "master"
    dd4hep.DetectorList = ["/world", "Magnet"]
    mgr.ExtSvc += [dd4hep]

TestFileDB.test_file_db['MiniBrunel_2018_MinBias_FTv4_DIGI'].run(
    configurable=app)
