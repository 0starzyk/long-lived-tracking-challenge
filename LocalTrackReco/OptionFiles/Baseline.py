###############################################################################
# (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################

from PyConf.application import make_data_with_FetchDataFromFile
from Moore import options, run_reconstruction
from Moore.config import Reconstruction
from RecoConf.event_filters import require_gec


from RecoConf.hlt2_tracking import (make_hlt2_tracks)
from RecoConf.mc_checking import (
    check_tracking_efficiency, make_links_lhcbids_mcparticles_tracking_system,
    make_links_tracks_mcparticles, mc_unpackers)
from RecoConf.mc_checking_categories import (get_mc_categories,
                                             get_hit_type_mask)


from PyConf.Algorithms import (PrimaryVertexChecker,PrTrackerDumper,PrTrackRecoDumper,PrDownstreamDumper,PrGhostDumper)

from PyConf.application import make_odin

from RecoConf.hlt1_tracking import (make_VPClus_hits,make_PrStoreFTHit_hits,make_PrStoreUTHit_hits)

from RecoConf.hlt1_tracking import (make_pvs)


input_files = [
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000021_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000043_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000026_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000029_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000042_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000025_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000049_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000054_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000058_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000063_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000055_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000076_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000086_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000088_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000007_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000031_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000027_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000036_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000047_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000066_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000065_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000094_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000005_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000019_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000015_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000035_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000039_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000060_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000052_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000062_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000071_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000099_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000004_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000008_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000003_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000016_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000030_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000034_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000024_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000046_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000053_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000057_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000056_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000061_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000084_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000083_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000085_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000095_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000091_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000006_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000009_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000001_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000040_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000028_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000044_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000048_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000064_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000069_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000080_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000082_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000090_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000010_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000020_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000023_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000011_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000017_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000033_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000041_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000070_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000077_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000078_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000081_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000002_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000013_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000018_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000014_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000012_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000032_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000089_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000022_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000037_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000038_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000045_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000051_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000072_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000059_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000073_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000067_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000068_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000074_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000075_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000079_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000087_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000096_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000097_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000093_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000092_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00104521/0000/00104521_00000098_2.xdst',
## LS  Files Below:
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000005_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000006_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000044_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000029_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000043_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000048_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000019_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000046_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000068_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000061_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000076_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000088_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000013_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000024_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000053_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000028_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000058_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000035_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000030_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000025_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000062_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000074_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000079_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000089_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000055_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000042_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000023_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000039_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000066_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000009_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000001_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000059_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000016_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000040_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000015_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000018_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000041_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000050_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000036_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000065_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000070_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000064_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000067_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000073_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000087_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000092_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000002_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000038_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000054_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000060_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000049_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000027_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000047_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000075_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000003_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000051_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000031_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000045_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000063_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000071_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000078_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000083_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000084_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000086_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000007_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000008_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000004_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000014_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000026_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000022_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000033_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000072_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000077_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000081_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000080_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000037_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000020_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000056_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000057_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000034_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000032_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000012_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000021_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000017_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000052_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000069_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000082_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000091_2.xdst',
'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/Upgrade/XDST/00102230/0000/00102230_00000090_2.xdst',
]

options.input_files = input_files
options.input_type = 'ROOT'
options.input_raw_format = 4.3
# When running from Upgrade MC, must use the post-juggling locations of the raw
# event

# options.evt_max = 100
options.evt_max = 1000
options.simulation = True
options.data_type = 'Upgrade'
options.dddb_tag = 'dddb-20190223'
options.conddb_tag = 'sim-20180530-vc-md100'
options.geometry_version = 'trunk'
options.conditions_version = 'master'

options.output_file = 'hlt2_example.dst'
options.output_type = 'ROOT'
options.output_manifest_file = "hlt2_example.tck.json"

# options.histo_file = "MCMatching_baseline_MiniBias.root"
# options.ntuple_file = "ntuple.root"

def upgradeFunction():
    track_type = "Downstream"
    hlt2_tracks = make_hlt2_tracks()
    links_to_hits = make_links_lhcbids_mcparticles_tracking_system()
    links_to_tracks = make_links_tracks_mcparticles(
        InputTracks=hlt2_tracks[track_type], LinksToLHCbIDs=links_to_hits)
    pr_checker = check_tracking_efficiency(
        TrackType=track_type,
        InputTracks=hlt2_tracks[track_type],
        LinksToTracks=links_to_tracks,
        LinksToLHCbIDs=links_to_hits,
        MCCategories=get_mc_categories(track_type),
        HitTypesToCheck=get_hit_type_mask(track_type),
    )

    RecDumper=PrTrackRecoDumper(
            TrackLocation = hlt2_tracks["Downstream"]["v1"],
            VPLightClusterLocation = make_VPClus_hits(),
            ODINLocation = make_odin(),
            LinksLocation = links_to_tracks,
            FTHitsLocation = make_PrStoreFTHit_hits(),
            UTHitsLocation = make_PrStoreUTHit_hits(),
            MCParticlesLocation=mc_unpackers()["MCParticles"],
        )

    
    

    # data = [match_debug, forward_debug, forward_mc_debug, param_data]
    data = [pr_checker,RecDumper]


    return Reconstruction('upgradeFunction', data, [require_gec()])


options.histo_file = "output_hist.root"
options.ntuple_file = "output_tuple.root"
run_reconstruction(options, upgradeFunction)
