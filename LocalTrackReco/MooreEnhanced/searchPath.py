from LbEnv.ProjectEnv.options import (
    SearchPath,
    SearchPathEntry,
    EnvSearchPathEntry,
    NightlyPathEntry,
    LHCbDevPathEntry,
)

path = SearchPath([EnvSearchPathEntry('CMAKE_PREFIX_PATH', '/cvmfs/lhcb.cern.ch/lib/lhcb:/cvmfs/lhcb.cern.ch/lib/lcg/releases:/cvmfs/lhcb.cern.ch/lib/lcg/app/releases:/cvmfs/lhcb.cern.ch/lib/lcg/external:/cvmfs/lhcb.cern.ch/lib/contrib:/cvmfs/lhcb.cern.ch/lib/var/lib/LbEnv/2930/stable/linux-64/lib/python3.9/site-packages/LbDevTools/data/cmake'), EnvSearchPathEntry('CMTPROJECTPATH', ''), EnvSearchPathEntry('LHCBPROJECTPATH', '')])
