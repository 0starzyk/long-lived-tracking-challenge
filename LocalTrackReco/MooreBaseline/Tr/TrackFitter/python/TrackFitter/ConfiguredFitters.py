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
# ====================================================================
#  Track fitting options
# ====================================================================

from future.utils import raise_
from TrackSys.Configuration import TrackSys
from Gaudi.Configuration import allConfigurables

from Configurables import (
    TrackEventFitter, TrackMasterFitter, TrackKalmanFilter,
    TrackProjectorSelector, TrackMasterExtrapolator, TrackSimpleExtraSelector,
    TrackDistanceExtraSelector, TrackHerabExtrapolator,
    SimplifiedMaterialLocator, DetailedMaterialLocator, MeasurementProvider,
    StateDetailedBetheBlochEnergyCorrectionTool, StateThickMSCorrectionTool)


def ConfiguredMasterFitter(Name,
                           FieldOff=None,
                           SimplifiedGeometry=None,
                           NoDriftTimes=None,
                           KalmanSmoother=None,
                           LiteClusters=False,
                           ApplyMaterialCorrections=None,
                           StateAtBeamLine=True,
                           MaxNumberOutliers=2,
                           FastOutlierIteration=False,
                           MSRossiAndGreisen=False):
    from GaudiKernel.Configurable import ConfigurableAlgTool

    # set the mutable default arguments
    if FieldOff is None: FieldOff = TrackSys().fieldOff()
    if SimplifiedGeometry is None:
        SimplifiedGeometry = TrackSys().simplifiedGeometry()
    if NoDriftTimes is None: NoDriftTimes = TrackSys().noDrifttimes()
    if KalmanSmoother is None: KalmanSmoother = TrackSys().kalmanSmoother()
    if ApplyMaterialCorrections is None:
        ApplyMaterialCorrections = not TrackSys().noMaterialCorrections()

    #if isinstance(Name, ConfigurableAlgTool) :
    fitter = Name
    #else :
    #    if allConfigurables.get( Name ) :
    #        raise ValueError, 'ConfiguredMasterFitter: instance with name '+Name+' already exists'
    #    fitter = TrackMasterFitter(Name)

    # apply material corrections
    if not ApplyMaterialCorrections:
        fitter.ApplyMaterialCorrections = False
        fitter.Extrapolator.ApplyMultScattCorr = False
        fitter.Extrapolator.ApplyEnergyLossCorr = False
        fitter.Extrapolator.ApplyElectronEnergyLossCorr = False

    # provide a state at the beamline
    fitter.StateAtBeamLine = StateAtBeamLine

    # set the maximum number of outliers removed in the fit
    fitter.MaxNumberOutliers = MaxNumberOutliers

    # set up the material locator
    # Using 'RossiAndGreisen' will omit the log term in the multiple scattering calculation (which should give better pulls).
    if SimplifiedGeometry:
        fitter.addTool(SimplifiedMaterialLocator, name="MaterialLocator")
        fitter.Extrapolator.addTool(
            SimplifiedMaterialLocator, name="MaterialLocator")
    else:
        fitter.addTool(DetailedMaterialLocator, name="MaterialLocator")
        fitter.Extrapolator.addTool(
            DetailedMaterialLocator, name="MaterialLocator")

    fitter.MaterialLocator.addTool(StateThickMSCorrectionTool,
                                   "StateMSCorrectionTool")
    fitter.MaterialLocator.StateMSCorrectionTool.UseRossiAndGreisen = MSRossiAndGreisen
    fitter.Extrapolator.MaterialLocator.addTool(StateThickMSCorrectionTool,
                                                "StateMSCorrectionTool")
    fitter.Extrapolator.MaterialLocator.StateMSCorrectionTool.UseRossiAndGreisen = MSRossiAndGreisen

    # not yet used for DC09 production
    # else:
    #    fitter.addTool(DetailedMaterialLocator(), name="MaterialLocator")
    #    fitter.MaterialLocator.GeneralDedxTool="StateDetailedBetheBlochEnergyCorrectionTool"
    #    fitter.Extrapolator.addTool(DetailedMaterialLocator(), name="MaterialLocator")
    #    fitter.Extrapolator.GeneralDedxTool= "StateDetailedBetheBlochEnergyCorrectionTool"
    #    fitter.Extrapolator.addTool(DetailedMaterialLocator(), name="MaterialLocator")
    #    fitter.Extrapolator.MaterialLocator.GeneralDedxTool= "StateDetailedBetheBlochEnergyCorrectionTool"

    # special settings for field-off
    if FieldOff:
        fitter.Extrapolator.addTool(
            TrackSimpleExtraSelector, name="ExtraSelector")
        fitter.Extrapolator.ExtraSelector.ExtrapolatorName = "TrackLinearExtrapolator"
        fitter.Extrapolator.ApplyEnergyLossCorr = False
        fitter.Extrapolator.ApplyElectronEnergyLossCorr = False
        fitter.ApplyEnergyLossCorr = False
        fitter.NodeFitter.DoF = 4

    # change the smoother
    if KalmanSmoother:
        fitter.NodeFitter.ForceBiDirectionalFit = False

    if FastOutlierIteration:
        fitter.UpdateReferenceInOutlierIterations = False

    return fitter


def ConfiguredEventFitter(Name,
                          TracksInContainer,
                          TracksOutContainer,
                          FieldOff=None,
                          SimplifiedGeometry=None,
                          NoDriftTimes=None,
                          KalmanSmoother=None,
                          LiteClusters=False,
                          ApplyMaterialCorrections=None,
                          StateAtBeamLine=True,
                          MaxNumberOutliers=2,
                          MSRossiAndGreisen=False):
    # make sure the name is unique
    if allConfigurables.get(Name):
        raise_(
            ValueError, 'ConfiguredEventFitter: instance with name ' + Name +
            ' already exists')
    # create the event fitter
    eventfitter = TrackEventFitter(Name)
    eventfitter.TracksInContainer = TracksInContainer
    eventfitter.TracksOutContainer = TracksOutContainer
    # this addTool should not be necessary but unfortunately there is a problem with the toolhandle configuration
    eventfitter.addTool(TrackMasterFitter, name="Fitter")
    # configure the fitter
    ConfiguredMasterFitter(
        eventfitter.Fitter,
        FieldOff=FieldOff,
        SimplifiedGeometry=SimplifiedGeometry,
        NoDriftTimes=NoDriftTimes,
        KalmanSmoother=KalmanSmoother,
        LiteClusters=LiteClusters,
        ApplyMaterialCorrections=ApplyMaterialCorrections,
        StateAtBeamLine=StateAtBeamLine,
        MaxNumberOutliers=MaxNumberOutliers,
        MSRossiAndGreisen=MSRossiAndGreisen)
    return eventfitter


def ConfiguredHltFitter(Name):
    fitter = ConfiguredMasterFitter(
        Name, FieldOff=None, SimplifiedGeometry=True, LiteClusters=True)
    fitter.NodeFitter.ForceBiDirectionalFit = True
    fitter.NumberFitIterations = 1
    fitter.MaxNumberOutliers = 2
    fitter.UpdateTransport = False
    fitter.AddDefaultReferenceNodes = True
    fitter.UpdateReferenceInOutlierIterations = False
    ## Force to use drift time
    fitter.MeasProvider.IgnoreMuon = True
    return fitter


def ConfiguredForwardFitter(Name,
                            FieldOff=None,
                            LiteClusters=True,
                            ForceUseDriftTime=True):
    fitter = ConfiguredMasterFitter(
        Name,
        FieldOff=FieldOff,
        SimplifiedGeometry=True,
        LiteClusters=LiteClusters)
    fitter.NumberFitIterations = 1
    fitter.MaxNumberOutliers = 0
    fitter.AddDefaultReferenceNodes = False
    fitter.NodeFitter.ForceBiDirectionalFit = False
    fitter.FillExtraInfo = False
    fitter.UpdateReferenceInOutlierIterations = False
    return fitter


def ConfiguredForwardEventFitter(Name,
                                 TracksInContainer,
                                 TracksOutContainer,
                                 ForceUseDriftTime=True):
    eventfitter = TrackEventFitter(Name)
    eventfitter.TracksInContainer = TracksInContainer
    eventfitter.TracksOutContainer = TracksOutContainer
    fittername = Name + ".Fitter"
    eventfitter.addTool(
        ConfiguredForwardFitter(
            Name=fittername, ForceUseDriftTime=ForceUseDriftTime),
        name="Fitter")
    return eventfitter


def ConfiguredHltEventFitter(Name, TracksInContainer, TracksOutContainer):
    eventfitter = TrackEventFitter(Name)
    eventfitter.TracksInContainer = TracksInContainer
    eventfitter.TracksOutContainer = TracksOutContainer
    fittername = Name + ".Fitter"
    eventfitter.addTool(ConfiguredHltFitter(Name=fittername), name="Fitter")
    return eventfitter


def ConfiguredVeloOnlyEventFitter(Name, TracksInContainer, TracksOutContainer):
    eventfitter = ConfiguredFitter(
        Name, TracksInContainer, TracksOutContainer, FieldOff=True)
    eventfitter.Fitter.MeasProvider.IgnoreIT = True
    eventfitter.Fitter.MeasProvider.IgnoreOT = True
    eventfitter.Fitter.MeasProvider.IgnoreTT = True
    eventfitter.Fitter.MeasProvider.IgnoreMuon = True
    return eventfitter


def ConfiguredStraightLineFitter(Name,
                                 TracksInContainer,
                                 TracksOutContainer,
                                 NoDriftTimes=None):
    fitter = ConfiguredFitter(
        Name,
        TracksInContainer,
        TracksOutContainer,
        FieldOff=True,
        NoDriftTimes=NoDriftTimes,
        StateAtBeamLine=False,
        ApplyMaterialCorrections=False)
    fitter.Fitter.AddDefaultReferenceNodes = False
    return eventfitter


def ConfiguredForwardStraightLineFitter(Name):
    fitter = ConfiguredForwardFitter(Name, FieldOff=True)
    fitter.ApplyMaterialCorrections = False
    return fitter


def ConfiguredForwardStraightLineEventFitter(Name, TracksInContainer,
                                             TracksOutContainer):
    eventfitter = TrackEventFitter(Name)
    eventfitter.TracksInContainer = TracksInContainer
    eventfitter.TracksOutContainer = TracksOutContainer
    fittername = Name + ".Fitter"
    eventfitter.addTool(
        ConfiguredForwardStraightLineFitter(Name=fittername), name="Fitter")
    return eventfitter


####################################################################
## Below are configurations for 'Fit' algorithms used in           #
## RecoTracking.py. Don't use these.                               #
####################################################################


def ConfiguredFitVelo(Name="FitVelo",
                      TracksInContainer="Rec/Track/PreparedVelo",
                      TracksOutContainer="Rec/Track/PreparedVelo"):
    # note that we ignore curvatue in velo. in the end that seems the
    # most sensible thing to do.
    eventfitter = ConfiguredEventFitter(Name, TracksInContainer,
                                        TracksOutContainer)
    return eventfitter


def ConfiguredFitVeloTT(Name="FitVeloTT",
                        TracksInContainer="Rec/Track/VeloTT",
                        TracksOutContainer="Rec/Track/VeloTT"):
    eventfitter = ConfiguredEventFitter(Name, TracksInContainer,
                                        TracksOutContainer)
    # Make the nodes even though this is a refit, to add StateAtBeamLine
    # (PatVeloTT also fits but does not make nodes)
    eventfitter.Fitter.MakeNodes = True
    return eventfitter


def ConfiguredFitSeed(Name="FitSeedForMatch",
                      TracksInContainer="Rec/Track/Seed",
                      TracksOutContainer="Rec/Track/Seed"):
    eventfitter = ConfiguredEventFitter(Name, TracksInContainer,
                                        TracksOutContainer)
    ## Such a small error is needed for TrackMatching to work properly,
    eventfitter.Fitter.ErrorQoP = [0.04, 5e-08]
    return eventfitter


def ConfiguredFitForward(Name="FitForward",
                         TracksInContainer="Rec/Track/Forward",
                         TracksOutContainer="Rec/Track/Forward"):
    eventfitter = ConfiguredEventFitter(Name, TracksInContainer,
                                        TracksOutContainer)
    return eventfitter


def ConfiguredFitMatch(Name="FitMatch",
                       TracksInContainer="Rec/Track/Match",
                       TracksOutContainer="Rec/Track/Match"):
    eventfitter = ConfiguredEventFitter(Name, TracksInContainer,
                                        TracksOutContainer)
    return eventfitter


def ConfiguredFitDownstream(Name="FitDownstream",
                            TracksInContainer="Rec/Track/Downstream",
                            TracksOutContainer="Rec/Track/Downstream"):
    eventfitter = ConfiguredEventFitter(Name, TracksInContainer,
                                        TracksOutContainer)
    return eventfitter


def ConfiguredFit(Name, TracksInContainer, TracksOutContainer):
    eventfitter = ConfiguredEventFitter(Name, TracksInContainer,
                                        TracksOutContainer)
    return eventfitter
