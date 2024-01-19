/*****************************************************************************\
* (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
// Class: ReadMLP_5
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : MLP::MLP
TMVA Release   : 4.2.1         [262657]
ROOT Release   : 6.06/08       [394760]
Creator        : mexu
Date           : Thu Jan 31 10:19:21 2019
Host           : Linux lcgapp-slc6-x86-64-2.cern.ch 2.6.32-504.1.3.el6.x86_64 #1 SMP Wed Nov 12 06:58:35 CET 2014 x86_64
x86_64 x86_64 GNU/Linux Dir            : /afs/cern.ch/work/m/mexu/workspace/Ghost_Study/Final_Verstion Training events:
52558 Analysis type  : [Classification]


#OPT -*-*-*-*-*-*-*-*-*-*-*-*- options -*-*-*-*-*-*-*-*-*-*-*-*-

# Set by User:
NCycles: "600" [Number of training cycles]
HiddenLayers: "N+5" [Specification of hidden layer architecture]
NeuronType: "ReLU" [Neuron activation function type]
EstimatorType: "CE" [MSE (Mean Square Estimator) for Gaussian Likelihood or CE(Cross-Entropy) for Bernoulli Likelihood]
V: "False" [Verbose output (short form of "VerbosityLevel" below - overrides the latter one)]
VarTransform: "N" [List of variable transformations performed before training, e.g.,
"D_Background,P_Signal,G,N_AllClasses" for: "Decorrelation, PCA-transformation, Gaussianisation, Normalisation, each for
the given class of events ('AllClasses' denotes all events of all classes, if no class indication is given, 'All' is
assumed)"] H: "False" [Print method-specific help message] CreateMVAPdfs: "True" [Create PDFs for classifier outputs
(signal and background)] TestRate: "5" [Test for overtraining performed at each #th epochs] UseRegulator: "False" [Use
regulator to avoid over-training] # Default: RandomSeed: "1" [Random seed for initial synapse weights (0 means unique
seed for each run; default value '1')] NeuronInputType: "sum" [Neuron input function type] VerbosityLevel: "Default"
[Verbosity level] IgnoreNegWeightsInTraining: "False" [Events with negative weights are ignored in the training (but are
included for testing and performance evaluation)] TrainingMethod: "BP" [Train with Back-Propagation (BP), BFGS Algorithm
(BFGS), or Genetic Algorithm (GA - slower and worse)] LearningRate: "2.000000e-02" [ANN learning rate parameter]
DecayRate: "1.000000e-02" [Decay rate for learning parameter]
EpochMonitoring: "False" [Provide epoch-wise monitoring plots according to TestRate (caution: causes big ROOT output
file!)] Sampling: "1.000000e+00" [Only 'Sampling' (randomly selected) events are trained each epoch] SamplingEpoch:
"1.000000e+00" [Sampling is used for the first 'SamplingEpoch' epochs, afterwards, all events are taken for training]
SamplingImportance: "1.000000e+00" [ The sampling weights of events in epochs which successful (worse estimator than
before) are multiplied with SamplingImportance, else they are divided.] SamplingTraining: "True" [The training sample is
sampled] SamplingTesting: "False" [The testing sample is sampled] ResetStep: "50" [How often BFGS should reset history]
Tau: "3.000000e+00" [LineSearch "size step"]
BPMode: "sequential" [Back-propagation learning mode: sequential or batch]
BatchSize: "-1" [Batch size: number of events/batch, only set if in Batch Mode, -1 for BatchSize=number_of_events]
ConvergenceImprove: "1.000000e-30" [Minimum improvement which counts as improvement (<0 means automatic convergence
check is turned off)] ConvergenceTests: "-1" [Number of steps (without improvement) required for convergence (<0 means
automatic convergence check is turned off)] UpdateLimit: "10000" [Maximum times of regulator update] CalculateErrors:
"False" [Calculates inverse Hessian matrix at the end of the training to be able to calculate the uncertainties of an
MVA value] WeightRange: "1.000000e+00" [Take the events for the estimator calculations from small deviations from the
desired value to large deviations only over the weight range]
##


#VAR -*-*-*-*-*-*-*-*-*-*-*-* variables *-*-*-*-*-*-*-*-*-*-*-*-

NVar 10
UpgradeGhostInfo_obsFT        UpgradeGhostInfo_obsFT        UpgradeGhostInfo_obsFT        UpgradeGhostInfo_obsFT 'F'
[9,12] UpgradeGhostInfo_FitTChi2     UpgradeGhostInfo_FitTChi2     UpgradeGhostInfo_FitTChi2 UpgradeGhostInfo_FitTChi2
'F'    [0.0107197267935,29.6306381226] UpgradeGhostInfo_FitTNDoF     UpgradeGhostInfo_FitTNDoF UpgradeGhostInfo_FitTNDoF
UpgradeGhostInfo_FitTNDoF                                       'F'    [2,7] UpgradeGhostInfo_obsUT
UpgradeGhostInfo_obsUT        UpgradeGhostInfo_obsUT        UpgradeGhostInfo_obsUT 'F'    [2,8]
UpgradeGhostInfo_UToutlier    UpgradeGhostInfo_UToutlier    UpgradeGhostInfo_UToutlier    UpgradeGhostInfo_UToutlier 'F'
[0,2] UpgradeGhostInfo_veloHits     UpgradeGhostInfo_veloHits     UpgradeGhostInfo_veloHits UpgradeGhostInfo_veloHits
'F'    [212,8281] UpgradeGhostInfo_utHits       UpgradeGhostInfo_utHits       UpgradeGhostInfo_utHits
UpgradeGhostInfo_utHits                                         'F'    [218,5250] TRACK_CHI2 TRACK_CHI2 TRACK_CHI2
TRACK_CHI2                                                      'F'    [0.579863786697,38.7384567261] TRACK_PT TRACK_PT
TRACK_PT                      TRACK_PT                                                        'F'
[1.11619496346,11829.7861328] TRACK_ETA                     TRACK_ETA                     TRACK_ETA TRACK_ETA 'F'
[1.52495276928,8.42633724213] NSpec 0


============================================================================ */
#include "Kernel/STLExtensions.h"
#include "Kernel/TMV_utils.h"
#include "vdt/exp.h"
#include <array>
#include <string_view>

namespace Data::ReadGhostProbabilityDownstream {
  namespace {
    constexpr auto ActivationFnc       = []( float x ) { return x > 0 ? x : 0; };
    constexpr auto OutputActivationFnc = []( float x ) {
      // sigmoid
      return 1.f / ( 1.f + vdt::fast_expf( -x ) );
    };
    // build network structure
    // weight matrix from layer 0 to 1
    constexpr auto fWeightMatrix0to1 = std::array<std::array<float, 11>, 15>{
        {{0.0564278702527005, 0.0736385664395001, 0.111533698925958, -1.31428313514404, 0.0704978919437891,
          -0.156728845697468, 0.109500018264481, 0.00377299516653829, -24.0739682065448, -0.427773252339787,
          -23.528906669402},
         {-0.41378742994724, 0.501039742840929, 1.0969194785462, 0.216930356862149, 1.54269908083731, -2.29106477648207,
          -0.630905421958254, -1.54467238810462, 0.606585139519856, 1.21349406419577, -1.70009776141255},
         {-0.475703414108781, 1.35622309884546, 0.886233642836136, 3.66287985114785, -3.60756070813725,
          -0.297026255360717, -0.999231588912256, -1.4899472359416, 8.02346708145781, -3.07789748742688,
          5.06299750780424},
         {0.0638605145835328, 0.180439170179029, 0.121169070538529, 0.599037209301457, -0.0134201468677977,
          0.0904094968116527, 0.00414618947077304, -0.300670285279403, 2.90089567591911, 5.18033555611658,
          3.14738506199541},
         {-2.02341099091201, 1.51883843868399, 0.214635916367808, 1.48492194941506, 1.81705243901495, 0.579924004502985,
          -0.747525742491653, 0.444289133603457, 0.796972487106634, 1.4611422586121, -0.623442252006784},
         {-1.06540506155234, 0.228027281928625, 0.629059537315245, 2.48737033049966, 0.373998234933337, 1.5955126824544,
          0.998947082346914, -2.02803366273659, 0.923113302860426, 1.83731117315979, -1.27054780808538},
         {1.95308987855022, -0.464087453795175, -3.76259443090128, 0.194147304149005, 0.96531696932321,
          1.67342425079388, -1.88620665146531, 1.55162442136255, -1.19505889597544, 0.193743228484704,
          -0.129343843235445},
         {0.296563697573794, 0.357475565173729, 1.10750586433397, -1.63096381036845, -1.0342548740511, 1.20456422715473,
          -2.05475364798695, 0.0412023199198812, 1.60396056370954, 1.45155496275357, -1.05007226997895},
         {-0.705406570398858, -0.375012561342689, -0.0679161446006385, -0.0963954653775027, 0.748845226144644,
          -1.55718464681044, -3.00814085217088, 1.64671233642228, 0.643030115278188, -0.595020189921331,
          -1.19957445502576},
         {0.0845256285038076, 1.93760771006849, 0.519354405297375, 2.77837559741746, -0.532563935449298,
          0.138213135231722, -0.55874805412034, -3.31949858887676, 4.93596713254006, -0.243541431404453,
          5.34072592545661},
         {-0.783348645745542, 1.8870223751757, -0.930808522286314, -1.41490164755772, -0.0439230637009331,
          0.294959139665001, 0.862483020169701, -0.51145783378522, 4.86773949762358, -0.10237931360217,
          4.3711435317074},
         {-1.4597697161777, 2.52115033863364, 0.129079758003422, -0.34112026769147, -0.27844674164863,
          0.640590260383645, -2.45584188972217, 0.754143613655884, 0.33070540645097, 1.8231084798916,
          -1.09014756573786},
         {0.266488181029207, -1.44557692855659, -0.458796414951242, -8.82227151206785, 1.56616269000572,
          0.079657192322938, 0.165486447972991, 1.47662467651742, 3.66242218200325, -0.124380940498492,
          1.75843962523415},
         {-0.0812294094337669, -0.66203818588675, -0.124180980552198, 3.98067843204028, 0.665033390209112,
          0.0461039681255906, 0.575446408467601, 0.804721302422991, -10.1376277475253, 0.428192458730899,
          -6.64965421150736},
         {-0.459471607715004, 1.10363051803478, 0.0591645659814928, 1.1310150728631, 2.01100218271877,
          -0.48362692736632, 0.190173201926637, -0.947094901236744, 1.1202334827207, -1.31536974846196,
          0.465214449432741}}};

    constexpr auto fWeightMatrix1to2 =
        std::array<float, 16>{{-2.25595011494829, 0.947123406322107, 0.718208541996511, -2.64340263708723,
                               -0.832819423197834, 2.44743166402018, -0.556999795728784, 0.618371520673036,
                               0.442753277419524, 0.718535077332206, -1.59623813389326, 2.10328463027449,
                               -0.85021866403389, -0.930993801262542, -1.39535195016149, 1.11889113591803}};

    constexpr auto fMin =
        std::array<float, 10>{{9, 0.0107197267935, 2, 2, 0, 212, 218, 0.579863786697, 1.11619496346, 1.52495276928}};

    constexpr auto fMax =
        std::array<float, 10>{{12, 29.6306381226, 7, 8, 2, 8281, 5250, 38.7384567261, 11829.7861328, 8.42633724213}};

    // Normalization transformation
    constexpr auto transformer = TMV::Utils::Transformer{fMin, fMax};

    // the training input variables
    constexpr auto validator = TMV::Utils::Validator{
        "ReadGhostProbabilityDownstream",
        std::tuple{"UpgradeGhostInfo_obsFT", "UpgradeGhostInfo_FitTChi2", "UpgradeGhostInfo_FitTNDoF",
                   "UpgradeGhostInfo_obsUT", "UpgradeGhostInfo_UToutlier", "UpgradeGhostInfo_veloHits",
                   "UpgradeGhostInfo_utHits", "TRACK_CHI2", "TRACK_PT", "TRACK_ETA"}};

    constexpr auto l0To1 = TMV::Utils::Layer{fWeightMatrix0to1, ActivationFnc};
    constexpr auto l1To2 = TMV::Utils::Layer{fWeightMatrix1to2, OutputActivationFnc};
    constexpr auto MVA   = TMV::Utils::MVA{validator, transformer, 0, l0To1, l1To2};
  } // namespace
} // namespace Data::ReadGhostProbabilityDownstream

//_______________________________________________________________________

struct ReadGhostProbabilityDownstream final {

  // constructor
  ReadGhostProbabilityDownstream( LHCb::span<const std::string_view, 10> theInputVars ) {
    Data::ReadGhostProbabilityDownstream::MVA.validate( theInputVars );
  }

  // the classifier response
  // "inputValues" is a vector of input values in the same order as the
  // variables given to the constructor
  static constexpr auto GetMvaValue( LHCb::span<const float, 10> input ) {
    return Data::ReadGhostProbabilityDownstream::MVA( input );
  }
};
