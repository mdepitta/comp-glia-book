/*----------------------------------------------------------------------------
  AstroSim: Simulation of astrocyte networks Ca2+ dynamics
  Copyright (c) 2016-2017 Jules Lallouette
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#include <stdio.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_errno.h>

#include "ResultSaver.h"
#include "AbstractFactory.h"
#include "ChISimulationManager.h"
#include "StimulationStrat.h"
#include "ChIModel.h"
#include "GridSearchSimulation.h"
#include "PropagationModels.h"
#include "ErrorCodes.h"
#include "Synapse.h"
#include "Neuron.h"

using namespace AstroModel;
using namespace std;

void InitializeFactories();

int main(int argc, char * argv[])
{
	InitializeFactories();
	gsl_set_error_handler_off();

	int returnVal = 0;

	///////////////////////////////////////////////////
	// Parameters                                    //
	///////////////////////////////////////////////////

	ParamHandler & handler = ParamHandler::GlobalParams;

	bool debug,  save,  load,  saveModel,  loadModel,  useSubDir,
		saveSim, loadSim, splitGrid, simulate, useParamsFromFile,
		notSpacialNetwork,      normLinkStrength,     showNbSims,
		onlySpikePoiss,   useModTransTimeSERS,  useToroidalSpace,
		caSpontRelPoiss,                      functTopoUseThresh, 
		functTopoDirectedLinks,                functTopoSpecMetr,
		cnstrFrmSimUseSameMetr,  voroUseGhostPoints,  makeDirNet,
		frmFileDirected,  fourierTransConcatSig,  fourTrUseSTFFT,
		netDimFullNet,     functTopoCompL,      functTopoCompHCC, 
		functTopoUseMinForBidir,  threshCaSpontRel,  preRunToEqu,
		frmFileUseSparse, frmFileUseValAsStrengths, propOnlyOneWave,
		waveDetectUPOPANS, poissIsolateNode, poissUseGluStim,
		starLikeMakeSquare;
	double tStart,    tEnd,    Step,   F,   IPBias,   savingStep, 
		poissPeriod,   poissLength,  propMinDelay,  propMaxDelay, 
		randStimLength,	     randPauseLength,     propInstWindow, 
		maxLinkDist,        linkRadius,        erdosRenyiMeanDeg, 
		spatialScaleFreeRc,      swRewireProb,      thresholdMod, 
		netDimRepFract,     modThreshModel,    mplMeanPathLength, 
		mpltStart,   mpltEnd,   ModThreshBval,   SersA,    SersB,
		SersPeriod,   SersActTime,   SersTransTime,  SERSStimRad, 
		correlMaxLag,   correlPauseTime,    transEntropTimeEmbed,
		threshDetCaThresh,      cellDistIncr,      activCaThresh, 
		commPrecision,        scaleFreeGamma,       aStructBuild, 
		varStructBuild,  minDistStructBuild,   randMDStructBuild,
		minModTransTimeFactSERS,          simpleThreshSERSThresh, 
		simpleThreshSERSSpontFire,    simpleThreshSERSRecovProba,
		freqEstTimeWin,                          functTopoThresh, 
		functTopoCommonMeanDegThresh,       functTopoStdDevCoeff,
		functTopoCommonStdDevCoeff,    voroParamMult,    gluStim, 
		randIsoCellsRatio,fourTrSTFFTWinSize, fourTrSTFFTWinStep,
		fourTrDomFreqThrRat,  modVoroMinDecile,  fluxComputDelay,
		preRunTime, minFreqToSave, maxFreqToSave, stepFreqToSave,
		wavltMinFreqToSave, wavltMaxFreqToSave, wavltStepFreqToSave, 
		wavltDromFreqThreRat,               wavltPowrThreshFract, 
		defStimCouplStrength, slbInterCellDist, slbInterCompDist,
		voroMaxLinkDist, KStimVal, KStimIncrRate, ThreshDetCouplStr,
		poissGluQuantalRel, poissOmegaC, somaCouplStr;
	int N, Nastr,   Nneur,   desiredNb,    seed,    swNeighbDist, 
		thresholdStimRad,   functTopoCommonMethod,  maxRadToSave, 
		DefaultCouplingMethod, shellScrambleStartNode;
	unsigned int dim,      regularDegree,     spatialScaleFreeNl, 
		hccStart,     hccEnd,     hccRepeat,     propModSimStart, 
		propModSimEnd,   repeatSim,   nbSinks,   threshDetDegree, 
		mplNbStims,          mpltExp,          transEntropNbBins, 
		threshDetnbStimulated,  threshDetSinks,  slbBranchLength,
		slbNbBraches, slbNbEndBall, somaRadius;
	string modelLoadingPath,     resultFileName,      subDirPath, 
		simLoadingPath,      simSavingPath,      modelSavingPath, 
		savingPath,         loadingPath,         defaultGridName, 
		gridScriptFileName,       mainPath,       paramsFilePath, 
		constructType,     solverClassName,    couplingClassName,
		simTypeName,      poissCouplName,       defStimCouplName, 
		randCouplName,       threshCouplName,       mplCouplName, 
		correlCouplName,  structBuildType,  cnstrFrmSimModelName,
		neuronClassName,   synapseClassName,  neuronNetClassName, 
		astroNetClassName,      frmFilePath,      FrmFileNeurNet, 
		FrmFileAstrNet,    loadFrmFileStruct,     frmFileAstrPos, 
		frmFileNeurPos,                    randIsoCellsStratName, 
		frmAstrToSynFilePath,   neurConstrType,   astrConstrType,
		dummyNeurSpTrainPath,                    frmFileListPath,
		stochResSigStimClass, stochResNoiseStimClass, shellScrambleSubSTrat;
	vector<string> dataToSave,     gridParamNames,    paramNames, 
		paramVals,  metricNames,  stimTypes,  gridSingleParNames, 
		gridSingleParValues,                condOnSingleParNames,
		condOnSingleParValues,               condSingleParamName, 
		condSingleParamValue,   condOnParNames,  condOnParValues, 
		condParamName,  functTopoMetrList,  cnstrFrmSimParamName,
		cnstrFrmSimParamVals,   fourierTrNames,   fourierTrElems,
		functTopoByNetStratNames,   wavltTrNames,   wavltTrElems, 
		gridCondReplPath, stimStratNames;
	vector<int> stimIndices,      stimExpand,      modThreshStim,
		SERSmodStim, cnstrFrmSimParamPos, netDimNodeList;
	vector<double> stimStart,      stimEnd,     paramStartValues, 
		paramEndValues,     paramStepValues,     condParStartVal,
		condParEndVal,       condParStepVal,      propStartQuant, 
		propEndQuant, propStepQuant, propStartRatio, propEndRatio, 
		propStepRatio;

	handler <= "-useParams", useParamsFromFile = false, paramsFilePath = "";

	handler <= "-Save", save = false, savingPath = "GridSearch.dat";
	handler <= "-Load", load = false, loadingPath = "GridSearch.dat";
	handler <= "-modelSave", saveModel = false, modelSavingPath = "ModelSave.dat";
	handler <= "-modelLoad", loadModel = false, modelLoadingPath = "ModelSave.dat";
	handler <= "-simSave", saveSim = false, simSavingPath = "Simulation.dat";
	handler <= "-simLoad", loadSim = false, simLoadingPath = "Simulation.dat";

	handler <= "-SolverClass", solverClassName = "ODERungeKuttaSolverDouble";
	handler.AddAllowedValsList("-SolverClass", 0, 
		AbstractFactory<ODE::ODESolver<double, double> >::GetFactoriesNames());

	handler <= "-Sim", simulate = false;
	handler <= "-showNbSims", showNbSims = false;
	handler <= "-SplitGrid", splitGrid = false, desiredNb = 20, defaultGridName = "GridSearch",
		gridScriptFileName = "LaunchGrids.job";
	handler <= "-DgridRange", gridParamNames, paramStartValues, paramEndValues, paramStepValues;
	handler <= "-DgridPar", gridSingleParNames, gridSingleParValues;
	handler <= "-DcondGridPar", condOnSingleParNames, condOnSingleParValues, condSingleParamName, condSingleParamValue;
	handler <= "-DcondGridRange", condOnParNames, condOnParValues, condParamName, condParStartVal, condParEndVal, condParStepVal;

	handler <= "-DParVal", paramNames, paramVals;

	handler <= "-PropagModelSim", propModSimStart = 0, propModSimEnd = 100;
	handler <= "-ThresholdMod", thresholdMod = 0.1, thresholdStimRad = 1;
	handler <= "-ModThreshModel", modThreshModel = 1.0, ModThreshBval = 0;
	handler <= "-ModThreshStim", modThreshStim;
	handler <= "-SERSModel", SersPeriod = 9, SersActTime = 5, SersTransTime = 3.5, SERSStimRad = 1;
	handler <= "-SERSModelStims", SERSmodStim;
	handler <= "-SecOrdSERSModel", SersA = 0.06, SersB = 0.4;
	handler <= "-SecOrdSERSModelModTransT", useModTransTimeSERS = false, minModTransTimeFactSERS = 0.1;
	handler <= "-SimpleThreshSERSModel", simpleThreshSERSThresh = 0.2, simpleThreshSERSSpontFire = 0.0, simpleThreshSERSRecovProba = 0.25;

	handler <= "-debug", debug = false;
	handler <= "-tStart", tStart = 0;
	handler <= "-tEnd", tEnd = 1;
	handler <= "-Step", Step = 0.01;
	handler <= "-T", tStart, tEnd, Step;
	handler <= "-Nastro", Nastr = 1;
	handler <= "-Nneur", Nneur = 1;
	handler <= "-N", N = 1;
	handler <= "-F", F = 0.001, DefaultCouplingMethod = 0;
	handler <= "-IPBias", IPBias = 0.001;
	handler <= "-XtraCellGluStim", gluStim = 0.2;
	handler <= "-XtraCellKStim", KStimVal = 3.0, KStimIncrRate = 1.0e-14;
	handler <= "-seed", seed = time(0);
	handler <= "-PreRunTimeToEq", preRunToEqu = false, preRunTime = 20;

	handler <= "-SaveResults", resultFileName = "AstroRes";
	handler <= "-SavingStep", savingStep = 0.1;
	handler <= "-SubDir", useSubDir = false, subDirPath;
	handler <= "-GridCondReplPath", gridCondReplPath;
	handler <= "-Path", mainPath = ".";
	handler <= "-repeat", repeatSim = 1;

	handler <= "-CorrelParams", correlMaxLag = 10;
	handler <= "-TransferEntropy", transEntropTimeEmbed = 1.0, transEntropNbBins = 2;
	handler <= "-FunctTopoByNetConstr", functTopoByNetStratNames;
		handler.AddAllowedValsList("-FunctTopoByNetConstr", 0, AbstractFactory<NetworkConstructStrat>::GetFactoriesNames());
	handler <= "-FTopoSpecifMetrs", functTopoSpecMetr = false, functTopoMetrList;
		handler.AddAllowedValsList("-FTopoSpecifMetrs", 1, AbstractFactory<Metric>::GetFactoriesNames());
	handler <= "-FTopoCommon", functTopoCommonMethod = 1, functTopoCommonMeanDegThresh = 6.0, functTopoCommonStdDevCoeff = 1.5, modVoroMinDecile = 0.5;
	handler <= "-FTopoComputeL", functTopoCompL = false;
	handler <= "-FTopoComputeHCC", functTopoCompHCC = false;
	handler <= "-FTopoUseMinForBidir", functTopoUseMinForBidir = false;
	handler <= "-FourierTransform", fourierTrNames, fourierTrElems, minFreqToSave = 0.025, maxFreqToSave = 0.3, stepFreqToSave = 0.025;
	handler <= "-FourierTransConcatSigs", fourierTransConcatSig = false;
	handler <= "-FourierTransUseShortTimeFFT", fourTrUseSTFFT = false, fourTrSTFFTWinSize = 50, fourTrSTFFTWinStep = 1, fourTrDomFreqThrRat = 0.05;
	handler <= "-WaveletTransformElems", wavltTrNames, wavltTrElems;
	handler <= "-WaveletTransform", wavltMinFreqToSave = 0.025, wavltMaxFreqToSave = 0.3, wavltStepFreqToSave = 0.025, wavltDromFreqThreRat = 0.05, wavltPowrThreshFract = 0.6;

	handler <= "-commDetect", commPrecision = 0.000001;

	handler <= "-StimType", stimTypes; // Default Poisson RandomMono ThreshDetermStim MeanPathStim CorrelDeterm
		handler.AddAllowedValsList("-StimType", 0, AbstractFactory<StimulationStrat>::GetFactoriesNames());

	handler <= "-defStim", stimStratNames, stimIndices, stimStart, stimEnd, stimExpand;
		handler.AddAllowedValsList("-defStim", 0, AbstractFactory<StimulationStrat>::GetFactoriesNames());
	handler <= "-defStimCoupl", defStimCouplName = SigmoidCoupling::ClassName, defStimCouplStrength = 0.002;
	handler <= "-poissStim", poissPeriod = 200, poissLength = 10, 
		poissCouplName = SigmoidCoupling::ClassName;
	handler <= "-poissOnlySpike", onlySpikePoiss = false;
	handler <= "-poissCaSpontRel", caSpontRelPoiss = false;
	handler <= "-poissIsolateNode", poissIsolateNode = false;
	handler <= "-poissGluStim", poissUseGluStim = false, poissGluQuantalRel = 2.4, poissOmegaC = 60;
	handler <= "-randStim", randStimLength = 50, randPauseLength = 20, 
		randCouplName = SigmoidCoupling::ClassName;
	handler <= "-ThreshDeterm", nbSinks = 1, threshDetnbStimulated = 1, 
		threshDetCaThresh = 0.0007, threshCouplName = SigmoidCoupling::ClassName, ThreshDetCouplStr = 0.002;
	handler <= "-ThreshDetermCaSpontRel", threshCaSpontRel = false;
	handler <= "-MeanPathStimParams", mplNbStims = 3, mplMeanPathLength = 3, 
		mpltStart = tStart, mpltEnd = tEnd, mpltExp = 0, mplCouplName = SigmoidCoupling::ClassName;
	handler <= "-CorrelStim", correlPauseTime = correlMaxLag, correlCouplName = SigmoidCoupling::ClassName;
		handler.AddAllowedValsList("-defStimCoupl", 0, AbstractFactory<NetworkEdge>::GetFactoriesNames());
		handler.AddAllowedValsList("-poissStim", 2, AbstractFactory<NetworkEdge>::GetFactoriesNames());
		handler.AddAllowedValsList("-randStim", 2, AbstractFactory<NetworkEdge>::GetFactoriesNames());
		handler.AddAllowedValsList("-ThreshDeterm", 3, AbstractFactory<NetworkEdge>::GetFactoriesNames());
		handler.AddAllowedValsList("-MeanPathStimParams", 5, AbstractFactory<NetworkEdge>::GetFactoriesNames());
		handler.AddAllowedValsList("-CorrelStim", 1, AbstractFactory<NetworkEdge>::GetFactoriesNames());

	handler <= "-StructBuilder", structBuildType = StandardStructureBuilder::ClassName;
		handler.AddAllowedValsList("-StructBuilder", 0, AbstractFactory<SpatialStructureBuilder>::GetFactoriesNames());
	handler <= "-StdStructBuild", aStructBuild = 0.000070, varStructBuild = 0.000055, minDistStructBuild = 0.000005;
	handler <= "-RandStructBuild", randMDStructBuild = 0.000050;
	handler <= "-StarLikeBuild", slbInterCellDist = 0.00005, slbInterCompDist = 0.000001, slbBranchLength = 25, slbNbBraches = 5, slbNbEndBall = 1;
	handler <= "-StarLikeMakeSquare", starLikeMakeSquare = false;
	handler <= "-StarLikeNet", somaRadius = 3, somaCouplStr = 100;
	handler <= "-ShellScrambler", shellScrambleSubSTrat = SpatialRegularStrat::ClassName, shellScrambleStartNode = -1;
		handler.AddAllowedValsList("-ShellScrambler", 0, AbstractFactory<NetworkConstructStrat>::GetFactoriesNames());

	handler <= "-Construct", constructType = SpatialRegularStrat::ClassName; // Regular LinkRadius ErdosRenyi SpatialScaleFree SmallWorld ThreshDetermNet
		handler.AddAllowedValsList("-Construct", 0, AbstractFactory<NetworkConstructStrat>::GetFactoriesNames());
	handler <= "-AstrNetConstruct", astrConstrType = SpatialRegularStrat::ClassName;
		handler.AddAllowedValsList("-AstrNetConstruct", 0, AbstractFactory<NetworkConstructStrat>::GetFactoriesNames());
	handler <= "-NeurNetConstruct", neurConstrType = SpatialRegularStrat::ClassName;
		handler.AddAllowedValsList("-NeurNetConstruct", 0, AbstractFactory<NetworkConstructStrat>::GetFactoriesNames());
	handler <= "-scaleFree", scaleFreeGamma = 2.0;
	handler <= "-regularConstr", regularDegree = 6, maxLinkDist = 0.00015;
	handler <= "-linkRadiusConstr", linkRadius = 0.00015;
	handler <= "-erdosRenyiMeanDeg", erdosRenyiMeanDeg = 6;
	handler <= "-spatialScaleFree", spatialScaleFreeRc = 0.0001, spatialScaleFreeNl = 5;
	handler <= "-SWParams", swRewireProb = 0, swNeighbDist = 1;
	handler <= "-ThreshDetermNet", threshDetDegree = 10, threshDetSinks = 2;
	handler <= "-functTopoUseThresh", functTopoUseThresh = false, functTopoThresh = 0.1;
	handler <= "-functTopoDirectedLinks", functTopoDirectedLinks = false;
	handler <= "-functTopoAutoThresh", functTopoStdDevCoeff = 1.5;
	handler <= "-constrFromSimParams", cnstrFrmSimParamName, cnstrFrmSimParamPos, cnstrFrmSimParamVals;
	handler <= "-constrFromSimModelName", cnstrFrmSimModelName = "ChIModel";
		handler.AddAllowedValsList("-constrFromSimModelName", 0, AbstractFactory<ChIModel>::GetFactoriesNames());
	handler <= "-constrFromSimUseSameMetr", cnstrFrmSimUseSameMetr = false;
	handler <= "-voronoiParam", voroParamMult = 10000000, voroMaxLinkDist = 1.0;
	handler <= "-voronoiParamGhostPoints", voroUseGhostPoints = false;
	handler <= "-frmFilePath", frmFilePath = "";
	handler <= "-frmFileDirected", frmFileDirected = false;
	handler <= "-frmFileUseSparse", frmFileUseSparse = false;
	handler <= "-frmFileUseValAsStrengths", frmFileUseValAsStrengths = false;
	handler <= "-frmFileListPath", frmFileListPath = "";
	handler <= "-LoadFrmFileStructBuild", loadFrmFileStruct = "";
	handler <= "-isolatedCells", randIsoCellsRatio = 0.1, randIsoCellsStratName = VoronoiConstructStrat::ClassName;
		handler.AddAllowedValsList("-isolatedCells", 1, AbstractFactory<NetworkConstructStrat>::GetFactoriesNames());

	handler <= "-normLink", normLinkStrength = false;
	handler <= "-makeDirNet", makeDirNet = false;

	handler <= "-Coupling", couplingClassName = "SigmoidFunction";
	handler.AddAllowedValsList("-Coupling", 0, AbstractFactory<NetworkEdge>::GetFactoriesNames());

	handler <= "-notSpacial", notSpacialNetwork = false;

	handler <= "-WaveDetectUPOPANS", waveDetectUPOPANS = false;
	handler <= "-activParams", activCaThresh = 0.7e-3, freqEstTimeWin = 10, fluxComputDelay = 7;
	handler <= "-propDelays", propMinDelay = 5.0, propMaxDelay = 20.0, propInstWindow = 5.0, cellDistIncr = 1.0;
	handler <= "-propOnlyOneWave", propOnlyOneWave = false;
	handler <= "-propSaveQuant", DefaultVal(propStartQuant, 0.0), DefaultVal(propEndQuant, 1.0), DefaultVal(propStepQuant, 0.1);
	handler <= "-propSaveRatio", DefaultVal(propStartRatio, 0.1), DefaultVal(propEndRatio, 1.0), DefaultVal(propStepRatio, 0.1);
	handler <= "-dim", dim = 2;
	handler <= "-toroidalSpace", useToroidalSpace = false;
	handler <= "-sD", dataToSave;
		handler.AddAllowedVal("-sD", 0, "FullAdjacencyMatrix");
		handler.AddAllowedVal("-sD", 0, "SparseAdjacencyMatrix");
		handler.AddAllowedVal("-sD", 0, "Positions");
		handler.AddAllowedVal("-sD", 0, "DegreeDistr");
		handler.AddAllowedVal("-sD", 0, "StrengthsDistr");
		handler.AddAllowedVal("-sD", 0, "DegreeDistrStats");
		handler.AddAllowedVal("-sD", 0, "ClustCoeffDistr");
		handler.AddAllowedVal("-sD", 0, "HierarchClustCoeff");
		handler.AddAllowedVal("-sD", 0, "FullAllPairDist");
		handler.AddAllowedVal("-sD", 0, "PathLengthInfos");
		handler.AddAllowedVal("-sD", 0, "CellToCellDist");
		handler.AddAllowedVal("-sD", 0, "NetworkDimensions");
		handler.AddAllowedVal("-sD", 0, "NetDimFullData");
		handler.AddAllowedVal("-sD", 0, "NetDimFinalValues");
		handler.AddAllowedVal("-sD", 0, "NetworkCommunityStats");
		handler.AddAllowedVal("-sD", 0, "NetworkFullCommunities");
		//
		handler.AddAllowedVal("-sD", 0, "Stimulations");
		handler.AddAllowedVal("-sD", 0, "MeanPathStimInfos");
		//
		handler.AddAllowedVal("-sD", 0, "ActivatedCells");
		handler.AddAllowedVal("-sD", 0, "OutFluxesFull");
		handler.AddAllowedVal("-sD", 0, "FullWaves");
		handler.AddAllowedVal("-sD", 0, "WaveMetrics");
		handler.AddAllowedVal("-sD", 0, "Concentrations");
		handler.AddAllowedVal("-sD", 0, "WaveSummary");
		handler.AddAllowedVal("-sD", 0, "SpikingCoherence");
		handler.AddAllowedVal("-sD", 0, "ThresholdDetermination");
		handler.AddAllowedVal("-sD", 0, "ThresholdDeterminationFull");
		handler.AddAllowedVal("-sD", 0, "FullFrontierCells");
		handler.AddAllowedVal("-sD", 0, "FrontierCellsStats");
		//
		handler.AddAllowedVal("-sD", 0, "PropagActivatedCells");
		handler.AddAllowedVal("-sD", 0, "DetailedActivatedCells");
		//
		handler.AddAllowedVal("-sD", 0, "SpatialDistByNbCell");
		handler.AddAllowedVal("-sD", 0, "MaxNbCellDistribution");
		handler.AddAllowedVal("-sD", 0, "SpatialDistByCellDist");
		handler.AddAllowedVal("-sD", 0, "MaxWaveSpeedDistribution");
		handler.AddAllowedVal("-sD", 0, "NbCellInWaveByClustCoeff");
		handler.AddAllowedVal("-sD", 0, "NbWavesDistribution");
		//
		handler.AddAllowedVal("-sD", 0, "CellToCellDistStats");
		handler.AddAllowedVal("-sD", 0, "ScalarStatistics");
		//
		handler.AddAllowedVal("-sD", 0, "GridFilePathIndexes");
		handler.AddAllowedVal("-sD", 0, "GridPaths");
		handler.AddAllowedVal("-sD", 0, "GridScalarStatistics");
		handler.AddAllowedVal("-sD", 0, "CondParamReplace");
		//
		handler.AddAllowedVal("-sD", 0, "ZeroLagCorrelationMatrix");
		handler.AddAllowedVal("-sD", 0, "MaxCorrelationValuesMatrix");
		handler.AddAllowedVal("-sD", 0, "MaxCorrelationLagsMatrix");
		handler.AddAllowedVal("-sD", 0, "CommonNeighbsZeroCorrel");
		handler.AddAllowedVal("-sD", 0, "CorrelDirectionalityMatrix");
		handler.AddAllowedVal("-sD", 0, "MaxCorrelationValuesNoZeroMatrix");
		handler.AddAllowedVal("-sD", 0, "MaxCorrelationNoZeroLagsMatrix");
		handler.AddAllowedVal("-sD", 0, CorrelationsMetric::ZeroCorrTopoName + "_ROCData");
		handler.AddAllowedVal("-sD", 0, CorrelationsMetric::MaxCorrTopoName + "_ROCData");
		handler.AddAllowedVal("-sD", 0, CorrelationsMetric::MaxCorrNoZeroTopoName + "_ROCData");
		handler.AddAllowedVal("-sD", 0, CorrelationsMetric::ZeroCorrTopoName + "_ConnecByThresh");
		handler.AddAllowedVal("-sD", 0, CorrelationsMetric::MaxCorrTopoName + "_ConnecByThresh");
		handler.AddAllowedVal("-sD", 0, CorrelationsMetric::MaxCorrNoZeroTopoName + "_ConnecByThresh");
		//
		handler.AddAllowedVal("-sD", 0, "TransferEntropyFullMat");
		handler.AddAllowedVal("-sD", 0, TransferEntropyMetric::TransEntrFunctTopoName + "_ROCData");
		handler.AddAllowedVal("-sD", 0, TransferEntropyMetric::TransEntrFunctTopoName + "_ConnecByThresh");
		//
		handler.AddAllowedVal("-sD", 0, "FullFourierTransform");
		handler.AddAllowedVal("-sD", 0, "FullFourierSpectrogram");
		handler.AddAllowedVal("-sD", 0, "FourierDominantFrequencies");
		//
		handler.AddAllowedVal("-sD", 0, "WaveletSpectroGrams");

	handler <= "-HCCParams", hccStart = 1, hccEnd = 3, hccRepeat = 20;
	handler <= "-NetDimParams", netDimRepFract = 1.0, maxRadToSave = 5;
	handler <= "-NetDimNodeList", netDimNodeList;
	handler <= "-NetDimFullNet", netDimFullNet = false;

	handler <= "-SimulationType", simTypeName = "RepeatSimulationChIModel";
	handler.AddAllowedValsList("-SimulationType", 0, AbstractFactory<SimulationManager>::GetFactoriesNames());

	// RepeatSimulationChIModel
	// RepeatSimulationThresholdModel
	// RepeatSimulationModifiedThresholdModel
	// RepeatSimulationSERSModel
	// RepeatSimulationChIModelThresholdDetermination

	handler <= "-NeuronClass", neuronClassName = SFALIFNeuron::ClassName;
	handler.AddAllowedValsList("-NeuronClass", 0, AbstractFactory<Neuron>::GetFactoriesNames());

	handler <= "-DummyNeuronSpikeTrain", dummyNeurSpTrainPath;

	handler <= "-SynapseClass", synapseClassName = TMSynapse::ClassName;
	handler.AddAllowedValsList("-SynapseClass", 0, AbstractFactory<Synapse>::GetFactoriesNames());

	handler <= "-NeuronNetClass", neuronNetClassName = NeuronNetModel::ClassName;
	handler.AddAllowedValsList("-NeuronNetClass", 0, AbstractFactory<NeuronNetModel>::GetFactoriesNames());

	handler <= "-AstroNetClass", astroNetClassName = ChIModel::ClassName;
	handler.AddAllowedValsList("-AstroNetClass", 0, AbstractFactory<ChIModel>::GetFactoriesNames());

	handler <= "-FrmFileNeurNet", FrmFileNeurNet;
	handler <= "-FrmFileAstrNet", FrmFileAstrNet;

	handler <= "-FrmFileAstrPos", frmFileAstrPos;
	handler <= "-FrmFileNeurPos", frmFileNeurPos;

	handler <= "-AstrToSynFile", frmAstrToSynFilePath;

	handler <= "-StochRes", stochResSigStimClass = DefaultStimStrat::ClassName, stochResNoiseStimClass = PoissonianStimStrat::ClassName;
	handler.AddAllowedValsList("-StochRes", 0, AbstractFactory<StimulationStrat>::GetFactoriesNames());
	handler.AddAllowedValsList("-StochRes", 1, AbstractFactory<StimulationStrat>::GetFactoriesNames());

	handler <= "-aM", metricNames;
	handler.AddAllowedValsList("-aM", 0, AbstractFactory<Metric>::GetFactoriesNames());

	// PropagationMetrics
		// BoolWrp_BoolLnk_CellStateSaverMetric
		// SERStates_BoolLnk_CellStateSaverMetric
	// ChIModelMetrics
		// ActivatedCellsMetric
		// ConcentrationsMetric
		// PropagationDistanceMetric
		// CorrelationsMetric
		// TransferEntropyMetric
		// ThresholdDeterminationMetric
		// WaveFrontDetectMetric
	// Network Metrics
		// DegreeDistributionNetworkMetric
		// ClusteringCoefficientDistributionNetworkMetric
		// AllPairDistancesNetworkMetric
		// AdjacencyMatrixNetworkMetric
		// SpatialPositionsNetworkMetric
		// DimensionsNetworkMetric
	// Simulation Metrics
		// PropagationDistanceStatMetric
		// NetworkTopologyAndSpatialStatMetric
		// GridSearchMetric
	// Stimulation Metrics
		// StimulatedCellsMetric
		// MeanPathStimulationMetric

	if (not handler.Parse(argc, argv))
		return PARAMETER_PARSING_FAILED;
		
	// If additional params are specified in external file, load the file
	if (useParamsFromFile)
	{
		if (not handler.LoadParams(paramsFilePath))
			return PARAMETER_LOADING_FROM_FILE_FAILED;
	}

	if (showNbSims)
		mainPath = ".";

	globalRNG.SetSeed(seed);

	stringstream path;
	path << mainPath << "/data";
	if (useSubDir)
		path << "/" << subDirPath;

	//Update saving paths
	if (savingPath[0] != '/')
		savingPath = path.str() + "/" + savingPath;
	if (loadingPath[0] != '/')
		loadingPath = path.str() + "/" + loadingPath;
	if (simSavingPath[0] != '/')
		simSavingPath = path.str() + "/" + simSavingPath;
	if (simLoadingPath[0] != '/')
		simLoadingPath = path.str() + "/" + simLoadingPath;
	if (modelSavingPath[0] != '/')
		modelSavingPath = path.str() + "/" + modelSavingPath;
	if (modelLoadingPath[0] != '/')
		modelLoadingPath = path.str() + "/" + modelLoadingPath;

	///////////////////////////////////////////////////
	// Structures loading / Creation                 //
	///////////////////////////////////////////////////

	GridSearchSimulation *simManag = 0;
	SimulationManager *subSim = 0;

	if (load)
	{
		TRACE_UP("Loading grid search simulation")

		ifstream stream(loadingPath.c_str());
		if (stream.good())
			simManag = new GridSearchSimulation(stream, handler);
		else
		{
			cerr << "Couldn't open grid search file : " << loadingPath << endl;
			return GRIDSEARCH_LOADING_FAILED;
		}
		TRACE_DOWN("Grid search simulation loaded")
	}
	else
	{
		TRACE_UP("Creating repeat simulation and grid search simulation")

		assert(AbstractFactory<SimulationManager>::Factories[simTypeName]);
		subSim = AbstractFactory<SimulationManager>::Factories[simTypeName]->Create();

		simManag = new GridSearchSimulation(subSim, false, handler, path.str());

		// Add gridsearch parameters by single values
		for (unsigned int i = 0 ; i < gridSingleParNames.size() ; ++i)
			simManag->AddGridParameterValue(gridSingleParNames[i], gridSingleParValues[i]);
		// Add gridsearch parameters and their ranges
		for (unsigned int i = 0 ; i < gridParamNames.size() ; ++i)
			simManag->AddGridParameterRange(gridParamNames[i], paramStartValues[i], 
				paramEndValues[i], paramStepValues[i]);
		// gridSize
		// gridJitter
		// gridMinDist

		// IP3Bias
		// DefaultCouplingStrength
		// DefaultCouplingStrengthStdDev
		// DefaultLinkThreshold
		// DefaultLinkThresholdStdDev
		// DefaultLinkScale
		// DefaultLinkScaleStdDev
		// InitCaVarRatio
		// InithVarRatio
		// InitIP3VarRatio
		// DefaultCaVal
		// DefaulthVal
		// DefaultIP3Val
		// c1
		// rC
		// rL
		// vER
		// a2
		// Ker
		// vd
		// Kplcd
		// kd
		// k3

		// SpatialRegularDegree
		// SpatialConnectionRadius
		// ErdosRenyiMeanDegree
		// SpatialScaleFreeInteractionRange
		// SpatialScaleFreeNewLinks
		// SmallWorldProba
		// SmallWorldNeighbDist
		// PropagationModelStart
		// PropagationModelEnd
		// ThresholdModelThresh
		// ThresholdModelStimRadius
		// [functTopoName]_FunctionalTopoConstrThresh (?)
		// [functTopoName]_FunctionalTopoStdDevCoeff (?)
		// FunctTopoThresholdEstimMethod
		// FunctTopoThresholdEstimMeanDeg
		// FunctTopoThresholdEstimStdDevCoeff
		// RandIsoCellsRatio
		// RandIsoCellsSubStratName

		// ActivationThreshold
		// FrequencyEstimTimeWin

		// DefaultStimInd_xxxxx
		// DefaultStimStart_xxxxx
		// DefaultStimEnd_xxxxx
		// DefaultStimRadius_xxxxx

		// ThresholdDeterminationDegree
		// ThresholdDeterminationDerivs
		// ThreshDetermNbStim
		// ThreshDetermNbSinks
		// ThreshDetermMaxCa
		// MeanPathStimNbStims
		// MeanPathStimMeanPathLength
		// MeanPathStimTStart
		// MeanPathStimTEnd
		// MeanPathStimExp
		// ModThreshPropConst
		// ModThreshBVal
		// SERSModelPeriod
		// SERSModelActTime
		// SERSModelTransTime
		// SERSModelStimRadius
		// NetworkConstructStrat
		// NetworkSize
		// C0
		// d1 
		// d2 
		// d3 
		// d5 
		// v3k
		// K3k
		// r5p
		// ChISimStartTime
		// ChISimEndTime
		// PoissonianCouplStrength
		// PoissonianIP3Bias
		// SimpleThreshSERSModelRelatThresh
		// SecOrdSERSModelA
		// SecOrdSERSModelB
		// SecOrdSERSModelModTransT
		// SecOrdSERSModelModTransTFact
		// SimpleThreshSERSModelRelatThresh
		// SimpleThreshSERSModelSpontFireRate
		// SimpleThreshSERSModelRecovProba
		// NetworkEdgeClassName

		// Add conditional grid parameter values
		for (unsigned int i = 0 ; i < condOnSingleParNames.size() ; ++i)
			simManag->AddCondGridParameterValue(condOnSingleParNames[i], 
				condOnSingleParValues[i], condSingleParamName[i], condSingleParamValue[i]);
		// Add conditional grid parameter ranges
		for (unsigned int i = 0 ; i < condOnParNames.size() ; ++i)
			simManag->AddCondGridParameterRange(condOnParNames[i], condOnParValues[i], 
				condParamName[i], condParStartVal[i], condParEndVal[i], condParStepVal[i]);

		// Add global parameters values
		assert(paramNames.size() == paramVals.size());
		for (unsigned int i = 0 ; i < paramNames.size() ; ++i)
			simManag->AddGlobalParameterValue(paramNames[i], paramVals[i]);

		// Add all metrics
		for (unsigned int i = 0 ; i < metricNames.size() ; ++i)
		{
			assert(AbstractFactory<Metric>::Factories[metricNames[i]]);
			if (not simManag->AddMetric(AbstractFactory<Metric>::Factories[metricNames[i]]->Create(), true))
				std::cerr << "Problem detected while adding " << metricNames[i] << " to " << simManag->GetClassName() << std::endl;
		}

		TRACE_DOWN("Repeat simulation and grid search simulation have been created")
	}

	///////////////////////////////////////////////////
	// Results saving                                //
	///////////////////////////////////////////////////
	if (showNbSims)
		std::cout 
			<< "Total number of simulations : " 
			<< simManag->GetFullGridSize() << std::endl;

	///////////////////////////////////////////////////
	// Results saving                                //
	///////////////////////////////////////////////////

	// Register for saving the files specified by the -SD parameter.
	// If simManag was loaded from file, it allows the addition of
	// -sD files and works with doublons.
	assert(simManag);
	for (unsigned int i = 0 ; i < dataToSave.size() ; ++i)
		simManag->AddFileToSave(dataToSave[i]);

	///////////////////////////////////////////////////
	// New structures creation                       //
	///////////////////////////////////////////////////

	const string defaultSimName = "Simulation";
	if (splitGrid)
	{
		if (not simManag->ComputeAndSaveGridSearches(desiredNb, 
			path.str() + "/" + gridScriptFileName, 
			path.str() + "/" + defaultGridName, 
			path.str() + "/" + defaultSimName + "_g"))
		{
			std::cerr << "Failed to save smaller grid searches : " << gridScriptFileName << std::endl;
			return SUBGRID_SAVING_FAILED;
		}
	}

	///////////////////////////////////////////////////
	// Structures Saving                             //
	///////////////////////////////////////////////////

	if (save)
	{
		ofstream stream(savingPath.c_str());
		if (not simManag->SaveToStream(stream))
		{
			cerr << "Failed to save grid search to file : " << savingPath << endl;
			return GRIDSEARCH_SAVING_FAILED;
		}
	}
	if (saveSim)
	{
		ofstream stream(simSavingPath.c_str());
		if (not subSim->SaveToStream(stream))
		{
			cerr << "Failed to save sub simulation to file : " << simSavingPath << endl;
			return SIMULATION_SAVING_FAILED;
		}
	}

	///////////////////////////////////////////////////
	// Problem solving                               //
	///////////////////////////////////////////////////

	if (simulate)
	{
		TRACE_UP("Launching Simulation.")
		returnVal |= simManag->LaunchSimulation(string("Sim_") + StringifyFixed(seed));
		TRACE_DOWN("Simulation Ended.")
	}

	///////////////////////////////////////////////////
	// Model saving                                  //
	///////////////////////////////////////////////////

	if (saveModel)
	{
		ofstream stream(modelSavingPath.c_str());
		if (load)
			subSim = simManag->GetSubSimulation();

		AbstractRepeatSimulation *subSimTemp;
		if ((subSimTemp = dynamic_cast<AbstractRepeatSimulation *>(subSim)))
		{
			if (not subSimTemp->SaveModelToStream(stream))
			{
				cerr << "Failed to save model to file : " << modelSavingPath << endl;
				return MODEL_LOADING_FAILED;
			}
		}

		if (load)
			subSim = 0;
	}

TRACE("Delete simManag")
	if (simManag)
		delete simManag;

TRACE("Delete subSim")
	if (subSim)
		delete subSim;

	return returnVal;
}
