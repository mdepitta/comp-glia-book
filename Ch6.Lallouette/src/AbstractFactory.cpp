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

#include "AbstractFactory.h"
#include "CouplingFunction.h"
#include "StimulationStrat.h"
#include "Network.h"
#include "SpatialNetwork.h"
#include "SimulationMetrics.h"
#include "SimulationManager.h"
#include "ChISimulationManager.h"
#include "GridSearchSimulation.h"
#include "PropagationModels.h"
#include "PropagationMetrics.h"
#include "Neuron.h"
#include "Synapse.h"
#include "AstroNeuroModel.h"
#include "ChIModel.h"

using namespace std;
using namespace AstroModel;

//********************************************************************//
//****************** A B S R A C T   F A C T O R Y *******************//
//********************************************************************//

template <> string Network<CouplingFunction>::ClassName("StandardCouplingFunctionNetwork");
template <> string Network<BooleanLink>::ClassName("StandardBooleanNetwork");
template <> string Network<Synapse>::ClassName("StandardSynapticNetwork");
template <> string SpatialNetwork<CouplingFunction>::ClassName("SpatialCouplingFunctionNetwork");
template <> string SpatialNetwork<BooleanLink>::ClassName("SpatialBooleanNetwork");
template <> string SpatialNetwork<Synapse>::ClassName("SpatialSynapticNetwork");

template <> string ODE::RungeKuttaSolver<double, double>::ClassName("ODERungeKuttaSolverDouble");
template <> string ODE::EulerSolver<double, double>::ClassName("ODEEulerSolverDouble");

// Network Links
template <> map<string, AbstractFactory<NetworkEdge>*> 
	AbstractFactory<NetworkEdge>::Factories = map<string, AbstractFactory<NetworkEdge>*>();
static DerivedFactory<NetworkEdge, SigmoidCoupling> SigmoFuncFact;
static DerivedFactory<NetworkEdge, LinearCoupling>  LinFuncFact;
static DerivedFactory<NetworkEdge, BooleanLink> BooleanLinkFact;
static DerivedFactory<NetworkEdge, TMSynapse> TMSynEdgeFact;
static DerivedFactory<NetworkEdge, TMSynapseOptim> TMSynOptEdgeFact;

// Neuron models
template <> map<string, AbstractFactory<Neuron>*> 
	AbstractFactory<Neuron>::Factories = map<string, AbstractFactory<Neuron>*>();
static DerivedFactory<Neuron, SFALIFNeuron> SFALIFNeuronFact;
static DerivedFactory<Neuron, DummyNeuron> DummyNeuronFact;

// Synapse models
template <> map<string, AbstractFactory<Synapse>*> 
	AbstractFactory<Synapse>::Factories = map<string, AbstractFactory<Synapse>*>();
static DerivedFactory<Synapse, TMSynapse> TMSynapseFact;
static DerivedFactory<Synapse, TMSynapseOptim> TMSynapseOptimFact;

// Stimulation Strategies
template <> map<string, AbstractFactory<StimulationStrat>*> 
	AbstractFactory<StimulationStrat>::Factories = map<string, AbstractFactory<StimulationStrat>*>();
static DerivedFactory<StimulationStrat, DefaultStimStrat> DefStimStratFact;
static DerivedFactory<StimulationStrat, PoissonianStimStrat> PoissStimStratFact;
static DerivedFactory<StimulationStrat, RandomMonoStim> RandomMonoStimStratFact;
static DerivedFactory<StimulationStrat, ThresholdDeterminationStimStrat> ThresholdDeterminationStimStratFact;
static DerivedFactory<StimulationStrat, MeanPathStimStrat> MeanPathStimStratFact;
static DerivedFactory<StimulationStrat, CorrelDetermStimStrat> CorrelDetermStimStratFact;
static DerivedFactory<StimulationStrat, ExtraCellGluStim> ExtraCellGluStimFact;

// Networks
template <> map<string, AbstractFactory<AbstractNetwork >*> 
	AbstractFactory<AbstractNetwork >::Factories = map<string, AbstractFactory<AbstractNetwork>*>();
static DerivedFactory<AbstractNetwork, Network<CouplingFunction> > StandardNetFact;
static DerivedFactory<AbstractNetwork, Network<BooleanLink> > BooleanStandardNetFact;
static DerivedFactory<AbstractNetwork, SpatialNetwork<CouplingFunction> > SpatialNetFact;
static DerivedFactory<AbstractNetwork, SpatialNetwork<BooleanLink> > BooleanSpatialNetFact;
static DerivedFactory<AbstractNetwork, SpatialNetwork<Synapse> > SpatialSynaptNetFact;

// Network building strategies
template <> map<string, AbstractFactory<NetworkConstructStrat >*> 
	AbstractFactory<NetworkConstructStrat >::Factories = map<string, AbstractFactory<NetworkConstructStrat >*>();
static DerivedFactory<NetworkConstructStrat, ScaleFreeStrat> ScaleFreeConstructStrat;
static DerivedFactory<NetworkConstructStrat, SpatialScaleFreeStrat> SpatialScaleFreeConstructStrat;
static DerivedFactory<NetworkConstructStrat, SpatialRegularStrat> SpatialRegularConstructStrat;
static DerivedFactory<NetworkConstructStrat, SpatialConnectionRadiusStrat> SpatialConnectionRadiusConstructStrat;
static DerivedFactory<NetworkConstructStrat, VoronoiConstructStrat> VoronoiConstructStratFact;
static DerivedFactory<NetworkConstructStrat, ErdosRenyiRandomStrat> ErdosRenyiRandomConstructStrat;
static DerivedFactory<NetworkConstructStrat, SmallWorldStrat> SmallWorldConstructStrat;
static DerivedFactory<NetworkConstructStrat, ThresholdDeterminationStrat> ThresholdDeterminationConstructStrat;
static DerivedFactory<NetworkConstructStrat, FunctionalTopoStrat> FunctionalTopoConstructStrat;
static DerivedFactory<NetworkConstructStrat, ConstrFromChISim> ConstrFromChISimStratFact;
static DerivedFactory<NetworkConstructStrat, FromFileConstrStrat> ConstrFromFiletratFact;
static DerivedFactory<NetworkConstructStrat, FileListConstrStrat> ConstrFromFileListStratFact;
static DerivedFactory<NetworkConstructStrat, RandomIsolatedCells> RandIsoCellsStrat;
static DerivedFactory<NetworkConstructStrat, StarLikeAstroNetConstrStrat> StartLikeAstrStratFact;
static DerivedFactory<NetworkConstructStrat, ShellStructScramblerStrat> ShellStructScramblerStratFact;

// Spatial structures building strategies
template <> map<string, AbstractFactory<SpatialStructureBuilder>*> 
	AbstractFactory<SpatialStructureBuilder>::Factories = map<string, AbstractFactory<SpatialStructureBuilder>*>();
static DerivedFactory<SpatialStructureBuilder, StandardStructureBuilder> StandardSructureStrat;
static DerivedFactory<SpatialStructureBuilder, RandomStructureBuilder> RandomSructureStrat;
static DerivedFactory<SpatialStructureBuilder, LoadFromFileStructureBuilder> LoadFromFileSructureStrat;
static DerivedFactory<SpatialStructureBuilder, StarLikeAstrocyteStruct> StarLikeSructureStrat;

// Metrics
template <> map<string, AbstractFactory<Metric>*> 
	AbstractFactory<Metric>::Factories = map<string, AbstractFactory<Metric>*>();
static DerivedFactory<Metric, DegreeDistrComp> DegreeDistCompFact;
static DerivedFactory<Metric, ClusteringCoeffs> ClustCoeffDistCompFact;
static DerivedFactory<Metric, AllPairDistances> AllPairDistancesCompFact;
static DerivedFactory<Metric, AdjMatrixComp> AdjMatrixCompFact;
static DerivedFactory<Metric, PositionsComp> PositionsCompFact;
static DerivedFactory<Metric, StimulatedCells> StimulatedCellsFact;
static DerivedFactory<Metric, ActivatedCells> ActivatedCellsFact;
static DerivedFactory<Metric, ConcentrationMetrics> ConcentrationMetricsFact;
static DerivedFactory<Metric, PropagationDistance> PropagationDistanceFact;
static DerivedFactory<Metric, PropagationDistanceStat> PropagationDistanceStatFact;
static DerivedFactory<Metric, NetworkTopologyStats> NetworkTopologyStatsFact;
static DerivedFactory<Metric, NetworkDimensions> NetworkDimensionsFact;
static DerivedFactory<Metric, GridSearchMetric> GridSearchMetricFact;
static DerivedFactory<Metric, ThresholdDetermination> ThresholdDeterminationFact;
static DerivedFactory<Metric, MeanPathStimMetric> MeanPathStimMetricFact;
static DerivedFactory<Metric, CorrelationsMetric> CorrelationsMetricFact;
static DerivedFactory<Metric, TransferEntropyMetric> TransferEntropyMetricFact;
static DerivedFactory<Metric, FunctTopoByNetConstrStrat> FunctTopoByNetConstrMetricFact;
static DerivedFactory<Metric, WaveFrontDetect> WaveFrontDetectFact;
static DerivedFactory<Metric, FourierTransform> FourierTransFact;
static DerivedFactory<Metric, WaveletTransform> WaveletTransFact;
static DerivedFactory<Metric, ScalarStatisticsMetric> ScalarStatisticsMetricFact;
static DerivedFactory<Metric, StochResMetric> StochResMetrFact;

static DerivedFactory<Metric, CellStateSaver<BooleanWrapper, BooleanLink> > CellStateSaverBoolStBoolLnkFact;
static DerivedFactory<Metric, CellStateSaver<SERStates, BooleanLink> > CellStateSaverSERSStBoolLnkFact;

// Simulations
template <> map<string, AbstractFactory<SimulationManager>*> 
	AbstractFactory<SimulationManager>::Factories = map<string, AbstractFactory<SimulationManager>*>();
static DerivedFactory<SimulationManager, GridSearchSimulation> GridSearchSimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<ChIModel> > RepeatChISimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<ThresholdModel> > RepeatThreshSimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<ModifiedThresholdModel> > RepeatModifThreshSimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<SecondOrdNeighbSERSModel> > RepeatSecondOrdSERSSimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<SimpleThresholdSERSModel> > RepeatSimpleThreshSERSSimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<ChIModelThresholdDetermination> > RepeatChIThreshDetermSimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<NeuronNetModel> > RepeatNeuronNetSimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<AstroNeuroNetModel> > RepeatAstroNeuroSimFact;
static DerivedFactory<SimulationManager, RepeatSimulation<ChIModelStochasticRes> > RepeatChIStochResFact;

// ODE Solvers
template <> map<string, AbstractFactory<ODE::ODESolver<double, double> >* >
	AbstractFactory<ODE::ODESolver<double, double> >::Factories = map<string, AbstractFactory<ODE::ODESolver<double, double> >* >();
static DerivedFactory<ODE::ODESolver<double, double>, ODE::RungeKuttaSolver<double> > ODERungeKuttaSolverFact;
static DerivedFactory<ODE::ODESolver<double, double>, ODE::EulerSolver<double> > ODEEulerSolverFact;

// ChI Models
template <> map<string, AbstractFactory<ChIModel>*> 
	AbstractFactory<ChIModel>::Factories = map<string, AbstractFactory<ChIModel>*>();
static DerivedFactory<ChIModel, ChIModel > StdChIModFact;

// Neuron net models
template <> map<string, AbstractFactory<NeuronNetModel>*> 
	AbstractFactory<NeuronNetModel>::Factories = map<string, AbstractFactory<NeuronNetModel>*>();
static DerivedFactory<NeuronNetModel, NeuronNetModel> StdNeurNetModFact;

void InitializeFactories()
{
	// Coupling functions
	AbstractFactory<NetworkEdge>::Factories.insert(make_pair(SigmoidCoupling::ClassName, &SigmoFuncFact));
	AbstractFactory<NetworkEdge>::Factories.insert(make_pair(LinearCoupling::ClassName, &LinFuncFact));
	AbstractFactory<NetworkEdge>::Factories.insert(make_pair(BooleanLink::ClassName, &BooleanLinkFact));
	AbstractFactory<NetworkEdge>::Factories.insert(make_pair(TMSynapse::ClassName, &TMSynEdgeFact));
	AbstractFactory<NetworkEdge>::Factories.insert(make_pair(TMSynapseOptim::ClassName, &TMSynOptEdgeFact));

	// Stimulation strategies
	AbstractFactory<StimulationStrat>::Factories.insert(make_pair(DefaultStimStrat::ClassName, &DefStimStratFact));
	AbstractFactory<StimulationStrat>::Factories.insert(make_pair(PoissonianStimStrat::ClassName, &PoissStimStratFact));
	AbstractFactory<StimulationStrat>::Factories.insert(make_pair(RandomMonoStim::ClassName, &RandomMonoStimStratFact));
	AbstractFactory<StimulationStrat>::Factories.insert(make_pair(ThresholdDeterminationStimStrat::ClassName, &ThresholdDeterminationStimStratFact));
	AbstractFactory<StimulationStrat>::Factories.insert(make_pair(MeanPathStimStrat::ClassName, &MeanPathStimStratFact));
	AbstractFactory<StimulationStrat>::Factories.insert(make_pair(CorrelDetermStimStrat::ClassName, &CorrelDetermStimStratFact));
	AbstractFactory<StimulationStrat>::Factories.insert(make_pair(ExtraCellGluStim::ClassName, &ExtraCellGluStimFact));

	// Networks
	AbstractFactory<AbstractNetwork >::Factories.insert(make_pair(Network<CouplingFunction>::ClassName, &StandardNetFact));
	AbstractFactory<AbstractNetwork >::Factories.insert(make_pair(SpatialNetwork<CouplingFunction>::ClassName, &SpatialNetFact));
	AbstractFactory<AbstractNetwork >::Factories.insert(make_pair(Network<BooleanLink>::ClassName, &BooleanStandardNetFact));
	AbstractFactory<AbstractNetwork >::Factories.insert(make_pair(SpatialNetwork<BooleanLink>::ClassName, &BooleanSpatialNetFact));
	AbstractFactory<AbstractNetwork >::Factories.insert(make_pair(SpatialNetwork<Synapse>::ClassName, &SpatialSynaptNetFact));

	// Neurons
	AbstractFactory<Neuron>::Factories.insert(make_pair(SFALIFNeuron::ClassName, &SFALIFNeuronFact));
	AbstractFactory<Neuron>::Factories.insert(make_pair(DummyNeuron::ClassName, &DummyNeuronFact));

	// Synapses
	AbstractFactory<Synapse>::Factories.insert(make_pair(TMSynapse::ClassName, &TMSynapseFact));
	AbstractFactory<Synapse>::Factories.insert(make_pair(TMSynapseOptim::ClassName, &TMSynapseOptimFact));
	
	// Network building strategies
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(ScaleFreeStrat::ClassName, &ScaleFreeConstructStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(SpatialScaleFreeStrat::ClassName, &SpatialScaleFreeConstructStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(SpatialRegularStrat::ClassName, &SpatialRegularConstructStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(SpatialConnectionRadiusStrat::ClassName, &SpatialConnectionRadiusConstructStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(VoronoiConstructStrat::ClassName, &VoronoiConstructStratFact));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(ErdosRenyiRandomStrat::ClassName, &ErdosRenyiRandomConstructStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(SmallWorldStrat::ClassName, &SmallWorldConstructStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(ThresholdDeterminationStrat::ClassName, &ThresholdDeterminationConstructStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(FunctionalTopoStrat::ClassName, &FunctionalTopoConstructStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(ConstrFromChISim::ClassName, &ConstrFromChISimStratFact));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(FromFileConstrStrat::ClassName, &ConstrFromFiletratFact));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(FileListConstrStrat::ClassName, &ConstrFromFileListStratFact));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(RandomIsolatedCells::ClassName, &RandIsoCellsStrat));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(StarLikeAstroNetConstrStrat::ClassName, &StartLikeAstrStratFact));
	AbstractFactory<NetworkConstructStrat>::Factories.insert(make_pair(ShellStructScramblerStrat::ClassName, &ShellStructScramblerStratFact));

	// Spatial structures building strategies
	AbstractFactory<SpatialStructureBuilder>::Factories.insert(make_pair(StandardStructureBuilder::ClassName, &StandardSructureStrat));
	AbstractFactory<SpatialStructureBuilder>::Factories.insert(make_pair(RandomStructureBuilder::ClassName, &RandomSructureStrat));
	AbstractFactory<SpatialStructureBuilder>::Factories.insert(make_pair(LoadFromFileStructureBuilder::ClassName, &LoadFromFileSructureStrat));
	AbstractFactory<SpatialStructureBuilder>::Factories.insert(make_pair(StarLikeAstrocyteStruct::ClassName, &StarLikeSructureStrat));

	// Metrics
	AbstractFactory<Metric>::Factories.insert(make_pair(DegreeDistrComp::ClassName, &DegreeDistCompFact));
	AbstractFactory<Metric>::Factories.insert(make_pair(ClusteringCoeffs::ClassName, &ClustCoeffDistCompFact));
	AbstractFactory<Metric>::Factories.insert(make_pair(AllPairDistances::ClassName, &AllPairDistancesCompFact));
	AbstractFactory<Metric>::Factories.insert(make_pair(AdjMatrixComp::ClassName, &AdjMatrixCompFact));
	AbstractFactory<Metric>::Factories.insert(make_pair(PositionsComp::ClassName, &PositionsCompFact));
	AbstractFactory<Metric>::Factories.insert(make_pair(StimulatedCells::ClassName, &StimulatedCellsFact));
	AbstractFactory<Metric>::Factories.insert(make_pair(ActivatedCells::ClassName, &ActivatedCellsFact));
	AbstractFactory<Metric>::Factories.insert(make_pair(ConcentrationMetrics::ClassName, &ConcentrationMetricsFact));
	AbstractFactory<Metric>::Factories.insert(make_pair(PropagationDistance::ClassName, &PropagationDistanceFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(PropagationDistanceStat::ClassName, &PropagationDistanceStatFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(NetworkTopologyStats::ClassName, &NetworkTopologyStatsFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(CellStateSaver<BooleanWrapper, BooleanLink>::ClassName, &CellStateSaverBoolStBoolLnkFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(NetworkDimensions::ClassName, &NetworkDimensionsFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(GridSearchMetric::ClassName, &GridSearchMetricFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(ThresholdDetermination::ClassName, &ThresholdDeterminationFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(MeanPathStimMetric::ClassName, &MeanPathStimMetricFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(CellStateSaver<SERStates, BooleanLink>::ClassName, &CellStateSaverSERSStBoolLnkFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(CorrelationsMetric::ClassName, &CorrelationsMetricFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(TransferEntropyMetric::ClassName, &TransferEntropyMetricFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(FunctTopoByNetConstrStrat::ClassName, &FunctTopoByNetConstrMetricFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(WaveFrontDetect::ClassName, &WaveFrontDetectFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(ScalarStatisticsMetric::ClassName, &ScalarStatisticsMetricFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(FourierTransform::ClassName, &FourierTransFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(WaveletTransform::ClassName, &WaveletTransFact));  
	AbstractFactory<Metric>::Factories.insert(make_pair(StochResMetric::ClassName, &StochResMetrFact));


	// Simulations
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(GridSearchSimulation::ClassName, &GridSearchSimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<ChIModel>::ClassName, &RepeatChISimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<ThresholdModel>::ClassName, &RepeatThreshSimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<ModifiedThresholdModel>::ClassName, &RepeatModifThreshSimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<ChIModelThresholdDetermination>::ClassName, &RepeatChIThreshDetermSimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<SecondOrdNeighbSERSModel>::ClassName, &RepeatSecondOrdSERSSimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<SimpleThresholdSERSModel>::ClassName, &RepeatSimpleThreshSERSSimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<NeuronNetModel>::ClassName, &RepeatNeuronNetSimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<AstroNeuroNetModel>::ClassName, &RepeatAstroNeuroSimFact));
	AbstractFactory<SimulationManager>::Factories.insert(make_pair(RepeatSimulation<ChIModelStochasticRes>::ClassName, &RepeatChIStochResFact));

	// ODE Solvers
	AbstractFactory<ODE::ODESolver<double, double> >::Factories.insert(make_pair(ODE::RungeKuttaSolver<double>::ClassName, &ODERungeKuttaSolverFact));
	AbstractFactory<ODE::ODESolver<double, double> >::Factories.insert(make_pair(ODE::EulerSolver<double>::ClassName, &ODEEulerSolverFact));

	// ChI Models
	AbstractFactory<ChIModel>::Factories.insert(make_pair(ChIModel::ClassName, &StdChIModFact));

	// Neuron net Models
	AbstractFactory<NeuronNetModel>::Factories.insert(make_pair(NeuronNetModel::ClassName, &StdNeurNetModFact));
}


