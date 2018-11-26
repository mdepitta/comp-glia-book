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

#include "ChIModel.h"
#include "ODESolvers.h"
#include "ChICell.h"
#include "utility.h"
#include "AbstractFactory.h"
#include "SpatialStructureBuilder.h"
#include "SpatialNetwork.h"
#include "ChIModelMetrics.h"
#include "StimulationMetrics.h"
#include "ErrorCodes.h"

#include <assert.h>

using namespace AstroModel;
using namespace std;

//********************************************************************//
//********************** C H I  M O D E L ****************************//
//********************************************************************//

string ChIModel::ClassName = "ChIModel";

//**********************************************************************
// Default Constructor
//**********************************************************************
ChIModel::ChIModel(ParamHandler & h) : 
	ODENetworkDynamicsModel<CouplingFunction, ChICell>::ODENetworkDynamicsModel(h),
	StimulableCellNetwork::StimulableCellNetwork(h),
	C0(2.0e-03), d1(0.13e-03), d2(1.049e-03), d3(0.9434e-03), 
	d5(0.08234e-03), v3k(4.5e-03), K3k(0.7e-03), r5p(0.21)
{
	TRACE("*** Initializing ChI Model ***")
	SetFunct(new ODE::ChINetworkFunct(*this), true);
}

//**********************************************************************
// Loading Constructor
//**********************************************************************
ChIModel::ChIModel(std::ifstream & stream, ParamHandler & h) : 
	ODENetworkDynamicsModel<CouplingFunction, ChICell>::ODENetworkDynamicsModel(h),
	StimulableCellNetwork::StimulableCellNetwork(h)
{
	if (not LoadFromStream(stream))
		cerr << "Failed to load the model !" << endl;
}

//**********************************************************************
// Destructor
//**********************************************************************
ChIModel::~ChIModel()
{
}

//**********************************************************************
// Sets up cells and allocate data for ODEs
//**********************************************************************
void ChIModel::SetUpCellsAndODEs(unsigned int nbCells)
{
	SetFunct(new ODE::ChINetworkFunct(*this), true);
	ODENetworkDynamicsModel<CouplingFunction, ChICell>::SetUpCellsAndODEs(nbCells);
}

//**********************************************************************
// Initializes the model
//**********************************************************************
void ChIModel::Initialize(ResultSaver saver)
{
	ODENetworkDynamicsModel<CouplingFunction, ChICell>::Initialize(saver);
	Stimulable::Initialize();

	// Initialize metrics
	metrics.InitializeMetricsDefault();
}

//**********************************************************************
// Method called before solving the ODE problem
//**********************************************************************
bool ChIModel::PreSimulationCall(ResultSaver saver)
{
	bool ok = ODENetworkDynamicsModel<CouplingFunction, ChICell>::PreSimulationCall(saver);

	if (preRunToEqu)
	{
		// Temporarily deactivate metrics and stimStrats
		std::vector<StimulationStrat *> tmpStimStrats = stimStrats;
		stimStrats.clear();
		isPreRunning = true;

		solver->Solve(*this, tStart, tStart + preRunTime);

		stimStrats = tmpStimStrats;
		isPreRunning = false;
	}
	ODE::ODEProblem<double, double>::UseCurrentValsAsInitVals();

	// Save wanted data
	ok &= metrics.ComputeMetrics<NeedFrequentUpdateMetric>(*this);

	return ok;
}

//**********************************************************************
// Method called before solving the ODE problem
//**********************************************************************
bool ChIModel::PostSimulationCall(ResultSaver saver)
{
	bool ok = true;

TRACE_UP("*** Computing and saving after simulation ChIMetrics ***")
	ok &= ODENetworkDynamicsModel<CouplingFunction, ChICell>::PostSimulationCall(saver);
	ok &= metrics.ComputeMetrics<ChIModelAfterSimMetric>(*this);
	ok &= metrics.SaveMetrics<ChIModelAfterSimMetric>(saver);
TRACE_DOWN("*** After simulation metric data saved ***")

	// Save computed dynamic data and metrics
TRACE_UP("*** Saving dynamic data ***")
	ok &= metrics.SaveMetrics<ChIModelDynMetric>(saver);
TRACE("ChIModelDynMetric saved.")

	Stimulable::SaveStimStratMetrics(saver);
TRACE_DOWN("*** Finished saving dynamic data ***")

	return ok;
}

//**********************************************************************
// Computes fluxes across cells
//**********************************************************************
void ChIModel::ComputeFluxes(double )
{
	const std::vector<std::vector<unsigned int> > & neighbors = network->GetAllNeighbors();
	for (unsigned int i = 0 ;  i < neighbors.size() ; ++i)
	{
		cells[i]->totFlux = 0;
		cells[i]->caSpontLeak = false;
		for (unsigned int j = 0 ; j < neighbors[i].size() ; ++j)
			cells[i]->totFlux += (*((*network)[i][neighbors[i][j]]))(
					cells[i]->dynVals[ChICell::IP3] - cells[neighbors[i][j]]->dynVals[ChICell::IP3]);
	}
}

//**********************************************************************
// Changes the flux of a cell
//**********************************************************************
void ChIModel::ModifFluxes(unsigned int i, double flux)
{
	//assert(i < cells.size());
	cells[i]->totFlux += flux;
}

//**********************************************************************
// Stimulates Ca2+ release from an astrocyte
//**********************************************************************
void ChIModel::StimSpontaneousCa(unsigned int i)
{
	cells[i]->caSpontLeak = true;
}

//**********************************************************************
// Dynamically dispatch a metric according to its type
//**********************************************************************
bool ChIModel::AddMetric(Metric *_m, bool _f)
{
	if (ODENetworkDynamicsModel<CouplingFunction, ChICell>::AddMetric(_m, _f))
		return true;
	else if (Stimulable::AddMetric(_m, _f))
		return true;
	else
		return ChIModel::metrics.AddMetricAndDependencies(_m, _f, this);
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool ChIModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = ODENetworkDynamicsModel<CouplingFunction, ChICell>::LoadFromStream(stream);
	ok &= Stimulable::LoadFromStream(stream);
	if (ok)
	{
		// Model Parameters
		stream >> C0;
		stream >> d1;
		stream >> d2;
		stream >> d3;
		stream >> d5;
		stream >> v3k;
		stream >> K3k;
		stream >> r5p;

		// metrics
		metrics.FreeAndClean();
		ok &= metrics.LoadFromStream(stream);

		return ok and (stream.good() or stream.eof());
	}
	else 
		return false;
}

//**********************************************************************
// Load funct from stream
//**********************************************************************
bool ChIModel::LoadFunctFromStream(std::ifstream & )
{
	// ODE Function
	SetFunct(new ODE::ChINetworkFunct(*this), true);
	return true;
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool ChIModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = ODENetworkDynamicsModel<CouplingFunction, ChICell>::SaveToStream(stream);
	ok &= Stimulable::SaveToStream(stream);
	// Model Parameters
	stream 
		<< C0     << endl
		<< d1     << endl
		<< d2     << endl
		<< d3     << endl
		<< d5     << endl
		<< v3k    << endl
		<< K3k    << endl
		<< r5p    << endl;

	// Metrics
	ok &= metrics.SaveToStream(stream);

	return ok;
}

//**********************************************************************
// Returns a ParamHandler object with references to internal parameters
//**********************************************************************
ParamHandler ChIModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += ODENetworkDynamicsModel<CouplingFunction, ChICell>::BuildModelParamHandler();
	params += Stimulable::BuildModelParamHandler();

	params <= "C0",  C0;
	params <= "d1",  d1; 
	params <= "d2",  d2; 
	params <= "d3",  d3; 
	params <= "d5",  d5; 
	params <= "v3k", v3k;
	params <= "K3k", K3k;
	params <= "r5p", r5p;

	params += metrics.BuildModelParamHandler();
	return params;
}

//**********************************************************************
// Returns all metrics and submetrics
//**********************************************************************
vector<Metric *> ChIModel::GetAllMetrics() const
{
	vector<Metric *> mTot, mTemp;

	mTemp = ODENetworkDynamicsModel<CouplingFunction, ChICell>::GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	mTemp = Stimulable::GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	for (unsigned int i = 0 ; i < metrics.GetMetricsRaw().size() ; ++i)
		mTot.push_back(metrics.GetMetricsRaw()[i]);

	return mTot;
}

//********************************************************************//
// Returns the total fluxes going out of the cell
//********************************************************************//
double ChIModel::GetTotalFlux(unsigned int i) const
{
	return cells[i]->totFlux;
}

//********************************************************************//
// Return the dynamic value that constitutes the excitable part of the system
//********************************************************************//
double ChIModel::GetExcDynVal(unsigned int cellNb) const
{
	return GetDynVal(cellNb, ChICell::Ca);
}

//********************************************************************//
//* C H I  M O D E L   T H R E S H O L D   D E T E R M I N A T I O N *//
//********************************************************************//

string ChIModelThresholdDetermination::ClassName = "ChIModelThresholdDetermination";

//**********************************************************************
// Default Constructor
//**********************************************************************
ChIModelThresholdDetermination::ChIModelThresholdDetermination(ParamHandler & h) :
	ChIModel(h)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ChIModelThresholdDetermination::ChIModelThresholdDetermination(
	std::ifstream & stream, ParamHandler & h) :
	ChIModel(stream, h)
{
	
}

//**********************************************************************
// Destructor
//**********************************************************************
ChIModelThresholdDetermination::~ChIModelThresholdDetermination()
{

}

//**********************************************************************
// Launch simulation
//**********************************************************************
int ChIModelThresholdDetermination::Simulate(ResultSaver saver)
{
TRACE(GetNetwork().GetConstructStrat().GetClassName())
	int returnVal = 0;

	// Save wanted data
TRACE_UP("*** Saving metric data ***")
	this->network->ComputeAndSaveMetrics(saver);
	this->metrics.ComputeMetrics<ChIModelStaticMetric>(*this);
TRACE_DOWN("*** Metric data saved ***")

	// There should be only one stimulation strategy for this model
	assert(stimStrats.size() == 1);

	ThresholdDeterminationStimStrat *specStimStrat = 
		dynamic_cast<ThresholdDeterminationStimStrat *>(this->stimStrats[0]);
	assert(specStimStrat);

	ThresholdDetermination *threshMetr = 
		GetSpecificMetric<Metric, ThresholdDetermination>(this->GetAllMetrics());
	assert(threshMetr);

	bool hasAlreadySpiked = false;

	// Get the total number of neighbors
	unsigned int k = this->GetNetwork().GetNodeDegree(0);
	// For each possible number of stimulated cells
	for (unsigned int i = 1 ; i <= k ; ++i)
	{
		specStimStrat->SetStimulatedCells(i);
		// Initialize cells
		for (unsigned int j = 0 ; j < cells.size() ; ++j)
			cells[j]->Initialize();
		this->metrics.InitializeMetrics<ChIModelDynMetric>();
		specStimStrat->Initialize();
		// Simulating
		assert(solver);
		solver->Solve(*this, tStart, tEnd);

		if ((not hasAlreadySpiked and threshMetr->HasSpiked()) or (i == k))
		{
			hasAlreadySpiked = true;
			if (not metrics.SaveMetrics<ChIModelDynMetric>(saver))
				returnVal |= MODEL_METRIC_SAVING_PROBLEM;
			if (not specStimStrat->SaveMetrics(saver))
				returnVal |= MODEL_METRIC_SAVING_PROBLEM;
			if (not ODE::ODEProblem<double, double>::SaveMetrics(saver))
				returnVal |= MODEL_METRIC_SAVING_PROBLEM;
		}
	}

	// Save computed dynamic data and metrics
TRACE_UP("*** Saving threshold determination data ***")
	if (not this->metrics.SaveMetrics<ChIModelThreshDetermMetric>(saver))
		returnVal |= MODEL_METRIC_SAVING_PROBLEM;

TRACE_DOWN("*** Threshold determination data saved ***")

	return returnVal;
}

//**********************************************************************
// Returns the network central node degree
//**********************************************************************
unsigned int ChIModelThresholdDetermination::GetDegree() const
{
	const ThresholdDeterminationStrat & threshConstrStrat = 
		dynamic_cast<const ThresholdDeterminationStrat &>(GetNetwork().GetConstructStrat());
	return threshConstrStrat.GetDegree();
}

//**********************************************************************
// Returns the number of derivations (sinks) for stim cells
//**********************************************************************
unsigned int ChIModelThresholdDetermination::GetNbDerivs() const
{
	const ThresholdDeterminationStrat & threshConstrStrat = 
		dynamic_cast<const ThresholdDeterminationStrat &>(GetNetwork().GetConstructStrat());
	return threshConstrStrat.GetNbDerivs();
}

//********************************************************************//
//* C H I  M O D E L   S A M E   N E T   P A R A M   C O M B         *//
//********************************************************************//

string ChIModelStochasticRes::ClassName = "ChIModelStochasticRes";

//**********************************************************************
// Default Constructor
//**********************************************************************
ChIModelStochasticRes::ChIModelStochasticRes(ParamHandler & h) :
	ChIModel(h)
{
	signalStimStratClassName = param.getParam<std::string>("-StochRes", 0);
	noiseStimStratClassName = param.getParam<std::string>("-StochRes", 1);
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ChIModelStochasticRes::ChIModelStochasticRes(
	std::ifstream & stream, ParamHandler & h) :
	ChIModel(h)
{
TRACE("Building stoch res")
	if (not LoadFromStream(stream))
		cerr << "Failed to load the model !" << endl;
}

//**********************************************************************
// Destructor
//**********************************************************************
ChIModelStochasticRes::~ChIModelStochasticRes()
{

}

//**********************************************************************
// Launch simulation
//**********************************************************************
int ChIModelStochasticRes::Simulate(ResultSaver saver)
{
	int returnVal = 0;

	this->metrics.InitializeMetrics<ChIModelDynMetric>();

	PreSimulationCall(saver);
	InitializeDynVals();

	ParamHandler params = BuildModelParamHandler();

	////////////////////
	// Signal only    //
	////////////////////

	params.SetVal(signalStimStratClassName + "_Activ", "1");
	params.SetVal(noiseStimStratClassName + "_Activ", "0");
	TRACE(signalStimStratClassName + "_Activ")
	TRACE(noiseStimStratClassName + "_Activ")
	// First stimulation strategies initialization
	InitializeStimStrats(true);

	solver->Solve(*this, tStart, tEnd);
	this->metrics.ComputeMetrics<StochResMetric>(*this);

	PostSimulationCall(saver("Signal_only"));
	InitializeDynVals();
	SetAllCellsToEquilibrium();
	this->metrics.InitializeMetrics<ChIModelDynMetric>();
	ODE::ODEProblem<double, double>::Initialize(saver);

	////////////////////
	// Noise only     //
	////////////////////

	params.SetVal(signalStimStratClassName + "_Activ", "0");
	params.SetVal(noiseStimStratClassName + "_Activ", "1");
	InitializeStimStrats();

	solver->Solve(*this, tStart, tEnd);
	this->metrics.ComputeMetrics<StochResMetric>(*this);

	PostSimulationCall(saver("Noise_only"));
	InitializeDynVals();
	SetAllCellsToEquilibrium();
	this->metrics.InitializeMetrics<ChIModelDynMetric>();
	ODE::ODEProblem<double, double>::Initialize(saver);

	////////////////////
	// Signal + noise //
	////////////////////

	params.SetVal(signalStimStratClassName + "_Activ", "1");
	params.SetVal(noiseStimStratClassName + "_Activ", "1");
	InitializeStimStrats();

	solver->Solve(*this, tStart, tEnd);
	this->metrics.ComputeMetrics<StochResMetric>(*this);

	PostSimulationCall(saver("Signal_and_noise"));
	InitializeDynVals();
	SetAllCellsToEquilibrium();
	this->metrics.InitializeMetrics<ChIModelDynMetric>();
	ODE::ODEProblem<double, double>::Initialize(saver);

	////////////////////
	// Other metrics  //
	////////////////////

	if (not PostSimulationCall(saver))
		returnVal |= MODEL_METRIC_SAVING_PROBLEM;

	return returnVal;
}

//**********************************************************************
// Initialize stimulation strategies
//**********************************************************************
void ChIModelStochasticRes::InitializeStimStrats(bool firstInit)
{
	RandomStimStrat *tmp = 0;
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
	{
		if (not firstInit and (tmp = dynamic_cast<RandomStimStrat*>(stimStrats[i])))
			tmp->InitializeWithSameRandomSequence();
		else
			stimStrats[i]->Initialize();
	}
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool ChIModelStochasticRes::LoadFromStream(std::ifstream & stream)
{
TRACE("here")
	bool ok = ChIModel::LoadFromStream(stream);
TRACE("there")

	stream >> signalStimStratClassName;
	stream >> noiseStimStratClassName;

	return ok and stream.good();
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool ChIModelStochasticRes::SaveToStream(std::ofstream & stream) const
{
	bool ok = ChIModel::SaveToStream(stream);

	stream << signalStimStratClassName << std::endl;
	stream << noiseStimStratClassName << std::endl;

	return ok and stream.good();
}

