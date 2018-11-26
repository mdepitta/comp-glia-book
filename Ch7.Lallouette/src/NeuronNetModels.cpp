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

#include "NeuronNetModels.h"
#include "ODESolvers.h"
#include "utility.h"
#include "AbstractFactory.h"
#include "ErrorCodes.h"
#include "ChIModelMetrics.h"

#include <assert.h>

using namespace AstroModel;
using namespace std;

//********************************************************************//
//*************** N E U R O N   N E T W O R K   M O D E L ************//
//********************************************************************//

std::string NeuronNetModel::ClassName = "NeuronNetModel";

//**********************************************************************
// Default Constructor
//**********************************************************************
NeuronNetModel::NeuronNetModel(ParamHandler & h) : 
	NetworkDynamicsModel::NetworkDynamicsModel(h),
	ODE::ODEProblem<double, double>::ODEProblem(0, h),
	solver(0)
{
	TRACE("*** Creating Neuron Net Model ***")

	neuronClassName = param.getParam<std::string>("-NeuronClass");
	synapseClassName = param.getParam<std::string>("-SynapseClass");
	network->SetEdgeClassName(synapseClassName);

	// Initialize ODEProblem
	Initialize();
	SetUpNeuronsAndODEs(network->size());
	InitializeDynVals();

	// ODE Solver
	string solverClassName = param.getParam<string>("-SolverClass");
	if (AbstractFactory<ODE::ODESolver<double, double> >::
		Factories[solverClassName])
	{
		solver = AbstractFactory<ODE::ODESolver<double, double> >::
			Factories[solverClassName]->Create();
		solver->SetStepSize(param.getParam<double>("-Step"));
	}
	else
		cerr << "Couldn't create the following ODE solver : " 
			<< solverClassName << endl;
}

//**********************************************************************
// Loading Constructor
//**********************************************************************
NeuronNetModel::NeuronNetModel(std::ifstream & stream, ParamHandler & h) : 
	NetworkDynamicsModel::NetworkDynamicsModel(h),
	ODE::ODEProblem<double, double>::ODEProblem(0),
	solver(0)
{
	if (not this->LoadFromStream(stream))
		cerr << "Failed to load the model !" << endl;
}

//**********************************************************************
// Destructor
//**********************************************************************
NeuronNetModel::~NeuronNetModel()
{
	// ODE related
	if (solver)
		delete solver;

	// Neurons
	ClearNeurons();
}

//**********************************************************************
// Sets up cells and allocate data for ODEs
//**********************************************************************
void NeuronNetModel::SetUpNeuronsAndODEs(unsigned int nbNeur)
{
	assert(AbstractFactory<Neuron>::Factories[neuronClassName]);

	for (unsigned int i = 0 ; i < neurons.size() ; ++i)
		delete neurons[i];
	neurons.clear();

	// Create neurons
	for (unsigned int i = 0 ; i < nbNeur ; ++i)
	{
		neurons.push_back(AbstractFactory<Neuron>::
			Factories[neuronClassName]->Create());
	}

	SetUpODEsAndSynapticStruct();
}

//**********************************************************************
//**********************************************************************
void NeuronNetModel::SetUpODEsAndSynapticStruct()
{
	// Should only be called after network initialization
	SetFunct(new ODE::NeuronNetworkFunc(*this), true);

	// Add Synapses from network
	synapses = network->GetAllocatedLinks();

	AllocateMemory(GetTotNbDynVals());

	// Calls SetVals with 0 / 0 arguments (default args)
	// So that it doesn't change the allocated vals.
	SetVals();
}

//**********************************************************************
// Gives the total number of desired dyn vals (not equivalent to 
// GetNbVals from OPEProblem<double, double>
//**********************************************************************
unsigned int NeuronNetModel::GetTotNbDynVals() const
{
	unsigned nbNeurDynVals = 0;
	unsigned nbSynDynVals = 0;

	for (unsigned int i = 0 ; i < neurons.size() ; ++i)
		nbNeurDynVals += neurons[i]->GetNbDynVal();

	for (unsigned int i = 0 ; i < synapses.size() ; ++i)
		nbSynDynVals += synapses[i]->GetNbDynVal();

	return nbNeurDynVals + nbSynDynVals;
}

//**********************************************************************
// Change ODE vals to given pointer (for external use)
//**********************************************************************
void NeuronNetModel::SetVals(double *_vals, unsigned int _nbVals)
{
	ODE::ODEProblem<double, double>::SetVals(_vals, _nbVals);

	// Update neurons for correct val pointers and synapses
	unsigned int nbNeurDynVals = 0;
	for (unsigned int i = 0 ; i < neurons.size() ; ++i)
		neurons[i]->ClearSynapses();
	for (unsigned int i = 0 ; i < neurons.size() ; ++i)
	{
		neurons[i]->SetDynVals(this->vals + nbNeurDynVals, false);
		neurons[i]->SetModel(this);
		neurons[i]->SetValNamesPostfix(neurons[i]->GetClassName() + "_" +
			StringifyFixed(i));
		// Add axonal and dendritic synapses and update synaptic
		// pointers to pre and post neurons
		for (unsigned int j = 0 ; j < (*network)[i].size() ; ++j)
			if ((*network)[i][j])
			{
				neurons[i]->AddAxonSyn((*network)[i][j]);
				neurons[j]->AddDendrSyn((*network)[i][j]);
				(*network)[i][j]->SetPreSynNeur(neurons[i]);
				(*network)[i][j]->SetPostSynNeur(neurons[j]);
			}
		nbNeurDynVals += neurons[i]->GetNbDynVal();
	}

	// Update synapses for correct val pointers
	unsigned int nbSynDynVals = nbNeurDynVals;
	for (unsigned int i = 0 ; i < synapses.size() ; ++i)
	{
		synapses[i]->SetDynVals(this->vals + nbSynDynVals, false);
		synapses[i]->SetModel(this);
		synapses[i]->SetValNamesPostfix(synapses[i]->GetClassName() + "_" +
			StringifyFixed(i));
		nbSynDynVals += synapses[i]->GetNbDynVal();
	}

	// Initializing init vals
	for (unsigned int i = 0 ; i < this->nbVals ; ++i)
		this->initVals[i] = this->vals[i];
}

//**********************************************************************
// Initializes the model
//**********************************************************************
void NeuronNetModel::Initialize(ResultSaver saver)
{
	NetworkDynamicsModel::Initialize(saver);
	ODE::ODEProblem<double, double>::Initialize(saver);

	SetUpODEsAndSynapticStruct();

	InitializeDynVals();
	// Initialize metrics
	metrics.InitializeMetricsDefault();
}

//**********************************************************************
// Initializes the dynamic values
//**********************************************************************
void NeuronNetModel::InitializeDynVals()
{
	// Initialize cells
	for (unsigned int i = 0 ; i < neurons.size() ; ++i)
		neurons[i]->Initialize();
	// Initialize synapses
	for (unsigned int i = 0 ; i < synapses.size() ; ++i)
		synapses[i]->Initialize();
}

//**********************************************************************
// Method called before solveing the ODE problem
//**********************************************************************
bool NeuronNetModel::PreSimulationCall(ResultSaver saver)
{
	bool ok = true;
	// Save wanted data
TRACE_UP("*** Saving metric data ***")
	ok &= this->network->ComputeAndSaveMetrics(saver);
	ok &= metrics.ComputeMetrics<NeedFrequentUpdateMetric>(*this);
TRACE_DOWN("*** metric data saved ***")
	return ok;
}

//**********************************************************************
// Method called before solving the ODE problem
//**********************************************************************
bool NeuronNetModel::PostSimulationCall(ResultSaver saver)
{
	bool ok = true;

TRACE_UP("*** Computing and saving after simulation metric ***")
	ok &= metrics.ComputeMetrics<AfterSimMetric>(*this);
	ok &= metrics.SaveMetrics<AfterSimMetric>(saver);
TRACE_DOWN("*** After simulation metric data saved ***")

	// Save computed dynamic data and metrics
TRACE("*** Saving dynamic data ***")
	ok &= metrics.SaveMetrics<DynMetric>(saver);
	ok &= ODE::ODEProblem<double, double>::SaveMetrics(saver);

	return ok;

}

//**********************************************************************
// Launch the model simulation with passed solver
//**********************************************************************
int NeuronNetModel::Simulate(ResultSaver saver)
{
	int returnVal = 0;

	PreSimulationCall(saver);

	// Simulate
TRACE_UP("*** Starting Simulation tStart = " << tStart << " tEnd = " << tEnd << " ***")
	assert(solver);
	solver->Solve(*this, tStart, tEnd);
TRACE_DOWN("*** Simulation Ended ***")

	if (not PostSimulationCall(saver))
		returnVal |= MODEL_METRIC_SAVING_PROBLEM;

	return returnVal;
}

//**********************************************************************
// Notifies the model that all values have been updated for timestep t
//**********************************************************************
void NeuronNetModel::UpdateVals(double t)
{
	tCurr = t; 
	// Check if neurons fired
	for (unsigned int i = 0 ; i < neurons.size() ; ++i)
	{
		if (neurons[i]->CheckForSpiking(t))
		{
			//TRACE("t = " << t << " // Spike ! // " << i)
		}
	}

	metrics.ComputeMetrics<NeedFrequentUpdateMetric>(*this);
	ODE::ODEProblem<double, double>::UpdateVals(t);
}

//**********************************************************************
// Set all cells to equilibrium
//**********************************************************************
void NeuronNetModel::SetAllCellsToEquilibrium()
{
	solver->ResetVals();
}

//**********************************************************************
// Dynamically dispatch a metric according to its type
//**********************************************************************
bool NeuronNetModel::AddMetric(Metric *_m, bool _f)
{
	if (NetworkDynamicsModel::AddMetric(_m, _f))
		return true;
	else if (ODE::ODEProblem<double, double>::AddMetric(_m, _f))
		return true;
	else if (metrics.AddMetricAndDependencies(_m, _f, this))
		return true;

	return false;
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool NeuronNetModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = NetworkDynamicsModel::LoadFromStream(stream);
	ok &= ODE::ODEProblem<double, double>::LoadFromStream(stream);
	if (ok)
	{
		stream >> neuronClassName;
		stream >> synapseClassName;

		// Prepare network for correct synapse creation after initialization
		network->SetEdgeClassName(synapseClassName);

		// Neurons
		ClearNeurons();

		unsigned int nbNeur;
		stream >> nbNeur;
		std::string neurClName;
		for (unsigned int i = 0 ; i < nbNeur ; ++i)
		{
			stream >> neurClName;
			assert(AbstractFactory<Neuron>::Factories[neurClName]);
			neurons.push_back(AbstractFactory<Neuron>::
				Factories[neurClName]->CreateFromStream(stream));
		}
		SetUpODEsAndSynapticStruct();

		// ODE Solver
		if (solver)
			delete solver;
		std::string solvName;
		stream >> solvName;

		assert((AbstractFactory<ODE::ODESolver<double, double> >::Factories[solvName]));
		solver = AbstractFactory<ODE::ODESolver<double, double> >::
			Factories[solvName]->CreateFromStream(stream);

		// metrics
		metrics.FreeAndClean();
		ok &= metrics.LoadFromStream(stream);

		return ok and (stream.good() or stream.eof());
	}
	else 
		return false;
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool NeuronNetModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = NetworkDynamicsModel::SaveToStream(stream);
	ok &= ODE::ODEProblem<double, double>::SaveToStream(stream);
	stream
		<< neuronClassName << endl
		<< synapseClassName << endl;
	// Neurons
	stream << neurons.size() << endl;
	for (unsigned int i = 0 ; i < neurons.size() ; ++i)
	{
		stream << neurons[i]->GetClassName() << endl;
		ok &= neurons[i]->SaveToStream(stream);
	}

	// Solver
	assert(solver);
	stream << solver->GetClassName() << std::endl;
	solver->SaveToStream(stream);

	// Metrics
	ok &= metrics.SaveToStream(stream);

	return ok;
}

//**********************************************************************
// Returns a ParamHandler object with references to internal parameters
//**********************************************************************
ParamHandler NeuronNetModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += NetworkDynamicsModel::BuildModelParamHandler();
	params += ODE::ODEProblem<double, double>::BuildModelParamHandler();

	if (not neurons.empty())
		params += neurons.back()->BuildModelParamHandler();
	if (not synapses.empty())
		params += synapses.back()->BuildModelParamHandler();

	params += metrics.BuildModelParamHandler();
	return params;
}

//**********************************************************************
// Returns all metrics and submetrics
//**********************************************************************
vector<Metric *> NeuronNetModel::GetAllMetrics() const
{
	vector<Metric *> mTot, mTemp;

	mTemp = NetworkDynamicsModel::GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	mTemp = ODE::ODEProblem<double, double>::GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	for (unsigned int i = 0 ; i < metrics.GetMetricsRaw().size() ; ++i)
		mTot.push_back(metrics.GetMetricsRaw()[i]);
	return mTot;
}

//**********************************************************************
// Deletes allocated neurons and clear vector
//**********************************************************************
void NeuronNetModel::ClearNeurons()
{
	for (unsigned int i = 0 ; i < neurons.size() ; ++i)
		delete neurons[i];
	neurons.clear();
}


//**********************************************************************
//**********************************************************************
Synapse * NeuronNetModel::GetSynapse(unsigned int ind)
{
	assert(ind < synapses.size());
	return synapses[ind];
}

//**********************************************************************
// Return the indice of the given synapse
//**********************************************************************
unsigned int NeuronNetModel::GetSynInd(Synapse * _s) const
{
	for (unsigned int i = 0 ; i < synapses.size() ; ++i)
		if (synapses[i] == _s)
			return i;
	return DEFAULT_MAX_VAL;
}

