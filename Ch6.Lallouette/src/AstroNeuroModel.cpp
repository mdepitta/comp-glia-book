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

#include "AstroNeuroModel.h"
#include "ODESolvers.h"
#include "utility.h"
#include "AbstractFactory.h"
#include "ErrorCodes.h"
#include "ChIModel.h"
#include "NeuronNetModels.h"

#include <assert.h>
#include <string>

using namespace AstroModel;
using namespace std;

//********************************************************************//
//*************** N E U R O N   N E T W O R K   M O D E L ************//
//********************************************************************//

std::string AstroNeuroNetModel::ClassName = "AstrocyteNeuronNetModel";

//**********************************************************************
// Default Constructor
//**********************************************************************
AstroNeuroNetModel::AstroNeuroNetModel(ParamHandler & h) : 
	Model::Model(h),
	ODE::ODEProblem<double, double>::ODEProblem(0, h),
	solver(0)
{
	TRACE("*** Creating Neuron Net Model ***")


	neuronNetClassName = param.getParam<std::string>("-NeuronNetClass");
	astroNetClassName = param.getParam<std::string>("-AstroNetClass");

	int Nneur = param.getParam<int>("-Nneur");
	param.SetVal("-N", Stringify(Nneur).c_str(), 0);
	std::string _pNeur = param.getParam<std::string>("-FrmFileNeurNet");
	param.SetVal("-frmFilePath", _pNeur.c_str());
	param.SetVal("-frmFileDirected", "1");
	std::string _posNeur = param.getParam<std::string>("-FrmFileNeurPos");
	param.SetVal("-LoadFrmFileStructBuild", _posNeur.c_str());
	param.SetVal("-Construct", param.getParam<std::string>("-NeurNetConstruct").c_str());
	assert(AbstractFactory<NeuronNetModel>::Factories[neuronNetClassName]);
	neuronNet = AbstractFactory<NeuronNetModel>::Factories[neuronNetClassName]->Create();


	int Nastr = param.getParam<int>("-Nastro");
	h.SetVal("-N", Stringify(Nastr).c_str(), 0);
	std::string _pAstr = param.getParam<std::string>("-FrmFileAstrNet");
	param.SetVal("-frmFilePath", _pAstr.c_str());
	param.SetVal("-frmFileDirected", "0");
	std::string _posAstr = param.getParam<std::string>("-FrmFileAstrPos");
	param.SetVal("-LoadFrmFileStructBuild", _posAstr.c_str());
	param.SetVal("-Construct", param.getParam<std::string>("-AstrNetConstruct").c_str());
	assert(AbstractFactory<ChIModel>::Factories[astroNetClassName]);
	astroNet = AbstractFactory<ChIModel>::Factories[astroNetClassName]->Create();

	// Initialize ODEProblem and both models
	SetUpModelsAndODEs();

	astrToSyn = std::vector<std::vector<Synapse *> >(astroNet->GetNbCells(), std::vector<Synapse *>());
	std::string tmpStr = ParamHandler::GlobalParams.getParam<std::string>("-AstrToSynFile", 0);
	std::ifstream tmpStream(tmpStr.c_str());
	unsigned int nbAstr, nbSyn, indSyn;
	tmpStream >> nbAstr;
	for (unsigned int i = 0 ; i < nbAstr ; ++i)
	{
		tmpStream >> nbSyn;
		for (unsigned int j = 0 ; j < nbSyn ; ++j)
		{
			tmpStream >> indSyn;
			astrToSyn[i].push_back(neuronNet->GetSynapse(indSyn));
		}
	}


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
AstroNeuroNetModel::AstroNeuroNetModel(std::ifstream & stream, ParamHandler & h) : 
	Model(h),
	ODE::ODEProblem<double, double>::ODEProblem(0),
	solver(0)
{
	if (not this->LoadFromStream(stream))
		cerr << "Failed to load the model !" << endl;
}

//**********************************************************************
// Destructor
//**********************************************************************
AstroNeuroNetModel::~AstroNeuroNetModel()
{
	// ODE related
	if (solver)
		delete solver;
	if (neuronNet)
		delete neuronNet;
	if (astroNet)
		delete astroNet;
}

//**********************************************************************
// Allocate memory and sets up model with corresponding dynVals
//**********************************************************************
void AstroNeuroNetModel::SetUpModelsAndODEs()
{
	SetFunct(new ODE::AstroNeuronNetFunc(*this), true);

	unsigned int nbNeurDynVals = neuronNet->GetTotNbDynVals();
	unsigned int nbAstrDynVals = astroNet->GetTotNbDynVals();

	AllocateMemory(nbNeurDynVals + nbAstrDynVals);

	neuronNet->SetVals(this->vals, nbNeurDynVals);
	neuronNet->InitializeDynVals();

	astroNet->SetVals(this->vals + nbNeurDynVals, nbAstrDynVals);
	astroNet->InitializeDynVals();


	std::vector<std::string> vnams;
	vnams = neuronNet->GetValNames();
	for (unsigned int i = 0 ; i < vnams.size() ; ++i)
		SetValName(i, vnams[i]);

	vnams = astroNet->GetValNames();
	for (unsigned int i = 0 ; i < vnams.size() ; ++i)
		SetValName(nbNeurDynVals + i, vnams[i]);
}

//**********************************************************************
// Initializes the model
//**********************************************************************
void AstroNeuroNetModel::Initialize(ResultSaver saver)
{
	ODE::ODEProblem<double, double>::Initialize(saver);

	// Initialize neuron network
	neuronNet->Initialize(saver);
	// Initialize astrocyte network
	astroNet->Initialize(saver);

	SetUpModelsAndODEs();

	// Initialize metrics
	metrics.InitializeMetricsDefault();
}

//**********************************************************************
// Method called before solving the ODE problem
//**********************************************************************
bool AstroNeuroNetModel::PreSimulationCall(ResultSaver saver)
{
	bool ok = true;
	assert(neuronNet);
	assert(astroNet);
TRACE("Pre-simulation call for " << neuronNet->GetClassName())
	ok &= neuronNet->PreSimulationCall(saver);
TRACE("Pre-simulation call for " << astroNet->GetClassName())
	ok &= astroNet->PreSimulationCall(saver);

	ok &= metrics.ComputeMetrics<NeedFrequentUpdateMetric>(*this);

	return ok;
}

//**********************************************************************
// Method called before solving the ODE problem
//**********************************************************************
bool AstroNeuroNetModel::PostSimulationCall(ResultSaver saver)
{
	bool ok = true;

	assert(neuronNet);
	assert(astroNet);

TRACE("Pre-simulation call for " << neuronNet->GetClassName())
	ok &= neuronNet->PostSimulationCall(saver);
TRACE("Pre-simulation call for " << astroNet->GetClassName())
	ok &= astroNet->PostSimulationCall(saver);

TRACE_UP("*** Computing and saving after simulation metric ***")
	ok &= ODE::ODEProblem<double, double>::ComputeAfterSimMetrics();
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
int AstroNeuroNetModel::Simulate(ResultSaver saver)
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
void AstroNeuroNetModel::UpdateVals(double t)
{
	tCurr = t; 
	neuronNet->UpdateVals(t);
	astroNet->UpdateVals(t);

	metrics.ComputeMetrics<NeedFrequentUpdateMetric>(*this);
	ODE::ODEProblem<double, double>::UpdateVals(t);
}

//**********************************************************************
// Set all cells to equilibrium
//**********************************************************************
void AstroNeuroNetModel::SetAllCellsToEquilibrium()
{
	solver->ResetVals();
}

//**********************************************************************
// Dynamically dispatch a metric according to its type
//**********************************************************************
bool AstroNeuroNetModel::AddMetric(Metric *_m, bool _f)
{
	if (Model::AddMetric(_m, _f))
		return true;
	else if (ODE::ODEProblem<double, double>::AddMetric(_m, _f))
		return true;
	else if (metrics.AddMetricAndDependencies(_m, _f, this))
		return true;
	else if (astroNet and astroNet->AddMetric(_m, _f))
		return true;
	else if (neuronNet and neuronNet->AddMetric(_m, _f))
		return true;

	return false;
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool AstroNeuroNetModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = Model::LoadFromStream(stream);
	ok &= ODE::ODEProblem<double, double>::LoadFromStream(stream);
	if (ok)
	{
		// Neuron Net
		stream >> neuronNetClassName;
		assert(AbstractFactory<NeuronNetModel>::Factories[neuronNetClassName]);
		if (neuronNet)
			delete neuronNet;
		neuronNet = AbstractFactory<NeuronNetModel>::Factories[neuronNetClassName]
			->CreateFromStream(stream);

		// Astrocyte Net
		stream >> astroNetClassName;
		assert(AbstractFactory<ChIModel>::Factories[astroNetClassName]);
		if (astroNet)
			delete astroNet;
		astroNet = AbstractFactory<ChIModel>::Factories[astroNetClassName]
			->CreateFromStream(stream);

		// AstrToSyn
		unsigned int nbAstrToSyn, nbSyn, indSyn;
		stream >> nbAstrToSyn;
		assert(nbAstrToSyn == astroNet->GetNbCells());
		astrToSyn = std::vector<std::vector<Synapse *> >(nbAstrToSyn, 
			std::vector<Synapse *>());
		for (unsigned int i = 0 ; i < nbAstrToSyn ; ++i)
		{
			stream >> nbSyn;
			for (unsigned int j = 0 ; j < nbSyn ; ++j)
			{
				stream >> indSyn;
TRACE(indSyn)
				astrToSyn[i].push_back(neuronNet->GetSynapse(indSyn));
			}
		}

		// ODE Solver
		if (solver)
			delete solver;
		std::string solvName;
		stream >> solvName;
		solver = AbstractFactory<ODE::ODESolver<double, double> >::Factories[solvName]
			->CreateFromStream(stream);

		// metrics
		metrics.FreeAndClean();
		ok &= metrics.LoadFromStream(stream);

		Initialize();
	}
	return ok and (stream.good() or stream.eof());
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool AstroNeuroNetModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = Model::SaveToStream(stream);
	ok &= ODE::ODEProblem<double, double>::SaveToStream(stream);

	// Neuron network
	assert(neuronNet);
	stream << neuronNet->GetClassName() << std::endl;
	ok &= neuronNet->SaveToStream(stream);

	// Astrocyte network
	assert(astroNet);
	stream << astroNet->GetClassName() << std::endl;
	ok &= astroNet->SaveToStream(stream);

	// Astro to synapses connectivity
	assert(astrToSyn.size() == astroNet->GetNbCells());
	stream << astrToSyn.size() << std::endl;
	for (unsigned int i = 0 ; i < astrToSyn.size() ; ++i)
	{
		stream << astrToSyn[i].size() << std::endl;
		for (unsigned int j = 0 ; j < astrToSyn[i].size() ; ++j)
			stream << neuronNet->GetSynInd(astrToSyn[i][j]) << std::endl;
	}

	// solver
	assert(solver);
	stream << solver->GetClassName() << std::endl;
	ok &= solver->SaveToStream(stream);

	// Metrics
	ok &= metrics.SaveToStream(stream);

	return ok and stream.good();
}

//**********************************************************************
// Returns a ParamHandler object with references to internal parameters
//**********************************************************************
ParamHandler AstroNeuroNetModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += ODE::ODEProblem<double, double>::BuildModelParamHandler();

	assert(neuronNet);
	assert(astroNet);
	params.AddAllWithPrefix(neuronNet->BuildModelParamHandler(), "NeuronNet_");
	params.AddAllWithPrefix(astroNet->BuildModelParamHandler(), "AstroNet_");

	params += metrics.BuildModelParamHandler();
	return params;
}

//**********************************************************************
// Returns all metrics and submetrics
//**********************************************************************
vector<Metric *> AstroNeuroNetModel::GetAllMetrics() const
{
	vector<Metric *> mTot, mTemp;

	for (unsigned int i = 0 ; i < metrics.GetMetricsRaw().size() ; ++i)
		mTot.push_back(metrics.GetMetricsRaw()[i]);
	
	mTemp = ODE::ODEProblem<double, double>::GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	assert(neuronNet);
	assert(astroNet);
	mTemp = neuronNet->GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	mTemp = astroNet->GetAllMetrics();
	for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
		mTot.push_back(mTemp[i]);

	return mTot;
}


