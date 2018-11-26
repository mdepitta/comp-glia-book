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

#include "Neuron.h"
#include "Synapse.h"
#include "NeuronNetModels.h"

using namespace AstroModel;
using namespace std;

//********************************************************************//
//********** N E U R O N   M O D E L   B A S E   C L A S S ***********//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
Neuron::Neuron(const NeuronNetModel * _model) : funct(0), model(_model)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
Neuron::Neuron(std::ifstream & )
{

}

//**********************************************************************
// Copy constructor
//**********************************************************************
Neuron::Neuron(const Neuron & n) : model(n.model)
{
}

//**********************************************************************
// Destructor
//**********************************************************************
Neuron::~Neuron()
{
	delete funct;
}

//**********************************************************************
//**********************************************************************
void Neuron::SetModel(const NeuronNetModel *_m)
{
	model = _m;
}

//**********************************************************************
// Returns true if the overall structure formed by the neuron and its
// dendritic and axonal synapses is coherent.
//**********************************************************************
bool Neuron::IsCorrectlyInitialized() const
{
	bool correct = true;
	for (unsigned int i = 0 ; i < dendrSyn.size() ; ++i)
		correct &= (dendrSyn[i]->postSynNeur == this);
	for (unsigned int i = 0 ; i < axonSyn.size() ; ++i)
		correct &= (axonSyn[i]->preSynNeur == this);
	return correct;
}

//**********************************************************************
//**********************************************************************
void Neuron::AddDendrSyn(Synapse *_s)
{
	dendrSyn.push_back(_s);
}

//**********************************************************************
//**********************************************************************
void Neuron::AddAxonSyn(Synapse *_s)
{
	axonSyn.push_back(_s);
}

//**********************************************************************
//**********************************************************************
void Neuron::ClearSynapses()
{
	dendrSyn.clear();
	axonSyn.clear();
}

//**********************************************************************
//**********************************************************************
const NeuronNetModel * Neuron::GetModel() const
{
	return model;
}

//********************************************************************//
//******* D U M M Y   N E U R O N   ( S P I K E   T R A I N ) ********//
//********************************************************************//

std::string DummyNeuron::ClassName = "DummyNeuron";

//**********************************************************************
// Default Constructor
//**********************************************************************
DummyNeuron::DummyNeuron(const NeuronNetModel * _model, double *_dv, 
		bool _fv) : Neuron::Neuron(_model), dynVals(_dv), freeDynVals(_fv), 
	nbDynVals(DUMMYNEURON_NBVALS_PER_NEURON)
{
	funct = new ODE::DummyNeuronFunct(*this);
	if (not _dv)
	{
		dynVals = new double[nbDynVals];
		freeDynVals = true;
	}

	// Load spike trains from file
	std::string spikeTrainPath = ParamHandler::GlobalParams.getParam<std::string>("-DummyNeuronSpikeTrain");
	ifstream spikeTrain(spikeTrainPath.c_str());
	assert(spikeTrain.good());
	double tmpTime;
	unsigned int nbSpikes;
	spikeTrain >> nbSpikes;
	for (unsigned int i = 0 ; i < nbSpikes ; ++i)
	{
		spikeTrain >> tmpTime;
		spikingTimes.push_back(tmpTime);
	}

	Initialize();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
DummyNeuron::DummyNeuron(std::ifstream & stream) : Neuron(stream)
{
	funct = new ODE::DummyNeuronFunct(*this);
	nbDynVals = DUMMYNEURON_NBVALS_PER_NEURON;
	dynVals = new double[nbDynVals];
	freeDynVals = true;
	DummyNeuron::LoadFromStream(stream);
}

//**********************************************************************
// Copy constructor
//**********************************************************************
DummyNeuron::DummyNeuron(const DummyNeuron & s) : Neuron(s), 
	dynVals(s.dynVals), freeDynVals(s.freeDynVals), nbDynVals(s.nbDynVals)
{
	funct = new ODE::DummyNeuronFunct(*this);
	if (freeDynVals)
	{
		dynVals = new double[nbDynVals];
		for (unsigned int i = 0 ; i < nbDynVals ; ++i)
			dynVals[i] = s.dynVals[i];
	}
	Initialize();
}

//**********************************************************************
// Destructor
//**********************************************************************
DummyNeuron::~DummyNeuron()
{
	if (freeDynVals)
		delete[] dynVals;
}

//**********************************************************************
//**********************************************************************
void DummyNeuron::Initialize()
{
	dynVals[0] = -1.0;
	nextSpike = 0;
}

//**********************************************************************
// Check for potential spiking
//**********************************************************************
bool DummyNeuron::CheckForSpiking(double t)
{
	// Get Neuron back to polarized state
	dynVals[0] = -1.0;

	// If the a spike is planned at this time
	if ((nextSpike < spikingTimes.size()) and (t >= spikingTimes[nextSpike]))
	{
		ForceSpiking(spikingTimes[nextSpike]);
		++nextSpike;
		return true;
	}
	return false;
}

//**********************************************************************
// Update value names in ODEProblem
//**********************************************************************
void DummyNeuron::SetValNamesPostfix(std::string pf)
{
	model->AddPostfixToValName(dynVals + V, pf + "_state");
}

//**********************************************************************
// Forces the neuron into emitting a spike, regardless of its state
//**********************************************************************
void DummyNeuron::ForceSpiking(double t)
{
	dynVals[0] = 0;
	// For each axonal synapse
	for (unsigned int i = 0 ; i < axonSyn.size() ; ++i)
		axonSyn[i]->PresynSpike(t);
}

//**********************************************************************
//**********************************************************************
bool DummyNeuron::LoadFromStream(std::ifstream & stream)
{
	unsigned int tmp;
	double time;
	spikingTimes.clear();
	stream >> tmp;
	for (unsigned int i = 0 ; i < tmp ; ++i)
	{
		stream >> time;
		spikingTimes.push_back(time);
	}

	Initialize();
	return not stream.eof();
}

//**********************************************************************
//**********************************************************************
bool DummyNeuron::SaveToStream(std::ofstream & stream) const
{
	stream << spikingTimes.size() << std::endl;
	for (unsigned int i = 0 ; i < spikingTimes.size() ; ++i)
		stream << spikingTimes[i] << std::endl;
	return stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler DummyNeuron::BuildModelParamHandler()
{
	ParamHandler params;
	return params;
}

//**********************************************************************
//**********************************************************************
unsigned int DummyNeuron::GetNbDynVal() const
{
	return DUMMYNEURON_NBVALS_PER_NEURON;
}

//**********************************************************************
// Change dynamic values to given pointer
//**********************************************************************
void DummyNeuron::SetDynVals(double *_dv, bool _f)
{
	if (freeDynVals)
		delete[] dynVals;
	dynVals = _dv;
	freeDynVals = _f;
}


//********************************************************************//
//** S P I K E   F R E Q   A D A P T   L E A K Y   I N T   F I R E ***//
//********************************************************************//

std::string SFALIFNeuron::ClassName = "SFALIFNeuron";

double SFALIFNeuron::DefaultC    = 0.01;                    // F.m^-2
double SFALIFNeuron::DefaultgL   = 0.5;                     // S.m^-2
double SFALIFNeuron::DefaultE0   = -0.060;                  // V
double SFALIFNeuron::DefaultVt   = -0.050;                  // V
double SFALIFNeuron::DefaultVr   = SFALIFNeuron::DefaultE0; //temp
double SFALIFNeuron::Defaulttauw = 0.6;                     // s 
double SFALIFNeuron::Defaultb    = 0.000000000005;          // A // 0.00005;//

double SFALIFNeuron::DefaultV    = SFALIFNeuron::DefaultE0;
double SFALIFNeuron::Defaultw    = 0.0;

//**********************************************************************
// Default Constructor
//**********************************************************************
SFALIFNeuron::SFALIFNeuron(const NeuronNetModel * _model, double *_dv, 
		bool _fv) : Neuron::Neuron(_model),
	defaultBiophysParams(true), C(DefaultC), gL(DefaultgL), E0(DefaultE0),
	Vt(DefaultVt), Vr(DefaultVr), tauw(Defaulttauw), b(Defaultb), inputCurr(0),
	dynVals(_dv), freeDynVals(_fv), 
	nbDynVals(SFALIFNEURON_NBVALS_PER_NEURON)
{
	funct = new ODE::SFALIFNeuronFunct(*this);
	if (not _dv)
	{
		dynVals = new double[nbDynVals];
		freeDynVals = true;
	}
	Initialize();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
SFALIFNeuron::SFALIFNeuron(std::ifstream & stream) : Neuron(stream)
{
	funct = new ODE::SFALIFNeuronFunct(*this);
	nbDynVals = SFALIFNEURON_NBVALS_PER_NEURON;
	dynVals = new double[nbDynVals];
	freeDynVals = true;
	SFALIFNeuron::LoadFromStream(stream);
}

//**********************************************************************
// Copy constructor
//**********************************************************************
SFALIFNeuron::SFALIFNeuron(const SFALIFNeuron & s) : Neuron(s), 
	defaultBiophysParams(s.defaultBiophysParams), C(s.C), gL(s.gL), 
	E0(s.E0), Vt(s.Vt), Vr(s.Vr), tauw(s.tauw), b(s.b), inputCurr(s.inputCurr),
	dynVals(s.dynVals), freeDynVals(s.freeDynVals), nbDynVals(s.nbDynVals)
{
	funct = new ODE::SFALIFNeuronFunct(*this);
	if (freeDynVals)
	{
		dynVals = new double[nbDynVals];
		for (unsigned int i = 0 ; i < nbDynVals ; ++i)
			dynVals[i] = s.dynVals[i];
	}
	Initialize();
}

//**********************************************************************
// Destructor
//**********************************************************************
SFALIFNeuron::~SFALIFNeuron()
{
	if (freeDynVals)
		delete[] dynVals;
}

//**********************************************************************
//**********************************************************************
void SFALIFNeuron::Initialize()
{
	dynVals[0] = DefaultV;
	dynVals[1] = Defaultw;
	if (defaultBiophysParams)
	{
		C    = DefaultC;   
		gL   = DefaultgL;  
		E0   = DefaultE0;  
		Vt   = DefaultVt;  
		Vr   = DefaultVr;  
		tauw = Defaulttauw;
		b    = Defaultb;   
	}
	inputCurr = 0.0;
}

//**********************************************************************
//**********************************************************************
void SFALIFNeuron::SetToEquilibrium()
{
	dynVals[V] = DefaultV;
	dynVals[w] = Defaultw;
}

//**********************************************************************
// Check for potential spiking
//**********************************************************************
bool SFALIFNeuron::CheckForSpiking(double t)
{
	// If the threshold potential for spiking is crossed
	if (dynVals[V] >= Vt)
	{
		ForceSpiking(t);
		return true;
	}
	return false;
}

//**********************************************************************
// Update value names in ODEProblem
//**********************************************************************
void SFALIFNeuron::SetValNamesPostfix(std::string pf)
{
	model->AddPostfixToValName(dynVals + V, pf + "_V");
	model->AddPostfixToValName(dynVals + w, pf + "_w");
}

//**********************************************************************
// Forces the neuron into emitting a spike, regardless of its state
//**********************************************************************
void SFALIFNeuron::ForceSpiking(double t)
{
	dynVals[V] = Vr;
	dynVals[w] += b;
	// For each axonal synapse
	for (unsigned int i = 0 ; i < axonSyn.size() ; ++i)
		axonSyn[i]->PresynSpike(t);
}

//**********************************************************************
//**********************************************************************
bool SFALIFNeuron::LoadFromStream(std::ifstream & stream)
{
	stream >> defaultBiophysParams;
	stream >> C;   
	stream >> gL;  
	stream >> E0;  
	stream >> Vt;  
	stream >> Vr;  
	stream >> tauw;
	stream >> b;   
	Initialize();
	return not stream.eof();
}

//**********************************************************************
//**********************************************************************
bool SFALIFNeuron::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< defaultBiophysParams << endl
		<< C    << endl
		<< gL   << endl
		<< E0   << endl
		<< Vt   << endl
		<< Vr   << endl
		<< tauw << endl
		<< b    << endl;
	return stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler SFALIFNeuron::BuildModelParamHandler()
{
	ParamHandler params;
	params <= "SFALIFDefaultVVal", DefaultV;
	params <= "SFALIFDefaultwVal", Defaultw;
	params <= "SFALIFC"   , DefaultC;   
	params <= "SFALIFgL"  , DefaultgL;  
	params <= "SFALIFE0"  , DefaultE0;  
	params <= "SFALIFVt"  , DefaultVt;  
	params <= "SFALIFVr"  , DefaultVr;  
	params <= "SFALIFtauw", Defaulttauw;
	params <= "SFALIFb"   , Defaultb;   
	return params;
}

//**********************************************************************
//**********************************************************************
unsigned int SFALIFNeuron::GetNbDynVal() const
{
	return SFALIFNEURON_NBVALS_PER_NEURON;
}

//**********************************************************************
// Change dynamic values to given pointer
//**********************************************************************
void SFALIFNeuron::SetDynVals(double *_dv, bool _f)
{
	if (freeDynVals)
		delete[] dynVals;
	dynVals = _dv;
	freeDynVals = _f;
}

