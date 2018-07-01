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

#include "Synapse.h"
#include "Neuron.h"
#include "NeuronNetModels.h"
#include <algorithm>
#include <gsl/gsl_sf_exp.h>

using namespace AstroModel;
using namespace std;

//********************************************************************//
//********** S Y N A P S E   M O D E L   B A S E   C L A S S *********//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
Synapse::Synapse(const NeuronNetModel * _model, double *_dv, bool _fv, 
	unsigned int _nbdv) : 
	funct(0), model(_model), preSynNeur(0), postSynNeur(0), 
	dynVals(_dv), freeDynVals(_fv), nbDynVals(_nbdv)
{
}

//**********************************************************************
// Cosntructor from stream
//**********************************************************************
Synapse::Synapse(std::ifstream & ) : funct(0), model(0), 
	preSynNeur(0), postSynNeur(0), dynVals(0), freeDynVals(false),
	nbDynVals(0)
{
}

//**********************************************************************
// Copy constructor
//**********************************************************************
Synapse::Synapse(const Synapse & s) : model(s.model), 
	dynVals(s.dynVals), freeDynVals(s.freeDynVals), nbDynVals(s.nbDynVals)
{
	if (freeDynVals)
	{
		dynVals = new double[nbDynVals];
		for (unsigned int i = 0 ; i < nbDynVals ; ++i)
			dynVals[i] = s.dynVals[i];
	}
}

//**********************************************************************
// Destructor
//**********************************************************************
Synapse::~Synapse()
{
	delete funct;
	if (freeDynVals)
		delete[] dynVals;
}

//**********************************************************************
//**********************************************************************
void Synapse::SetModel(const NeuronNetModel *_m)
{
	model = _m;
}

//**********************************************************************
// Returns true if the overall structure formed by the synapse and its
// pre and post synaptic neurons is coherent.
//**********************************************************************
bool Synapse::IsCorrectlyInitialized() const
{
	bool correct = true;
	correct &= (preSynNeur != 0);
	correct &= (postSynNeur != 0);
	if (correct)
	{
		correct &= (std::find(preSynNeur->axonSyn.begin(), 
			preSynNeur->axonSyn.end(), this) != 
			preSynNeur->axonSyn.end());
		correct &= (std::find(postSynNeur->dendrSyn.begin(), 
			postSynNeur->dendrSyn.end(), this) != 
			postSynNeur->dendrSyn.end());
	}
	return correct;
}

//**********************************************************************
//**********************************************************************
void Synapse::SetPreSynNeur(Neuron *_n)
{
	preSynNeur = _n;
}

//**********************************************************************
//**********************************************************************
void Synapse::SetPostSynNeur(Neuron *_n)
{
	postSynNeur = _n;
}

//**********************************************************************
//**********************************************************************
const NeuronNetModel * Synapse::GetModel() const
{
	return model;
}

//**********************************************************************
// Change dynamic values to given pointer
//**********************************************************************
void Synapse::SetDynVals(double *_dv, bool _f)
{
	if (freeDynVals)
		delete[] dynVals;
	dynVals = _dv;
	freeDynVals = _f;
}

//********************************************************************//
//* G L U T A M A T E R G I C   S Y N A P S E   B A S E   C L A S S **//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
GlutamatergicSynapse::GlutamatergicSynapse(const NeuronNetModel * _model, double *_dv, bool _fv, 
	unsigned int _nbdv) : Synapse(_model, _dv, _fv, _nbdv)
{
}

//**********************************************************************
// Cosntructor from stream
//**********************************************************************
GlutamatergicSynapse::GlutamatergicSynapse(std::ifstream & stream) :
	Synapse(stream)
{
}


//********************************************************************//
//********** T S O D Y K S   M A R K R A M   S Y N A P S E ***********//
//********************************************************************//

std::string TMSynapse::ClassName = "TMSynapse";

double TMSynapse::Defaultx     = 1.0;
double TMSynapse::Defaultu     = 0.0;
double TMSynapse::Defaultgamma = 0.0;

double TMSynapse::DefaultOd   = 1.0;//1.0 // 50   // temp, for no depletion : 10000.0
double TMSynapse::DefaultOf   = 2.0;//2.0;     // temp, for no depletion : 0.0
double TMSynapse::DefaultU0   = 0.25;    // not definitive
double TMSynapse::DefaultrhoA = 0.00065; // taken from ratio for astro -> syn
double TMSynapse::Defaultnv   = 4.0;     // 
double TMSynapse::DefaultGv   = 50000.0;    // taken for astro vesicles, check
double TMSynapse::DefaultOc   = 60.0; //(was 60.0 before the 28/02/2013)    // taken from ratio for astro -> syn

double TMSynapse::DefaultSpillOvFract = 0.025;

//**********************************************************************
// Default Constructor
//**********************************************************************
TMSynapse::TMSynapse(const NeuronNetModel * _model, double *_dv, 
		bool _fv) : GlutamatergicSynapse::GlutamatergicSynapse(_model, 
	_dv, _fv, TMSYNAPSE_NBVALS_PER_SYN), defaultBiophysParams(true), 
	Od(DefaultOd), Of(DefaultOf), U0(DefaultU0), rhoA(DefaultrhoA), 
	nv(Defaultnv), Gv(DefaultGv), Oc(DefaultOc), 
	spillOvFract(DefaultSpillOvFract)
{
	funct = new ODE::TMSynapseFunct(*this);
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
TMSynapse::TMSynapse(std::ifstream & stream) : GlutamatergicSynapse(stream)
{
	funct = new ODE::TMSynapseFunct(*this);
	nbDynVals = TMSYNAPSE_NBVALS_PER_SYN;
	dynVals = new double[nbDynVals];
	freeDynVals = true;
	TMSynapse::LoadFromStream(stream);
}

//**********************************************************************
// Copy constructor
//**********************************************************************
TMSynapse::TMSynapse(const TMSynapse & s) : GlutamatergicSynapse(s), 
	defaultBiophysParams(s.defaultBiophysParams), Od(s.Od), Of(s.Of),
	U0(s.U0), rhoA(s.rhoA), nv(s.nv), Gv(s.Gv), Oc(s.Oc), 
	spillOvFract(s.spillOvFract)
{
	funct = new ODE::TMSynapseFunct(*this);
	Initialize();
}

//**********************************************************************
// Destructor
//**********************************************************************
TMSynapse::~TMSynapse()
{
}

//**********************************************************************
//**********************************************************************
void TMSynapse::Initialize()
{
	dynVals[0] = Defaultx;
	dynVals[1] = Defaultu;
	dynVals[2] = Defaultgamma;
	if (defaultBiophysParams)
	{
		Od   = DefaultOd;
		Of   = DefaultOf;  
		U0   = DefaultU0;  
		rhoA = DefaultrhoA;
		nv   = Defaultnv;  
		Gv   = DefaultGv;  
		Oc   = DefaultOc;  
		spillOvFract = DefaultSpillOvFract;
	}
}

//**********************************************************************
//**********************************************************************
void TMSynapse::SetToEquilibrium()
{
	dynVals[x] = Defaultx;
	dynVals[u] = Defaultu;
	dynVals[gamma] = Defaultgamma;
}

//**********************************************************************
// Treats a presynaptic spike
//**********************************************************************
void TMSynapse::PresynSpike(double)
{
	dynVals[u] += U0 * (1.0 - dynVals[u]);
	// use u(t+) and x(t-) to compute RRs
	dynVals[gamma] += dynVals[x] * dynVals[u] * rhoA * Gv * nv;
	dynVals[x] -= dynVals[u]*dynVals[x];
}

//**********************************************************************
// Update value names in ODEProblem
//**********************************************************************
void TMSynapse::SetValNamesPostfix(std::string pf)
{
	model->AddPostfixToValName(dynVals + x, pf + "_x");
	model->AddPostfixToValName(dynVals + u, pf + "_u");
	model->AddPostfixToValName(dynVals + gamma, pf + "_gamma");
}

//**********************************************************************
// Update edge values in case parameters were changed
//**********************************************************************
void TMSynapse::UpdateEdge()
{
	Initialize();
}

//**********************************************************************
//**********************************************************************
bool TMSynapse::LoadFromStream(std::ifstream & stream)
{
	stream >> defaultBiophysParams;
	stream >> Od;  
	stream >> Of;  
	stream >> U0;  
	stream >> rhoA;
	stream >> nv;  
	stream >> Gv;  
	stream >> Oc;  
	stream >> spillOvFract;
	Initialize();
	return not stream.eof();
}

//**********************************************************************
//**********************************************************************
bool TMSynapse::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< defaultBiophysParams << endl
		<< Od   << endl
		<< Of   << endl
		<< U0   << endl
		<< rhoA << endl
		<< nv   << endl
		<< Gv   << endl
		<< Oc   << endl
		<< spillOvFract << endl;
	return stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler TMSynapse::BuildModelParamHandler()
{
	ParamHandler params;
	params <= "DefaultxVal", Defaultx;
	params <= "DefaultuVal", Defaultu;
	params <= "DefaultgammaVal", Defaultgamma;
	params <= "Od", DefaultOd;
	params <= "Of", DefaultOf;
	params <= "U0", DefaultU0;
	params <= "rhoA", DefaultrhoA;   
	params <= "nv", Defaultnv;
	params <= "Gv", DefaultGv;   
	params <= "Oc", DefaultOc;   
	params <= "spillOverFract", DefaultSpillOvFract;
	return params;
}

//**********************************************************************
//**********************************************************************
double TMSynapse::GetGluVal() const
{
	return dynVals[TMSynapse::gamma];
}

//**********************************************************************
//**********************************************************************
unsigned int TMSynapse::GetNbDynVal() const
{
	return TMSYNAPSE_NBVALS_PER_SYN;
}

//********************************************************************//
//**** O P T I M   T S O D Y K S   M A R K R A M   S Y N A P S E *****//
//********************************************************************//

std::string TMSynapseOptim::ClassName = "TMSynapseOptim";

//**********************************************************************
// Default Constructor
//**********************************************************************
TMSynapseOptim::TMSynapseOptim(const NeuronNetModel * _model, double *_dv, 
		bool _fv) : TMSynapse(_model, _dv, _fv), lastSpikeTime(0),
		lastU(TMSynapse::Defaultu), lastX(TMSynapse::Defaultx)
{
	nbDynVals = TMSYNAPSE_OPTIM_NBVALS_PER_SYN;
	delete funct;
	funct = new ODE::TMSynapseFunctOptim(*this);
	if (not _dv)
	{
		delete[] dynVals;
		dynVals = new double[nbDynVals];
		freeDynVals = true;
	}
	Initialize();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
TMSynapseOptim::TMSynapseOptim(std::ifstream & stream) : TMSynapse(stream)
{
	lastSpikeTime = 0;
	lastU = TMSynapse::Defaultu;
	lastX = TMSynapse::Defaultx;
	delete funct;
	funct = new ODE::TMSynapseFunctOptim(*this);
	nbDynVals = TMSYNAPSE_OPTIM_NBVALS_PER_SYN;
	dynVals = new double[nbDynVals];
	freeDynVals = true;
}

//**********************************************************************
// Copy constructor
//**********************************************************************
TMSynapseOptim::TMSynapseOptim(const TMSynapseOptim & s) : TMSynapse(s),
	lastSpikeTime(0), lastU(s.lastU), lastX(s.lastX)
{
	delete funct;
	funct = new ODE::TMSynapseFunctOptim(*this);
	Initialize();
}

//**********************************************************************
// Destructor
//**********************************************************************
TMSynapseOptim::~TMSynapseOptim()
{
}

//**********************************************************************
//**********************************************************************
void TMSynapseOptim::Initialize()
{
	dynVals[0] = Defaultgamma;
	lastSpikeTime = 0;
	lastU = TMSynapse::Defaultu;
	lastX = TMSynapse::Defaultx;
	if (defaultBiophysParams)
	{
		Od   = DefaultOd;
		Of   = DefaultOf;  
		U0   = DefaultU0;  
		rhoA = DefaultrhoA;
		nv   = Defaultnv;  
		Gv   = DefaultGv;  
		Oc   = DefaultOc;  
		spillOvFract = DefaultSpillOvFract;
	}
}

//**********************************************************************
//**********************************************************************
void TMSynapseOptim::SetToEquilibrium()
{
	dynVals[gammaOpt] = Defaultgamma;
}

//**********************************************************************
// Treats a presynaptic spike
//**********************************************************************
void TMSynapseOptim::PresynSpike(double t)
{
	lastX = 1.0 - (1.0 - lastX) * gsl_sf_exp(-Od * (t - lastSpikeTime));
	lastU = lastU * gsl_sf_exp(-Of * (t - lastSpikeTime));

	lastU += U0 * (1.0 - lastU);

	// Use u(t+) and x(t-) to compute RRs
	dynVals[gammaOpt] += lastX * lastU * rhoA * Gv * nv;

	lastX -= lastU * lastX;

	lastSpikeTime = t;
}

//**********************************************************************
// Update value names in ODEProblem
//**********************************************************************
void TMSynapseOptim::SetValNamesPostfix(std::string pf)
{
	model->AddPostfixToValName(dynVals + gammaOpt, pf + "_gammaOpt");
}

//**********************************************************************
//**********************************************************************
double TMSynapseOptim::GetGluVal() const
{
	return dynVals[TMSynapseOptim::gammaOpt];
}

//**********************************************************************
//**********************************************************************
unsigned int TMSynapseOptim::GetNbDynVal() const
{
	return TMSYNAPSE_OPTIM_NBVALS_PER_SYN;
}

