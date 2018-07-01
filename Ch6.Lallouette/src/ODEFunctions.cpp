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

#include "ODEFunctions.h"
#include "ChICell.h"
#include "ChIModel.h"
#include "Synapse.h"
#include "Neuron.h"
#include "NeuronNetModels.h"
#include "AstroNeuroModel.h"

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include <math.h>

using namespace ODE;
using namespace AstroModel;
using namespace std;

std::string ChICellFunct::ClassName = "ChICellFunct";
std::string TMSynapseFunct::ClassName = "TMSynapseFunct";
std::string TMSynapseFunctOptim::ClassName = "TMSynapseFunctOptim";
std::string SFALIFNeuronFunct::ClassName = "SFALIFNeuronFunct";
std::string DummyNeuronFunct::ClassName = "DummyNeuronFunct";
std::string ChINetworkFunct::ClassName = "ChINetworkFunct";
std::string NeuronNetworkFunc::ClassName = "NeuronNetworkFunc";
std::string AstroNeuronNetFunc::ClassName = "AstroNeuronNetFunc";

//********************************************************************//
//********************** C H I  M O D E L ****************************//
//********************************************************************//


//**********************************************************************
// Constructor
//**********************************************************************
ChICellFunct::ChICellFunct(AstroModel::ChICell & _cell) : cell(_cell), model(*_cell.model)
{

}

//**********************************************************************
// Constructor
//**********************************************************************
ChINetworkFunct::ChINetworkFunct(AstroModel::ChIModel & _model) : model(_model)
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void ChICellFunct::CompFunc(const double &, const double *v, double *f) const
{
	static double Ca2;
	static double Ca4;
	static double Q2;
	static double hInf; 
	static double tauH; 
	static double mInf; 
	static double nInf; 
	static double chanProb;
	static double Jchan;
	static double Jleak;
	static double Jpump;
	static double Jspont;
	static double Pplcd;
	static double D5p;  
	static double D3k;  
	static double dIP3i;

	double & diffCa    = f[0];
	double & diffh     = f[1];
	double & diffIP3   = f[2];
	const double & Ca  = v[0];
	const double & h   = v[1];
	const double & IP3 = v[2];

	Q2     = model.d2 * (IP3 + model.d1) / (IP3 + model.d3);

	hInf   = Q2  / (Q2 + Ca);
	tauH   = 1.0 / (cell.a2 * (Q2 + Ca));
	mInf   = IP3  / (IP3 + model.d1);
	nInf   = Ca  / (Ca + model.d5);

	Ca2 = INTPOW2(Ca);
	Ca4 = INTPOW2(Ca2);
	chanProb = mInf * nInf * h;
	Jchan  = cell.rC  * (model.C0 - (1.0 + cell.c1) * Ca) * INTPOW3(chanProb);
	Jleak  = cell.rL  * (model.C0 - (1.0 + cell.c1) * Ca);
	Jpump  = cell.vER * Ca2 / (INTPOW2(cell.Ker) + Ca2);
	// Modified version to match Osama's formula
	Jspont = cell.caSpontLeak ? cell.rL : 0;

	Pplcd  = cell.vd * cell.kd / (cell.kd + IP3) * Ca2 / (Ca2 + INTPOW2(cell.Kplcd));
	D5p    = model.r5p * IP3;
	D3k    = model.v3k * Ca4 / (Ca4 + INTPOW4(model.K3k)) * IP3 / (IP3 + cell.k3);

	dIP3i  = Pplcd - D3k - D5p;

	diffCa  = Jchan + Jleak - Jpump + Jspont;
	diffh   = (hInf - h) / tauH;
	diffIP3 = dIP3i - cell.totFlux + cell.gluIP3Prod;
}

//**********************************************************************
// Function network operator, it returns the derivatives of every dynVal
//**********************************************************************
void ChINetworkFunct::CompFunc(const double & t, const double *v, double *f) const
{
	static unsigned int dynValDim = model.GetNbDynVal(0);

	model.ComputeFluxes(t);
	model.Stimulate(t);

	for (unsigned int i = 0 ; i < model.cells.size() ; ++i)
		model.cells[i]->funct->CompFunc(t, v + i * dynValDim, f + i * dynValDim);
}

//********************************************************************//
//********* T S O D Y K S   M A R K R A M   S Y N A P S E S **********//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
TMSynapseFunct::TMSynapseFunct(AstroModel::TMSynapse & _syn) : 
	synapse(_syn), model(*_syn.GetModel())
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void TMSynapseFunct::CompFunc(const double & , const double *v, 
	double *f) const
{
	double & diffx       = f[0];
	double & diffu       = f[1];
	double & diffgamma   = f[2];
	const double & x     = v[0];
	const double & u     = v[1];
	const double & gamma = v[2];

	diffx     = synapse.Od * (1.0 - x);
	diffu     = synapse.Of * (- u);
	diffgamma = -synapse.Oc * gamma;
}

//********************************************************************//
//**** O P T I M   T S O D Y K S   M A R K R A M   S Y N A P S E S ***//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
TMSynapseFunctOptim::TMSynapseFunctOptim(AstroModel::TMSynapseOptim & _syn) : 
	TMSynapseFunct(_syn)
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void TMSynapseFunctOptim::CompFunc(const double & , const double *v, 
	double *f) const
{
	double & diffgamma   = f[0];
	const double & gamma = v[0];

	diffgamma = -synapse.Oc * gamma;
}

//********************************************************************//
//************** S F A L I F   N E U R O N S   M O D E L *************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
SFALIFNeuronFunct::SFALIFNeuronFunct(AstroModel::SFALIFNeuron & _neur) : 
	neuron(_neur), model(*_neur.GetModel())
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void SFALIFNeuronFunct::CompFunc(const double & , const double *v, 
	double *f) const
{
	static double leak;
	static double tmpAMPA;

	double & diffV   = f[0];
	double & diffw   = f[1];
	const double & V = v[0];
	const double & w = v[1];

	leak = neuron.gL / neuron.C * (V - neuron.E0);
	tmpAMPA = 0.0;
	for (unsigned int i = 0 ; i < neuron.dendrSyn.size() ; ++i)
		tmpAMPA += neuron.dendrSyn[i]->GetDynVal(2) * 0.1;
	tmpAMPA = tmpAMPA / neuron.C;

	diffV = - leak - w / neuron.C + neuron.inputCurr / neuron.C;
	diffw = - w / neuron.tauw;
}

//********************************************************************//
//**************** D U M M Y   N E U R O N S   M O D E L *************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
DummyNeuronFunct::DummyNeuronFunct(AstroModel::DummyNeuron & )
{

}

//**********************************************************************
// Function cell operator, it returns the derivative of the dynamical
// values of a cell.
//**********************************************************************
void DummyNeuronFunct::CompFunc(const double & , const double *, 
	double *f) const
{
	// Do nothing you dummy !
	f[0] = 0.0;
}

//********************************************************************//
//************* N E U R O N   N E T W O R K  M O D E L ***************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
NeuronNetworkFunc::NeuronNetworkFunc(AstroModel::NeuronNetModel & _model) : model(_model)
{

}

//**********************************************************************
// Function network operator, it returns the derivatives of every dynVal
//**********************************************************************
void NeuronNetworkFunc::CompFunc(const double & t, const double *v, double *f) const
{
	unsigned int nbValDyn = 0;
	for (unsigned int i = 0 ; i < model.neurons.size() ; ++i)
	{
		model.neurons[i]->funct->CompFunc(t, v + nbValDyn, f + nbValDyn);
		nbValDyn += model.neurons[i]->GetNbDynVal();
	}

	for (unsigned int i = 0 ; i < model.synapses.size() ; ++i)
	{
		model.synapses[i]->funct->CompFunc(t, v + nbValDyn, f + nbValDyn);
		nbValDyn += model.synapses[i]->GetNbDynVal();
	}
}

//**********************************************************************
// Constructor
//**********************************************************************
AstroNeuronNetFunc::AstroNeuronNetFunc(AstroModel::AstroNeuroNetModel & _model) : model(_model)
{

}

//**********************************************************************
// Function network operator, it returns the derivatives of every dynVal
//**********************************************************************
void AstroNeuronNetFunc::CompFunc(const double & t, const double *v, double *f) const
{
	model.neuronNet->function->CompFunc(t, v, f);
	unsigned int nbDynVal = model.neuronNet->GetTotNbDynVals();

	double DefaultSpillOvFract = 0.025;
	GlutamatergicSynapse *gluSyn = 0;
	TMSynapse *tmSyn = 0;
	double glu = 0;
	double Calc = 0.0;
	for (unsigned int i = 0 ; i < model.astrToSyn.size() ; ++i)
	{
		model.astroNet->cells[i]->gluIP3Prod = 0;
		Calc = v[model.astroNet->cells[i]->dynVals - model.vals + ChICell::Ca];
		for (unsigned int j = 0 ; j < model.astrToSyn[i].size() ; ++j)
		{
			gluSyn = dynamic_cast<GlutamatergicSynapse *>(model.astrToSyn[i][j]);
			if (gluSyn)
			{
				tmSyn = dynamic_cast<TMSynapse *>(gluSyn);
				glu = (tmSyn ? tmSyn->spillOvFract : DefaultSpillOvFract) * gluSyn->GetGluVal();
			}	
			model.astroNet->cells[i]->gluIP3Prod += 
				model.astroNet->cells[i]->vbeta * (pow(glu, 0.7) / 
				(pow(glu, 0.7) + pow(model.astroNet->cells[i]->kR + 
					model.astroNet->cells[i]->kP*(Calc / (Calc + model.astroNet->cells[i]->kpi)), 0.7)));
		}
	}

	model.astroNet->function->CompFunc(t, v + nbDynVal, f + nbDynVal);
}

