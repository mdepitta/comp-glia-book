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

#include "ChICell.h"
#include "ChIModel.h"

using namespace AstroModel;
using namespace std;

//********************************************************************//
//*********************** C H I   C E L L ****************************//
//********************************************************************//

double ChICell::DefaultCa  = 0.0351e-03;
double ChICell::Defaulth   = 0.9122;
double ChICell::DefaultIP3 = 0.3046e-03;

double ChICell::InitCaVarRatio = 0.0;
double ChICell::InithVarRatio = 0.0;
double ChICell::InitIP3VarRatio = 0.0;
double ChICell::Defaultc1 = 0.185;
double ChICell::DefaultrC = 6.0;
double ChICell::DefaultrL = 0.11;
double ChICell::DefaultvER = 0.9e-03;
double ChICell::Defaulta2 = 0.2e03;
double ChICell::DefaultKer = 0.05e-03;  
double ChICell::Defaultvd = 0.7e-03;   
double ChICell::DefaultKplcd = 0.1e-03;
double ChICell::Defaultkd = 1.5e-03;   
double ChICell::Defaultk3 = 1.0e-03;   
double ChICell::Defaultvbeta = 0.2e-03;
double ChICell::DefaultkR = 1.3e-03;
double ChICell::DefaultkP = 10e-03;
double ChICell::Defaultkpi = 0.6e-03;

std::string ChICell::ClassName = "ChICell";
unsigned int ChICell::NbValsPerCell = CHIMODEL_NBVALS_PER_CELL;

//**********************************************************************
// Default Constructor
//**********************************************************************
ChICell::ChICell(const ODENetworkDynamicsModel<CouplingFunction, ChICell> * _model, 
	double *_dv, bool _fv) : model(0),
	defaultBiophysParams(true), c1(Defaultc1), rC(DefaultrC), 
	rL(DefaultrL), vER(DefaultvER), a2(Defaulta2), Ker(DefaultKer), 
	vd(Defaultvd), Kplcd(DefaultKplcd), kd(Defaultkd), k3(Defaultk3), 
	vbeta(Defaultvbeta), kR(DefaultkR), kP(DefaultkP), kpi(Defaultkpi), 
	dynVals(_dv), freeDynVals(_fv), nbDynVals(CHIMODEL_NBVALS_PER_CELL), 
	totFlux(0), caSpontLeak(false)
{
	model = static_cast<const ChIModel*>(_model);
	funct = new ODE::ChICellFunct(*this);
	if (not _dv)
	{
		dynVals = new double[nbDynVals];
		freeDynVals = true;
	}
	ChICell::Initialize();
}

//**********************************************************************
// Copy constructor
//**********************************************************************
ChICell::ChICell(const ChICell & c) : model(c.model), 
	defaultBiophysParams(c.defaultBiophysParams), c1(c.c1), 
	rC(c.rC), rL(c.rL), vER(c.vER), a2(c.a2), Ker(c.Ker), vd(c.vd), 
	Kplcd(c.Kplcd), kd(c.kd), k3(c.k3), vbeta(c.vbeta), kR(c.kR), kP(c.kP), 
	kpi(c.kpi), dynVals(c.dynVals), freeDynVals(c.freeDynVals), 
	nbDynVals(c.nbDynVals), totFlux(c.totFlux),	caSpontLeak(c.caSpontLeak)
{
	funct = new ODE::ChICellFunct(*this);
	if (freeDynVals)
	{
		dynVals = new double[nbDynVals];
		for (unsigned int i = 0 ; i < nbDynVals ; ++i)
			dynVals[i] = c.dynVals[i];
	}
	Initialize();
}

//**********************************************************************
// Destructor
//**********************************************************************
ChICell::~ChICell()
{
	delete funct;
	if (freeDynVals)
		delete[] dynVals;
}

//**********************************************************************
//**********************************************************************
void ChICell::Initialize()
{
	dynVals[0] = DefaultCa  * (1.0 + InitCaVarRatio  * (2.0 * UnifRand() - 1.0));
	dynVals[1] = Defaulth   * (1.0 + InithVarRatio   * (2.0 * UnifRand() - 1.0));
	dynVals[2] = DefaultIP3 * (1.0 + InitIP3VarRatio * (2.0 * UnifRand() - 1.0));
	totFlux = 0;
	caSpontLeak = false;
	gluIP3Prod = 0.0;
	if (defaultBiophysParams)
	{
		c1    = Defaultc1;
		rC    = DefaultrC;
		rL    = DefaultrL;
		vER   = DefaultvER;
		a2    = Defaulta2;
		Ker   = DefaultKer;
		vd    = Defaultvd;
		Kplcd = DefaultKplcd;
		kd    = Defaultkd;
		k3    = Defaultk3;
		vbeta = Defaultvbeta;
		kR    = DefaultkR;
		kP    = DefaultkP;
		kpi   = Defaultkpi;
	}
}

//**********************************************************************
//**********************************************************************
void ChICell::SetToEquilibrium()
{
	dynVals[Ca] = DefaultCa;
	dynVals[h] = Defaulth;
	dynVals[IP3] = DefaultIP3;
	totFlux = 0;
	caSpontLeak = false;
	gluIP3Prod = 0;
}

//**********************************************************************
//**********************************************************************
bool ChICell::LoadFromStream(std::ifstream & stream)
{
	stream >> c1;
	stream >> rC;   
	stream >> rL;   
	stream >> vER;  
	stream >> a2;   
	stream >> Ker;  
	stream >> vd;   
	stream >> Kplcd;
	stream >> kd;   
	stream >> k3;   
	stream >> vbeta;
	stream >> kR;
	stream >> kP;
	stream >> kpi;
	return not stream.eof();
}

//**********************************************************************
//**********************************************************************
bool ChICell::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< c1    << endl
		<< rC    << endl
		<< rL    << endl
		<< vER   << endl
		<< a2    << endl
		<< Ker   << endl
		<< vd    << endl
		<< Kplcd << endl
		<< kd    << endl
		<< k3    << endl
		<< vbeta << endl
		<< kR    << endl
		<< kP    << endl
		<< kpi   << endl;
	return stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler ChICell::BuildModelParamHandler()
{
	ParamHandler params;
	params <= "InitCaVarRatio", InitCaVarRatio;
	params <= "InithVarRatio", InithVarRatio;
	params <= "InitIP3VarRatio", InitIP3VarRatio;
	params <= "DefaultCaVal", DefaultCa;
	params <= "DefaulthVal", Defaulth;
	params <= "DefaultIP3Val", DefaultIP3;
	params <= "c1", Defaultc1;
	params <= "rC", DefaultrC;
	params <= "rL", DefaultrL;
	params <= "vER", DefaultvER;
	params <= "a2", Defaulta2;
	params <= "Ker", DefaultKer;  
	params <= "vd", Defaultvd;   
	params <= "Kplcd", DefaultKplcd;
	params <= "kd", Defaultkd;   
	params <= "k3", Defaultk3;   
	params <= "vbeta", Defaultvbeta;   
	params <= "kR", DefaultkR;   
	params <= "kP", DefaultkP;   
	params <= "kpi", Defaultkpi;   
	return params;
}

//**********************************************************************
// Change dynamic values to given pointer
//**********************************************************************
void ChICell::SetDynVals(double *_dv, bool _f)
{
	if (freeDynVals)
		delete[] dynVals;
	dynVals = _dv;
	freeDynVals = _f;
}

//**********************************************************************
// Return the number of dyn vals per cell
//**********************************************************************
unsigned int ChICell::GetNbDynVals()
{
	return nbDynVals;
}

//**********************************************************************
// Update value names in ODEProblem
//**********************************************************************
void ChICell::SetValNamesPostfix(std::string pf)
{
	model->AddPostfixToValName(dynVals + Ca, pf + "_Ca");
	model->AddPostfixToValName(dynVals + h, pf + "_h");
	model->AddPostfixToValName(dynVals + IP3, pf + "_IP3");
}

