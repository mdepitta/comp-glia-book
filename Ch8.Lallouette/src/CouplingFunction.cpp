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

#include "CouplingFunction.h"

#include "Network.h"
#include <math.h>
#include "utility.h"

using namespace AstroModel;
using namespace std;

string SigmoidCoupling::ClassName("SigmoidFunction");
string LinearCoupling::ClassName("LinearFunction");
string BooleanLink::ClassName("BooleanLink");

//********************************************************************//
//**************** C O U P L I N G   F U N C T I O N *****************//
//********************************************************************//

double CouplingFunction::DefaultCouplingStrength = 2.0;
double CouplingFunction::DefaultCouplStrengthStdDev = 0;
int CouplingFunction::DefaultCouplDistrMethod = CouplingFunction::TruncGauss;

//**********************************************************************
// Default Constructor
//**********************************************************************
CouplingFunction::CouplingFunction(double _F, double _sigF, 
	int _meth) : F(0), constStrength(false)
{
	static bool firstCreate = true;
	if (firstCreate)
	{
		DefaultCouplingStrength = 
			ParamHandler::GlobalParams.getParam<double>("-F", 0);
		DefaultCouplDistrMethod = 
			ParamHandler::GlobalParams.getParam<int>("-F", 1);
		firstCreate = false;
	}

	// If modifications here, also apply the same in UpdateEdge()
	redrawStrength(_F, _sigF, _meth);
}

//**********************************************************************
// Normalizes link strengths according to node degrees and 
// current strength values
//**********************************************************************
void CouplingFunction::normalizeLinkStrength(AbstractNetwork & network)
{
	Network<CouplingFunction> & net = dynamic_cast<Network<
		CouplingFunction>& >(network);

	std::vector<double> nodeDegrees;
	std::vector<double> nodeDegreesSquared;
	std::vector<double> linkStr;
	for (unsigned int i = 0 ; i < net.size() ; ++i)
	{
		nodeDegrees.push_back(net.GetNodeDegree(i));
		nodeDegreesSquared.push_back(pow(net.GetNodeDegree(i),2.0));
		for (unsigned int j = i + 1 ; j < net[i].size() ; ++j)
			if (net[i][j])
				linkStr.push_back(net[i][j]->GetStrength());
	}

	double meanDeg = ComputeMean(nodeDegrees);
	double meanLinkStr = ComputeMean(linkStr);

	double m2 = ComputeMean(nodeDegreesSquared);
	if (meanDeg != 0)
	{
		double rho = m2 / meanDeg;
		double maxK = ComputeMax(nodeDegrees);

		if (maxK != rho)
		{
			for (unsigned int i = 0 ; i < net.size() ; ++i)
				for (unsigned int j = 0 ; j < net[i].size() ; ++j)
				{
					if (net[i][j])
					{
						double ki = net.GetNodeDegree(i);
						double kj = net.GetNodeDegree(j);
						net[i][j]->ChangeStrength(meanLinkStr * 
							(1 + (2 * rho - (ki + kj))/
							(2.0 * (maxK - rho))));
						net[i][j]->SetConstStrength(true);
					}
				}
		}
	}
}

//**********************************************************************
// Update edge values in case parameters were changed
//**********************************************************************
void CouplingFunction::UpdateEdge()
{
	// If modifications here, also apply to the constructor
	if (not constStrength)
	{
		F = 0;
		redrawStrength(DefaultCouplingStrength, DefaultCouplStrengthStdDev, 
			DefaultCouplDistrMethod);
	}
}

//**********************************************************************
// Set Strength according to specified distribution
//**********************************************************************
void CouplingFunction::redrawStrength(double _F, double _sigF, int _meth)
{
	if (_sigF == 0)
		F = _F;
	else
	{
		if (_meth == TruncGauss)
		{
			F = GaussianRand(_F, _sigF);
			if (F < 0)
				F = 0;
		}
		else if (_meth == SymTruncGauss)
		{
			do
			{
				F = GaussianRand(_F, _sigF);
			}
			while ((F < 0) or (F > 2*_F));
		}
		else if (_meth == StickSymTruncGauss)
		{
			F = GaussianRand(_F, _sigF);
			if (F < 0)
				F = 0;
			else if (F > 2*_F)
				F = 2*_F;
		}
		else if (_meth == GammaLaw)
		{
			F = GammaRand(_F * _F / (_sigF * _sigF), _sigF * _sigF / _F);
		}
	}
}

//**********************************************************************
// Build parameter handling object, handles static default params
//**********************************************************************
ParamHandler CouplingFunction::BuildModelParamHandler()
{
	ParamHandler params;
	params <= "DefaultCouplingStrength", DefaultCouplingStrength;
	params <= "DefaultCouplingStrengthStdDev", DefaultCouplStrengthStdDev;
	params <= "DefaultCouplingDistrMethod", DefaultCouplDistrMethod;
	return params;
}

//********************************************************************//
//***************** S I G M O I D   C O U P L I N G ******************//
//********************************************************************//

double SigmoidCoupling::DefaultLinkThreshold = 0.3e-03;
double SigmoidCoupling::DefaultLnkThreshStdDev = 0;
double SigmoidCoupling::DefaultLinkScale = 0.05e-03;
double SigmoidCoupling::DefaultLnkScaleStdDev = 0;

//**********************************************************************
// Default Constructor
//**********************************************************************
SigmoidCoupling::SigmoidCoupling(double _F, double _sigF, 
	double _IP3Scale, double _IP3ScaleSD, double _IP3Thresh,
	double _IP3ThreshSD) : CouplingFunction::CouplingFunction(_F, _sigF),
		IP3Scale(0), IP3Thresh(0)
{
	if (_IP3ScaleSD == 0)
		IP3Scale = _IP3Scale;
	else
		while (IP3Scale <= 0)
			IP3Scale = GaussianRand(_IP3Scale, _IP3ScaleSD);

	if (_IP3ThreshSD == 0)
		IP3Thresh = _IP3Thresh;
	else
		while (IP3Thresh <= 0)
			IP3Thresh = GaussianRand(_IP3Thresh, _IP3ThreshSD);
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
SigmoidCoupling::SigmoidCoupling(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Sigmoid function operator
//**********************************************************************
double SigmoidCoupling::operator()(double DIP3) const
{
	return F / 2.0 * (1.0 + fast_tanh((fabs(DIP3) - IP3Thresh) / IP3Scale))
		* ((DIP3 > 0) ? 1.0 : -1.0);
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool SigmoidCoupling::LoadFromStream(std::ifstream & stream)
{
	stream >> F;
	stream >> constStrength;
	stream >> IP3Scale;
	stream >> IP3Thresh;
	return not stream.eof();
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool SigmoidCoupling::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< F         << endl
		<< constStrength << endl
		<< IP3Scale  << endl
		<< IP3Thresh << endl;
	return true;
}

//**********************************************************************
// Build parameter handling object, handles static default params
//**********************************************************************
ParamHandler SigmoidCoupling::BuildModelParamHandler()
{
	ParamHandler params;
	params += CouplingFunction::BuildModelParamHandler();
	params <= "DefaultLinkThreshold", DefaultLinkThreshold;
	params <= "DefaultLinkThresholdStdDev", DefaultLnkThreshStdDev;
	params <= "DefaultLinkScale", DefaultLinkScale;
	params <= "DefaultLinkScaleStdDev", DefaultLnkScaleStdDev;
	return params;
}

//**********************************************************************
// Update edge values in case parameters were changed
//**********************************************************************
void SigmoidCoupling::UpdateEdge()
{
	CouplingFunction::UpdateEdge();
	if (not constStrength)
	{
		IP3Scale = 0;
		if (DefaultLnkScaleStdDev == 0)
			IP3Scale = DefaultLinkScale;
		else
			while (IP3Scale <= 0)
				IP3Scale = GaussianRand(DefaultLinkScale, DefaultLnkScaleStdDev);

		IP3Thresh = 0;
		if (DefaultLnkThreshStdDev == 0)
			IP3Thresh = DefaultLinkThreshold;
		else
			while (IP3Thresh <= 0)
				IP3Thresh = GaussianRand(DefaultLinkThreshold, DefaultLnkThreshStdDev);
	}
}

//********************************************************************//
//****************** L I N E A R   C O U P L I N G *******************//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
LinearCoupling::LinearCoupling(double _F) : 
	CouplingFunction::CouplingFunction(_F)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
LinearCoupling::LinearCoupling(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Linear function operator
//**********************************************************************
double LinearCoupling::operator()(double DIP3) const
{
	return F * DIP3;
}

//**********************************************************************
// Update edge values in case parameters were changed
//**********************************************************************
void LinearCoupling::UpdateEdge()
{
	CouplingFunction::UpdateEdge();
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool LinearCoupling::LoadFromStream(std::ifstream & stream)
{
	stream >> F;
	stream >> constStrength;
	return not stream.eof();
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool LinearCoupling::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< F << endl
		<< constStrength << endl;
	return true;
}

