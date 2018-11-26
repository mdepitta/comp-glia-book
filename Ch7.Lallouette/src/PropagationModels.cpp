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

#include "PropagationModels.h"

#include "utility.h"

using namespace std;
using namespace AstroModel;

string ThresholdModel::ClassName = "ThresholdPropagationModel";
string ModifiedThresholdModel::ClassName = "ModifiedThresholdPropagationModel";
string SecondOrdNeighbSERSModel::ClassName = "SecondOrdNeighbSERSModel";
string SimpleThresholdSERSModel::ClassName = "SimpleThresholdSERSModel";

//********************************************************************//
//******************* U T I L I T Y   C L A S S E S ******************//
//********************************************************************//

static std::string boolValNames[2] = {"ON_STATE", "OFF_STATE"};
static bool boolVals[2] = {true, false};
static BooleanWrapper boolDecl(true, boolValNames, boolVals);

static std::string SERValNames[3] = {"REFRACTORY_STATE", "SUSCEPTIBLE_STATE", "EXCITED_STATE"};
static SERStatesType SERVals[3] = {Refractory, Susceptible, Excited};
static SERStates SERDecl(Susceptible, SERValNames, SERVals);

//********************************************************************//
//******************* T H R E S H O L D   M O D E L ******************//
//********************************************************************//

//**********************************************************************
// Initializes the model
//**********************************************************************
void ThresholdModel::Initialize(ResultSaver)
{
	PropagationModel<BooleanWrapper, BooleanLink>::Initialize();

	std::set<unsigned int> expInds, tempInds;

	bool stimPointWasEmpty = false;
	if (stimPoints.empty())
	{
		stimPoints.push_back(floor(UnifRand() * states.size()));
		stimPointWasEmpty = true;
	}

	std::set<unsigned int> actNodes;
	// for each stimulated node
	for (unsigned int i = 0 ; i < stimPoints.size() ; ++i)
	{
		expInds.insert(stimPoints[i]);
		// Get neighborhood at radius stimRadius
		for (int radius = stimRadius ; radius > 0 ; --radius)
		{
			for (std::set<unsigned int>::iterator it = expInds.begin() ; it != expInds.end() ; ++it)
			{
				std::set<unsigned int> neighbs(this->GetNetwork().GetNeighbors(*it).begin(), this->GetNetwork().GetNeighbors(*it).end());
				tempInds.insert(neighbs.begin(), neighbs.end());
			}
			expInds.insert(tempInds.begin(), tempInds.end());
			tempInds.clear();
		}
		// Activate resulting nodes
		for (std::set<unsigned int>::iterator it = expInds.begin() ; it != expInds.end() ; ++it)
		{
			states[*it] = true;
			actNodes.insert(*it);
		}
		expInds.clear();
	}
	nbStimNodes = actNodes.size();

	if (stimPointWasEmpty)
		stimPoints.clear();
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool ThresholdModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = PropagationModel<BooleanWrapper, BooleanLink>::LoadFromStream(stream);
	stream >> threshold;
	stream >> stimRadius;
	stimPoints.clear();
	unsigned int tmpSize = 0;
	stream >> tmpSize;
	unsigned int point = 0;
	for (unsigned int i = 0 ; i < tmpSize ; ++i)
	{
		stream >> point;
		stimPoints.push_back(point);
	}
	ThresholdModel::Initialize();
	return ok and stream.good();
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool ThresholdModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = PropagationModel<BooleanWrapper, BooleanLink>::SaveToStream(stream);
	stream 
		<< threshold << std::endl
		<< stimRadius << std::endl;
	stream << stimPoints.size() << std::endl;
	for (unsigned int i = 0 ; i < stimPoints.size() ; ++i)
		stream << stimPoints[i] << std::endl;
	return ok and stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler ThresholdModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += PropagationModel<BooleanWrapper, BooleanLink>::BuildModelParamHandler();
	params <= "ThresholdModelThresh", threshold;
	params <= "ThresholdModelStimRadius", stimRadius;
	return params;
}

//**********************************************************************
//**********************************************************************
unsigned int ThresholdModel::chooseNode() const
{
	return floor(this->network->size() * UnifRand());
}

//**********************************************************************
//**********************************************************************
BooleanWrapper ThresholdModel::getNewState(unsigned int ind) const
{
	std::set<unsigned int> neighbs(this->network->GetNeighbors(ind).begin(), this->network->GetNeighbors(ind).end());

	double nbActNeighb = 0;
	for (std::set<unsigned int>::iterator it = neighbs.begin() ; it != neighbs.end() ; ++it)
		if (this->states[*it])
			++nbActNeighb;
	return (neighbs.size() > 0) ? (((nbActNeighb / ((double) neighbs.size())) > threshold) ? BooleanWrapper(true) : this->states[ind]) : this->states[ind];
}

//**********************************************************************
//**********************************************************************
bool ThresholdModel::isFinished() const
{
	return false;
}

//********************************************************************//
//******** M O D I F I E D   T H R E S H O L D   M O D E L ***********//
//********************************************************************//

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool ModifiedThresholdModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = ThresholdModel::LoadFromStream(stream);
	stream >> propConst;
	stream >> Bval;
	ThresholdModel::Initialize();
	return ok and stream.good();
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool ModifiedThresholdModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = ThresholdModel::SaveToStream(stream);
	stream << propConst << std::endl;
	stream << Bval << std::endl;
	return ok and stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler ModifiedThresholdModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += ThresholdModel::BuildModelParamHandler();
	params <= "ModThreshPropConst", propConst;
	params <= "ModThreshBVal", Bval;
	return params;
}

//**********************************************************************
//**********************************************************************
BooleanWrapper ModifiedThresholdModel::getNewState(unsigned int ind) const
{
	if (not this->states[ind])
	{
		std::set<unsigned int> neighbs(this->network->GetNeighbors(ind).begin(), this->network->GetNeighbors(ind).end());
		std::set<unsigned int> tempNeighbs;

		double tempUnAct = 0;
		double nbActNeighb = 0;
		for (std::set<unsigned int>::iterator it = neighbs.begin() ; it != neighbs.end() ; ++it)
			if (this->states[*it])
			{
				// Count how many un activated neighbors the neighbhoring node has
				tempUnAct = 0;
				tempNeighbs = std::set<unsigned int>(this->network->GetNeighbors(*it).begin(), this->network->GetNeighbors(*it).end());
				for (std::set<unsigned int>::iterator it2 = tempNeighbs.begin() ; it2 != tempNeighbs.end() ; ++it2)
					if (not this->states[*it2])
						++tempUnAct;
				// set activating potential of current neighbor
				nbActNeighb += propConst / tempUnAct;
			}
		double k = ((double) neighbs.size());
		return (neighbs.size() > 0) ? (((nbActNeighb) >= (threshold * k + Bval)) ? BooleanWrapper(true) : this->states[ind]) : this->states[ind];
	}
	return this->states[ind];
}

//********************************************************************//
//************************ S E R S   M O D E L ***********************//
//********************************************************************//

//**********************************************************************
// Constructor
//**********************************************************************
SERSModel::SERSModel(ParamHandler & h) : 
	PropagationModel<SERStates, BooleanLink>::PropagationModel(h)
{
	period     = h.getParam<double>("-SERSModel", 0);
	actTime    = h.getParam<double>("-SERSModel", 1);
	transTime  = h.getParam<double>("-SERSModel", 2);
	stimRadius = h.getParam<double>("-SERSModel", 3);
	stimPoints = param.getParam<std::vector<int> >("-SERSModelStims", 0);
	stimOrderInd = 0;
}

//**********************************************************************
// Initializes the model
//**********************************************************************
void SERSModel::Initialize(ResultSaver)
{
	PropagationModel<SERStates, BooleanLink>::Initialize();
	std::set<unsigned int> expInds, tempInds;

	states = std::vector<SERStates>(this->GetNetwork().size(), Susceptible);

	bool stimPointWasEmpty = false;
	if (stimPoints.empty())
	{
		stimPoints.push_back(floor(UnifRand() * states.size()));
		stimPointWasEmpty = true;
	}

	allStimNodes.clear();
	// for each stimulated node
	for (unsigned int i = 0 ; i < stimPoints.size() ; ++i)
	{
		expInds.insert(stimPoints[i]);
		// Get neighborhood at radius stimRadius
		for (int radius = stimRadius ; radius > 0 ; --radius)
		{
			for (std::set<unsigned int>::iterator it = expInds.begin() ; it != expInds.end() ; ++it)
			{
				std::set<unsigned int> neighbs(this->GetNetwork().GetNeighbors(*it).begin(), this->GetNetwork().GetNeighbors(*it).end());
				tempInds.insert(neighbs.begin(), neighbs.end());
			}
			expInds.insert(tempInds.begin(), tempInds.end());
			tempInds.clear();
		}
		// add nodes to all stimulated nodes
		for (std::set<unsigned int>::iterator it = expInds.begin() ; it != expInds.end() ; ++it)
		{
			this->states[*it] = SERStates(Excited);
			allStimNodes.insert(*it);
		}
		expInds.clear();
	}

	if (stimPointWasEmpty)
		stimPoints.clear();

	stimOrder = std::vector<unsigned int>(this->network->size(), 0);
	for (unsigned int i = 0 ; i < stimOrder.size() ; ++i)
		stimOrder[i] = i;
	stimOrderInd = 0;
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool SERSModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = PropagationModel<SERStates, BooleanLink>::LoadFromStream(stream);
	stream >> period;   
	stream >> actTime;  
	stream >> transTime;
	stream >> stimRadius;
	stimPoints.clear();
	unsigned int tmpSize = 0;
	stream >> tmpSize;
	unsigned int point = 0;
	for (unsigned int i = 0 ; i < tmpSize ; ++i)
	{
		stream >> point;
		stimPoints.push_back(point);
	}
	SERSModel::Initialize();
	return ok and stream.good();
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool SERSModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = PropagationModel<SERStates, BooleanLink>::SaveToStream(stream);
	stream 
		<< period << std::endl
		<< actTime << std::endl
		<< transTime << std::endl
		<< stimRadius << std::endl;
	stream << stimPoints.size() << std::endl;
	for (unsigned int i = 0 ; i < stimPoints.size() ; ++i)
		stream << stimPoints[i] << std::endl;
	return ok and stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler SERSModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += PropagationModel<SERStates, BooleanLink>::BuildModelParamHandler();
	params <= "SERSModelPeriod", period;
	params <= "SERSModelActTime", actTime;
	params <= "SERSModelTransTime", transTime;
	params <= "SERSModelStimRadius", stimRadius;
	return params;
}

//**********************************************************************
//**********************************************************************
unsigned int SERSModel::chooseNode() const
{
	// Randomize the update sequence
	if (stimOrderInd == 0)
	{
		unsigned int tmpInd;
		unsigned int tmpVal;
		for (unsigned int i = 0 ; i < stimOrder.size() ; ++i)
		{
			tmpInd = floor(this->network->size() * UnifRand());
			tmpVal = stimOrder[tmpInd];
			stimOrder[tmpInd] = stimOrder[i];
			stimOrder[i] = tmpVal;
		}
	}
	unsigned int chosenNode = stimOrder[stimOrderInd];
	stimOrderInd = (stimOrderInd + 1) % stimOrder.size();
	return chosenNode;
}

//**********************************************************************
//**********************************************************************
bool SERSModel::isFinished() const
{
	return false;
}

//**********************************************************************
//**********************************************************************
SERStates SERSModel::getNewState(unsigned int ind) const
{
	if (allStimNodes.find(ind) == allStimNodes.end())
	{
		if (this->states[ind] == Susceptible)
		{
			return TrueWithProba(this->probaSusceptibleToExcited(ind)) ?
				SERStates(Excited) : this->states[ind];
		}
		else if (this->states[ind] == Excited)
			return TrueWithProba(this->probaExcitedToRefractory(ind)) ? 
				SERStates(Refractory) : this->states[ind];
		else if (this->states[ind] == Refractory)
			return TrueWithProba(this->probaRefractoryToSusceptible(ind)) ? 
				SERStates(Susceptible) : this->states[ind];
		else
			return this->states[ind];
	}
	else
		return SERStates(Excited);
}

//**********************************************************************
//**********************************************************************
double SERSModel::probaSusceptibleToExcited(unsigned int ) const
{
	return 1.0 / transTime;
}

//**********************************************************************
//**********************************************************************
double SERSModel::probaExcitedToRefractory(unsigned int ) const
{
	return 1.0 / actTime;
}

//**********************************************************************
//**********************************************************************
double SERSModel::probaRefractoryToSusceptible(unsigned int ) const
{
	return 1.0 / (period - actTime);
}

//********************************************************************//
//*********** S E C O N D   O R D E R   S E R S   M O D E L **********//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
SecondOrdNeighbSERSModel::SecondOrdNeighbSERSModel(ParamHandler & h) : 
	SERSModel::SERSModel(h)
{
	A = h.getParam<double>("-SecOrdSERSModel", 0);
	B = h.getParam<double>("-SecOrdSERSModel", 1);
	useModTransTime = h.getParam<bool>("-SecOrdSERSModelModTransT", 0);
	minTransTimeFact = h.getParam<double>("-SecOrdSERSModelModTransT", 1);
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool SecondOrdNeighbSERSModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = SERSModel::LoadFromStream(stream);
	stream >> A;        
	stream >> B;
	stream >> useModTransTime;
	stream >> minTransTimeFact;
	return ok and stream.good();
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool SecondOrdNeighbSERSModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = SERSModel::SaveToStream(stream);
	stream 
		<< A << std::endl
		<< B << std::endl
		<< useModTransTime << std::endl
		<< minTransTimeFact << std::endl;
	return ok and stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler SecondOrdNeighbSERSModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += SERSModel::BuildModelParamHandler();
	params <= "SecOrdSERSModelA", A;
	params <= "SecOrdSERSModelB", B;
	params <= "SecOrdSERSModelModTransT", useModTransTime;
	params <= "SecOrdSERSModelModTransTFact", minTransTimeFact;
	return params;
}

//**********************************************************************
//**********************************************************************
double SecondOrdNeighbSERSModel::probaSusceptibleToExcited(
	unsigned int ind) const
{
	std::set<unsigned int> neighbs(this->network->GetNeighbors(ind).begin(), this->network->GetNeighbors(ind).end());
	std::set<unsigned int> tempNeighbs;
	double gamma = 0;

	double tempUnAct = 0;
	double sumAlphas = 0;
	for (std::set<unsigned int>::iterator it = neighbs.begin() ; 
			it != neighbs.end() ; ++it)
		if (this->states[*it] == Excited)
		{
			// Count how many unactivated neighbors the neighbhoring node has
			tempUnAct = 0;
			tempNeighbs = std::set<unsigned int>(this->network->GetNeighbors(*it).begin(), this->network->GetNeighbors(*it).end());
			for (std::set<unsigned int>::iterator it2 = tempNeighbs.begin() ;
					it2 != tempNeighbs.end() ; ++it2)
				if (this->states[*it2] != Excited)
					++tempUnAct;
			// add activating potential of current neighbor to total
			// activating potential
			sumAlphas += 1.0 / tempUnAct;
		}
	double k = ((double) neighbs.size());
	double chi = A * k + B;

	if (useModTransTime)
	{
		double tmax = transTime;
		double tmin = minTransTimeFact * transTime;
		double specifTransTime = (sumAlphas - chi) / (k - chi) * 
			(tmin - tmax) + tmax;
		
		gamma = (sumAlphas >= chi) ? 1.0 / specifTransTime : 0;
	}
	else
		gamma = (sumAlphas >= chi) ? 1.0 / transTime : 0;

	return gamma;
}

//********************************************************************//
//******* S I M P L E   T H R E S H O L D   S E R S   M O D E L ******//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
SimpleThresholdSERSModel::SimpleThresholdSERSModel(ParamHandler & h) : 
	SERSModel::SERSModel(h)
{
	relativeThreshold = h.getParam<double>("-SimpleThreshSERSModel", 0);
	spontaneousFiring = h.getParam<double>("-SimpleThreshSERSModel", 1);
	recoveryProba = h.getParam<double>("-SimpleThreshSERSModel", 2);
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool SimpleThresholdSERSModel::LoadFromStream(std::ifstream & stream)
{
	bool ok = SERSModel::LoadFromStream(stream);
	stream >> relativeThreshold;
	stream >> spontaneousFiring;
	stream >> recoveryProba;
	return ok and stream.good();
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool SimpleThresholdSERSModel::SaveToStream(std::ofstream & stream) const
{
	bool ok = SERSModel::SaveToStream(stream);
	stream 
		<< relativeThreshold << std::endl
		<< spontaneousFiring << std::endl
		<< recoveryProba << std::endl;
	return ok and stream.good();
}

//**********************************************************************
//**********************************************************************
ParamHandler SimpleThresholdSERSModel::BuildModelParamHandler()
{
	ParamHandler params;
	params += SERSModel::BuildModelParamHandler();
	params <= "SimpleThreshSERSModelRelatThresh", relativeThreshold;
	params <= "SimpleThreshSERSModelSpontFireRate", spontaneousFiring;
	params <= "SimpleThreshSERSModelRecovProba", recoveryProba;
	return params;
}

//**********************************************************************
//**********************************************************************
double SimpleThresholdSERSModel::probaSusceptibleToExcited(
	unsigned int ind) const
{
	std::set<unsigned int> neighbs(this->network->GetNeighbors(ind).begin(), this->network->GetNeighbors(ind).end());
	std::set<unsigned int> tempNeighbs;

	// Counting the number of excited neighbors
	double nbActNeighbs = 0;
	for (std::set<unsigned int>::iterator it = neighbs.begin() ; 
			it != neighbs.end() ; ++it)
		if (this->states[*it] == Excited)
			++nbActNeighbs;

	double k = ((double) neighbs.size());

	return (nbActNeighbs / k >= relativeThreshold) ? 1.0 : spontaneousFiring;
}

//**********************************************************************
//**********************************************************************
double SimpleThresholdSERSModel::probaExcitedToRefractory(
	unsigned int ) const
{
	return 1.0;
}

//**********************************************************************
//**********************************************************************
double SimpleThresholdSERSModel::probaRefractoryToSusceptible(
	unsigned int ) const
{
	return recoveryProba;
}
