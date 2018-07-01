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

#include "GridSearchSimulation.h"
#include "AbstractFactory.h"
#include "ErrorCodes.h"

using namespace AstroModel;
using namespace std;

//********************************************************************//
//************ G R I D   S E A R C H   S I M U L A T I O N ***********//
//********************************************************************//

string GridSearchSimulation::ClassName("GridSearchSimulation");

//**********************************************************************
// Default Constructor
//**********************************************************************
GridSearchSimulation::GridSearchSimulation(SimulationManager *_ss, bool _f,
	ParamHandler & _h, std::string path) : 
	subSimulation(_ss), freeSubSimulation(_f), handler(_h), saver(path)
{
	std::vector<std::string> condReplPath = 
		handler.getParam<std::vector<std::string> >("-GridCondReplPath", 0);
	// Load cond replace
	for (unsigned int i = 0 ; i < condReplPath.size() ; ++i)
		LoadCondReplaceFromFile(condReplPath[i]);
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
GridSearchSimulation::GridSearchSimulation(std::ifstream & stream,
	ParamHandler & _h) :
	subSimulation(0), freeSubSimulation(false), handler(_h)
{
	stringstream path;
	std::string mainPath = handler.getParam<std::string>("-Path");
	bool useSubDir = handler.getParam<bool>("-SubDir", 0);
	std::string subDirPath = handler.getParam<std::string>("-SubDir", 1);

	path << mainPath << "/data";
	if (useSubDir)
		path << "/" << subDirPath;
	saver = ResultSaver(path.str());

	LoadFromStream(stream);
}

//**********************************************************************
// Default Destructor
//**********************************************************************
GridSearchSimulation::~GridSearchSimulation()
{
	if (freeSubSimulation)
		delete subSimulation;
}

//**********************************************************************
// Only calls LaunchSimulation with a suitable saver
//**********************************************************************
int GridSearchSimulation::LaunchSimulation(std::string simName)
{
	ResultSaver localSaver = saver(simName);
	return LaunchSimulation(localSaver);
}

//**********************************************************************
// Launch a grid search simulation with registered parameters
//**********************************************************************
int GridSearchSimulation::LaunchSimulation(ResultSaver localSaver)
{
	int returnVal = 0;
	ParamHandler params = subSimulation->BuildParamHandler();
	// Setting global params
	for (std::map<std::string, std::string>::iterator it = 
			globalParameters.begin() ; it != globalParameters.end() ; ++it)
		for (unsigned int i = 0 ; i < params.GetNbParamsForName(it->first) 
				; ++i)
			params.SetVal(it->first, it->second.c_str(), i);


	currIndices = vector<unsigned int>(gridParameters.size(), 0);
	currCondInd.clear();
	bool firstTime = true;
	bool paramCombModified = true;
	while (firstTime or ((not currIndices.empty()) and 
		(currIndices[0] < gridParameters[0].second.size())))
	{
		// Update currCondGrid and currCondInd when the parameter 
		// combinaison in gridParameters has changed.
		if (paramCombModified)
		{
			currCondGrid.clear();
			std::map<std::pair<std::string, std::string>, 
				std::vector<std::pair<std::string, 
				std::vector<std::string> > > >::iterator it;
			// Update currCondGrid first
			// Fir each grid parameter
			for (unsigned int i = 0 ; i < gridParameters.size() ; ++i)
			{
				// If the current value "triggers" a conditional 
				// grid value
				if ((it = condGridParameters.find(std::make_pair(
					gridParameters[i].first, 
					gridParameters[i].second[currIndices[i]]))) != 
					condGridParameters.end())
				{
					// For each attached values
					for (unsigned int j = 0 ; j < it->second.size() ; ++j)
					{
						// Add the values to currCondGrid
						unsigned int tempInd = 
							getParamInd(it->second[j].first, currCondGrid);
						if (tempInd == currCondGrid.size())
							currCondGrid.push_back(
								std::make_pair(it->second[j].first, 
								it->second[j].second));
						else
							currCondGrid[tempInd].second.insert(
								currCondGrid[tempInd].second.end(), 
								it->second[j].second.begin(), 
								it->second[j].second.end());
					}
				}
			}
			// Update currCondInd
			currCondInd = std::vector<unsigned int>(currCondGrid.size(), 0);
			paramCombModified = false;
		}

		string savPath = updateGridValuesAndGetPath();
		updateGridValuesAndGetPath();
	
//==============================//
	TRACE_UP("Launching sub simulation with new param combinaison.")
	for (unsigned int i = 0 ; i < currIndices.size() ; ++i)
		TRACE(gridParameters[i].first << " : " 
			<< gridParameters[i].second[currIndices[i]])
	for (unsigned int i = 0 ; i < currCondInd.size() ; ++i)
		TRACE(currCondGrid[i].first << " : " 
			<< currCondGrid[i].second[currCondInd[i]])
//==============================//

		returnVal |= subSimulation->LaunchSimulation(localSaver(savPath));

//==============================//
	TRACE_DOWN("Sub simulation ended.")
	TRACE_UP("Computing metrics...")
//==============================//

		if (not metrics.ComputeMetricsDefault(*this))
			returnVal |= GRIDSEARCH_METRIC_COMPUTATION_PROBLEM;

//==============================//
	TRACE_DOWN("Metrics successfully computed.")
//==============================//

		// Vary currCondInd first
		if (not currCondInd.empty())
		{
			++currCondInd.back();
			for (unsigned int i = currCondInd.size() - 1 ; i > 0 ; --i)
				if (currCondInd[i] >= currCondGrid[i].second.size())
				{
					currCondInd[i] = 0;
					++currCondInd[i - 1];
				}
		}
		// If the parameters of currCondInd have been fully explored
		if (currCondInd.empty() or 
			(currCondInd[0] >= currCondGrid[0].second.size()))
		{
			if (not currIndices.empty())
			{
				++currIndices.back();
				for (unsigned int i = currIndices.size() - 1 ; i > 0 ; --i)
					if (currIndices[i] >= gridParameters[i].second.size())
					{
						currIndices[i] = 0;
						++currIndices[i - 1];
					}
			}
			paramCombModified = true;
		}

		firstTime = false;
	}

//==============================//
	TRACE_UP("End of grid search, saving metrics.")
//==============================//

	if (not metrics.SaveMetricsDefault(localSaver))
		returnVal |= GRIDSEARCH_METRIC_SAVING_PROBLEM;

//==============================//
	TRACE_DOWN("Grid search metrics saved.")
	TRACE_UP("Compressing the results.")
//==============================//


//==============================//
	TRACE_DOWN("Results successfully compressed.")
//==============================//

	return returnVal;
}


//**********************************************************************
// Returns its own metrics and the metric tree of all its "children" objects
//**********************************************************************
std::vector<Metric *> GridSearchSimulation::GetAllMetrics() const
{
	std::vector<Metric *> temp = metrics.GetMetricsRaw();
	if (subSimulation)
	{
		std::vector<Metric *> tempSubSim = subSimulation->GetAllMetrics();
		temp.insert(temp.begin(), tempSubSim.begin(), tempSubSim.end());
	}
	return temp;
}

//**********************************************************************
// Cut a large grid search into smaller subGridsearches in order to save them
//**********************************************************************
vector<GridSearchSimulation *> GridSearchSimulation::
	CreateSmallerGridSearches(unsigned int desiredSize, int tmp) const
{
	vector<GridSearchSimulation*> smallerGrids;
	unsigned int totalGridSize = GetFullGridSize();

	if (desiredSize < totalGridSize)
	{
		std::pair<GridSearchSimulation*, GridSearchSimulation*> subs = 
			separateinTwo();

		std::vector<GridSearchSimulation *> sub1 = 
			subs.first->CreateSmallerGridSearches(desiredSize, tmp + 1);
		std::vector<GridSearchSimulation *> sub2 = 
			subs.second->CreateSmallerGridSearches(desiredSize, tmp + 1);

		smallerGrids.insert(smallerGrids.end(), sub1.begin(), sub1.end());
		smallerGrids.insert(smallerGrids.end(), sub2.begin(), sub2.end());

		delete subs.first;
		delete subs.second;
	}
	else
		smallerGrids.push_back(getCopy());

	return smallerGrids;
}

bool GridSearchSimulation::ComputeAndSaveGridSearches(unsigned int nb, 
	std::string scriptPath, std::string gridPath, std::string )
{
IF_DEBUG("Splitting Search Grid")
	ofstream scriptStream(scriptPath.c_str());
	vector<GridSearchSimulation*> smallGrids = CreateSmallerGridSearches(nb);
	for (unsigned int i = 0 ; i < smallGrids.size() ; ++i)
	{
		string smallGridPath = gridPath + StringifyFixed(i) + ".gsim";
		ofstream stream(smallGridPath.c_str());
		if (not smallGrids[i]->SaveToStream(stream))
		{
			std::cerr << "Couldn't save" << smallGridPath << std::endl;
			return false;
		}
		scriptStream << "-Load " << smallGridPath << endl;

		delete smallGrids[i];
	}
	return true;
}

//**********************************************************************
// Add a file to be saved
//**********************************************************************
void GridSearchSimulation::AddFileToSave(std::string name)
{
	saver <= name;
}

//**********************************************************************
// Load conditional replace from file (saved by GridSearchMetric)
//**********************************************************************
bool GridSearchSimulation::LoadCondReplaceFromFile(std::string path)
{
	std::ifstream stream(path.c_str());
	unsigned int befNamesSize, aftNamesSize;
	std::string tmpStr;

	assert(stream.good());
	int oldBefInd = condReplaceBeforeNames.size() - 1;
	int oldAftInd = condReplaceAfterNames.size() - 1;
	bool found;
	int tmpInd;

	stream >> befNamesSize;
	std::vector<unsigned int> befInds(befNamesSize, 0);
	for (unsigned int i = 0 ; i < befNamesSize ; ++i)
	{
		stream >> tmpStr;
		found = false;
		for (tmpInd = 0 ; not found and (tmpInd <= oldBefInd) ; ++tmpInd)
			found |= (condReplaceBeforeNames[tmpInd] == tmpStr);
		if (not found)
		{
			condReplaceBeforeNames.push_back(tmpStr);
			befInds[i] = condReplaceBeforeNames.size() - 1;
		}
		else
			befInds[i] = tmpInd - 1;
	}
	stream >> aftNamesSize;
	std::vector<unsigned int> aftInds(aftNamesSize, 0);
	for (unsigned int i = 0 ; i < aftNamesSize ; ++i)
	{
		stream >> tmpStr;
		found = false;
		for (tmpInd = 0 ; not found and (tmpInd <= oldAftInd) ; ++tmpInd)
			found |= (condReplaceAfterNames[tmpInd] == tmpStr);
		if (not found)
		{
			condReplaceAfterNames.push_back(tmpStr);
			aftInds[i] = condReplaceAfterNames.size() - 1;
		}
		else
			aftInds[i] = tmpInd - 1;
	}

	// Update old condReplaceVals
	for (unsigned int i = 0 ; i < condReplaceVals.size() ; ++i)
	{
		for (unsigned int j = oldBefInd + 1 ; j < condReplaceBeforeNames.size() ; ++j)
			condReplaceVals[i].first.push_back("NULL");
		for (unsigned int j = oldAftInd + 1 ; j < condReplaceAfterNames.size() ; ++j)
			condReplaceVals[i].second.push_back("NULL");
	}

	unsigned int nbConds;
	stream >> nbConds;
	for (unsigned int i = 0 ; i < nbConds ; ++i)
	{
		condReplaceVals.push_back(std::make_pair(
			std::vector<std::string>(condReplaceBeforeNames.size(), "NULL"), 
			std::vector<std::string>(condReplaceAfterNames.size(), "NULL")));
		for (unsigned int j = 0 ; j < befNamesSize ; ++j)
		{
			stream >> tmpStr;
			condReplaceVals.back().first[befInds[j]] = tmpStr;
		}
		for (unsigned int j = 0 ; j < aftNamesSize ; ++j)
		{
			stream >> tmpStr;
			condReplaceVals.back().second[aftInds[j]] = tmpStr;
		}
	}

	return stream.good();
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool GridSearchSimulation::LoadFromStream(std::ifstream & stream)
{
	bool ok = true;

	totalParamNames.clear();
	// grid parameters
	gridParameters.clear();
	unsigned int gridSize, paramSize;
	string paramName, paramVal;
	stream >> gridSize;
	for (unsigned int i = 0 ; i < gridSize ; ++i)
	{
		stream >> paramName;
		stream >> paramSize;
		if (std::find(totalParamNames.begin(), totalParamNames.end(), paramName) ==
				totalParamNames.end())
			totalParamNames.push_back(paramName);
		gridParameters.push_back(make_pair(paramName, vector<string>()));
		for (unsigned int j = 0 ; j < paramSize ; ++j)
		{
			stream >> paramVal;
			gridParameters.back().second.push_back(paramVal);
		}
	}
	// conditional grid parameters loading
	condGridParameters.clear();
	unsigned int condGridSize, condParamGroupSize, condParamSize;
	string str1, str2;
	std::pair<std::string, std::string> key;
	stream >> condGridSize;
	for (unsigned int i = 0 ; i < condGridSize ; ++i)
	{
		stream >> str1;
		stream >> str2;
		key = make_pair(str1, str2);
		stream >> condParamGroupSize;
		for (unsigned int j = 0 ; j < condParamGroupSize ; ++j)
		{
			stream >> str1;
			if (std::find(totalParamNames.begin(), totalParamNames.end(), str1) ==
					totalParamNames.end())
				totalParamNames.push_back(str1);
			condGridParameters[key].push_back(make_pair(str1, std::vector<std::string>()));
			stream >> condParamSize;
			for (unsigned int k = 0 ; k < condParamSize ; ++k)
			{
				stream >> str1;
				condGridParameters[key].back().second.push_back(str1);
			}
		}
	}

	// Global parameters loading
	globalParameters.clear();
	unsigned int globParamSize;
	std::string parName, parVal;
	stream >> globParamSize;
	for (unsigned int i = 0 ; i < globParamSize ; ++i)
	{
		stream >> parName;
		stream >> parVal;
		globalParameters[parName] = parVal;
	}

	// cond replace params
	condReplaceBeforeNames.clear();
	condReplaceAfterNames.clear();
	condReplaceVals.clear();
	unsigned int nbBefNames, nbAftNames, nbcondVals;
	std::string tmpStr;
	stream >> nbBefNames;
	for (unsigned int i = 0 ; i < nbBefNames ; ++i)
	{
		stream >> tmpStr;
		condReplaceBeforeNames.push_back(tmpStr);
	}
	stream >> nbAftNames;
	for (unsigned int i = 0 ; i < nbAftNames ; ++i)
	{
		stream >> tmpStr;
		condReplaceAfterNames.push_back(tmpStr);
	}
	stream >> nbcondVals;
	for (unsigned int i = 0 ; i < nbcondVals ; ++i)
	{
		condReplaceVals.push_back(std::make_pair(std::vector<std::string>(), 
			std::vector<std::string>()));
		for (unsigned int j = 0 ; j < condReplaceBeforeNames.size() ; ++j)
		{
			stream >> tmpStr;
			condReplaceVals.back().first.push_back(tmpStr);
		}
		for (unsigned int j = 0 ; j < condReplaceAfterNames.size() ; ++j)
		{
			stream >> tmpStr;
			condReplaceVals.back().second.push_back(tmpStr);
		}
	}

	// Sub simulation loading
	string subSimType;
	stream >> subSimType;
TRACE(subSimType)
	if (handler.getParam<bool>("-simLoad"))
	{
		string loadPath = handler.getParam<string>("-simLoad", 1);
		ifstream loadStream(loadPath.c_str());
TRACE(loadPath)
		assert(AbstractFactory<SimulationManager>::Factories[subSimType]);
		SimulationManager *tempSimManag = AbstractFactory<SimulationManager>::Factories[subSimType]
					->CreateFromStream(loadStream);
		//
		subSimulation = tempSimManag;
		//
	}
	metrics.FreeAndClean();
TRACE("Loading GridSearch Metrics ...")
	ok &= metrics.LoadFromStream(stream);
TRACE("Finished loading Gridsearch Metrics !")
	ok &= saver.LoadFromStream(stream);

	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool GridSearchSimulation::SaveToStream(std::ofstream & stream) const
{
	bool ok = true;

	// grid parameters
	stream << gridParameters.size() << endl;
	for (unsigned int i = 0 ; i < gridParameters.size() ; ++i)
	{
		stream << gridParameters[i].first << endl;
		stream << gridParameters[i].second.size() << endl;
		for (unsigned int j = 0 ; j < gridParameters[i].second.size() ; ++j)
			stream << gridParameters[i].second[j] << endl;
	}
	// conditional grid parameters
	stream << condGridParameters.size() << endl;
	for (std::map<std::pair<std::string, std::string>, 
		std::vector<std::pair<std::string, std::vector<std::string> > > >::
		const_iterator it = condGridParameters.begin() ; 
		it != condGridParameters.end() ; ++it)
	{
		stream << it->first.first << std::endl;
		stream << it->first.second << std::endl;
		stream << it->second.size() << std::endl;
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
		{
			stream << it->second[i].first << std::endl;
			stream << it->second[i].second.size() << std::endl;
			for (unsigned j = 0 ; j < it->second[i].second.size() ; ++j)
				stream << it->second[i].second[j] << std::endl;
		}
	}
	// global parameters
	stream << globalParameters.size() << std::endl;
	for (std::map<std::string, std::string>::const_iterator it = globalParameters.begin() ; it != globalParameters.end() ; ++it)
		stream << it->first << std::endl << it->second << std::endl;
	// cond replace params
	stream << condReplaceBeforeNames.size() << std::endl;
	for (unsigned int i = 0 ; i < condReplaceBeforeNames.size() ; ++i)
		stream << condReplaceBeforeNames[i] << std::endl;
	stream << condReplaceAfterNames.size() << std::endl;
	for (unsigned int i = 0 ; i < condReplaceAfterNames.size() ; ++i)
		stream << condReplaceAfterNames[i] << std::endl;
	stream << condReplaceVals.size() << std::endl;
	for (unsigned int i = 0 ; i < condReplaceVals.size() ; ++i)
	{
		for (unsigned int j = 0 ; j < condReplaceVals[i].first.size() ; ++j)
			stream << condReplaceVals[i].first[j] << std::endl;
		for (unsigned int j = 0 ; j < condReplaceVals[i].second.size() ; ++j)
			stream << condReplaceVals[i].second[j] << std::endl;
	}
	// save sub simulation
	assert(subSimulation);
	stream << subSimulation->GetClassName() << endl;
	// metrics
	ok &= metrics.SaveToStream(stream);
	// saver
	ok &= saver.SaveToStream(stream);

	return ok and stream.good();
}

//**********************************************************************
// Returns a vector sontaining the values of the current parameter 
// combination.
//**********************************************************************
vector<string> GridSearchSimulation::GetCurrParamComb() const
{
	vector<string> temp;
	unsigned int ind;
	for (unsigned int i = 0 ; i < totalParamNames.size() ; ++i)
	{
		ind = getParamInd(totalParamNames[i], gridParameters);
		if (ind < gridParameters.size())
			temp.push_back(gridParameters[ind].second[currIndices[ind]]);
		else
		{
			ind = getParamInd(totalParamNames[i], currCondGrid);
			if (ind < currCondGrid.size())
				temp.push_back(currCondGrid[ind].second[currCondInd[ind]]);
			else
				temp.push_back("NULL");
		}
	}
	return temp;
}

//**********************************************************************
// Returns the name of the parameters
//**********************************************************************
vector<string> GridSearchSimulation::GetParamNames() const
{
	return totalParamNames;
}

//**********************************************************************
// Returns the total grid size
//**********************************************************************
unsigned int GridSearchSimulation::GetFullGridSize() const
{
	unsigned int totalGridSize = 0;
	unsigned int tempSize = 0;
	std::map<std::pair<std::string, std::string>, 
		std::vector<std::pair<std::string, 
		std::vector<std::string> > > >::const_iterator it;
	std::vector<unsigned int> tempInd(gridParameters.size(), 0);
	while (not gridParameters.empty() and
		(tempInd[0] < gridParameters[0].second.size()))
	{
		tempSize = 1;
		for (unsigned int i = 0 ; i < tempInd.size() ; ++i)
		{
			it = condGridParameters.find(std::make_pair(gridParameters[i].first, 
				gridParameters[i].second[tempInd[i]]));
			if (it != condGridParameters.end())
			{
				for (unsigned int j = 0 ; j < it->second.size() ; ++j)
					tempSize *= it->second[j].second.size();
			}
		}
		totalGridSize += tempSize;

		// Increase combination
		tempInd.back()++;
		for (unsigned int i = tempInd.size() - 1 ; i > 0 ; --i)
			if (tempInd[i] >= gridParameters[i].second.size())
			{
				tempInd[i] = 0;
				++tempInd[i-1];
			}
	}
	return totalGridSize;
}

//**********************************************************************
// Returns the indice of a parameter if it exists, otherwise returns
// gridParam size
//**********************************************************************
unsigned int GridSearchSimulation::getParamInd(string name, 
	const std::vector<std::pair<std::string, 
	std::vector<std::string> > > & paramStruct) const
{
	unsigned int ind = 0;
	while ((ind < paramStruct.size()) and (paramStruct[ind].first != name))
		++ind;
	return ind;
}

//**********************************************************************
// Update the grid values through paramHandler
// Returns the path to be used by the resultSaver during simulations
//**********************************************************************
std::string GridSearchSimulation::updateGridValuesAndGetPath()
{
	// Check whether a subset of the current comb is to be replaced
	std::vector<std::string> tmpSubComb;
TRACE(condReplaceBeforeNames.size())
	bool found = false;
	for (unsigned int i = 0 ; i < condReplaceBeforeNames.size() ; ++i)
	{
		found = false;
		// Check in gridParameters
		for (unsigned int j = 0 ; not found and j < gridParameters.size() ; ++j)
			if (gridParameters[j].first == condReplaceBeforeNames[i])
			{
				tmpSubComb.push_back(gridParameters[j].second[currIndices[j]]);
				found = true;
			}
		// Check in currCondGrid
		for (unsigned int j = 0 ; not found and j < currCondGrid.size() ; ++j)
			if (currCondGrid[j].first == condReplaceBeforeNames[i])
			{
				tmpSubComb.push_back(currCondGrid[j].second[currCondInd[j]]);
				found = true;
			}
		if (not found)
			tmpSubComb.push_back("NULL");
	}
	// Check if curr comb is in condReplaceVals
	found = false;
	unsigned int tmpInd;
	bool tmpOk;
	for (tmpInd = 0 ; not found and  tmpInd < condReplaceVals.size() ; ++tmpInd)
	{
		tmpOk = true;
		for (unsigned int j = 0 ; tmpOk and (j < condReplaceVals[tmpInd].first.size()) ; ++j)
		{
			tmpOk &= (condReplaceVals[tmpInd].first[j] == tmpSubComb[j]);
		}
		found |= tmpOk;
	}
	--tmpInd;


	ParamHandler params = subSimulation->BuildParamHandler();
	string tempStr = "";
	for (int i = gridParameters.size() - 1 ; i >= 0 ; --i)
	{
		// If hasn't been replaced
		if (not found or (std::find(condReplaceBeforeNames.begin(), 
			condReplaceBeforeNames.end(), gridParameters[i].first) == 
			condReplaceBeforeNames.end()))
		{
			for (unsigned int j = 0 ; 
					j < params.GetNbParamsForName(gridParameters[i].first) 
					; ++j)
				params.SetVal(gridParameters[i].first, 
					gridParameters[i].second[currIndices[i]].c_str(), j);
		}

		tempStr += gridParameters[i].first + "_" + 
			gridParameters[i].second[currIndices[i]] + "/";
	}
	for (int i = currCondGrid.size() - 1 ; i >= 0 ; --i)
	{
		if (not found or (std::find(condReplaceBeforeNames.begin(), 
			condReplaceBeforeNames.end(), currCondGrid[i].first) == 
			condReplaceBeforeNames.end()))
		{
			for (unsigned int j = 0 ; 
					j < params.GetNbParamsForName(currCondGrid[i].first) 
					; ++j)
				params.SetVal(currCondGrid[i].first, 
					currCondGrid[i].second[currCondInd[i]].c_str(), j);
		}

		tempStr += currCondGrid[i].first + "_" + 
			currCondGrid[i].second[currCondInd[i]] + 
			((i == 0) ? "" : "/");
	}
	// Update to be replaced values
	if (found)
	{
		for (unsigned int i = 0 ; i < condReplaceAfterNames.size() ; ++i)
		{
			if (condReplaceVals[tmpInd].second[i] != "NULL")
			{
				for (unsigned int j = 0 ; 
						j < params.GetNbParamsForName(condReplaceAfterNames[i]) ; ++j)
					params.SetVal(condReplaceAfterNames[i], 
						condReplaceVals[tmpInd].second[j].c_str(), j);
			}
		}
	}
	return tempStr;
}

//**********************************************************************
// Separate the current grid search in two smaller gridsearches
//**********************************************************************
std::pair<GridSearchSimulation*, GridSearchSimulation*> 
	GridSearchSimulation::separateinTwo() const
{
	GridSearchSimulation *newSim1 = new GridSearchSimulation(subSimulation, false);
	GridSearchSimulation *newSim2 = new GridSearchSimulation(subSimulation, false);
	// copy the saver
	newSim1->saver = saver;
	newSim2->saver = saver;
	// Fill globalParameters
	newSim1->globalParameters = globalParameters;
	newSim2->globalParameters = globalParameters;
	// Fill metrics
	for (unsigned int i = 0 ; i < metrics.size() ; ++i)
	{
		newSim1->metrics.push_back(make_pair(metrics[i].first, false));
		newSim2->metrics.push_back(make_pair(metrics[i].first, false));
	}
	unsigned int divideInd = 0;
	// Find a parameter with more than one value
	for (divideInd = 0 ; (divideInd < gridParameters.size()) and
		(gridParameters[divideInd].second.size() == 1) ; ++divideInd)
	{
		newSim1->gridParameters.push_back(gridParameters[divideInd]);
		newSim2->gridParameters.push_back(gridParameters[divideInd]);
	}

	if (divideInd < gridParameters.size())
	{
		newSim1->gridParameters.push_back(std::make_pair(
			gridParameters[divideInd].first, std::vector<std::string>()));
		newSim2->gridParameters.push_back(std::make_pair(
			gridParameters[divideInd].first, std::vector<std::string>()));

		unsigned int maxInd1 = gridParameters[divideInd].second.size() / 2;
		std::pair<std::string, std::string> key;
		for (unsigned int i = 0 ; i < maxInd1 ; ++i)
		{
			newSim1->gridParameters.back().second.push_back(
				gridParameters[divideInd].second[i]);
			key = std::make_pair(gridParameters[divideInd].first, 
				gridParameters[divideInd].second[i]);
			if (condGridParameters.find(key) != condGridParameters.end())
				newSim1->condGridParameters[key] = condGridParameters[key];
		}
		for (unsigned int i = maxInd1 ; 
				i < gridParameters[divideInd].second.size() ; ++i)
		{
			newSim2->gridParameters.back().second.push_back(
				gridParameters[divideInd].second[i]);
			key = std::make_pair(gridParameters[divideInd].first, 
				gridParameters[divideInd].second[i]);
			if (condGridParameters.find(key) != condGridParameters.end())
				newSim2->condGridParameters[key] = condGridParameters[key];
		}

		// Add other parameters
		for (unsigned int i = divideInd + 1 ; i < gridParameters.size() ; ++i)
		{
			newSim1->gridParameters.push_back(gridParameters[i]);
			newSim2->gridParameters.push_back(gridParameters[i]);
		}

		// Add other condGridParams
		for (std::map<std::pair<std::string, std::string>, 
			std::vector<std::pair<std::string, std::vector<std::string> > > >::
			const_iterator it = condGridParameters.begin() ; 
			it != condGridParameters.end() ; ++it)
		{
			if (it->first.first != gridParameters[divideInd].first)
			{
				newSim1->condGridParameters.insert(*it);
				newSim2->condGridParameters.insert(*it);
			}
		}
	}
	else
	{
		std::vector<std::pair<std::string, 
			std::vector<std::string> > > currCondGridTmp;

		std::pair<std::string, std::string> key;
		std::vector<std::pair<std::string, std::string> > refPars;
		// fill currCondGrid
		for (unsigned int i = 0 ; i < gridParameters.size() ; ++i)
		{
			key = std::make_pair(gridParameters[i].first, gridParameters[i].second[0]);
			if (condGridParameters.find(key) != condGridParameters.end())
			{
				for (unsigned int j = 0 ; j < condGridParameters[key].size() ; ++j)
				{
					currCondGridTmp.push_back(condGridParameters[key][j]);
					refPars.push_back(key);
				}
			}
		}

		unsigned int condDivInd = 0;
		// Find a cond parameter with more than one value
		for (condDivInd = 0 ; (condDivInd < currCondGridTmp.size()) and 
			(currCondGridTmp[condDivInd].second.size() == 1) ; ++condDivInd)
		{
			newSim1->condGridParameters[refPars[condDivInd]].push_back(
				currCondGridTmp[condDivInd]);
			newSim2->condGridParameters[refPars[condDivInd]].push_back(
				currCondGridTmp[condDivInd]);
		}

		if (condDivInd < currCondGridTmp.size())
		{
			newSim1->condGridParameters[refPars[condDivInd]].push_back(
				std::make_pair(currCondGridTmp[condDivInd].first, 
				std::vector<std::string>()));
			newSim2->condGridParameters[refPars[condDivInd]].push_back(
				std::make_pair(currCondGridTmp[condDivInd].first, 
				std::vector<std::string>()));

			unsigned int nbDiv = currCondGridTmp[condDivInd].second.size() / 2;
			for (unsigned int i = 0 ; i < nbDiv ; ++i)
				newSim1->condGridParameters[refPars[condDivInd]].back().
					second.push_back(currCondGridTmp[condDivInd].second[i]);
			for (unsigned int i = nbDiv ; i < currCondGridTmp[condDivInd].second.size() ; ++i)
				newSim2->condGridParameters[refPars[condDivInd]].back().
					second.push_back(currCondGridTmp[condDivInd].second[i]);


			for (unsigned int i = condDivInd + 1 ; 
				(i < currCondGridTmp.size()) ; ++i)
			{
				newSim1->condGridParameters[refPars[i]].push_back(
					currCondGridTmp[i]);
				newSim2->condGridParameters[refPars[i]].push_back(
					currCondGridTmp[i]);
			}
		}
		else
		{
			TRACE("ERROR")
			assert(false);
		}
	}

	return std::make_pair(newSim1, newSim2);
}

//**********************************************************************
// Copy the current GridSearch
//**********************************************************************
GridSearchSimulation* GridSearchSimulation::getCopy() const
{
	GridSearchSimulation *newSim1 = new GridSearchSimulation(subSimulation, false);
	// copy the saver
	newSim1->saver = saver;
	// Fill globalParameters
	newSim1->globalParameters = globalParameters;
	// Fill gridparameters
	newSim1->gridParameters = gridParameters;
	// Fill cond parameters
	newSim1->condGridParameters = condGridParameters;
	// Fill metrics
	for (unsigned int i = 0 ; i < metrics.size() ; ++i)
		newSim1->metrics.push_back(make_pair(metrics[i].first, false));

	return newSim1;
}

