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

#ifndef GRIDSEARCHSIMULATION_H
#define GRIDSEARCHSIMULATION_H

#include "SimulationManager.h"
#include "ODESolvers.h"
#include "ChISimulationManager.h"

#include <vector>

#define DEFAULT_SIMULATION_PATH "./data/default"
#define RANGE_PERMISSIVITY 0.001

namespace AstroModel
{
	typedef SpecificMetric<GridSearchSimulation> GenericGridSearchMetric;
	
/**********************************************************************/
/* Grid search simulation                                             */
/**********************************************************************/
	class GridSearchSimulation : public SimulationManager
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Standard Constructor
		GridSearchSimulation(SimulationManager *_ss = 0, bool _f = false,
			ParamHandler & _h = ParamHandler::GlobalParams, 
			std::string path = DEFAULT_SIMULATION_PATH);
		// Constructor from stream
		GridSearchSimulation(std::ifstream & stream, 
			ParamHandler & _h = ParamHandler::GlobalParams);
		// Standard Destructor
		virtual ~GridSearchSimulation();

		//===========================================================||
		// Standard Simulation Manager methods                       ||
		//===========================================================||
		// Only calls LaunchSimulation with a suitable saver
		virtual int LaunchSimulation(std::string simName);
		// Launches the simulation and saves data with saver
		virtual int LaunchSimulation(ResultSaver saver);
		virtual ParamHandler BuildParamHandler() 
			{ return ParamHandler(); }
		virtual const std::vector<Metric *> & GetMetrics() const
			{ static std::vector<Metric *> temp; return temp; }
		// Returns its own metrics and the metric tree of all its "children" objects
		virtual std::vector<Metric *> GetAllMetrics() const;

		//===========================================================||
		// Grid search specific methods                              ||
		//===========================================================||
		template <typename T> void AddGridParameterValue(
			std::string name, const T & v)
		{
			unsigned int ind;
			if ((ind = getParamInd(name, gridParameters)) ==
				gridParameters.size())
			{
				gridParameters.push_back(make_pair(name, 
					std::vector<std::string>()));
				ind = gridParameters.size() - 1;
			}
			gridParameters[ind].second.push_back(Stringify(v));
			if (std::find(totalParamNames.begin(), totalParamNames.end(), name) ==
					totalParamNames.end())
				totalParamNames.push_back(name);
		}

		template <typename T> void AddGridParameterRange(std::string name, 
			const T & start, const T & end, const T & step)
		{
			for (T val = start ; val <= (end + step * RANGE_PERMISSIVITY) ; val += step)
				AddGridParameterValue(name, val);
		}

		template <typename T> void AddCondGridParameterValue(
			std::string condName, std::string condVal, std::string name, 
			const T & val)
		{
			std::pair<std::string, std::string> key = 
				std::make_pair(condName, Stringify(condVal));
			unsigned int ind = getParamInd(name, condGridParameters[key]);
			if (ind == condGridParameters[key].size())
			{
				condGridParameters[key].push_back(std::make_pair(name,
					std::vector<std::string>()));
				ind = condGridParameters[key].size() - 1;
			}
			condGridParameters[key][ind].second.push_back(Stringify(val));
			if (std::find(totalParamNames.begin(), totalParamNames.end(), name) ==
					totalParamNames.end())
				totalParamNames.push_back(name);
		}

		template <typename T> void AddCondGridParameterRange(
			std::string condName, std::string condVal, std::string name, 
			const T & start, const T & end, const T & step)
		{
			for (T val = start ; val <= (end + step * RANGE_PERMISSIVITY) ; val += step)
				AddCondGridParameterValue(condName, condVal, name, val);
		}

		template <typename T> void AddGlobalParameterValue(std::string name, const T & v)
		{
			globalParameters[name] = Stringify(v);
		}

		virtual bool AddMetric(Metric *_m, bool _f = false)
		{
			if (metrics.AddMetricAndDependencies(_m, _f, this))
				return true;
			else
				return subSimulation ? subSimulation->AddMetric(_m, _f) : false;
		}
		// Cut a large grid search into smaller subGridsearches in order to save them
		std::vector<GridSearchSimulation *> CreateSmallerGridSearches(
			unsigned int desiredSize, int tmp = 1) const;
		// Compute smaller grid searches and save them
		bool ComputeAndSaveGridSearches(unsigned int nb, std::string scriptPath, 
			std::string gridPath, std::string simPath);
		// Add a file to be saved
		void AddFileToSave(std::string name);
		// Load conditional replace from file (saved by GridSearchMetric)
		bool LoadCondReplaceFromFile(std::string path);

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the cell from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the cell to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		virtual SimulationManager * GetSubSimulation() const
		{ return subSimulation; }
		virtual std::vector<std::string> GetCurrParamComb() const;
		virtual std::vector<std::string> GetParamNames() const;
		virtual const std::map<std::string, std::string> & GetGlobalParams() const
		{ return globalParameters; }
		virtual unsigned int GetFullGridSize() const;

	protected:
		std::vector<std::pair<std::string, std::vector<std::string> > > gridParameters;
		//
		mutable std::map<std::pair<std::string, std::string>, std::vector<std::pair<std::string, std::vector<std::string> > > > condGridParameters;
		std::vector<std::pair<std::string, std::vector<std::string> > > currCondGrid;
		//
		std::map<std::string, std::string> globalParameters;
		std::vector<unsigned int> currIndices;
		//
		std::vector<unsigned int> currCondInd;

		// All grid parameter names 
		// (including conditional ones but excluding global ones)
		std::vector<std::string> totalParamNames;
		//

		// Conditional replacment of values
		std::vector<std::string> condReplaceBeforeNames;
		std::vector<std::string> condReplaceAfterNames;
		std::vector<std::pair<std::vector<std::string>, 
			std::vector<std::string> > > condReplaceVals;

		SimulationManager *subSimulation;
		bool freeSubSimulation;

		SortedMetrics<GridSearchSimulation> metrics;

		ParamHandler & handler;
		ResultSaver saver;

		//===========================================================||
		// Protected methods                                         ||
		//===========================================================||
		unsigned int getParamInd(std::string name, const std::vector<std::pair<std::string, std::vector<std::string> > > & paramStruct) const;
		std::string updateGridValuesAndGetPath();
		// Returns true if debug mode is on
		inline bool DebugMode() const 
		{ return handler.getParam<bool>("-debug"); }
		virtual std::pair<GridSearchSimulation*, GridSearchSimulation*> separateinTwo() const;
		virtual GridSearchSimulation* getCopy() const;
	};
}

#endif
