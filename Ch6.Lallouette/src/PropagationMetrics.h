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

#ifndef PROPAGATIONMETRICS_H
#define PROPAGATIONMETRICS_H

#include "MetricComputeStrat.h"

#include <vector>

namespace AstroModel
{
	// Forward declarations
	template <typename StateType, typename LinkType> class PropagationModel;

/**********************************************************************/
/* Abstract Base Class for cell state savers                          */
/**********************************************************************/
	class AbstractCellStateSaver
	{
	public:
		virtual const std::map<std::string, std::set<unsigned int> > & 
			GetLastCumulStates() const = 0;
		virtual unsigned int GetNbStimulatedNodes() const = 0;
	};

/**********************************************************************/
/* Evolution of cell states over time                                 */
/**********************************************************************/
	template <typename StateType, typename LinkType>
	class CellStateSaver : 
		public SpecificMetric<PropagationModel<StateType, LinkType> >, 
		public virtual AbstractCellStateSaver
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		CellStateSaver(ParamHandler & h = ParamHandler::GlobalParams) :
			currTime(0)
		{
			timeStep = h.getParam<double>("-SavingStep", 0);
			Initialize();
		}
		// Full constructor
		CellStateSaver(double _ts) : timeStep(_ts), currTime(0) 
		{ 
			Initialize();
		}
		// Constructor from stream
		CellStateSaver(std::ifstream & stream) : currTime(0)
		{ LoadFromStream(stream); }

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const PropagationModel<StateType, LinkType> & model)
		{
			double modTime = model.GetTime();
			if (modTime - currTime >= timeStep)
			{
				currTime = modTime;
				allStates.push_back(
					std::make_pair(currTime, std::vector<StateType>()));
				for (unsigned int i = 0 ; i < model.GetNbNodes() ; ++i)
				{
					allStates.back().second.push_back(model.GetNodeState(i));
					cumulStates[model.GetNodeState(i).GetName()].insert(i);
				}
			}
			nbStimNodes = model.GetNbStimulatedNodes();
			return true;
		}

		virtual bool SaveMetric(ResultSaver saver) const
		{
			bool ok = true;
			std::string detailedActivCells("DetailedActivatedCells");
			if (saver.isSaving(detailedActivCells) and not allStates.empty())
			{
				std::ofstream & stream = saver.getStream();
				this->AddSavedFile(detailedActivCells, saver.getCurrFile());

				stream << "Time";
				for (unsigned int i = 0 ; i < allStates[0].second.size() ; ++i)
					stream << "\tCell_" << StringifyFixed(i);
				stream << std::endl;

				for (unsigned int t = 0 ; t < allStates.size() ; ++t)
				{
					stream << allStates[t].first;
					for (unsigned int i = 0 ; i < allStates[t].second.size() ; ++i)
						stream << "\t" << allStates[t].second[i];
					stream << std::endl;
				}

				ok &= stream.good();
			}
			std::string propagActivCells("PropagActivatedCells");
			if (saver.isSaving(propagActivCells) and not allStates.empty())
			{
				std::vector<std::set<unsigned int> > cumulCounts;
				std::vector<unsigned int> nbEachState;

				std::ofstream & stream = saver.getStream();
				this->AddSavedFile(propagActivCells, saver.getCurrFile());

				stream << "Time";
				const std::vector<std::string> & stateNames = StateType::GetStateNames();
				const std::vector<StateType> & stateVals = StateType::GetStateVals();
				for (unsigned int i = 0 ; i != stateNames.size() ; ++i)
				{
					cumulCounts.push_back(std::set<unsigned int>());
					stream << "\tNbCurr_" << stateNames[i] << "\tNbCumul_" << stateNames[i];
				}
				stream << std::endl;

				for (unsigned int i = 0 ; i < allStates.size() ; ++i)
				{
					stream << allStates[i].first;
					nbEachState = std::vector<unsigned int>(stateVals.size(), 0);
					for (unsigned int j = 0 ; j < allStates[i].second.size() ; ++j)
					{
						for (unsigned int k = 0 ; k < stateVals.size() ; ++k)
							if (allStates[i].second[j] == stateVals[k])
							{
								++nbEachState[k];
								cumulCounts[k].insert(j);
								break;
							}
					}
					for (unsigned int j = 0 ; j != stateNames.size() ; ++j)
					{
						stream << "\t" << nbEachState[j] << "\t" << cumulCounts[j].size();
					}
					stream << std::endl;
				}

				ok &= stream.good();
			}
			return ok;
		}
		
		virtual void Initialize()
		{
			Metric::Initialize();
			currTime = 0;
			allStates.clear();
			cumulStates.clear();
			for (std::vector<std::string>::const_iterator it = 
					StateType::GetStateNames().begin() ;
					it != StateType::GetStateNames().end() ; ++it)
				cumulStates.insert(std::make_pair(*it, std::set<unsigned int>()));
		}

		virtual Metric * BuildCopy() const
		{
			return new CellStateSaver<StateType, LinkType>(timeStep); 
		}

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			stream >> timeStep;
			Initialize();
			return stream.good() and not stream.eof();
		}

		virtual bool SaveToStream(std::ofstream & stream) const
		{
			stream << timeStep << std::endl;
			return stream.good();
		}

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual const std::map<std::string, std::set<unsigned int> > & 
			GetLastCumulStates() const
		{
			return cumulStates;
		}
		virtual unsigned int GetNbStimulatedNodes() const { return nbStimNodes; }

	protected:
		double timeStep;

		double currTime;
		std::vector<std::pair<double, std::vector<StateType> > > allStates;
		std::map<std::string, std::set< unsigned int> > cumulStates;
		unsigned int nbStimNodes;
	};

}

#endif
