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

#ifndef SIMULATIONMANAGER_H
#define SIMULATIONMANAGER_H

#include "Savable.h"
#include "ParamHandler.h"
#include "ResultSaver.h"


namespace AstroModel
{

	class Metric;

/**********************************************************************/
/* Base class                                                         */
/**********************************************************************/
	class SimulationManager : public SaveAndLoadFromStream
	{
	public:
		// Returns the class name
		virtual std::string GetClassName() const = 0;
		virtual ~SimulationManager() {}

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the cell from a stream
		virtual bool LoadFromStream(std::ifstream & stream) = 0;
		// Saves the cell to a stream
		virtual bool SaveToStream(std::ofstream & stream) const = 0;

		//===========================================================||
		// Standard Simulation Manager methods                       ||
		//===========================================================||
		virtual int LaunchSimulation(std::string )
		{
			std::cerr << "LaunchSimulation(std::string simName)"
				<< " is not implemented in this class" << std::endl;
			return 0;
		}
		virtual int LaunchSimulation(ResultSaver saver) = 0;
		virtual ParamHandler BuildParamHandler() = 0;
		virtual const std::vector<Metric *> & GetMetrics() const = 0;
		virtual std::vector<Metric *> GetAllMetrics() const = 0;
		virtual bool AddMetric(Metric *_m, bool _f = false) = 0;
	};
}

#endif

