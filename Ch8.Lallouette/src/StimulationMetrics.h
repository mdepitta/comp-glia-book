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

#ifndef STIMULATIONMETRICS_H
#define STIMULATIONMETRICS_H

#include "StimulationStrat.h"
#include "MetricComputeStrat.h"

#include <vector>

namespace AstroModel
{
	// Forward declarations
	class StimulationStrat;
	class StimulationInd;

	// Typedefs
	typedef SpecificMetric<StimulationStrat> StimulationMetric;

/**********************************************************************/
/* Stimulation Metrics : Stimulated Cells                             */
/**********************************************************************/
	class StimulatedCells : public StimulationMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		StimulatedCells();
		// Constructor from stream
		StimulatedCells(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const StimulationStrat & stimStrat);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual const std::set<unsigned int> & GetCumulStimCells() const
			{ return stimCells; }
		virtual const std::string & GetCurrOrigin() const
			{ return currOrigin; }

	protected:
		double lastT;
		std::string currOrigin;
		std::vector<std::pair<bool, double> > stimulated;
		std::set<unsigned int> stimCells;

		mutable std::vector<StimulationInd> stimulations;

		static bool stimFileHasBeenSaved;
	};

/**********************************************************************/
/* Stimulation Metrics : Multiple stimulation mean path metric        */
/**********************************************************************/
	class MeanPathStimMetric : public StimulatedCells
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		MeanPathStimMetric() {}
		// Constructor from stream
		MeanPathStimMetric(std::ifstream & stream) { LoadFromStream(stream); }

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const StimulationStrat & stimStrat);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		//virtual bool LoadFromStream(std::ifstream & stream);
		//virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		double effMeanPathLength;
	};
}
#endif
