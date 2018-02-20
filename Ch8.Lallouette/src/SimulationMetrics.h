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

#ifndef SIMULATIONMETRICS_H
#define SIMULATIONMETRICS_H

#include "MetricComputeStrat.h"
#include "ChIModelMetrics.h"
#include "NetworkMetrics.h"

#include <vector>

namespace AstroModel
{
	// Forward declaration
	class AbstractRepeatSimulation;
	template <typename ModelType> class RepeatSimulation;
	class GridSearchSimulation;

/**********************************************************************/
/* Simulation Metrics : Repeat simulation Metric                      */
/**********************************************************************/
	class RepeatSimulationMetric : public SpecificMetric<AbstractRepeatSimulation>
	{
	};

/**********************************************************************/
/* Simulation Metrics : Propagation Distance Statistics               */
/**********************************************************************/
	class PropagationDistanceStat : public RepeatSimulationMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		PropagationDistanceStat();
		// Constructor from stream
		PropagationDistanceStat(std::ifstream & stream);
		virtual ~PropagationDistanceStat();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const AbstractRepeatSimulation & manag);
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

	protected:
		typedef std::map<unsigned int, std::vector<double> > SDistByCellMapType;
		typedef std::map<double, std::vector<double> > SDistByCDistT;
		std::vector<std::pair<double, double> > nbCellsByClustCoeff;

		SDistByCellMapType spatialDistByCellNb;
		SDistByCDistT spatialDistByCellDist;
		std::vector<double> maxNbCells;
		std::vector<double> maxWaveSpeed;
		std::vector<int> nbWaves;
	};

/**********************************************************************/
/* Simulation Metrics : Network topologies and spatial metrics        */
/**********************************************************************/
	class NetworkTopologyStats : public RepeatSimulationMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		NetworkTopologyStats();
		// Constructor from stream
		NetworkTopologyStats(std::ifstream & stream);
		virtual ~NetworkTopologyStats();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const AbstractRepeatSimulation & manag);
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

	protected:
		class CellToCellDistData
		{
		public:
			CellToCellDistData(double m, double s, double mi) : 
				meanDist(m), stdDevDist(s), minDist(mi) {}
			double meanDist;
			double stdDevDist;
			double minDist;
		};
		std::vector<CellToCellDistData> cellToCellDist;
	};

/**********************************************************************/
/* Simulation Metrics : Scalar Statistics across all repetitions      */
/**********************************************************************/
	class ScalarStatisticsMetric : public RepeatSimulationMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		ScalarStatisticsMetric();
		// Constructor from stream
		ScalarStatisticsMetric(std::ifstream & stream);
		virtual ~ScalarStatisticsMetric();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const AbstractRepeatSimulation & manag);
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
		const std::map<std::string, double> & GetScalarStats() const
			{ return scalarStats; }

	protected:
		std::map<std::string, std::vector<double> > fullScalars;
		std::map<std::string, 
			std::vector<std::vector<double> > > fullDistribs;
		std::map<std::string, double> scalarStats;

		//===========================================================||
		// Specific ScalarStatisticsMetric methods                   ||
		//===========================================================||
		void updateScalarStats();
		std::string modifStatName(const std::string & str) const;
	};

/**********************************************************************/
/* Simulation Metrics : Grid Search Simulation Metric                 */
/**********************************************************************/
	class GridSearchMetric : public SpecificMetric<GridSearchSimulation>
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		GridSearchMetric() {}
		// Constructor from stream
		GridSearchMetric(std::ifstream & stream) { LoadFromStream(stream); }

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const GridSearchSimulation & manag);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		typedef std::vector<std::string> ParamComb;

		std::vector<std::string> paramNames;
		// ParamComb -> (FileName -> [filePathInd])
		std::map<ParamComb, 
			std::map<std::string, std::vector<int> > > savedFilesByParam;
		std::map<ParamComb, 
			std::map<std::string, double> > scalarStatsByParam;

		std::vector<std::string> paths;
		std::set<std::string> totalFileNames;
		std::set<std::string> totalScalarNames;
	};
}

#endif
