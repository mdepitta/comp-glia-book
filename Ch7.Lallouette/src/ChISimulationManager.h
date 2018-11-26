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

#ifndef CHISIMULATIONMANAGER_H
#define CHISIMULATIONMANAGER_H

#include "SimulationManager.h"
#include "ParamHandler.h"
#include "ResultSaver.h"
#include "ODESolvers.h"
#include "SimulationMetrics.h"
#include "AbstractFactory.h"
#include "ErrorCodes.h"
#include "NeuronNetModels.h"

namespace AstroModel
{
/**********************************************************************/
/* Abstract Repeat Model Simulation Manager                           */
/**********************************************************************/
	class AbstractRepeatSimulation : public SimulationManager
	{
	public:
		virtual std::vector<AbstractRepeatSimulation *> CreateSmallerSimulations(unsigned int desired) = 0;
		virtual bool SaveModelToStream(std::ofstream & stream) const = 0;
		virtual ~AbstractRepeatSimulation() {}
	};


/**********************************************************************/
/* Repeat Model Simulation Manager                                    */
/**********************************************************************/
	template <typename ModelType>
	class RepeatSimulation : public AbstractRepeatSimulation
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		RepeatSimulation(ParamHandler & h = ParamHandler::GlobalParams) : model(0), freeModel(false), handler(h)
		{
			runStart = 0;
			runEnd = handler.getParam<unsigned int>("-repeat", 0) - 1;
			InitializeStandard();
		}
		// Standard Constructor
		RepeatSimulation(unsigned int _rs, unsigned int _re,
			ParamHandler & _h = ParamHandler::GlobalParams) :
			model(0), freeModel(false), runStart(_rs), runEnd(_re), handler(_h)
		{
			InitializeStandard();
		}
		// Constructor from stream
		RepeatSimulation(std::ifstream & stream, 
			ParamHandler & _h = ParamHandler::GlobalParams) :
			model(0), freeModel(false), runStart(0), runEnd(0), handler(_h)
		{
			InitializeStandard();
			LoadFromStream(stream);
		}
		// Copy Constructor
		RepeatSimulation(const RepeatSimulation<ModelType> & _m) :
			model(_m.model), freeModel(false), runStart(_m.runStart), 
			runEnd(_m.runEnd), handler(_m.handler)
		{
			for (unsigned int i = 0 ; i < _m.metrics.size() ; ++i)
				metrics.push_back(std::make_pair(_m.metrics[i].first, false));
		}
		// Default Destructor
		~RepeatSimulation()
		{
TRACE("Delete model")
			if (model and freeModel)
				delete model;
TRACE("Model deleted")
		}

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the cell from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = true;
			if (not stream.good())
				return false;
			metrics.FreeAndClean();

			// Simulation Metrics
			ok &= metrics.LoadFromStream(stream);

			stream >> runStart;
			stream >> runEnd;

			return ok and stream.good() and not stream.eof();
		}

		// Saves the cell to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = true;
			ok &= metrics.SaveToStream(stream);

			stream << runStart << std::endl;
			stream << runEnd   << std::endl;

			return ok and stream.good();
		}

		//===========================================================||
		// Standard Simulation Manager methods                       ||
		//===========================================================||
		virtual int LaunchSimulation(ResultSaver saver)
		{
			int returnVal = 0;
		//==============================//
			TRACE_UP("Initializing the "<< ModelType::ClassName << " model for repeat simulation.")
		//==============================//

			// Metric initialization
			metrics.InitializeMetricsDefault();

			for (unsigned int run = runStart ; run <= runEnd ; ++run)
			{
				TRACE_UP("Starting run " << run << ".")

				// Initializing model
				model->Initialize(saver("Run" + StringifyFixed(run)));
				// Simulation
				returnVal |= model->Simulate(saver("Run" + StringifyFixed(run)));
				// Computing metrics
				if (not metrics.ComputeMetricsDefault(*this))
					returnVal |= SIMULATION_METRIC_COMPUTATION_PROBLEM;

				TRACE_DOWN("Run " << run << " ended.")
			}

		//==============================//
			TRACE_DOWN("All runs are done.")
			TRACE_UP("Saving " << ClassName << " metrics.")
		//==============================//

			// Saving metrics
			if (not metrics.SaveMetricsDefault(
				saver("Summary_Runs_" + StringifyFixed(runStart) + "_to_" + StringifyFixed(runEnd))))
				returnVal |= SIMULATION_METRIC_SAVING_PROBLEM;

		//==============================//
			TRACE_DOWN(ClassName << " metrics have been saved.")
		//==============================//

			return returnVal;
		}

		virtual bool SaveModelToStream(std::ofstream & stream) const
		{
			return model ? model->SaveToStream(stream) : false;
		}

		bool LoadModelFromStream(std::ifstream & stream)
		{
			return model ? model->LoadFromStream(stream) : false;
		}

		virtual std::vector<AbstractRepeatSimulation *> CreateSmallerSimulations(unsigned int desired)
		{
			// First check if the simulation has a statistics metric on it.
			// If it does, repeats cannot be separated in different simulations
			if (GetSpecificMetric<Metric, ScalarStatisticsMetric>(GetMetrics()))
				desired = 1;

			std::vector<AbstractRepeatSimulation *> simulations;
			RepeatSimulation<ModelType> *tempPtr;
			unsigned int step = ceil((double)(runEnd - runStart + 1) / (double)desired);
			for (unsigned int i = runStart ; i <= runEnd ; i += step)
			{
				tempPtr = new RepeatSimulation<ModelType>(*this);
				tempPtr->runStart = i;
				tempPtr->runEnd   = std::min(i + step - 1, runEnd);
				simulations.push_back(tempPtr);
			}
			return simulations;
		}

		//===========================================================||
		// Standard Repeat Simulation Manager methods                ||
		//===========================================================||
		virtual ParamHandler BuildParamHandler()
		{
			ParamHandler params;
			assert(model);
			params += model->BuildModelParamHandler();
			params += metrics.BuildModelParamHandler();
			return params;
		}

		virtual bool AddMetric(Metric *_m, bool _f = false)
		{
			if (metrics.AddMetricAndDependencies(_m, _f, this))
				return true;
			else
				return model ? model->AddMetric(_m, _f) : false;
		}

		virtual const std::vector<Metric*> & GetMetrics() const
		{
			return metrics.GetMetricsRaw();
		}

		// Returns its own metrics and the metric tree of all its "children" objects
		virtual std::vector<Metric *> GetAllMetrics() const
		{
			std::vector<Metric *> temp = metrics.GetMetricsRaw();
			if (model)
			{
				std::vector<Metric *> tempModel = model->GetAllMetrics();
				temp.insert(temp.begin(), tempModel.begin(), tempModel.end());
			}
			return temp;
		}

	protected:
		ModelType *model;
		bool freeModel;

		unsigned int runStart;
		unsigned int runEnd;

		ParamHandler &handler;

		SortedMetrics<AbstractRepeatSimulation> metrics;

		//===========================================================||
		// Protected methods                                         ||
		//===========================================================||
		virtual void InitializeStandard()
		{
			if (model and freeModel)
				delete model;
			if (handler.getParam<bool>("-modelLoad"))
			{
				std::ifstream loadStream(handler.getParam<std::string>("-modelLoad", 1).c_str());
TRACE(handler.getParam<std::string>("-modelLoad", 1))
TRACE(ModelType::ClassName)
				model = new ModelType(loadStream, handler);
				freeModel = true;
			}
			else
			{
				model = new ModelType(handler);
				freeModel = true;
			}
		}

		// Returns true if debug mode is on
		inline bool DebugMode() const { return handler.getParam<bool>("-debug"); }
	};

}

#endif
