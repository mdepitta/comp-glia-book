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

#ifndef ODEPROBLEMS_H
#define ODEPROBLEMS_H

#include <vector>
#include "ODEFunctions.h"
#include "ResultSaver.h"
#include "MetricComputeStrat.h"
#include "AbstractFactory.h"

namespace ODE
{
	// Forward declaration
	class ChINetworkFunct;
	template <typename Val, typename StepT> class ODESolver;

/**********************************************************************/
/* Base class                                                         */
/**********************************************************************/
	template<typename Val = double, typename TStep = double>
	class ODEProblem : public SaveAndLoadFromStream
	{
	public:
		Function<Val, TStep> * function;
		bool freeFunct;
		unsigned int nbVals;
		Val *initVals;
		Val *vals;
		Val *derivs;
		bool freeVals;

		mutable std::vector<std::string> valNames;

		ODEProblem(Function<Val, TStep> * _f = 0, ParamHandler & param = 
			ParamHandler::GlobalParams, bool _ff = false) :
			function(_f), freeFunct(_ff), nbVals(0),
			initVals(0), vals(0), derivs(0), freeVals(false),
			preRunToEqu(false), preRunTime(0), isPreRunning(false), solver(0)
		{
			integrStep = param.getParam<double>("-Step");
			nOut = param.getParam<double>("-SavingStep") / integrStep;
			tStart = param.getParam<double>("-tStart");
			tEnd = param.getParam<double>("-tEnd");

			// Simulation Parameters
			preRunToEqu = param.getParam<bool>("-PreRunTimeToEq", 0);
			preRunTime = param.getParam<double>("-PreRunTimeToEq", 1);
			
			// ODE Solver
			std::string solverClassName = param.getParam<std::string>("-SolverClass");
			if (AstroModel::AbstractFactory<ODE::ODESolver<double, double> >::Factories[solverClassName])
			{
				solver = AstroModel::AbstractFactory<ODE::ODESolver<double, double> >::Factories[solverClassName]->Create();
				solver->SetStepSize(param.getParam<double>("-Step"));
			}
			else
				std::cerr << "Couldn't create the following ODE solver : " << solverClassName << std::endl;
		}

		virtual ~ODEProblem()
		{
			FreeMemory();
			if (freeFunct)
				delete function;
			if (solver)
				delete solver;
		}

		virtual void AllocateMemory(unsigned int _nbVals) 
		{
			FreeMemory();
			assert(_nbVals > 0);
			nbVals = _nbVals;
			initVals = new Val[nbVals];
			vals = new Val[nbVals];
			derivs = new Val[nbVals];
			valNames = std::vector<std::string>(nbVals, "");
			freeVals = true;
		}

		// If _vals is not null, free actual vals and change the pointer
		virtual void SetVals(Val *_vals = 0, unsigned int _nbVals = 0)
		{
			if (_vals)
			{
				FreeMemory();
				vals = _vals;
				if (_nbVals > 0)
					nbVals = _nbVals;
				initVals = new Val[nbVals];
				derivs = new Val[nbVals];
				valNames = std::vector<std::string>(nbVals, "");
				freeVals = false;
			}
		}

		virtual void FreeMemory()
		{
			if (initVals)
				delete[] initVals;
			if (vals and freeVals)
				delete[] vals;
			if (derivs)
				delete[] derivs;
			valNames.clear();
			freeVals = false;
		}
		virtual void UpdateVals(double) 
		{
			metrics.template ComputeMetrics<AstroModel::NeedFrequentUpdateMetric>(*this);
		}

		virtual void SetFunct(Function<Val, TStep> * _f, bool _ff = false)
		{
			if (function and freeFunct)
				delete function;
			function = _f;
			freeFunct = _ff;
		}

		virtual void SetToInitVals()
		{
			assert(nbVals > 0 and initVals and vals);
			for (unsigned int i = 0 ; i < nbVals ; ++i)
				vals[i] = initVals[i];
		}

		virtual void UseCurrentValsAsInitVals()
		{
			for (unsigned int i = 0 ; i < nbVals ; ++i)
				this->initVals[i] = this->vals[i];
		}

		inline unsigned int GetNbVals() const { return nbVals; }
		inline Val GetVal(unsigned int i) const { return vals[i]; }

		virtual void AddPostfixToValName(Val * v, std::string pf) const
		{
			assert(v >= vals);
			if (valNames.empty())
				valNames[(v - vals)] += "_";
			valNames[(v - vals)] += pf;
		}

		// Loads from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = true;

			stream >> tStart;
			stream >> tEnd;
			stream >> tCurr;
			stream >> nOut;
			stream >> integrStep;
			// Simulation Parameters
			stream >> preRunToEqu;
			stream >> preRunTime;

			// Solver
			if (solver)
				delete solver;
			std::string solvName;
			stream >> solvName;
			assert((AstroModel::AbstractFactory<ODE::ODESolver<Val, TStep> >::Factories[solvName]));
			solver = AstroModel::AbstractFactory<ODE::ODESolver<Val, TStep> >::
				Factories[solvName]->CreateFromStream(stream);
			// init vals filling
			for (unsigned int i = 0 ; i < this->nbVals ; ++i)
				this->initVals[i] = this->vals[i];

			// metrics
			metrics.FreeAndClean();
			ok &= metrics.LoadFromStream(stream);

			return ok and stream.good();
		}
		// Saves to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = true;

			stream 
				<< tStart << std::endl
				<< tEnd << std::endl
				<< tCurr << std::endl
				<< nOut << std::endl
				<< integrStep << std::endl;

			// Simulation Parameters
			stream 
				<< preRunToEqu << std::endl
				<< preRunTime << std::endl;

			// Solver
			assert(solver);
			stream << solver->GetClassName() << std::endl;
			solver->SaveToStream(stream);

			// Metrics
			ok &= metrics.SaveToStream(stream);

			return ok and stream.good();
		}
		// Initializes the model
		virtual void Initialize(ResultSaver )
		{
			tCurr = tStart;
			// Initialize metrics
			metrics.InitializeMetricsDefault();
		}

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a metric
		virtual bool AddMetric(AstroModel::Metric *_m, bool _f = false)
		{
			return metrics.AddMetricAndDependencies(_m, _f, this);
		}
		// Only compute metrics
		virtual bool ComputeAfterSimMetrics()
		{
			return metrics.template ComputeMetrics<AstroModel::AfterSimMetric>(*this);
		}
		// Only save metrics
		virtual bool SaveMetrics(ResultSaver saver)
		{
			return metrics.SaveMetricsDefault(saver);
		}

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Returns tStart
		virtual const double & GetTStart() const { return tStart; }
		// Return tEnd
		virtual const double & GetTEnd() const { return tEnd; }
		// Returns current simulation time
		virtual double GetTime() const { return tCurr; }
		// Return the integration step
		virtual double GetIntegrStep() const { return integrStep; }
		// Get dyn val names
		const std::vector<std::string> & GetValNames() const 
		{ return valNames; }
		// Change a dyn val name
		void SetValName(unsigned int ind, std::string _str)
		{
			assert(ind < valNames.size());
			valNames[ind] = _str;
		}
		//
		inline bool IsPreRunning() const { return isPreRunning; }

		virtual std::string GetClassName() const 
		{ return std::string("ODEProblem_") + typeid(Val).name() + 
			"_"  + typeid(TStep).name(); }

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler()
		{
			ParamHandler params;

			params <= "StartTime", tStart;
			params <= "EndTime", tEnd;

			return params;
		}

		// Returns all metrics and submetrics
		virtual std::vector<AstroModel::Metric *> GetAllMetrics() const
		{
			return metrics.GetMetricsRaw();
		}

	protected:

		//===========================================================||
		// Time and integration parameters                           ||
		//===========================================================||
		double tStart;     // Starting time of the simulation
		double tEnd;       // Ending time of the simulation
		double tCurr;      // Current time in the simulation
		unsigned int nOut; // Save data every nOut step
		double integrStep; // Integration step

		// Do we need to run a pre simulation to start from equilibrium values ?
		bool preRunToEqu;
		double preRunTime;
		bool isPreRunning;

		ODESolver<Val, TStep> * solver;

		//===========================================================||
		// Metrics                                                   ||
		//===========================================================||
		AstroModel::SortedMetrics<ODEProblem<Val, TStep> > metrics;

	};

}

#endif
