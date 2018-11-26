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

#ifndef ASTRONEUROMODEL_H
#define ASTRONEUROMODEL_H

#include <vector>
#include "Model.h"
#include "ODEProblems.h"
#include "ParamHandler.h"
#include "MetricComputeStrat.h"
#include "Synapse.h"
#include "Neuron.h"

namespace ODE
{
	// Forward declarations
}

namespace AstroModel
{
	// Forward declarations

/**********************************************************************/
/* Neuron Network Model (abstract base class)                         */
/**********************************************************************/
	class AstroNeuroNetModel : 
		public Model,
		public ODE::ODEProblem<double, double>
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::AstroNeuronNetFunc;

	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		AstroNeuroNetModel(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		AstroNeuroNetModel(std::ifstream & stream, 
			ParamHandler & h = ParamHandler::GlobalParams);
		// Destructor
		virtual ~AstroNeuroNetModel();

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Allocate memory and sets up model with corresponding dynVals
		virtual void SetUpModelsAndODEs();
		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver);
		// Launch simulation
		virtual int Simulate(ResultSaver saver);
		// Method called before solving the ODE problem
		virtual bool PreSimulationCall(ResultSaver saver);
		// Method called after solving the ODE problem
		virtual bool PostSimulationCall(ResultSaver saver);
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model computations                                ||
		//===========================================================||

		//===========================================================||
		// Setters and callbacks                                     ||
		//===========================================================||
		// Notifies the model that all values have been updated for 
		// timestep t
		virtual void UpdateVals(double t);
		// Set all cells to equilibrium
		void SetAllCellsToEquilibrium();

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a network metric to the network
		virtual bool AddMetric(Metric *_m, bool _f = false);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const;
		// Is the given cell currently stimulated ?
		bool IsStimulated(unsigned int ind) const;

	protected:

		virtual double getModelVersionNum() const { return 0.1; }

		std::string neuronNetClassName;
		std::string astroNetClassName;

		NeuronNetModel *neuronNet;
		ChIModel *astroNet;

		std::vector<std::vector<Synapse *> > astrToSyn;

		//===========================================================||
		// Network parameters                                        ||
		//===========================================================||

		//===========================================================||
		// ODE Solver                                                ||
		//===========================================================||
		ODE::ODESolver<double, double> * solver;

		//===========================================================||
		// Metrics                                                   ||
		//===========================================================||
		SortedMetrics<AstroNeuroNetModel> metrics;
	};

}

#endif

