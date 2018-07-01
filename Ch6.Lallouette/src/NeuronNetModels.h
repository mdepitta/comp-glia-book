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

#ifndef NEURONNETMODELS_H
#define NEURONNETMODELS_H

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
	class SFALIFNeuronFunct;
	class TMSynapseFunct;
}

namespace AstroModel
{
	// Forward declarations

/**********************************************************************/
/* Neuron Network Model (abstract base class)                         */
/**********************************************************************/
	class NeuronNetModel : 
		public NetworkDynamicsModel<Synapse>,
		public ODE::ODEProblem<double, double>
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::NeuronNetworkFunc;

	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		NeuronNetModel(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		NeuronNetModel(std::ifstream & stream, 
			ParamHandler & h = ParamHandler::GlobalParams);
		// Destructor
		virtual ~NeuronNetModel();

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Sets up cells and allocate data for ODEs
		virtual void SetUpNeuronsAndODEs(unsigned int nbNeur);
		// Sets up ODE data and synaptic structures
		virtual void SetUpODEsAndSynapticStruct();
		// Change ODE vals to given pointer (for external use)
		virtual void SetVals(double *_vals = 0, unsigned int _nbVals = 0);
		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver);
		// Initializes the dynamic values
		virtual void InitializeDynVals();
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
		// Returns the number of neurons in the model
		inline unsigned int GetNbNeurons() const { return neurons.size(); }
		// Retursn the number of synapses in the model
		inline unsigned int GetNbSynapses() const { return synapses.size(); }
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const;
		// Is the given cell currently stimulated ?
		bool IsStimulated(unsigned int ind) const;
		// Deletes allocated neurons and clear vector
		void ClearNeurons();
		// Gives the total number of desired dyn vals (not equivalent to GetNbVals
		// from OPEProblem<double, double>
		virtual unsigned int GetTotNbDynVals() const;
		virtual Synapse * GetSynapse(unsigned int ind);
		// Return the indice of the given synapse
		virtual unsigned int GetSynInd(Synapse *) const;

	protected:

		virtual double getModelVersionNum() const { return 0.1; }

		std::string neuronClassName;
		std::string synapseClassName;

		std::vector<Neuron *> neurons; // Neurons, managed by this class
		std::vector<Synapse *> synapses; // Synapses, managed by Network
		
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
		SortedMetrics<NeuronNetModel> metrics;
	};

}

#endif
