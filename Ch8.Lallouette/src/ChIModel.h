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

#ifndef CHIMODEL_H
#define CHIMODEL_H

#include <vector>
#include "Model.h"
#include "ODEProblems.h"
#include "ParamHandler.h"
#include "CouplingFunction.h"
#include "StimulationStrat.h"
#include "ChICell.h"
#include "Network.h"
#include "ChIModelMetrics.h"
#include "StimulationMetrics.h"


namespace ODE
{
	// Forward declarations
	class ChICellFunct;
	class ChINetworkFunct;
}

namespace AstroModel
{
	// Forward declarations
	class ChICell;
	class StimulationStrat;

/**********************************************************************/
/* ChI Model                                                          */
/**********************************************************************/
	class ChIModel : 
		public ODENetworkDynamicsModel<CouplingFunction, ChICell>, 
		public StimulableCellNetwork
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::ChICellFunct;
		friend class ODE::ChINetworkFunct;
		friend class ODE::KChICellFunct;
		friend class ODE::KChINetworkFunct;
		friend class ODE::AstroNeuronNetFunc;
		friend class KChICell;

		friend class PoissonianStimStrat;

	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		ChIModel(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		ChIModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams);
		// Destructor
		virtual ~ChIModel();

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Sets up cells and allocate data for ODEs
		virtual void SetUpCellsAndODEs(unsigned int nbCells);
		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver);
		// Method called before solving the ODE problem
		virtual bool PreSimulationCall(ResultSaver saver);
		// Method called after solving the ODE problem
		virtual bool PostSimulationCall(ResultSaver saver);
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Load funct from stream
		virtual bool LoadFunctFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model computations                                ||
		//===========================================================||
		// Computes fluxes across cells
		virtual void ComputeFluxes(double t);

		//===========================================================||
		// Setters and callbacks                                     ||
		//===========================================================||
		// Changes the flux of a cell
		void ModifFluxes(unsigned int i, double flux);
		// Stimulates Ca2+ release from an astrocyte
		void StimSpontaneousCa(unsigned int i);
		// Notifies the model that all values have been updated for timestep t
		virtual void UpdateVals(double t)
		{ 
			tCurr = t; 
			if (not isPreRunning)
			{
				metrics.ComputeMetrics<ChIModelNeedFrequentUpdateMetric>(*this);
				metrics.ComputeMetrics<ChIModelDynMetric>(*this);
				ODENetworkDynamicsModel<CouplingFunction, ChICell>::UpdateVals(t);
			}
		}

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
		// Getters                                                   ||
		//===========================================================||
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const;
		// Returns the total fluxes going out of the cell
		virtual double GetTotalFlux(unsigned int i) const;
		// Return the dynamic value that constitutes the excitable part of the system
		virtual double GetExcDynVal(unsigned int cellNb) const;

		// Returns a const ref on the network
		virtual const AbstractNetwork & GetNetwork() const 
			{return NetworkDynamicsModel<CouplingFunction>::GetNetwork();}
		// Returns tStart
		virtual const double & GetTStart() const 
			{ return ODE::ODEProblem<double, double>::GetTStart(); }
		// Return tEnd
		virtual const double & GetTEnd() const 
			{ return ODE::ODEProblem<double, double>::GetTEnd(); }
		// Returns the number of cells in the model
		inline unsigned int GetNbCells() const 
			{ return ODENetworkDynamicsModel<CouplingFunction, ChICell>::GetNbCells(); }
		// Set all cells to equilibrium
		virtual void SetAllCellsToEquilibrium()
			{ return ODENetworkDynamicsModel<CouplingFunction, ChICell>::SetAllCellsToEquilibrium(); }
		// Returns the neighbors of cell i
		virtual const std::vector<unsigned int> & GetNeighbors(unsigned int i) const
			{ return ODENetworkDynamicsModel<CouplingFunction, ChICell>::GetNeighbors(i); }
		// Is the given cell currently stimulated ?
		virtual bool IsStimulated(unsigned int ind) const
			{ return Stimulable::IsStimulated(ind); }


	protected:

		virtual double getModelVersionNum() const { return 0.2; }

		//===========================================================||
		// Network parameters                                        ||
		//===========================================================||
		double C0;  // Total cell free [Ca2+] per cytosolic volume
		double d1;  // IP3 dissociation constant
		double d2;  // Ca2+ inactivation dissociation constant
		double d3;  // IP3 dissociation constant
		double d5;  // Ca2+ activation dissociation constant
		double v3k; // Rate of IP3 degradation by 3K
		double K3k; // Half maximal degradation rate of IP3 by IP3-3K 
		double r5p; // Rate of IP3 degradation by IP-5P 

		//===========================================================||
		// Associated ODE Problem                                    ||
		//===========================================================||

		//===========================================================||
		// Metrics                                                   ||
		//===========================================================||
		SortedMetrics<ChIModel> metrics;
	};

/**********************************************************************/
/* ChI Model Threshold Determination                                  */
/**********************************************************************/
	class ChIModelThresholdDetermination : public ChIModel
	{
	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		ChIModelThresholdDetermination(ParamHandler & h);
		// Constructor from stream
		ChIModelThresholdDetermination(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams);
		// Destructor
		virtual ~ChIModelThresholdDetermination();

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Launch simulation
		virtual int Simulate(ResultSaver saver);

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Returns the network central node degree
		virtual unsigned int GetDegree() const;
		// Returns the number of derivations (sinks) for stim cells
		virtual unsigned int GetNbDerivs() const;

	protected:

		virtual double getModelVersionNum() const { return 1.2; }
	};

/**********************************************************************/
/* ChI Model Same Net Param Comb                                      */
/**********************************************************************/
	class ChIModelStochasticRes : public ChIModel
	{
	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		ChIModelStochasticRes(ParamHandler & h);
		// Constructor from stream
		ChIModelStochasticRes(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams);
		// Destructor
		virtual ~ChIModelStochasticRes();

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Launch simulation
		virtual int Simulate(ResultSaver saver);
		// Initialize Stimulation strategies 
		void InitializeStimStrats(bool firstInit = false);
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		std::string getSignalStimStartName() const { return signalStimStratClassName; }
		std::string getNoiseStimStartName() const { return noiseStimStratClassName; }

	protected:

		virtual double getModelVersionNum() const { return 1.3; }

		std::string signalStimStratClassName;
		std::string noiseStimStratClassName;
	};
}

#endif
