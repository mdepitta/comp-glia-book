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

#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <vector>
#include "ResultSaver.h"
#include "Savable.h"
#include "Network.h"
#include "SpatialStructureBuilder.h"
#include "SpatialNetwork.h"
#include "ParamHandler.h"
#include "MetricComputeStrat.h"
#include "ODESolvers.h"
#include "ErrorCodes.h"

namespace ODE
{
	// Forward declarations
	template <typename Val, typename TStep> class ODESolver;
}

namespace AstroModel
{
/**********************************************************************/
/* Abstract model class                                               */
/**********************************************************************/
	class Model : public SaveAndLoadFromStream
	{
	public:
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		Model(ParamHandler & h) : param(h) {}
		virtual ~Model() {}

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const = 0;
		// Initializes the model
		virtual void Initialize(ResultSaver saver) = 0;
		// Launch the model simulation
		virtual int Simulate(ResultSaver saver) = 0;
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			if (DebugMode())
				std::cout << "*** Loading Model From Stream ***" << std::endl;

			// File type version
			double versionNum;
			stream >> versionNum;
			if (versionNum != this->getModelVersionNum())
			{
				std::cout << "Current model version is " << this->getModelVersionNum()
					 << " but the file has been saved by version " << versionNum << std::endl;
			}
			return true;
		}
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			if (DebugMode())
				std::cout << "*** Saving Model To Stream ***" << std::endl;

			// File type version
			stream << this->getModelVersionNum() << std::endl;
			return true;
		}

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a network metric to the network
		virtual bool AddMetric(Metric *, bool)
		{
			return false;
		}

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler() = 0;

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Returns true if debug mode is on
		inline bool DebugMode() const { return param.getParam<bool>("-debug"); }
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const = 0;

	protected:
		virtual double getModelVersionNum() const { return 1.0; }

		// Param Handler
		ParamHandler & param;
	};

/**********************************************************************/
/* Network dynamics model                                             */
/**********************************************************************/
	template <typename NetworkEdges>
	class NetworkDynamicsModel : public Model
	{
	public:
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		NetworkDynamicsModel(ParamHandler & h = ParamHandler::GlobalParams) : Model::Model(h), network(0)
		{
			// Network
			if (not network)
			{
				if (param.getParam<bool>("-notSpacial"))
					network = new Network<NetworkEdges>(h);
				else
					network = new SpatialNetwork<NetworkEdges>(h);
			}
		}
		// Constructor from stream
		NetworkDynamicsModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams) :
			Model::Model(h), network(0)
		{
			if (not NetworkDynamicsModel<NetworkEdges>::LoadFromStream(stream))
				std::cerr << "Failed to load the model !" << std::endl;
		}
		// Destructor
		virtual ~NetworkDynamicsModel()
		{
			delete network;
		}
		
		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver)
		{
			// Initialize network
			network->BuildNetwork(saver);
			network->Initialize();
		}
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = Model::LoadFromStream(stream);
			if (ok)
			{
				std::string tempName;
				// Network topology
				if (network)
					delete network;
				stream >> tempName;
		TRACE(tempName)
				assert(AbstractFactory<AbstractNetwork>::Factories[tempName]);
				network = dynamic_cast<Network<NetworkEdges>* >(AbstractFactory<AbstractNetwork>::Factories[tempName]
								->CreateFromStream(stream));

				return true;
			}
			else
				return false;
		}
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = Model::SaveToStream(stream);
			// Network topology
			stream << network->GetClassName() << std::endl;
			ok &= network->SaveToStream(stream);
			return ok;
		}

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a network metric to the network
		virtual bool AddMetric(Metric *_m, bool _f = false)
		{
			if (Model::AddMetric(_m, _f))
				return true;
			else
				return (network ? network->AddMetric(_m, _f) : false);
		}

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler()
		{
			ParamHandler params;
			if (network)
				params += network->BuildModelParamHandler();
			return params;
		}

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Returns a const ref on the network
		inline const Network<NetworkEdges> & GetNetwork() const 
		{ return *network; }
		// Is the network topology defined ?
		virtual bool IsNetworkTopoDefined() const
		{
			return network and network->HasBeenBuilt();
		}
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const
		{
			std::vector<Metric *> mTot;
			if (network)
				mTot = network->GetAllMetrics();
			return mTot;
		}
		// Get spatial distance between two cells
		double GetDistance(unsigned int i, unsigned int j) const
		{
			assert((i < network->size()) and (j < network->size()));
			SpatialNetwork<NetworkEdges> *spatial = 
				dynamic_cast<SpatialNetwork<NetworkEdges>*>(network);
			if (spatial)
				return spatial->Pos(i)->DistanceTo(*(spatial->Pos(j)));
			else
				return 0;
		}
		// Returns the neighbors of cell i
		virtual const std::vector<unsigned int> & GetNeighbors(unsigned int i) const
			{ return network->GetNeighbors(i); }

	protected:
		// Network topology
		Network<NetworkEdges> *network; // Network structure
	};


	class AbstractODENetDynProblem :
		public ODE::ODEProblem<double, double>
	{
	public:
		AbstractODENetDynProblem(ParamHandler & h = ParamHandler::GlobalParams) :
			ODE::ODEProblem<double, double>::ODEProblem(0, h) {}

		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver)
		{
			ODE::ODEProblem<double, double>::Initialize(saver);
			// Initialize metrics
			metrics.InitializeMetricsDefault();
		}
		// Method called before solving the ODE problem
		virtual bool PreSimulationCall(ResultSaver )
		{
			bool ok = metrics.template ComputeMetrics<NeedFrequentUpdateMetric>(*this);
			return ok;
		}
		// Method called after solving the ODE problem
		virtual bool PostSimulationCall(ResultSaver saver)
		{
			bool ok = ODE::ODEProblem<double, double>::ComputeAfterSimMetrics();
			ok &= ODE::ODEProblem<double, double>::SaveMetrics(saver);

			ok &= metrics.template SaveMetrics<DynMetric>(saver);
			TRACE("AbstractODENetDynMetric saved.")

			ok &= metrics.template ComputeMetrics<AfterSimMetric>(*this);
			ok &= metrics.template SaveMetrics<AfterSimMetric>(saver);
			return ok;
		}
		virtual void UpdateVals(double t) 
		{
			ODE::ODEProblem<double, double>::UpdateVals(t);
			metrics.ComputeMetrics<NeedFrequentUpdateMetric>(*this);
		}
		// Only save metrics
		virtual bool SaveMetrics(ResultSaver saver)
		{
			return metrics.SaveMetricsDefault(saver);
		}
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = ODE::ODEProblem<double, double>::LoadFromStream(stream);
			// metrics
			metrics.FreeAndClean();
			ok &= metrics.LoadFromStream(stream);
			return ok and (stream.good() or stream.eof());
		}
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = ODE::ODEProblem<double, double>::SaveToStream(stream);
			// Metrics
			ok &= metrics.SaveToStream(stream);
			return ok;
		}

		// Add a network metric to the network
		virtual bool AddMetric(Metric *_m, bool _f = false)
		{
			if (ODE::ODEProblem<double, double>::AddMetric(_m, _f))
				return true;
			else
				return metrics.AddMetricAndDependencies(_m, _f, this);
		}
		virtual ParamHandler BuildModelParamHandler()
		{
			ParamHandler params;
			params += ODE::ODEProblem<double, double>::BuildModelParamHandler();
			params += metrics.BuildModelParamHandler();
			return params;
		}
		// Returns the total fluxes going out of the cell
		virtual double GetTotalFlux(unsigned int i) const = 0;
		// Return the dynamic value that constitutes the excitable part of the system
		virtual double GetExcDynVal(unsigned int cellNb) const = 0;
		// Is the given cell currently stimulated ?
		virtual bool IsStimulated(unsigned int ind) const = 0;
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const
		{
			std::vector<Metric *> mTot, mTemp;
			mTemp = ODE::ODEProblem<double, double>::GetAllMetrics();
			for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
				mTot.push_back(mTemp[i]);

			for (unsigned int i = 0 ; i < metrics.GetMetricsRaw().size() ; ++i)
				mTot.push_back(metrics.GetMetricsRaw()[i]);

			return mTot;
		}
		// Returns the neighbors of cell i
		virtual const std::vector<unsigned int> & GetNeighbors(unsigned int i) const = 0;
		// Returns the number of cells in the model
		virtual unsigned int GetNbCells() const = 0;

	protected:
		//===========================================================||
		// Metrics                                                   ||
		//===========================================================||
		SortedMetrics<AbstractODENetDynProblem> metrics;
	};

/**********************************************************************/
/* Network ODE dynamics model                                         */
/**********************************************************************/
	template <typename NetworkEdges, typename NodeType>
	class ODENetworkDynamicsModel : 
		public NetworkDynamicsModel<NetworkEdges>,
		public AbstractODENetDynProblem
	{
	public:
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		ODENetworkDynamicsModel(ParamHandler & h = ParamHandler::GlobalParams) : 
			NetworkDynamicsModel<NetworkEdges>::NetworkDynamicsModel(h), AbstractODENetDynProblem(h)//ODE::ODEProblem<double, double>::ODEProblem(0, h)
		{
			// Get Net size
			int nbCellsTemp = h.getParam<int>("-N");
			assert(nbCellsTemp > 0);
			unsigned int nbCells = nbCellsTemp;

			// Initialize ODEProblem
			SetUpCellsAndODEs(nbCells);
		}
		// Constructor from stream
		ODENetworkDynamicsModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams) :
			NetworkDynamicsModel<NetworkEdges>::NetworkDynamicsModel(h), AbstractODENetDynProblem(h)//ODE::ODEProblem<double, double>(h)
		{
			//if (not this->LoadFromStream(stream))
			if (not ODENetworkDynamicsModel<NetworkEdges, NodeType>::LoadFromStream(stream))
				std::cerr << "Failed to load the model !" << std::endl;
		}
		// Destructor
		virtual ~ODENetworkDynamicsModel()
		{
			// Cells
			for (unsigned int i = 0 ; i < cells.size() ; ++i)
				delete cells[i];
		}
		
		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const 
			{ return ODE::ODEProblem<double, double>::GetClassName(); }
		// Sets up cells and allocate data for ODEs
		virtual void SetUpCellsAndODEs(unsigned int nbCells)
		{
			//SetFunct(new ODE::ChINetworkFunct(*this), true);
			AllocateMemory(nbCells * NodeType::NbValsPerCell);

			for (unsigned int i = 0 ; i < cells.size() ; ++i)
				delete cells[i];
			cells.clear();
			for (unsigned int i = 0 ; i < nbCells ; ++i)
			{
				cells.push_back(new NodeType(this, this->vals + 
					i * NodeType::NbValsPerCell, false));
			}

			// Calls SetVals with 0 / 0 arguments (default args)
			// So that it doesn't change the allocated vals.
			SetVals();
		}
		// Change ODE vals to given pointer (for external use)
		virtual void SetVals(double *_vals = 0, unsigned int _nbVals = 0)
		{
			ODE::ODEProblem<double, double>::SetVals(_vals, _nbVals);

			for (unsigned int i = 0 ; i < cells.size() ; ++i)
			{
				cells[i]->SetDynVals(this->vals + i * NodeType::NbValsPerCell, false);
				cells[i]->SetValNamesPostfix(cells[i]->GetClassName() + "_" +
					StringifyFixed(i));
			}

			ODE::ODEProblem<double, double>::UseCurrentValsAsInitVals();
		}

		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver)
		{
			NetworkDynamicsModel<NetworkEdges>::Initialize(saver);
			AbstractODENetDynProblem::Initialize(saver);

			if (cells.size() != NetworkDynamicsModel<NetworkEdges>::network->size())
				SetUpCellsAndODEs(NetworkDynamicsModel<NetworkEdges>::network->size());

			InitializeDynVals();

			// Initialize metrics
			metrics.InitializeMetricsDefault();
		}
		// Initializes the dynamic values
		virtual void InitializeDynVals()
		{
			// Initialize cells
			for (unsigned int i = 0 ; i < cells.size() ; ++i)
				cells[i]->Initialize();
		}
		// Method called before solving the ODE problem
		virtual bool PreSimulationCall(ResultSaver saver)
		{
			bool ok = NetworkDynamicsModel<NetworkEdges>::network->ComputeAndSaveMetrics(saver);
			ok &= AbstractODENetDynProblem::PreSimulationCall(saver);
			ok &= metrics.template ComputeMetrics<NeedFrequentUpdateMetric>(*this);
			return ok;
		}
		// Method called after solving the ODE problem
		virtual bool PostSimulationCall(ResultSaver saver)
		{
		TRACE_UP("*** Computing and saving after simulation ODENetworkDynamicsModelMetrics ***")
			bool ok = AbstractODENetDynProblem::PostSimulationCall(saver);
			ok &= metrics.template ComputeMetrics<AfterSimMetric>(*this);
			ok &= metrics.template SaveMetrics<AfterSimMetric>(saver);
		TRACE_DOWN("*** After simulation metric data saved ***")
			return ok;
		}
		// Launch simulation
		virtual int Simulate(ResultSaver saver)
		{
			int returnVal = 0;

			PreSimulationCall(saver);

			// Simulate
		TRACE_UP("*** Starting Simulation tStart = " << tStart << " tEnd = " << tEnd << " ***")
			assert(solver);
			solver->Solve(*this, tStart, tEnd);
		TRACE_DOWN("*** Simulation Ended ***")

			if (not PostSimulationCall(saver))
				returnVal |= MODEL_METRIC_SAVING_PROBLEM;

			return returnVal;
		}

		virtual void UpdateVals(double t) 
		{
			AbstractODENetDynProblem::UpdateVals(t);
			metrics.template ComputeMetrics<NeedFrequentUpdateMetric>(*this);
		}

		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = NetworkDynamicsModel<NetworkEdges>::LoadFromStream(stream);
			ok &= AbstractODENetDynProblem::LoadFromStream(stream);
			ok &= LoadFunctFromStream(stream);
			ok &= LoadCellsFromStream(stream);
			// metrics
			metrics.FreeAndClean();
			ok &= metrics.LoadFromStream(stream);
			// Initialize ODEProblem
			SetUpCellsAndODEs(cells.size());
			return ok and (stream.good() or stream.eof());
		}
		// Load funct from stream
		virtual bool LoadFunctFromStream(std::ifstream & stream) = 0;
		// Load cells from stream
		virtual bool LoadCellsFromStream(std::ifstream & stream)
		{
			bool ok = true;
			// Cells
			for (unsigned int i = 0 ; i < cells.size() ; ++i)
				delete cells[i];
			cells.clear();

			unsigned int nbCells;
			stream >> nbCells;
			//
			AllocateMemory(nbCells * NodeType::NbValsPerCell);
			for (unsigned int i = 0 ; i < nbCells ; ++i)
			{
				cells.push_back(new NodeType(this, this->vals + i * NodeType::NbValsPerCell, false));
				ok &= cells.back()->LoadFromStream(stream);
			}
			return ok;
		}

		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = NetworkDynamicsModel<NetworkEdges>::SaveToStream(stream);
			ok &= AbstractODENetDynProblem::SaveToStream(stream);
			// Cells
			stream << cells.size() << std::endl;
			for (unsigned int i = 0 ; i < cells.size() ; ++i)
				ok &= cells[i]->SaveToStream(stream);
			// Metrics
			ok &= metrics.SaveToStream(stream);
			return ok;
		}


		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a network metric to the network
		virtual bool AddMetric(Metric *_m, bool _f = false)
		{
			if (NetworkDynamicsModel<NetworkEdges>::AddMetric(_m, _f))
				return true;
			else if (AbstractODENetDynProblem::AddMetric(_m, _f))
				return true;
			else
				return metrics.AddMetricAndDependencies(_m, _f, this);
		}

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler()
		{
			ParamHandler params;
			params += NetworkDynamicsModel<NetworkEdges>::BuildModelParamHandler();
			params += AbstractODENetDynProblem::BuildModelParamHandler();
			assert(not cells.empty());
			params += cells.back()->BuildModelParamHandler();
			params += metrics.BuildModelParamHandler();
			return params;
		}

		//===========================================================||
		// Setters and callbacks                                     ||
		//===========================================================||
		// Set all cells to equilibrium
		virtual void SetAllCellsToEquilibrium()
			{ solver->ResetVals(); }
		// Sets a specific cell to equilibrium
		virtual void SetCellToEquilibrium(unsigned int i)
			{ cells[i]->SetToEquilibrium(); }

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Returns a dynamic value of a cell
		virtual double GetDynVal(unsigned int cellNb, int name) const
			{ return cells[cellNb]->GetDynVal(name); }
		// Returns the number of dynamical values
		inline unsigned int GetNbDynVal(unsigned int ind) const 
			{ return cells[ind]->GetNbDynVals();}
		// Returns the number of cells in the model
		virtual unsigned int GetNbCells() const
			{ return cells.size(); }
		// Gives the total number of desired dyn vals (not equivalent to GetNbVals
		// from OPEProblem<double, double>
		virtual unsigned int GetTotNbDynVals() const
			{ return cells.size() * NodeType::NbValsPerCell; }
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const
		{
			std::vector<Metric *> mTot, mTemp;
			mTemp = NetworkDynamicsModel<NetworkEdges>::GetAllMetrics();
			for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
				mTot.push_back(mTemp[i]);

			mTemp = AbstractODENetDynProblem::GetAllMetrics();
			for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
				mTot.push_back(mTemp[i]);

			for (unsigned int i = 0 ; i < metrics.GetMetricsRaw().size() ; ++i)
				mTot.push_back(metrics.GetMetricsRaw()[i]);

			return mTot;
		}
		// Returns the total fluxes going out of the cell (should return 0 if irrelevant)
		virtual double GetTotalFlux(unsigned int i) const = 0;
		// Returns a constant ref on cell i
		virtual const NodeType & GetCell(unsigned int i) const
			{ return *cells[i]; }
		// Return the integration step
		virtual double GetIntegrStep() const
			{ return ODE::ODEProblem<double, double>::GetIntegrStep(); }
		// Returns current simulation time
		virtual double GetTime() const
			{ return ODE::ODEProblem<double, double>::GetTime(); }
		// Returns the neighbors of cell i
		virtual const std::vector<unsigned int> & GetNeighbors(unsigned int i) const
			{ return NetworkDynamicsModel<NetworkEdges>::GetNeighbors(i); }

	protected:
		// Cells
		std::vector<NodeType*> cells;

		//===========================================================||
		// Metrics                                                   ||
		//===========================================================||
		SortedMetrics<ODENetworkDynamicsModel<NetworkEdges, NodeType> > metrics;
	};
}

#endif
