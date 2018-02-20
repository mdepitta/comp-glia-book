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

#ifndef PROPAGATIONMODEL_H
#define PROPAGATIONMODEL_H

#include <set>
#include "Model.h"
#include "CouplingFunction.h"
#include "ErrorCodes.h"

namespace AstroModel
{
/**********************************************************************/
/* Propagation model : Utility structures                             */
/**********************************************************************/
	template <typename T, int nbVals>
	class NamedStatesWrapper
	{
	public:
		NamedStatesWrapper(T v = T(), std::string names[nbVals] = 0, T values[nbVals] = 0) : val(v)
		{
			if (not namesInitialized and names and values and (nbVals > 0))
			{
				for (unsigned int i = 0 ; i < nbVals ; ++i)
				{
					stateNames.push_back(std::string(names[i]));
					stateVals.push_back(NamedStatesWrapper<T, nbVals>(values[i]));
				}
				namesInitialized = true;
			}
			name = GetName();
		}
		operator T() const { return val; }
		NamedStatesWrapper<T, nbVals> & operator = (T v) 
			{ val = v; name = ""; name = GetName(); return *this; }
		bool operator == (const NamedStatesWrapper<T, nbVals> & other)
			{ return val == other.val; }
		const std::string & GetName() const
		{
			if (name == "")
			{
				for (unsigned int i = 0 ; i < stateVals.size() ; ++i)
					if (stateVals[i] == *this)
						name = stateNames[i];
			}
			return name;
		}

		static const std::vector<NamedStatesWrapper<T, nbVals> > & GetStateVals() 
		{ return stateVals; }
		static const std::vector<std::string> & GetStateNames()
		{ return stateNames; }

	protected:
		static bool namesInitialized;
		static std::vector<std::string> stateNames;
		static std::vector<NamedStatesWrapper<T, nbVals> > stateVals;

		T val;
		mutable std::string name; // not guaranteed to have a value, call GetName().
	};

	template <typename T, int nbVals> bool NamedStatesWrapper<T, nbVals>::namesInitialized = false;
	template <typename T, int nbVals> std::vector<std::string> NamedStatesWrapper<T, nbVals>::stateNames;
	template <typename T, int nbVals> std::vector<NamedStatesWrapper<T, nbVals> > NamedStatesWrapper<T, nbVals>::stateVals;

	typedef NamedStatesWrapper<bool, 2> BooleanWrapper;

	enum SERStatesType {Refractory = 0, Susceptible = 1, Excited = 2};
	typedef NamedStatesWrapper<SERStatesType, 3> SERStates;
	
/**********************************************************************/
/* Propagation model abstract class                                   */
/**********************************************************************/
	template <typename StateType, typename LinkType>
	class PropagationModel : public NetworkDynamicsModel<LinkType>
	{
	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		PropagationModel(ParamHandler & h) : NetworkDynamicsModel<LinkType>::NetworkDynamicsModel(h)
		{
			start = h.getParam<unsigned int>("-PropagModelSim", 0);
			end   = h.getParam<unsigned int>("-PropagModelSim", 1);
		}
		// Constructor from stream
		PropagationModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams) :
			NetworkDynamicsModel<LinkType>::NetworkDynamicsModel(h)
		{
			PropagationModel<StateType, LinkType>::LoadFromStream(stream);
		}
		// Destructor
		virtual ~PropagationModel()
		{
			metrics.FreeAndClean();
		}

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Initializes the model
		virtual void Initialize(ResultSaver = ResultSaver::NullSaver)
		{
			NetworkDynamicsModel<LinkType>::Initialize();
			currStep = start;
			states = std::vector<StateType>(this->GetNetwork().size(), StateType());

			metrics.InitializeMetricsDefault();
		}

		// Launch simulation
		virtual int Simulate(ResultSaver saver)
		{
			int returnVal = 0;

			// Save wanted data
			if (this->DebugMode())
				std::cout << "*** Saving metric data ***" << std::endl;
			this->network->ComputeAndSaveMetrics(saver);

			// Simulate
			if (this->DebugMode())
				std::cout << "*** Starting Simulation ***" << std::endl;

			for (currStep = start ; (currStep <= end) and not this->isFinished() ; ++currStep)
			{
				unsigned int ind = this->chooseNode();
				states[ind] = this->getNewState(ind);
				// Metrics computations
				metrics.ComputeMetricsDefault(*this);
			}

			// Save computed dynamic data and metrics
			if (this->DebugMode())
				std::cout << "*** Saving dynamic data ***" << std::endl;
			if (not metrics.SaveMetricsDefault(saver))
				returnVal |= MODEL_METRIC_SAVING_PROBLEM;

			return returnVal;
		}
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = NetworkDynamicsModel<LinkType>::LoadFromStream(stream);
			stream >> start;
			stream >> end;
			metrics.FreeAndClean();
			ok &= metrics.LoadFromStream(stream);
			PropagationModel<StateType, LinkType>::Initialize();
			return ok and stream.good();
		}
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = NetworkDynamicsModel<LinkType>::SaveToStream(stream);
			stream 
				<< start << std::endl
				<< end << std::endl;
			ok &= metrics.SaveToStream(stream);
			return ok and stream.good();
		}

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler()
		{
			ParamHandler params;
			params += NetworkDynamicsModel<LinkType>::BuildModelParamHandler();
			params <= "PropagationModelStart", start;
			params <= "PropagationModelEnd", end;
			params += metrics.BuildModelParamHandler();
			return params;
		}

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a network metric to the network
		virtual bool AddMetric(Metric *_m, bool _f = false)
		{
			if (NetworkDynamicsModel<LinkType>::AddMetric(_m, _f))
				return true;
			else if (metrics.AddMetricAndDependencies(_m, _f, this))
				return true;
			else
				return false;
		}

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		inline unsigned int GetNbNodes() const { return states.size(); }
		inline const StateType & GetNodeState(unsigned int i) const
		{ return states[i]; }

		inline unsigned int GetCurrStep() const 
		{ return currStep; }

		virtual double GetTime() const 
		{ return ((double)(currStep - start)) / ((double) states.size()); }

		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const
		{
			std::vector<Metric *> mTot = NetworkDynamicsModel<LinkType>::GetAllMetrics();
			for (unsigned int i = 0 ; i < metrics.GetMetricsRaw().size() ; ++i)
				mTot.push_back(metrics.GetMetricsRaw()[i]);
			return mTot;
		}

		virtual unsigned int GetNbStimulatedNodes() const = 0;

	protected:
		std::vector<StateType> states;

		unsigned int start;
		unsigned int end;
		unsigned int currStep;

		SortedMetrics<PropagationModel<StateType, LinkType> > metrics;

		virtual unsigned int chooseNode() const = 0;
		virtual StateType getNewState(unsigned int ind) const = 0;
		virtual bool isFinished() const = 0;
	};

/**********************************************************************/
/* Threshold Model                                                    */
/**********************************************************************/
	class ThresholdModel : public PropagationModel<BooleanWrapper, BooleanLink>
	{
	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		ThresholdModel(ParamHandler & h) : PropagationModel<BooleanWrapper, BooleanLink>::PropagationModel(h)
		{
			threshold  = h.getParam<double>("-ThresholdMod", 0);
			stimRadius = h.getParam<int>("-ThresholdMod", 1);
			stimPoints = param.getParam<std::vector<int> >("-ModThreshStim", 0);
		}
		// Constructor from stream
		ThresholdModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams) :
			PropagationModel<BooleanWrapper, BooleanLink>::PropagationModel(h)
		{
			ThresholdModel::LoadFromStream(stream);
		}
		// Constructor from stream with global parameters
		ThresholdModel(std::ifstream & stream) : PropagationModel<BooleanWrapper, BooleanLink>::PropagationModel(stream) { }

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver);
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual unsigned int GetNbStimulatedNodes() const 
			{ return nbStimNodes; }

	protected:
		double threshold;
		int stimRadius;
		std::vector<int> stimPoints;
		unsigned int nbStimNodes;

		virtual unsigned int chooseNode() const;
		virtual BooleanWrapper getNewState(unsigned int ind) const;
		virtual bool isFinished() const;
	};

/**********************************************************************/
/* Modified Threshold Model                                           */
/**********************************************************************/
	class ModifiedThresholdModel : public ThresholdModel
	{
	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		ModifiedThresholdModel(ParamHandler & h) : ThresholdModel::ThresholdModel(h) 
		{
			propConst = h.getParam<double>("-ModThreshModel", 0);
			Bval = h.getParam<double>("-ModThreshModel", 1);
		}
		// Constructor from stream
		ModifiedThresholdModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams) :
			ThresholdModel::ThresholdModel(h) 
		{ 
			// Modify if a LoadFromStream method is added to ModifiedThresholdModel
			LoadFromStream(stream);
		}

		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

	protected:
		double propConst; // Proportionality constant
		double Bval;      // Constant threshold value

		virtual BooleanWrapper getNewState(unsigned int ind) const;
	};

/**********************************************************************/
/* SERS Model                                                         */
/**********************************************************************/
	class SERSModel : public PropagationModel<SERStates, BooleanLink>
	{
	public:
		//static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		SERSModel(ParamHandler & h);
		// Constructor from stream
		SERSModel(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams) :
			PropagationModel<SERStates, BooleanLink>::PropagationModel(h)
		{
			SERSModel::LoadFromStream(stream);
		}

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Initializes the model
		virtual void Initialize(ResultSaver saver = ResultSaver::NullSaver);
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual unsigned int GetNbStimulatedNodes() const 
			{ return allStimNodes.size(); }

	protected:
		double period;    // period of oscillation 
						  //  (time between successive excitations)
		double actTime;   // mean time during which a node stays in 
						  //  the excited state
		double transTime; // mean time needed in order to transmit 
						  //  an excitation

		int stimRadius;
		std::vector<int> stimPoints;

		std::set<unsigned int> allStimNodes;

		mutable std::vector<unsigned int> stimOrder;
		mutable unsigned int stimOrderInd;

		virtual unsigned int chooseNode() const;
		virtual bool isFinished() const;
		virtual SERStates getNewState(unsigned int ind) const;
		virtual double probaSusceptibleToExcited(unsigned int ind) const;
		virtual double probaExcitedToRefractory(unsigned int ind) const;
		virtual double probaRefractoryToSusceptible(unsigned int ind) const;
	};

/**********************************************************************/
/* Second order neighbor SERS Model                                   */
/**********************************************************************/
	class SecondOrdNeighbSERSModel : public SERSModel
	{
	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		SecondOrdNeighbSERSModel(ParamHandler & h);
		// Constructor from stream
		SecondOrdNeighbSERSModel(std::ifstream & stream, ParamHandler & h =
			ParamHandler::GlobalParams) : SERSModel::SERSModel(h)
		{
			SecondOrdNeighbSERSModel::LoadFromStream(stream);
		}

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

	protected:
		double A;         // prop coeff in chi = A*k + B
		double B;         // bias in same formula
		bool useModTransTime;
		double minTransTimeFact;

		virtual double probaSusceptibleToExcited(unsigned int ind) const;
	};

/**********************************************************************/
/* Simple threshold SERS Model                                        */
/*    Same model that Hütt et al. used in "Stochastic resonance in    */
/*    discrete excitable dynamics on graphs".                         */
/**********************************************************************/
	class SimpleThresholdSERSModel : public SERSModel
	{
	public:
		static std::string ClassName;
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		SimpleThresholdSERSModel(ParamHandler & h);
		// Constructor from stream
		SimpleThresholdSERSModel(std::ifstream & stream, ParamHandler & h =
			ParamHandler::GlobalParams) : SERSModel::SERSModel(h)
		{
			SimpleThresholdSERSModel::LoadFromStream(stream);
		}
		// Constructor from stream with global parameters
		SimpleThresholdSERSModel(std::ifstream & stream) : 
			SERSModel::SERSModel(stream) { }

		//===========================================================||
		// Standard model methods                                    ||
		//===========================================================||
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

	protected:
		double relativeThreshold;
		double spontaneousFiring;
		double recoveryProba;

		virtual double probaSusceptibleToExcited(unsigned int ind) const;
		virtual double probaExcitedToRefractory(unsigned int ind) const;
		virtual double probaRefractoryToSusceptible(unsigned int ind) const;
	};
}

#endif

