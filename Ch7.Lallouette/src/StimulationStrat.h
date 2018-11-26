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

#ifndef STIMULATIONSTRAT_H
#define STIMULATIONSTRAT_H

#include "Savable.h"
#include "ParamHandler.h"
#include "MetricComputeStrat.h"
#include "Network.h"
//#include "StimulationMetrics.h"

#include <string>

namespace AstroModel
{
	// Forward declarations
	class StimulationStrat;
	class NetworkEdges;

/**********************************************************************/
/* Stimulable base class                                              */
/**********************************************************************/
	class Stimulable
	{
	public:
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Constructor
		Stimulable(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		Stimulable(std::ifstream & stream, ParamHandler & h = ParamHandler::GlobalParams);
		// Destructor
		virtual ~Stimulable();

		//===========================================================||
		// Standard stimulable methods                               ||
		//===========================================================||
		// Initializes the stimulable object
		virtual void Initialize();
		// Stimulates cells using the stimulation strategy
		void Stimulate(double t);
		// Save Stim strat metrics
		bool SaveStimStratMetrics(ResultSaver saver) const;
		// Loads the model from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the model to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a network metric to the network
		virtual bool AddMetric(Metric *_m, bool _f = false);

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const;
		// Is the given cell currently stimulated ?
		virtual bool IsStimulated(unsigned int ind) const;
		// Returns true if the given stimulation strat is added to the model and activated
		virtual bool IsStimStratActivated(std::string) const;
		// Return tEnd
		virtual const double & GetTEnd() const = 0;
		// Return tStart
		virtual const double & GetTStart() const = 0;

	protected:
		// Activation mask for stimulation strategies
		std::vector<int> stimStratActMask;
		// Stimulation strategies
		std::vector<StimulationStrat *> stimStrats;
	};

/**********************************************************************/
/* Stimulable with cells                                              */
/**********************************************************************/
	class StimulableCellNetwork : public Stimulable
	{
	public:
		// Constructor
		StimulableCellNetwork(ParamHandler & h = ParamHandler::GlobalParams);

		// Returns the number of cells in the model
		virtual unsigned int GetNbCells() const = 0;
		// Returns a const ref on the network
		virtual const AbstractNetwork & GetNetwork() const = 0;
		// Returns the neighbors of cell i
		virtual const std::vector<unsigned int> & GetNeighbors(unsigned int i) const = 0;
		// Set all cells to equilibrium
		virtual void SetAllCellsToEquilibrium() = 0;
	};

/**********************************************************************/
/* Utility Stimulation Structure                                      */
/**********************************************************************/
	class StimulationInd
	{
	public:
		std::string stratName;
		unsigned int ind; // Cell indice
		double start;     // Stimulation start time
		double end;       // Stimulation end time
		unsigned int expand;
		unsigned int expandToDo;
		std::string origin;
		StimulationInd(std::string _n, unsigned int _i, double _s, double _e, 
			unsigned int _ex = 0, std::string _o = "Unknown") :
			stratName(_n), ind(_i), start(_s), end(_e), expand(_ex), 
			expandToDo(_ex), origin(_o) {}
	};

/**********************************************************************/
/* Abstract class                                                     */
/**********************************************************************/
	class StimulationStrat : public SaveAndLoadFromStream
	{
	public:
		// Destructor
		virtual ~StimulationStrat();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Strategy main method, calls the specialized protected 
		// stimulate method
		virtual void Stimulate(Stimulable & model, double t);
		// Returns the class name
		virtual std::string GetClassName() const = 0;
		// Initializes the stimulation strategy
		virtual void Initialize();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the cell from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the cell to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Add a metric
		virtual bool AddMetric(Metric *_m, bool _f = false) 
		{ return metrics.AddMetricAndDependencies(_m, _f, this); }
		// Save all metrics
		bool SaveMetrics(ResultSaver saver) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		inline bool IsStimulated(unsigned int ind) const 
			{ return (ind < stimulated.size()) ? stimulated[ind] : false; };
		inline unsigned int GetSize() const { return stimulated.size(); }
		inline double GetTime() const { return tCurr; }
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const;

	protected:
		double tCurr; // Current time
		// State of each cell (being stimulated or not)
		std::vector<bool> stimulated; 

		SortedMetrics<StimulationStrat> metrics;

		// Stimulates the given model
		virtual void stimulate(Stimulable & model) = 0;
	};

/**********************************************************************/
/* Network of cells simulation class                                  */
/**********************************************************************/
	class NetworkStimulationStrat : public StimulationStrat
	{

	protected:
		// Actually stimulate with appropriate method
		virtual void stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const = 0;
		// Stimulates the given model 
		virtual void stimulate(Stimulable & model);
		// Stimulates the given network model 
		virtual void stimulateNet(StimulableCellNetwork & model) = 0;
	};

/**********************************************************************/
/* Abstract Random stimulation class                                  */
/**********************************************************************/
	class RandomStimStrat : public NetworkStimulationStrat
	{
	public:
		// Initializes the stimulation strategy with the same random 
		// sequence that it already has.
		virtual void InitializeWithSameRandomSequence() = 0;
	};

/**********************************************************************/
/* Default Stimulation strategy                                       */
/* (list of defined cells and periodic square stimulations)           */
/**********************************************************************/
	class DefaultStimStrat : public NetworkStimulationStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		DefaultStimStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		DefaultStimStrat(double bias, CouplingFunction *_f);
		// Special constructor
		DefaultStimStrat(unsigned int ind, double _s, double _e,
			double bias = 1.0, CouplingFunction *_f = 0);
		// Constructor from stream
		DefaultStimStrat(std::ifstream & stream);
		// Destructor
		virtual ~DefaultStimStrat();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Loads the strategy from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the strategy to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;
		// Initializes the stimulation strategy
		virtual void Initialize();

		//===========================================================||
		// Specific Default Stimulation Strategy methods             ||
		//===========================================================||
		// Adds a stimulation order on cell 'CellNbr' 
		void AddCellStim(std::string name, unsigned int cellNbr, double tStart, 
			double tEnd, unsigned int exp = 0);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

	protected:
		//===========================================================||
		// Parameters and local variables                            ||
		//===========================================================||
		// Coupling function to bias cell for stimulation
		CouplingFunction *funct; 
		// Need to free the coupling function on death ?
		bool freeFunct;
		// IP3 bias concentration in bias cell
		double IP3Bias;
		// vector of predicted stimulations
		std::vector<StimulationInd> cellsToStim; 
		// unused stimulation points were removed
		bool unusedRemoved;

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Strategy main method, stimulates the given model
		virtual void stimulateNet(StimulableCellNetwork & model);
		// Actually stimulate with appropriate method
		virtual void stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const;
	};

/**********************************************************************/
/* Poissonian stimulation strategy                                    */
/* (random stimulation according to a poissonian distribution)        */
/**********************************************************************/
	class PoissonianStimStrat : public RandomStimStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		//Default Constructor from paramHandler
		PoissonianStimStrat(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Default Constructor
		PoissonianStimStrat(double _T, double _s, double bias, 
			CouplingFunction *_f, bool _oS, bool _cSR, double _ct,
			bool _in, bool _ugs, double _gqr, double _poc);
		// Constructor from stream
		PoissonianStimStrat(std::ifstream & stream);
		// Destructor
		virtual ~PoissonianStimStrat();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Loads the strategy from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the strategy to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;
		// Initializes the stimulation strategy
		virtual void Initialize();

		//===========================================================||
		// Specific Stimulation Strategy methods                     ||
		//===========================================================||
		// Initialize stimulation strategy but don't delete spikeTrain
		virtual void InitializeWithSameRandomSequence();

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

	protected:
		//===========================================================||
		// Parameters and local variables                            ||
		//===========================================================||
		double period;   // Mean period between two stimulations
		double stimTime; // Mean time of stimulation
		double IP3Bias;  // IP3 bias concentration in bias cell
		// If true, stimulates a minimum amount of time 
		// to make the cell spike
		bool onlyDriveToSpike; 
		bool useCaSpontRelease;
		double caThresh;
		bool isolateNodes;
		bool useGluStim;
		double gluQuantalRelease;
		double poissOmegaC;

		CouplingFunction *funct; // Coupling function to bias cell for stimulation
		bool freeFunct;          // Need to free the coupling function on death ?

		// Position in spike train for each cell [cell ind](spikeTrain ind)
		std::vector<unsigned int> spikeInd;

		// Spike Train [cell ind][spike Ind](start Time, end time)
		std::vector<std::vector<std::pair<double, double> > > spikeTrain;

		// counts the number of time each GJC has been shut down
		std::vector<std::vector<int> > gjcShutDown;

		// Check cell tags to determine whether a cell must be stimulated
		std::vector<unsigned long int> cellTags;

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Strategy main method, stimulates the given model
		virtual void stimulateNet(StimulableCellNetwork & model);
		// Actually stimulate with appropriate method
		virtual void stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const;
		// Get to next spike for a given cell
		virtual void stopStimAndSetNextSpike(unsigned int i);

		//===========================================================||
		// Utility mathematical functions                            ||
		//===========================================================||
		// Returns a random value from an exponential distribution
		double expDelay(double mean);
		// Returns a time during which a cell is going to stay activated
		double activationTime(double mean);
		// Generate a spike train
		void generateSpikeTrain(unsigned int nbCells, double tMax);

		// Isolate a node
		void isolateANode(unsigned int i, StimulableCellNetwork & model);
		// UnIsolate a node
		void unIsolateANode(unsigned int i, StimulableCellNetwork & model);

		// Parse cell tags from network
		void parseCellTags(StimulableCellNetwork & model);
	};

/**********************************************************************/
/* Random mono stimulation strategy                                   */
/* (One randomly chosen stimulated cell at a time)                    */
/**********************************************************************/
	class RandomMonoStim : public RandomStimStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		//Default Constructor from paramHandler
		RandomMonoStim(ParamHandler & h = ParamHandler::GlobalParams);
		// Default Constructor
		RandomMonoStim(double _sl, double _pl, double bias, 
			CouplingFunction *_f);
		// Constructor from stream
		RandomMonoStim(std::ifstream & stream);
		// Destructor
		virtual ~RandomMonoStim();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Loads the strategy from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the strategy to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;
		// Initializes the stimulation strategy
		virtual void Initialize();

		//===========================================================||
		// Specific Stimulation Strategy methods                     ||
		//===========================================================||
		// Initialize stimulation strategy but don't delete spikeTrain
		virtual void InitializeWithSameRandomSequence();

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

	protected:
		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Strategy main method, stimulates the given model
		virtual void stimulateNet(StimulableCellNetwork & model);
		// Actually stimulate with appropriate method
		virtual void stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const;
		// Generates stimulation train
		virtual void generateStimTrain(unsigned int nbCells, double tMax);

		//===========================================================||
		// Parameters and local variables                            ||
		//===========================================================||
		double stimLength;
		double pauseLength;
		double currStimStart;
		double currPauseStart;
		unsigned int currStimInd;

		double IP3Bias;  // IP3 bias concentration in bias cell

		CouplingFunction *funct; // Coupling function to bias cell for stimulation
		bool freeFunct;          // Need to free the coupling function on death ?

		std::vector<unsigned int> stimTrain;
	};

/**********************************************************************/
/* Threshold determination stimulation strategy                       */
/**********************************************************************/
	class ThresholdDeterminationStimStrat : public NetworkStimulationStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		//Default Constructor from paramHandler
		ThresholdDeterminationStimStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Default Constructor
		ThresholdDeterminationStimStrat(unsigned int ns, 
			unsigned int nb, double bias, 
			CouplingFunction *_f, double _Cat, bool _spontCa);
		// Constructor from stream
		ThresholdDeterminationStimStrat(std::ifstream & stream);
		// Destructor
		virtual ~ThresholdDeterminationStimStrat();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Loads the strategy from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the strategy to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;
		// Initializes the stimulation strategy
		virtual void Initialize();

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Specific Default Stimulation Strategy methods             ||
		//===========================================================||
		virtual void SetStimulatedCells(unsigned int nb)
		{ nbStimulated = nb; }

	protected:
		//===========================================================||
		// Parameters and local variables                            ||
		//===========================================================||
		unsigned int nbStimulated;
		unsigned int nbNeighbSinks;

		double IP3Bias;  // IP3 bias concentration in bias cell
		double CaThresh; // Calcium concentration at which stimulation is stopped
		bool useCaSpontRelease; // Use outflux from ER instead of IP3 bias

		std::vector<bool> neverSpiked;

		CouplingFunction *funct; // Coupling function to bias cell for stimulation
		bool freeFunct;          // Need to free the coupling function on death ?
		CouplingFunction *sinkFunct; // Coupling function to sink cell 

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Strategy main method, stimulates the given model
		virtual void stimulateNet(StimulableCellNetwork & model);
		// Useless here
		virtual void stimulateSpecific(StimulableCellNetwork *, unsigned int ) const {};
	};

/**********************************************************************/
/* Multiple Stimulation strategy                                      */
/*     List of randomly chosen stimulated cells separated by a given  */
/*     mean path length                                               */
/**********************************************************************/
	class MeanPathStimStrat : public DefaultStimStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		//Default Constructor from paramHandler
		MeanPathStimStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Default Constructor
		MeanPathStimStrat(unsigned int _ns, double _mpl, 
			double _ts,	double _te, unsigned int _se, 
			double bias, CouplingFunction *_f);
		// Constructor from stream
		MeanPathStimStrat(std::ifstream & stream);
		// Destructor
		virtual ~MeanPathStimStrat();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Loads the strategy from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the strategy to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;
		// Initializes the stimulation strategy
		virtual void Initialize();

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		double GetCurrEffMeanPathLength() const { return currEffMeanPath; }

	protected:
		//===========================================================||
		// Parameters and local variables                            ||
		//===========================================================||
		unsigned int nbStims;
		double meanPathLength;
		double tStart;
		double tEnd;
		unsigned int stimExp;

		bool stimCellsChosen;
		double currEffMeanPath;

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Strategy main method, stimulates the given model
		virtual void stimulateNet(StimulableCellNetwork & model);
	};

/**********************************************************************/
/* Correlation determination stimulation strategy                     */
/*     Stimulates all cells successively during a given amount of     */
/*     time and resets all cells to equilibrium when changing         */
/*     stimulation point                                              */
/**********************************************************************/
	class CorrelDetermStimStrat : public DefaultStimStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		//Default Constructor from paramHandler
		CorrelDetermStimStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Default Constructor
		CorrelDetermStimStrat(double _ps, double bias, CouplingFunction *_f);
		// Constructor from stream
		CorrelDetermStimStrat(std::ifstream & stream);
		// Destructor
		virtual ~CorrelDetermStimStrat();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Loads the strategy from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the strategy to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;
		// Initializes the stimulation strategy
		virtual void Initialize();

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||

	protected:
		//===========================================================||
		// Parameters and local variables                            ||
		//===========================================================||
		std::vector<double> resetTimes;
		unsigned int nextResetTimeInd;
		double pauseTime;

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Strategy main method, stimulates the given model
		virtual void stimulateNet(StimulableCellNetwork & model);
	};

/**********************************************************************/
/* Fixed extracellular glutamate stimulation                          */
/**********************************************************************/
	class ExtraCellGluStim : public DefaultStimStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		//Default Constructor from paramHandler
		ExtraCellGluStim(ParamHandler & h = ParamHandler::GlobalParams);
		// Default Constructor
		ExtraCellGluStim(double _glu, CouplingFunction *_f);
		// Constructor from stream
		ExtraCellGluStim(std::ifstream & stream);
		// Destructor
		virtual ~ExtraCellGluStim();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Loads the strategy from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the strategy to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||

	protected:
		//===========================================================||
		// Parameters and local variables                            ||
		//===========================================================||
		double glu;

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Actually stimulate with appropriate method
		virtual void stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const;
	};

/**********************************************************************/
/* Fixed extracellular glutamate stimulation                          */
/**********************************************************************/
	class ExtraCellKStim : public DefaultStimStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		//Default Constructor from paramHandler
		ExtraCellKStim(ParamHandler & h = ParamHandler::GlobalParams);
		// Default Constructor
		ExtraCellKStim(double _K, CouplingFunction *_f);
		// Constructor from stream
		ExtraCellKStim(std::ifstream & stream);
		// Destructor
		virtual ~ExtraCellKStim();

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Loads the strategy from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the strategy to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||

	protected:
		//===========================================================||
		// Parameters and local variables                            ||
		//===========================================================||
		double K; // Potassium level
		double incrRate; // Increase rate

		//===========================================================||
		// Standard Stimulation Strategy methods                     ||
		//===========================================================||
		// Actually stimulate with appropriate method
		virtual void stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const;
	};
}

#endif
