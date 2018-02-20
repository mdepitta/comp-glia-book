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

#ifndef CHIMODELMETRICS_H
#define CHIMODELMETRICS_H

#include "MetricComputeStrat.h"
#include "Network.h"
#include "ODEProblems.h"

#include <vector>
#include <queue>

// Needed to avoid conflicting declarations between wavelet and fft
#define GSL_DISABLE_DEPRECATED
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_fft_complex.h>

namespace AstroModel
{
	// Forward declarations
	class ChIModel;
	class AbstractNetwork;
	class AbstractODENetDynProblem;

	// Typedefs
	typedef SpecificMetric<ChIModel> ChIModelMetric;
	class ChIModelNeedFrequentUpdateMetric : 
		public ChIModelMetric, public NeedFrequentUpdateMetric {};
	class ChIModelDynMetric : public ChIModelMetric, public DynMetric {};
	class ODENetDynModelDynMetric : public SpecificMetric<AbstractODENetDynProblem>, public DynMetric {};
	class ChIModelStaticMetric : 
		public ChIModelMetric, public StaticMetric {};
	class ChIModelAfterSimMetric : 
		public ChIModelMetric, public AfterSimMetric {};
	class ChIModelThreshDetermMetric : public ChIModelNeedFrequentUpdateMetric {};

/**********************************************************************/
/* ChI Model Metrics : Utility structures                             */
/**********************************************************************/
	class SpikingCoherence
	{
	public:
		SpikingCoherence() : stdDevActTimes(-1), maxMinActTimesRatio(-1),
			minTime(__LONG_MAX__), maxTime(0) {}
		void AddNewSpike(unsigned int ind, double time)
		{
			if ((time > 0) and (std::find(alreadySpiked.begin(), alreadySpiked.end(), ind) == alreadySpiked.end()))
			{
				alreadySpiked.push_back(ind);
				spikeTimes.push_back(time);
				minTime = std::min(minTime, time);
				maxTime = std::max(maxTime, time);
				stdDevActTimes = ComputeStdDev(spikeTimes);
				maxMinActTimesRatio = minTime ? maxTime / minTime : -1;
			}
		}
		unsigned int GetNbActCels() const { return alreadySpiked.size(); }

		double stdDevActTimes;
		double maxMinActTimesRatio;
		double minTime;
		double maxTime;
	protected:
		std::vector<unsigned int> alreadySpiked;
		std::vector<double> spikeTimes;
	};

	class TimedActivations
	{
	public:
		TimedActivations(double _t) : time(_t) {}
		double time;
		std::vector<unsigned int> cells;
	};

	class TimedPropagation
	{
	public:
		TimedPropagation(double t, double sd, double isd, double cd, double cn) : 
			time(t), spatialDist(sd), instSpatialDist(isd), cellDist(cd), cellNb(cn) {}
		double time;
		double spatialDist;
		double instSpatialDist;
		double cellDist;
		unsigned int cellNb;
	};

	class CellToCellPropagation
	{
	public:
		CellToCellPropagation(unsigned int i, double t) : 
			ind(i), lastActiv(t), distToInitiator(0), localTopoSpeed(0), localSpatialSpeed(0) {}
		CellToCellPropagation(unsigned int i, unsigned int p, double d, double t) : 
			ind(i), lastActiv(t), distToInitiator(d), localTopoSpeed(0), localSpatialSpeed(0), 
			parents(std::vector<unsigned int>(1,p)) {}
		unsigned int ind;
		double lastActiv;
		double distToInitiator;
		double localTopoSpeed;
		double localSpatialSpeed;
		std::vector<unsigned int> parents;
		std::vector<unsigned int> children;
	};

	class TimedWave
	{
	public:
		TimedWave(unsigned int i, double t) : 
			tStart(t), initiator(i), 
			propagation(std::vector<TimedPropagation>(1, TimedPropagation(t, 0, 0, 0, 1))),
			structure(std::vector<CellToCellPropagation>(1, CellToCellPropagation(i, t)))
		{ involvedCells.insert(i); }
		void AddCell(CellToCellPropagation & prop)
		{
			structure.push_back(prop);
			involvedCells.insert(prop.ind);
			while (firstActCoherence.size() < (prop.distToInitiator + 1))
				firstActCoherence.push_back(SpikingCoherence());
			firstActCoherence[prop.distToInitiator].AddNewSpike(prop.ind, prop.lastActiv);
		}
		unsigned int GetCellNb() const { return involvedCells.size(); }

		double tStart;
		double tEnd;
		unsigned int initiator;
		std::vector<TimedPropagation> propagation;
		std::vector<CellToCellPropagation> structure;
		std::set<unsigned int> involvedCells;
		std::vector<SpikingCoherence> firstActCoherence;
	};

/**********************************************************************/
/* ChI Model Metrics : Activated cells                                */
/**********************************************************************/
	class ActivatedCells : public ODENetDynModelDynMetric//ChIModelDynMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		ActivatedCells(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		ActivatedCells(double _at, double _fetw, double _fcd);
		// Constructor from stream
		ActivatedCells(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		//virtual bool ComputeMetric(const ChIModel & model);
		virtual bool ComputeMetric(const AbstractODENetDynProblem & model);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		virtual ParamHandler BuildModelParamHandler();
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Compute and access data                                   ||
		//===========================================================||
		// Compute and return the time dependent frequency estimations
		virtual const std::vector<std::pair<double, double> > & 
			CompAndGetFreqEst() const;
		// Return the mean frequency estimation
		virtual double GetMeanFreqEst() const;
		// Return the mean influx in cells before firing 
		virtual double ComputeMeanInflux(unsigned int i) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		inline const std::vector<TimedActivations> & GetActivations() const
			{ return activations; }
		inline int GetCumulActCells() const 
			{ return cumulActivCells.size(); }
		inline double GetMaxInfluxInTimeWin(unsigned int i) const;
		inline double GetMeanInfluxAfterSpike(unsigned int i) const;

	protected:
		double activThresh;
		double freqEstTimeWin;
		double fluxComputingDelay;
		std::vector<TimedActivations> activations;
		// pair : time, frequency estimate
		mutable std::vector<std::pair<double, double> > freqEstim;
		mutable double meanFreqEstim;
		unsigned int totNbCells;

		std::set<unsigned int> cumulActivCells;
		std::vector<double> lastActivatedTime; // -1 for not activated recently
		std::vector<double> tmpOutFluxes;
		std::vector<std::vector<double> > tmpInFluxes;
		std::vector<std::pair<double, double> > totalFluxes; // (k, flux)
		std::vector<std::vector<double> > totalInFluxes; // [ind][spikeNbr]
		std::vector<std::pair<double, double> > maxInFluxInTimeWin; // (current, max)
		std::vector<std::queue<double> > neighborSpikeTimes; // [ind][FIFO spike times from neighbors]
		std::vector<std::vector<double> > influxAfterSpike; // [ind][event nbr]
		double lastTime;
		double tStart;
		double tEnd;
	};

/**********************************************************************/
/* ChI Model Metrics : Propagation distance                           */
/**********************************************************************/
	class PropagationDistance : public ChIModelDynMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		PropagationDistance(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		PropagationDistance(double _miD, double _maD, double _cdi, double _idtw, bool _amio);
		// Constructor from stream
		PropagationDistance(std::ifstream & stream);
		virtual ~PropagationDistance();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ChIModel & model);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		inline unsigned int size() const { return waves.size(); }
		const TimedWave & operator[](unsigned int ind) const;
		double GetMaxDelay() const { return maxDelay; }
	
	protected:
		unsigned int indActiv;
		double maxDelay;
		double minDelay;
		double instDistTimeWindow;
		double cellDistIncr;
		bool allMergeInOne;
		std::vector<std::pair<bool, double> > isInWaveSince;

		std::vector<TimedWave> waves;

		std::vector<double> startQuant;
		std::vector<double> endQuant;
		std::vector<double> stepQuant;

		std::vector<double> startRatio;
		std::vector<double> endRatio;
		std::vector<double> stepRatio;

		double lastTime;
	};

/**********************************************************************/
/* ChI Model Metrics : Concentrations of Ca2+, IP3 and h              */
/**********************************************************************/
	class ConcentrationMetrics : public DynMetric, public SpecificMetric<ODE::ODEProblem<double, double> >//public ChIModelDynMetric
	{
	protected:
		//typedef std::vector<std::vector<double> > CellsConcentrations;
		typedef std::pair<double, std::vector<double> > TimedConcentrations;

	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		ConcentrationMetrics(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		ConcentrationMetrics(double _ss);
		// Constructor from stream
		ConcentrationMetrics(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		//virtual bool ComputeMetric(const ChIModel & model);
		virtual bool ComputeMetric(const ODE::ODEProblem<double, double> & prob);
		virtual bool SaveMetric(ResultSaver saver) const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;
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
		const std::vector<TimedConcentrations> & GetConcentrations() const
		{ return concentrations; }
		const std::vector<std::string> & GetValNames() const
		{ return valNames; }
		double GetSavingStep() const { return savingStep; }

	protected:
		unsigned int nOut;
		double savingStep;
		mutable std::vector<TimedConcentrations> concentrations;
		mutable std::vector<std::string> valNames;
	};

/**********************************************************************/
/* Functional topology metric                                         */
/**********************************************************************/
	class FunctionalTopoMetric : public ChIModelAfterSimMetric
	{
		class ThreshInfo
		{
		public:
			ThreshInfo() : thresh(0), connec(0), 
				meanPath(0), ratUnco(0), falsePosRat(0), 
				truePosRat(0), hcc(0) {}
			double thresh;
			double connec;
			double meanPath;
			double ratUnco;
			double falsePosRat;
			double truePosRat;
			double hcc;
			double efficiency;
		};

		class simplifROCData
		{
		public:
			double falsePosRat;
			double truePosRat;
			double connec;
			double onLinksNb;
		};

	public:
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		FunctionalTopoMetric(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		FunctionalTopoMetric(int _m, double _t, double _sdC, double _mvmd, 
			bool _cl, bool _chcc, bool _umfb, bool _uSM, 
			const std::vector<std::string> & _fTML, VoronoiConstructStrat *_vcs = 0);
		// Constructor from stream
		FunctionalTopoMetric(std::ifstream & stream);
		// Destructor
		virtual ~FunctionalTopoMetric();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Functional topo specific methods                          ||
		//===========================================================||
		void ComputeFunctionalTopo(const ChIModel & model, 
			const std::string & name,
			const std::vector<std::vector<double> > & mat);
		const AbstractNetwork * GetFunctTopo(const std::string & name) const;

	protected:
		//===========================================================||
		// Parameters                                                ||
		//===========================================================||
		// 1 = same mean deg // 2 = specif mean Deg 
		// 3 = opt ROC       // 4 = auto thresh (mean + A * stddev)
		// 5 = based on highest Lfunc value
		// 6 = same mean deg as voronoi diagram
		// 7 = based on strongest decrease of L when graph is connected
		// 8 = based on voronoi and modulated by funct topo
		int threshEstimMethod;
		double threshSpecifMeanDeg; // meanDeg thresh
		double stdDevCoeff;
		double modVoroMinDecile;
		bool computeL;
		bool computeHCC;
		bool useMinForBidir;
		//
		bool functTopoUseSpecmetr;
		std::vector<std::string> functTopoMetrList;

		//===========================================================||
		// Utility objects                                           ||
		//===========================================================||
		VoronoiConstructStrat *voroConstrStrat;

		//===========================================================||
		// Data                                                      ||
		//===========================================================||
		// Functional topology
		std::map<std::string, Network<BooleanLink> *> functTopos;
		// <name> -> [thresh, connectivity, etc.]
		std::map<std::string, std::vector<ThreshInfo> > threshInf;
		// <name> -> indice in vector corersponding to chosen threshold
		std::map<std::string, unsigned int> chosenThresh;

		// threshold value for optimal point on ROC curve (first)
		// and distance to diagonal (second)
		std::map<std::string, std::pair<double, double> > optRocPoint;

		// Returns the indice of std::vector<ThreshInfo> corresponding
		// to given connectivity
		unsigned int getIndFromConnec(const std::string & name, 
			double connec) const;
		// Computes the ROC score of an adj mat
		simplifROCData computeROCPoint(const ChIModel & model, 
			const std::vector<std::vector<double> > &adjMat) const;

	};

/**********************************************************************/
/* Correlations metric                                                */
/**********************************************************************/
	class CorrelationsMetric : public FunctionalTopoMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		const static std::string ZeroCorrTopoName;
		const static std::string MaxCorrTopoName;
		const static std::string MaxCorrNoZeroTopoName;

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		CorrelationsMetric(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		CorrelationsMetric(int _m, double _t, double _sdC, double _mvmd, 
			bool _cl, bool _chcc, bool _umfb, bool _uSM, 
			const std::vector<std::string> & _fTML, double _mL);
		// Constructor from stream
		CorrelationsMetric(std::ifstream & stream);
		virtual ~CorrelationsMetric();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ChIModel & model);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		double GetMaxCorrMeanAbsLag() const
			{ return getMeanAbsLag(GetFunctTopo(MaxCorrTopoName),
				maxCrossCorrLags); }
		double GetMaxCorrNoZeroMeanAbsLag() const
			{ return getMeanAbsLag(GetFunctTopo(MaxCorrNoZeroTopoName),
				maxCrossCorrLagsNoZero); }

	protected:
		std::vector<std::vector<double> > zeroLagCorr;
		std::vector<std::vector<double> > maxCrossCorrVals;
		std::vector<std::vector<double> > maxCrossCorrLags;
		std::vector<std::vector<double> > maxCrossCorrValsNoZero;
		std::vector<std::vector<double> > maxCrossCorrLagsNoZero;

		std::vector<std::vector<double> > nbCommonNeighbs;
		std::vector<std::vector<double> > directionnality;

		double maxLag;

		// Get threshold to match connectivity
		double findThreholdForConnectivity(
			const std::vector<std::vector<double> > & corr, 
			unsigned int nbEdges) const;
		// Return the mean of the absolute values of lags for a given cross
		// correlation topology
		double getMeanAbsLag(
			const AbstractNetwork * topo, 
			const std::vector<std::vector<double> > & lags) const;
	};

/**********************************************************************/
/* Transfer Entropy Metric                                            */
/**********************************************************************/
	class TransferEntropyMetric : public FunctionalTopoMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		const static std::string TransEntrFunctTopoName;

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		TransferEntropyMetric(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		TransferEntropyMetric(int _m, double _t, double _sdC, double _mvmd, 
			bool _cl, bool _chcc, bool _umfb, bool _uSM, 
			const std::vector<std::string> & _fTML, double _te, 
			unsigned int _nb);
		// Constructor from stream
		TransferEntropyMetric(std::ifstream & stream);
		virtual ~TransferEntropyMetric();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ChIModel & model);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		const std::vector<std::vector<double> > GetValues() const;

	protected:
		double timeEmbed;
		unsigned int nbBins;
		unsigned int nbCarBins;

		std::vector<double> binsLimits; // Inclusive upper limits of bins
		
		// signal[t][cell](ind sigtable_k, ind_sigtable_kp1)
		std::vector<std::vector<std::pair<int, int> > > discrSig;
		// pairSig_for_i_and_j[t]
		//	(ind sigJointTable_xkt_yktp1, ind sigJointTable_xkp1tp1_yktp1)
		std::vector<std::pair<int, int> > tmpDiscrPairSig;

		std::vector<std::vector<double> > signalJointDistrib_kp1;
		std::vector<std::vector<double> > signalJointDistrib_k;
		std::vector<std::string> sigTable_kp1;
		std::vector<std::string> sigTable_k;

		std::vector<std::vector<std::vector<double> > > sigPairJointDistrib_xkp1tp1_yktp1;
		std::vector<std::vector<std::vector<double> > > sigPairJointDistrib_xkt_yktp1;
		// [ind sig pair](ind sigtable_kp1, ind sigtable_k)
		std::vector<std::pair<int, int> > sigJointTable_xkp1tp1_yktp1;
		// [ind sig pair](ind sigtable_k, ind sigtable_k)
		std::vector<std::pair<int, int> > sigJointTable_xkt_yktp1;


		std::vector<std::vector<double> > transferEntropy;

		void computeSignalBins(const std::vector<std::vector<double> > & vals);
		unsigned int binSignalValue(double val) const;
		inline std::string createSigStr(unsigned int k) const
		{ return std::string(k * nbCarBins, ' '); }
		void writeSigToStr(std::string & dest, unsigned int ind, double val) const;
		void writeSigToStr(std::string & dest, unsigned int ind, unsigned int bin) const;

		// Returns the indice of the first instance of val in vect
		template <typename T> unsigned int findInd(const std::vector<T> & vect, 
			const T & val) const
		{
			unsigned int ind;
			for (ind = 0 ; ind < vect.size() ; ++ind)
				if (vect[ind] == val)
					return ind;
			return ind;
		}
	};

/**********************************************************************/
/* Functional reconstruct by Net construct strat                      */
/**********************************************************************/
	class FunctTopoByNetConstrStrat : public FunctionalTopoMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		const static std::string FunctTopoByNetName;

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		FunctTopoByNetConstrStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		FunctTopoByNetConstrStrat(int _m, double _t, double _sdC, 
			double _mvmd, bool _cl, bool _chcc, bool _umfb, bool _uSM, 
			const std::vector<std::string> & _fTML, 
			std::vector<NetworkConstructStrat*> _ncs);
		// Constructor from stream
		FunctTopoByNetConstrStrat(std::ifstream & stream);
		virtual ~FunctTopoByNetConstrStrat();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ChIModel & model);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||

	protected:
		std::vector<NetworkConstructStrat *> netConstrStrats;
		std::vector<std::vector<std::vector<double> > > adjMats;

		std::vector<std::vector<double> > functTopoByNet;
	};

/**********************************************************************/
/* ChI Model Threshold determination Metrics : ThresholdDetermination */
/**********************************************************************/
	class ThresholdDetermination : public ChIModelThreshDetermMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		ThresholdDetermination();
		// Constructor from stream
		ThresholdDetermination(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ChIModel & model);
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
		bool HasSpiked() { return spiked.empty() ? false : spiked.back(); }
		double GetSecondNeighbThresh() const;
		double GetSecondNeighbThreshUnderTime(double time) const;

	protected:
		std::vector<bool> spiked;
		std::vector<double> timeSpiked;
		std::vector<double> nbStims;
		std::vector<double> nbNeighbs;
		std::vector<double> nbDerivs;

		double lastCompTime;
	};

/**********************************************************************/
/* ChI Model Metrics : Wave Front Detection and Stats Metric          */
/**********************************************************************/
	class WaveFrontDetect : public ChIModelAfterSimMetric
	{
		class NodeValues
		{
		public:
			NodeValues(unsigned int i, double d, double in, double ain, double dti) :
				ind(i), degree(d), influx(in), actualInflux(ain), distToInit(dti) {}
			unsigned int ind;
			double degree;
			double influx;
			double actualInflux;
			double distToInit;
		};

	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		WaveFrontDetect(ParamHandler & h = ParamHandler::GlobalParams);
		WaveFrontDetect(bool _upopans);
		// Constructor from stream
		WaveFrontDetect(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ChIModel & model);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		std::vector<double> GetFrontInluxes() const;
		std::vector<double> GetNonActInfluxes() const;

	protected:
		bool useParentsOfParentsAsNonSinks;

		std::vector<std::vector<NodeValues> > frontiers;
		std::vector<std::vector<NodeValues> > nonActNeighbs;
		std::map<double, std::vector<NodeValues> > shellStructure;
		std::vector<std::vector<double> > nodeNbInShell; // [shell][init point]
		std::vector<std::vector<double> > maxInfluxByShell; // [shell][cell] 

		double meanDegFront;
		double meanDegNonAct;
		double meanInfluxFront;
		double meanInfluxNonAct;
		double meanActualInfluxFront;
		double meanActualInfluxNonAct;
		double meanDistToInitFront;
		double meanDistToInitNonAct;

		double stdDevDegFront;
		double stdDevDegNonAct;
		double stdDevInfluxFront;
		double stdDevInfluxNonAct;
		double stdDevActualInfluxFront;
		double stdDevActualInfluxNonAct;
		double stdDevDistToInitFront;
		double stdDevDistToInitNonAct;

		// Compute the maximum theoretical influx received by the cell 
		// during a given wave.
		double computeMaxInflux(unsigned int ind, 
			const TimedWave & wave, const AbstractNetwork & net, 
			double maxDelay) const;
	};

/**********************************************************************/
/* ChI Model Metrics : Fourier Transform of all cell signals          */
/**********************************************************************/
	class FourierTransform : public AfterSimMetric, 
		public SpecificMetric<ODE::ODEProblem<double, double> >
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		FourierTransform(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		FourierTransform(std::ifstream & stream);
		// Destructor
		virtual ~FourierTransform();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ODE::ODEProblem<double, double> & prob);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||

	protected:
		bool concatenateSigs;
		bool useShortTimeFFT;
		double stfftWinSize;
		double stfftWinStep;
		double domFreqThreshRat;

		double minFreqToSave;
		double maxFreqToSave;
		double stepFreqToSave;

		// String elements allowing the categorization of signals
		std::map<std::string, std::vector<std::string> > elems;
		// Spectrums associated with signals (pairs (freq, weight))
		std::map<std::string, std::vector<std::pair<double, double> > > spectrums;

		// name -> [time][freq]
		std::map<std::string, std::vector<std::vector<double> > > spectrograms;

		std::map<std::string, std::vector<double> > dominantFrequencies;

		// Internal fourier transform data
		gsl_fft_complex_wavetable * wavetable;
		gsl_fft_complex_workspace * workspace;
		double *sigData_comp;
		unsigned int n;

		// Initilize spectrum with frequency values and put spectrum vals to 0
		void initializeSpectrum(const std::string & specName, unsigned int len, double dt);
		// Compute spectrum on the given signal delimited by iStart and iEnd, 
		// add the spectrum to spectrums
		void computeAndAddSpectrum(	
			const std::vector<double> & sig, unsigned int iStart, unsigned int iEnd, 
			const std::string & specName);
	};

/**********************************************************************/
/* ChI Model Metrics : Wavelet Transform of all cell signals          */
/**********************************************************************/
	class WaveletTransform : public AfterSimMetric, 
		public SpecificMetric<ODE::ODEProblem<double, double> >
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		WaveletTransform(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		WaveletTransform(std::ifstream & stream);
		// Destructor
		virtual ~WaveletTransform();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ODE::ODEProblem<double, double> & prob);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||

	protected:
		double domFreqThreshRat;

		double minFreqToSave;
		double maxFreqToSave;
		double stepFreqToSave;
		double powrThreshFract;

		// String elements allowing the categorization of signals
		std::map<std::string, std::vector<std::string> > elems;

		// name -> [cell][level](pseudoFreq, [time])
		std::map<std::string, std::vector<std::vector<
			std::pair<double, std::vector<double> > > > > spectrograms;

		std::map<std::string, std::vector<double> > dominantFrequencies;

		std::map<std::string, unsigned int> totNbSignals;

		// Internal fourier transform data
		gsl_wavelet *wavelet;
		gsl_wavelet_workspace *workspace;
		double *sigData;
		unsigned int n;
		unsigned int newN;
		unsigned int nbLevels;

		// Initilize spectrum with frequency values and put spectrum vals to 0
		//void initializeSpectrum(const std::string & specName, unsigned int len, double dt);

		// Compute spectrogram on the given signal delimited by iStart and iEnd and add it to spectrograms 
		void computeAndAddSpectrogram(	
			const std::vector<double> & sig, unsigned int iStart, unsigned int iEnd, 
			const std::string & specName, double samplingFreq);

		// Compute dominant frequencies on the given signal and indice in spectrogram
		// (to be called after computeAndAddSpectrogram)
		void computeAndAddDominantFrequencies(const std::string & specName, unsigned int ind);
	};

/**********************************************************************/
/* ChI Model Metrics : Propagation distance                           */
/**********************************************************************/
	class StochResMetric : public ChIModelMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		StochResMetric(ParamHandler & h = ParamHandler::GlobalParams);
		// Constructor from stream
		StochResMetric(std::ifstream & stream);
		virtual ~StochResMetric();

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const ChIModel & model);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;
	
	protected:
		double sigAndNoise_freq;
		double sigAndNoise_act;
		double noiseOnly_freq;
		double noiseOnly_act;
		double sigOnly_freq;
		double sigOnly_act;
	};
}

#endif
