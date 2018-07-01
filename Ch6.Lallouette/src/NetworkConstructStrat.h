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

#ifndef NETWORKCONSTRUCTSTRAT_H
#define NETWORKCONSTRUCTSTRAT_H

#include "Savable.h"
#include "AbstractFactory.h"
#include "utility.h"
#include "ParamHandler.h"
#include "CouplingFunction.h"

#include <math.h>
#include <stdlib.h>

namespace AstroModel
{
	// Forward declarations
	class AbstractNetwork;
	class AbstractSpatialNetwork;
	class ChIModel;
	class Metric;

	// Abstract base class for construct strat that may not need
	// to be called for each run
	class StaticConstrStrat{
	public:
		virtual bool NeedsToBeconstructed() const = 0;
	};

/**********************************************************************/
/* Abstract class                                                     */
/**********************************************************************/
	class NetworkConstructStrat : public SaveAndLoadFromStream
	{
	public:
		// Returns the class name
		virtual std::string GetClassName() const = 0;

		// Default constructor
		NetworkConstructStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		NetworkConstructStrat(bool _norm) : normLkStr(_norm) {}

		//===========================================================||
		// Standard Network construction strategy method             ||
		//===========================================================||
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & fact, ResultSaver saver) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler() = 0;

		// Loads the construction strat from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the construction strat to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;
		
		// Returns a copy of the network construct strat
		virtual NetworkConstructStrat* BuildCopy() const;

		virtual std::vector<Metric *> GetAllMetrics() const;

		inline bool GetMakeDir() const { return makeDir; }

	protected:
		bool normLkStr;
		bool makeDir;
	};

/**********************************************************************/
/* Abstract Spatial class                                             */
/**********************************************************************/
	class SpatialNetConstrStrat : public NetworkConstructStrat
	{
	public:
		// Default constructor
		SpatialNetConstrStrat(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		SpatialNetConstrStrat(bool _n);

		//===========================================================||
		// Standard Network construction strategy method             ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;

	protected:
		//===========================================================||
		// Standard Spatial Network construction strategy methods    ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool buildSpatialNetwork(AbstractSpatialNetwork & network,
			AbstractFactory<NetworkEdge> & factory) const = 0;
	};

/**********************************************************************/
/* Scale free network construction strategy                           */
/**********************************************************************/
// Deprecated, use SpatialScaleFree with large rc instead
	class ScaleFreeStrat : public NetworkConstructStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		ScaleFreeStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		ScaleFreeStrat(double _g, bool _n);
		// Constructor from stream
		ScaleFreeStrat(std::ifstream & stream);

		//===========================================================||
		// Standard Network construction strategy methods            ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		double gamma; // Power law exponent
	};

/**********************************************************************/
/* Spatial Scale free network construction strategy                   */
/**********************************************************************/
	class SpatialScaleFreeStrat : public SpatialNetConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		SpatialScaleFreeStrat(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		SpatialScaleFreeStrat(double _rc, unsigned int _nl, bool _n);
		// Constructor from stream
		SpatialScaleFreeStrat(std::ifstream & stream);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		double rc;             // Spatial influence
		unsigned int newLinks; // Maximum number of new links
							   // for each added node

		//===========================================================||
		// Standard Spatial Network construction strategy methods    ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool buildSpatialNetwork(AbstractSpatialNetwork & network,
			AbstractFactory<NetworkEdge> & factory) const;
		// Compute raw link probabilities
		double LinkProbaRaw(AbstractSpatialNetwork & network, 
			unsigned int add, unsigned int exist, 
			bool considerUnconnected = false) const;
		// Normalize probabilities
		void NormalizeProbas(std::vector<double> & probas) const;
	};

/**********************************************************************/
/* Spatial regular network construction strategy                      */
/**********************************************************************/
	class SpatialRegularStrat : public SpatialNetConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		SpatialRegularStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		SpatialRegularStrat(unsigned int _d, double _mld, bool _n);
		// Constructor from stream
		SpatialRegularStrat(std::ifstream & stream);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		unsigned int degree; // Degree of each node
		double maxLinkDist;

		//===========================================================||
		// Standard Spatial Network construction strategy methods    ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool buildSpatialNetwork(AbstractSpatialNetwork & network,
			AbstractFactory<NetworkEdge> & factory) const;
	};

/**********************************************************************/
/* Spatial connection radius network construction strategy            */
/**********************************************************************/
	class SpatialConnectionRadiusStrat : public SpatialNetConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		SpatialConnectionRadiusStrat(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		SpatialConnectionRadiusStrat(double _r, bool _n);
		// Constructor from stream
		SpatialConnectionRadiusStrat(std::ifstream & stream);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		double radius; // Connection radius

		//===========================================================||
		// Standard Spatial Network construction strategy methods    ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool buildSpatialNetwork(AbstractSpatialNetwork & network,
			AbstractFactory<NetworkEdge> & factory) const;
	};

/**********************************************************************/
/* Construct the network from the voronoi diagram defined by the nodes*/
/**********************************************************************/
	class VoronoiConstructStrat : public SpatialNetConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		VoronoiConstructStrat(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		VoronoiConstructStrat(double _dmF, bool _ugp, double _mld);
		// Constructor from stream
		VoronoiConstructStrat(std::ifstream & stream);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		// Multiplicative factor for coordinated because the voronoi 
		// computation library is rounding of all coordinates
		double distMultFact;
		// If true, adds points on the sides to avoid longrange links
		// between nodes of the borders
		bool useGhostPoints;
		// Remove links if they are longer than this threshold
		double maxLinkDist;

		//===========================================================||
		// Standard Spatial Network construction strategy methods    ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool buildSpatialNetwork(AbstractSpatialNetwork & network,
			AbstractFactory<NetworkEdge> & factory) const;
	};

/**********************************************************************/
/* Erdos Renyi Random graph construction strategy                     */
/**********************************************************************/
	class ErdosRenyiRandomStrat : public NetworkConstructStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		ErdosRenyiRandomStrat(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		ErdosRenyiRandomStrat(double _k, bool _n);
		// Constructor from stream
		ErdosRenyiRandomStrat(std::ifstream & stream);

		//===========================================================||
		// Standard Network construction strategy methods            ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		double meanDegree; // Target mean degree of the network
	};

/**********************************************************************/
/* Small World network construction strategy                          */
/**********************************************************************/
	class SmallWorldStrat : public SpatialNetConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		SmallWorldStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		SmallWorldStrat(double _p, int _nD, bool _n);
		// Constructor from stream
		SmallWorldStrat(std::ifstream & stream);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		double prob; // Rewiring probability
		unsigned int neighbDist; // grid distance

		//===========================================================||
		// Standard Spatial Network construction strategy methods    ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool buildSpatialNetwork(AbstractSpatialNetwork & network,
			AbstractFactory<NetworkEdge> & factory) const;
	};

/**********************************************************************/
/* Erdos Renyi Random graph construction strategy                     */
/**********************************************************************/
	class ThresholdDeterminationStrat : public NetworkConstructStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		ThresholdDeterminationStrat(ParamHandler & h =
			ParamHandler::GlobalParams);
		// Special constructor
		ThresholdDeterminationStrat(unsigned int _k, unsigned int _d, 
			bool _n);
		// Constructor from stream
		ThresholdDeterminationStrat(std::ifstream & stream);

		//===========================================================||
		// Standard Network construction strategy methods            ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		// Returns the network central node degree
		virtual unsigned int GetDegree() const;
		// Returns the number of derivations (sinks) for stim cells
		virtual unsigned int GetNbDerivs() const;

	protected:
		unsigned int degree;  // degree of the observed node
		unsigned int nbDeriv; // number of outter cell per neighbor cells
	};

/**********************************************************************/
/* Fonctionnal topology construction strategy                         */
/**********************************************************************/
	class FunctionalTopoStrat : public NetworkConstructStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		FunctionalTopoStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		FunctionalTopoStrat(
			const std::vector<std::vector<double> > & _vals, bool _dL,
			bool _uT, double _t, double _sdC, bool _n = false, bool _uat = false,
			double _at = 0);
		// Constructor from stream
		FunctionalTopoStrat(std::ifstream & stream);

		//===========================================================||
		// Standard Network construction strategy methods            ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual void SetValues(const std::vector<std::vector<double> > & vals) const;
		virtual void SetThreshold(double thr, bool use = true);

	protected:
		mutable std::vector<std::vector<double> > values;

		bool dirLinks;
		bool useFixedThresh;
		bool useAbsoluteThresh;
		double threshold;
		double stdDevCoeff;
		double absoluteThresh;
	};

/**********************************************************************/
/* Fonctionnal topology construction strategy                         */
/**********************************************************************/
	class ConstrFromChISim : public NetworkConstructStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		ConstrFromChISim(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		ConstrFromChISim(const std::vector<std::string> _metr, bool _n = false);
		// Constructor from stream
		ConstrFromChISim(std::ifstream & stream);

		//===========================================================||
		// Standard Network construction strategy methods            ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;

		virtual std::vector<Metric *> GetAllMetrics() const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||

	protected:
		mutable ChIModel *model;

		std::string modelClassName;
		std::map<std::string, std::map<int, std::string> > params;

		std::vector<std::string> metricNames;
	};

/**********************************************************************/
/* Load from file construction strategy                               */
/**********************************************************************/
	class FromFileConstrStrat : public NetworkConstructStrat, public StaticConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		FromFileConstrStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		FromFileConstrStrat(std::string _p, bool _us, bool _d, bool _n = false, bool _uvas = false);
		// Constructor from stream
		FromFileConstrStrat(std::ifstream & stream);

		//===========================================================||
		// Standard Network construction strategy methods            ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;
		// StaticConstrStrat method, returns whether the network needs to be rebuilt
		virtual bool NeedsToBeconstructed() const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		mutable std::string path;
		mutable std::string newPath;
		mutable bool directed;
		bool newDirected;
		bool useSparse;
		bool useValuesAsStrengths;
	};

/**********************************************************************/
/* Load from file using file lists                                    */
/**********************************************************************/
	class FileListConstrStrat : public FromFileConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		FileListConstrStrat(ParamHandler & h = ParamHandler::GlobalParams);
		// Special constructor
		FileListConstrStrat(std::string _lp, bool _us, bool _d, bool _n = false);
		// Constructor from stream
		FileListConstrStrat(std::ifstream & stream);

		//===========================================================||
		// Standard Network construction strategy methods            ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;
		bool NeedsToBeconstructed() const;

		//===========================================================||
		// Specific methods                                          ||
		//===========================================================||
		// Load File List from given file
		void LoadFileListFromPath(std::string _lp) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		mutable std::string listPath;
		mutable std::string newListPath;

		mutable std::vector<std::string> fileList;
		mutable unsigned int posInList;
	};

/**********************************************************************/
/* Random Isolated cells construction strategy                        */
/*    Builds the network with another strat and isolates randomly     */
/*    chosen cells.                                                   */
/**********************************************************************/
	class RandomIsolatedCells : public NetworkConstructStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		RandomIsolatedCells(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		RandomIsolatedCells(double _rat, std::string _sN, bool _n);
		// Constructor from stream
		RandomIsolatedCells(std::ifstream & stream);

		//===========================================================||
		// Standard Network construction strategy methods            ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool BuildNetwork(AbstractNetwork & network, 
			AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		double ratio; // Mean ratio of isolated cells
		std::string subStratName; // Substrat class name
		mutable NetworkConstructStrat *subStrat; // Substrat
	};

/**********************************************************************/
/* Star-like astrocyte network wiring construction strategy           */
/**********************************************************************/
	class StarLikeAstroNetConstrStrat : public SpatialNetConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		StarLikeAstroNetConstrStrat(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		StarLikeAstroNetConstrStrat(unsigned int _sr, unsigned int _nbb, 
			double _scs, bool _n);
		// Constructor from stream
		StarLikeAstroNetConstrStrat(std::ifstream & stream);

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		unsigned int somaRadius; // Connection radius
		unsigned int nbBranches;
		double somaCouplStr;

		//===========================================================||
		// Standard Spatial Network construction strategy methods    ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool buildSpatialNetwork(AbstractSpatialNetwork & network,
			AbstractFactory<NetworkEdge> & factory) const;

		void setEdgeToSpecificCoupling(NetworkEdge *edge, double cpl) const;
	};

/**********************************************************************/
/* Shell structure scrambler network construct strategy               */
/*    Randomizes the shell structure and conserves the degree distrib */
/**********************************************************************/
	class ShellStructScramblerStrat : public SpatialNetConstrStrat
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		ShellStructScramblerStrat(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		ShellStructScramblerStrat(std::string _ss, int _sn, bool _n);
		// Constructor from stream
		ShellStructScramblerStrat(std::ifstream & stream);
		// Destructor
		virtual ~ShellStructScramblerStrat();

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

	protected:
		std::string subConstrStratName;
		int startNode;

		SpatialNetConstrStrat *strat;

		//===========================================================||
		// Standard Spatial Network construction strategy methods    ||
		//===========================================================||
		// Builds the network topology with given factory
		virtual bool buildSpatialNetwork(AbstractSpatialNetwork & network,
			AbstractFactory<NetworkEdge> & factory) const;

		virtual void randomizeSubSquare(
			const std::set<unsigned int> & pop1, 
			const std::set<unsigned int> & pop2, 
			std::vector<std::vector<double> > & adjMat ) const;

		// Returns the next shell
		virtual std::set<unsigned int> getNextShell(
			const std::set<unsigned int> & currShell, 
			const std::vector<std::vector<double> > & adjMat) const;

		// Return neighbors of a node
		virtual std::set<unsigned int> getNeighbs(unsigned int i, 
			const std::vector<std::vector<double> > & adjMat) const;

		// CLear all the links between two node populations
		void clearLinksBetween(const std::set<unsigned int> & pop1, 
			const std::set<unsigned int> & pop2, 
			std::vector<std::vector<double> > & adjMat) const;

		virtual void switchLines(unsigned int l1, unsigned int l2, 
			const std::set<unsigned int> & subcols, 
			std::vector<std::vector<double> > & adjMat) const;

		virtual void switchColumns(unsigned int c1, unsigned int c2, 
			const std::set<unsigned int> & sublines, 
			std::vector<std::vector<double> > & adjMat) const;

		virtual void switchCases(unsigned int l1, unsigned int c1, 
			unsigned int l2, unsigned int c2, 
			std::vector<std::vector<double> > & adjMat) const;
	};
}

#endif
