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

#ifndef NETWORKMETRICS_H
#define NETWORKMETRICS_H

#include "MetricComputeStrat.h"
#include "utility.h"

namespace AstroModel
{
	// Forward declarations
	class AbstractNetwork;
	class AbstractSpatialNetwork;
	class CouplingFunction;

	// Typedefs
	typedef SpecificMetric<AbstractNetwork> NetworkMetric;
	typedef SpecificMetric<AbstractSpatialNetwork> SpatialNetMetric;

/**********************************************************************/
/* Network metrics : Degree Distribution                              */
/**********************************************************************/
	class DegreeDistrComp : public NetworkMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		DegreeDistrComp();
		// Constructor from stream
		DegreeDistrComp(std::ifstream & stream) { LoadFromStream(stream); }

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const AbstractNetwork & network);
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
		double operator[](unsigned int i) const;
		double GetMeanDegree() const { return meanDegree; }
		double GetStdDevDegree() const { return stdDevDegree; }
	
	protected:
		std::vector<double> degrees; // Cells degrees
		std::vector<double> totLinkStrengths; // Total link strength per node
		double meanDegree;
		double stdDevDegree;
		double noZeroMeanDegree;
		double noZeroStdDevDegree;
		double meanStrengths;
		double stdDevStrengths;
	};

/**********************************************************************/
/* Clustering coefficients                                            */
/**********************************************************************/
	class ClusteringCoeffs : public NetworkMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		ClusteringCoeffs(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		ClusteringCoeffs(unsigned int _s, unsigned int _e, unsigned int _r);
		// Constructor from stream
		ClusteringCoeffs(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const AbstractNetwork & network);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream &);
		virtual bool SaveToStream(std::ofstream &) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		double operator[](unsigned int i) const;
	
	protected:
		std::vector<double> clustCoeffs; // Cells degrees
		double meanClustCoeff;
		double stdDevClustCoeff;

		double estMeanHierarchClustCoeff; // Estimated Mean Hierarchical Clustering Coefficient
		unsigned int hccStart;
		unsigned int hccEnd;
		unsigned int hccRepeat;
	};

/**********************************************************************/
/* All pair distances                                                 */
/**********************************************************************/
	class AllPairDistances : public NetworkMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		AllPairDistances();
		// Constructor from stream
		AllPairDistances(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const AbstractNetwork & network);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;
		// Returns the scalar values to be saved by the scalar stats metric
		virtual std::map<std::string, double> GetScalarStatsToSave() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream &) { return true; }
		virtual bool SaveToStream(std::ofstream &) const { return true; }

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		const std::vector<double> & operator[](unsigned int i) const;
		inline const std::vector<std::vector<double> > GetDistances() const
		{ return distances; }
		double GetRatioUnconnected() const { return ratioUnconnectedPaths; }
	
	protected:
		std::vector<std::vector<double> > distances; // Cells degrees

		double meanShortestPath;
		double stdDevShortestPath;
		double ratioUnconnectedPaths;
	};

/**********************************************************************/
/* Network metrics : Adjacency Matrix                                 */
/**********************************************************************/
	class AdjMatrixComp : public NetworkMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		AdjMatrixComp();
		// Constructor from stream
		AdjMatrixComp(std::ifstream & stream);
		// Default destructor
		~AdjMatrixComp() {}

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const AbstractNetwork & network);
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
		// Get the adjacency matrix
		virtual const std::vector<std::vector<double> > & GetAdjMat() const;

	protected:
		std::vector<std::vector<double> > adjMat; // Adjacency matrix
	};

/**********************************************************************/
/* Spatial Network metrics : Positions                                */
/**********************************************************************/
	class PositionsComp : public SpatialNetMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		PositionsComp();
		// Constructor from stream
		PositionsComp(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		virtual bool ComputeMetric(const AbstractSpatialNetwork & network);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream & stream);
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Getters                                                   ||
		//===========================================================||
		inline double GetMeanDist() const { return meanDist; }
		inline double GetStdDevDist() const { return stdDevDist; }
		inline double GetMinDist() const { return minDist; }

	protected:
		std::vector<Position> positions;
		std::vector<double> distances;
		double meanDist;   // Mean min distances
		double stdDevDist; // Std dev min distances
		double minDist;    // Minimum min distance
	};

/**********************************************************************/
/* Network Dimensions                                                 */
/**********************************************************************/
	class NetworkDimensions : public NetworkMetric
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const {return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		NetworkDimensions(ParamHandler & h = ParamHandler::GlobalParams);
		// Full constructor
		NetworkDimensions(double repfract, std::vector<int> _nl, 
			bool _fne, int _mrtsas);
		// Constructor from stream
		NetworkDimensions(std::ifstream & stream);

		//===========================================================||
		// Standard Metric computing methods                         ||
		//===========================================================||
		std::map<std::string, double> GetScalarStatsToSave() const;
		virtual bool ComputeMetric(const AbstractNetwork & network);
		virtual bool SaveMetric(ResultSaver saver) const;
		virtual void Initialize();
		virtual Metric * BuildCopy() const;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		virtual bool LoadFromStream(std::ifstream &);
		virtual bool SaveToStream(std::ofstream &) const;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
	
	protected:
		double predN(double r) const;
		double predL(double r) const;
		double predC(double r) const;
		double predNRect(double r) const;
		double predLRect(double r) const;
		double predCRect(double r) const;

		double fractCompDim;
		std::vector<int> nodeList;
		bool fullNetExp;
		unsigned int maxRadiusToSaveAsStats;

		// Computed data
		std::vector<std::pair<double, double> > rDist;      // rStart and rEnd
		std::vector<std::pair<double, double> > stdDim;     // d  and sigma(d)
		std::vector<std::pair<double, double> > stdPrefact; // C  and sigma(C)
		std::vector<std::pair<double, double> > LDim;       // d' and sigma(d')
		std::vector<std::pair<double, double> > LPrefact;   // F  and sigma(D)
		std::vector<std::pair<double, double> > CDim;       // d' and sigma(d')
		std::vector<std::pair<double, double> > CPrefact;   // E  and sigma(D)
		std::vector<std::pair<double, double> > lgDim;      // d' and sigma(d')
		std::vector<std::pair<double, double> > lgPrefact;  // D  and sigma(D)
		std::vector<std::pair<double, double> > zetaEst;    // zeta = 2F / (C(k-2))

		// Averaged data
		std::vector<double> nbNodeTot;
		std::vector<double> nbOutLinksTot;
		std::vector<double> nbIntraLinksTot;
		std::vector<double> nbPoints;

		// Final Values for dimensions and prefactors
		double d;
		double C;
		double dprime;
		double E;
		double F;
		double D;
		double zeta;

		// Parameters for local scaling effects correction
		// log(V(r) / (C*r^d)) = a * r^-b
		double ad;
		double bd;
		double aL;
		double bL;
		double aC;
		double bC;
	};

}

#endif
