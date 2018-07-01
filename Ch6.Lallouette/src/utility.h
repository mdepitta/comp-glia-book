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

#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <assert.h>
#include <set>
#include <gsl/gsl_rng.h>

#include "Savable.h"
#include "ParamHandler.h"

#define SYSTEM_NO_ERROR_CODE 0
#define NB_SPACE_TRACE 2
#define DEFAULT_MAX_PATH     999999
#define DEFAULT_MAX_VAL      999999
#define DEFAULT_SMALL_VAL    0.00000000000001
#define MEAN_PATH_PRECISION  0.0001

#define IF_DEBUG(msg) if (ParamHandler::GlobalParams.getParam<bool>("-debug", 0)) std::cout << msg << std::endl;

#define SPECIAL_TREE_FUNC(nb) RepeatStr(std::string("|") + std::string(NB_SPACE_TRACE, ' '), nb)

#define TRACE(msg) IF_DEBUG(SPECIAL_TREE_FUNC(NTRACE_GLOB) << __FILE__ << ":" << __LINE__ << " // " << msg)
#define TRACE_UP(msg) IF_DEBUG(SPECIAL_TREE_FUNC(NTRACE_GLOB++) << __FILE__ << ":" << __LINE__ << " // " << msg)
#define TRACE_DOWN(msg) IF_DEBUG(SPECIAL_TREE_FUNC(--NTRACE_GLOB) << __FILE__ << ":" << __LINE__ << " // " << msg)

extern unsigned int NTRACE_GLOB;

std::string RepeatStr(std::string str, unsigned int nb);

// Random number generator
class RngWrapper
{
public:
	RngWrapper();
	~RngWrapper();
	double RandUnif();
	void SetSeed(unsigned long int seed);
	gsl_rng *rng;
};

extern RngWrapper globalRNG;

class Position : public SaveAndLoadFromStream
{
public:
	//===========================================================||
	// Constructors / Destructor                                 ||
	//===========================================================||
	Position(unsigned int _dim = 1, double init = 0) : dim(_dim)
	{
		vals = new double[dim];
		for (unsigned int i = 0 ; i < dim ; ++i)
			vals[i] = init;
	}

	Position(const Position & pos) : dim(pos.dim)
	{
		vals = new double[dim];
		for (unsigned int i = 0 ; i < dim ; ++i)
			vals[i] = pos.vals[i];
	}

	Position(std::ifstream & stream) : dim(0), vals(0)
	{
		LoadFromStream(stream);
	}

	virtual ~Position()
	{
		if (vals)
			delete[] vals;
	}

	//===========================================================||
	// Standard Save and Load methods                            ||
	//===========================================================||
	virtual bool LoadFromStream(std::ifstream & stream)
	{
		stream >> dim;
		if (vals)
			delete[] vals;
		vals = new double [dim];
		for (unsigned int i = 0 ; i < dim ; ++i)
			stream >> vals[i];
		return not stream.eof();
	}

	virtual bool SaveToStream(std::ofstream & stream) const
	{
		stream << dim << std::endl;
		for (unsigned int i = 0 ; i < dim ; ++i)
			stream << vals[i] << std::endl;
		return stream.good();
	}

	//===========================================================||
	// Accessors                                                 ||
	//===========================================================||
	double & operator[](unsigned int i) 
	{
		static double voidVal;
		if (i < dim)
			return  vals[i];
		else
			return voidVal;
	}
	// Const version of operator[]
	const double & operator[](unsigned int i) const 
	{
		static double voidVal;
		if (i < dim)
			return  vals[i];
		else
			return voidVal;
	}

	unsigned int size() const { return dim; }

	//===========================================================||
	// Operators                                                 ||
	//===========================================================||
	Position operator+(const Position & pos) const
	{
		assert(dim == pos.dim);

		Position temp(dim, 0);
		for (unsigned int i = 0 ; i < dim ; ++i)
			temp.vals[i] = vals[i] + pos.vals[i];
		return temp;
	}

	Position operator-(const Position & pos) const
	{
		assert(dim == pos.dim);

		Position temp(dim, 0);
		for (unsigned int i = 0 ; i < dim ; ++i)
			temp.vals[i] = vals[i] - pos.vals[i];
		return temp;
	}

	double operator*(const Position & pos) const
	{
		assert(dim == pos.dim);

		double valTemp = 0;
		for (unsigned int i = 0 ; i < dim ; ++i)
			valTemp += vals[i] * pos.vals[i];
		return valTemp;
	}

	Position operator*(double scal) const
	{
		Position tmp(dim, 0);
		for (unsigned int i = 0 ; i < dim ; ++i)
			tmp.vals[i] = vals[i] * scal;
		return tmp;
	}

	Position & operator=(const Position & pos)
	{
		if (vals)
			delete[] vals;
		dim = pos.dim;
		vals = new double[dim];
		for (unsigned int i = 0 ; i < dim ; ++i)
			vals[i] = pos.vals[i];
		return *this;
	}

	Position abs() const
	{
		Position temp(dim, 0);
		for (unsigned int i = 0 ; i < dim ; ++i)
			temp[i] = fabs(vals[i]);
		return temp;
	}

	double norm() const { return sqrt(*this * *this); }

	double DistanceTo(const Position & pos) { return (*this - pos).norm(); }

protected:
	unsigned int dim;
	double *vals;
};

// Converts a value handled by std::ostringstream in a string
template <typename T> std::string Stringify(T o)
{
	std::ostringstream os;
	os << o;
	return(os.str());
}

// Converts a value handled by std::ostringstream in a string
template <typename T> std::string StringifyFixed(T o, unsigned int w = 5, char c = '0')
{
	std::ostringstream os;
	os << std::setfill(c) << std::setw(w) << o;
	return(os.str());
}

// Push back a value in a vector (std::vector<T>::push_back(const T&) is void)
template <typename T> std::vector<T> & DefaultVal(std::vector<T> &vect, T val)
{
	vect.push_back(val);
	return vect;
}

//template<typename BType, typename RType> RType * GetSpecificMetricNoFree(const std::vector<BType *> & metrics)
template<typename BType, typename RType> RType * GetSpecificMetric(std::vector<BType *> metrics)
{
	RType *tempPtr = 0;
	for (unsigned int i = 0 ; (i < metrics.size()) and not tempPtr ; ++i)
		tempPtr = dynamic_cast<RType*>(metrics[i]);
	return tempPtr;
}

template<typename BType, typename RType> std::vector<RType *> GetAllSpecificMetrics(std::vector<BType *> metrics)
{
	std::vector<RType *> result;
	RType *tempPtr = 0;
	for (unsigned int i = 0 ; i < metrics.size() ; ++i)
	{
		tempPtr = dynamic_cast<RType*>(metrics[i]);
		if (tempPtr)
			result.push_back(tempPtr);
	}
	return result;
}

double GaussianRand(double mean = 0.0, double sd = 1.0);
double UnifRand();
double ExpRand(double mean);
double GammaRand(double a, double b);
bool TrueWithProba(double p);
double lgamma(double x);
float fast_tanh(float x);

double ComputeMean(const std::vector<double> & vals);
double ComputeSum(const std::vector<double> & vals);
double ComputeVar(const std::vector<double> & vals, bool corrected = false);
double ComputeStdDev(const std::vector<double> & vals);
double ComputeMin(const std::vector<double> & vals);
double ComputeMin(const std::vector<double> & vals, unsigned int iStart, unsigned int iEnd);
double ComputeMax(const std::vector<double> & vals);
double ComputeMax(const std::vector<double> & vals, unsigned int iStart, unsigned int iEnd);

unsigned int GetIndiceMin(const std::vector<double> & vals);
unsigned int GetIndiceMin(const std::vector<double> & vals, unsigned int iStart, unsigned int iEnd);
unsigned int GetIndiceMax(const std::vector<double> & vals);
unsigned int GetIndiceMax(const std::vector<double> & vals, unsigned int iStart, unsigned int iEnd);

// Compute quantiles of a distribution
void ComputeAndAddQuantiles(std::map<std::string, double>& scalStats, std::string name, std::vector<double> & vect, const std::vector<double> & quants);
// Compute the frequency of events at given ratios
void ComputeAndAddRatioFreq(std::map<std::string, double>& scalStats, std::string name, std::vector<double> & vect, std::vector<double> ratios, double totalTime);

std::vector<std::vector<double> > LoadSparseMatrix(std::ifstream & stream, double fillVal = 1.0);
bool SaveSparseMatrix(std::ofstream & stream, const std::vector<std::vector<double> > & mat);

std::vector<std::vector<double> > LoadFullMatrix(std::ifstream & stream);
bool SaveFullMatrix(std::ofstream & stream, const std::vector<std::vector<double> > & mat);

void RandomSwapsOnVector(std::vector<unsigned int> & vect);

void ComputeAllPairDistsFromBidirAdjMat(const std::vector<std::vector<double> > & mat, std::vector<std::vector<double> > & dist, bool useW = false);

// Floyd Warschall
template <typename T> void ComputeAllPairDistances(const T & network, std::vector<std::vector<double> > & distances)
{
	distances = std::vector<std::vector<double> >(network.size(), std::vector<double>(network.size(), 0));
	for (unsigned int i = 0 ; i < network.size() ; ++i)
		for (unsigned int j = i + 1 ; j < network.size() ; ++j)
		{
			distances[i][j] = (network.AreConnected(i, j) ? 1.0 : DEFAULT_MAX_PATH);
			distances[j][i] = distances[i][j];
		}

	for (unsigned int k = 0 ; k < network.size() ; ++k)
		for (unsigned int i = 0 ; i < network.size() ; ++i)
			for (unsigned int j = i + 1 ; j < network.size() ; ++j)
			{
				distances[i][j] = std::min(distances[i][j], 
					distances[i][k]+distances[k][j] );
				distances[j][i] = distances[i][j];
			}
}

// BFS from each node
template <typename T> void ComputeAllPairDistancesSparse(const T & network, std::vector<std::vector<double> > & distances)
{
	distances = std::vector<std::vector<double> >(network.size(), std::vector<double>(network.size(), DEFAULT_MAX_PATH));
	std::set<unsigned int> border, old_border;
	unsigned int dist;
	for (unsigned int i = 0 ; i < distances.size() ; ++i)
	{
		old_border.clear();
		old_border.insert(i);
		dist = 1;
		do
		{
			border.clear();
			for (std::set<unsigned int>::const_iterator it = old_border.begin() ; it != old_border.end() ; ++it)
				for (std::vector<unsigned int>::const_iterator neighb = network.GetNeighbors(*it).begin() ; 
						neighb != network.GetNeighbors(*it).end() ; ++neighb)
					if (distances[i][*neighb] == DEFAULT_MAX_PATH)
					{
						border.insert(*neighb);
						distances[i][*neighb] = dist;
					}
			++dist;
			old_border = border;
		} while (border.size() > 0);
	}
}


// COmpute efficiency given a distance matrix
double ComputeEfficiency(const std::vector<std::vector<double> > & dist);

template <typename T> void AddNeighbors(const T & network, unsigned int ind, std::set<unsigned int> & neighb)
{
	for (unsigned int i = 0 ; i < network.size(ind) ; ++i)
	{
		if (network.AreConnected(ind, i))
			neighb.insert(i);
	}
}

template <> void AddNeighbors(const std::vector<std::vector<double> > & network, 
	unsigned int ind, std::set<unsigned int> & neighb);

template <typename T> unsigned int CountEdgesInSubGraph(const T & network, const std::set<unsigned int> & inds)
{
	unsigned int edges = 0;
	for (std::set<unsigned int>::const_iterator it = inds.begin() ; it != inds.end() ; ++it)
	{
		for (std::set<unsigned int>::const_iterator it2 = inds.begin() ; it2 != inds.end() ; ++it2)
			if (network.AreConnected(*it, *it2))
				++edges;
	}
	return edges;
}

template <> unsigned int CountEdgesInSubGraph(const std::vector<std::vector<double> > & network, const std::set<unsigned int> & inds);

template <typename T> double ComputeHierarchClustCoeff(const T & network, unsigned int s, unsigned int e, unsigned int repeat)
{
	std::vector<double> clustCoeffVals;
	for (unsigned int r = 0 ; r < repeat ; ++r)
	{
		unsigned int ind = floor(UnifRand() * (double)network.size());
		std::set<unsigned int> tempNeighb, intNeighb;
		tempNeighb.clear();
		intNeighb.clear();
		intNeighb.insert(ind);
		for (unsigned int d = 0 ; d < s - 1 ; ++d)
		{
			tempNeighb = intNeighb;
			for (std::set<unsigned int>::iterator it = tempNeighb.begin() ; it != tempNeighb.end() ; ++it)
				AddNeighbors(network, *it, intNeighb);
		}

		std::set<unsigned int> realNeighb = intNeighb;
		for (unsigned int d = s - 1 ; d < e ; ++d)
		{
			tempNeighb = realNeighb;
			for (std::set<unsigned int>::iterator it = tempNeighb.begin() ; it != tempNeighb.end() ; ++it)
				AddNeighbors(network, *it, realNeighb);
		}
			
		for (std::set<unsigned int>::iterator it = intNeighb.begin() ; it != intNeighb.end() ; ++it)
			realNeighb.erase(*it);

		unsigned int nbEdges = CountEdgesInSubGraph(network, realNeighb) / 2;

		if (realNeighb.size() > 1)
			clustCoeffVals.push_back(2 * (double)nbEdges / (double)(realNeighb.size() * (realNeighb.size() - 1)));
	}
	if (not clustCoeffVals.empty())
		return ComputeMean(clustCoeffVals);
	else
		return -1;
}

// Rank values in matrix
std::vector<double> GetSortedValuesFromMat(const std::vector<std::vector<double> > & mat, bool avoidDiag = true);

// Compute thresholded version of given matrix
void ApplyThresholdOnMat(const std::vector<std::vector<double> > & mat, 
	std::vector<std::vector<double> > & thresh, double threshold);

// Save a vector to a stream
template <typename T> bool SaveVectorToStream(const std::vector<T> & vect, std::ofstream & stream)
{
	stream << vect.size() << std::endl;
	for (unsigned int i = 0 ; i < vect.size() ; ++i)
		stream << vect[i] << std::endl;
	return stream.good();
}

// Load a vector from a stream
template <typename T> bool LoadVectorFromStream(std::vector<T> & vect, std::ifstream & stream)
{
	vect.clear();
	unsigned int tmpSize;
	stream >> tmpSize;
	T val;
	for (unsigned int i = 0 ; i < tmpSize ; ++i)
	{
		stream >> val;
		vect.push_back(val);
	}
	return stream.good();
}

#endif
