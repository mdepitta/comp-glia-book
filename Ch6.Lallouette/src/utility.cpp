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

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include "utility.h"

using namespace std;

unsigned int NTRACE_GLOB = 0;
RngWrapper globalRNG;


RngWrapper::RngWrapper()
{
	static bool firstCall = true;
	if (firstCall)
	{
		gsl_rng_env_setup();
		firstCall = false;
	}
	rng = gsl_rng_alloc(gsl_rng_mt19937);
}

RngWrapper::~RngWrapper()
{
	gsl_rng_free(rng);
}

double RngWrapper::RandUnif()
{
	return gsl_rng_uniform(rng);
}

void RngWrapper::SetSeed(unsigned long int seed)
{
	gsl_rng_set(rng, seed);
}

std::string RepeatStr(std::string str, unsigned int nb)
{
	std::string tempStr;
	for (unsigned int i = 0 ; i < nb ; ++i)
		tempStr += str;
	return tempStr;
}

double GaussianRand(double mean, double sd)
{
	double fac, sqrRad, v1, v2;
	do
	{
		v1 = 2.0 * UnifRand() - 1.0;
		v2 = 2.0 * UnifRand() - 1.0;
		sqrRad = v1 * v1 + v2 * v2;
	} while ((sqrRad >= 1.0) or (sqrRad == 0.0));
	fac = sqrt(-2.0 * log(sqrRad) / sqrRad);
	return sd * v2 * fac + mean;
}

double UnifRand()
{
	return globalRNG.RandUnif();
}

double ExpRand(double mean)
{
	return gsl_ran_exponential(globalRNG.rng, mean);
}

double GammaRand(double a, double b)
{
	return gsl_ran_gamma(globalRNG.rng, a, b);
}

bool TrueWithProba(double p)
{
	return UnifRand() < p;
}

double ComputeMean(const std::vector<double> & vals)
{
	double temp = 0;
	for (unsigned int i = 0 ; i < vals.size() ; ++i)
		temp += vals[i];
	return temp / (double) vals.size();
}

double ComputeSum(const std::vector<double> & vals)
{
	double temp = 0;
	for (unsigned int i = 0 ; i < vals.size() ; ++i)
		temp += vals[i];
	return temp;
}

double ComputeVar(const std::vector<double> & vals, bool corrected)
{
	double temp = 0;
	double mean = ComputeMean(vals);
	for (unsigned int i = 0 ; i < vals.size() ; ++i)
		temp += pow(vals[i] - mean, 2.0);
	return temp / (double) (corrected ? vals.size() - 1 : vals.size());
}

double ComputeStdDev(const std::vector<double> & vals)
{
	return sqrt(ComputeVar(vals));
}

double ComputeMin(const std::vector<double> & vals)
{
	return vals.empty() ? 0 : vals[GetIndiceMin(vals)];
}

double ComputeMin(const std::vector<double> & vals, unsigned int iStart, unsigned int iEnd)
{
	return vals.empty() ? 0 : vals[GetIndiceMin(vals, iStart, iEnd)];
}

double ComputeMax(const std::vector<double> & vals)
{
	return vals.empty() ? 0 : vals[GetIndiceMax(vals)];
}

double ComputeMax(const std::vector<double> & vals, unsigned int iStart, unsigned int iEnd)
{
	return vals.empty() ? 0 : vals[GetIndiceMax(vals, iStart, iEnd)];
}

unsigned int GetIndiceMin(const std::vector<double> & vals)
{
	return GetIndiceMin(vals, 0, vals.size());
}

unsigned int GetIndiceMin(const std::vector<double> & vals, unsigned int iStart, unsigned int iEnd)
{
	if (vals.empty())
		return 0;
	double minVal = vals[iStart];
	double minInd = iStart;
	for (unsigned int i = (iStart+1) ; i < iEnd ; ++i)
	{
		if (vals[i] < minVal)
		{
			minVal = vals[i];
			minInd = i;
		}
	}
	return minInd;
}

unsigned int GetIndiceMax(const std::vector<double> & vals)
{
	return GetIndiceMax(vals, 0, vals.size());
}

unsigned int GetIndiceMax(const std::vector<double> & vals, unsigned int iStart, unsigned int iEnd)
{
	if (vals.empty())
		return 0;
	double maxVal = vals[iStart];
	double maxInd = iStart;
	for (unsigned int i = (iStart+1) ; i < iEnd ; ++i)
	{
		if (vals[i] > maxVal)
		{
			maxVal = vals[i];
			maxInd = i;
		}
	}
	return maxInd;
}

// Compute quantiles of a distribution
void ComputeAndAddQuantiles(std::map<std::string, double>& scalStats, std::string name, std::vector<double> & vect, const std::vector<double> & quants)
{
	std::sort(vect.begin(), vect.end());
	double *vals = new double[vect.size()];
	for (unsigned int i = 0 ; i < vect.size() ; ++i)
		vals[i] = vect[i];
	for (unsigned int i = 0 ; i < quants.size() ; ++i)
		scalStats[name + "_Quant_" + Stringify(quants[i])] = 
			gsl_stats_quantile_from_sorted_data (vals, 1, vect.size(), quants[i]);
	delete[] vals;
}

// Compute the frequency of events at given ratios
void ComputeAndAddRatioFreq(std::map<std::string, double>& scalStats, std::string name, std::vector<double> & vect, std::vector<double> ratios, double totalTime)
{
	double min = 0;
    double max = 0;
	if (not vect.empty())
	{
		std::sort(vect.begin(), vect.end());
		std::sort(ratios.begin(), ratios.end());
		min = vect[0];
		max = vect.back();
	}
	unsigned int currInd = 0;
	for(unsigned int i = 0 ; i < ratios.size() ; ++i)
	{
		while ((currInd < vect.size()) and (vect[currInd] < (ratios[i] * (max - min) + min)))
			++currInd;
		scalStats[name + "_Frequency_less_than_" + Stringify(ratios[i])] = ((double)currInd) / totalTime;
	}
	scalStats[name + "_MinVal"] = min;
	scalStats[name + "_MaxVal"] = max;
}

vector<vector<double> > LoadSparseMatrix(std::ifstream & stream, double fillVal)
{
	unsigned int m, n, ind;
	stream >> m;

	vector<vector<double> > mat(m, vector<double>(m, 0));

	for (unsigned int i = 0 ; i < m ; ++i)
	{
		stream >> n;
		for (unsigned int j = 0 ; j < n ; ++j)
		{
			stream >> ind;
			mat[i][ind] = fillVal;
			mat[ind][i] = mat[i][ind];
		}
	}

	assert(stream.good() and not stream.eof());
	return mat;
}

bool SaveSparseMatrix(std::ofstream & stream, const vector<vector<double> > & mat)
{
	stream << mat.size() << endl;
	for (unsigned int i = 0 ; i < mat.size() ; ++i)
	{
		vector<unsigned int> indices;
		for (unsigned int j = i ; j < mat[i].size() ; ++j)
		{
			if (mat[i][j])
				indices.push_back(j);
		}
		stream << indices.size() << endl;
		for (unsigned int j = 0 ; j < indices.size() ; ++j)
			stream << indices[j] << "\t";
		stream << endl;
	}
	return stream.good();
}

bool SaveFullMatrix(ofstream & stream, const vector<vector<double> > & mat)
{
	for (unsigned int i = 0 ; i < mat.size() ; ++i)
	{
		for (unsigned int j = 0 ; j < mat[i].size() ; ++j)
			stream << mat[i][j] << "\t";
		stream << endl;
	}
	return stream.good();
}

vector<vector<double> > LoadFullMatrix(ifstream & stream)
{
	vector<double> data;
	double val;
	stream >> val;
	do
	{
		data.push_back(val);
		stream >> val;
	} while (not stream.eof());

	unsigned int m = floor(sqrt(data.size()));
	vector<vector<double> > mat(m, std::vector<double>(m, 0));
	unsigned int ind = 0;
	for (unsigned int i = 0 ; i < m ; ++i)
		for (unsigned int j = 0 ; j < m ; ++j,++ind)
			mat[i][j] = data[ind];

	return mat;
}

void RandomSwapsOnVector(std::vector<unsigned int> & vect)
{
	for (unsigned int i = 0 ; i < vect.size() ; ++i)
	{
		unsigned int tempInd = floor(UnifRand() * (double)vect.size());
		unsigned int lastVal = vect[i];
		vect[i] = vect[tempInd];
		vect[tempInd] = lastVal;
	}
}

void ComputeAllPairDistsFromBidirAdjMat(const std::vector<std::vector<double> > & mat, std::vector<std::vector<double> > & dist, bool useW)
{
	dist = std::vector<std::vector<double> >(
		mat.size(), std::vector<double>(mat.size(), 0));
	for (unsigned int i = 0 ; i < mat.size() ; ++i)
		for (unsigned int j = i + 1 ; j < mat[i].size() ; ++j)
		{
			dist[i][j] = ((mat[i][j]) > 0 ? 
				(useW ? mat[i][j] : 1.0) : DEFAULT_MAX_PATH);
			dist[j][i] = dist[i][j];
		}

	for (unsigned int k = 0 ; k < mat.size() ; ++k)
		for (unsigned int i = 0 ; i < mat.size() ; ++i)
			for (unsigned int j = i + 1 ; j < mat.size() ; ++j)
			{
				dist[i][j] = std::min(dist[i][j], 
					dist[i][k]+dist[k][j]);
				dist[j][i] = dist[i][j];
			}
}

// COmpute efficiency given a distance matrix
double ComputeEfficiency(const std::vector<std::vector<double> > & dist)
{
	double tmpEff = 0;
	for (unsigned int i = 0 ; i < dist.size() ; ++i)
		for (unsigned int j = 0 ; j < dist[i].size() ; ++j)
			if (i != j)
			{
				if (dist[i][j] != DEFAULT_MAX_PATH)
					tmpEff += 1.0 / dist[i][j];
			}

	return tmpEff / ((double)dist.size() * ((double) dist.size() - 1.0));
}


//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns log(gamma) of real argument.
//  NOTE: Returns 1e308 if argument is 0 or negative.
//
double lgamma(double x)
{
    double x0,x2,xp,gl,gl0;
    int n, k;
	n = 0;
    static double a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590};
    
    x0 = x;
    if (x <= 0.0) return 1e308;
    else if ((x == 1.0) || (x == 2.0)) return 0.0;
    else if (x <= 7.0) {
		n = (int)(7-x);
        x0 = x+n;
    }
    x2 = 1.0/(x0*x0);
    xp = 2.0*M_PI;
    gl0 = a[9];
    for (k=8;k>=0;k--) {
        gl0 = gl0*x2 + a[k];
    }
    gl = gl0/x0+0.5*log(xp)+(x0-0.5)*log(x0)-x0;
    if (x <= 7.0) {
        for (k=1;k<=n;k++) {
            gl -= log(x0-1.0);
            x0 -= 1.0;
        }
    }
    return gl;
}


// Compute quick approximation to tanh
float fast_tanh(float x)
{
  float x2 = x * x;
  float a = x * (135135.0f + x2 * (17325.0f + x2 * (378.0f + x2)));
  float b = 135135.0f + x2 * (62370.0f + x2 * (3150.0f + x2 * 28.0f));
  float th = a / b;
  return (th < -1.0f) ? -1.0f : ((th > 1.0f) ? 1.0f : th);
}

//**********************************************************************
// Rank values in matrix
//**********************************************************************
std::vector<double> GetSortedValuesFromMat(
	const std::vector<std::vector<double> > & mat, bool avoidDiag)
{
	std::vector<double> sortVals;
	for (unsigned int i = 0 ; i < mat.size() ; ++i)
		for (unsigned int j = 0 ; j < mat[i].size() ; ++j)
			if (i != j or not avoidDiag)
				sortVals.push_back(mat[i][j]);
	std::sort(sortVals.begin(), sortVals.end());

	return sortVals;
}

//**********************************************************************
// Compute thresholded version of given matrix
//**********************************************************************
void ApplyThresholdOnMat(
	const std::vector<std::vector<double> > & mat, 
	std::vector<std::vector<double> > & thresh, double threshold)
{
	// Fill the thresholded topology
	thresh = std::vector<std::vector<double> >(
		mat.size(), std::vector<double>(mat.size(), 0));	
	for (unsigned int i = 0 ; i < mat.size() ; ++i)
		for (unsigned int j = 0 ; j < mat[i].size() ; ++j)
		{
			if (mat[i][j] >= threshold)
				thresh[i][j] = 1;
			else
				thresh[i][j] = 0;
		}
}

// Template specialization for std adj mat
template <> unsigned int CountEdgesInSubGraph(const std::vector<std::vector<double> > & network, const std::set<unsigned int> & inds)
{
	unsigned int edges = 0;
	for (std::set<unsigned int>::const_iterator it = inds.begin() ; it != inds.end() ; ++it)
	{
		for (std::set<unsigned int>::const_iterator it2 = inds.begin() ; it2 != inds.end() ; ++it2)
			if (network[*it][*it2] > 0)
				++edges;
	}
	return edges;
}

// Template specialization for std adj mat
template <> void AddNeighbors(const std::vector<std::vector<double> > & network, 
	unsigned int ind, std::set<unsigned int> & neighb)
{
	for (unsigned int i = 0 ; i < network[ind].size() ; ++i)
	{
		if (network[ind][i] > 0)
			neighb.insert(i);
	}
}

