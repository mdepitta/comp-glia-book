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

#include "NetworkMetrics.h"
#include "Network.h"
#include "SpatialNetwork.h"
#include "CouplingFunction.h"

#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_sf_exp.h>
#include <algorithm>
#include "MetricNames.h"

using namespace AstroModel;
using namespace std;

//********************************************************************//
//************** D E G R E E   D I S T R I B U T I O N ***************//
//********************************************************************//

string DegreeDistrComp::ClassName(DEGREE_DISTRIBUTION_NETWORK_METRIC);
string ClusteringCoeffs::ClassName(CLUSTERING_COEFFICIENT_DISTRIBUTION_NETWORK_METRIC);
string AllPairDistances::ClassName(ALL_PAIR_DISTANCES_NETWORK_METRIC);
string AdjMatrixComp::ClassName(ADJACENCY_MATRIX_NETWORK_METRIC);
string PositionsComp::ClassName(SPATIAL_POSITION_NETWORK_METRIC);
string NetworkDimensions::ClassName(DIMENSION_NETWORK_METRIC);

//**********************************************************************
// Default Constructor
//**********************************************************************
DegreeDistrComp::DegreeDistrComp()
{

}

//**********************************************************************
// Compute the degree distribution
//**********************************************************************
bool DegreeDistrComp::ComputeMetric(const AbstractNetwork & network)
{
	degrees = vector<double>(network.size(), 0);
	totLinkStrengths = vector<double>(network.size(), 0);
	CouplingFunction *cplFct = 0;

	for(unsigned int i = 0 ; i < network.size() ; ++i)
		for (unsigned int j = i + 1 ; j < network.size(i) ; ++j)
			if (network.AreConnected(i,j))
			{
				++degrees[i];
				++degrees[j];
				// If edges have strengths
				if ((cplFct = dynamic_cast<CouplingFunction*>(network.GetAbstractEdge(i, j))))
				{
					totLinkStrengths[i] += cplFct->GetStrength();
					totLinkStrengths[j] += cplFct->GetStrength();
				}
			}
	meanDegree = ComputeMean(degrees);
	stdDevDegree = ComputeStdDev(degrees);

	meanStrengths = ComputeMean(totLinkStrengths);
	stdDevStrengths = ComputeStdDev(totLinkStrengths);

	// Computing mean and stdev without considering unconnected nodes 
	std::vector<double> noZeroDegrees;
	for (unsigned int i = 0 ; i < network.size() ; ++i)
		if (degrees[i] > 0)
			noZeroDegrees.push_back(degrees[i]);
	noZeroMeanDegree = ComputeMean(noZeroDegrees);
	noZeroStdDevDegree = ComputeStdDev(noZeroDegrees);

	return true;
}

//**********************************************************************
// Save the degree distribution
//**********************************************************************
bool DegreeDistrComp::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string degreeDistrName("DegreeDistr");
	std::string strengthsDistrName("StrengthsDistr");
	std::string degreeDistrStatsName("DegreeDistrStats");

	if (saver.isSaving(degreeDistrName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(degreeDistrName, saver.getCurrFile());

		for (unsigned int i = 0 ; i < degrees.size() ; ++i)
			stream << degrees[i] << endl;

		allSaved &= stream.good();
	}
	if (saver.isSaving(strengthsDistrName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(strengthsDistrName, saver.getCurrFile());

		for (unsigned int i = 0 ; i < totLinkStrengths.size() ; ++i)
			stream << totLinkStrengths[i] << endl;

		allSaved &= stream.good();
	}
	if (saver.isSaving(degreeDistrStatsName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(degreeDistrStatsName, saver.getCurrFile());

		stream << "MeanDegree\tStdDevDegree\t"
			<< "NoZeroMeanDegree\tNoZeroStdDevDegree\t"
			<< "MeanStrengths\tStdDevStrengths" << std::endl;
		stream 
			<< meanDegree << "\t" 
			<< stdDevDegree << "\t" 
			<< noZeroMeanDegree << "\t" 
			<< noZeroStdDevDegree << "\t"
			<< meanStrengths << "\t"
			<< stdDevStrengths << std::endl;

		allSaved &= stream.good();
	}
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void DegreeDistrComp::Initialize()
{
	Metric::Initialize();
	degrees.clear();
	totLinkStrengths.clear();
	meanDegree = 0;
	stdDevDegree = 0;
	noZeroMeanDegree = 0;
	noZeroStdDevDegree = 0;
	meanStrengths = 0;
	stdDevStrengths = 0;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * DegreeDistrComp::BuildCopy() const
{
	return new DegreeDistrComp(); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> DegreeDistrComp::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	res["MeanDegree"] = meanDegree;
	res["StdDevDegree"] = stdDevDegree;
	res["NoZeroMeanDegree"] = noZeroMeanDegree;
	res["NoZeroStdDevDegree"] = noZeroStdDevDegree;
	res["MeanStrengths"] = meanStrengths;
	res["StdDevStrengths"] = stdDevStrengths;
	return res;
}

//**********************************************************************
// Returns the degree of the ith cell
//**********************************************************************
double DegreeDistrComp::operator[](unsigned int i) const
{
	assert(i < degrees.size());
	return degrees[i];
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool DegreeDistrComp::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool DegreeDistrComp::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//********************************************************************//
//*********** C L U S T E R I N G   C O E F F I C I E N T ************//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
ClusteringCoeffs::ClusteringCoeffs(ParamHandler & h)
{
	hccStart = h.getParam<unsigned int>("-HCCParams", 0);
	hccEnd = h.getParam<unsigned int>("-HCCParams", 1);
	hccRepeat = h.getParam<unsigned int>("-HCCParams", 2);
}

//**********************************************************************
// Full constructor
//**********************************************************************
ClusteringCoeffs::ClusteringCoeffs(unsigned int _s, unsigned int _e, 
	unsigned int _r) :
		hccStart(_s), hccEnd(_e), hccRepeat(_r)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ClusteringCoeffs::ClusteringCoeffs(std::ifstream & stream) 
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Compute the clustering coefficient distribution
//**********************************************************************
bool ClusteringCoeffs::ComputeMetric(const AbstractNetwork & network)
{
	clustCoeffs = vector<double>(network.size(), 0);
	for(unsigned int i = 0 ; i < network.size() ; ++i)
	{
		vector<unsigned int> neighbors;
		for (unsigned int j = 0 ; j < network.size(i) ; ++j)
		{
			if (network.AreConnected(i, j))
				neighbors.push_back(j);
		}
		double nbLinks = 0;
		for (unsigned int k = 0 ; k < neighbors.size() ; ++k)
		{
			for (unsigned m = k + 1 ; m < neighbors.size() ; ++m)
				if (network.AreConnected(neighbors[k], neighbors[m]))
					++nbLinks;
		}
		if (nbLinks > 0)
			clustCoeffs[i] = 2.0 * nbLinks / 
				((double)neighbors.size() * (neighbors.size() - 1.0));
		else
			clustCoeffs[i] = 0.0;
	}

	// Hierarchical Clust Coeffs :
	estMeanHierarchClustCoeff = ComputeHierarchClustCoeff(
		network, hccStart, hccEnd, hccRepeat);

	meanClustCoeff = ComputeMean(clustCoeffs);
	stdDevClustCoeff = ComputeStdDev(clustCoeffs);

	return true;
}

//**********************************************************************
// Save the clustering coefficients distribution
//**********************************************************************
bool ClusteringCoeffs::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string clustCoeffDistrName("ClustCoeffDistr");
	std::string hierarchClustCoeffName("HierarchClustCoeff");

	if (saver.isSaving(clustCoeffDistrName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(clustCoeffDistrName, saver.getCurrFile());

		for (unsigned int i = 0 ; i < clustCoeffs.size() ; ++i)
			stream << clustCoeffs[i] << endl;

		allSaved &= stream.good();
	}
	if (saver.isSaving(hierarchClustCoeffName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(hierarchClustCoeffName, saver.getCurrFile());

		stream << "Start\tEnd\tHierarchClustCoeff" << endl;
		stream 
			<< hccStart << "\t"
			<< hccEnd   << "\t"
			<< estMeanHierarchClustCoeff << endl;

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void ClusteringCoeffs::Initialize()
{
	Metric::Initialize();
	clustCoeffs.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * ClusteringCoeffs::BuildCopy() const
{
	return new ClusteringCoeffs(hccStart, hccEnd, hccRepeat); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> ClusteringCoeffs::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	res["MeanClustCoeff"] = meanClustCoeff;
	res["StdDevClustCoeff"] = stdDevClustCoeff;
	std::string nameHierarchClust = "EstHierarchClustCoeff_";
	nameHierarchClust += Stringify(hccStart) + "_" + Stringify(hccEnd);
	res[nameHierarchClust] = estMeanHierarchClustCoeff;
	return res;
}

//**********************************************************************
// Returns the degree of the ith cell
//**********************************************************************
double ClusteringCoeffs::operator[](unsigned int i) const
{
	assert(i < clustCoeffs.size());
	return clustCoeffs[i];
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool ClusteringCoeffs::LoadFromStream(std::ifstream & stream)
{
	stream >> hccStart;
	stream >> hccEnd;
	stream >> hccRepeat;
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool ClusteringCoeffs::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< hccStart  << endl
		<< hccEnd    << endl
		<< hccRepeat << endl;
	return stream.good();
}

//********************************************************************//
//*************** A L L   P A I R   D I S T A N C E S ****************//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
AllPairDistances::AllPairDistances() 
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
AllPairDistances::AllPairDistances(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Compute the all pair distances
//**********************************************************************
bool AllPairDistances::ComputeMetric(const AbstractNetwork & network)
{
	//ComputeAllPairDistances(network, distances);
	ComputeAllPairDistancesSparse(network, distances);

	std::vector<double> allPairDists;
	ratioUnconnectedPaths = 0;
	for (unsigned int i = 0 ; i < distances.size() ; ++i)
		for (unsigned int j = i + 1 ; j < distances.size() ; ++j)
		{
			if (distances[i][j] != DEFAULT_MAX_PATH)
				allPairDists.push_back(distances[i][j]);
			else
				++ratioUnconnectedPaths;
		}

	ratioUnconnectedPaths /= (distances.size() - 1) * distances.size() / 2.0;
	meanShortestPath = ComputeMean(allPairDists);
	stdDevShortestPath = ComputeStdDev(allPairDists);

	return true;
}

//**********************************************************************
// Save the clustering coefficients distribution
//**********************************************************************
bool AllPairDistances::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string fullAllPairDistName("FullAllPairDist");
	std::string pathLengthInfosName("PathLengthInfos");

	if (saver.isSaving(fullAllPairDistName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(fullAllPairDistName, saver.getCurrFile());

		SaveFullMatrix(stream, distances);

		allSaved &= stream.good();
	}
	if (saver.isSaving(pathLengthInfosName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(pathLengthInfosName, saver.getCurrFile());

		stream << "AvgPathLength\tStdDevPathLength\tRatioUnconnected" << std::endl;

		stream 
			<< meanShortestPath << "\t" 
			<< stdDevShortestPath << "\t"
			<< ratioUnconnectedPaths << std::endl;

		allSaved &= stream.good();
	}
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void AllPairDistances::Initialize()
{
	Metric::Initialize();
	distances.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * AllPairDistances::BuildCopy() const
{
	return new AllPairDistances(); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> AllPairDistances::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	res["MeanShortestPath"] = meanShortestPath;
	res["StdDevShortestPath"] = stdDevShortestPath;
	res["RatioUnconnectedPaths"] = ratioUnconnectedPaths;
	return res;
}

//**********************************************************************
// Returns the distances from the ith cell
//**********************************************************************
const std::vector<double> & AllPairDistances::operator[](unsigned int i) const
{
	assert(i < distances.size());
	return distances[i];
}


//********************************************************************//
//***************** A D J A C E N C Y   M A T R I X ******************//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
AdjMatrixComp::AdjMatrixComp() 
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
AdjMatrixComp::AdjMatrixComp(std::ifstream & stream)
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Compute the degree distribution
//**********************************************************************
bool AdjMatrixComp::ComputeMetric(const AbstractNetwork & network)
{
	adjMat = vector<vector<double> >(network.size(), vector<double>(network.size(), 0));
	for (unsigned int i = 0 ; i < network.size() ; ++i)
		for (unsigned int j = 0 ; j < network.size(i) ; ++j)
			if (network.AreConnected(i, j))
				adjMat[i][j] = 1;
	return true;
}

//**********************************************************************
// Save the Adjacency matrix   
//**********************************************************************
bool AdjMatrixComp::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string fulladjacencyMatrixName("FullAdjacencyMatrix");
	std::string sparseAdjacencyMatrixname("SparseAdjacencyMatrix");

	if (saver.isSaving(fulladjacencyMatrixName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(fulladjacencyMatrixName, saver.getCurrFile());

		SaveFullMatrix(stream, adjMat);

		allSaved &= stream.good();
	}

	if (saver.isSaving(sparseAdjacencyMatrixname))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(sparseAdjacencyMatrixname, saver.getCurrFile());

		SaveSparseMatrix(stream, adjMat);

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void AdjMatrixComp::Initialize()
{
	Metric::Initialize();
	adjMat.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * AdjMatrixComp::BuildCopy() const
{
	return new AdjMatrixComp(); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool AdjMatrixComp::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool AdjMatrixComp::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//**********************************************************************
// Get the adjacency matrix
//**********************************************************************
const std::vector<std::vector<double> > & AdjMatrixComp::GetAdjMat() const
{
	return adjMat;
}

//********************************************************************//
//**************** S P A T I A L   P O S I T I O N S *****************//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
PositionsComp::PositionsComp() 
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
PositionsComp::PositionsComp(std::ifstream & stream) 
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Compute positions, mean distances from cell to cell
//**********************************************************************
bool PositionsComp::ComputeMetric(const AbstractSpatialNetwork & network)
{
	// Positions
	positions = vector<Position>(network.size(), Position());
	for (unsigned int i = 0 ; i < network.size() ; ++i)
		positions[i] = *(network.Pos(i));
	// Distances
	distances.clear();
	double minD;
	for (unsigned int i = 0 ; i < positions.size() ; ++i)
	{
		minD = 999999;
		for (unsigned int j = 1 ; j < positions.size() ; ++j)
			if ((i != j) and (minD > positions[i].DistanceTo(positions[j])))
				minD = positions[i].DistanceTo(positions[j]);
		distances.push_back(minD);
	}
	// Statistics
	meanDist = ComputeMean(distances);
	stdDevDist = ComputeStdDev(distances);
	minDist = ComputeMin(distances);

	return true;
}

//**********************************************************************
// Save positions and cell to cell distances
//**********************************************************************
bool PositionsComp::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string positionsName("Positions");
	std::string cellToCellDistName("CellToCellDist");

	if (saver.isSaving(positionsName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(positionsName, saver.getCurrFile());

		for (unsigned int i = 0 ; i < positions.size() ; ++i)
		{
			for (unsigned int j = 0 ; j < positions[i].size() ; ++j)
			{
				stream << positions[i][j] << "\t";
			}
			stream << endl;
		}

		allSaved &= stream.good();
	}
	if (saver.isSaving(cellToCellDistName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(cellToCellDistName, saver.getCurrFile());

		for (unsigned int i = 0 ; i < distances.size() ; ++i)
			stream << distances[i] << endl;

		allSaved &= stream.good();
	}
	saver.Param("MeanDist", meanDist);
	saver.Param("StdDevDist", stdDevDist);
	saver.Param("MinDist", minDist);
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void PositionsComp::Initialize()
{
	Metric::Initialize();
	positions.clear();
	distances.clear();
	meanDist = 0;
	stdDevDist = 0;
	minDist = 0;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * PositionsComp::BuildCopy() const
{
	return new PositionsComp(); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool PositionsComp::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool PositionsComp::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//********************************************************************//
//**************** N E T W O R K   D I M E N S I O N S ***************//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(DIMENSION_NETWORK_METRIC, DEGREE_DISTRIBUTION_NETWORK_METRIC)

//**********************************************************************
// Default constructor
//**********************************************************************
NetworkDimensions::NetworkDimensions(ParamHandler & h)
{
	fractCompDim = h.getParam<double>("-NetDimParams", 0);
	maxRadiusToSaveAsStats = h.getParam<int>("-NetDimParams", 1);
	nodeList = h.getParam<std::vector<int> >("-NetDimNodeList", 0);
	fullNetExp = h.getParam<bool>("-NetDimFullNet", 0);
}

//**********************************************************************
// Full constructor
//**********************************************************************
NetworkDimensions::NetworkDimensions(double repfract, std::vector<int> _nl,
	bool _fne, int _mrtsas) : fractCompDim(repfract), nodeList(_nl), 
	fullNetExp(_fne), maxRadiusToSaveAsStats(_mrtsas)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
NetworkDimensions::NetworkDimensions(std::ifstream & stream) 
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> NetworkDimensions::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	for (unsigned int r = 0 ; r < maxRadiusToSaveAsStats ; ++r)
	{
		res[std::string("NetDimStats_NbNodes_r") + StringifyFixed(r+1)] = 
			(nbNodeTot.size() > r) ? (nbNodeTot[r] - ((r > 0) ? nbNodeTot[r-1] : 0)) : 0;
		res[std::string("NetDimStats_NbOutLinks_r") + StringifyFixed(r+1)] = 
			(nbOutLinksTot.size() > r) ? (nbOutLinksTot[r] - 
			((r > 0) ? nbOutLinksTot[r-1] : 0)) : 0;
		res[std::string("NetDimStats_NbIntraLinks_r") + StringifyFixed(r+1)] = 
			(nbIntraLinksTot.size() > r) ? (nbIntraLinksTot[r] - 
			((r > 0) ? nbIntraLinksTot[r-1] : 0)) : 0;
		res[std::string("NetDimStats_NbPoints_r") + StringifyFixed(r+1)] = 
			(nbPoints.size() > r) ? nbPoints[r] : 0;
	}
	return res;
}

//**********************************************************************
// Compute the all pair distances
//**********************************************************************
bool NetworkDimensions::ComputeMetric(const AbstractNetwork & network)
{
	std::vector<std::vector<unsigned int> > dist, nbNodes, nbIntraLinks, nbOutLinks;
	unsigned int minNbPoint = 3;

	// Maximum radius
	unsigned int rMax;
	if (network.IsSpatial() or fullNetExp)
		rMax = network.size(); //Rely on edge indications if the network is embeded in a space
	else
		rMax = gsl_sf_log(network.size()); // Assume that the network is small world

	unsigned int Nmax = network.size();

	// Fill node list if it's empty
	bool listWasEmpty = false;
	if (nodeList.empty())
	{
		listWasEmpty = true;
		for (unsigned int i = 0 ; i < network.size() ; (i += floor(1.0 / fractCompDim)))
			nodeList.push_back(i);
	}

	unsigned int ind = nodeList[0];
	for (unsigned int i = 0 ; i < nodeList.size() ; ++i, (ind = nodeList[i]))
	{
		dist.push_back(std::vector<unsigned int>());
		nbNodes.push_back(std::vector<unsigned int>());
		nbIntraLinks.push_back(std::vector<unsigned int>());
		nbOutLinks.push_back(std::vector<unsigned int>());

		std::set<unsigned int> visitedNodes, currShell, newShell, tempDiff;

		currShell.insert(ind);
		bool edgeTouched = false;
		for (unsigned int r = 0 ; (fullNetExp or ((r < rMax) and (not edgeTouched))) and 
			(visitedNodes.size() < Nmax) and not currShell.empty() ; ++r)
		{
			dist[i].push_back(r);
			// NbNodes (r)
			visitedNodes.insert(currShell.begin(), currShell.end());
			nbNodes[i].push_back(visitedNodes.size());
			// Stop if connected component is fully visited
			edgeTouched |= (r > 0) ? (nbNodes[i][r] == nbNodes[i][r-1]) : false;

			// IntraLinks (r)
			unsigned int nbIntraTemp = 0;
			for (std::set<unsigned int>::iterator j = currShell.begin() ; j != currShell.end() ; ++j)
				for (std::set<unsigned int>::iterator k = j ; k != currShell.end() ; ++k)
					if (network.AreConnected(*j, *k))
						++nbIntraTemp;
			nbIntraLinks[i].push_back(nbIntraTemp + ((nbIntraLinks[i].size() > 0) ? nbIntraLinks[i].back() : 0));

			// New shell (r + 1)
			newShell.clear();
			for (std::set<unsigned int>::iterator j = currShell.begin() ; j != currShell.end() ; ++j)
				for (unsigned int k = 0 ; k < network.size() ; ++k)
					if (network.AreConnected(*j, k))
					{
						newShell.insert(k);
						edgeTouched |= network.IsNodeOnEdge(k);
					}
			tempDiff.clear();
			std::set_difference(newShell.begin(), newShell.end(), visitedNodes.begin(), visitedNodes.end(), std::inserter(tempDiff, tempDiff.end()));
			newShell = tempDiff;

			// OutLinks (r)
			unsigned int nbOutTemp = 0;
			for (std::set<unsigned int>::iterator j = currShell.begin() ; j != currShell.end() ; ++j)
				for (std::set<unsigned int>::iterator k = newShell.begin() ; k != newShell.end() ; ++k)
					if (network.AreConnected(*j, *k))
						++nbOutTemp;
			nbOutLinks[i].push_back(nbOutTemp + ((nbOutLinks[i].size() > 0) ? nbOutLinks[i].back() : 0));
			
			currShell = newShell;
		}
	}

	std::vector<unsigned int> rmaxes;
	for (unsigned int i = 0 ; i < dist.size() ; ++i)
		rmaxes.push_back(*std::max_element(dist[i].begin(), dist[i].end()));
	unsigned int rMaxEff = *std::max_element(rmaxes.begin(), rmaxes.end());
	for (unsigned int r = 1 ; r < rMaxEff + 1 ; ++r)
	{
		nbNodeTot.push_back(0);
		nbIntraLinksTot.push_back(0);
		nbOutLinksTot.push_back(0);
		nbPoints.push_back(0);
		for (unsigned int i = 0 ; i < nbNodes.size() ; ++i)
			if (nbNodes[i].size() > r)
			{
				nbNodeTot.back() += nbNodes[i][r];
				nbIntraLinksTot.back() += nbIntraLinks[i][r];
				nbOutLinksTot.back() += nbOutLinks[i][r];
				nbPoints.back()++;
			}
		nbNodeTot.back() /= nbPoints.back();
		nbIntraLinksTot.back() /= nbPoints.back();
		nbOutLinksTot.back() /= nbPoints.back();
		// lookout for finite size effects
	}
	if (not nbOutLinksTot.empty())
		nbOutLinksTot.pop_back(); // we don't consider the last r because the L value isn't accurate

	// Get some needed metrics from the network
	DegreeDistrComp *degDistMetr = GetSpecificMetric<Metric, DegreeDistrComp>(network.GetAllMetrics());
	assert(degDistMetr);
	double meanDegree; //, stdDevDegree;
	meanDegree = degDistMetr->GetMeanDegree();

	double *datx;
	double *datyN;
	double *datyL;
	double *datyC;
	double *datyLC;
	for (unsigned int r = 1 ; (rMaxEff >= minNbPoint) and (r <= (rMaxEff - minNbPoint + 1)) ; ++r)
	{
		datx = new double [nbNodeTot.size() - r + 1];
		datyN = new double [nbNodeTot.size() - r + 1];
		datyL = new double [nbOutLinksTot.size() - r + 1];
		datyC = new double [nbIntraLinksTot.size() - r + 1];
		datyLC = new double [nbOutLinksTot.size() - r + 1];
		for (unsigned int r2 = 0 ; r2 < (nbNodeTot.size() - r + 1) ; ++r2)
		{
			datx[r2] = gsl_sf_log(r2 + r);
			datyN[r2] = gsl_sf_log(nbNodeTot[r2 + r - 1]);
			datyC[r2] = (nbIntraLinksTot[r2 + r - 1] > 0) ? gsl_sf_log(nbIntraLinksTot[r2 + r - 1]) : 0;
		}
		for (unsigned int r2 = 0 ; r2 < (nbOutLinksTot.size() - r + 1) ; ++r2)
		{
			datyL[r2] = gsl_sf_log(nbOutLinksTot[r2 + r - 1]);
			datyLC[r2] = gsl_sf_log(nbOutLinksTot[r2 + r - 1] + nbIntraLinksTot[r2 + r - 1]);
		}
		double c0N, c1N, cov00N, cov01N, cov11N, sumsqN,
			c0L, c1L, cov00L, cov01L, cov11L, sumsqL,
			c0C, c1C, cov00C, cov01C, cov11C, sumsqC,
			c0LC, c1LC, cov00LC, cov01LC, cov11LC, sumsqLC;
		gsl_fit_linear(datx, 1, datyN, 1, nbNodeTot.size() - r + 1, &c0N, &c1N, &cov00N, &cov01N, &cov11N, &sumsqN);
		gsl_fit_linear(datx, 1, datyL, 1, nbOutLinksTot.size() - r + 1, &c0L, &c1L, &cov00L, &cov01L, &cov11L, &sumsqL);
		gsl_fit_linear(datx, 1, datyC, 1, nbNodeTot.size() - r + 1, &c0C, &c1C, &cov00C, &cov01C, &cov11C, &sumsqC);
		gsl_fit_linear(datx, 1, datyLC, 1, nbOutLinksTot.size() - r + 1, &c0LC, &c1LC, &cov00LC, &cov01LC, &cov11LC, &sumsqLC);
		delete[] datx;
		delete[] datyN;
		delete[] datyL;
		delete[] datyC;
		delete[] datyLC;
		rDist.push_back(std::make_pair(r, rMaxEff));
		stdDim.push_back(std::make_pair(c1N, cov11N));
		stdPrefact.push_back(std::make_pair(gsl_sf_exp(c0N), cov00N));
		LDim.push_back(std::make_pair(c1L, cov11L));
		LPrefact.push_back(std::make_pair(gsl_sf_exp(c0L), cov00L));
		CDim.push_back(std::make_pair(c1C, cov11C));
		CPrefact.push_back(std::make_pair((c0C != 0) ? gsl_sf_exp(c0C) : 0, cov00C));
		lgDim.push_back(std::make_pair(c1LC, cov11LC));
		lgPrefact.push_back(std::make_pair(gsl_sf_exp(c0LC), cov00LC));
		zetaEst.push_back(std::make_pair(2.0 * CPrefact.back().first / (stdPrefact.back().first * (meanDegree - 2.0)), 0));
	}

	// Computation of final values for dimensions and prefactors
	d = 0; C = 0;
	for (unsigned int i = 0 ; i < stdDim.size() ; ++i)
		if (stdDim[i].first > d)
		{
			d = stdDim[i].first;
			C = stdPrefact[i].first;
		}
	dprime = 0;	E = 0; F = 0; D = 0;
	for (unsigned int i = 0 ; i < lgDim.size() ; ++i)
		if (lgDim[i].first > dprime)
		{
			dprime = lgDim[i].first;
			D = lgPrefact[i].first;
			E = LPrefact[i].first;
			F = CPrefact[i].first;
			zeta = zetaEst[i].first;
		}

	// Computation of parameters for local scaling effects correction
	double *datxN;
	double *datxL;
	double *datxC;
	std::vector<double> tempdatN, tempdatL, tempdatC;
	bool dirN, dirL, dirC;

	dirN = (nbNodeTot.size() > 0) and (nbNodeTot[0] > predN(1));
	for (unsigned int r = 0 ; (r < nbNodeTot.size()) and (dirN ? (nbNodeTot[r] > predN(r+1)) : (nbNodeTot[r] < predN(r+1))) ; ++r)
		tempdatN.push_back(gsl_sf_log(gsl_sf_log(dirN ? (nbNodeTot[r] / predN(r+1)) : (predN(r+1) / nbNodeTot[r]))));
	dirL = (nbOutLinksTot.size() > 0) and (nbOutLinksTot[0] > predL(1));
	for (unsigned int r = 0 ; (r < nbOutLinksTot.size()) and (dirL ? (nbOutLinksTot[r] > predL(r+1)) : (nbOutLinksTot[r] < predL(r+1))) ; ++r)
		tempdatL.push_back(gsl_sf_log(gsl_sf_log(dirL ? (nbOutLinksTot[r] / predL(r+1)) : (predL(r+1) / nbOutLinksTot[r]))));
	dirC = (nbIntraLinksTot.size() > 0) and (nbIntraLinksTot[0] > predC(1));
	for (unsigned int r = 0 ; (r < nbIntraLinksTot.size()) and (dirC ? (nbIntraLinksTot[r] > predC(r+1)) : (nbIntraLinksTot[r] < predC(r+1))) ; ++r)
		tempdatC.push_back(gsl_sf_log(gsl_sf_log(dirC ? (nbIntraLinksTot[r] / predC(r+1)) : (predC(r+1) / nbIntraLinksTot[r]))));

	double c0, c1, cov00, cov01, cov11, sumsq;
	ad = 0; bd = 0; aL = 0; bL = 0; aC = 0; bC = 0;
	if (tempdatN.size() > 1)
	{
		datxN = new double [tempdatN.size()];
		datyN = new double [tempdatN.size()];
		for (unsigned int r = 0 ; r < tempdatN.size() ; ++r)
		{
			datxN[r] = gsl_sf_log(r + 1);
			datyN[r] = tempdatN[r];
		}
		gsl_fit_linear(datxN, 1, datyN, 1, tempdatN.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		ad = gsl_sf_exp(c0);
		bd = - c1;
		delete[] datxN;
		delete[] datyN;
	}
	if (tempdatL.size() > 1)
	{
		datxL = new double [tempdatL.size()];
		datyL = new double [tempdatL.size()];
		for (unsigned int r = 0 ; r < tempdatL.size() ; ++r)
		{
			datxL[r] = gsl_sf_log(r + 1);
			datyL[r] = tempdatL[r];
		}
		gsl_fit_linear(datxL, 1, datyL, 1, tempdatL.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		aL = gsl_sf_exp(c0);
		bL = - c1;
		delete[] datxL;
		delete[] datyL;
	}
	if (tempdatC.size() > 1)
	{
		datxC = new double [tempdatC.size()];
		datyC = new double [tempdatC.size()];
		for (unsigned int r = 0 ; r < tempdatC.size() ; ++r)
		{
			datxC[r] = gsl_sf_log(r + 1);
			datyC[r] = tempdatC[r];
		}
		gsl_fit_linear(datxC, 1, datyC, 1, tempdatC.size(), &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
		aC = gsl_sf_exp(c0);
		bC = - c1;
		delete[] datxC;
		delete[] datyC;
	}

	if (listWasEmpty)
		nodeList.clear();

	return true;
}

//**********************************************************************
// Save the clustering coefficients distribution
//**********************************************************************
bool NetworkDimensions::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string networkDimensionsName("NetworkDimensions");
	std::string netDimFullDataName("NetDimFullData");
	std::string netDimFinalValuesName("NetDimFinalValues");

	if (saver.isSaving(networkDimensionsName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(networkDimensionsName, saver.getCurrFile());

		stream << "rStart\trEnd\tstdDim\tstdDevStdDim\tstdPrefact\tstdDevStdPrefact\t"
			<< "LDim\tstdDevLDim\tLPrefact\tstdDevLPrefact\t" 
			<< "CDim\tstdDevCDim\tCPrefact\tstdDevCPrefact\t" 
			<< "lgDim\tstdDevLgDim\tlgPrefact\tstdDevLgPrefact\t" 
			<< "ZetaEst\tstdDevZetaEst\t"
			<< std::endl;
		for (unsigned int i = 0 ; i < rDist.size() ; ++i)
		{
			stream 
				<< rDist[i].first      << "\t" << rDist[i].second      << "\t"
				<< stdDim[i].first     << "\t" << stdDim[i].second     << "\t"
				<< stdPrefact[i].first << "\t" << stdPrefact[i].second << "\t"
				<< LDim[i].first       << "\t" << LDim[i].second       << "\t"
				<< LPrefact[i].first   << "\t" << LPrefact[i].second   << "\t"
				<< CDim[i].first       << "\t" << CDim[i].second       << "\t"
				<< CPrefact[i].first   << "\t" << CPrefact[i].second   << "\t"
				<< lgDim[i].first      << "\t" << lgDim[i].second      << "\t"
				<< lgPrefact[i].first  << "\t" << lgPrefact[i].second  << "\t"
				<< zetaEst[i].first    << "\t" << zetaEst[i].second    << "\t"
				<< std::endl;
		}

		allSaved &= stream.good();
	}
	if (saver.isSaving(netDimFullDataName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(netDimFullDataName, saver.getCurrFile());

		stream << "r\tNbNodes\tNbIntra\tNbOut\tShellNodes\tShellIntra\tShellOut\tNbPoints\tNbNodesPrev\tNbIntraPrev\tNbOutPrev\t" << std::endl;
		for (unsigned int i  = 0 ; i < nbPoints.size() ; ++i)
		{
			stream 
				<< i + 1 << "\t"
				<< ((i < nbNodeTot.size()) ? nbNodeTot[i] : -1) << "\t"
				<< ((i < nbIntraLinksTot.size()) ? nbIntraLinksTot[i] : -1) << "\t"
				<< ((i < nbOutLinksTot.size()) ? nbOutLinksTot[i] : -1) << "\t";
			if (i > 0)
				stream
					<< ((i < nbNodeTot.size()) ? nbNodeTot[i] - nbNodeTot[i-1] : -1) << "\t"
					<< ((i < nbIntraLinksTot.size()) ? nbIntraLinksTot[i] - nbIntraLinksTot[i-1] : -1) << "\t"
					<< ((i < nbOutLinksTot.size()) ? nbOutLinksTot[i] - nbOutLinksTot[i-1]: -1) << "\t";
			else
				stream 
					<< ((i < nbNodeTot.size()) ? nbNodeTot[i] : -1) << "\t"
					<< ((i < nbIntraLinksTot.size()) ? nbIntraLinksTot[i] : -1) << "\t"
					<< ((i < nbOutLinksTot.size()) ? nbOutLinksTot[i] : -1) << "\t";

			stream
				<< ((i < nbPoints.size()) ? nbPoints[i] : -1) << "\t"
				<< predNRect(i+1) << "\t"
				<< predCRect(i+1) << "\t"
				<< predLRect(i+1) << std::endl;
		}

		allSaved &= stream.good();
	}
	if (saver.isSaving(netDimFinalValuesName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(netDimFinalValuesName, saver.getCurrFile());

		stream << "d\tC\tdprime\tE\tF\tD\tad\tbd\taL\tbL\taC\tbC\tzetaEst\t" << std::endl;
		stream
			<< d << "\t" << C << "\t" 
			<< dprime << "\t" << E << "\t" << F << "\t" << D << "\t"
			<< ad << "\t" << bd << "\t"
			<< aL << "\t" << bL << "\t"
			<< aC << "\t" << bC << "\t" 
			<< zeta << std::endl;

		allSaved &= stream.good();
	}
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void NetworkDimensions::Initialize()
{
	Metric::Initialize();
	rDist.clear();
	stdDim.clear();
	stdPrefact.clear();
	LDim.clear();
	LPrefact.clear();
	CDim.clear();
	CPrefact.clear();
	lgDim.clear();
	lgPrefact.clear();
	zetaEst.clear();
	nbNodeTot.clear();
	nbOutLinksTot.clear();
	nbIntraLinksTot.clear();
	nbPoints.clear();
	d = 0; C = 0;
	dprime = 0;	E = 0; F = 0; D = 0;
	zeta = 0;
	ad = 0; bd = 0; aL = 0; bL = 0; aC = 0; bC = 0;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * NetworkDimensions::BuildCopy() const
{
	return new NetworkDimensions(fractCompDim, nodeList, fullNetExp, 
		maxRadiusToSaveAsStats); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool NetworkDimensions::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();
	stream >> fractCompDim;
	unsigned int tmpSize;
	stream >> tmpSize;
	nodeList = std::vector<int>(tmpSize, 0);
	for (unsigned int i = 0 ; i < tmpSize ; ++i)
		stream >> nodeList[i];
	stream >> fullNetExp;
	stream >> maxRadiusToSaveAsStats;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool NetworkDimensions::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< fractCompDim << std::endl
		<< nodeList.size() << std::endl;
	for (unsigned int i = 0 ; i < nodeList.size() ; ++i)
		stream << nodeList[i] << std::endl;
	stream 
		<< fullNetExp << std::endl
		<< maxRadiusToSaveAsStats << std::endl;
	return stream.good();
}

//**********************************************************************
// Quick estimations for different values
//**********************************************************************
double NetworkDimensions::predN(double r) const { return (C * pow(r, d)); }
double NetworkDimensions::predL(double r) const { return (E * pow(r, dprime)); }
double NetworkDimensions::predC(double r) const { return (F * pow(r, dprime)); }
double NetworkDimensions::predNRect(double r) const { return predN(r) * gsl_sf_exp(ad * pow(r, -bd)); }
double NetworkDimensions::predLRect(double r) const { return predL(r) * gsl_sf_exp(aL * pow(r, -bL)); }
double NetworkDimensions::predCRect(double r) const { return predC(r) * gsl_sf_exp(aC * pow(r, -bC)); }

