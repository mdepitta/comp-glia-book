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

#include "NetworkConstructStrat.h"

#include "Network.h"
#include "SpatialNetwork.h"
#include "ChIModel.h"
#include <hull.h>

#include <algorithm>

using namespace AstroModel;
using namespace std;

string ScaleFreeStrat::ClassName(
	"ScaleFreeBuildingStrategy");
string SpatialScaleFreeStrat::ClassName(
	"SpatialScaleFreeBuildingStrategy");
string SpatialRegularStrat::ClassName(
	"SpatialRegularBuildingStrategy");
string SpatialConnectionRadiusStrat::ClassName(
	"SpatialConnectionRadiusBuildingStrategy");
string VoronoiConstructStrat::ClassName(
	"VoronoiDiagramBuildingStrategy");
string ErdosRenyiRandomStrat::ClassName(
	"ErdosRenyiRandomBuildingStrategy");
string SmallWorldStrat::ClassName(
	"SmallWorldBuildingStrategy");
string ThresholdDeterminationStrat::ClassName(
	"ThresholdDeterminationStrategy");
string FunctionalTopoStrat::ClassName(
	"FunctionalTopoStrategy");
string ConstrFromChISim::ClassName(
	"ConstructFromChISimulationStrategy");
string FromFileConstrStrat::ClassName(
	"LoadFromFileBuildingStrategy");
string FileListConstrStrat::ClassName(
	"LoadFromFileListBuildingStrategy");
string RandomIsolatedCells::ClassName(
	"RandomIsolatedCellsBuildingStrategy");
string StarLikeAstroNetConstrStrat::ClassName(
	"StarLikeAstroNetBuildingStrategy");
string ShellStructScramblerStrat::ClassName(
	"ShellStructureScramblerBuildingStrategy");

//********************************************************************//
//*********** N E T W O R K   C O N S T R U C T   S T R A T **********//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
NetworkConstructStrat::NetworkConstructStrat(ParamHandler & h)
{
	normLkStr = h.getParam<bool>("-normLink");
	makeDir = h.getParam<bool>("-makeDirNet");
}

//**********************************************************************
// Build the network
//**********************************************************************
bool NetworkConstructStrat::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & fact, ResultSaver) const
{
	if (normLkStr)
	{
		NetworkEdge *temp = fact.Create();
		if (dynamic_cast<CouplingFunction*>(temp))
		{
			CouplingFunction::normalizeLinkStrength(network);
		}
		delete temp;
	}
	return true;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool NetworkConstructStrat::LoadFromStream(std::ifstream & stream)
{
	stream >> normLkStr;
	stream >> makeDir;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool NetworkConstructStrat::SaveToStream(std::ofstream & stream) const
{
	stream << normLkStr << std::endl;
	stream << makeDir << std::endl;
	return stream.good();
}

//**********************************************************************
// Returns a copy of the network construct strat
//**********************************************************************
NetworkConstructStrat* NetworkConstructStrat::BuildCopy() const
{
	std::string tmpFileName = "TmpBuildCopyConstrStrat";
	std::string clName = this->GetClassName();

	std::ofstream stream1(tmpFileName.c_str());
	this->SaveToStream(stream1);
	stream1.close();

	std::ifstream stream2(tmpFileName.c_str());
	assert(AbstractFactory<NetworkConstructStrat>::Factories[clName]);
	return AbstractFactory<NetworkConstructStrat>::
		Factories[clName]->CreateFromStream(stream2);
}

//**********************************************************************
// 
//**********************************************************************
std::vector<Metric *> NetworkConstructStrat::GetAllMetrics() const
{
	return std::vector<Metric *>();
}

//********************************************************************//
//** S P A T I A L   N E T W O R K   C O N S T R U C T   S T R A T ***//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
SpatialNetConstrStrat::SpatialNetConstrStrat(ParamHandler & h) : 
	NetworkConstructStrat::NetworkConstructStrat(h) 
{

}

//**********************************************************************
// Special constructor
//**********************************************************************
SpatialNetConstrStrat::SpatialNetConstrStrat(bool _n) :
	NetworkConstructStrat::NetworkConstructStrat(_n) 
{

}

//**********************************************************************
// Build the network
//**********************************************************************
bool SpatialNetConstrStrat::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
	AbstractSpatialNetwork *sptlNetwrk = 
		dynamic_cast<AbstractSpatialNetwork*>(&network);
	if (sptlNetwrk)
	{
		bool ok = true;
		ok &= this->buildSpatialNetwork(*sptlNetwrk, factory);
		return ok and NetworkConstructStrat::BuildNetwork(network, factory, saver);
	}
	else
		return false;
}

//********************************************************************//
//****************** S C A L E   F R E E   S T R A T *****************//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
ScaleFreeStrat::ScaleFreeStrat(ParamHandler & h) :
	NetworkConstructStrat::NetworkConstructStrat(h)
{
	h.getParam<double>("-scaleFree", 0);
}

//**********************************************************************
// Special constructor
//**********************************************************************
ScaleFreeStrat::ScaleFreeStrat(double _g, bool _n) : 
NetworkConstructStrat::NetworkConstructStrat(_n), gamma(_g) 
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ScaleFreeStrat::ScaleFreeStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool ScaleFreeStrat::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
	// Deprecated, use SpatialScaleFree with large rc instead
	network.SetAbstractEdge(0, 1, factory.Create());
	network.SetAbstractEdge(1, 0, network.GetAbstractEdge(0, 1));
	double sumDeg = 1;
	for (unsigned int i = 1 ; i < network.size() ; ++i)
	{
		for (unsigned int j = 0 ; j < i ; ++j)
		{
			if (UnifRand() < network.GetNodeDegree(j) / sumDeg)
			{
				network.SetAbstractEdge(i, j, factory.Create());
				network.SetAbstractEdge(j, i, network.GetAbstractEdge(i, j));
				sumDeg += 2;
			}
		}
	}
	network.SetDirected(false);
	return NetworkConstructStrat::BuildNetwork(network, factory, saver);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler ScaleFreeStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "ScaleFreeGamma", gamma;
	return params;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool ScaleFreeStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> gamma;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool ScaleFreeStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << gamma << std::endl;
	return stream.good();
}

//********************************************************************//
//********** S P A T I A L   S C A L E   F R E E   S T R A T *********//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
SpatialScaleFreeStrat::SpatialScaleFreeStrat(ParamHandler & h) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(h)
{
	rc =  h.getParam<double>("-spatialScaleFree", 0);
	newLinks =  h.getParam<unsigned int>("-spatialScaleFree", 1);
}

//**********************************************************************
// Special constructor
//**********************************************************************
SpatialScaleFreeStrat::SpatialScaleFreeStrat(double _rc, 
	unsigned int _nl, bool _n) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(_n), rc(_rc), 
	newLinks(_nl)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
SpatialScaleFreeStrat::SpatialScaleFreeStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler SpatialScaleFreeStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "SpatialScaleFreeInteractionRange", rc;
	params <= "SpatialScaleFreeNewLinks", newLinks;
	return params;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool SpatialScaleFreeStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> rc;
	stream >> newLinks;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool SpatialScaleFreeStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream 
		<< rc       << std::endl
		<< newLinks << std::endl;
	return stream.good();
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool SpatialScaleFreeStrat::buildSpatialNetwork(AbstractSpatialNetwork &
	network, AbstractFactory<NetworkEdge> & factory) const
{
	std::vector<unsigned int> addedNodes;

	unsigned int tempInd = floor((double) network.size() * UnifRand());
	std::vector<double> tempProbas(network.size(), 0);
	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		if (i != tempInd)
			tempProbas[i] = LinkProbaRaw(network, i, tempInd, true);
	}
	NormalizeProbas(tempProbas);

	std::vector<unsigned int> checkOrder(network.size(), 0);
	for (unsigned int i = 0 ; i < checkOrder.size() ; ++i)
		checkOrder[i] = i;
	RandomSwapsOnVector(checkOrder);

	bool linked = false;
	unsigned int ind = 0;
	while (not linked)
	{
		if ((linked = TrueWithProba(tempProbas[checkOrder[ind]])))
		{
			addedNodes.push_back(tempInd);
			addedNodes.push_back(checkOrder[ind]);
			network.SetAbstractEdge(tempInd, checkOrder[ind], 
				factory.Create());
			network.SetAbstractEdge(checkOrder[ind], tempInd, 
				network.GetAbstractEdge(tempInd, checkOrder[ind]));
		}
		++ind;
		ind = ind % checkOrder.size();
	}

	// Choosing the nodes addition order
	std::vector<unsigned int> addOrder;
	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		bool alreadyAdded = false;
		for (unsigned int j = 0 ; j < addedNodes.size() ; ++j)
			alreadyAdded |= (i == addedNodes[j]);
		if (not alreadyAdded)
			addOrder.push_back(i);
	}
	RandomSwapsOnVector(addOrder);

	// For each node, preferentially attach it to existent nodes
	for (unsigned int i = 0 ; i < addOrder.size() ; ++i)
	{
		std::vector<double> probas(addedNodes.size(), 0);
		for (unsigned int j = 0 ; j < addedNodes.size() ; ++j)
			probas[j] = LinkProbaRaw(network, addOrder[i], addedNodes[j]);
		NormalizeProbas(probas);

		unsigned int currNbLinks = 0;
		std::vector<double> tempProbas = probas;
		bool allNull = false;
		for (unsigned int j = 0 ; (currNbLinks < newLinks) and 
			(not allNull) ; ++j, j = j % probas.size())
		{
			if (TrueWithProba(tempProbas[j]))
			{
				probas[j] = 0;
				NormalizeProbas(probas);
				tempProbas = probas;
				
				network.SetAbstractEdge(addOrder[i], addedNodes[j], 
					factory.Create());
				network.SetAbstractEdge(addedNodes[j], addOrder[i],
					network.GetAbstractEdge(addOrder[i], addedNodes[j]));
				++currNbLinks;

				unsigned int nbNull = 0;
				for (unsigned int k = 0 ; k < probas.size() ; ++k)
					nbNull = nbNull + (probas[k] ? 0 : 1);
				allNull = (nbNull == probas.size());
			}
			else
			{
				tempProbas[j] = 0;
				NormalizeProbas(tempProbas);
			}
		}
		addedNodes.push_back(addOrder[i]);
	}
	return true;
}

//**********************************************************************
// Compute raw link probabilities
//**********************************************************************
double SpatialScaleFreeStrat::LinkProbaRaw(AbstractSpatialNetwork & network, 
	unsigned int add, unsigned int exist, bool considerUnconnected) const
{
	return ((double) network.GetNodeDegree(exist) + 
		(considerUnconnected ? 1.0 : 0.0)) *
		exp(-1.0 * network.GetDistanceBetween(add, exist) / rc);
}

//**********************************************************************
// Normalize probabilities
//**********************************************************************
void SpatialScaleFreeStrat::NormalizeProbas(std::vector<double> & probas) const
{
	double totVal = 0;
	for (unsigned int i = 0 ; i < probas.size() ; ++i)
		totVal += probas[i];
	if (totVal != 0)
	{
		for (unsigned int i = 0 ; i < probas.size() ; ++i)
			probas[i] /= totVal;
	}
	else
	{
		for (unsigned int i = 0 ; i < probas.size() ; ++i)
			probas[i] /= 1.0 / ((double)(probas.size()));
	}
}

//********************************************************************//
//************* S P A T I A L   R E G U L A R   S T R A T ************//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
SpatialRegularStrat::SpatialRegularStrat(ParamHandler & h) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(h)
{
	degree = h.getParam<unsigned int>("-regularConstr", 0);
	maxLinkDist = h.getParam<double>("-regularConstr", 1);
}

//**********************************************************************
// Special constructor
//**********************************************************************
SpatialRegularStrat::SpatialRegularStrat(unsigned int _d, double _mld, 
	bool _n) : 	SpatialNetConstrStrat::SpatialNetConstrStrat(_n), 
	degree(_d),	maxLinkDist(_mld) 
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
SpatialRegularStrat::SpatialRegularStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler SpatialRegularStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "SpatialRegularDegree", degree;
	params <= "SpatialRegularMaxLinkDist", maxLinkDist;
	return params;
}

bool IndComparator ( const std::pair<double, unsigned int> & l, const std::pair<double, unsigned int> & r)
{ return l.first < r.first; }

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool SpatialRegularStrat::buildSpatialNetwork(AbstractSpatialNetwork & network,
	AbstractFactory<NetworkEdge> & factory) const
{
	// Compute all spatial distances first
	std::vector<std::vector<unsigned int> > sortedSpatialDistances(network.size(), std::vector<unsigned int>());
	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		std::vector<std::pair<double, unsigned int> > tempDistances(network.size(), std::make_pair(DEFAULT_MAX_PATH, 0));
		for (unsigned int j = 0 ; j < network.size() ; ++j)
		{
			tempDistances[j].first = network.GetDistanceBetween(i, j);
			tempDistances[j].second = j;
		}
		std::sort(tempDistances.begin(), tempDistances.end(), IndComparator);
		for (unsigned int k =  0 ; k < tempDistances.size() ; ++k)
			sortedSpatialDistances[i].push_back(tempDistances[k].second);
	}

	for (unsigned int k = 1 ; k <= degree ; ++k)
	{
		std::vector<unsigned int> checkOrder(network.size(), 0);
		for (unsigned int i = 0 ; i < checkOrder.size() ; ++i)
			checkOrder[i] = i;
		for (unsigned int i = 0 ; i < checkOrder.size() ; ++i)
		{
			unsigned int tempInd = floor(UnifRand() * 
				(double)checkOrder.size());
			unsigned int lastVal = checkOrder[i];
			checkOrder[i] = checkOrder[tempInd];
			checkOrder[tempInd] = lastVal;
		}

		unsigned int i = checkOrder[0];
		for (unsigned int o = 0 ; o < checkOrder.size() ; ++o)
		{
			i = checkOrder[o];
			if (network.GetNodeDegree(i) < k)
			{
				unsigned int minCell = sortedSpatialDistances[i][0];
				unsigned int m = 0;
				for (m = 0 ; (m < sortedSpatialDistances[i].size()) and 
					not ((i != minCell) and (not network.GetAbstractEdge(i, minCell)) and (network.GetNodeDegree(minCell) < k)) ;
						minCell = sortedSpatialDistances[i][m++]);
				if ((minCell < network.size()) and (network.GetDistanceBetween(i, minCell) <= maxLinkDist))
				{
					network.SetAbstractEdge(i, minCell, factory.Create());
					network.SetAbstractEdge(minCell, i, network.GetAbstractEdge(i, minCell));
				}
			}
		}
	}
	return true;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool SpatialRegularStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> degree;
	stream >> maxLinkDist;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool SpatialRegularStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << degree << std::endl;
	stream << maxLinkDist << std::endl;
	return stream.good();
}

//********************************************************************//
//*** S P A T I A L   C O N N E C T I O N   R A D I U S   S T R A T **//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
SpatialConnectionRadiusStrat::SpatialConnectionRadiusStrat(
	ParamHandler & h) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(h)
{
	radius = h.getParam<double>("-linkRadiusConstr", 0);
}

//**********************************************************************
// Special constructor
//**********************************************************************
SpatialConnectionRadiusStrat::SpatialConnectionRadiusStrat(double _r, 
	bool _n) : SpatialNetConstrStrat::SpatialNetConstrStrat(_n), 
	radius(_r) 
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
SpatialConnectionRadiusStrat::SpatialConnectionRadiusStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler SpatialConnectionRadiusStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "SpatialConnectionRadius", radius;
	return params;
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool SpatialConnectionRadiusStrat::buildSpatialNetwork(AbstractSpatialNetwork & network,
	AbstractFactory<NetworkEdge> & factory) const
{
	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		for (unsigned int j = i + 1 ; j < network.size(i) ; ++j)
		{
			if (network.GetDistanceBetween(i, j) < radius)
			{
				network.SetAbstractEdge(i, j, factory.Create());
				network.SetAbstractEdge(j, i, network.GetAbstractEdge(i, j));
			}
		}
	}
	return true;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool SpatialConnectionRadiusStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> radius;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool SpatialConnectionRadiusStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << radius << std::endl;
	return stream.good();
}

//********************************************************************//
//*** V O R O N O I   D I A G R A M   C O N S T R U C T   S T R A T **//
//********************************************************************//

//**********************************************************************
// C-like functions and data to be used by libhull
//**********************************************************************
FILE *INFILE, *OUTFILE, *TFILE;
FILE *DFILE = stderr;

static std::vector<site> voroSites;
static unsigned int voroCurrInd;
static std::vector<std::vector<bool> > tmpAdj;

site get_next_site(void) 
{
	return (voroCurrInd < voroSites.size()) ?
		voroSites[voroCurrInd++] : NULL;
}

long site_numm(site p) 
{
	if (p==hull_infinity) 
	{
		return -1;
	}
	if (!p) 
	{
		return -2;
	}
	for (unsigned int i = 0 ; i < voroSites.size() ; ++i)
		if (voroSites[i] == p)
		{
			return i;
		}
	return -3;
}

void* visit_test(simplex *sim, void *)
{
	if (tmpAdj.empty())
		tmpAdj = std::vector<std::vector<bool> >(voroSites.size(), 
			std::vector<bool>(voroSites.size(), false));

	std::vector<unsigned int> tmpNums;

	if (sim->neigh and not sim->peak.vert)
	{
		// Check that neighb don't contain infinite vertex
		bool containsInfi = false;
		for (int i = 0 ; (not containsInfi) and (i < cdim) ; ++i)
			containsInfi |= (site_numm(sim->neigh[i].vert) < 0);

		if (not containsInfi)
		{
			long ind1, ind2;
			for (int i = 0 ; i < cdim ; ++i)
			{
				ind1 = site_numm(sim->neigh[i].vert);
				if (ind1 >= 0)
					for (int j = i + 1 ; j < cdim ; ++j)
					{
						ind2 = site_numm(sim->neigh[j].vert);
						if (ind2 >= 0)
						{
							tmpAdj[ind1][ind2] = true;
							tmpAdj[ind2][ind1] = true;
						}
					}
			}
		}
	}
	return 0;
}

//**********************************************************************
// Default constructor
//**********************************************************************
VoronoiConstructStrat::VoronoiConstructStrat(
	ParamHandler & h) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(h)
{
	distMultFact = h.getParam<double>("-voronoiParam", 0);
	useGhostPoints = h.getParam<bool>("-voronoiParamGhostPoints", 0);
	maxLinkDist = h.getParam<double>("-voronoiParam", 1);
}

//**********************************************************************
// Special constructor
//**********************************************************************
VoronoiConstructStrat::VoronoiConstructStrat(double _dmF, bool _ugp, double _mld) :
	distMultFact(_dmF), useGhostPoints(_ugp), maxLinkDist(_mld)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
VoronoiConstructStrat::VoronoiConstructStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool VoronoiConstructStrat::buildSpatialNetwork(
	AbstractSpatialNetwork & network,
	AbstractFactory<NetworkEdge> & factory) const
{
	if (voroSites.empty())
	{
		voroSites.resize(network.size() + (useGhostPoints ? (2*network.GetDim()) : 0));
		for (unsigned int i = 0 ; i < voroSites.size() ; ++i)
			voroSites[i] = (new double[network.GetDim()]);
	}

	voroCurrInd = 0;

	std::vector<std::pair<double, double> > minMaxD = 
		std::vector<std::pair<double, double> >(network.GetDim(), 
		std::make_pair(DEFAULT_MAX_VAL, -DEFAULT_MAX_VAL));
	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		for (unsigned int d = 0 ; d < network.GetDim() ; ++d)
		{
			voroSites[voroCurrInd][d] = (*network.Pos(i))[d] * distMultFact;
			minMaxD[d].first = std::min(minMaxD[d].first, voroSites[voroCurrInd][d]);
			minMaxD[d].second = std::max(minMaxD[d].second, voroSites[voroCurrInd][d]);
		}
		++voroCurrInd;
	}
	// Add ghost points to avoid border effects
	if (useGhostPoints)
		for (unsigned int d = 0 ; d < network.GetDim() ; ++d)
		{
			voroSites[voroCurrInd][d] = minMaxD[d].first - 0.5 * 
				(minMaxD[d].second - minMaxD[d].first);
			for (unsigned int d2 = 0 ; d2 < network.GetDim() ; ++d2)
				if (d2 != d)
					voroSites[voroCurrInd][d2] = 0.5 * (minMaxD[d2].second + minMaxD[d2].first);
			++voroCurrInd;

			voroSites[voroCurrInd][d] = minMaxD[d].second + 0.5 * 
				(minMaxD[d].second - minMaxD[d].first);
			for (unsigned int d2 = 0 ; d2 < network.GetDim() ; ++d2)
				if (d2 != d)
					voroSites[voroCurrInd][d2] = 0.5 * (minMaxD[d2].second + minMaxD[d2].first);
			++voroCurrInd;
		}

	voroCurrInd = 0;

	simplex *root;
	root = build_convex_hull(get_next_site, site_numm, network.GetDim(), 1);

	tmpAdj.clear();
	visit_triang(root, visit_test);
	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		for (unsigned int j = 0 ; j < network.size() ; ++j)
		{
			if (tmpAdj[i][j] and not network.AreConnected(i, j) and
				(network.GetDistanceBetween(i, j) <= maxLinkDist))
			{
				network.SetAbstractEdge(i, j, factory.Create());
				network.SetAbstractEdge(j, i, 
					network.GetAbstractEdge(i, j));
			}

		}
	}
	tmpAdj.clear();

	free_hull_storage();

	for (unsigned int i = 0 ; i < voroSites.size() ; ++i)
		if (voroSites[i])
			delete[] voroSites[i];
	voroSites.clear();

	return true;
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler VoronoiConstructStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "VoronoiConstrDistMultFact", distMultFact;
	params <= "VoronoiConstrUseGhostPoints", useGhostPoints;
	params <= "VoronoiConstrMaxLinkDist", maxLinkDist;
	return params;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool VoronoiConstructStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> distMultFact;
	stream >> useGhostPoints;
	stream >> maxLinkDist;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool VoronoiConstructStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << distMultFact << std::endl
		<< useGhostPoints << std::endl
		<< maxLinkDist << std::endl;
	return stream.good();
}

//********************************************************************//
//********** E R D O S   R E N Y I   R A N D O M   S T R A T *********//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
ErdosRenyiRandomStrat::ErdosRenyiRandomStrat(ParamHandler & h) :
	NetworkConstructStrat::NetworkConstructStrat(h)
{
	meanDegree = h.getParam<double>("-erdosRenyiMeanDeg", 0);
}

//**********************************************************************
// Special constructor
//**********************************************************************
ErdosRenyiRandomStrat::ErdosRenyiRandomStrat(double _k, bool _n) :
	NetworkConstructStrat::NetworkConstructStrat(_n), meanDegree(_k)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ErdosRenyiRandomStrat::ErdosRenyiRandomStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler ErdosRenyiRandomStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "ErdosRenyiMeanDegree", meanDegree;
	return params;
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool ErdosRenyiRandomStrat::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
	assert(network.size() > 1);
	double p = meanDegree / ((double) network.size() - 1.0);

	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		for (unsigned int j = (makeDir ? 0 : (i + 1)) ; j < network.size(i) ; ++j)
		{
			if((i != j) and TrueWithProba(p))
			{
				network.SetAbstractEdge(i, j, factory.Create());
				if (not makeDir)
					network.SetAbstractEdge(j, i, network.GetAbstractEdge(i, j));
			}
		}
	}
	network.SetDirected(makeDir);
	return NetworkConstructStrat::BuildNetwork(network, factory, saver);
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool ErdosRenyiRandomStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> meanDegree;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool ErdosRenyiRandomStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << meanDegree << std::endl;
	return stream.good();
}

//********************************************************************//
//********* S P A T I A L   S M A L L   W O R L D   S T R A T ********//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
SmallWorldStrat::SmallWorldStrat(ParamHandler & h) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(h)
{
	prob = h.getParam<double>("-SWParams", 0);
	neighbDist = h.getParam<int>("-SWParams", 1);
}

//**********************************************************************
// Special constructor
//**********************************************************************
SmallWorldStrat::SmallWorldStrat(double _p, int _nD, bool _n) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(_n), prob(_p), 
	neighbDist(_nD) 
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
SmallWorldStrat::SmallWorldStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler SmallWorldStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "SmallWorldProba", prob;
	params <= "SmallWorldNeighbDist", neighbDist;
	return params;
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool SmallWorldStrat::buildSpatialNetwork(AbstractSpatialNetwork & network,
	AbstractFactory<NetworkEdge> & factory) const
{
	unsigned int dim = network.GetDim();
	unsigned int nbSide = ceil(pow(network.size(), 1.0 / (double)dim)); 
	bool shortcutCreated = false;
	bool isToroidal = network.IsToroidal();
	std::vector<unsigned int> coords;
	std::vector<unsigned int> toLink;
	unsigned int tempInd, randInd;

	// Check that the network is (hyper)cubic
	assert(pow(nbSide, dim) == network.size());

	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		coords = std::vector<unsigned int>(dim, 0);
		toLink.clear();
		
		for (unsigned int d = 0 ; d < dim ; ++d)
			coords[d] = ((unsigned int)floor((double)i / pow((double)nbSide, d))) % nbSide;

		// List links to be made
		for (unsigned int d = 0 ; d < dim ; ++d)
			for (unsigned int m = 1 ; m <= neighbDist ; ++m)
				if (isToroidal or (coords[d] + m < nbSide))
				{
					tempInd = 0;
					for (unsigned int dp = 0 ; dp < dim ; ++dp)
						tempInd += 
							((coords[dp] + ((dp == d) ? m : 0)) % nbSide) *
							pow((double)nbSide, dp);
					toLink.push_back(tempInd);
				}
		
		// Make links and rewire with proba prob
		for (unsigned int j = 0 ; j < toLink.size() ; ++j)
		{
			if (not network.GetAbstractEdge(i, toLink[j]))
			{
				if(TrueWithProba(prob))
				{
					for (randInd = toLink[j] ; 
						randInd == toLink[j] ; randInd = 
						floor(UnifRand() * (double)network.size()));
					toLink[j] = randInd;
					shortcutCreated |= true;
				}
				network.SetAbstractEdge(i, toLink[j], factory.Create());
				network.SetAbstractEdge(toLink[j], i, 
					network.GetAbstractEdge(i, toLink[j]));
			}
		}
		network.SetIsStrictlySpatial(not shortcutCreated);
	}

	return true;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool SmallWorldStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> prob;
	stream >> neighbDist;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool SmallWorldStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << 
		prob << std::endl << 
		neighbDist << std::endl;
	return stream.good();
}

//********************************************************************//
//***** T H R E S H O L D   D E T E R M I N A T I O N   S T R A T ****//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
ThresholdDeterminationStrat::ThresholdDeterminationStrat(ParamHandler & h)
	: NetworkConstructStrat::NetworkConstructStrat(h)
{
	degree = h.getParam<unsigned int >("-ThreshDetermNet", 0);
	nbDeriv = h.getParam<unsigned int >("-ThreshDetermNet", 1);
}

//**********************************************************************
// Special constructor
//**********************************************************************
ThresholdDeterminationStrat::ThresholdDeterminationStrat(unsigned int _k,
	unsigned int _d, bool _n) : 
	NetworkConstructStrat::NetworkConstructStrat(_n), degree(_k), 
	nbDeriv(_d) 
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ThresholdDeterminationStrat::ThresholdDeterminationStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler ThresholdDeterminationStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "ThresholdDeterminationDegree", degree;
	params <= "ThresholdDeterminationDerivs", nbDeriv;
	return params;
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool ThresholdDeterminationStrat::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
	assert(network.size() > degree * (nbDeriv + 1));
	for (unsigned int i = 1 ; i <= degree ; ++i)
	{
		network.SetAbstractEdge(0, degree + i, factory.Create());
		network.SetAbstractEdge(degree + i, 0, network.GetAbstractEdge(0, degree + i));

		network.SetAbstractEdge(2*degree + i, degree + i, factory.Create());
		network.SetAbstractEdge(degree + i, 2*degree + i, network.GetAbstractEdge(2*degree + i, degree + i));

		network.SetAbstractEdge(2*degree + i, i, factory.Create());
		network.SetAbstractEdge(i, 2*degree + i, network.GetAbstractEdge(2*degree + i, i));

		for (unsigned int j = 1 ; j <= nbDeriv ; ++j)
		{
			network.SetAbstractEdge(degree + i, (j+2) * degree + i, factory.Create());
			network.SetAbstractEdge((j+2) * degree + i, degree + i, network.GetAbstractEdge(degree + i, (j+2) * degree + i));

			network.SetAbstractEdge((j+2) * degree + i, (j+nbDeriv+2) * degree + i, factory.Create());
			network.SetAbstractEdge((j+nbDeriv+2) * degree + i, (j+2) * degree + i, network.GetAbstractEdge((j+2) * degree + i, (j+nbDeriv+2) * degree + i));
		}
	}
	network.SetDirected(false);
	return NetworkConstructStrat::BuildNetwork(network, factory, saver);
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool ThresholdDeterminationStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> degree;
	stream >> nbDeriv;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool ThresholdDeterminationStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << degree << std::endl;
	stream << nbDeriv << std::endl;
	return stream.good();
}

//**********************************************************************
// Returns the network central node degree
//**********************************************************************
unsigned int ThresholdDeterminationStrat::GetDegree() const
{
	return degree;
}

//**********************************************************************
// Returns the number of derivations (sinks) for stim cells
//**********************************************************************
unsigned int ThresholdDeterminationStrat::GetNbDerivs() const
{
	return nbDeriv;
}

//********************************************************************//
//***** F U N C T I O N A L   T O P O   C O N S T R   S T R A T ******//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
FunctionalTopoStrat::FunctionalTopoStrat(ParamHandler & h) :
	NetworkConstructStrat::NetworkConstructStrat(h), useAbsoluteThresh(false),
	absoluteThresh(0)
{
	dirLinks = h.getParam<bool>("-functTopoDirectedLinks", 0);
	useFixedThresh = h.getParam<bool>("-functTopoUseThresh", 0);
	threshold = h.getParam<double>("-functTopoUseThresh", 1);
	stdDevCoeff = h.getParam<double>("-functTopoAutoThresh", 0);
TRACE(stdDevCoeff)
}

//**********************************************************************
// Special constructor
//**********************************************************************
FunctionalTopoStrat::FunctionalTopoStrat(
	const std::vector<std::vector<double> > & _vals, bool _dL, bool _uT, 
			double _t, double _sdC, bool _n, bool _uat, double _at) : 
	NetworkConstructStrat::NetworkConstructStrat(_n), dirLinks(_dL), 
	useFixedThresh(_uT), useAbsoluteThresh(_uat), threshold(_t), 
	stdDevCoeff(_sdC), absoluteThresh(_at)
{
	SetValues(_vals);
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
FunctionalTopoStrat::FunctionalTopoStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool FunctionalTopoStrat::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
	assert(threshold <= 1.0 and threshold >= 0.0);
	assert(network.size() == values.size());
	
	std::vector<double> sortVals = GetSortedValuesFromMat(values);
	double threshVal = 0;
	
	if (useAbsoluteThresh)
		threshVal = absoluteThresh;
	else if (useFixedThresh)
		threshVal = sortVals[floor((1.0 - threshold) * sortVals.size())];
	else
	{
		// automatic threshval detection
		double nbBinsTmp = 40.0;
		double binSize = (sortVals.back() - sortVals[0]) / nbBinsTmp;
		std::vector<double> bins(nbBinsTmp, 0);
		std::vector<double> vals(nbBinsTmp, 0);
		unsigned int currInd = 0;
		double tmpMean = 0;
		for (unsigned int i = 0 ; i < sortVals.size() ; ++i)
		{
			tmpMean += sortVals[i];
			++bins[currInd];
			if (currInd < floor((sortVals[i] - sortVals[0]) / binSize))
			{
				vals[currInd] = tmpMean / bins[currInd];
				tmpMean = 0;
				currInd = floor((sortVals[i] - sortVals[0]) / binSize);
			}
		}
		unsigned int maxInd = 0;
		for (unsigned int i = 0 ; i < bins.size() ; ++i)
			if (bins[maxInd] < bins[i])
				maxInd = i;


		std::vector<double> tmpBins(bins);
		unsigned int binInd;
		for (unsigned int i = 0 ; i <= maxInd ; ++i)
		{
			tmpBins[i] = 0;
			binInd = std::min((unsigned int)tmpBins.size()-1, 2*maxInd - i);
			tmpBins[binInd] = std::min(tmpBins[binInd], 
				std::max(0.0, bins[binInd] - bins[i]));
		}
		// Compute new mean
		double newMean = 0;
		double totVals = 0;
		for (unsigned int i = 0 ; i < tmpBins.size() ; ++i)
		{
			totVals += tmpBins[i];
			newMean += tmpBins[i] * vals[i];
		}
		newMean /= totVals;

		threshVal = newMean - ComputeMean(sortVals) + 
			stdDevCoeff * ComputeStdDev(sortVals);
	}

	// Apply threshold
	std::vector<std::vector<double> > threshTopo(network.size(), 
		std::vector<double>(network.size(), 0));
	ApplyThresholdOnMat(values, threshTopo, threshVal);

	// Construct the actual network
	for (unsigned int i = 0 ; i < network.size() ; ++i)
	{
		for (unsigned int j = (dirLinks ? 0 : (i + 1)) ; j < network.size() ; ++j)
		{
			if (threshTopo[i][j] > 0)
			{
				network.SetAbstractEdge(i, j, factory.Create());
				if (not dirLinks)
					network.SetAbstractEdge(j, i, network.GetAbstractEdge(i, j));
			}
		}
	}

	network.SetDirected(dirLinks);
	return NetworkConstructStrat::BuildNetwork(network, factory, saver);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler FunctionalTopoStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	if (useFixedThresh)
		params <= "FunctionalTopoConstrThresh", threshold;
	params <= "FunctionalTopoStdDevCoeff", stdDevCoeff;
	return params;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool FunctionalTopoStrat::LoadFromStream(std::ifstream & stream)
{
	bool ok = NetworkConstructStrat::LoadFromStream(stream);
	stream >> dirLinks;
	stream >> useFixedThresh;
	stream >> useAbsoluteThresh;
	stream >> threshold;
	stream >> stdDevCoeff;
	stream >> absoluteThresh;
	return ok and stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool FunctionalTopoStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = NetworkConstructStrat::SaveToStream(stream);
	stream << dirLinks << std::endl
		<< useFixedThresh << std::endl
		<< useAbsoluteThresh << std::endl
		<< threshold << std::endl
		<< stdDevCoeff << std::endl
		<< absoluteThresh << std::endl;
	return ok and stream.good();
}

//**********************************************************************
// Set values to the given matrix
//**********************************************************************
void FunctionalTopoStrat::SetValues(const std::vector<std::vector<double> > & vals) const
{
	if (not dirLinks)
	{
		// use a bidirectional version of the values matrix
		for (unsigned int i = 0 ; i < values.size() ; ++i)
			for (unsigned int j = i + 1 ; j < values[i].size() ; ++j)
			{
				values[i][j] = max(values[i][j], values[j][i]);
				values[j][i] = values[i][j];
			}
	}

	values = vals;
}

//**********************************************************************
//**********************************************************************
void FunctionalTopoStrat::SetThreshold(double thr, bool use)
{
	threshold = thr;
	useFixedThresh = use;
}

//********************************************************************//
//**** C O N S T R U C T   F R O M   C H I   S I M U L A T I O N *****//
//********************************************************************//

static std::vector<std::vector<double> > emptyMat;

//**********************************************************************
// Default constructor
//**********************************************************************
ConstrFromChISim::ConstrFromChISim(ParamHandler & h) :
	NetworkConstructStrat::NetworkConstructStrat(h), model(0)
{
	std::vector<std::string> paramNames = 
		h.getParam<std::vector<std::string> >("-constrFromSimParams", 0);
	std::vector<int> paramPos = 
		h.getParam<std::vector<int> >("-constrFromSimParams", 1);
	std::vector<std::string> paramVals = 
		h.getParam<std::vector<std::string> >("-constrFromSimParams", 2);
	
	modelClassName = h.getParam<std::string>("-constrFromSimModelName", 0);
	bool useSameMetr = h.getParam<bool>("-constrFromSimUseSameMetr", 0);
	if (useSameMetr)
		metricNames = h.getParam<std::vector<std::string> >("-aM");

	assert(paramNames.size() == paramPos.size());
	assert(paramPos.size() == paramVals.size());
	for (unsigned int i = 0 ; i < paramNames.size() ; ++i)
		params[paramNames[i].substr(1, paramNames[i].length())][paramPos[i]] = paramVals[i];
}

//**********************************************************************
// Special constructor
//**********************************************************************
ConstrFromChISim::ConstrFromChISim(const std::vector<std::string> _metr, bool _n) : 
	NetworkConstructStrat::NetworkConstructStrat(_n), 
	model(0), metricNames(_metr)
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ConstrFromChISim::ConstrFromChISim(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool ConstrFromChISim::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
	std::map<std::string, std::map<int, std::string> > paramsOldVal;
	// Create ChIModel each time
	if (model)
		delete model;

	ParamHandler & handler = ParamHandler::GlobalParams;
	for (std::map<std::string, std::map<int, std::string> >::const_iterator it = 
			params.begin() ; it != params.end() ; ++it)
	{
		for (std::map<int, std::string>::const_iterator it2 = it->second.begin() ;
				it2 != it->second.end() ; ++it2)
		{
			paramsOldVal[it->first][it2->first] = handler.getStringParam(it->first, it2->first);
			handler.SetVal(it->first, it2->second.c_str(), it2->first);
		}
	}

	// Create ChIModel
	if (AbstractFactory<ChIModel>::Factories[modelClassName])
		model = AbstractFactory<ChIModel>::Factories[modelClassName]->Create();

	// Add topo funct metrics
	TransferEntropyMetric *TEMetr = new TransferEntropyMetric();
	model->AddMetric(TEMetr, true);

	// Add additional metric if needed
	Metric *tmpMetr = 0;
	for (unsigned int i = 0 ; i < metricNames.size() ; ++i)
	{
		if (AbstractFactory<Metric>::Factories[metricNames[i]])
		{
			tmpMetr = AbstractFactory<Metric>::Factories[metricNames[i]]->Create();
			tmpMetr->AddOptStatParamName("ConstrFromSim_");
			model->AddMetric(tmpMetr, true);
		}
	}

	// Initializing model
	model->Initialize();
	// Simulation
	model->Simulate(saver("ConstrFromSimData"));

	// Get functionnal metric values (TE or CC)

	const AbstractNetwork * tmpNet = TEMetr->GetFunctTopo(
		TransferEntropyMetric::TransEntrFunctTopoName);
	bool ok = false;
	if (tmpNet)
	{
		bool isDir = tmpNet->IsDirected();
		network.SetDirected(isDir);
		// Construct the actual network
		for (unsigned int i = 0 ; i < network.size() ; ++i)
		{
			for (unsigned int j = (isDir ? 0 : i + 1) ; j < network.size() ; ++j)
			{
				if (tmpNet->AreConnected(i, j))
				{
					network.SetAbstractEdge(i, j, factory.Create());
					if (not isDir)
						network.SetAbstractEdge(j, i, network.GetAbstractEdge(i, j));
				}
			}
		}
		ok = true;
	}

	// Restore old parameters
	for (std::map<std::string, std::map<int, std::string> >::const_iterator it = 
			paramsOldVal.begin() ; it != paramsOldVal.end() ; ++it)
		for (std::map<int, std::string>::const_iterator it2 = it->second.begin() ;
				it2 != it->second.end() ; ++it2)
			handler.SetVal(it->first, it2->second.c_str(), it2->first);


	return ok and NetworkConstructStrat::BuildNetwork(network, factory, saver);
}

//**********************************************************************
// 
//**********************************************************************
std::vector<Metric *> ConstrFromChISim::GetAllMetrics() const
{
	return model ? model->GetAllMetrics() : std::vector<Metric *>();
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler ConstrFromChISim::BuildModelParamHandler() 
{
	ParamHandler params;
	return params;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool ConstrFromChISim::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool ConstrFromChISim::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	return stream.good();
}

//********************************************************************//
//************* L O A D   F R O M   F I L E   S T R A T **************//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
FromFileConstrStrat::FromFileConstrStrat(ParamHandler & h) :
	NetworkConstructStrat::NetworkConstructStrat(h)
{
	path = h.getParam<std::string>("-frmFilePath", 0);
	newPath = path;
	directed = h.getParam<bool>("-frmFileDirected", 0);
	newDirected = directed;
	useSparse = h.getParam<bool>("-frmFileUseSparse", 0);
	useValuesAsStrengths = h.getParam<bool>("-frmFileUseValAsStrengths", 0);
}

//**********************************************************************
// Special constructor
//**********************************************************************
FromFileConstrStrat::FromFileConstrStrat(std::string _p, bool _us, bool _d, bool _n, bool _uvas) : 
NetworkConstructStrat::NetworkConstructStrat(_n), path(_p), newPath(_p), 
directed(_d), newDirected(_d), useSparse(_us), useValuesAsStrengths(_uvas)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
FromFileConstrStrat::FromFileConstrStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool FromFileConstrStrat::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
	path = newPath;
	directed = newDirected;
	std::ifstream stream(path.c_str());
	assert(stream.good());
	std::vector<std::vector<double> > adjMat;
	if (useSparse)
		adjMat = LoadSparseMatrix(stream);
	else
		adjMat = LoadFullMatrix(stream);

	assert(network.size() == adjMat.size());
	for (unsigned int i = 0 ; i < adjMat.size() ; ++i)
		for (unsigned int j = (directed ? 0 : i + 1) ; j < adjMat[i].size() ; ++j)
		{
			if (adjMat[i][j] > 0)
			{
				network.SetAbstractEdge(i, j, factory.Create());
				// Allows to set weighted links
				if (useValuesAsStrengths)
				{
					CouplingFunction *tmp = dynamic_cast<CouplingFunction*>(network.GetAbstractEdge(i, j));
					assert(tmp);
					tmp->ChangeStrength(adjMat[i][j]);
					tmp->SetConstStrength(true);
				}
				if (not directed)
					network.SetAbstractEdge(j, i, network.GetAbstractEdge(i, j));
			}
		}
	network.SetDirected(directed);
	return NetworkConstructStrat::BuildNetwork(network, factory, saver);
}

//**********************************************************************
// StaticConstrStrat method, returns whether the network needs to be 
// rebuilt
//**********************************************************************
bool FromFileConstrStrat::NeedsToBeconstructed() const
{
	return ((path != newPath) or (directed != newDirected));
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler FromFileConstrStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "ConstrFromFilePath", newPath;
	params <= "ConstrFromFileDirected", newDirected;
	params <= "ConstrFromFileUseSparse", useSparse;
	params <= "ConstrFromFileUseValuesAsStrengths", useValuesAsStrengths;
	return params;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool FromFileConstrStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> path;
	if (path == "NULL")
		path = "";
	newPath = path;
	stream >> directed;
	stream >> useSparse;
	stream >> useValuesAsStrengths;
	newDirected = directed;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool FromFileConstrStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	if (not newPath.empty())
		stream << newPath << std::endl;
	else
		stream << "NULL" << std::endl;
	stream << directed << std::endl;
	stream << useSparse << std::endl;
	stream << useValuesAsStrengths << std::endl;
	return stream.good();
}

//********************************************************************//
//******** L O A D   F R O M   F I L E   L I S T   S T R A T *********//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
FileListConstrStrat::FileListConstrStrat(ParamHandler & h) :
	FromFileConstrStrat::FromFileConstrStrat(h), posInList(0)
{
	listPath = h.getParam<std::string>("-frmFileListPath", 0);
	newListPath = path;
	if (not listPath.empty())
		LoadFileListFromPath(listPath);
}

//**********************************************************************
// Special constructor
//**********************************************************************
FileListConstrStrat::FileListConstrStrat(std::string _lp, bool _us, bool _d, bool _n) : 
	FromFileConstrStrat::FromFileConstrStrat("", _us, _d, _n), listPath(_lp), newListPath(_lp), 
	posInList(0)
{
	if (not listPath.empty())
		LoadFileListFromPath(listPath);
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
FileListConstrStrat::FileListConstrStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool FileListConstrStrat::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
TRACE(newListPath)
	if (listPath != newListPath)
		LoadFileListFromPath(newListPath);
	assert(posInList < fileList.size());
	this->newPath = fileList[posInList++];

	return FromFileConstrStrat::BuildNetwork(network, factory, saver);
}

//**********************************************************************
// StaticConstrStrat method, returns whether the network needs to be 
// rebuilt
//**********************************************************************
bool FileListConstrStrat::NeedsToBeconstructed() const
{
	return (listPath != newListPath) or (posInList < fileList.size());
}

//**********************************************************************
// Load File List from given file
//**********************************************************************
void FileListConstrStrat::LoadFileListFromPath(std::string _lp) const
{
	fileList.clear();
	posInList = 0;
	listPath = _lp;
	newListPath = listPath;

	std::ifstream stream(listPath.c_str());
	assert(stream.good());
	unsigned int nbLines;
	std::string tmpStr;

	stream >> nbLines;
	for (unsigned int i = 0 ; i < nbLines ; ++i)
	{
		stream >> tmpStr;
		fileList.push_back(tmpStr);
	}
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler FileListConstrStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params += FromFileConstrStrat::BuildModelParamHandler();
	params <= "ConstrFromListFilePath", newListPath;
	params <= "FullAdjacencyMatrix", newListPath;
	return params;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool FileListConstrStrat::LoadFromStream(std::ifstream & stream)
{
	bool ok = true;
	ok &= FromFileConstrStrat::LoadFromStream(stream);
	stream >> listPath;
	if (listPath == "NULL")
		listPath = "";
	newListPath = listPath;
	if (not listPath.empty())
		LoadFileListFromPath(listPath);
	return ok and stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool FileListConstrStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = true;
	ok &= FromFileConstrStrat::SaveToStream(stream);
	if (not listPath.empty())
		stream << listPath << std::endl;
	else
		stream << "NULL" << std::endl;
	return ok and stream.good();
}

//********************************************************************//
//******* R A N D O M   I S O L A T E D   C E L L S   S T R A T ******//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
RandomIsolatedCells::RandomIsolatedCells(ParamHandler & h) :
	NetworkConstructStrat::NetworkConstructStrat(h)
{
	ratio = h.getParam<double>("-isolatedCells", 0);
	subStratName = h.getParam<std::string>("-isolatedCells", 1);
	assert(AbstractFactory<NetworkConstructStrat>::Factories[subStratName]);
	subStrat = AbstractFactory<NetworkConstructStrat>::
		Factories[subStratName]->Create();
}

//**********************************************************************
// Special constructor
//**********************************************************************
RandomIsolatedCells::RandomIsolatedCells(double _rat, std::string _sN, bool _n) :
	NetworkConstructStrat::NetworkConstructStrat(_n), ratio(_rat), subStratName(_sN)
{
	assert(AbstractFactory<NetworkConstructStrat>::Factories[_sN]);
	subStrat = AbstractFactory<NetworkConstructStrat>::
		Factories[_sN]->Create();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
RandomIsolatedCells::RandomIsolatedCells(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler RandomIsolatedCells::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "RandIsoCellsRatio", ratio;
	params <= "RandIsoCellsSubStratName", subStratName;
	return params;
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool RandomIsolatedCells::BuildNetwork(AbstractNetwork & network, 
	AbstractFactory<NetworkEdge> & factory, ResultSaver saver) const
{
	// If substrategy class name has changed
	if (not subStrat or (subStrat->GetClassName() != subStratName))
	{
		delete subStrat;
		assert(AbstractFactory<NetworkConstructStrat>::Factories[subStratName]);
		subStrat = AbstractFactory<NetworkConstructStrat>::
			Factories[subStratName]->Create();
	}
	assert(subStrat);
	assert(network.size() > 1);

	// Build the network with substrat
	assert(subStrat->BuildNetwork(network, factory, saver));

	// Delete links to and from randomly chosen nodes
	for (unsigned int i = 0 ; i < network.size() ; ++i)
		if (TrueWithProba(ratio))
			for (unsigned int j = 0 ; j < network.size() ; ++j)
				if ((i != j) and (network.GetAbstractEdge(i, j) or 
					network.GetAbstractEdge(j, i)))
				{
					delete network.GetAbstractEdge(i, j);
					if (network.GetAbstractEdge(i, j) != network.GetAbstractEdge(j, i))
						delete network.GetAbstractEdge(j, i);

					network.SetAbstractEdge(i, j, 0);
					network.SetAbstractEdge(j, i, 0);
				}

	network.SetDirected(subStrat->GetMakeDir());
	return NetworkConstructStrat::BuildNetwork(network, factory, saver);
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool RandomIsolatedCells::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> ratio;
	stream >> subStratName;
	assert(AbstractFactory<NetworkConstructStrat>::Factories[subStratName]);
	subStrat = AbstractFactory<NetworkConstructStrat>::
		Factories[subStratName]->CreateFromStream(stream);

	return subStrat and stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool RandomIsolatedCells::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << ratio << std::endl;
	stream << subStrat->GetClassName() << std::endl;
	assert(subStrat);
	subStrat->SaveToStream(stream);
	return stream.good();
}

//********************************************************************//
//******* S T A R   L I K E   A S T R O C Y T E   N E T W O R K ******//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
StarLikeAstroNetConstrStrat::StarLikeAstroNetConstrStrat(
	ParamHandler & h) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(h)
{
	somaRadius = h.getParam<unsigned int>("-StarLikeNet", 0);
	somaCouplStr = h.getParam<double>("-StarLikeNet", 1);
	nbBranches = h.getParam<unsigned int>("-StarLikeBuild", 3);
}

//**********************************************************************
// Special constructor
//**********************************************************************
StarLikeAstroNetConstrStrat::StarLikeAstroNetConstrStrat(unsigned int _sr, 
	unsigned int _nbb, double _scs, bool _n) : 
	SpatialNetConstrStrat::SpatialNetConstrStrat(_n), somaRadius(_sr), 
	nbBranches(_nbb), somaCouplStr(_scs)
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
StarLikeAstroNetConstrStrat::StarLikeAstroNetConstrStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler StarLikeAstroNetConstrStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "StarLikeSomaRadius", somaRadius;
	params <= "StarLikeNbBranches", nbBranches;
	params <= "StarLikeSomaCouplingStrength", somaCouplStr;
	return params;
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool StarLikeAstroNetConstrStrat::buildSpatialNetwork(AbstractSpatialNetwork & network,
	AbstractFactory<NetworkEdge> & factory) const
{
	std::vector<std::vector<std::vector<unsigned int> > > subElemStruct;
	std::vector<unsigned int> cellCenters;

	// Parse sub-element structure from tags
	for(unsigned int i = 0 ; i < network.size() ; ++i)
	{
		unsigned int tmp = network.GetNodeTag(i);
		unsigned int astrNum = (tmp & 0x7F0000) >> 16;
		unsigned int procNum = (tmp & 0xFF00) >> 8;
		unsigned int subNum  = (tmp & 0xFF);
		if (astrNum > 0)
		{
			for(unsigned int j = subElemStruct.size() ; j < astrNum ; ++j)
			{
				subElemStruct.push_back(std::vector<std::vector<unsigned int> >(
							nbBranches, std::vector<unsigned int>()));
				cellCenters.push_back(0);
			}
			// If the subelement is part of a process
			if (procNum > 0)
			{
				for(unsigned int j = subElemStruct[astrNum-1][procNum-1].size() ; j < (subNum + 1) ; ++j)
				{
					subElemStruct[astrNum-1][procNum-1].push_back(0);
				}
				subElemStruct[astrNum-1][procNum-1][subNum] = i;
			}
			else // if it is the cell center
			{
				cellCenters[astrNum - 1] = i;
			}
		}
	}

	// Add links 
	// For each Astro
	for(unsigned int i = 0 ; i < subElemStruct.size() ; ++i)
	{
		// Tag center as part of the soma
		network.SetNodeTag(cellCenters[i], 0x1000000 + 
			network.GetNodeTag(cellCenters[i]));

		// For each process
		for(unsigned int j = 0 ; j < subElemStruct[i].size() ; ++j)
		{
			// Add links from center to begining of each process
			network.SetAbstractEdge(cellCenters[i], subElemStruct[i][j][0], factory.Create());
			network.SetAbstractEdge(subElemStruct[i][j][0], cellCenters[i], 
					network.GetAbstractEdge(cellCenters[i], subElemStruct[i][j][0]));
			setEdgeToSpecificCoupling(network.GetAbstractEdge(cellCenters[i], 
				subElemStruct[i][j][0]), somaCouplStr);
			// Add links inside processes
			for(unsigned int k = 0 ; k < (subElemStruct[i][j].size()-1) ; ++k)
			{
				network.SetAbstractEdge(subElemStruct[i][j][k], subElemStruct[i][j][k+1], factory.Create());
				network.SetAbstractEdge(subElemStruct[i][j][k+1], subElemStruct[i][j][k], 
						network.GetAbstractEdge(subElemStruct[i][j][k], subElemStruct[i][j][k+1]));
				// If the edge is part of the soma
				if ((k + 1) < somaRadius)
					setEdgeToSpecificCoupling(network.GetAbstractEdge(
						subElemStruct[i][j][k], subElemStruct[i][j][k+1]), somaCouplStr);
			}
		}
		// Add links inside soma
		for(unsigned int r = 0 ; r < somaRadius ; ++r)
		{
			for(unsigned int j = 0 ; j < subElemStruct[i].size() ; ++j)
			{
				network.SetAbstractEdge(subElemStruct[i][j][r], 
						subElemStruct[i][(j+1)%nbBranches][r], factory.Create());
				network.SetAbstractEdge(subElemStruct[i][(j+1)%nbBranches][r],
						subElemStruct[i][j][r], network.GetAbstractEdge(subElemStruct[i][j][r],
							subElemStruct[i][(j+1)%nbBranches][r]));
				setEdgeToSpecificCoupling(network.GetAbstractEdge(subElemStruct[i][j][r],
						subElemStruct[i][(j+1)%nbBranches][r]), somaCouplStr);
				network.SetNodeTag(subElemStruct[i][j][r], 0x1000000 + 
					network.GetNodeTag(subElemStruct[i][j][r]));
			}
		}
	}
	return true;
}

//**********************************************************************
//**********************************************************************
void StarLikeAstroNetConstrStrat::setEdgeToSpecificCoupling(
	NetworkEdge *edge, double cpl) const
{
	CouplingFunction *cplfct = dynamic_cast<CouplingFunction *>(edge);
	assert(cplfct);
	cplfct->ChangeStrength(cpl);
	cplfct->SetConstStrength(true);
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool StarLikeAstroNetConstrStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> somaRadius;
	stream >> nbBranches;
	stream >> somaCouplStr;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool StarLikeAstroNetConstrStrat::SaveToStream(std::ofstream & stream) const
{
	NetworkConstructStrat::SaveToStream(stream);
	stream << somaRadius << std::endl
		<< nbBranches << std::endl
		<< somaCouplStr << std::endl;
	return stream.good();
}

//********************************************************************//
//***** S H E L L   S T R U C T   S C R A M B L E R   S T R A T ******//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
ShellStructScramblerStrat::ShellStructScramblerStrat(
	ParamHandler & h) :
	SpatialNetConstrStrat::SpatialNetConstrStrat(h)
{
	subConstrStratName = h.getParam<std::string>("-ShellScrambler", 0);
	startNode = h.getParam<int>("-ShellScrambler", 1);
	assert(AbstractFactory<NetworkConstructStrat>::
		Factories[subConstrStratName]);
	NetworkConstructStrat *tmp = AbstractFactory<NetworkConstructStrat>::
		Factories[subConstrStratName]->Create();
	strat = dynamic_cast<SpatialNetConstrStrat *>(tmp);
	assert(strat);
}

//**********************************************************************
// Special constructor
//**********************************************************************
ShellStructScramblerStrat::ShellStructScramblerStrat(std::string _ss, 
	int _sn, bool _n) : 
	SpatialNetConstrStrat::SpatialNetConstrStrat(_n), subConstrStratName(_ss), 
	startNode(_sn)
{
	assert(AbstractFactory<NetworkConstructStrat>::
		Factories[subConstrStratName]);
	NetworkConstructStrat *tmp = AbstractFactory<NetworkConstructStrat>::
		Factories[subConstrStratName]->Create();
	strat = dynamic_cast<SpatialNetConstrStrat *>(tmp);
	assert(strat);
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ShellStructScramblerStrat::ShellStructScramblerStrat(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
ShellStructScramblerStrat::~ShellStructScramblerStrat()
{
	delete strat;
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler ShellStructScramblerStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "ShellScramblerSubStratName", subConstrStratName;
	params <= "ShellScramblerStartNode", startNode;
	// If we need to reload the sub constr strat
	if (strat->GetClassName() != subConstrStratName)
	{
		delete strat;
		assert(AbstractFactory<NetworkConstructStrat>::
			Factories[subConstrStratName]);
		NetworkConstructStrat *tmp = AbstractFactory<NetworkConstructStrat>::
			Factories[subConstrStratName]->Create();
		strat = dynamic_cast<SpatialNetConstrStrat *>(tmp);
		assert(strat);
	}
	params += strat->BuildModelParamHandler();
	return params;
}

//**********************************************************************
// Builds the network topology with given factory
//**********************************************************************
bool ShellStructScramblerStrat::buildSpatialNetwork(AbstractSpatialNetwork & network,
	AbstractFactory<NetworkEdge> & factory) const
{
	strat->BuildNetwork(network, factory, ResultSaver::NullSaver);

	std::set<unsigned int> currShell, nextShell;
	currShell.insert(startNode);
	std::vector<std::vector<double> > adjMat = network.GetAdjMat();
	std::vector<std::vector<double> > newAdjMat = adjMat; 

	// Randomize shell structure
	do
	{
		nextShell = getNextShell(currShell, adjMat);
		// Treat inter-shell links
		randomizeSubSquare(currShell, nextShell, newAdjMat);
		// Treat intra-shell links
		randomizeSubSquare(nextShell, nextShell, newAdjMat);
		// Remove inter and intra-shell links in original adjacency matrix
		clearLinksBetween(currShell, nextShell, adjMat);
		clearLinksBetween(nextShell, nextShell, adjMat);
		currShell = nextShell;
	} while (! nextShell.empty());

	// Clear old net
	for(unsigned int i = 0 ; i < network.size() ; ++i)
	{
		for (unsigned int j = (i+1) ; j < network.size() ; ++j)
		{
			delete network.GetAbstractEdge(i, j);
			network.SetAbstractEdge(i, j, 0);
			network.SetAbstractEdge(j, i, 0);
		}
	}

	// Copy new net
	for(unsigned int i = 0 ; i < newAdjMat.size() ; ++i)
	{
		for (unsigned int j = (i+1) ; j < newAdjMat[i].size() ; ++j)
		{
			if (newAdjMat[i][j] > 0)
			{
				network.SetAbstractEdge(i, j, factory.Create());
				network.SetAbstractEdge(j, i, network.GetAbstractEdge(i, j));
			}
		}
	}
	return true;
}

//**********************************************************************
//**********************************************************************
void ShellStructScramblerStrat::randomizeSubSquare(
	const std::set<unsigned int> & pop1, const std::set<unsigned int> & pop2,
	std::vector<std::vector<double> > & adjMat) const
{
	std::vector<unsigned int> pop1v(pop1.begin(), pop1.end());
	std::vector<unsigned int> pop2v(pop2.begin(), pop2.end());

	// Classify degrees of pop1
	std::vector<std::vector<unsigned int> > classifByInDegPop1(pop2.size() + 1, 
		std::vector<unsigned int>());
	std::vector<unsigned int> pop1Deg;
	for(unsigned int i = 0 ; i < pop1v.size() ; ++i)
	{
		unsigned int tmpVal = 0;
		for (unsigned int j = 0 ; j < pop2.size() ; ++j)
			if (adjMat[pop1v[i]][pop2v[j]] > 0)
				++tmpVal;
		classifByInDegPop1[tmpVal].push_back(pop1v[i]);
		pop1Deg.push_back(tmpVal);
	}

	// Classify degrees of pop2
	std::vector<std::vector<unsigned int> > classifByInDegPop2(pop1.size() + 1, 
		std::vector<unsigned int>());
	std::vector<unsigned int> pop2Deg;
	for(unsigned int i = 0 ; i < pop2v.size() ; ++i)
	{
		unsigned int tmpVal = 0;
		for (unsigned int j = 0 ; j < pop1.size() ; ++j)
			if (adjMat[pop2v[i]][pop1v[j]] > 0)
				++tmpVal;
		classifByInDegPop2[tmpVal].push_back(pop2v[i]);
		pop2Deg.push_back(tmpVal);
	}
	

	// Switch lines
	for(unsigned int i = 0 ; i < pop1v.size() ; ++i)
	{
		unsigned int destLineInd = floor(UnifRand() * (double) classifByInDegPop1[pop1Deg[i]].size());
		switchLines(pop1v[i], classifByInDegPop1[pop1Deg[i]][destLineInd], pop2, adjMat);
	}
	// Switch columns
	for(unsigned int i = 0 ; i < pop2v.size() ; ++i)
	{
		unsigned int destColumnInd = floor(UnifRand() * (double) classifByInDegPop2[pop2Deg[i]].size());
		switchColumns(pop2v[i], classifByInDegPop2[pop2Deg[i]][destColumnInd], pop1, adjMat);
	}
}

//**********************************************************************
//**********************************************************************
void ShellStructScramblerStrat::switchLines(unsigned int l1, unsigned int l2,
	const std::set<unsigned int> & subcols, 
	std::vector<std::vector<double> > & adjMat) const
{
	for (std::set<unsigned int>::const_iterator it = subcols.begin() ; it != subcols.end() ; ++it)
	{
		switchCases(l1, *it, l2, *it, adjMat);
	}
}

//**********************************************************************
//**********************************************************************
void ShellStructScramblerStrat::switchColumns(unsigned int c1, unsigned int c2, 
	const std::set<unsigned int> & sublines, 
	std::vector<std::vector<double> > & adjMat) const
{
	for (std::set<unsigned int>::const_iterator it = sublines.begin() ; it != sublines.end() ; ++it)
	{
		switchCases(*it, c1, *it, c2, adjMat);
	}
}

//**********************************************************************
//**********************************************************************
void ShellStructScramblerStrat::switchCases(
	unsigned int l1, unsigned int c1,
	unsigned int l2, unsigned int c2, 
	std::vector<std::vector<double> > & adjMat) const
{
	double tmpVal = adjMat[l2][c2];
	adjMat[l2][c2] = adjMat[l1][c1];
	adjMat[l1][c1] = tmpVal;
	tmpVal = adjMat[c2][l2];
	adjMat[c2][l2] = adjMat[c1][l1];
	adjMat[c1][l1] = tmpVal;
}

//**********************************************************************
// Returns the next shell provided that all links to previous and 
// intra-links were removed
//**********************************************************************
std::set<unsigned int> ShellStructScramblerStrat::getNextShell(
	const std::set<unsigned int> & currShell, 
	const std::vector<std::vector<double> > & adjMat) const
{
	std::set<unsigned int> nextShell;
	for (std::set<unsigned int>::const_iterator it = currShell.begin() ; it != currShell.end() ; ++it)
	{
		std::set<unsigned int> neighbs = getNeighbs(*it, adjMat);
		nextShell.insert(neighbs.begin(), neighbs.end());
	}
	return nextShell;
}

//**********************************************************************
// Return neighbors of a node
//**********************************************************************
std::set<unsigned int> ShellStructScramblerStrat::getNeighbs(unsigned int i, 
	const std::vector<std::vector<double> > & adjMat) const
{
	std::set<unsigned int> neighbs;
	for(unsigned int j = 0 ; j < adjMat[i].size() ; ++j)
	{
		if (adjMat[i][j] > 0)
			neighbs.insert(j);
	}
	return neighbs;
}

//**********************************************************************
// Clear all the links between two node populations
//**********************************************************************
void ShellStructScramblerStrat::clearLinksBetween(
	const std::set<unsigned int> & pop1, const std::set<unsigned int> & pop2, 
	std::vector<std::vector<double> > & adjMat) const
{
	for (std::set<unsigned int>::const_iterator it = pop1.begin() ; it != pop1.end() ; ++it)
	{
		for (std::set<unsigned int>::const_iterator it2 = pop2.begin() ; it2 != pop2.end() ; ++it2)
		{
			adjMat[*it][*it2] = 0;
			adjMat[*it2][*it] = 0;
		}
	}
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool ShellStructScramblerStrat::LoadFromStream(std::ifstream & stream)
{
	NetworkConstructStrat::LoadFromStream(stream);
	stream >> subConstrStratName;
	assert(AbstractFactory<NetworkConstructStrat>::
		Factories[subConstrStratName]);
	NetworkConstructStrat *tmp = AbstractFactory<NetworkConstructStrat>::
		Factories[subConstrStratName]->CreateFromStream(stream);
	strat = dynamic_cast<SpatialNetConstrStrat *>(tmp);
	assert(strat);

	stream >> startNode;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool ShellStructScramblerStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = NetworkConstructStrat::SaveToStream(stream);
	stream << subConstrStratName << std::endl;
	ok &= strat->SaveToStream(stream);
	stream << startNode << std::endl;
	return ok and stream.good();
}
