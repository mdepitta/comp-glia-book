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

#ifndef NETWORK_H
#define NETWORK_H

#include <vector>

#include "Savable.h"
#include "NetworkConstructStrat.h"
#include "AbstractFactory.h"
#include "NetworkMetrics.h"

#define SPECIAL_END_INDEX -1

namespace AstroModel
{
/**********************************************************************/
/* Abstract Base class (interface)                                    */
/**********************************************************************/
	class AbstractNetwork
	{
	public:
		virtual ~AbstractNetwork() {}
		virtual std::string GetClassName() const = 0;
		virtual unsigned int size() const = 0;
		virtual unsigned int size(unsigned int i) const = 0;
		virtual unsigned int GetNeededSize() const = 0;
		virtual bool AreConnected(unsigned int i, unsigned int j) const = 0;
		virtual double GetNodeDegree(unsigned int ind) const = 0;
		virtual const std::vector<unsigned int> & GetNeighbors(unsigned int i) const = 0;
		virtual NetworkEdge * GetAbstractEdge(unsigned int i, unsigned int j) const = 0;
		virtual void SetAbstractEdge(unsigned int i, unsigned int j, NetworkEdge * edge) = 0;
		virtual void SetDirected(bool _b) = 0;
		virtual bool IsNodeOnEdge(unsigned int i) const = 0;
		virtual bool IsSpatial() const = 0;
		virtual bool IsDirected() const = 0;
		virtual bool AddMetric(Metric *, bool) = 0; 
		virtual std::vector<Metric *> GetAllMetrics() const = 0;
		virtual const NetworkConstructStrat & GetConstructStrat() const = 0;
		virtual bool BuildNetwork(ResultSaver) = 0;
		virtual unsigned int GetNbEdges(bool dir = true) const = 0;
		virtual std::vector<std::vector<double> > GetAdjMat() const = 0;
		virtual unsigned long int GetNodeTag(unsigned int i) const = 0;
		virtual void SetNodeTag(unsigned int i, unsigned long int tag) = 0;
	};

/**********************************************************************/
/* Base class                                                         */
/**********************************************************************/
	template <typename LinkType>
	class Network : public SaveAndLoadFromStream, public virtual AbstractNetwork
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		Network(ParamHandler & h = ParamHandler::GlobalParams) : 
			hasBeenBuilt(false), isDirected(false)
		{
			netSize = (unsigned int)h.getParam<int>("-N");
			assert(netSize > 0);
			W = std::vector<std::vector<LinkType *> >(netSize, 
				std::vector<LinkType *>(netSize, 0));
			constructStratName = h.getParam<std::string>("-Construct", 0);
			assert(AbstractFactory<NetworkConstructStrat>::
				Factories[constructStratName]);
			construct = AbstractFactory<NetworkConstructStrat>::
				Factories[constructStratName]->Create();
			assert(construct);
			freeConstruct = true;
			netEdgeClassName = h.getParam<std::string>("-Coupling");
			nodeTags = std::vector<unsigned long int>(netSize, 0);
		}
		// Standard Constructor
		Network(unsigned int nbCells, NetworkConstructStrat *_c, 
			bool free = false, std::string _necn = "SigmoidFunction")
			: W(std::vector<std::vector<LinkType *> >(nbCells, 
				std::vector<LinkType *>(nbCells, 0))), 
			construct(_c), freeConstruct(free), netSize(nbCells), 
			netEdgeClassName(_necn), hasBeenBuilt(false), isDirected(false),
			nodeTags(std::vector<unsigned long int>(nbCells, 0))
		{
			assert(construct);
			constructStratName = _c->GetClassName();
		}

		// Constructor from stream
		Network(std::ifstream & stream):
			construct(0), freeConstruct(false), hasBeenBuilt(false), isDirected(false)
		{
			LoadFromStream(stream);
		}

		// Destructor
		virtual ~Network()
		{
			if (freeConstruct and construct)
				delete construct;
			ClearNetwork();
		}

		//===========================================================||
		// Network modification methods (create, clear, ...)         ||
		//===========================================================||
		// Builds a network using the network building strategy
		virtual bool BuildNetwork(ResultSaver saver = ResultSaver::NullSaver)
		{
			// Only rebuild the network if it hasn't been built or if it needs to be rebuilt
			StaticConstrStrat *tmpstrat = dynamic_cast<StaticConstrStrat*>(construct);
			if (not (hasBeenBuilt and (tmpstrat and not tmpstrat->NeedsToBeconstructed())))
			{
				ClearNetwork();
				if (W.size() != netSize)
				{
					assert(netSize);
					W = std::vector<std::vector<LinkType*> >(netSize, 
						std::vector<LinkType *>(netSize, 0));
					nodeTags = std::vector<unsigned long int>(netSize, 0);
				}

				assert(AbstractFactory<NetworkEdge>::Factories[netEdgeClassName]);
				if (construct)
				{
					hasBeenBuilt = construct->BuildNetwork(*this, 
						*AbstractFactory<NetworkEdge>::Factories[netEdgeClassName], saver);
					buildNeighbors();
				}
			}
			else // Just update the links in case parameters were changed
			{
				for (unsigned int i = 0 ; i < W.size() ; ++i)
					for (unsigned int j = i ; j < W[i].size() ; ++j)
						if (W[i][j])
						{
							if (W[j][i] and (W[i][j] != W[j][i]))
							{
								W[i][j]->UpdateEdge();
								W[j][i]->UpdateEdge();
							}
							else
								W[i][j]->UpdateEdge();
						}
			}
			return hasBeenBuilt;
		}

		// Clears the network and free its contents
		virtual void ClearNetwork()
		{
			for (unsigned int i = 0 ; i < W.size() ; ++i)
				for (unsigned int j = i ; j < W[i].size() ; ++j)
					if (W[i][j])
					{
						delete W[i][j];
						if (W[j][i] != W[i][j])
							delete W[j][i];
						W[i][j] = 0;
						W[j][i] = 0;
					}

			neighbors.clear();

			hasBeenBuilt = false;
		}

		// Initializes the network
		virtual void Initialize()
		{
			metrics.InitializeMetricsDefault();
		}

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the cell from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = true;
			ClearNetwork();
			W.clear();
			// Network
			unsigned int w;
			int tempInd;
			std::string tempName;

			stream >> netSize;
			for (unsigned int i = 0 ; i < netSize ; ++i)
			{
				W.push_back(std::vector<LinkType *>());
				stream >> w;
				stream >> tempInd;
				for (unsigned int j = 0 ; j < w ; ++j)
				{
					if (static_cast<int>(j) != tempInd)
						W.back().push_back(0);
					else
					{
						stream >> tempName;
						LinkType *tempEdge = dynamic_cast<LinkType*>(
							AbstractFactory<NetworkEdge>::Factories[tempName]->CreateFromStream(stream));
						W.back().push_back(tempEdge);
						if (tempInd != SPECIAL_END_INDEX)
							stream >> tempInd;
					}
				}
			}
			// Construction strategy
			if (construct and freeConstruct)
				delete construct;
			stream >> tempName;
			constructStratName = tempName;
			construct = AbstractFactory<NetworkConstructStrat>::Factories[tempName]
							->CreateFromStream(stream);
			freeConstruct = true;

			stream >> netEdgeClassName;

			// Metrics
			metrics.FreeAndClean();
			ok &= metrics.LoadFromStream(stream);

			stream >> hasBeenBuilt;

			ok &= LoadVectorFromStream(nodeTags, stream);

			return ok and stream.good() and not stream.eof();
		}

		// Saves the cell to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = true;
			// Network
			stream << W.size() << std::endl;
			for (unsigned int i = 0 ; i < W.size() ; ++i)
			{
				stream << W[i].size() << std::endl;
				for (unsigned int j = 0 ; j < W[i].size() ; ++j)
				{
					if (W[i][j])
					{
						stream 
							<< j << std::endl
							<< W[i][j]->GetClassName() << std::endl;
						ok &= W[i][j]->SaveToStream(stream);
					}
				}
				stream << SPECIAL_END_INDEX << std::endl;
			}
			// Construction strategy
			stream << construct->GetClassName() << std::endl;
			ok &= construct->SaveToStream(stream);

			stream << netEdgeClassName << std::endl;

			// Metrics
			ok &= metrics.SaveToStream(stream);

			stream << hasBeenBuilt << std::endl;

			ok &= SaveVectorToStream(nodeTags, stream);

			return ok and stream.good();
		}

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Only compute metrics
		virtual bool ComputeMetrics()
		{
			return metrics.ComputeMetricsDefault(*this);
		}
		// Only save metrics
		virtual bool SaveMetrics(ResultSaver saver)
		{
			return metrics.SaveMetricsDefault(saver);
		}
		// Compute and Save desired metrics
		virtual bool ComputeAndSaveMetrics(ResultSaver saver)
		{
			return ComputeMetrics() and SaveMetrics(saver);
		}

		// Add a metric
		virtual bool AddMetric(Metric *_m, bool _f = false) 
		{
			return metrics.AddMetricAndDependencies(_m, _f, this);
		}

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler() 
		{
			ParamHandler params;
			params <= "NetworkConstructStrat", constructStratName;
			params <= "NetworkSize", netSize;
			if (construct)
			{
				if (constructStratName != construct->GetClassName())
				{
					if (freeConstruct)
						delete construct;
					assert(AbstractFactory<NetworkConstructStrat>::
						Factories[constructStratName]);
					construct = AbstractFactory<NetworkConstructStrat>::
						Factories[constructStratName]->Create();
					freeConstruct = true;
				}
				params += construct->BuildModelParamHandler();
			}
			// Network edge name and params
			params <= "NetworkEdgeClassName", netEdgeClassName;
			assert(AbstractFactory<NetworkEdge>::Factories[netEdgeClassName]);
			LinkType *tempLink = dynamic_cast<LinkType*>(AbstractFactory<NetworkEdge>::Factories[
				netEdgeClassName]->Create());
			assert(tempLink);
			params += tempLink->BuildModelParamHandler();
			delete tempLink;

			params += metrics.BuildModelParamHandler();
			return params;
		}

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		// Access operator to the link matrix
		virtual std::vector<LinkType*> & operator[](unsigned int i)
			{ return W[i]; }
		// Const version of acces operator
		virtual const std::vector<LinkType*> & operator[](unsigned int i) const
			{ return W[i]; }
		// Returns a const ref on neighbors
		inline const std::vector<std::vector<unsigned int> > & GetAllNeighbors() const
			{ return neighbors; }
		// Returns the neighbors of cell i
		inline const std::vector<unsigned int> & GetNeighbors(unsigned int i) const
			{ return neighbors[i]; }

		// Returns the size of the first vector
		virtual unsigned int size() const { return W.size(); }
		// Returns the size of the ith second vector
		virtual unsigned int size(unsigned int i) const { return W[i].size(); }
		// Returns the needed size for the network (can be different from the actual one)
		virtual unsigned int GetNeededSize() const { return netSize; }
		// Returns true of two cells are connected
		virtual bool AreConnected(unsigned int i, unsigned int j) const
		{ return W[i][j];}
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const
		{
			std::vector<Metric *> mTot, mTemp;
			if (construct)
				mTot = construct->GetAllMetrics();
			mTemp = metrics.GetMetricsRaw();
			for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
				mTot.push_back(mTemp[i]);
			return mTot;
		}
		// Returns the degree of a node
		virtual double GetNodeDegree(unsigned int ind) const
		{
			int tempDegree = 0;
			for (unsigned int i = 0 ; i < W.size() ; ++i)
				if (W[ind][i])
					++tempDegree;
			return tempDegree;
		}
		// Returns an abstract pointer to the edge between i and j
		virtual NetworkEdge * GetAbstractEdge(unsigned int i, unsigned int j) const
		{
			return W[i][j];
		}
		// Sets the edge between i and j to edge
		virtual void SetAbstractEdge(unsigned int i, unsigned int j, NetworkEdge * edge)
		{
			if (not edge)
				W[i][j] = 0;
			else
			{
				LinkType *tempEdge = dynamic_cast<LinkType*>(edge);
				if (tempEdge and (i < W.size()) and (j < W[i].size()))
				{
					if (W[i][j] and (W[i][j] != W[j][i]))
						delete W[i][j];
					W[i][j] = tempEdge;
				}
			}
		}
		virtual bool IsNodeOnEdge(unsigned int) const { return false; }
		virtual bool IsSpatial() const { return false; }
		virtual bool IsDirected() const { return isDirected; }
		virtual const NetworkConstructStrat & GetConstructStrat() const 
		{
			assert(construct);
			return *construct; 
		}
		virtual bool HasBeenBuilt() const { return hasBeenBuilt; }
		virtual void SetDirected(bool _b) { isDirected = _b; }
		virtual unsigned int GetNbEdges(bool dir = true) const
		{
			unsigned int res = 0;
			for (unsigned int i = 0 ; i < W.size() ; ++i)
				for (unsigned int j = (dir ? 0 : i + 1) ; j < W[i].size() ; ++j)
					if (i !=j)
						++res;
			return res;
		}
		virtual std::vector<LinkType*> GetAllocatedLinks(bool dir = true) const
		{
			std::vector<LinkType*> res;
			for (unsigned int i = 0 ; i < W.size() ; ++i)
				for (unsigned int j = (dir ? 0 : i + 1) ; j < W[i].size() ; ++j)
					if ((i !=j) and W[i][j])
						res.push_back(W[i][j]);
			return res;
		}
		virtual void SetEdgeClassName(std::string _nam)
		{
			netEdgeClassName = _nam;
		}
		virtual std::vector<std::vector<double> > GetAdjMat() const
		{
			std::vector<std::vector<double> > tmpAdj(W.size(), std::vector<double>(W.size(), 0));
			for (unsigned int i = 0 ; i < W.size() ; ++i)
				for (unsigned int j = 0 ; j < W[i].size() ; ++j)
					if (AreConnected(i, j))
						tmpAdj[i][j] = 1;
			return tmpAdj;
		}

		virtual unsigned long int GetNodeTag(unsigned int i) const
		{
			assert(i < nodeTags.size());
			return nodeTags[i];
		}
		virtual void SetNodeTag(unsigned int i, unsigned long int tag) 
		{
			assert(i < nodeTags.size());
			nodeTags[i] = tag;
		}

	protected:
		std::vector<std::vector<LinkType*> > W; // Link matrix
		std::string constructStratName;         // Name of the network construct strategy to use
		NetworkConstructStrat *construct;       // Network construction strategy
		bool freeConstruct;                     // Delete construct strategy ?
		unsigned int netSize;                   // Size of the network
		std::string netEdgeClassName;           // Network edge class name
		bool hasBeenBuilt;
		bool isDirected;

		std::vector<unsigned long int> nodeTags; // Gives the possibility to tag nodes differently

		// Links stored in an efficient structure
		std::vector<std::vector<unsigned int> > neighbors;

		//===========================================================||
		// Metrics                                                   ||
		//===========================================================||
		SortedMetrics<AbstractNetwork> metrics;

		// Build the neighbors index (efficient when the network is sparse)
		virtual void buildNeighbors()
		{
			neighbors.clear();
			for (unsigned int i = 0 ;  i < W.size() ; ++i)
			{
				neighbors.push_back(std::vector<unsigned int>());
				for (unsigned int j = 0 ; j < W[i].size() ; ++j)
					if (W[i][j])
						neighbors.back().push_back(j);
			}
		}

	};
}

#endif
