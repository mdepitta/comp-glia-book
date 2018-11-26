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

#ifndef SPATIALNETWORK_H
#define SPATIALNETWORK_H

#include "Network.h"
#include "utility.h"
#include "SpatialStructureBuilder.h"
#include "ParamHandler.h"

#include <vector>
#include <algorithm>

namespace AstroModel
{
	// Forward declarations
	class AbstractNetwork;
	template <typename LinkType> class Network;

/**********************************************************************/
/* Abstract Base class (interface)                                    */
/**********************************************************************/
	class AbstractSpatialNetwork : public virtual AbstractNetwork
	{
	public:
		virtual ~AbstractSpatialNetwork() {}
		virtual Position* & Pos(unsigned int i) = 0;
		virtual Position* const & Pos(unsigned int i) const = 0;
		virtual double GetDistanceBetween(unsigned int i, unsigned int j) const = 0;
		virtual unsigned int GetDim() const = 0;
		virtual void SetDim(unsigned int _d) = 0;
		virtual bool BuildSpatialStructure() = 0;
		virtual void SetBoundingSpaceRect(Position start, Position end) = 0;
		virtual void SetOnEdge(unsigned int i, bool v) = 0;
		virtual void SetIsStrictlySpatial(bool v) = 0;
		virtual bool IsToroidal() const = 0;
	};

/**********************************************************************/
/* Spatial Network                                                    */
/**********************************************************************/
	template <typename LinkType> 
	class SpatialNetwork : public Network<LinkType>, public AbstractSpatialNetwork
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		SpatialNetwork(ParamHandler & h = ParamHandler::GlobalParams) :
			Network<LinkType>::Network(h), isStrictlySpatial(true)
		{
			positions = std::vector<Position *>(this->size(), 0);
			onEdge = std::vector<bool>(this->size(), false);
			dim = h.getParam<unsigned int>("-dim");
			toroidalSpace = h.getParam<bool>("-toroidalSpace", 0);

			std::string structBuildName = h.getParam<std::string>("-StructBuilder", 0);
			assert(AbstractFactory<SpatialStructureBuilder>::
				Factories[structBuildName]);
			structBuilder = AbstractFactory<SpatialStructureBuilder>::
				Factories[structBuildName]->Create();
			assert(structBuilder);
			freeStruct = true;
		}
		// Special constructor
		SpatialNetwork(unsigned int nbCells, NetworkConstructStrat *_c, std::string _necn,
			SpatialStructureBuilder *_s, unsigned int _dim,
			bool freeConstr = false, bool _freeStruct = false) : Network<LinkType>::Network(nbCells, _c, freeConstr, _necn), 
			dim(_dim), positions(nbCells, 0), onEdge(nbCells, false), structBuilder(_s), freeStruct(_freeStruct), isStrictlySpatial(true)
		{
			assert(structBuilder);
		}
		// Constructor from stream
		SpatialNetwork(std::ifstream & stream) : structBuilder(0), freeStruct(false), isStrictlySpatial(true)
		{
			LoadFromStream(stream);
		}

		// Destructor
		virtual ~SpatialNetwork()
		{
			ClearSpatialData();
			if (freeStruct and structBuilder)
				delete structBuilder;
		}

		//===========================================================||
		// Network modification methods (create, clear, ...)         ||
		//===========================================================||
		// Builds the spatial structure and the network
		virtual bool BuildNetwork(ResultSaver saver = ResultSaver::NullSaver)
		{
			bool ok = true;
			// rebuild positions if needed
			if (positions.size() != this->netSize)
			{
				assert(this->netSize);
				ClearSpatialData();
				positions = std::vector<Position *>(this->netSize, 0);
				onEdge = std::vector<bool>(this->netSize, false);
			}
			ok &= BuildSpatialStructure();
			return ok and Network<LinkType>::BuildNetwork(saver);
		}

		// Builds the spatial structure according to the builder strategy
		virtual bool BuildSpatialStructure()
		{
			if (structBuilder)
			{
				ConstStructBuild *tmpBuild = dynamic_cast<ConstStructBuild*>(structBuilder);
				if ((not tmpBuild) or tmpBuild->NeedsToBeBuilt())
				{
					ClearSpatialData();
					return structBuilder->BuildSpatialStructure(*this);
				}
				else
					return true;
			}
			else
				return false;
		}

		// Clears the network and free its contents
		virtual void ClearNetwork()
		{
			Network<LinkType>::ClearNetwork();
		}

		// Delete the positions
		virtual void ClearSpatialData()
		{
			for (unsigned int i = 0 ; i < positions.size() ; ++i)
				if (positions[i])
				{
					delete positions[i];
					positions[i] = 0;
				}
			onEdge = std::vector<bool>(positions.size(), false);
		}

		// Initializes the network
		virtual void Initialize()
		{
			Network<LinkType>::Initialize();
			spatialMetrics.InitializeMetricsDefault();
		}

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the cell from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = true;
			// Standard Network
			ok &= Network<LinkType>::LoadFromStream(stream);
			// dimension
			stream >> dim;
			// Cell positions
			ClearSpatialData();
			positions.clear();
			unsigned int nbCells;
			stream >> nbCells;
			for (unsigned int i = 0 ; (i < nbCells) and ok ; ++i)
				positions.push_back(new Position(stream));
			// Cell edge statuses
			stream >> nbCells;
			onEdge.clear();
			bool boolTemp;
			for (unsigned int i = 0 ; (i < nbCells) and ok ; ++i)
			{
				stream >> boolTemp;
				onEdge.push_back(boolTemp);
			}
			// Construction strategy
			std::string tempName;
			if (structBuilder and freeStruct)
				delete structBuilder;
			stream >> tempName;
			structBuilder = AbstractFactory<SpatialStructureBuilder>::Factories[tempName]
							->CreateFromStream(stream);
			freeStruct = true;
			// Spatial metrics
			spatialMetrics.FreeAndClean();
			ok &= spatialMetrics.LoadFromStream(stream);
			// strictly spatial
			stream >> isStrictlySpatial;
			// toroidal space ?
			stream >> toroidalSpace;
			// Bounding rect
			bStart = Position(stream);
			bEnd = Position(stream);

			return ok and stream.good();
		}

		// Saves the cell to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = true;
			// Network
			ok &= Network<LinkType>::SaveToStream(stream);
			// dimension
			stream << dim << std::endl;
			// Cell positions
			stream << positions.size() << std::endl;
			Position tempPos(GetDim());
			for (unsigned int i = 0 ; i < positions.size() ; ++i)
			{
				if (positions[i])
					ok &= positions[i]->SaveToStream(stream);
				else
					tempPos.SaveToStream(stream);
			}
			// Cell edge statuses
			stream << onEdge.size() << std::endl;
			for (unsigned int i = 0 ; i < onEdge.size() ; ++i)
				stream << onEdge[i] << std::endl;
			// Spatial Structure building strategy
			stream << structBuilder->GetClassName() << std::endl;
			ok &= structBuilder->SaveToStream(stream);
			// Spatial Metrics
			ok &= spatialMetrics.SaveToStream(stream);
			// strictly spatial
			stream << isStrictlySpatial << std::endl;
			// toroidal space ?
			stream << toroidalSpace << std::endl;
			// Bounding rect
			ok &= bStart.SaveToStream(stream);
			ok &= bEnd.SaveToStream(stream);

			return ok and stream.good();
		}

		//===========================================================||
		// Special metrics handling functions                        ||
		//===========================================================||
		// Compute and Save desired metrics
		virtual bool ComputeAndSaveMetrics(ResultSaver saver)
		{
			bool ok = true;

			ok &= Network<LinkType>::ComputeAndSaveMetrics(saver);
			ok &= spatialMetrics.ComputeMetricsDefault(*this);
			ok &= spatialMetrics.SaveMetricsDefault(saver);
			return ok;
		}

		// Add a metric
		virtual bool AddMetric(Metric *_m, bool _f = false) 
		{
			if (Network<LinkType>::AddMetric(_m, _f))
				return true;
			else
				return spatialMetrics.AddMetricAndDependencies(_m, _f, this); 
		}

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler() 
		{
			ParamHandler params;
			params += Network<LinkType>::BuildModelParamHandler();
			if (structBuilder)
				params += structBuilder->BuildModelParamHandler();
			params += spatialMetrics.BuildModelParamHandler();
			return params;
		}

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		// Returns a reference to the pointer to the position of the ith cell
		virtual Position* & Pos(unsigned int i)
		{
			static Position* voidPos;
			if (i < positions.size())
				return positions[i];
			else
				return voidPos;
		}
		// Const version of Pos
		virtual Position* const & Pos(unsigned int i) const
		{
			static Position* voidPos;
			if (i < positions.size())
				return positions[i];
			else
				return voidPos;
		}
		// Returns the distance between two point according to 
		// the chosen spatial structure (euclidean or toroidal)
		virtual double GetDistanceBetween(unsigned int i, unsigned int j) const
		{
			assert(i < positions.size() and j < positions.size());
			Position & a = *positions[i];
			Position & b = *positions[j];
			Position euclidDist = (a - b).abs();
			if (toroidalSpace)
			{
				Position toroidDist(a.size(), 0);
				for (unsigned d = 0 ; d < toroidDist.size() ; ++d)
				{
					toroidDist[d] = std::min(euclidDist[d], 
						std::min((a[d] - bStart[d]) + (bEnd[d] - b[d]), 
							(b[d] - bStart[d]) + (bEnd[d] - a[d])));
				}
				return toroidDist.norm();
			}
			else
				return euclidDist.norm();
		}

		// Returns the dimension number
		virtual unsigned int GetDim() const { return dim; }
		// Returns the dimension number
		virtual void SetDim(unsigned int _d) { dim = _d; }
		// Returns all metrics and submetrics
		virtual std::vector<Metric *> GetAllMetrics() const
		{
			std::vector<Metric *> mTot, mTemp;
			mTemp = Network<LinkType>::GetAllMetrics();
			for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
				mTot.push_back(mTemp[i]);
			for (unsigned int i = 0 ; i < spatialMetrics.GetMetricsRaw().size() ; ++i)
				mTot.push_back(spatialMetrics.GetMetricsRaw()[i]);
			return mTot;
		}

		virtual bool IsNodeOnEdge(unsigned int i) const 
		{
			return (i < onEdge.size() and not toroidalSpace) ? onEdge[i] : false;
		}

		virtual void SetBoundingSpaceRect(Position start, Position end)
		{
			bStart = start;
			bEnd = end;
		}

		virtual void SetOnEdge(unsigned int i, bool v)
		{
			if (i < onEdge.size())
				onEdge[i] = v;
		}
		virtual bool IsSpatial() const { return isStrictlySpatial; }
		virtual void SetIsStrictlySpatial(bool v) { isStrictlySpatial = v;}
		virtual bool IsToroidal() const { return toroidalSpace; };

	protected:
		// Dimension of the embeding space
		unsigned int dim;
		// Positions of each cells
		std::vector<Position *> positions;
		// Edge status of nodes
		std::vector<bool> onEdge;
		// Spatial structure builder strategy
		SpatialStructureBuilder *structBuilder;
		// Delete spatial structure building strategy ?
		bool freeStruct;
		// Tells wether mean path length scales as a polynom of the number of nodes
		bool isStrictlySpatial;
		// Bounding rect encompassing the network (used to simulate toroids)
		Position bStart;
		Position bEnd;
		// Defines if the space is to be considered toroidal (periodic boundary conditions)
		bool toroidalSpace;

		SortedMetrics<AbstractSpatialNetwork> spatialMetrics;
	};
}

#endif
