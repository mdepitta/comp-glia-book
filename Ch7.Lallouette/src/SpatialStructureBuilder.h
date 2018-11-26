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

#ifndef SPATIALSTRUCTUREBUILDER_H
#define SPATIALSTRUCTUREBUILDER_H

#include "Savable.h"
#include "ParamHandler.h"
#include "utility.h"

#include <math.h>
#include <assert.h>

namespace AstroModel
{
	// Forward declarations
	class AbstractSpatialNetwork;

	// Abstract base class for structures that don't need to get rebuilt
	class ConstStructBuild
	{
	public:
		virtual bool NeedsToBeBuilt() const = 0;
	};

/**********************************************************************/
/* Abstract class                                                     */
/**********************************************************************/
	class SpatialStructureBuilder : public SaveAndLoadFromStream
	{
	public:
		// Returns the class name
		virtual std::string GetClassName() const = 0;
		//===========================================================||
		// Standard Spatial Structure Builder Strategy method        ||
		//===========================================================||
		virtual bool BuildSpatialStructure(AbstractSpatialNetwork &
			network) const = 0;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler() = 0;
	};

/**********************************************************************/
/* Standard Regular                                                   */
/**********************************************************************/
	class StandardStructureBuilder : public SpatialStructureBuilder
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		StandardStructureBuilder(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		StandardStructureBuilder(double _a, double _var, double _md);
		// Constructor from stream
		StandardStructureBuilder(std::ifstream & stream);

		//===========================================================||
		// Standard Spatial Structure Builder Strategy method        ||
		//===========================================================||
		virtual bool BuildSpatialStructure(AbstractSpatialNetwork & 
			network) const;

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
		double a;       // Initial lattice grid size
		double var;     // Variance of jitter
		double minDist; // Minimum distance between two cells

		//===========================================================||
		// Recursive lattice filling method                          ||
		//===========================================================||
		virtual unsigned int recursiveFill(AbstractSpatialNetwork & network, 
			Position currPos, unsigned int dimLvl, unsigned int nbEachDim,	
			unsigned int currInd, unsigned int dim, bool oE) const;
	};

/**********************************************************************/
/* Copy structure builder                                             */
/**********************************************************************/
	class CopyStructureBuilder : public SpatialStructureBuilder
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		CopyStructureBuilder(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		CopyStructureBuilder(const AbstractSpatialNetwork & net);
		// Constructor from stream
		CopyStructureBuilder(std::ifstream & stream);

		//===========================================================||
		// Standard Spatial Structure Builder Strategy method        ||
		//===========================================================||
		virtual bool BuildSpatialStructure(AbstractSpatialNetwork & 
			network) const;

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
		// Accesors                                                  ||
		//===========================================================||
		void ChangeRefNet(const AbstractSpatialNetwork & net);

	protected:
		const AbstractSpatialNetwork* refNet;
	};

/**********************************************************************/
/* Uniform random structure builder                                   */
/**********************************************************************/
	class RandomStructureBuilder : public SpatialStructureBuilder
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		RandomStructureBuilder(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		RandomStructureBuilder(double _md);
		// Constructor from stream
		RandomStructureBuilder(std::ifstream & stream);

		//===========================================================||
		// Standard Spatial Structure Builder Strategy method        ||
		//===========================================================||
		virtual bool BuildSpatialStructure(AbstractSpatialNetwork & 
			network) const;

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
		double meanDist; // Cube side size
	};

/**********************************************************************/
/* Load From File Builder                                             */
/**********************************************************************/
	class LoadFromFileStructureBuilder : public SpatialStructureBuilder,
		public ConstStructBuild
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		LoadFromFileStructureBuilder(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		LoadFromFileStructureBuilder(std::string _p);
		// Constructor from stream
		LoadFromFileStructureBuilder(std::ifstream & stream);

		//===========================================================||
		// Standard Spatial Structure Builder Strategy method        ||
		//===========================================================||
		virtual bool BuildSpatialStructure(AbstractSpatialNetwork & 
			network) const;
		// Returns whether the file has to be reloaded
		virtual bool NeedsToBeBuilt() const;

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
		std::string newPath;
	};

/**********************************************************************/
/* Tree-like astrocytes structure builder                             */
/**********************************************************************/
	class StarLikeAstrocyteStruct : public SpatialStructureBuilder
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default constructor
		StarLikeAstrocyteStruct(ParamHandler & h = 
			ParamHandler::GlobalParams);
		// Special constructor
		StarLikeAstrocyteStruct(double icd, double icod, unsigned int bl, 
			unsigned int nbB, unsigned int nbeb, bool ms);
		// Constructor from stream
		StarLikeAstrocyteStruct(std::ifstream & stream);

		//===========================================================||
		// Standard Spatial Structure Builder Strategy method        ||
		//===========================================================||
		virtual bool BuildSpatialStructure(AbstractSpatialNetwork & 
			network) const;

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
		double interCellDist; // Distance between two cell body centers
		double interCompartDist; // Distance between two subcompartments
		unsigned int branchLength; // Number of subcompartments in a single branch (both sides)
		unsigned int nbBranches; // Number of branches for each cell body
		unsigned int nbEndBall; // Number of compartment in bigger process end
		bool makeSquare;

		unsigned int makeGrid(AbstractSpatialNetwork & network, Position botLeft, 
			unsigned int currInd, unsigned int nbRows, unsigned int nbCols) const;
		unsigned int makeStar(AbstractSpatialNetwork & network, unsigned int astroNum,
			Position center, unsigned int currInd) const;
		unsigned int makeBranch(AbstractSpatialNetwork & network, 
			unsigned int astroNum, unsigned int processNum, Position center,
			Position direction, unsigned int currInd) const;
	};
}

#endif
