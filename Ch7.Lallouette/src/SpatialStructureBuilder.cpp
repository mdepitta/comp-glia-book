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

#include "SpatialStructureBuilder.h"

#include "SpatialNetwork.h"

#include <cmath>

using namespace AstroModel;
using namespace std;

string StandardStructureBuilder::ClassName(
	"StandardSpatialStructureBuilder");
string CopyStructureBuilder::ClassName(
	"CopySpatialStructureBuilder");
string RandomStructureBuilder::ClassName(
	"RandomSpatialStructureBuilder");
string LoadFromFileStructureBuilder::ClassName(
	"LoadFromFileStructureBuilder");
string StarLikeAstrocyteStruct::ClassName(
	"StarLikeAstrocyteStructureBuilder");

//********************************************************************//
//******** S T A N D A R D   S T R U C T U R E   B U I L D E R *******//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
StandardStructureBuilder::StandardStructureBuilder(ParamHandler & h)
{
	a = h.getParam<double>("-StdStructBuild", 0);
	var = h.getParam<double>("-StdStructBuild", 1);
	minDist = h.getParam<double>("-StdStructBuild", 2);
}

//**********************************************************************
// Special constructor
//**********************************************************************
StandardStructureBuilder::StandardStructureBuilder(double _a,
	double _var, double _md) : 
a(_a), var(_var), minDist(_md)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
StandardStructureBuilder::StandardStructureBuilder(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler StandardStructureBuilder::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "gridSize"   , a;
	params <= "gridJitter" , var;
	params <= "gridMinDist", minDist;
	return params;
}

//**********************************************************************
// Builds the spatial structure
//**********************************************************************
bool StandardStructureBuilder::BuildSpatialStructure(AbstractSpatialNetwork & network) const
{
	unsigned int dim = network.GetDim();
	int nbEachDim = ceil(pow((double)network.GetNeededSize(), 1.0 / (double)dim));

	Position tempPos(dim);
	for (unsigned int dimNb = 0 ; dimNb < dim ; ++dimNb)
		tempPos[dimNb] = 0;

	this->recursiveFill(network, tempPos, 0, nbEachDim, 0, dim, false);
	network.SetBoundingSpaceRect(Position(dim, -(a - minDist) / 2.0),
		Position(dim, ((double)nbEachDim - 1.0) * a + (a - minDist) / 2.0));
	return true;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool StandardStructureBuilder::LoadFromStream(std::ifstream & stream)
{
	stream >> a;
	stream >> var;
	stream >> minDist;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool StandardStructureBuilder::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< a       << std::endl
		<< var     << std::endl
		<< minDist << std::endl;
	return stream.good();
}

//**********************************************************************
// Recursive lattice filling method
//**********************************************************************
unsigned int StandardStructureBuilder::recursiveFill(
	AbstractSpatialNetwork & network, Position currPos, 
	unsigned int dimLvl, unsigned int nbEachDim,	
	unsigned int currInd, unsigned int dim, bool oE) const
{
	assert(a - minDist > 0);
	assert(dim > 0);
	Position tempPos(currPos);
	bool localOECond = false;
	for (unsigned int j = 0 ; j < nbEachDim ; ++j)
	{
		localOECond = (j == 0) or (j == (nbEachDim - 1));
		currPos[dimLvl] = (double)j * a;
		if (dimLvl < dim - 1)
			currInd = recursiveFill(network, currPos, dimLvl + 1, nbEachDim, currInd, dim, oE or localOECond);
		else
		{
			do
			{
				for (unsigned int d = 0 ; d < dim ; ++d)
					tempPos[d] = currPos[d] + GaussianRand(0, var);
			} while (tempPos.DistanceTo(currPos) > ((a - minDist) / 2.0));
			network.Pos(currInd) = new Position(tempPos);
			network.SetOnEdge(currInd, oE or localOECond);
			++currInd;
		}
	}
	return currInd;
}

//********************************************************************//
//************ C O P Y   S T R U C T U R E   B U I L D E R ***********//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
CopyStructureBuilder::CopyStructureBuilder(ParamHandler &) : refNet(0)
{
}

//**********************************************************************
// Special constructor
//**********************************************************************
CopyStructureBuilder::CopyStructureBuilder(
	const AbstractSpatialNetwork & net) : refNet(&net)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
CopyStructureBuilder::CopyStructureBuilder(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler CopyStructureBuilder::BuildModelParamHandler() 
{
	ParamHandler params;
	return params;
}

//**********************************************************************
// Builds the spatial structure
//**********************************************************************
bool CopyStructureBuilder::BuildSpatialStructure(
	AbstractSpatialNetwork & network) const
{
	assert(network.size() == refNet->size());
	assert(network.GetDim() == refNet->GetDim());
	for (unsigned int i = 0 ; i < network.size() ; ++i)
		network.Pos(i) = new Position(*refNet->Pos(i));
	return true;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool CopyStructureBuilder::LoadFromStream(std::ifstream & stream)
{
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool CopyStructureBuilder::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//********************************************************************//
//********* R A N D O M   S T R U C T U R E   B U I L D E R **********//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
RandomStructureBuilder::RandomStructureBuilder(ParamHandler & h)
{
	meanDist = h.getParam<double>("-RandStructBuild", 0);
}

//**********************************************************************
// Special constructor
//**********************************************************************
RandomStructureBuilder::RandomStructureBuilder(double _md) :
	meanDist(_md)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
RandomStructureBuilder::RandomStructureBuilder(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler RandomStructureBuilder::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "RandomStructMeanDist", meanDist;
	return params;
}

//**********************************************************************
// Builds the spatial structure
//**********************************************************************
bool RandomStructureBuilder::BuildSpatialStructure(AbstractSpatialNetwork & network) const
{
	unsigned int dim = network.GetDim();

	double N = pow(
		exp(lgamma((double) dim / 2.0 + 1.0) + lgamma(1.0 + 1.0 / (double) dim) / (double) dim) 
		/ (meanDist * sqrt(M_PI)), (double) dim);

	double cubeSide = pow((double) network.GetNeededSize() / N, 1.0 / (double) dim);
	std::cout << cubeSide << std::endl;

	Position tempPos(dim, 0);
	for (unsigned int i = 0 ; i < network.GetNeededSize() ; ++i)
	{
		for (unsigned int d = 0 ; d < tempPos.size() ; ++d)
			tempPos[d] = UnifRand() * cubeSide;
		network.Pos(i) = new Position(tempPos);
	}
	Position bStart(dim, 0);
	Position bEnd(dim, cubeSide);
	network.SetBoundingSpaceRect(bStart, bEnd);

	return true;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool RandomStructureBuilder::LoadFromStream(std::ifstream & stream)
{
	stream >> meanDist;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool RandomStructureBuilder::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< meanDist << std::endl;
	return stream.good();
}


//********************************************************************//
//** L O A D   F R O M   F I L E   S T R U C T U R E   B U I L D E R *//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
LoadFromFileStructureBuilder::LoadFromFileStructureBuilder(ParamHandler & h)
{
	newPath = h.getParam<std::string>("-LoadFrmFileStructBuild", 0);
}

//**********************************************************************
// Special constructor
//**********************************************************************
LoadFromFileStructureBuilder::LoadFromFileStructureBuilder(std::string _p) :
	newPath(_p)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
LoadFromFileStructureBuilder::LoadFromFileStructureBuilder(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler LoadFromFileStructureBuilder::BuildModelParamHandler() 
{
	ParamHandler params;
	params <= "LoadFromFileStructureBuilderPath", newPath;
	return params;
}

//**********************************************************************
// Builds the spatial structure
//**********************************************************************
bool LoadFromFileStructureBuilder::BuildSpatialStructure(AbstractSpatialNetwork & network) const
{
	path = newPath;
	std::ifstream stream(path.c_str());
	assert(stream.good());
	unsigned int dim;
	stream >> dim;
	network.SetDim(dim);
	unsigned int N;
	stream >> N;
	assert(N == network.GetNeededSize());

	Position minPos(dim, DEFAULT_MAX_VAL);
	Position maxPos(dim, -DEFAULT_MAX_VAL);
	Position tempPos(dim);
	double tmpVal;
	for (unsigned int i = 0 ; i < N ; ++i)
	{
		for (unsigned int dimNb = 0 ; dimNb < dim ; ++dimNb)
		{
			stream >> tmpVal;
			tempPos[dimNb] = tmpVal;
			minPos[dimNb] = (minPos[dimNb] > tmpVal) ? tmpVal : minPos[dimNb];
			maxPos[dimNb] = (maxPos[dimNb] < tmpVal) ? tmpVal : maxPos[dimNb];
		}
		network.Pos(i) = new Position(tempPos);
	}

	network.SetBoundingSpaceRect(minPos, maxPos);
	return true;
}

//**********************************************************************
// Returns whether the file has to be reloaded
//**********************************************************************
bool LoadFromFileStructureBuilder::NeedsToBeBuilt() const
{
	return path != newPath;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool LoadFromFileStructureBuilder::LoadFromStream(std::ifstream & stream)
{
	stream >> newPath;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool LoadFromFileStructureBuilder::SaveToStream(std::ofstream & stream) const
{
	stream << newPath << std::endl;
	return stream.good();
}

//********************************************************************//
//******* S T A R   L I K E   A S T R O C Y T E   B U I L D E R ******//
//********************************************************************//

//**********************************************************************
// Default constructor
//**********************************************************************
StarLikeAstrocyteStruct::StarLikeAstrocyteStruct(ParamHandler & h)
{
	interCellDist = h.getParam<double>("-StarLikeBuild", 0);
	interCompartDist = h.getParam<double>("-StarLikeBuild", 1);
	branchLength = h.getParam<unsigned int>("-StarLikeBuild", 2);
	nbBranches = h.getParam<unsigned int>("-StarLikeBuild", 3);
	nbEndBall = h.getParam<unsigned int>("-StarLikeBuild", 4);
	makeSquare = h.getParam<bool>("-StarLikeMakeSquare", 0);
}

//**********************************************************************
// Special constructor
//**********************************************************************
StarLikeAstrocyteStruct::StarLikeAstrocyteStruct(double icd, double icod, 
	unsigned int bl, unsigned int nbB, unsigned int nbeb, bool ms) : 
	interCellDist(icd), interCompartDist(icod),	branchLength(bl), 
	nbBranches(nbB), nbEndBall(nbeb), makeSquare(ms)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
StarLikeAstrocyteStruct::StarLikeAstrocyteStruct(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Builds the param handler
//**********************************************************************
ParamHandler StarLikeAstrocyteStruct::BuildModelParamHandler() 
{
	ParamHandler params;

	params <= "StarLikeInterCellDist", interCellDist;
	params <= "StarLikeInterCompartDist", interCompartDist;
	params <= "StarLikeBranchLength", branchLength;
	params <= "StarLikeNbBranches", nbBranches;
	params <= "StarLikeNbEndBall", nbEndBall;
	params <= "StarLikeMakeSquare", makeSquare;

	return params;
}

//**********************************************************************
// Builds the spatial structure
//**********************************************************************
bool StarLikeAstrocyteStruct::BuildSpatialStructure(AbstractSpatialNetwork & network) const
{
	unsigned int dim = network.GetDim();
	assert(dim == 2);

	// Create the cell bodies
	unsigned int nbCells = network.GetNeededSize() / (branchLength * nbBranches);
	assert(nbCells * (branchLength * nbBranches + 1) == network.GetNeededSize());

	unsigned int nbRows = floor(sqrt((double)nbCells));
	unsigned int nbCols = ceil((double) nbCells / (double) nbRows);
	unsigned int currInd = 0;
	currInd = makeGrid(network, Position(dim, 0), currInd, nbRows, nbCols);
	assert(currInd == network.GetNeededSize());

	network.SetBoundingSpaceRect(Position(dim, 0), Position(dim, nbCols * interCellDist));

	return true;
}

//**********************************************************************
// Loads the strategy from a stream
//**********************************************************************
bool StarLikeAstrocyteStruct::LoadFromStream(std::ifstream & stream)
{
	stream >> interCellDist;
	stream >> interCompartDist;
	stream >> branchLength;
	stream >> nbBranches;
	stream >> makeSquare;
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the strategy to a stream
//**********************************************************************
bool StarLikeAstrocyteStruct::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< interCellDist << std::endl
		<< interCompartDist << std::endl
		<< branchLength << std::endl
		<< nbBranches << std::endl
		<< makeSquare << std::endl;
	return stream.good();
}


//**********************************************************************
// Build a grid of stars
//**********************************************************************
unsigned int StarLikeAstrocyteStruct::makeGrid(AbstractSpatialNetwork & network,
	Position botLeft, unsigned int currInd, unsigned int nbRows, unsigned int nbCols) const
{
	Position cent = botLeft;
	unsigned int astroNum = 0;
	for (unsigned int i = 0 ; i < nbRows and currInd < network.GetNeededSize() ; ++i)
	{
		for (unsigned int j = 0 ; j < nbCols and currInd < network.GetNeededSize() ; ++j) 
		{
			currInd = makeStar(network, astroNum, cent, currInd);
			cent[0] += interCellDist;
			++astroNum;
		}
		cent[1] += interCellDist;
	}
	return currInd;
}

//**********************************************************************
// Builds a full star
//**********************************************************************
unsigned int StarLikeAstrocyteStruct::makeStar(AbstractSpatialNetwork & network, 
	unsigned int astroNum, Position center, unsigned int currInd) const
{
	Position dir(network.GetDim(), 0);
	double angle = 0;
	// Add center
	network.Pos(currInd) = new Position(center);
	network.SetNodeTag(currInd, (0x7F0000 & (astroNum << 16)) + 
								  (0xFF00 & 0) + 
									(0xFF & 0));
	++currInd;
	// Add branches
	for (unsigned int i = 0 ; i < nbBranches ; ++i)
	{
		angle = 2.0 * 3.14159265 * ((double) i) / ((double) nbBranches);
		dir[0] = cos(angle);
		dir[1] = sin(angle);
		currInd = makeBranch(network, astroNum, i, center, dir, currInd);
	}
	return currInd;
}

//**********************************************************************
// Builds a star branch
//**********************************************************************
unsigned int StarLikeAstrocyteStruct::makeBranch(AbstractSpatialNetwork & network, 
	unsigned int astroNum, unsigned int processNum, Position center, 
	Position direction, unsigned int currInd) const
{
	double scalPos = interCompartDist;
	bool isOk = true;
	bool wasOk;
	for (unsigned int i = 0 ; i < branchLength ; ++i)
	{
		if (isOk)
		{
			network.Pos(currInd) = new Position(center + (direction * scalPos));
			network.SetOnEdge(currInd, i == (branchLength-1));
			scalPos += interCompartDist;
			// Tag identification code:
			// subelement number S on 8 bits (starts at 0)
			// process number P on on 8 bits (starts at 1)
			// astrocyte number A on 7 bits (starts at 1)
			// border labeling on 1 bit
			// 0xBA A P P P S S S
			network.SetNodeTag(currInd, (0x7F0000 & ((astroNum + 1) << 16)) + 
										  (0xFF00 & ((processNum + 1) << 8)) + 
											(0xFF & i));
		}
		else
		{
			network.Pos(currInd) = new Position(network.GetDim(),-interCellDist);
			network.SetOnEdge(currInd, true);
		}

		++currInd;
		wasOk = isOk;
		isOk &= !makeSquare or ((fabs((direction*scalPos)[0]) <= interCellDist / 2.0) and 
			(fabs((direction*scalPos)[1]) <= interCellDist / 2.0));
		// If we reached the end of the process, make the end ball
		if (wasOk and not isOk)
		{
			// Tag the last node as the cell border
			network.SetNodeTag(currInd, network.GetNodeTag(currInd) + 0x800000);
			Position tmpPos;
			for (unsigned int j = 0 ; (j < nbEndBall) and (i < (branchLength-1)) ; ++j,++i)
			{
				tmpPos = center + (direction * (scalPos - interCompartDist));
				tmpPos[0] += UnifRand() * interCompartDist;
				tmpPos[1] += UnifRand() * interCompartDist;
				network.Pos(currInd) = new Position(tmpPos);
				network.SetOnEdge(currInd, true);
				++currInd;
			}
		}
	}
	return currInd;
}

