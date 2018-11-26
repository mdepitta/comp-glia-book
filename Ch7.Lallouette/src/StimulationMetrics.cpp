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

#include "StimulationMetrics.h"
#include "StimulationStrat.h"
#include "MetricNames.h"

using namespace AstroModel;
using namespace std;

string StimulatedCells::ClassName(STIMULATED_CELLS_METRIC);
string MeanPathStimMetric::ClassName(MEAN_PATH_STIM_METRIC);

//********************************************************************//
//***************** S T I M U L A T E D   C E L L S ******************//
//********************************************************************//

bool StimulatedCells::stimFileHasBeenSaved = false;

//**********************************************************************
// Default Constructor
//**********************************************************************
StimulatedCells::StimulatedCells() : currOrigin("Unknown")
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
StimulatedCells::StimulatedCells(std::ifstream & stream) 
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Compute the activated cells
//**********************************************************************
bool StimulatedCells::ComputeMetric(const StimulationStrat & stimStrat)
{
	double t = stimStrat.GetTime();
	currOrigin = stimStrat.GetClassName();

	if (stimulated.size() != stimStrat.GetSize())
		stimulated = vector<pair<bool, double> >(stimStrat.GetSize(), 
			make_pair(false, t));

	for (unsigned int i = 0 ; i < stimStrat.GetSize() ; ++i)
	{
		if ((not stimulated[i].first) and stimStrat.IsStimulated(i))
		{
			stimulated[i].first = true;
			stimulated[i].second = t;
			stimCells.insert(i);
		}
		else if (stimulated[i].first and (not stimStrat.IsStimulated(i)))
		{
			stimulated[i].first = false;
			stimulations.push_back( StimulationInd(stimStrat.GetClassName(),
				i, stimulated[i].second, t, 0, currOrigin));
		}
	}

	lastT = t;
	return true;
}

//**********************************************************************
// Save the activated cells
//**********************************************************************
bool StimulatedCells::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	// End Current started stimulations
	for (unsigned int i = 0 ; i < stimulated.size() ; ++i)
		if (stimulated[i].first)
			stimulations.push_back(
				StimulationInd("", i, stimulated[i].second, lastT, 0, currOrigin));

	std::string stimulationsName("Stimulations");

	// Saving
	if (saver.isSaving(stimulationsName) and not stimulations.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(stimulationsName, saver.getCurrFile());
		
		if (not stimFileHasBeenSaved)
			stream << "CellNb\tTStart\tTEnd" << endl; 

		for (unsigned int i = 0 ; i < stimulations.size() ; ++i)
			 stream 
			 	<< stimulations[i].ind   << "\t"
			 	<< stimulations[i].start << "\t"
			 	<< stimulations[i].end   << endl; 

		allSaved &= stream.good();
		stimFileHasBeenSaved = true;
	}
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void StimulatedCells::Initialize()
{
	Metric::Initialize();
	lastT = 0;
	stimulated.clear();
	stimulations.clear();
	stimCells.clear();
	stimFileHasBeenSaved = false;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * StimulatedCells::BuildCopy() const
{
	return new StimulatedCells(); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool StimulatedCells::LoadFromStream(std::ifstream & stream)
{
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool StimulatedCells::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//********************************************************************//
//****** M E A N   P A T H   S T I M U L A T I O N   M E T R I C *****//
//********************************************************************//

//**********************************************************************
// Compute the activated cells
//**********************************************************************
bool MeanPathStimMetric::ComputeMetric(const StimulationStrat & stimStrat)
{
	bool computeOk = true; 
	computeOk &= StimulatedCells::ComputeMetric(stimStrat);
	const MeanPathStimStrat *mps = 
		dynamic_cast<const MeanPathStimStrat *>(&stimStrat);
	assert(mps);
	
	effMeanPathLength = mps->GetCurrEffMeanPathLength();

	return computeOk;
}

//**********************************************************************
// Save the activated cells
//**********************************************************************
bool MeanPathStimMetric::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	allSaved &= StimulatedCells::SaveMetric(saver);

	std::string mpsName("MeanPathStimInfos");
	if (saver.isSaving(mpsName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(mpsName, saver.getCurrFile());
		
		stream << "EffectiveMeanPathLength" << endl;
		stream << effMeanPathLength << std::endl;

		allSaved &= stream.good();
	}
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void MeanPathStimMetric::Initialize()
{
	StimulatedCells::Initialize();
	effMeanPathLength = -1;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * MeanPathStimMetric::BuildCopy() const
{
	return new MeanPathStimMetric(); 
}

