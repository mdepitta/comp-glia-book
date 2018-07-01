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

#include "SimulationMetrics.h"
#include "GridSearchSimulation.h"
#include "ChISimulationManager.h"
#include "MetricNames.h"
#include "PropagationMetrics.h"
#include "StimulationMetrics.h"

#include "libalglib/statistics.h"

#include <algorithm>

using namespace AstroModel;
using namespace std;

//********************************************************************//
//************ C H I   S I M U L A T I O N   M E T R I C *************//
//********************************************************************//

string PropagationDistanceStat::ClassName(PROPAGATION_DISTANCE_STAT_METRIC);
string NetworkTopologyStats::ClassName(NETWORK_TOPOLOGY_AND_SPATIAL_STAT_METRIC);
string ScalarStatisticsMetric::ClassName(SCALAR_STATISTICS_METRIC);
string GridSearchMetric::ClassName(GRID_SEARCH_METRIC);

//********************************************************************//
//******** P R O P A G A T I O N   D I S T A N C E   S T A T *********//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(PROPAGATION_DISTANCE_STAT_METRIC, PROPAGATION_DISTANCE_METRIC)
ADD_METRIC_DEPENDENCY(PROPAGATION_DISTANCE_STAT_METRIC, CLUSTERING_COEFFICIENT_DISTRIBUTION_NETWORK_METRIC)

//**********************************************************************
// Default Constructor
//**********************************************************************
PropagationDistanceStat::PropagationDistanceStat()
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
PropagationDistanceStat::PropagationDistanceStat(ifstream & stream)
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Default Destructor
//**********************************************************************
PropagationDistanceStat::~PropagationDistanceStat()
{
}

//**********************************************************************
// Compute the metric
//**********************************************************************
bool PropagationDistanceStat::ComputeMetric(const AbstractRepeatSimulation & manag)
{
	// For each wave
	double maxNb = 0;
	double maxSpeed = 0;

	ClusteringCoeffs *netClustCoeff = 
		GetSpecificMetric<Metric, ClusteringCoeffs>(manag.GetAllMetrics());
	PropagationDistance *propDist = 
		GetSpecificMetric<Metric, PropagationDistance>(manag.GetAllMetrics());

	nbWaves.push_back(propDist->size());
	for (unsigned int w = 0 ; w < propDist->size() ; ++w)
	{
		// nb Cells in wave by clustering coeff of initiator
		assert(netClustCoeff);
		nbCellsByClustCoeff.push_back(make_pair((*netClustCoeff)[(*propDist)[w].initiator], 
			(*propDist)[w].involvedCells.size()));

		double startTime = (*propDist)[w].propagation[0].time;
		for (unsigned int p = 0 ; p < (*propDist)[w].propagation.size() ; ++p)
		{
			maxNb = std::max(maxNb, (double)(*propDist)[w].propagation[p].cellNb);
			if (p != 0)
				maxSpeed = std::max(maxSpeed, (*propDist)[w].propagation[p].spatialDist / ((*propDist)[w].propagation[p].time - startTime));
			
			spatialDistByCellNb[(*propDist)[w].propagation[p].cellNb].push_back((*propDist)[w].propagation[p].spatialDist);
			spatialDistByCellDist[(*propDist)[w].propagation[p].cellDist].push_back((*propDist)[w].propagation[p].spatialDist);

		}
	}
	maxNbCells.push_back(maxNb);
	maxWaveSpeed.push_back(maxSpeed);

	return true;
}


//**********************************************************************
// Save the metric
//**********************************************************************
bool PropagationDistanceStat::SaveMetric(ResultSaver saver) const
{
	static const string spatialDistByNbCellName("SpatialDistByNbCell");
	static const string spatialDistByCellDistName("SpatialDistByCellDist");
	static const string maxNbCellDistributionName("MaxNbCellDistribution");
	static const string NbCellByClustCoeffName("NbCellInWaveByClustCoeff");
	static const string maxWaveSpeedDistributionName("MaxWaveSpeedDistribution");
	static const string nbWavesDistribName("NbWavesDistribution");

	bool allSaved = true;

	if (saver.isSaving(spatialDistByNbCellName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(spatialDistByNbCellName, saver.getCurrFile());

		stream << "CellNb\tSpatialDistMean\tSpatialDistStdDev" << endl;
		for (SDistByCellMapType::const_iterator it = spatialDistByCellNb.begin() ; it != spatialDistByCellNb.end() ; ++it)
		{
			stream
				<< it->first << "\t"
				<< ComputeMean(it->second) << "\t"
				<< ComputeStdDev(it->second) << endl;
		}

		allSaved &= stream.good();
	}
	if (saver.isSaving(spatialDistByCellDistName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(spatialDistByCellDistName, saver.getCurrFile());

		stream << "CellDist\tSpatialDistMean\tSpatialDistStdDev" << endl;
		for (SDistByCDistT::const_iterator it = spatialDistByCellDist.begin() ; it != spatialDistByCellDist.end() ; ++it)
		{
			stream
				<< it->first << "\t"
				<< ComputeMean(it->second) << "\t"
				<< ComputeStdDev(it->second) << endl;
		}

		allSaved &= stream.good();
	}
	if (saver.isSaving(maxNbCellDistributionName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(maxNbCellDistributionName, saver.getCurrFile());

		for (unsigned int i = 0 ; i < maxNbCells.size() ; ++i)
			stream << maxNbCells[i] << endl;

		allSaved &= stream.good();
	}
	if (saver.isSaving(NbCellByClustCoeffName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(NbCellByClustCoeffName, saver.getCurrFile());

		stream << "ClusteringCoefficient\tNumberOfInvolvedCells" << endl;
		for (unsigned int i = 0 ; i < nbCellsByClustCoeff.size() ; ++i)
			stream 
				<< nbCellsByClustCoeff[i].first << "\t"
				<< nbCellsByClustCoeff[i].second << endl;

		allSaved &= stream.good();
	}
	if (saver.isSaving(maxWaveSpeedDistributionName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(maxWaveSpeedDistributionName, saver.getCurrFile());

		for (unsigned int i = 0 ; i < maxWaveSpeed.size() ; ++i)
			stream << maxWaveSpeed[i] << endl;

		allSaved &= stream.good();
	}
	if (saver.isSaving(nbWavesDistribName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(nbWavesDistribName, saver.getCurrFile());

		for (unsigned int i = 0 ; i < nbWaves.size() ; ++i)
			stream << nbWaves[i] << endl;

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void PropagationDistanceStat::Initialize()
{
	Metric::Initialize();
	spatialDistByCellNb.clear();
	spatialDistByCellDist.clear();
	nbCellsByClustCoeff.clear();
	maxNbCells.clear();
	maxWaveSpeed.clear();
	nbWaves.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * PropagationDistanceStat::BuildCopy() const
{
	return new PropagationDistanceStat(); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool PropagationDistanceStat::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool PropagationDistanceStat::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//********************************************************************//
//****** N E T W O R K   T O P O L O G Y   S T A T   M E T R I C *****//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(NETWORK_TOPOLOGY_AND_SPATIAL_STAT_METRIC, SPATIAL_POSITION_NETWORK_METRIC)


//**********************************************************************
// Default Constructor
//**********************************************************************
NetworkTopologyStats::NetworkTopologyStats()
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
NetworkTopologyStats::NetworkTopologyStats(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Default Destructor
//**********************************************************************
NetworkTopologyStats::~NetworkTopologyStats()
{
}

//**********************************************************************
// Compute the metric
//**********************************************************************
bool NetworkTopologyStats::ComputeMetric(const AbstractRepeatSimulation & sim)
{
	PositionsComp *positionsMetric = GetSpecificMetric<Metric, PositionsComp>(sim.GetAllMetrics());
	assert(positionsMetric);
	cellToCellDist.push_back(CellToCellDistData(
		positionsMetric->GetMeanDist(),
		positionsMetric->GetStdDevDist(),
		positionsMetric->GetMinDist()));
	return true;
}

//**********************************************************************
// Save the metric
//**********************************************************************
bool NetworkTopologyStats::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	static const string CellToCellDistName("CellToCellDistStats");
	if (saver.isSaving(CellToCellDistName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(CellToCellDistName, saver.getCurrFile());

		stream << "MeanDist\tStdDevDist\tMinDist" << endl;
		for (unsigned int i = 0 ; i < cellToCellDist.size() ; ++i)
			stream 
				<< cellToCellDist[i].meanDist << "\t"
				<< cellToCellDist[i].stdDevDist << "\t"
				<< cellToCellDist[i].minDist << endl;

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void NetworkTopologyStats::Initialize()
{
	Metric::Initialize();
	cellToCellDist.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * NetworkTopologyStats::BuildCopy() const
{
	return new NetworkTopologyStats(); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool NetworkTopologyStats::LoadFromStream(std::ifstream & stream)
{
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool NetworkTopologyStats::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//********************************************************************//
//********* S C A L A R   S T A T I S T I C S   M E T R I C **********//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************

//**********************************************************************
// Default Constructor
//**********************************************************************
ScalarStatisticsMetric::ScalarStatisticsMetric()
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ScalarStatisticsMetric::ScalarStatisticsMetric(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Default Destructor
//**********************************************************************
ScalarStatisticsMetric::~ScalarStatisticsMetric()
{
}

//**********************************************************************
// Compute the metric
//**********************************************************************
bool ScalarStatisticsMetric::ComputeMetric(
	const AbstractRepeatSimulation & sim)
{
	bool scalarsChanged = false;

	// First scan all metrics for easily retrievable scalars
	std::vector<Metric *> allMetrics = sim.GetAllMetrics();
	std::map<std::string, double> tmpMap;
	for (unsigned int i = 0 ; i < allMetrics.size() ; ++i)
	{
		if (allMetrics[i])
		{
			tmpMap = allMetrics[i]->GetScalarStatsToSave();
			for (std::map<std::string, double>::iterator it = tmpMap.begin() ;
					it != tmpMap.end() ; ++it)
			{
				if (it->second < DEFAULT_MAX_VAL-1)
				{
					fullScalars[modifStatName(allMetrics[i]->ModifStatName(it->first))].push_back(it->second);
					scalarsChanged = true;
				}
			}
		}
	}

	// Stimulated nodes
	std::string totStimParamName = "ChITotalStimulatedCells";
	std::set<unsigned int> defaultStimNodes;
	std::vector<StimulatedCells *> stimMetrs =
		GetAllSpecificMetrics<Metric,StimulatedCells>(sim.GetAllMetrics());
	if (not stimMetrs.empty())
	{
		std::set<unsigned int> cellInds;
		for (unsigned int i = 0 ; i < stimMetrs.size() ; ++i)
		{
			assert(stimMetrs[i]);
			// Track nodes stimulated by default stim strat (for stochastic resonance)
			if (stimMetrs[i]->GetCurrOrigin() == DefaultStimStrat::ClassName)
				defaultStimNodes.insert(stimMetrs[i]->GetCumulStimCells().begin(),
				stimMetrs[i]->GetCumulStimCells().end());

			cellInds.insert(stimMetrs[i]->GetCumulStimCells().begin(),
				stimMetrs[i]->GetCumulStimCells().end());
		}
		fullScalars[totStimParamName].push_back(cellInds.size());

		scalarsChanged = true;
	}

	// ChIModel activated cells
	ActivatedCells *actCells = GetSpecificMetric<Metric,
		ActivatedCells>(sim.GetAllMetrics());
	if (actCells)
	{
		if (fullScalars.find(totStimParamName) != fullScalars.end())
			fullScalars["TotalNonStimActivatedCells"].push_back(
				actCells->GetCumulActCells() - 
				fullScalars[totStimParamName].back());

		scalarsChanged = true;
	}

	// ChIModel propagation metric
	PropagationDistance *propDist = GetSpecificMetric<Metric,
		PropagationDistance>(sim.GetAllMetrics());
	if (propDist)
	{
		// Cells activated by central stimulation (for stochastic resonance)
		PropagationDistance refOnProp = *propDist;
		std::set<unsigned int> signalActCells;
		for (unsigned int i = 0 ; i < refOnProp.size() ; ++i)
		{
			if (defaultStimNodes.find(refOnProp[i].initiator) != defaultStimNodes.end())
				signalActCells.insert(refOnProp[i].involvedCells.begin(), 
					refOnProp[i].involvedCells.end());
		}
		fullScalars["ChIDefaultStimActivatedCells"].push_back(signalActCells.size());

		scalarsChanged = true;
	}

	// Propagation models activated cells
	AbstractCellStateSaver *stateCellSav = GetSpecificMetric<Metric, 
		AbstractCellStateSaver>(sim.GetAllMetrics());
	if (stateCellSav)
	{
		const std::map<std::string, std::set<unsigned int> > & cumStates = 
			stateCellSav->GetLastCumulStates();

		for (std::map<std::string, std::set<unsigned int> >::const_iterator
				it = cumStates.begin() ; it != cumStates.end() ; ++it)
			fullScalars[it->first + "_CUMUL"].push_back(it->second.size());

		fullScalars["PropagModelNbStimNodes"].push_back(
			stateCellSav->GetNbStimulatedNodes());

		scalarsChanged = true;
	}

	// Wave front metrics
	WaveFrontDetect *wavefrntmetr = GetSpecificMetric<Metric, 
		WaveFrontDetect>(sim.GetAllMetrics());
	if (wavefrntmetr)
	{
		fullDistribs["WaveFrontInfluxes"].push_back(
			wavefrntmetr->GetFrontInluxes());

		fullDistribs["NonActCellsInfluxes"].push_back(
			wavefrntmetr->GetNonActInfluxes());

		// size of distribs
		unsigned int nbPointFront = 0;
		unsigned int nbPointNonAct = 0;
		for (unsigned int i = 0 ;
				i < fullDistribs["WaveFrontInfluxes"].size() ; ++i)
			nbPointFront += fullDistribs["WaveFrontInfluxes"][i].size();
		for (unsigned int i = 0 ; 
				i < fullDistribs["NonActCellsInfluxes"].size() ; ++i)
			nbPointNonAct += fullDistribs["NonActCellsInfluxes"][i].size();
		// distribs declaration and initialization
		alglib::real_1d_array distrFront, distrNonAct;
		double *frontPtr  = new double[nbPointFront];
		double *nonActPtr = new double[nbPointNonAct];
		unsigned int tempInd = 0;
		for (unsigned int i = 0 ;
				i < fullDistribs["WaveFrontInfluxes"].size() ; ++i)
			for (unsigned int j = 0 ;
					j < fullDistribs["WaveFrontInfluxes"][i].size() ; ++j)
				frontPtr[tempInd++] = fullDistribs[
					"WaveFrontInfluxes"][i][j];
		tempInd = 0;
		for (unsigned int i = 0 ;
				i < fullDistribs["NonActCellsInfluxes"].size() ; ++i)
			for (unsigned int j = 0 ;
					j < fullDistribs["NonActCellsInfluxes"][i].size() ; ++j)
				nonActPtr[tempInd++] = fullDistribs[
					"NonActCellsInfluxes"][i][j];
		distrFront.setcontent(nbPointFront, frontPtr);
		distrNonAct.setcontent(nbPointNonAct, nonActPtr);
		// Free temp vectors
		delete [] frontPtr;
		delete [] nonActPtr;
		// Mann-whitney U test
		double bothtails, lefttail, righttail;
		alglib::mannwhitneyutest(distrFront, nbPointFront, 
			distrNonAct, nbPointNonAct, bothtails, lefttail, righttail);
		// bothtails null hyp : medians are equal
		// leftail null hyp : median of first sample is >= than second
		// righttail null hyp : median of first sample is <= than second
		scalarStats["WaveInfluxBothTailsPVal"] = bothtails;
		scalarStats["WaveInfluxLeftTailPVal"] = lefttail;
		scalarStats["WaveInfluxRightTailPVal"] = righttail;
		
		scalarsChanged = true;
	}

	// Threshold determination

	ThresholdDetermination *threshDeterMetr = GetSpecificMetric<Metric, 
		ThresholdDetermination>(sim.GetAllMetrics());
	if (threshDeterMetr)
	{
		fullScalars["ThreshDetermSecondNeighbThresh"].push_back(
			threshDeterMetr->GetSecondNeighbThresh());
		fullScalars["ThreshDetermSecondNeighbThreshUnder10s"].push_back(
			threshDeterMetr->GetSecondNeighbThreshUnderTime(10.0));
		fullScalars["ThreshDetermSecondNeighbThreshUnder15s"].push_back(
			threshDeterMetr->GetSecondNeighbThreshUnderTime(15.0));
		fullScalars["ThreshDetermSecondNeighbThreshUnder20s"].push_back(
			threshDeterMetr->GetSecondNeighbThreshUnderTime(20.0));
	}

	if (scalarsChanged)
		updateScalarStats();

	return true;
}

//**********************************************************************
// Save the metric
//**********************************************************************
bool ScalarStatisticsMetric::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	static const string scalarStatsName("ScalarStatistics");
	if (saver.isSaving(scalarStatsName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(scalarStatsName, saver.getCurrFile());

		for (std::map<std::string, double>::const_iterator it = 
				scalarStats.begin() ; it != scalarStats.end() ; ++it)
			stream << it->first << "\t";
		stream << std::endl;

		for (std::map<std::string, double>::const_iterator it = 
				scalarStats.begin() ; it != scalarStats.end() ; ++it)
			stream << it->second << "\t";
		stream << std::endl;

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void ScalarStatisticsMetric::Initialize()
{
	Metric::Initialize();
	fullScalars.clear();
	scalarStats.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * ScalarStatisticsMetric::BuildCopy() const
{
	return new ScalarStatisticsMetric(); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool ScalarStatisticsMetric::LoadFromStream(std::ifstream & stream)
{
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool ScalarStatisticsMetric::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//**********************************************************************
// Updates the scalar statistics from the full scalar distributions
//**********************************************************************
void ScalarStatisticsMetric::updateScalarStats()
{
	for (std::map<std::string, std::vector<double> >::iterator it = 
		fullScalars.begin() ; it != fullScalars.end() ; ++it)
	{
		scalarStats[it->first + "_Mean"] = ComputeMean(it->second);
		scalarStats[it->first + "_StdDev"] = ComputeStdDev(it->second);
	}
}

//**********************************************************************
// Make sure that the stat name is compatible with MySQL
//**********************************************************************
std::string ScalarStatisticsMetric::modifStatName(const std::string & str) const
{
	std::string tmpStr = str;
	size_t pos;
	while ((pos = tmpStr.find_first_of('.')) != std::string::npos)
		tmpStr[pos] = 'p';

	return tmpStr;
}

//********************************************************************//
//**** G R I D   S E A R C H   S I M U L A T I O N   M E T R I C *****//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
//ADD_METRIC_DEPENDENCY(GRID_SEARCH_METRIC, SCALAR_STATISTICS_METRIC)

//**********************************************************************
// Compute the metric
//**********************************************************************
bool GridSearchMetric::ComputeMetric(const GridSearchSimulation & manag)
{
	if (paramNames.size() != (manag.GetParamNames().size() + 
		manag.GetGlobalParams().size()))
	{
		paramNames = manag.GetParamNames();
		for (std::map<std::string, std::string>::const_iterator it =
				manag.GetGlobalParams().begin() ; 
				it != manag.GetGlobalParams().end() ; ++it)
			paramNames.push_back(it->first);
	}

	const vector<Metric*> subMetrics = manag.GetSubSimulation()->
		GetAllMetrics();

	ParamComb comb = manag.GetCurrParamComb();
	// Add global params
	for (std::map<std::string, std::string>::const_iterator it = 
			manag.GetGlobalParams().begin() ; 
			it != manag.GetGlobalParams().end() ; ++it)
		comb.push_back(it->second);

	vector<pair<string, string> > tempVect;
	for (unsigned int i = 0 ; i < subMetrics.size() ; ++i)
	{
		// Get saved file paths from submetrics
		tempVect = subMetrics[i]->PullSavedFilePaths();
		for (unsigned int j = 0 ; j < tempVect.size() ; ++j)
		{
			paths.push_back(tempVect[j].second);
			totalFileNames.insert(tempVect[j].first);
			savedFilesByParam[comb][tempVect[j].first].push_back(paths.size() - 1);
		}
	}

	// Scalar statistics
	ScalarStatisticsMetric *scalStats = GetSpecificMetric<Metric, ScalarStatisticsMetric>(manag.GetAllMetrics());

	if(scalStats)
	{
		for (std::map<std::string, double>::const_iterator it = 
			scalStats->GetScalarStats().begin() ; 
			it != scalStats->GetScalarStats().end() ; ++it)
		{
			scalarStatsByParam[comb].insert(*it);
			totalScalarNames.insert(it->first);
		}
	}

	return true;
}

//**********************************************************************
// Save the metric
//**********************************************************************
bool GridSearchMetric::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string indexesFile("GridFilePathIndexes");
	std::string pathFile("GridPaths");
	std::string scalarStatsFile("GridScalarStatistics");

	std::string condParamReplace("CondParamReplace");


	if (saver.isSaving(indexesFile))
	{
		std::ofstream & stream = saver.getStream();
		this->AddSavedFile(indexesFile, saver.getCurrFile());
		
		for (unsigned int i = 0 ; i < paramNames.size() ; ++i)
			stream << paramNames[i] << "\t";

		for (set<string>::iterator it = totalFileNames.begin() ; 
				it != totalFileNames.end() ; ++it)
			stream << *it << "\t";
		stream << endl;

		for (map<ParamComb, map<string, vector<int> > >::const_iterator it = 
			savedFilesByParam.begin() ; it != savedFilesByParam.end() ; ++it)
		{
			for (unsigned int i = 0 ; i < it->first.size() ; ++i)
				stream << it->first[i] << "\t";

			for (set<string>::iterator it2 = totalFileNames.begin() ; 
				it2 != totalFileNames.end() ; ++it2)
			{
				map<string, vector<int> >::const_iterator foundName;
				if ((foundName = it->second.find(*it2)) != it->second.end())
					stream << foundName->second[0] + 1 << "\t";
				else
					stream << 0 << "\t";
			}
			stream << endl;
		}

		allSaved &= stream.good();
	}

	if (saver.isSaving(pathFile))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(pathFile, saver.getCurrFile());

		for (unsigned int i = 0 ; i < paths.size() ; ++i)
			stream << paths[i] << endl;

		allSaved &= stream.good();
	}

	if (saver.isSaving(scalarStatsFile) and not scalarStatsByParam.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(scalarStatsFile, saver.getCurrFile());

		for (unsigned int i = 0 ; i < paramNames.size() ; ++i)
			stream << paramNames[i] << "\t";

		for (set<string>::iterator it = totalScalarNames.begin() ; 
				it != totalScalarNames.end() ; ++it)
			stream << *it << "\t";
		stream << endl;

		for (map<ParamComb, map<string, double> >::const_iterator it = 
			scalarStatsByParam.begin() ; 
			it != scalarStatsByParam.end() ; ++it)
		{
			for (unsigned int i = 0 ; i < it->first.size() ; ++i)
				stream << it->first[i] << "\t";

			for (set<string>::iterator it2 = totalScalarNames.begin() ; 
				it2 != totalScalarNames.end() ; ++it2)
			{
				map<string, double>::const_iterator foundName;
				if ((foundName = it->second.find(*it2)) != it->second.end())
					stream << foundName->second << "\t";
				else
					stream << "NaN" << "\t";
			}
			stream << endl;
		}

		allSaved &= stream.good();
	}

	if (saver.isSaving(condParamReplace))
	{
		std::ofstream & stream = saver.getStream();
		this->AddSavedFile(indexesFile, saver.getCurrFile());

		unsigned int strW = ceil(log(savedFilesByParam.size() *
			totalFileNames.size()) / log(10));
		unsigned int filePathInd = 0;

		stream << paramNames.size() << std::endl;
		for (unsigned int i = 0 ; i < paramNames.size() ; ++i)
			stream << paramNames[i] << std::endl;
		stream << totalFileNames.size() << std::endl;
		for (std::set<std::string>::const_iterator it = totalFileNames.begin() ; 
				it != totalFileNames.end() ; ++it)
			stream << *it << std::endl;

		stream << savedFilesByParam.size() << std::endl;
		for (std::map<ParamComb, std::map<std::string, std::vector<int> > >::const_iterator 
			it = savedFilesByParam.begin() ; it != savedFilesByParam.end() ; ++it)
		{
			// Write paramComb
			for (unsigned int i = 0 ; i < it->first.size() ; ++i)
				stream << it->first[i] << std::endl;
			// Write file paths if needed
			for (std::set<std::string>::const_iterator it2 = totalFileNames.begin() ; 
				it2 != totalFileNames.end() ; ++it2)
			{
				std::map<std::string, std::vector<int> >::const_iterator foundName;
				if ((foundName = it->second.find(*it2)) != it->second.end())
				{
					if (foundName->second.size() > 1)
					{
						// Save a separate file containing filepath list
						ResultSaver tmpSav = saver("condReplacePathLists");
						std::string tmpFileName = std::string("FilePathList") + 
							StringifyFixed(filePathInd++, strW);
						tmpSav <= tmpFileName;
						if (tmpSav.isSaving(tmpFileName))
						{
							std::ofstream & subStream = tmpSav.getStream();
							this->AddSavedFile(tmpFileName, tmpSav.getCurrFile());

							subStream << foundName->second.size() << std::endl;
							for (unsigned int i = 0 ; i < foundName->second.size() ; ++i)
								subStream << paths[foundName->second[i]] << std::endl;

							stream << tmpSav.getCurrFile() << std::endl;
						}
						else
							stream << "NULL" << std::endl;
					}
					else
						stream << paths[foundName->second[0]] << std::endl;
				}
				else
					stream << "NULL" << std::endl;
			}

		}
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void GridSearchMetric::Initialize()
{
	Metric::Initialize();
	paramNames.clear();
	savedFilesByParam.clear();
	scalarStatsByParam.clear();
	paths.clear();
	totalFileNames.clear();
	totalScalarNames.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * GridSearchMetric::BuildCopy() const
{
	return new GridSearchMetric(); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool GridSearchMetric::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool GridSearchMetric::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

