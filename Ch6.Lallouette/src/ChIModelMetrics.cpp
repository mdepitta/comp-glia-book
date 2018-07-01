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

#include "ChIModelMetrics.h"
#include "ChIModel.h"
#include "ChICell.h"
#include "MetricNames.h"

#include <set>
#include <cmath>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_interp.h>

using namespace AstroModel;
using namespace std;

string ActivatedCells::ClassName(ACTIVATED_CELLS_METRICS);
string ConcentrationMetrics::ClassName(CONCENTRATIONS_METRIC);
string PropagationDistance::ClassName(PROPAGATION_DISTANCE_METRIC);
string CorrelationsMetric::ClassName(CORRELATIONS_METRIC);
string TransferEntropyMetric::ClassName(TRANSFER_ENTROPY_METRIC);
string FunctTopoByNetConstrStrat::ClassName(FUNCT_TOPO_BY_NET_STRAT_METRIC);
string ThresholdDetermination::ClassName(THRESHOLD_DETERMINATION_METRIC);
string WaveFrontDetect::ClassName(WAVE_FRONT_DETECT_METRIC);
string FourierTransform::ClassName(FOURIER_TRANSFORM_METRIC);
string WaveletTransform::ClassName(WAVELET_TRANSFORM_METRIC);
string StochResMetric::ClassName(STOCH_RES_METRIC);

//********************************************************************//
//****************** A C T I V A T E D   C E L L S *******************//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(ACTIVATED_CELLS_METRICS, DEGREE_DISTRIBUTION_NETWORK_METRIC)

//**********************************************************************
// Default Constructor
//**********************************************************************
ActivatedCells::ActivatedCells(ParamHandler & h) : lastTime(0), 
	tStart(-DEFAULT_MAX_VAL), tEnd(-DEFAULT_MAX_VAL)
{
	activThresh = h.getParam<double>("-activParams", 0);
	freqEstTimeWin = h.getParam<double>("-activParams", 1);
	fluxComputingDelay = h.getParam<double>("-activParams", 2);
}

//**********************************************************************
// Full constructor
//**********************************************************************
ActivatedCells::ActivatedCells(double _at, double _fetw, double _fcd) : 
	activThresh(_at), freqEstTimeWin(_fetw), fluxComputingDelay(_fcd), 
	totNbCells(0), lastTime(0), tStart(-DEFAULT_MAX_VAL), tEnd(-DEFAULT_MAX_VAL)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ActivatedCells::ActivatedCells(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Compute the activated cells
//**********************************************************************
//bool ActivatedCells::ComputeMetric(const ChIModel & model)
bool ActivatedCells::ComputeMetric(const AbstractODENetDynProblem & model)
{
	if (lastActivatedTime.size() != model.GetNbCells())
		lastActivatedTime = std::vector<double>(model.GetNbCells(), -1.0); 
	if (tmpOutFluxes.size() != model.GetNbCells())
		tmpOutFluxes = std::vector<double>(model.GetNbCells(), 0.0); 
	if (tmpInFluxes.size() != model.GetNbCells())
		tmpInFluxes = std::vector<std::vector<double> >(model.GetNbCells(), 
			std::vector<double>((int)(fluxComputingDelay / model.GetIntegrStep()), 0.0));

	if (totalInFluxes.size() != model.GetNbCells())
		totalInFluxes = std::vector<std::vector<double> >(
			model.GetNbCells(), std::vector<double>());
	if (maxInFluxInTimeWin.size() != model.GetNbCells())
		maxInFluxInTimeWin = std::vector<std::pair<double, double> >(
			model.GetNbCells(), std::make_pair(0.0, 0.0));
	if (neighborSpikeTimes.size() != model.GetNbCells())
		neighborSpikeTimes = std::vector<std::queue<double> >(
			model.GetNbCells(), std::queue<double>());
	if (influxAfterSpike.size() != model.GetNbCells())
		influxAfterSpike = std::vector<std::vector<double> >(
			model.GetNbCells(), std::vector<double>());

	if (tStart == -DEFAULT_MAX_VAL)
		tStart = model.GetTStart();
	if (tEnd == -DEFAULT_MAX_VAL)
		tEnd = model.GetTEnd();

	double t = model.GetTime();
	double tmpFlux;
	int nOut = 1;

	if (((int)(t / model.GetIntegrStep())) % nOut == 0)
	{
		totNbCells = model.GetNbCells();
		for (unsigned int i = 0 ; i < model.GetNbCells() ; ++i)
		{
			tmpFlux = std::max(0.0, - model.GetTotalFlux(i)) * (t - lastTime);
			// Compute maxInFluxInTimeWin
			maxInFluxInTimeWin[i].first += tmpFlux - 
				tmpInFluxes[i][((int)(t / model.GetIntegrStep())) % tmpInFluxes[i].size()];
			if (maxInFluxInTimeWin[i].first > maxInFluxInTimeWin[i].second)
				maxInFluxInTimeWin[i].second = maxInFluxInTimeWin[i].first;

			// In any case compute influxes
			tmpInFluxes[i][((int)(t / model.GetIntegrStep())) % tmpInFluxes[i].size()] = tmpFlux;

			// If a neighbor has just spiked exactly fluxComputingDelay seconds ago
			while (not neighborSpikeTimes[i].empty() and neighborSpikeTimes[i].front() + fluxComputingDelay <= t) {
				influxAfterSpike[i].push_back(ComputeSum(tmpInFluxes[i]));
				neighborSpikeTimes[i].pop();
			}

			if (model.GetExcDynVal(i) > activThresh)
			{
				// Check if it has already been marked as activated
				if (lastActivatedTime[i] < 0)
				{
					lastActivatedTime[i] = t;
					tmpOutFluxes[i] = 0;
					// Compute cumulated influxes
					totalInFluxes[i].push_back(ComputeSum(tmpInFluxes[i]));
					// Add Spike time to neighbors
					const std::vector<unsigned int> & neighbs = model.GetNeighbors(i);
					for (unsigned int j = 0 ; j < neighbs.size() ; ++j)
						neighborSpikeTimes[neighbs[j]].push(t);
				}
				// if last added time is different
				if (activations.empty() or (activations.back().time != model.GetTime()))
					activations.push_back(TimedActivations(model.GetTime()));
				activations.back().cells.push_back(i);
				cumulActivCells.insert(i);
			}


			// If cell is marked as recently activated, compute its fluxes
			if (lastActivatedTime[i] >= 0)
			{
				// If delay of fluxes computing is over, un mark the cell and store the flux and degree
				if ((t - lastActivatedTime[i]) > fluxComputingDelay)
				{
					// Check that the cell is not stimulated and not a neighbor of a stimulated cell
					bool discard = model.IsStimulated(i);
					const std::vector<unsigned int> & neighbs = model.GetNeighbors(i);
					for (unsigned int j = 0 ; (not discard) and (j < neighbs.size()) ; ++j)
						discard |= model.IsStimulated(neighbs[j]);
					if (not discard)
					{
						DegreeDistrComp  & degreeMetr = 
							*GetSpecificMetric<Metric, DegreeDistrComp>(model.GetAllMetrics());
						totalFluxes.push_back(std::make_pair(degreeMetr[i], tmpOutFluxes[i]));
					}
					lastActivatedTime[i] = -1.0;
					tmpOutFluxes[i] = 0;
				}
				else
					tmpOutFluxes[i] += std::max(0.0, model.GetTotalFlux(i)) * (t - lastTime);

			}
		}
	}

	lastTime = t;
	return true;
}

//**********************************************************************
// Save the activated cells
//**********************************************************************
bool ActivatedCells::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string activCellsName("ActivatedCells");
	if (saver.isSaving(activCellsName) and not activations.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(activCellsName, saver.getCurrFile());

		stream << "Time\tNbCurrActiv\tNbCumulActiv" << endl;

		set<unsigned int> activatedCells;
		unsigned int lastNbCurr = 0;
		unsigned int lastNbCumul = 0;
		for (unsigned int i = 0 ; i < activations.size() ; ++i)
		{
			for (unsigned int j = 0 ; j < activations[i].cells.size() ; ++j)
				activatedCells.insert(activations[i].cells[j]);

			if ((activations[i].cells.size() != lastNbCurr) or
				(activatedCells.size() != lastNbCumul))
			{
				lastNbCurr = activations[i].cells.size();
				lastNbCumul = activatedCells.size();
				stream 
					<< activations[i].time << "\t"
					<< activations[i].cells.size() << "\t"
					<< activatedCells.size() << endl;
			}
		}

		allSaved &= stream.good();
	}

	std::string outFluxesName("OutFluxesFull");
	if (saver.isSaving(outFluxesName) and not totalFluxes.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(outFluxesName, saver.getCurrFile());

		stream << "k\tOutFlux" << endl;

		for (unsigned int i = 0 ; i < totalFluxes.size() ; ++i)
			stream << totalFluxes[i].first << "\t" << totalFluxes[i].second << std::endl;

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void ActivatedCells::Initialize()
{
	Metric::Initialize();
	activations.clear();
	freqEstim.clear();
	cumulActivCells.clear();
	lastActivatedTime.clear();
	tmpOutFluxes.clear();
	tmpInFluxes.clear();
	totalFluxes.clear();
	totalInFluxes.clear();
	maxInFluxInTimeWin.clear();
	neighborSpikeTimes.clear();
	influxAfterSpike.clear();
	totNbCells = 0;
	meanFreqEstim = 0;
	lastTime = 0;
	tStart = -DEFAULT_MAX_VAL;
	tEnd = -DEFAULT_MAX_VAL;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * ActivatedCells::BuildCopy() const
{
	return new ActivatedCells(activThresh, freqEstTimeWin, fluxComputingDelay); 
}

//**********************************************************************
//**********************************************************************
ParamHandler ActivatedCells::BuildModelParamHandler()
{
	ParamHandler params;

	params <= "ActivationThreshold", activThresh;
	params <= "FrequencyEstimTimeWin", freqEstTimeWin;

	return params;
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> ActivatedCells::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	res["TotalNbCellsActivated"] = GetCumulActCells();
	res["MeanActivationFrequency"] = GetMeanFreqEst();
	return res;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool ActivatedCells::LoadFromStream(std::ifstream & stream)
{
	stream >> activThresh;
	stream >> freqEstTimeWin;
	stream >> fluxComputingDelay;
	Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool ActivatedCells::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< activThresh << endl
		<< freqEstTimeWin << endl
		<< fluxComputingDelay << endl;
	return stream.good();
}

//**********************************************************************
// Compute and return the time dependent frequency estimations
//**********************************************************************
const std::vector<std::pair<double, double> > & 
	ActivatedCells::CompAndGetFreqEst() const
{
	if (freqEstim.empty())
	{
		double nbSpikes = 0;
		std::map<unsigned int, bool> actCells;
		// Compute Mean Freq Estim
		nbSpikes = 0;
		actCells.clear();
		for (unsigned int j = 0 ; j < activations.size() ; ++j)
		{
				for (unsigned int k = 0 ; k < activations[j].cells.size() ; ++k)
				{
					if ((actCells.find(activations[j].cells[k]) == 
						actCells.end())
						or not actCells[activations[j].cells[k]])
					{
						actCells[activations[j].cells[k]] = true;
						nbSpikes += 1.0;
					}
				}
				for (std::map<unsigned int, bool>::iterator it = 
						actCells.begin(); it != actCells.end() ; ++it)
					if (std::find(activations[j].cells.begin(), 
							activations[j].cells.end(), it->first) == 
							activations[j].cells.end())
						it->second = false;
		}
		meanFreqEstim = nbSpikes / 
			((double)totNbCells * (activations.empty() ? 1 : (tEnd - tStart)));
	}
	return freqEstim;
}

//**********************************************************************
// Return the mean frequency estimation
//**********************************************************************
double ActivatedCells::GetMeanFreqEst() const
{
	CompAndGetFreqEst();
	return meanFreqEstim;
}

//**********************************************************************
// Return the mean influx in cells before firing 
//**********************************************************************
double ActivatedCells::ComputeMeanInflux(unsigned int i) const
{
	assert(i < totalInFluxes.size());
	return ComputeMean(totalInFluxes[i]);
}

//**********************************************************************
//**********************************************************************
double ActivatedCells::GetMaxInfluxInTimeWin(unsigned int i) const
{
	assert(i < maxInFluxInTimeWin.size());
	return maxInFluxInTimeWin[i].second;
}

//**********************************************************************
//**********************************************************************
double ActivatedCells::GetMeanInfluxAfterSpike(unsigned int i) const
{
	assert(i < influxAfterSpike.size());
	return ComputeMean(influxAfterSpike[i]);
}

//********************************************************************//
//************* P R O P A G A T I O N   D I S T A N C E **************//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(PROPAGATION_DISTANCE_METRIC, ACTIVATED_CELLS_METRICS)

//**********************************************************************
// Default constructor
//**********************************************************************
PropagationDistance::PropagationDistance(ParamHandler & h) : 
	indActiv(0), lastTime(0)
{
	minDelay = h.getParam<double>("-propDelays", 0), 
	maxDelay = h.getParam<double>("-propDelays", 1),
	instDistTimeWindow = h.getParam<double>("-propDelays", 2);
	cellDistIncr = h.getParam<double>("-propDelays", 3);

	allMergeInOne = h.getParam<bool>("-propOnlyOneWave", 0);

	startQuant = h.getParam<std::vector<double> >("-propSaveQuant", 0);
	endQuant   = h.getParam<std::vector<double> >("-propSaveQuant", 1);
	stepQuant  = h.getParam<std::vector<double> >("-propSaveQuant", 2);

	startRatio = h.getParam<std::vector<double> >("-propSaveRatio", 0);
	endRatio   = h.getParam<std::vector<double> >("-propSaveRatio", 1);
	stepRatio  = h.getParam<std::vector<double> >("-propSaveRatio", 2);
}

//**********************************************************************
// Full Constructor
//**********************************************************************
PropagationDistance::PropagationDistance(double _miD, double _maD, 
	double _cdi, double _idtw, bool _amio) :
	indActiv(0), maxDelay(_maD), minDelay(_miD), instDistTimeWindow(_idtw),
	cellDistIncr(_cdi), allMergeInOne(_amio), lastTime(0)
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
PropagationDistance::PropagationDistance(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Default destructor
//**********************************************************************
PropagationDistance::~PropagationDistance()
{
}

//**********************************************************************
// Compute propagation distances
//**********************************************************************
bool PropagationDistance::ComputeMetric(const ChIModel & model)
{
	ActivatedCells *activCells = 
		GetSpecificMetric<Metric, ActivatedCells>(model.GetAllMetrics());
	assert(activCells);

	if (isInWaveSince.size() != model.GetNbCells())
		isInWaveSince = vector<pair<bool, double> >(model.GetNbCells(),
			make_pair<bool, double>(false, 0));

	const std::vector<TimedActivations> & activ = activCells->GetActivations();
	unsigned int currCell;
	double t = model.GetTime();
	lastTime = t;
	// For each non analysed activation steps
	for (; indActiv < activ.size() ; ++indActiv)
	{
		// For each activated cells in the step
		for (unsigned int i = 0 ; i < activ[indActiv].cells.size() ; ++i)
		{
			currCell = activ[indActiv].cells[i];
			// If the cell has not already been counted
			if ((not isInWaveSince[currCell].first) or (isInWaveSince[currCell].second + minDelay < t))
			{
				isInWaveSince[currCell].first  = false;
				// For each wave
				for (unsigned int w = 0 ; (w < waves.size()) and not isInWaveSince[currCell].first ; ++w)
				{
					// If the currCell was the initiator of the wave and the wave isn't finished
					if ((currCell == waves[w].initiator) and (waves[w].structure.back().lastActiv + maxDelay >=t))
					{
						waves[w].structure.push_back(CellToCellPropagation(currCell, t));
						isInWaveSince[currCell].first  = true;
					}
					else
					{
						// For each cells in the wave w
						unsigned int wCell = waves[w].structure[0].ind;
						for (unsigned int f = 0 ; f < waves[w].structure.size() ; ++f)
						{
							wCell = waves[w].structure[f].ind;
							// If the current cell is linked to a cell in the wave w that was activated not too long ago
							if (allMergeInOne or ((waves[w].structure[f].lastActiv + maxDelay >= t) and 
								model.GetNetwork().AreConnected(currCell, wCell)))
							{
								// If the current cell has already been added
								if (waves[w].structure.back().ind == currCell)
								{
									// Add a parent to the current cell in the wave
									// if it is not already a parent
									if (find(waves[w].structure.back().parents.begin(),
											waves[w].structure.back().parents.end(), wCell) == 
											waves[w].structure.back().parents.end())
										waves[w].structure.back().parents.push_back(wCell);
									// Add currCell as a children to the parent
									waves[w].structure[f].children.push_back(currCell);
									// Update cell metric measure (cellDist)
									waves[w].structure.back().distToInitiator = 
										min(waves[w].structure.back().distToInitiator,
										waves[w].structure[f].distToInitiator + cellDistIncr);
									waves[w].propagation.back().cellDist = 
										max(waves[w].propagation.back().cellDist,
										waves[w].structure.back().distToInitiator);
								}
								else
								{
									// Add the current cell to the wave w
									double cellDistParent = waves[w].structure[f].distToInitiator;
									CellToCellPropagation tempCell(currCell, wCell, cellDistParent + cellDistIncr, t);
									waves[w].AddCell(tempCell);

									double spatialDist = max(waves[w].propagation.back().spatialDist, 
										model.GetDistance(waves[w].initiator, currCell));

									double instSpatialDist = model.GetDistance(waves[w].initiator, currCell);
									for (unsigned int k = 0 ; k < waves[w].propagation.size() ; ++k)
									{
										if (t - waves[w].propagation[k].time < instDistTimeWindow)
											instSpatialDist = max(waves[w].propagation[k].instSpatialDist, instSpatialDist);
									}

									double cellDist = max(waves[w].propagation.back().cellDist, 
										waves[w].structure.back().distToInitiator);

									// Add propagation mark
									waves[w].propagation.push_back(TimedPropagation(t, spatialDist, instSpatialDist, cellDist, waves[w].GetCellNb()));
								}
								// Add the cell as a child of it's latest added parent
								waves[w].structure[f].children.push_back(currCell);
								isInWaveSince[currCell].first  = true;
							}
						}
					}
					// if the cell was added to the current wave and all its parents have been added
					if (isInWaveSince[currCell].first)
					{
						double averageDelay = 0;
						double averageSpatialDist = 0;
						unsigned int nbParentsConsidered = 0;
						// For each parent cell
						for (unsigned int p = 0 ; p < waves[w].structure.back().parents.size() ; ++p)
						{
							double parentActTime = -1;
							// Find corresponding entry in waves[w].structure, starting from the end
							for (int pp = waves[w].structure.size() - 1 ; (pp >= 0) and (parentActTime == -1) ; --pp)
							{
								// Consider only parents which are closer to the initiator
								if ((waves[w].structure[pp].ind == waves[w].structure.back().parents[p]) 
									and (waves[w].structure[pp].distToInitiator < waves[w].structure.back().distToInitiator))
								{
									parentActTime = waves[w].structure[pp].lastActiv;
									averageDelay += t - parentActTime;
									averageSpatialDist += model.GetDistance(waves[w].structure[pp].ind, currCell);
									++nbParentsConsidered;
								}
							}
						}
						if ((nbParentsConsidered > 0) and (averageDelay > 0))
						{
							averageDelay /= (double) nbParentsConsidered;
							averageSpatialDist /= (double) nbParentsConsidered;
							waves[w].structure.back().localTopoSpeed = 1.0 / averageDelay;
							waves[w].structure.back().localSpatialSpeed = averageSpatialDist / averageDelay;
						}
					}
				}
				// If the cell couldn't be linked to any wave
				if (not isInWaveSince[currCell].first)
				{
					// Create a new wave with the cell as an initiator
					waves.push_back(TimedWave(currCell, t));

					isInWaveSince[currCell].first  = true;
				}
				isInWaveSince[currCell].second = t;
			}
		}
	}

	return true;
}

//**********************************************************************
// Save propagation distances
//**********************************************************************
bool PropagationDistance::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string fullWavesName("FullWaves");
	std::string waveMetricsName("WaveMetrics");
	std::string waveSummaryName("WaveSummary");
	std::string spikingCoherenceName("SpikingCoherence");
	if (saver.isSaving(fullWavesName) and not waves.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(fullWavesName, saver.getCurrFile());

		stream << "WaveId\tCellId\tTime\tParents\tChildren" << endl;
		
		for (unsigned int w = 0 ; w < waves.size() ; ++w)
			for (unsigned int c = 0 ; c < waves[w].structure.size() ; ++c)
				stream
					<< w << "\t"
					<< waves[w].structure[c].ind << "\t"
					<< waves[w].structure[c].lastActiv << "\t"
					<< waves[w].structure[c].parents.size() << "\t"
					<< waves[w].structure[c].children.size() << endl;

		allSaved &= stream.good();
	}
	if (saver.isSaving(waveMetricsName) and not waves.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(waveMetricsName, saver.getCurrFile());

		stream << "WaveId\tTime\tSpatialDist\tInstSpatialDist\tCellDist\tCellNb\tWaveSpeed" << endl;
		for (unsigned int w = 0 ; w < waves.size() ; ++w)
			for (unsigned int p = 0 ; p < waves[w].propagation.size() ; ++p)
				stream
					<<  w << "\t"
					<<  waves[w].propagation[p].time << "\t"
					<<  waves[w].propagation[p].spatialDist << "\t"
					<<  waves[w].propagation[p].instSpatialDist << "\t"
					<<  waves[w].propagation[p].cellDist << "\t"
					<<  waves[w].propagation[p].cellNb << "\t"
					<< (((p == 0) or (waves[w].propagation[p].time - waves[w].propagation[0].time == 0)) ? 0 :
						(waves[w].propagation[p].spatialDist / 
						(waves[w].propagation[p].time - waves[w].propagation[0].time))) << endl;

		allSaved &= stream.good();
	}
	if (saver.isSaving(waveSummaryName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(waveSummaryName, saver.getCurrFile());

		stream << "WaveId\tSize\tTotalSize\tMeanTopoSpeed\tMeanSpatialSpeed\tTopoExtent\tSpatialExtent" << endl;
		for (unsigned int w = 0 ; w < waves.size() ; ++w)
		{
			double meanTopoSpeed = 0;
			double meanSpatialSpeed = 0;
			double maxSpatialDist = 0;
			double maxCellDist = 0;
			double nbData = 0;
			for (unsigned int i = 0 ; i < waves[w].structure.size() ; ++i)
				if (waves[w].structure[i].localTopoSpeed != 0)
				{
					meanTopoSpeed += waves[w].structure[i].localTopoSpeed;
					meanSpatialSpeed += waves[w].structure[i].localSpatialSpeed;
					++nbData;
				}
			meanTopoSpeed /= nbData ? nbData : 1;
			meanSpatialSpeed /= nbData ? nbData : 1;
			maxSpatialDist = waves[w].propagation.back().spatialDist;
			maxCellDist = waves[w].propagation.back().cellDist;
			stream 
				<< w << "\t"
				<< waves[w].involvedCells.size() << "\t"
				<< waves[w].propagation.size() << "\t"
				<< meanTopoSpeed << "\t"
				<< meanSpatialSpeed << "\t"
				<< maxCellDist << "\t"
				<< maxSpatialDist << endl;
		}

		allSaved &= stream.good();
	}
	if (saver.isSaving(spikingCoherenceName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(spikingCoherenceName, saver.getCurrFile());

		stream << "WaveId\tTopoDist\tNbActCells\tStdDevFirstActTimes\ttMin\ttMax\tMaxMinFirstActTimesRatio" << endl;
		for (unsigned int w = 0 ; w < waves.size() ; ++w)
			for (unsigned int r = 0 ; r < waves[w].firstActCoherence.size() ; ++r)
				stream 
					<< w << "\t"
					<< r << "\t"
					<< waves[w].firstActCoherence[r].GetNbActCels() << "\t"
					<< waves[w].firstActCoherence[r].stdDevActTimes << "\t"
					<< waves[w].firstActCoherence[r].minTime << "\t"
					<< waves[w].firstActCoherence[r].maxTime << "\t"
					<< waves[w].firstActCoherence[r].maxMinActTimesRatio << std::endl;

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void PropagationDistance::Initialize()
{
	Metric::Initialize();
	indActiv = 0;
	isInWaveSince.clear();
	waves.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * PropagationDistance::BuildCopy() const
{
	return new PropagationDistance(minDelay, maxDelay, 
		cellDistIncr, instDistTimeWindow, allMergeInOne); 
}

//**********************************************************************
// Allow access to waves given their indices
//**********************************************************************
const TimedWave & PropagationDistance::operator[](unsigned int ind) const
{
	assert(ind < waves.size());
	return waves[ind];
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> PropagationDistance::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;

	std::vector<double> topoSpeeds, spatialSpeeds, spatialDists, 
		cellDists, distinctCells, totalCells;

	for (unsigned int w = 0 ; w < waves.size() ; ++w)
	{
		double meanTopoSpeed = 0;
		double meanSpatialSpeed = 0;
		double maxSpatialDist = 0;
		double maxCellDist = 0;
		double nbData = 0;
		for (unsigned int i = 0 ; i < waves[w].structure.size() ; ++i)
			if (waves[w].structure[i].localTopoSpeed != 0)
			{
				meanTopoSpeed += waves[w].structure[i].localTopoSpeed;
				meanSpatialSpeed += waves[w].structure[i].localSpatialSpeed;
				++nbData;
			}
		meanTopoSpeed /= nbData ? nbData : 1;
		meanSpatialSpeed /= nbData ? nbData : 1;
		maxSpatialDist = waves[w].propagation.back().spatialDist;
		maxCellDist = waves[w].propagation.back().cellDist;
		
		topoSpeeds.push_back(meanTopoSpeed);
		spatialSpeeds.push_back(meanSpatialSpeed);
		spatialDists.push_back(maxSpatialDist);
		cellDists.push_back(maxCellDist);
		distinctCells.push_back(waves[w].involvedCells.size());
		totalCells.push_back(waves[w].propagation.size());
	}

	std::vector<double> quants;
	std::vector<double> ratios;
	assert(startQuant.size() == endQuant.size());
	assert(endQuant.size() == stepQuant.size());
	assert(startRatio.size() == endRatio.size());
	assert(endRatio.size() == stepRatio.size());
	for (unsigned int i = 0 ; i < startQuant.size() ; ++i)
	{
		assert(stepQuant[i] > 0);
		double q = startQuant[i];
		do {
			quants.push_back(q);
			q += stepQuant[i];
		} while (q <= (endQuant[i] + DEFAULT_SMALL_VAL));
	}
	for (unsigned int i = 0 ; i < startRatio.size() ; ++i)
	{
		assert(stepRatio[i] > 0);
		double q = startRatio[i];
		do {
			ratios.push_back(q);
			q += stepRatio[i];
		} while (q <= (endRatio[i] + DEFAULT_SMALL_VAL));
	}

	ComputeAndAddQuantiles(res, "PropDist_TopoSpeed", topoSpeeds, quants);
	ComputeAndAddQuantiles(res, "PropDist_SpatialSpeed", spatialSpeeds, quants);
	ComputeAndAddQuantiles(res, "PropDist_SpatialDist", spatialDists, quants);
	ComputeAndAddQuantiles(res, "PropDist_CellDist", cellDists, quants);
	ComputeAndAddQuantiles(res, "PropDist_DistinctCellNb", distinctCells, quants);
	ComputeAndAddQuantiles(res, "PropDist_TotalCellNb", totalCells, quants);

	ComputeAndAddRatioFreq(res, "PropDist_TopoSpeed", topoSpeeds, ratios, lastTime);
	ComputeAndAddRatioFreq(res, "PropDist_SpatialSpeed", spatialSpeeds, ratios, lastTime);
	ComputeAndAddRatioFreq(res, "PropDist_SpatialDist", spatialDists, ratios, lastTime);
	ComputeAndAddRatioFreq(res, "PropDist_CellDist", cellDists, ratios, lastTime);
	ComputeAndAddRatioFreq(res, "PropDist_DistinctCellNb", distinctCells, ratios, lastTime);
	ComputeAndAddRatioFreq(res, "PropDist_TotalCellNb", totalCells, ratios, lastTime);

	res["PropDist_TotalNbEvents"] = waves.size();
	return res;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool PropagationDistance::LoadFromStream(std::ifstream & stream)
{
	bool ok = true;
	stream >> maxDelay;
	stream >> minDelay;
	stream >> instDistTimeWindow;
	stream >> cellDistIncr;
	stream >> allMergeInOne;
	LoadVectorFromStream(startQuant, stream);
	LoadVectorFromStream(endQuant, stream);
	LoadVectorFromStream(stepQuant, stream);
	LoadVectorFromStream(startRatio, stream);
	LoadVectorFromStream(endRatio, stream);
	LoadVectorFromStream(stepRatio, stream);
	this->Initialize();
	return ok and stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool PropagationDistance::SaveToStream(std::ofstream & stream) const
{
	bool ok = true;
	stream 
		<< maxDelay           << endl
		<< minDelay           << endl
		<< instDistTimeWindow << endl
		<< cellDistIncr       << endl
		<< allMergeInOne      << endl;
	SaveVectorToStream(startQuant, stream);
	SaveVectorToStream(endQuant, stream);
	SaveVectorToStream(stepQuant, stream);
	SaveVectorToStream(startRatio, stream);
	SaveVectorToStream(endRatio, stream);
	SaveVectorToStream(stepRatio, stream);
	return ok and stream.good();
}

//********************************************************************//
//************ C O N C E N T R A T I O N   M E T R I C S *************//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
ConcentrationMetrics::ConcentrationMetrics(ParamHandler & h) : nOut(0)
{
    savingStep = h.getParam<double>("-SavingStep", 0);
}

//**********************************************************************
// Full constructor
//**********************************************************************
ConcentrationMetrics::ConcentrationMetrics(double _ss) : 
	nOut(0), savingStep(_ss)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ConcentrationMetrics::ConcentrationMetrics(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Compute propagation distances
//**********************************************************************
//bool ConcentrationMetrics::ComputeMetric(const ChIModel & model)
bool ConcentrationMetrics::ComputeMetric(const ODE::ODEProblem<double, 
	double> & prob)
{
	if (valNames.empty())
		valNames = prob.valNames;
	double t = prob.GetTime();
	if (nOut == 0)
		nOut = max(1, (int) floor(savingStep / prob.GetIntegrStep()));

	if (((int)(t / prob.GetIntegrStep())) % nOut == 0)
	{
		concentrations.push_back(make_pair(prob.GetTime(), 
			std::vector<double>(prob.GetNbVals(), 0)));
		for (unsigned int c = 0 ; c < prob.GetNbVals() ; ++c)
			concentrations.back().second[c] = prob.GetVal(c);
	}

	return true;
}

//**********************************************************************
// Save propagation distances
//**********************************************************************
bool ConcentrationMetrics::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string concentrationsName("Concentrations");
	if (saver.isSaving(concentrationsName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(concentrationsName, saver.getCurrFile());

		stream << "Time";
		for (unsigned int i = 0 ; i < valNames.size() ; ++i)
			stream << "\t" << valNames[i];
		stream << endl;

		for (unsigned int i = 0 ; i < concentrations.size() ; ++i)
		{
			stream << concentrations[i].first;
			for (unsigned int c = 0 ; c < concentrations[i].second.size() ; ++c)
				stream << "\t" << concentrations[i].second[c];
			stream << endl;
		}

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void ConcentrationMetrics::Initialize()
{
	Metric::Initialize();
	nOut = 0;
	concentrations.clear();
	valNames.clear();
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
map<std::string, double> ConcentrationMetrics::GetScalarStatsToSave() const
{
	map<string, double> res;

	map<string, vector<unsigned int> > valTypes;
	for (unsigned int i = 0 ; i < valNames.size() ; ++i)
		valTypes[valNames[i].substr(valNames[i].find_last_of('_')+1)].push_back(i);

	for (map<string, vector<unsigned int> >::iterator it = valTypes.begin() ; it != valTypes.end() ; ++it)
	{
		res[it->first + "_Vals"] = 0;
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			for (unsigned int t = 0 ; t < concentrations.size() ; ++t)
			{
				res[it->first + "_Vals"] += concentrations[t].second[it->second[i]];
			}
		res[it->first + "_Vals"] /= (double)(it->second.size() * concentrations.size());
	}
	return res;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * ConcentrationMetrics::BuildCopy() const
{
	return new ConcentrationMetrics(savingStep); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool ConcentrationMetrics::LoadFromStream(std::ifstream & stream)
{
	stream >> savingStep;
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool ConcentrationMetrics::SaveToStream(std::ofstream & stream) const
{
	stream << savingStep << std::endl;
	return stream.good();
}

//********************************************************************//
//********** F U N C T I O N A L   T O P O   M E T R I C *************//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(TRANSFER_ENTROPY_METRIC, DEGREE_DISTRIBUTION_NETWORK_METRIC)
ADD_METRIC_DEPENDENCY(TRANSFER_ENTROPY_METRIC, ACTIVATED_CELLS_METRICS)
ADD_METRIC_DEPENDENCY(CORRELATIONS_METRIC, ACTIVATED_CELLS_METRICS)

//**********************************************************************
// Default Constructor
//**********************************************************************
FunctionalTopoMetric::FunctionalTopoMetric(ParamHandler & h)
{
	functTopoUseSpecmetr = h.getParam<bool>("-FTopoSpecifMetrs", 0);
	functTopoMetrList = h.getParam<std::vector<std::string> >(
		"-FTopoSpecifMetrs", 1);
	threshEstimMethod = h.getParam<int>("-FTopoCommon", 0);
	threshSpecifMeanDeg = h.getParam<double>("-FTopoCommon", 1);
	stdDevCoeff = h.getParam<double>("-FTopoCommon", 2);
	modVoroMinDecile = h.getParam<double>("-FTopoCommon", 3);
	computeL = h.getParam<bool>("-FTopoComputeL", 0);
	computeHCC = h.getParam<bool>("-FTopoComputeHCC", 0);
	useMinForBidir = h.getParam<bool>("-FTopoUseMinForBidir", 0);
	voroConstrStrat = new VoronoiConstructStrat(h);
}

//**********************************************************************
// Full constructor
//**********************************************************************
FunctionalTopoMetric::FunctionalTopoMetric(int _m, double _t, double _sdC, 
	double _mvmd, bool _cl, bool _chcc, bool _umfb, bool _uSM, 
	const std::vector<std::string> & _fTML, VoronoiConstructStrat *_vcs) : 
	threshEstimMethod(_m), threshSpecifMeanDeg(_t), stdDevCoeff(_sdC), 
	modVoroMinDecile(_mvmd), computeL(_cl), computeHCC(_chcc), useMinForBidir(_umfb),
	functTopoUseSpecmetr(_uSM), functTopoMetrList(_fTML), 
	voroConstrStrat(_vcs)
{
	if (not voroConstrStrat)
		voroConstrStrat = new VoronoiConstructStrat();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
FunctionalTopoMetric::FunctionalTopoMetric(std::ifstream & stream)
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Destructor
//**********************************************************************
FunctionalTopoMetric::~FunctionalTopoMetric()
{
	for (std::map<std::string, Network<BooleanLink> *>::iterator it = 
		functTopos.begin() ; it != functTopos.end() ; ++it)
	{
		if (it->second)
			delete it->second;
	}
	delete voroConstrStrat;
}

//**********************************************************************
//**********************************************************************
bool FunctionalTopoMetric::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	for (std::map<std::string, Network<BooleanLink> *>::const_iterator it = 
		functTopos.begin() ; it != functTopos.end() ; ++it)
	{
		if (it->second)
			it->second->SaveMetrics(saver(it->first));
	}

	std::string connecByThreshName;
	for (std::map<std::string, 	
		std::vector<ThreshInfo> >::const_iterator
		it = threshInf.begin() ; it != threshInf.end() ; ++it)
	{
		connecByThreshName = it->first + "_ConnecByThresh";
		if (saver.isSaving(connecByThreshName))
		{
			ofstream & stream = saver.getStream();
			this->AddSavedFile(connecByThreshName, saver.getCurrFile());

			stream << "Threshold\tMeanDeg\tMeanPathLength\tUnconRatio\tHCCEstim\tEfficiency" 
				<< std::endl;
			for (unsigned int i = 0 ; i < it->second.size() ; ++i)
				stream 
					<< it->second[i].thresh << "\t" 
					<< it->second[i].connec << "\t" 
					<< it->second[i].meanPath << "\t" 
					<< it->second[i].ratUnco << "\t"
					<< it->second[i].hcc << "\t"
					<< it->second[i].efficiency << std::endl;

			allSaved &= stream.good();
		}
	}

	std::string rocDataName;
	for (std::map<std::string, 
		std::vector<ThreshInfo> >::const_iterator
		it = threshInf.begin() ; it != threshInf.end() ; ++it)
	{
		rocDataName = it->first + "_ROCData";
		if (saver.isSaving(rocDataName))
		{
			ofstream & stream = saver.getStream();
			this->AddSavedFile(rocDataName, saver.getCurrFile());

			stream << "FalsePositiveRatio\tTruePositiveRatio" << std::endl;
			for (unsigned int i = 0 ; i < it->second.size() ; ++i)
				stream << it->second[i].falsePosRat << "\t" 
					<< it->second[i].truePosRat << std::endl;

			allSaved &= stream.good();
		}
	}

	return allSaved;
}

//**********************************************************************
//**********************************************************************
void FunctionalTopoMetric::Initialize()
{
	for (std::map<std::string, Network<BooleanLink> *>::iterator it = 
		functTopos.begin() ; it != functTopos.end() ; ++it)
	{
		if (it->second)
			delete it->second;
	}
	functTopos.clear();
	threshInf.clear();
	chosenThresh.clear();
	optRocPoint.clear();
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//**********************************************************************
std::map<std::string, double> FunctionalTopoMetric::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	std::map<std::string, double> tmpScals;
	std::map<std::string, std::pair<double, double> >::const_iterator ROCit;
	std::map<std::string, unsigned int>::const_iterator CTit;
	std::map<std::string, std::vector<ThreshInfo> >::const_iterator TIit;
	for (std::map<std::string, Network<BooleanLink> *>::const_iterator it = 
		functTopos.begin() ; it != functTopos.end() ; ++it)
	{
		if (it->second)
		{
			std::vector<Metric*> functMetr = it->second->GetAllMetrics();
			for (unsigned int i = 0 ; i < functMetr.size() ; ++i)
			{
				assert(functMetr[i]);
				tmpScals = functMetr[i]->GetScalarStatsToSave();
				for (std::map<std::string, double>::const_iterator it2 =
						tmpScals.begin() ; it2 != tmpScals.end() ; ++it2)
					res[it->first + "_" + it2->first] = it2->second;
			}
			if ((ROCit = optRocPoint.find(it->first)) != optRocPoint.end())
			{
				res[it->first + "_FunctTopoOptROCThresh"] = 
					ROCit->second.first;
				res[it->first + "_FunctTopoOptROCDist"] = 
					ROCit->second.second;
			}
			if (((CTit = chosenThresh.find(it->first)) != chosenThresh.end()) and
				((TIit = threshInf.find(it->first)) != threshInf.end()))
			{
				res[it->first + "_ChosenThreshROC_FalsePos"] = 
					TIit->second[CTit->second].falsePosRat;
				res[it->first + "_ChosenThreshROC_TruePos"] = 
					TIit->second[CTit->second].truePosRat;
				res[it->first + "_ChosenThresh_ThreshVal"] = 
					TIit->second[CTit->second].thresh;
			}
		}
	}
	return res;
}

//**********************************************************************
//**********************************************************************
ParamHandler FunctionalTopoMetric::BuildModelParamHandler()
{
	ParamHandler params;

	for (std::map<std::string, Network<BooleanLink> *>::const_iterator it = 
		functTopos.begin() ; it != functTopos.end() ; ++it)
	{
		if (it->second)
		{
			params += it->second->BuildModelParamHandler().AddPrefix(
				it->first + "_");
		}
	}
	params <= "FunctTopoThresholdEstimMethod", threshEstimMethod;
	params <= "FunctTopoThresholdEstimMeanDeg", threshSpecifMeanDeg;
	params <= "FunctTopoThresholdEstimStdDevCoeff", stdDevCoeff;
	params <= "FunctTopoThresholdEstimModVoroMinDec", modVoroMinDecile;

	return params;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool FunctionalTopoMetric::LoadFromStream(std::ifstream & stream)
{
	stream >> threshEstimMethod;
	stream >> threshSpecifMeanDeg;
	stream >> stdDevCoeff;
	stream >> modVoroMinDecile;
	stream >> computeL;
	stream >> computeHCC;
	stream >> useMinForBidir;
	stream >> functTopoUseSpecmetr;

	unsigned int nbMetr;
	std::string strTmp;
	functTopoMetrList.clear();

	stream >> nbMetr;
	for (unsigned int i = 0 ; i < nbMetr ; ++i)
	{
		stream >> strTmp;
		functTopoMetrList.push_back(strTmp);
	}

	if (voroConstrStrat)
		delete voroConstrStrat;
	voroConstrStrat = new VoronoiConstructStrat(stream);

	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool FunctionalTopoMetric::SaveToStream(std::ofstream & stream) const
{
	stream << threshEstimMethod << std::endl;
	stream << threshSpecifMeanDeg << std::endl;
	stream << stdDevCoeff << std::endl;
	stream << modVoroMinDecile << std::endl;
	stream << computeL << std::endl;
	stream << computeHCC << std::endl;
	stream << useMinForBidir << std::endl;
	stream << functTopoUseSpecmetr << std::endl;
	stream << functTopoMetrList.size() << std::endl;
	for (unsigned int i = 0 ; i < functTopoMetrList.size() ; ++i)
		 stream << functTopoMetrList[i] << std::endl;
	voroConstrStrat->SaveToStream(stream);
	return stream.good();
}

//**********************************************************************
// Compute the functionnal topology using the Network class
//**********************************************************************
void FunctionalTopoMetric::ComputeFunctionalTopo(const ChIModel & model, 
	const std::string & name, const std::vector<std::vector<double> > & mat)
{
	unsigned int nbCells = model.GetNbCells();

	// Connec by thresh and ROC curve computations
	double optRocThresh = 0;
	double sameDegThresh = 0;
	double optRocDist = 0;
	std::vector<std::vector<double> > bidirMat = 
		std::vector<std::vector<double> >(nbCells, std::vector<double>(nbCells, 0));
	std::vector<std::vector<double> > tmpThreshMat = bidirMat;

	for (unsigned int i = 0 ; i < bidirMat.size() ; ++i)
		for (unsigned int j = i + 1 ; j < bidirMat[i].size() ; ++j)
		{
			if (useMinForBidir)
				bidirMat[i][j] = min(mat[i][j], mat[j][i]);
			else
				bidirMat[i][j] = max(mat[i][j], mat[j][i]);
			bidirMat[j][i] = bidirMat[i][j];
		}

	std::vector<double> sortBidirVals = GetSortedValuesFromMat(bidirMat);
	std::vector<std::vector<double> > distances;

	threshInf[name].clear();

	double onLinksNb = 0;
	double nbConnected;
	double meanPath;
	unsigned int maxLInd = 0;
	// Params HCC
	unsigned int hccStart = ParamHandler::GlobalParams.
		getParam<unsigned int>("-HCCParams", 0);
	unsigned int hccEnd = ParamHandler::GlobalParams.
		getParam<unsigned int>("-HCCParams", 1);
	unsigned int hccRepeat = ParamHandler::GlobalParams.
		getParam<unsigned int>("-HCCParams", 2);

	if ((threshEstimMethod == 5) or (threshEstimMethod == 7))
		assert(computeL);

	for (unsigned int i = 0 ; i < sortBidirVals.size() ; ++i)
		if ((i == 0) or (sortBidirVals[i] != sortBidirVals[i - 1]))
		{
			ApplyThresholdOnMat(bidirMat, tmpThreshMat, sortBidirVals[i]);
			threshInf[name].push_back(ThreshInfo());
			threshInf[name].back().thresh = sortBidirVals[i];

			// L by thresh
			if (computeL)
			{
				ComputeAllPairDistsFromBidirAdjMat(tmpThreshMat, distances);
				nbConnected = 0;
				meanPath = 0;
				for (unsigned int k = 0 ; k < distances.size() ; ++k)
					for (unsigned int m = k+1 ; m < distances[k].size() ; ++m)
						if (distances[k][m] < DEFAULT_MAX_PATH)
						{
							meanPath += distances[k][m];
							nbConnected++;
						}
				threshInf[name].back().meanPath = meanPath / nbConnected;
				threshInf[name].back().ratUnco = 1.0 - 2.0 * nbConnected / 
					((double)distances.size()*((double)distances.size()-1.0));
				if (threshInf[name][maxLInd].meanPath < 
						threshInf[name].back().meanPath)
				{
					maxLInd = threshInf[name].size() - 1;
				}
				threshInf[name].back().efficiency = ComputeEfficiency(distances);
			}
			// HCC by thresh
			if (computeHCC)
			{
				threshInf[name].back().hcc = ComputeHierarchClustCoeff(
					tmpThreshMat, hccStart, hccEnd, hccRepeat);
			}

			// ROC Data and Connec by thresh
			simplifROCData rocDat = computeROCPoint(model, tmpThreshMat);
			threshInf[name].back().connec = rocDat.connec;
			if (model.IsNetworkTopoDefined())
			{
				threshInf[name].back().falsePosRat = rocDat.falsePosRat;
				threshInf[name].back().truePosRat  = rocDat.truePosRat;
				if (rocDat.truePosRat * rocDat.onLinksNb - rocDat.falsePosRat * (nbCells * (nbCells - 1.0) 
					/ 2.0 - rocDat.onLinksNb) > optRocDist)
				{
					optRocDist = rocDat.truePosRat * rocDat.onLinksNb - rocDat.falsePosRat * 
						(nbCells * (nbCells - 1.0) / 2.0 - rocDat.onLinksNb);
					optRocThresh = threshInf[name].back().connec;
				}
			}
			onLinksNb = rocDat.onLinksNb;
		}

	sameDegThresh = onLinksNb / (nbCells * (nbCells - 1.0) / 2.0);
	if (model.IsNetworkTopoDefined())
		optRocPoint[name] = make_pair(optRocThresh, optRocDist);

	threshInf[name].push_back(ThreshInfo());
	threshInf[name].back().thresh = sortBidirVals.back()
		+ 1.0/((double)DEFAULT_MAX_VAL);
	threshInf[name].back().ratUnco = 1;

	bool useThr = false;
	// case 5 vars
	double ratOfMax = 0.7;
	unsigned int halfMaxLInd = 0;
	// case 6 vars
	const AbstractSpatialNetwork * net = 0;
	CopyStructureBuilder *build = 0;
	SpatialNetwork<BooleanLink> *tmpNet = 0;
	DegreeDistrComp *deg = 0;
	// case 7 vars
	AllPairDistances *apd;
	double actNetRat;
	double kStep = 0.5;
	double connecStep = kStep / ((double)nbCells - 1.0);
	unsigned int currInd;
	unsigned int maxDropInd;
	double maxDropVal;
	// case 8 vars
	std::vector<std::vector<double> > modVoroAdjMat;
	simplifROCData rocDat;

	switch (threshEstimMethod)
	{
		case 1:
			assert(model.IsNetworkTopoDefined());
			useThr = true;
			chosenThresh[name] = getIndFromConnec(name, sameDegThresh);
		break;
		case 2:
			useThr = true;
			chosenThresh[name] = getIndFromConnec(name, 
				threshSpecifMeanDeg / (nbCells - 1.0));
		break;
		case 3:
			assert(model.IsNetworkTopoDefined());
			useThr = true;
			chosenThresh[name] = getIndFromConnec(name, optRocThresh);
		break;
		case 4:
			useThr = false;
			chosenThresh[name] = 0;
		break;
		case 5:
			useThr = true;
			halfMaxLInd = 0;
			for (unsigned int i = 0 ; i < maxLInd ; ++i)
				if (std::abs(threshInf[name][i].meanPath - 
						threshInf[name][maxLInd].meanPath * ratOfMax) < 
						std::abs(threshInf[name][halfMaxLInd].meanPath - 
						threshInf[name][maxLInd].meanPath * ratOfMax))
					halfMaxLInd = i;
			chosenThresh[name] = halfMaxLInd;
		break;
		// Degree from Voronoi diagram
		case 6:
			net = dynamic_cast<const AbstractSpatialNetwork *>(&(model.GetNetwork()));

			build = new CopyStructureBuilder(*net);
			tmpNet = 
				new SpatialNetwork<BooleanLink>(net->size(), voroConstrStrat, BooleanLink::ClassName,
				build, net->GetDim(), false, true);
			deg = new DegreeDistrComp();
			tmpNet->AddMetric(deg, true);
			// Build voronoi diagram
			tmpNet->BuildNetwork();

			// Compute and get mean degree
			tmpNet->ComputeMetrics();
			chosenThresh[name] = getIndFromConnec(name, 
				deg->GetMeanDegree() / ((double)net->size() - 1.0));
			useThr = true;
			delete tmpNet;
		break;
		case 7:
			apd = GetSpecificMetric<Metric, AllPairDistances>(model.GetNetwork().GetAllMetrics());
			assert(apd);
			actNetRat = apd->GetRatioUnconnected();
			kStep = 0.5;
			connecStep = kStep / ((double)nbCells - 1.0);
			currInd = threshInf[name].size() - 1;
			maxDropInd = 0;
			maxDropVal = 0;
			for (int i = threshInf[name].size() - 1 ; i > 0 ; --i) 
				if (threshInf[name][i].connec >= (threshInf[name][currInd].connec + connecStep))
				{
					if ((threshInf[name][currInd].ratUnco <= (actNetRat + 0.01)) and 
						(maxDropVal > (threshInf[name][currInd].meanPath 
							- threshInf[name][i].meanPath)))
					{
						maxDropInd = currInd;
						maxDropVal = threshInf[name][currInd].meanPath 
							- threshInf[name][i].meanPath;
					}
					currInd = i;
				}
			chosenThresh[name] = maxDropInd;
		break;
		case 8:
			// Get voronoi topology
			net = dynamic_cast<const AbstractSpatialNetwork *>(&(model.GetNetwork()));
			build = new CopyStructureBuilder(*net);
			tmpNet = 
				new SpatialNetwork<BooleanLink>(net->size(), voroConstrStrat, BooleanLink::ClassName,
				build, net->GetDim(), false, true);
			// Build voronoi diagram
			tmpNet->BuildNetwork();
			modVoroAdjMat = tmpNet->GetAdjMat();
			// Get decile value
			double tmpThresh = threshInf[name][floor(modVoroMinDecile*(threshInf[name].size()-1))].thresh;
TRACE(threshInf[name][floor(modVoroMinDecile*(threshInf[name].size()-1))].connec)
TRACE(tmpThresh)
			// Remove all links in voro whose value is below decile value
			for (unsigned int i = 0 ; i < modVoroAdjMat.size() ; ++i)
				for (unsigned int j = i + 1 ; j < modVoroAdjMat[i].size() ; ++j)
					if (bidirMat[i][j] < tmpThresh)
					{
						modVoroAdjMat[i][j] = 0;
						modVoroAdjMat[j][i] = 0;
					}
			delete tmpNet;

			rocDat = computeROCPoint(model, modVoroAdjMat);
			threshInf[name].push_back(ThreshInfo());
			threshInf[name].back().thresh = -1;
			threshInf[name].back().falsePosRat = rocDat.falsePosRat;
			threshInf[name].back().truePosRat = rocDat.truePosRat;
			threshInf[name].back().connec = rocDat.connec;
			chosenThresh[name] = threshInf[name].size() - 1;
		break;
	}

	// Functionnal topology computations
	if ((functTopos.find(name) != functTopos.end()) and functTopos[name])
		delete functTopos[name];
	FunctionalTopoStrat *strat;
	// Specific behavior if using modified voronoi diagram
	if (threshEstimMethod == 8)
		strat = new FunctionalTopoStrat(modVoroAdjMat, false, false, 
			0.0, stdDevCoeff, false, true, 0.5);
	else
		strat = new FunctionalTopoStrat(mat, false, useThr, 
			threshInf[name][chosenThresh[name]].connec, stdDevCoeff);
	
	functTopos[name] = new Network<BooleanLink>(nbCells, strat,
		true, "BooleanLink");
	functTopos[name]->BuildNetwork();

	// Add metrics
	Metric *tmpMetrPtr;
	if (functTopoUseSpecmetr)
	{
		for (unsigned int i = 0 ; i < functTopoMetrList.size() ; ++i)
		{
			assert(AbstractFactory<Metric>::Factories[functTopoMetrList[i]]);
			tmpMetrPtr = AbstractFactory<Metric>::
				Factories[functTopoMetrList[i]]->Create();
			if (not functTopos[name]->AddMetric(tmpMetrPtr, true))
				delete tmpMetrPtr;
		}
	}
	else
	{
		std::vector<Metric*> tmpMetr = model.GetNetwork().GetAllMetrics();
		for (unsigned int i = 0 ; i < tmpMetr.size() ; ++i)
		{
			tmpMetrPtr = tmpMetr[i]->BuildCopy();
			if (not functTopos[name]->AddMetric(tmpMetrPtr, true))
				delete tmpMetrPtr;
		}
	}

	functTopos[name]->ComputeMetrics();
}

//**********************************************************************
//**********************************************************************
const AbstractNetwork * FunctionalTopoMetric::GetFunctTopo(
	const std::string & name) const
{
	std::map<std::string, Network<BooleanLink> *>::const_iterator it;
	if ((it = functTopos.find(name)) != functTopos.end())
		return it->second;
	else
		return 0;
}

//**********************************************************************
// Returns the indice of std::vector<ThreshInfo> corresponding
// to given connectivity
//**********************************************************************
unsigned int FunctionalTopoMetric::getIndFromConnec(
	const std::string & name, double connec) const
{
	unsigned int res = 0;
	std::map<std::string, std::vector<ThreshInfo> >::const_iterator it;
	if ((it = threshInf.find(name)) != threshInf.end())
	{
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			if (std::abs(it->second[i].connec - connec) < 
					std::abs(it->second[res].connec - connec))
				res = i;
		return res;
	}
	else
		return DEFAULT_MAX_VAL;
}

//**********************************************************************
// Computes the ROC score of an adj mat
//**********************************************************************
FunctionalTopoMetric::simplifROCData FunctionalTopoMetric::computeROCPoint(
	const ChIModel & model, const std::vector<std::vector<double> > &adjMat) const
{
	unsigned int nbCells = model.GetNbCells();
	simplifROCData ret;
	ret.falsePosRat = 0;
	ret.truePosRat = 0;
	ret.onLinksNb = 0;
	double nbLinks = 0;
	for (unsigned int j = 0 ; j < adjMat.size() ; ++j)
		for (unsigned int k = j + 1 ; k < adjMat[j].size() ; ++k)
		{
			if (adjMat[j][k] > 0)
				++nbLinks;

			if (model.IsNetworkTopoDefined())
			{
				if (model.GetNetwork().AreConnected(j,k))
				{
					++ret.onLinksNb;
					if (adjMat[j][k] > 0)
						++ret.truePosRat;
				}
				else if (adjMat[j][k] > 0)
					++ret.falsePosRat;
			}
		}

	if (model.IsNetworkTopoDefined())
	{
		ret.falsePosRat /= (nbCells * (nbCells - 1.0) / 2.0 - ret.onLinksNb);
		ret.truePosRat /= ret.onLinksNb;
	}

	ret.connec = nbLinks / (nbCells * (nbCells - 1.0) / 2.0);
	return ret;
}

//********************************************************************//
//*************** C O R R E L A T I O N   M E T R I C ****************//
//********************************************************************//

const std::string CorrelationsMetric::ZeroCorrTopoName =
	"Corr_ZeroLag_FunctTopo";
const std::string CorrelationsMetric::MaxCorrTopoName = 
	"Corr_maxCross_FunctTopo";
const std::string CorrelationsMetric::MaxCorrNoZeroTopoName =
	"Corr_maxCross_noZero_FunctTopo";

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(CORRELATIONS_METRIC, CONCENTRATIONS_METRIC)
ADD_METRIC_DEPENDENCY(CORRELATIONS_METRIC, DEGREE_DISTRIBUTION_NETWORK_METRIC)

//**********************************************************************
// Default Constructor
//**********************************************************************
CorrelationsMetric::CorrelationsMetric(ParamHandler & h) : 
	FunctionalTopoMetric(h)
{
	maxLag = h.getParam<double>("-CorrelParams", 0);
}

//**********************************************************************
// Full constructor
//**********************************************************************
CorrelationsMetric::CorrelationsMetric(int _m, double _t, double _sdC, 
	double _mvmd, bool _cl, bool _chcc, bool _umfb,	bool _uSM, 
	const std::vector<std::string> & _fTML, double _mL) : 
	FunctionalTopoMetric(_m, _t, _sdC, _mvmd, _cl, _chcc, _umfb, _uSM, _fTML), maxLag(_mL)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
CorrelationsMetric::CorrelationsMetric(ifstream & stream)
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Default Destructor
//**********************************************************************
CorrelationsMetric::~CorrelationsMetric()
{

}

//**********************************************************************
// Compute the metric
//**********************************************************************
bool CorrelationsMetric::ComputeMetric(const ChIModel & model)
{
	ConcentrationMetrics *concentrMetric = 
		GetSpecificMetric<Metric, ConcentrationMetrics>(model.GetAllMetrics());
	assert(concentrMetric);

	DegreeDistrComp *degreeDistrMetr = 
		GetSpecificMetric<Metric, DegreeDistrComp>(model.GetAllMetrics());
	assert(degreeDistrMetr);

	const std::vector<std::pair<double, std::vector<double> > > & concentrRaw =
		concentrMetric->GetConcentrations();
	const std::vector<std::string> & valNames = concentrMetric->GetValNames();
	// First, tranform concentr to a more usable format
	std::vector<std::vector<double> > concentrations(model.GetNbCells(), 
		std::vector<double>(concentrRaw.size(), 0));
	unsigned int cellInd = 0;
	for (unsigned int i = 0 ; i < valNames.size() ; ++i)
		if (valNames[i].find(std::string("Ca")) != std::string::npos)
		{
			for (unsigned int t = 0 ; t < concentrRaw.size() ; ++t)
				concentrations[cellInd][t] = concentrRaw[t].second[i];
			++cellInd;
		}

	unsigned int nbCells = concentrations.size();
	double savingStep = concentrMetric->GetSavingStep();
	int maxLagStep = maxLag / savingStep;

	assert(nbCells > 0);
	assert(maxLagStep > 0);

	zeroLagCorr = std::vector<std::vector<double> >(nbCells, std::vector<double>(nbCells, 0));
	maxCrossCorrVals = std::vector<std::vector<double> >(nbCells, std::vector<double>(nbCells, 0));
	maxCrossCorrLags = std::vector<std::vector<double> >(nbCells, std::vector<double>(nbCells, 0));
	maxCrossCorrValsNoZero = std::vector<std::vector<double> >(nbCells, std::vector<double>(nbCells, 0));
	maxCrossCorrLagsNoZero = std::vector<std::vector<double> >(nbCells, std::vector<double>(nbCells, 0));
	directionnality = std::vector<std::vector<double> >(nbCells, std::vector<double>(nbCells, 0));

	double totCorr = 0;
	double normFact1 = 0;
	double normFact2 = 0;
	double firstTierMax, lastTierMax, firstTierLag, lastTierLag;
	// For each pair of cells
	for (unsigned int i = 0 ; i < nbCells ; ++i)
		for (unsigned int j = i + 1 ; j < nbCells ; ++j)
		{
			firstTierMax = 0;
			lastTierMax = 0;
			firstTierLag = ((double)(2*maxLagStep + 1)) / 3.0 - maxLagStep;
			lastTierLag = ((double)(2*maxLagStep + 1)) * 2.0 / 3.0 - maxLagStep;
			// For each timelag
			for (int lagStep = -maxLagStep ; lagStep <= maxLagStep ; ++lagStep)
			{
				totCorr = 0;
				normFact1 = 0;
				normFact2 = 0;
				for (unsigned int tStep = ((lagStep < 0) ? -lagStep : 0) ; 
				    tStep < (concentrations[i].size() - ((lagStep > 0) ? lagStep : 0)) ; ++tStep)
				{
					totCorr += concentrations[i][tStep] * concentrations[j][tStep + lagStep];
					normFact1 += pow(concentrations[i][tStep], 2.0);
					normFact2 += pow(concentrations[j][tStep + lagStep], 2.0);
				}
				totCorr /= sqrt(normFact1 * normFact2);

				if (lagStep == 0)
				{
				    zeroLagCorr[i][j] = totCorr;
				    zeroLagCorr[j][i] = totCorr;
				}
				if (std::abs(maxCrossCorrVals[i][j]) < std::abs(totCorr))
				{
				    maxCrossCorrVals[i][j] = totCorr;
				    maxCrossCorrLags[i][j] = lagStep * savingStep;
				}
				if (std::abs(maxCrossCorrVals[j][i]) < std::abs(totCorr))
				{
				    maxCrossCorrVals[j][i] = totCorr;
				    maxCrossCorrLags[j][i] = - lagStep * savingStep;
				}
				if ((lagStep <= firstTierLag) or (lagStep >= lastTierLag))
				{
					if (std::abs(maxCrossCorrValsNoZero[i][j]) < std::abs(totCorr))
					{
						maxCrossCorrValsNoZero[i][j] = totCorr;
						maxCrossCorrLagsNoZero[i][j] = lagStep * savingStep;
					}
					if (std::abs(maxCrossCorrValsNoZero[j][i]) < std::abs(totCorr))
					{
						maxCrossCorrValsNoZero[j][i] = totCorr;
						maxCrossCorrLagsNoZero[j][i] = - lagStep * savingStep;
					}
				}
				// If we are in the first tier
				if (lagStep <= firstTierLag)
				{
					firstTierMax = std::max(firstTierMax, std::abs(totCorr));
				}
				// If we are in the last tier
				else if (lagStep >= lastTierLag)
				{
					lastTierMax = std::max(lastTierMax, std::abs(totCorr));
				}
			}
			// Filling directionnality data
			directionnality[i][j] = lastTierMax / firstTierMax;
			directionnality[j][i] = firstTierMax / lastTierMax;
		}

	// Compute common neighbors matrix
	nbCommonNeighbs = std::vector<std::vector<double> >(nbCells, 
		std::vector<double>(nbCells, 0));
	for (unsigned int i = 0 ; i < nbCells ; ++i)
		for (unsigned int j = i + 1 ; j < nbCells ; ++j)
		{
			for (unsigned int k = 0 ; k < nbCells ; ++k)
				if (model.GetNetwork().AreConnected(i, k) and 
						model.GetNetwork().AreConnected(j, k))
					nbCommonNeighbs[i][j] += 1;
			nbCommonNeighbs[j][i] = nbCommonNeighbs[i][j];
		}

	FunctionalTopoMetric::ComputeFunctionalTopo(model,
		ZeroCorrTopoName, zeroLagCorr);
	FunctionalTopoMetric::ComputeFunctionalTopo(model, 
		MaxCorrTopoName, maxCrossCorrVals);
	FunctionalTopoMetric::ComputeFunctionalTopo(model, 
		MaxCorrNoZeroTopoName, maxCrossCorrValsNoZero);

	return true;
}

//**********************************************************************
// Save the metric
//**********************************************************************
bool CorrelationsMetric::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	allSaved &= FunctionalTopoMetric::SaveMetric(saver);

	std::string zeroLagCorrelMatName("ZeroLagCorrelationMatrix");
	if (saver.isSaving(zeroLagCorrelMatName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(zeroLagCorrelMatName, saver.getCurrFile());

		allSaved &= SaveFullMatrix(stream, zeroLagCorr);
	}

	std::string maxCorrelValsMatName("MaxCorrelationValuesMatrix");
	if (saver.isSaving(maxCorrelValsMatName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(maxCorrelValsMatName, saver.getCurrFile());

		allSaved &= SaveFullMatrix(stream, maxCrossCorrVals);
	}

	std::string maxCorrelLagsMatName("MaxCorrelationLagsMatrix");
	if (saver.isSaving(maxCorrelLagsMatName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(maxCorrelLagsMatName, saver.getCurrFile());

		allSaved &= SaveFullMatrix(stream, maxCrossCorrLags);
	}

	std::string maxCorrelValsNoZeroMatName("MaxCorrelationValuesNoZeroMatrix");
	if (saver.isSaving(maxCorrelValsNoZeroMatName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(maxCorrelValsNoZeroMatName, saver.getCurrFile());

		allSaved &= SaveFullMatrix(stream, maxCrossCorrValsNoZero);
	}

	std::string maxCorrelLagsNoZeroMatName("MaxCorrelationNoZeroLagsMatrix");
	if (saver.isSaving(maxCorrelLagsNoZeroMatName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(maxCorrelLagsNoZeroMatName, saver.getCurrFile());

		allSaved &= SaveFullMatrix(stream, maxCrossCorrLagsNoZero);
	}

	std::string zeroCorrelCommonNeighbName("CommonNeighbsZeroCorrel");
	if (saver.isSaving(zeroCorrelCommonNeighbName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(zeroCorrelCommonNeighbName, saver.getCurrFile());

		stream << "NbCommonNeighbs\tZeroLagCorrel" << std::endl;
		for (unsigned int i = 0 ; i < nbCommonNeighbs.size() ; ++i)
			for (unsigned int j = i + 1 ; j < nbCommonNeighbs[i].size() ; ++j)
				stream 
					<< nbCommonNeighbs[i][j] << "\t"
					<< zeroLagCorr[i][j] << std::endl;

		allSaved &= true;
	}

	std::string directionalityMatName("CorrelDirectionalityMatrix");
	if (saver.isSaving(directionalityMatName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(directionalityMatName, saver.getCurrFile());

		allSaved &= SaveFullMatrix(stream, directionnality);
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void CorrelationsMetric::Initialize()
{
	FunctionalTopoMetric::Initialize();
	zeroLagCorr.clear();
	maxCrossCorrVals.clear();
	maxCrossCorrLags.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * CorrelationsMetric::BuildCopy() const
{
	return new CorrelationsMetric(this->threshEstimMethod, 
		this->threshSpecifMeanDeg, this->stdDevCoeff, this->modVoroMinDecile, 
		this->computeL, this->computeHCC, this->useMinForBidir, this->functTopoUseSpecmetr, 
		this->functTopoMetrList, maxLag); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> CorrelationsMetric::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	res = FunctionalTopoMetric::GetScalarStatsToSave();
	return res;
}

//**********************************************************************
// Builds a parameter handling object
//**********************************************************************
ParamHandler CorrelationsMetric::BuildModelParamHandler()
{
	ParamHandler params;
	params += FunctionalTopoMetric::BuildModelParamHandler();
	return params;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool CorrelationsMetric::LoadFromStream(std::ifstream & stream)
{
	bool ok = FunctionalTopoMetric::LoadFromStream(stream);
	stream >> maxLag;

	return ok and stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool CorrelationsMetric::SaveToStream(std::ofstream & stream) const
{
	bool ok = FunctionalTopoMetric::SaveToStream(stream);
	stream << maxLag << std::endl;

	return ok and stream.good();
}

//**********************************************************************
// Get threshold to match connectivity
//**********************************************************************
double CorrelationsMetric::findThreholdForConnectivity(
	const std::vector<std::vector<double> > & corr, unsigned int nbEdges) const
{
	std::vector<double> corrVals = GetSortedValuesFromMat(corr);
	// Find the approximate value of correlation under which
	// all correlations are thresholded such that the mean degree
	// is conserved
	double threshCorr = 0;
	for (unsigned int i = 0 ; (i < corrVals.size()) and
			(i <= 2*nbEdges) ; ++i)
		threshCorr = corrVals[corrVals.size() - i - 1];
	
	return threshCorr;
}

//**********************************************************************
// Return the mean of the absolute values of lags for a given cross
// correlation topology
//**********************************************************************
double CorrelationsMetric::getMeanAbsLag(
	const AbstractNetwork * topo, 
	const std::vector<std::vector<double> > & lags) const
{
	double tempVal = 0;
	double nbVals = 0;
	assert(topo->size() == lags.size());
	for (unsigned int i = 0 ; i < topo->size() ; ++i)
		for (unsigned int j = 0 ; j < topo->size(i) ; ++j)
			if (topo->AreConnected(i, j))
			{
				tempVal += std::abs(lags[i][j]);
				++nbVals;
			}
	
	return tempVal / nbVals;
}

//********************************************************************//
//********** T R A N S F E R   E N T R O P Y   M E T R I C ***********//
//********************************************************************//

const std::string TransferEntropyMetric::TransEntrFunctTopoName = "TE_FunctTopo";

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(TRANSFER_ENTROPY_METRIC, CONCENTRATIONS_METRIC)

//**********************************************************************
// Default Constructor
//**********************************************************************
TransferEntropyMetric::TransferEntropyMetric(ParamHandler & h) : FunctionalTopoMetric(h)
{
	timeEmbed = h.getParam<double>("-TransferEntropy", 0);
	nbBins = h.getParam<unsigned int>("-TransferEntropy", 1);
	nbCarBins = Stringify(nbBins).length();
}

//**********************************************************************
// Full constructor
//**********************************************************************
TransferEntropyMetric::TransferEntropyMetric(int _m, double _t, 
	double _sdC, double _mvmd, bool _cl, bool _chcc, bool _umfb, bool _uSM, 
	const std::vector<std::string> & _fTML, double _te, unsigned int _nb) :
	FunctionalTopoMetric(_m, _t, _sdC, _mvmd, _cl, _chcc, _umfb, _uSM, _fTML), timeEmbed(_te),
	nbBins(_nb)
{
	nbCarBins = Stringify(nbBins).length();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
TransferEntropyMetric::TransferEntropyMetric(ifstream & stream)
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Default Destructor
//**********************************************************************
TransferEntropyMetric::~TransferEntropyMetric()
{
}

//**********************************************************************
// Compute the metric
//**********************************************************************
bool TransferEntropyMetric::ComputeMetric(const ChIModel & model)
{
	ConcentrationMetrics *concentrMetric = GetSpecificMetric<Metric, ConcentrationMetrics>(model.GetAllMetrics());
	assert(concentrMetric);

	const std::vector<std::pair<double, std::vector<double> > > & concentrRaw =
		concentrMetric->GetConcentrations();
	const std::vector<std::string> & valNames = concentrMetric->GetValNames();
	// First, tranform concentr to a more usable format
	std::vector<std::vector<double> > concentrations(model.GetNbCells(), 
		std::vector<double>(concentrRaw.size(), 0));
	unsigned int cellInd = 0;
	for (unsigned int i = 0 ; i < valNames.size() ; ++i)
		if (valNames[i].find(std::string("Ca")) != std::string::npos)
		{
			for (unsigned int t = 0 ; t < concentrRaw.size() ; ++t)
				concentrations[cellInd][t] = concentrRaw[t].second[i];
			++cellInd;
		}

	unsigned int nbCells = concentrations.size();
	unsigned int nbSteps = concentrations[0].size();
	double savingStep = concentrMetric->GetSavingStep();
	unsigned int k = ceil(timeEmbed / savingStep);

	nbCarBins = Stringify(nbBins).length();

	assert(nbCells > 0);
	assert(k > 0);

TRACE("Displaying parameters : NbBins " << nbBins << " // timeEmbed " << timeEmbed)
	computeSignalBins(concentrations);
TRACE("Bins computed")

	// Compute joint distributions for each cell
	signalJointDistrib_kp1 = 
		std::vector<std::vector<double> >(nbCells, std::vector<double>());
	signalJointDistrib_k = 
		std::vector<std::vector<double> >(nbCells, std::vector<double>());

	discrSig = std::vector<std::vector<std::pair<int, int> > >(
		nbSteps, std::vector<std::pair<int, int> >(nbCells, 
		std::make_pair<int, int>(-1, -1)));

	double totNb;
	std::string tempAddr_k, tempAddr_kp1;
	unsigned int ind_k, ind_kp1;
	unsigned int maxSigTabSize_k = 0;
	unsigned int maxSigTabSize_kp1 = 0;

	// For each cell
	for (unsigned int i = 0 ; i < nbCells ; ++i)
	{
		signalJointDistrib_k[i] = std::vector<double>(maxSigTabSize_k, 0);
		signalJointDistrib_kp1[i] = std::vector<double>(maxSigTabSize_kp1, 0);

		totNb = 0;
		// Span the whole time range to create the histogram
		for (unsigned int t = k - 1 ; t < nbSteps ; ++t)
		{
			tempAddr_kp1 = createSigStr(k + 1);
			tempAddr_k = createSigStr(k);
			// Get the corresponding address for the x_t^{k+1} vector
			if (t >= k)
				writeSigToStr(tempAddr_kp1, 0, concentrations[i][t - k]);
			for (unsigned int n = t - (k - 1) ; n <= t ; ++n)
			{
			    writeSigToStr(tempAddr_kp1, n - (t - k), concentrations[i][n]);
			    writeSigToStr(tempAddr_k, n - (t - (k - 1)), concentrations[i][n]);
			}

			ind_k = findInd(sigTable_k, tempAddr_k);
			if (ind_k == sigTable_k.size())
			{
				sigTable_k.push_back(tempAddr_k);
				signalJointDistrib_k[i].push_back(0);
			}
			// update nb of occurence
			signalJointDistrib_k[i][ind_k] += 1;
			// discr signal
			discrSig[t][i].first = ind_k;

			if (t >= k)
			{
				ind_kp1 = findInd(sigTable_kp1, tempAddr_kp1);
				if (ind_kp1 == sigTable_kp1.size())
				{
					sigTable_kp1.push_back(tempAddr_kp1);
					signalJointDistrib_kp1[i].push_back(0);
				}
				// update nb of occurence
				signalJointDistrib_kp1[i][ind_kp1] += 1;
				// discr signal
				discrSig[t][i].second = ind_kp1;
			}

			++totNb;
		}
		// Normalize the histogram
		for (unsigned int j = 0 ; j < signalJointDistrib_k[i].size() ; ++j)
		{
			signalJointDistrib_k[i][j] /= totNb;
		}
		for (unsigned int j = 0 ; j < signalJointDistrib_kp1[i].size() ; ++j)
		{
			signalJointDistrib_kp1[i][j] /= (totNb - 1.0);
		}

		// Update sig table sizes
		if (maxSigTabSize_k < signalJointDistrib_k[i].size())
			maxSigTabSize_k = signalJointDistrib_k[i].size();
		if (maxSigTabSize_kp1 < signalJointDistrib_kp1[i].size())
			maxSigTabSize_kp1 = signalJointDistrib_kp1[i].size();
	}

	// Compute transfer entropy values
	transferEntropy = std::vector<std::vector<double> >(nbCells, std::vector<double>(nbCells, 0));

	double jointProb;     // p(x_{t+1},x_t^k,y_{t+1}^k)
	double jointCondProb; // p(x_{t+1}|x_t^k,y_{t+1}^k)
	double autoCondProb;  // p(x_{t+1}|x_t^k)
	double pxtk;          // p(x_t^k)
	double pxtp1kp1;      // p(x_{t+1}^{k+1})
	double pxtkytk;       // p(x_t^k,y_{t+1}^k)

	sigPairJointDistrib_xkp1tp1_yktp1 = 
		std::vector<std::vector<std::vector<double> > >(nbCells, 
			std::vector<std::vector<double> >(nbCells, std::vector<double>()));
	sigPairJointDistrib_xkt_yktp1 = 
		std::vector<std::vector<std::vector<double> > >(nbCells, 
			std::vector<std::vector<double> >(nbCells, std::vector<double>()));

	std::pair<int, int> tempAddr_xkp1tp1_yktp1;
	std::pair<int, int> tempAddr_xkt_yktp1;

	unsigned int ind_xkt_yktp1;
	unsigned int ind_xkp1tp1_yktp1;
	unsigned int maxSigPairTabSize_xkt_yktp1 = 0;
	unsigned int maxSigPairTabSize_xkp1tp1_yktp1 = 0;
	std::vector<unsigned int> diffValsTaken;

	for (unsigned int i = 0 ; i < nbCells ; ++i) // Y
	{
		for (unsigned int j = 0 ; j < nbCells ; ++j) // X
		{
			if (i != j)
			{
				sigPairJointDistrib_xkt_yktp1[i][j] = 
					std::vector<double>(maxSigPairTabSize_xkt_yktp1, 0);
				sigPairJointDistrib_xkp1tp1_yktp1[i][j] = 
					std::vector<double>(maxSigPairTabSize_xkp1tp1_yktp1, 0);

				tmpDiscrPairSig = std::vector<std::pair<int, int> >(nbSteps, 
					std::make_pair<int, int>(-1, -1));

				totNb = 0;
				diffValsTaken.clear();

				// Compute full joint distrib between cells i and j
				// Span the whole time range to create the histogram
				for (unsigned int t = k - 1 ; t < (nbSteps - 1) ; ++t)
				{
					tempAddr_xkt_yktp1.first = discrSig[t][j].first;
					//
					tempAddr_xkt_yktp1.second = discrSig[t][i].first;
					//
					if (t >= k)
					{
						tempAddr_xkp1tp1_yktp1.first = discrSig[t+1][j].second;
						//
						tempAddr_xkp1tp1_yktp1.second = discrSig[t][i].first;
						//
					}

					ind_xkt_yktp1 = findInd(sigJointTable_xkt_yktp1, tempAddr_xkt_yktp1);
					if (ind_xkt_yktp1 == sigJointTable_xkt_yktp1.size())
					{
						sigJointTable_xkt_yktp1.push_back(tempAddr_xkt_yktp1);
						sigPairJointDistrib_xkt_yktp1[i][j].push_back(0);
					}
					// update nb of occurence
					sigPairJointDistrib_xkt_yktp1[i][j][ind_xkt_yktp1] += 1.0;
					// discr signal
					tmpDiscrPairSig[t].first = ind_xkt_yktp1;

					if (t >= k)
					{
						ind_xkp1tp1_yktp1 = findInd(sigJointTable_xkp1tp1_yktp1, 
							tempAddr_xkp1tp1_yktp1);
						if (ind_xkp1tp1_yktp1 == sigJointTable_xkp1tp1_yktp1.size())
						{
							sigJointTable_xkp1tp1_yktp1.push_back(tempAddr_xkp1tp1_yktp1);
							sigPairJointDistrib_xkp1tp1_yktp1[i][j].push_back(0);
						}
						if (sigPairJointDistrib_xkp1tp1_yktp1[i][j][ind_xkp1tp1_yktp1] == 0)
							diffValsTaken.push_back(t);
						//
						// update nb of occurence
						sigPairJointDistrib_xkp1tp1_yktp1[i][j][ind_xkp1tp1_yktp1] += 1.0;
						// discr signal
						tmpDiscrPairSig[t].second = ind_xkp1tp1_yktp1;
					}

					++totNb;
				}
				// Normalize the histogram
				for (unsigned int m = 0 ; m < sigPairJointDistrib_xkt_yktp1[i][j].size() ; ++m)
					sigPairJointDistrib_xkt_yktp1[i][j][m] /= totNb;
				for (unsigned int m = 0 ; m < sigPairJointDistrib_xkp1tp1_yktp1[i][j].size() ; ++m)
					sigPairJointDistrib_xkp1tp1_yktp1[i][j][m] /= (totNb - 1.0);

				// update max sig pair tab sizes
				if (maxSigPairTabSize_xkt_yktp1 < sigPairJointDistrib_xkt_yktp1[i][j].size())
					maxSigPairTabSize_xkt_yktp1 = sigPairJointDistrib_xkt_yktp1[i][j].size();
				if (maxSigPairTabSize_xkp1tp1_yktp1 < sigPairJointDistrib_xkp1tp1_yktp1[i][j].size())
					maxSigPairTabSize_xkp1tp1_yktp1 = sigPairJointDistrib_xkp1tp1_yktp1[i][j].size();

				// Span the whole time range
				unsigned int t;
				for (unsigned int tmpInd = 0 ; tmpInd < diffValsTaken.size() ; ++tmpInd)
				{
					t = diffValsTaken[tmpInd];

					// p(x_{t+1},x_t^k,y_{t+1}^k)
					jointProb = sigPairJointDistrib_xkp1tp1_yktp1[i][j]
						[tmpDiscrPairSig[t].second];

					// p(x_t^k,y_{t+1}^k)
					pxtkytk = sigPairJointDistrib_xkt_yktp1[i][j][tmpDiscrPairSig[t].first];

					// p(x_{t+1}|x_t^k,y_{t+1}^k)
					jointCondProb = min(jointProb / pxtkytk, 1.0);

					// p(x_{t+1}^{k+1})
					pxtp1kp1 = signalJointDistrib_kp1[j][discrSig[t + 1][j].second];

					// p(x_t^k)
					pxtk = signalJointDistrib_k[j][discrSig[t][j].first];

					// p(x_{t+1}|x_t^k)
					autoCondProb = min(pxtp1kp1 / pxtk, 1.0);

					// Add the partial entropy to the transfer entropy
					transferEntropy[i][j] += 
						jointProb * log(jointCondProb / autoCondProb);
				}
			}
		}
	}

	FunctionalTopoMetric::ComputeFunctionalTopo(model, TransEntrFunctTopoName, transferEntropy);

	return true;
}

//**********************************************************************
// Save the metric
//**********************************************************************
bool TransferEntropyMetric::SaveMetric(ResultSaver saver) const
{

	bool allSaved = true;

	allSaved &= FunctionalTopoMetric::SaveMetric(saver);

	std::string transferEntropyFullMatName("TransferEntropyFullMat");
	if (saver.isSaving(transferEntropyFullMatName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(transferEntropyFullMatName, saver.getCurrFile());

		allSaved &= SaveFullMatrix(stream, transferEntropy);
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void TransferEntropyMetric::Initialize()
{
	FunctionalTopoMetric::Initialize();
	binsLimits.clear();
	discrSig.clear();
	tmpDiscrPairSig.clear();
	signalJointDistrib_kp1.clear();
	signalJointDistrib_k.clear();
	sigTable_kp1.clear();
	sigTable_k.clear();
	sigPairJointDistrib_xkp1tp1_yktp1.clear();
	sigPairJointDistrib_xkt_yktp1.clear();
	sigJointTable_xkp1tp1_yktp1.clear();
	sigJointTable_xkt_yktp1.clear();

	transferEntropy.clear();
	nbCarBins = Stringify(nbBins).length();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * TransferEntropyMetric::BuildCopy() const
{
	return new TransferEntropyMetric(this->threshEstimMethod, 
		this->threshSpecifMeanDeg, this->stdDevCoeff, this->modVoroMinDecile,  
		this->computeL, this->computeHCC, this->useMinForBidir, this->functTopoUseSpecmetr, 
		this->functTopoMetrList, timeEmbed, nbBins); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> TransferEntropyMetric::GetScalarStatsToSave() const
{
	return FunctionalTopoMetric::GetScalarStatsToSave();
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool TransferEntropyMetric::LoadFromStream(std::ifstream & stream)
{
	bool ok = FunctionalTopoMetric::LoadFromStream(stream);
	stream >> timeEmbed;
	stream >> nbBins;

	return ok and stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool TransferEntropyMetric::SaveToStream(std::ofstream & stream) const
{
	bool ok = FunctionalTopoMetric::SaveToStream(stream);
	stream << timeEmbed << std::endl
		<< nbBins << std::endl;

	return ok and stream.good();
}


//**********************************************************************
// Computes the values necessary to bin raw signals
//**********************************************************************
void TransferEntropyMetric::computeSignalBins(
	const std::vector<std::vector<double> > & vals)
{
	binsLimits.clear();
	// Compute min and max values for the signal
	double minVal = vals[0][0];
	double maxVal = vals[0][0];
	for (unsigned int i = 0 ; i < vals.size() ; ++i)
		for (unsigned int t = 0 ; t < vals[i].size() ; ++t)
		{
			minVal = std::min(minVal, vals[i][t]);
			maxVal = std::max(maxVal, vals[i][t]);
		}
	// Divide in nbBins bins of same length
	for (double r = 0 ; r < nbBins ; ++r)
	{
		binsLimits.push_back(
			minVal + (maxVal - minVal)*((r + 1.0) / (double)nbBins));
	}
}

//**********************************************************************
// Returns the bin for the given value
//**********************************************************************
unsigned int TransferEntropyMetric::binSignalValue(double val) const
{
	unsigned int bin = 0;
	for (; (bin < binsLimits.size()) and (val > binsLimits[bin]) ; ++bin);
	return bin;
}

//**********************************************************************
//**********************************************************************
void TransferEntropyMetric::writeSigToStr(std::string & dest, 
    unsigned int ind, double val) const
{
    writeSigToStr(dest, ind, binSignalValue(val));
}

//**********************************************************************
//**********************************************************************
void TransferEntropyMetric::writeSigToStr(std::string & dest, 
    unsigned int ind, unsigned int bin) const
{
    std::string str(StringifyFixed(bin, nbCarBins));
    for (unsigned int i = 0 ; i < str.length() ; ++i)
		dest[ind * nbCarBins + i] = str[i];
}

//**********************************************************************
//**********************************************************************
const std::vector<std::vector<double> > TransferEntropyMetric::GetValues() const
{
	return transferEntropy;
}

//********************************************************************//
//********** R E C O N S T R   B Y   N E T   C O N S T R *************//
//********************************************************************//

const std::string FunctTopoByNetConstrStrat::FunctTopoByNetName =
	"FTopo_FromNet";

//**********************************************************************
// Dependencies
//**********************************************************************

//**********************************************************************
// Default Constructor
//**********************************************************************
FunctTopoByNetConstrStrat::FunctTopoByNetConstrStrat(ParamHandler & h) : 
	FunctionalTopoMetric(h)
{
	std::vector<std::string> stratNames = 
		h.getParam<std::vector<std::string> >("-FunctTopoByNetConstr", 0);

	for (unsigned int i = 0 ; i < stratNames.size() ; ++i)
	{
		assert(AbstractFactory<NetworkConstructStrat>::Factories[stratNames[i]]);
		netConstrStrats.push_back(AbstractFactory<NetworkConstructStrat>::
			Factories[stratNames[i]]->Create());
	}
}

//**********************************************************************
// Full constructor
//**********************************************************************
FunctTopoByNetConstrStrat::FunctTopoByNetConstrStrat(int _m, double _t, 
	double _sdC, double _mvmd, bool _cl, bool _chcc, bool _umfb, bool _uSM, 
	const std::vector<std::string> & _fTML, std::vector<NetworkConstructStrat*> _ncs) : 
	FunctionalTopoMetric(_m, _t, _sdC, _mvmd, _cl, _chcc, _umfb, _uSM, _fTML), netConstrStrats(_ncs)
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
FunctTopoByNetConstrStrat::FunctTopoByNetConstrStrat(ifstream & stream)
{
	LoadFromStream(stream); 
}

//**********************************************************************
// Default Destructor
//**********************************************************************
FunctTopoByNetConstrStrat::~FunctTopoByNetConstrStrat()
{
	for (unsigned int i = 0 ; i < netConstrStrats.size() ; ++i)
		delete netConstrStrats[i];
}

//**********************************************************************
// Compute the metric
//**********************************************************************
bool FunctTopoByNetConstrStrat::ComputeMetric(const ChIModel & model)
{
	const AbstractSpatialNetwork * net = 0;
	CopyStructureBuilder *build = 0;
	SpatialNetwork<BooleanLink> *tmpNet = 0;

	net = dynamic_cast<const AbstractSpatialNetwork *>(&(model.GetNetwork()));
	build = new CopyStructureBuilder(*net);
	for (unsigned int i = 0 ; i < netConstrStrats.size() ; ++i)
	{
		tmpNet = 
			new SpatialNetwork<BooleanLink>(net->size(), netConstrStrats[i], 
				BooleanLink::ClassName, build, net->GetDim(), false, false);
		// Build network from constr strat
		tmpNet->BuildNetwork();

		// Get Adj Mat
		adjMats.push_back(tmpNet->GetAdjMat());
		FunctionalTopoMetric::ComputeFunctionalTopo(model,
			FunctTopoByNetName + "_" + StringifyFixed(i, 3), 
			adjMats.back());

		delete tmpNet;
	}
	delete build;

	return true;
}

//**********************************************************************
// Save the metric
//**********************************************************************
bool FunctTopoByNetConstrStrat::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	allSaved &= FunctionalTopoMetric::SaveMetric(saver);

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void FunctTopoByNetConstrStrat::Initialize()
{
	FunctionalTopoMetric::Initialize();

	adjMats.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * FunctTopoByNetConstrStrat::BuildCopy() const
{
	std::vector<NetworkConstructStrat *> ncsDupl;
	for (unsigned int i = 0 ; i < netConstrStrats.size() ; ++i)
		ncsDupl.push_back(netConstrStrats[i]->BuildCopy());

	return new FunctTopoByNetConstrStrat(this->threshEstimMethod, 
		this->threshSpecifMeanDeg, this->stdDevCoeff, this->modVoroMinDecile,  
		this->computeL, this->computeHCC, this->useMinForBidir, this->functTopoUseSpecmetr, 
		this->functTopoMetrList, ncsDupl); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> FunctTopoByNetConstrStrat::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	res = FunctionalTopoMetric::GetScalarStatsToSave();
	return res;
}

//**********************************************************************
// Builds a parameter handling object
//**********************************************************************
ParamHandler FunctTopoByNetConstrStrat::BuildModelParamHandler()
{
	ParamHandler params;
	params += FunctionalTopoMetric::BuildModelParamHandler();
	return params;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool FunctTopoByNetConstrStrat::LoadFromStream(std::ifstream & stream)
{
	bool ok = FunctionalTopoMetric::LoadFromStream(stream);

	unsigned int nbStrats;
	std::string netConstrName;
	stream >> nbStrats;
	for (unsigned int i = 0 ; i < nbStrats ; ++i)
	{
		stream >> netConstrName;
		assert(AbstractFactory<NetworkConstructStrat>::Factories[netConstrName]);
		netConstrStrats.push_back(AbstractFactory<NetworkConstructStrat>::
			Factories[netConstrName]->CreateFromStream(stream));
	}

	return ok and stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool FunctTopoByNetConstrStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = FunctionalTopoMetric::SaveToStream(stream);

	stream << netConstrStrats.size() << std::endl;
	for (unsigned int i = 0 ; ok and (i < netConstrStrats.size()) ; ++i)
	{
		stream << netConstrStrats[i]->GetClassName() << std::endl;
		ok &= netConstrStrats[i]->SaveToStream(stream);
	}

	return ok and stream.good();
}

//********************************************************************//
//*********** T H R E S H O L D   D E T E R M I N A T I O N **********//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(THRESHOLD_DETERMINATION_METRIC, ACTIVATED_CELLS_METRICS)

//**********************************************************************
// Default Constructor
//**********************************************************************
ThresholdDetermination::ThresholdDetermination() 
{

}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ThresholdDetermination::ThresholdDetermination(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Compute the activated cells
//**********************************************************************
bool ThresholdDetermination::ComputeMetric(const ChIModel & model)
{
	const ChIModelThresholdDetermination & threshDetModel = 
		dynamic_cast<const ChIModelThresholdDetermination &>(model);

	ActivatedCells *activCells = 
		GetSpecificMetric<Metric, ActivatedCells>(model.GetAllMetrics());
	const std::vector<TimedActivations> & activ = activCells->GetActivations();

	// If a new simulation has begun
	if (model.GetTime() < lastCompTime)
	{
		double nbStim = 0;
		for (unsigned int i = 0 ; i < model.GetNbCells() ; ++i)
			if (model.IsStimulated(i))
				++nbStim;
		spiked.push_back(false);
		timeSpiked.push_back(-1);
		nbStims.push_back(nbStim);
		nbNeighbs.push_back(threshDetModel.GetDegree());
		nbDerivs.push_back(threshDetModel.GetNbDerivs());
	}
	// If reference cell has just spiked
	if (not activ.empty() and not spiked.back() and 
		(std::find(activ.back().cells.begin(), 
		activ.back().cells.end(), 0) != activ.back().cells.end()))
	{
		spiked.back() = true;
		timeSpiked.back() = model.GetTime();
	}

	lastCompTime = model.GetTime();
	
	return true;
}

//**********************************************************************
// Save the activated cells
//**********************************************************************
bool ThresholdDetermination::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string ThreshDetermName("ThresholdDetermination");
	if (saver.isSaving(ThreshDetermName) and not spiked.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(ThreshDetermName, saver.getCurrFile());

		stream << "HasSpiked\tSpikeTime\tThreshold\tSecondNeighbThresh" << endl;
		unsigned int ind;
		for (ind = 0 ; not spiked[ind] and (ind < (spiked.size() - 1)) ; ++ind);

		stream
			<< spiked[ind] << "\t"
			<< timeSpiked[ind] << "\t"
			<< (spiked[ind] ? 
				(nbStims[ind] / nbNeighbs[ind]) : DEFAULT_MAX_VAL) << "\t"
			<< (spiked[ind] ? 
				(nbStims[ind] / (nbDerivs[ind] + 1.0)) : DEFAULT_MAX_VAL)
			<< std::endl;

		allSaved &= stream.good();
	}
	std::string ThreshDetermFullName("ThresholdDeterminationFull");
	if (saver.isSaving(ThreshDetermFullName) and not spiked.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(ThreshDetermFullName, saver.getCurrFile());

		stream << "HasSpiked\tSpikeTime\tnbStims\tnbNeighbs" << endl;
		for (unsigned int j = 0 ; j < spiked.size() ; ++j)
			stream
				<< spiked[j] << "\t"
				<< timeSpiked[j] << "\t"
				<< nbStims[j] << "\t"
				<< nbNeighbs[j] << std::endl;

		allSaved &= stream.good();
	}

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void ThresholdDetermination::Initialize()
{
	Metric::Initialize();
	spiked.clear();
	timeSpiked.clear();
	nbStims.clear();
	nbNeighbs.clear();
	nbDerivs.clear();
	lastCompTime = DEFAULT_MAX_VAL;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * ThresholdDetermination::BuildCopy() const
{
	return new ThresholdDetermination(); 
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool ThresholdDetermination::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool ThresholdDetermination::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

//**********************************************************************
// Get the minimum thresh value for which the central cell has spiked
//**********************************************************************
double ThresholdDetermination::GetSecondNeighbThresh() const
{
	unsigned int ind;
	for (ind = 0 ; not spiked[ind] and (ind < (spiked.size() - 1)) ; ++ind);

	return (spiked[ind] ? (nbStims[ind] / (nbDerivs[ind] + 1.0)) : DEFAULT_MAX_VAL);
}

//**********************************************************************
//**********************************************************************
double ThresholdDetermination::GetSecondNeighbThreshUnderTime(double time) const
{
	unsigned int ind;
	for (ind = 0 ; (not spiked[ind] or timeSpiked[ind] > time) and (ind < (spiked.size() - 1)) ; ++ind);

	return ((spiked[ind] and timeSpiked[ind] <= time) ? 
		(nbStims[ind] / (nbDerivs[ind] + 1.0)) : DEFAULT_MAX_VAL);
}

//********************************************************************//
//**************** W A V E   F R O N T   M E T R I C *****************//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(WAVE_FRONT_DETECT_METRIC, PROPAGATION_DISTANCE_METRIC)
ADD_METRIC_DEPENDENCY(WAVE_FRONT_DETECT_METRIC, ALL_PAIR_DISTANCES_NETWORK_METRIC)
ADD_METRIC_DEPENDENCY(WAVE_FRONT_DETECT_METRIC, DEGREE_DISTRIBUTION_NETWORK_METRIC)
ADD_METRIC_DEPENDENCY(WAVE_FRONT_DETECT_METRIC, ACTIVATED_CELLS_METRICS)

//**********************************************************************
// Default Constructor
//**********************************************************************
WaveFrontDetect::WaveFrontDetect(ParamHandler & h)
{
	useParentsOfParentsAsNonSinks = h.getParam<bool>("-WaveDetectUPOPANS", 0);
}

//**********************************************************************
// Default Constructor
//**********************************************************************
WaveFrontDetect::WaveFrontDetect(bool _upopans) : 
	useParentsOfParentsAsNonSinks(_upopans) 
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
WaveFrontDetect::WaveFrontDetect(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Compute the activated cells
//**********************************************************************
bool WaveFrontDetect::ComputeMetric(const ChIModel & model)
{
	// Obtaining needed metrics
	PropagationDistance *activCells = 
		GetSpecificMetric<Metric, PropagationDistance>(model.GetAllMetrics());
	AllPairDistances *distMetr = 
		GetSpecificMetric<Metric, AllPairDistances>(model.GetAllMetrics());
	DegreeDistrComp  & degreeMetr = 
		*GetSpecificMetric<Metric, DegreeDistrComp>(model.GetAllMetrics());
	ActivatedCells *actCells = 
		GetSpecificMetric<Metric, ActivatedCells>(model.GetAllMetrics());

	double maxDelay = activCells->GetMaxDelay();
	
	// Initializing saved values
	frontiers = std::vector<std::vector<NodeValues> >(
		activCells->size(), std::vector<NodeValues>());
	nonActNeighbs = std::vector<std::vector<NodeValues> >(
		activCells->size(), std::vector<NodeValues>());

	// Temp variables
	std::set<unsigned int> invCellOnFront, nonActTemp;
	bool isOnFront;


	// For each wave
	for (unsigned int w = 0 ; w < activCells->size() ; ++w)
	{
		// All involved cells are potentially on the front
		invCellOnFront = std::set<unsigned int>(
			(*activCells)[w].involvedCells.begin(), 
			(*activCells)[w].involvedCells.end());
		unsigned int init = (*activCells)[w].initiator;

		// Fill nodeNbInShell
		std::vector<double> tmpShell;
		for (unsigned int i = 0 ; i < model.GetNbCells() ; ++i)
		{
			while (tmpShell.size() <= (*distMetr)[i][init])
				tmpShell.push_back(0);

			tmpShell[round((*distMetr)[i][init])] += 1.0;
			// Fill maxInfluxByShell
			while (maxInfluxByShell.size() <= (*distMetr)[i][init])
				maxInfluxByShell.push_back(std::vector<double>());

			maxInfluxByShell[round((*distMetr)[i][init])].push_back(
				actCells->GetMaxInfluxInTimeWin(i));

		}
		while (nodeNbInShell.size() < tmpShell.size())
			nodeNbInShell.push_back(std::vector<double>());
		for (unsigned int i = 0 ; i < tmpShell.size() ; i++)
			nodeNbInShell[i].push_back(tmpShell[i]);


		// For each involved cell A
		for (std::set<unsigned int>::iterator it = invCellOnFront.begin() ;
			it != invCellOnFront.end() ; ++it)
		{
			// Compute shell influxes
			shellStructure[(*distMetr)[*it][init]].push_back(NodeValues(*it, 
				degreeMetr[*it], computeMaxInflux(*it, (*activCells)[w], model.GetNetwork(), maxDelay),
				actCells->ComputeMeanInflux(*it), (*distMetr)[*it][init]));

			isOnFront = true;
			// For each other cell B
			for (std::set<unsigned int>::iterator it2 = invCellOnFront.begin() ;
				it2 != invCellOnFront.end() and isOnFront ; ++it2)
			{
				// If A is connected to B and B is further away,
				// A is not on the wave frontier
				if (model.GetNetwork().AreConnected(*it, *it2) and 
					((*distMetr)[*it][init] < (*distMetr)[*it2][init]))
					isOnFront = false;
			}
			// Add the cell to the wave frontier if needed
			if (isOnFront)
				frontiers[w].push_back(NodeValues(*it, degreeMetr[*it], 
					computeMaxInflux(*it, (*activCells)[w], model.GetNetwork(), maxDelay),
					actCells->ComputeMeanInflux(*it),
					//actCells->ComputeMeanInflux(*it),
					(*distMetr)[*it][init]));
		}

		// For each cell, look at unactivated cells that have at least one
		// activated neighbor.
		nonActTemp.clear();
		std::set<unsigned int> tempNeighbs;
		for (std::set<unsigned int>::iterator it = invCellOnFront.begin() ;
			it != invCellOnFront.end() ; ++it)
		{
			tempNeighbs = std::set<unsigned int>(model.GetNetwork().GetNeighbors(*it).begin(),
				model.GetNetwork().GetNeighbors(*it).end());
			// For each neighbor
			for (std::set<unsigned int>::iterator it2 = tempNeighbs.begin() ;
				it2 != tempNeighbs.end() ; ++it2)
			{
				if (invCellOnFront.find(*it2) == invCellOnFront.end())
					nonActTemp.insert(*it2);
			}
		}
		// Add results to nonActNeighbs
		for (std::set<unsigned int>::iterator it = nonActTemp.begin() ;
				it != nonActTemp.end() ; ++it)
			nonActNeighbs[w].push_back(NodeValues(*it, degreeMetr[*it], 
					computeMaxInflux(*it, (*activCells)[w], model.GetNetwork(), maxDelay),
					actCells->GetMeanInfluxAfterSpike(*it),
					(*distMetr)[*it][init]));
	}

	// Computing mean values and StdDevValues
	// For front cells
	std::vector<double> tempDeg, tempInflux, tempActualInflux, tempDistToInit;
	for (unsigned int w = 0 ; w < frontiers.size() ; ++w)
		for (unsigned int i = 0 ; i < frontiers[w].size() ; ++i)
		{
			tempDeg.push_back(frontiers[w][i].degree);
			tempInflux.push_back(frontiers[w][i].influx);
			tempActualInflux.push_back(frontiers[w][i].actualInflux);
			tempDistToInit.push_back(frontiers[w][i].distToInit);
		}
	meanDegFront = ComputeMean(tempDeg);
	meanInfluxFront = ComputeMean(tempInflux);
	meanActualInfluxFront = ComputeMean(tempActualInflux);
	meanDistToInitFront = ComputeMean(tempDistToInit);
	stdDevDegFront = ComputeStdDev(tempDeg);
	stdDevInfluxFront = ComputeStdDev(tempInflux);
	stdDevActualInfluxFront = ComputeStdDev(tempActualInflux);
	stdDevDistToInitFront = ComputeStdDev(tempDistToInit);

	// for non act neighbs cells
	tempDeg.clear();
	tempInflux.clear();
	tempActualInflux.clear();
	tempDistToInit.clear();
	for (unsigned int w = 0 ; w < nonActNeighbs.size() ; ++w)
		for (unsigned int i = 0 ; i < nonActNeighbs[w].size() ; ++i)
		{
			tempDeg.push_back(nonActNeighbs[w][i].degree);
			tempInflux.push_back(nonActNeighbs[w][i].influx);
			tempActualInflux.push_back(nonActNeighbs[w][i].actualInflux);
			tempDistToInit.push_back(nonActNeighbs[w][i].distToInit);
		}
	meanDegNonAct = ComputeMean(tempDeg);
	meanInfluxNonAct = ComputeMean(tempInflux);
	meanActualInfluxNonAct = ComputeMean(tempActualInflux);
	meanDistToInitNonAct = ComputeMean(tempDistToInit);
	stdDevDegNonAct = ComputeStdDev(tempDeg);
	stdDevInfluxNonAct = ComputeStdDev(tempInflux);
	stdDevActualInfluxNonAct = ComputeStdDev(tempActualInflux);
	stdDevDistToInitNonAct = ComputeStdDev(tempDistToInit);
	return true;
}

//**********************************************************************
// Save the activated cells
//**********************************************************************
bool WaveFrontDetect::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string fullFrontierCellsName("FullFrontierCells");
	if (saver.isSaving(fullFrontierCellsName) and not frontiers.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(fullFrontierCellsName, saver.getCurrFile());
		
		stream << "IdWave\tIdCell\tIsOnFront\tDegree\tIncomingFlux\tActualInFlux\tDistToInit" << std::endl;

		for (unsigned int w = 0 ; w < frontiers.size() ; ++w)
		{
			for (unsigned int j = 0 ; j < frontiers[w].size() ; ++j)
			{
				stream 
					<< w << "\t" 
					<< frontiers[w][j].ind << "\t"
					<< 1 << "\t"
					<< frontiers[w][j].degree << "\t"
					<< frontiers[w][j].influx << "\t"
					<< frontiers[w][j].actualInflux << "\t"
					<< frontiers[w][j].distToInit << std::endl;
			}
			for (unsigned int j = 0 ; j < nonActNeighbs[w].size() ; ++j)
			{
				stream 
					<< w << "\t" 
					<< nonActNeighbs[w][j].ind << "\t"
					<< 0 << "\t"
					<< nonActNeighbs[w][j].degree << "\t"
					<< nonActNeighbs[w][j].influx << "\t"
					<< nonActNeighbs[w][j].actualInflux << "\t"
					<< nonActNeighbs[w][j].distToInit << std::endl;
			}
		}

		allSaved &= stream.good();
	}

	std::string frontierCellsStatsName("FrontierCellsStats");
	if (saver.isSaving(frontierCellsStatsName) and not frontiers.empty())
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(frontierCellsStatsName, saver.getCurrFile());
		
		stream 
			<< "MeanDegFront\tMeanDegNonAct\t"
			<< "MeanInfluxFront\tMeanInfluxNonAct\t"
			<< "MeanActualInfluxFront\tMeanActualInfluxNonAct\t"
			<< "MeanDistToInitFront\tMeanDistToInitNonAct\t"
			<< "StdDevDegFront\tStdDevDegNonAct\t"
			<< "StdDevInfluxFront\tStdDevInfluxNonAct\t"
			<< "StdDevActualInfluxFront\tStdDevActualInfluxNonAct\t"
			<< "StdDevDistToInitFront\tStdDevDistToInitFrontNonAct" << std::endl;
		stream 
			<< meanDegFront << "\t"
			<< meanDegNonAct << "\t"
			<< meanInfluxFront << "\t"
			<< meanInfluxNonAct << "\t"
			<< meanActualInfluxFront << "\t"
			<< meanActualInfluxNonAct << "\t"
			<< meanDistToInitFront << "\t"
			<< meanDistToInitNonAct << "\t"
			<< stdDevDegFront << "\t"
			<< stdDevDegNonAct << "\t"
			<< stdDevInfluxFront << "\t"
			<< stdDevInfluxNonAct << "\t"
			<< stdDevActualInfluxFront << "\t"
			<< stdDevActualInfluxNonAct << "\t"
			<< stdDevDistToInitFront << "\t"
			<< stdDevDistToInitNonAct << std::endl;

		allSaved &= stream.good();
	}
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void WaveFrontDetect::Initialize()
{
	Metric::Initialize();
	frontiers.clear();
	nonActNeighbs.clear();
	shellStructure.clear();
	nodeNbInShell.clear();
	maxInfluxByShell.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * WaveFrontDetect::BuildCopy() const
{
	return new WaveFrontDetect(useParentsOfParentsAsNonSinks); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> WaveFrontDetect::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	res["WaveFrontInfluxes"] = 0;
	res["WaveFrontActualInfluxes"] = 0;
	res["WaveFrontDegree"] = 0;
	res["WaveFrontDistToInit"] = 0;
	double totalNb = 0;
	for (unsigned int w = 0 ; w < frontiers.size() ; ++w)
		for (unsigned int i = 0 ; i < frontiers[w].size() ; ++i)
		{
			res["WaveFrontInfluxes"] += frontiers[w][i].influx;
			res["WaveFrontActualInfluxes"] += frontiers[w][i].actualInflux;
			res["WaveFrontDegree"] += frontiers[w][i].degree;
			res["WaveFrontDistToInit"] += frontiers[w][i].distToInit;
			totalNb += 1;
		}

	res["WaveFrontInfluxes"] /= totalNb;
	res["WaveFrontActualInfluxes"] /= totalNb;
	res["WaveFrontDegree"] /= totalNb;
	res["WaveFrontDistToInit"] /= totalNb;


	res["NonActCellsInfluxes"] = 0;
	res["NonActCellsActualInfluxes"] = 0;
	res["NonActCellsDegree"] = 0;
	res["NonActCellsDistToInit"] = 0;
	totalNb = 0;
	for (unsigned int w = 0 ; w < nonActNeighbs.size() ; ++w)
		for (unsigned int i = 0 ; i < nonActNeighbs[w].size() ; ++i)
		{
			res["NonActCellsInfluxes"] += nonActNeighbs[w][i].influx;
			res["NonActCellsActualInfluxes"] += nonActNeighbs[w][i].actualInflux;
			res["NonActCellsDegree"] += nonActNeighbs[w][i].degree;
			res["NonActCellsDistToInit"] += nonActNeighbs[w][i].distToInit;
			totalNb += 1;
		}

	res["NonActCellsInfluxes"] /= totalNb;
	res["NonActCellsActualInfluxes"] /= totalNb;
	res["NonActCellsDegree"] /= totalNb;
	res["NonActCellsDistToInit"] /= totalNb;

	// Shell structure influxes
	double tmp = 0;
	std::set<unsigned int> cellNb;
	for (std::map<double, std::vector<NodeValues> >::const_iterator it = 
		shellStructure.begin() ; it != shellStructure.end() ; ++it)
	{
		tmp = 0.0;
		cellNb.clear();
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
		{
			tmp += it->second[i].actualInflux;
			cellNb.insert(it->second[i].ind);
		}
		res[std::string("InfluxInActivatedCellsOfShell_") + StringifyFixed(it->first)] = 
			tmp / (double) it->second.size();
		res[std::string("NbActCellsInShell_") + StringifyFixed(it->first)] = cellNb.size();
		res[std::string("NbNodesInShell_") + StringifyFixed(it->first)] = 
			ComputeMean(nodeNbInShell[round(it->first)]);
		res[std::string("ActivDensityInShell_") + StringifyFixed(it->first)] = 
			cellNb.size() / ComputeMean(nodeNbInShell[round(it->first)]);
	}

	// Max influx by shell
	for (unsigned int i = 0 ; i < maxInfluxByShell.size() ; ++i)
	{
		res[std::string("MaxInfluxInShell_") + StringifyFixed(i)] = 
			ComputeMean(maxInfluxByShell[i]);
	}

	return res;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool WaveFrontDetect::LoadFromStream(std::ifstream & stream)
{
	stream >> useParentsOfParentsAsNonSinks;
	this->Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool WaveFrontDetect::SaveToStream(std::ofstream & stream) const
{
	stream << useParentsOfParentsAsNonSinks << std::endl;
	return stream.good();
}

//**********************************************************************
// Returns the distribution of influxes for cells that are on the wave 
// front (the last cells that have been activated)
//**********************************************************************
std::vector<double> WaveFrontDetect::GetFrontInluxes() const
{
	std::vector<double> result;
	for (unsigned int w = 0 ; w < frontiers.size() ; ++w)
		for (unsigned int i = 0 ; i < frontiers[w].size() ; ++i)
			result.push_back(frontiers[w][i].influx);
	return result;
}

//**********************************************************************
// Returns the distribution of influxes for cells that are linked to
// activated cells but which didn't get activated
//**********************************************************************
std::vector<double> WaveFrontDetect::GetNonActInfluxes() const
{
	std::vector<double> result;
	for (unsigned int w = 0 ; w < nonActNeighbs.size() ; ++w)
		for (unsigned int i = 0 ; i < nonActNeighbs[w].size() ; ++i)
			result.push_back(nonActNeighbs[w][i].influx);
	return result;
}

//**********************************************************************
// Compute the maximum theoretical influx received by the cell 
// during a given wave.
//**********************************************************************
double WaveFrontDetect::computeMaxInflux(unsigned int ind, 
	const TimedWave & wave, const AbstractNetwork & net, double maxDelay) const
{
	double result, tempExc, sumExc;
	// If the cell has been activated in the wave
	if (wave.involvedCells.find(ind) != wave.involvedCells.end())
	{
		result = 0;
		// For each activation
		for (unsigned int i = 0 ; i < wave.structure.size() ; ++i)
		{
			if (wave.structure[i].ind == ind)
			{
				sumExc = 0;
				// For each parent cell
				for (unsigned int j = 0 ; 
					j < wave.structure[i].parents.size() ; ++j)
				{
					tempExc = 0;
					// Find the parent activation
					unsigned int k;
					for (k = 0 ; (k < wave.structure.size()) and 
						(wave.propagation[k].time < wave.propagation[i].time) ; ++k)
					{
						if (wave.structure[k].ind == wave.structure[i].parents[j])
						{
							if (useParentsOfParentsAsNonSinks)
							{
								tempExc = 1.0 / std::max(1.0, (net.GetNodeDegree(wave.structure[k].ind) - 
									(double) wave.structure[k].parents.size()));
							}
							else
							{
								double nbActNeighbs = 0;
								for (unsigned int m = 0 ; (m < wave.structure.size()) and 
									(wave.propagation[m].time < wave.propagation[i].time) ; ++m)
								{
									if ((wave.propagation[m].time >= wave.propagation[k].time) and 
											(net.AreConnected(wave.structure[m].ind, wave.structure[k].ind)))
										++nbActNeighbs;
								}
								tempExc = 1.0 / std::max(1.0, net.GetNodeDegree(wave.structure[k].ind) - 
									nbActNeighbs);
							}
						}
					}
					// Add the closest exc (in time) to the sum
					sumExc += tempExc;
				}
				// Take the maximum influx
				if (result < sumExc)
					result = sumExc;
			}
		}
	}
	// If the cell hasn't been activated in the wave
	else
	{
		sumExc = 0;
		std::set<unsigned int> neighbs(net.GetNeighbors(ind).begin(), net.GetNeighbors(ind).end());
		// For each neighbor
		for (std::set<unsigned int>::iterator it = neighbs.begin() ; 
			it != neighbs.end() ; ++it)
		{
			tempExc = 0;
			// If the neighbor has been activated
			if (wave.involvedCells.find(*it) != wave.involvedCells.end())
			{
				// For each activation of this neighbor
				for (unsigned int i = 0 ; i < wave.structure.size() ; ++i)
				{
					if (wave.structure[i].ind == *it)
					{
						if (useParentsOfParentsAsNonSinks)
						{
							tempExc = std::max(tempExc, 1.0 / std::max(1.0, (net.GetNodeDegree(*it) - 
								(double)wave.structure[i].parents.size())));
						}
						else
						{
							double nbActNeighbs = 0;
							for (unsigned int m = 0 ; (m < wave.structure.size()) and 
								(wave.propagation[m].time < (wave.propagation[i].time + maxDelay)) ; ++m)
							{
								if ((wave.propagation[m].time >= wave.propagation[i].time) and 
										(net.AreConnected(wave.structure[m].ind, *it)))
									++nbActNeighbs;
							}
							tempExc = 1.0 / std::max(1.0, net.GetNodeDegree(*it) - 
								nbActNeighbs);
						}
					}
				}
			}
			sumExc += tempExc;
		}
		result = sumExc;
	}
	return result;
}

//********************************************************************//
//********* F O U R I E R   T R A N S F O R M   M E T R I C **********//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(FOURIER_TRANSFORM_METRIC, CONCENTRATIONS_METRIC)

//**********************************************************************
// Default Constructor
//**********************************************************************
FourierTransform::FourierTransform(ParamHandler & h) :
	wavetable(0), workspace(0), sigData_comp(0)
{
	std::vector<std::string> names = 
		h.getParam<std::vector<std::string> >("-FourierTransform", 0);
	std::vector<std::string> totElems =
		h.getParam<std::vector<std::string> >("-FourierTransform", 1);
	assert(names.size() == totElems.size());
	minFreqToSave = h.getParam<double>("-FourierTransform", 2);
	maxFreqToSave = h.getParam<double>("-FourierTransform", 3);
	stepFreqToSave = h.getParam<double>("-FourierTransform", 4);

	for (unsigned int i = 0 ; i < names.size() ; ++i)
		elems[names[i]].push_back(totElems[i]);

	concatenateSigs = h.getParam<bool>("-FourierTransConcatSigs");
	useShortTimeFFT = h.getParam<bool>("-FourierTransUseShortTimeFFT", 0);
	stfftWinSize = h.getParam<double>("-FourierTransUseShortTimeFFT", 1);
	stfftWinStep = h.getParam<double>("-FourierTransUseShortTimeFFT", 2);
	domFreqThreshRat = h.getParam<double>("-FourierTransUseShortTimeFFT", 3);
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
FourierTransform::FourierTransform(std::ifstream & stream) :
	wavetable(0), workspace(0), sigData_comp(0)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
FourierTransform::~FourierTransform()
{
	if (sigData_comp)
		delete[] sigData_comp;
	if (wavetable)
		gsl_fft_complex_wavetable_free(wavetable);
	if (workspace)
		gsl_fft_complex_workspace_free(workspace);
}

//**********************************************************************
// Compute the activated cells
//**********************************************************************
bool FourierTransform::ComputeMetric(const ODE::ODEProblem<double, double> & prob)
{
	// Obtaining needed metrics
	ConcentrationMetrics *concentrMetr = 
		GetSpecificMetric<Metric, ConcentrationMetrics>(prob.GetAllMetrics());

	assert(concentrMetr);
	const std::vector<std::pair<double, std::vector<double> > > & concentrRaw =
		concentrMetr->GetConcentrations();
	const std::vector<std::string> & valNames = concentrMetr->GetValNames();
	double deltaT = concentrMetr->GetSavingStep();

TRACE(stfftWinSize << " // " << deltaT)
	unsigned int stfftDiscrWinSize = ceil(stfftWinSize / deltaT);
	unsigned int stfftDiscrWinStep = ceil(stfftWinStep / deltaT);
	assert(stfftDiscrWinSize > 1);

	double tmpMin = DEFAULT_MAX_VAL;
	double tmpMax = -DEFAULT_MAX_VAL;
	// First, tranform concentr to a more usable format
	unsigned int timeLen = concentrRaw.size();
	std::vector<std::vector<double> > concentr(valNames.size(), 
		std::vector<double>(concentrRaw.size(), 0));
	for (unsigned int i = 0 ; i < valNames.size() ; ++i)
		for (unsigned int t = 0 ; t < concentrRaw.size() ; ++t)
		{
			concentr[i][t] = concentrRaw[t].second[i];
			tmpMin = std::min(tmpMin, concentr[i][t]);
			tmpMax = std::max(tmpMax, concentr[i][t]);
		}
	double totAmplAcrossSigs = tmpMax - tmpMin;

	assert(stfftDiscrWinSize <= timeLen);

	// For each signal to be analyzed
	for (std::map<std::string, std::vector<std::string> >::iterator it = elems.begin() ; 
		it != elems.end() ; ++it)
	{
		// Construct signal mask
		std::vector<bool> sigMask(valNames.size(), true);
		for (unsigned int e = 0 ; e < elems[it->first].size() ; ++e)
			for (unsigned int i = 0 ; i < valNames.size() ; ++i)
				sigMask[i] = sigMask[i] and 
					(valNames[i].find(elems[it->first][e]) != std::string::npos);

		// Compute nbSubSigs
		std::vector<unsigned int> sigInds;
		for (unsigned int i = 0 ; i < sigMask.size() ; ++i)
			if (sigMask[i])
				sigInds.push_back(i);

		// Compute the actual tranforms
		if (concatenateSigs)
		{
			// Concatenate signals
			unsigned int tmpInd = 0;
			std::vector<double> totSig(sigInds.size() * timeLen, 0);
			for (unsigned int i = 0 ; i < sigInds.size() ; ++i)
				for (unsigned int t = 0 ; t < concentr[i].size() ; ++t)
					totSig[tmpInd++] = concentr[i][t];

			//centerAndNormSig(totSig);

			// If short time FFT is needed
			if (useShortTimeFFT)
			{
				initializeSpectrum(it->first, stfftDiscrWinSize, deltaT);
				double nbwins = 0;
				for (unsigned int t = 0 ; t <= (totSig.size() - stfftDiscrWinSize) ;
					t += stfftDiscrWinStep)
				{
					computeAndAddSpectrum(totSig, t, t + stfftDiscrWinSize, it->first);
					++nbwins;
				}
				// Normalize spectrum
				for (unsigned int i = 0 ; i < spectrums[it->first].size() ; ++i)
					spectrums[it->first][i].second /= nbwins;
			}
			// Otherwise perform full FFT on signal
			else
			{
				initializeSpectrum(it->first, totSig.size(), deltaT);
				computeAndAddSpectrum(totSig, 0, totSig.size(), it->first);
			}
		}
		else // Separated sub signals, resulting spectrum is mean spectrum
		{
			// Merged signal :
			if (useShortTimeFFT)
				initializeSpectrum(it->first, stfftDiscrWinSize, deltaT);
			else
				initializeSpectrum(it->first, timeLen, deltaT);

			// Compute spectrum sur each sub signal
			double nbComps = 0;
			for (unsigned int i = 0 ; i < sigInds.size() ; ++i)
			{
				if (useShortTimeFFT)
				{
					for (unsigned int t = 0 ; t <= (timeLen - stfftDiscrWinSize) ; 
						t += stfftDiscrWinStep)
					{
						spectrograms[it->first].push_back(
							std::vector<double>(spectrums[it->first].size(), 0));
						for (unsigned int j = 0 ; j < spectrums[it->first].size() ; ++j)
							spectrograms[it->first].back()[j] = spectrums[it->first][j].second;
						computeAndAddSpectrum(concentr[sigInds[i]], t, 
							t + stfftDiscrWinSize, it->first);
						++nbComps;
						for (unsigned int j = 0 ; j < spectrums[it->first].size() ; ++j)
							spectrograms[it->first].back()[j] = 
								spectrums[it->first][j].second - spectrograms[it->first].back()[j];
						if ((ComputeMax(concentr[sigInds[i]], t, t + stfftDiscrWinSize) - 
							ComputeMin(concentr[sigInds[i]], t, t + stfftDiscrWinSize)) > 
							domFreqThreshRat * totAmplAcrossSigs)
						{
							unsigned int indMax = GetIndiceMax(spectrograms[it->first].back(), 
								0, spectrograms[it->first].back().size() / 2);
							dominantFrequencies[it->first].push_back(
								spectrums[it->first][indMax].first);
						}
						else
							dominantFrequencies[it->first].push_back(0);
					}
				}
				else
				{
					computeAndAddSpectrum(concentr[sigInds[i]], 0, timeLen, it->first);
					++nbComps;
				}
			}
			// Normalize to mean
			TRACE(nbComps)
			for (unsigned int i = 0 ; i < spectrums[it->first].size() ; ++i)
			{
				spectrums[it->first][i].second /= nbComps;
			}

		}
	}
	return true;
}

//**********************************************************************
// Initilize spectrum with frequency values and put spectrum vals to 0
//**********************************************************************
void FourierTransform::initializeSpectrum(const std::string & specName, 
	unsigned int len, double dt)
{
	double freq;
	for (unsigned int i = 0 ; i < len ; ++i)
	{
		if (2*i <= len)
			freq = ((double)i) / ((double)len*dt);
		else
			freq = -((double)(len-i)) / ((double)len*dt);
		spectrums[specName].push_back(make_pair(freq, 0));
	}
}

//**********************************************************************
// Compute spectrum on the given signal delimited by iStart and iEnd, 
// add the spectrum to spectrums
//**********************************************************************
void FourierTransform::computeAndAddSpectrum(const std::vector<double> & sig, 
	unsigned int iStart, unsigned int iEnd, const std::string & specName)
{
	// Check that subsig length is correct
	if ((iEnd - iStart) != n)
	{
		n = iEnd - iStart;

		if (sigData_comp)
			delete[] sigData_comp;
		if (wavetable)
			gsl_fft_complex_wavetable_free(wavetable);
		if (workspace)
			gsl_fft_complex_workspace_free(workspace);

		sigData_comp = new double[2 * n];
		wavetable = gsl_fft_complex_wavetable_alloc(n);
		workspace = gsl_fft_complex_workspace_alloc(n);
	}

	double meanVal = 0;
	double maxVal = sig[0];
	double minVal = maxVal;
	// Compute min, max and mean values
	for (unsigned int t = iStart ; t < iEnd ; ++t)
	{
		meanVal += sig[t];
		if (maxVal < sig[t])
			maxVal = sig[t];
		if (minVal > sig[t])
			minVal = sig[t];
	}
	meanVal /= (double)(iEnd - iStart);

	// Fill complex table
	for (unsigned int t = 0 ; t < n ; ++t)
	{
		// Re Part
		sigData_comp[2*t] = (sig[iStart + t] - meanVal) / (maxVal - minVal);
		// Im Part
		sigData_comp[2*t+1] = 0;
	}

	// Compute actual fourier transform
	gsl_fft_complex_forward(sigData_comp, 1, n, wavetable, workspace);

	// Copy FFT to vector
	for (unsigned int j = 0 ; j < n ; ++j)
		spectrums[specName][j].second += 
			sqrt(sigData_comp[2*j]*sigData_comp[2*j] +
			sigData_comp[2*j+1]*sigData_comp[2*j+1]);
}

//**********************************************************************
// Save the activated cells
//**********************************************************************
bool FourierTransform::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;

	std::string FourDomFreqName("FourierDominantFrequencies");
	if (saver.isSaving(FourDomFreqName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(FourDomFreqName, saver.getCurrFile());

		for (std::map<std::string, std::vector<double> >::const_iterator it = 
				dominantFrequencies.begin() ; it != dominantFrequencies.end() ; ++it)
			stream << it->first << "\t";
		stream << std::endl;
		bool oneWritten = false;
		unsigned int currInd = 0;
		do
		{
			oneWritten = false;
			for (std::map<std::string, std::vector<double> >::const_iterator it = 
				dominantFrequencies.begin() ; it != dominantFrequencies.end() ; ++it)
			{
				if (currInd < it->second.size())
				{
					stream << it->second[currInd] << "\t";
					oneWritten |= true;
				}
				else
					stream << "NaN\t";
			}
			stream << std::endl;
			++currInd;
		} while (oneWritten);

		allSaved &= stream.good();
	}

	std::string FullSpectrogramName("FullFourierSpectrogram");
	if (saver.isSaving(FullSpectrogramName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(FullSpectrogramName, saver.getCurrFile());
	
		for (std::map<std::string, std::vector<std::vector<double> > >::const_iterator it = 
				spectrograms.begin() ; it != spectrograms.end() ; ++it)
		{
			allSaved &= SaveFullMatrix(stream, it->second);
			stream << std::endl;
		}
	}
	std::string fullFourTransName("FullFourierTransform");
	if (saver.isSaving(fullFourTransName))
	{
		ofstream & stream = saver.getStream();
		this->AddSavedFile(fullFourTransName, saver.getCurrFile());

		for (std::map<std::string, std::vector<std::pair<double, double> > >::
				const_iterator it = spectrums.begin() ; it != spectrums.end() ; ++it)
			stream << it->first << "_freq\t" << it->first << "_val\t";
		stream << std::endl;
		bool oneWritten = false;
		unsigned int currInd = 0;
		do
		{
			oneWritten = false;
			for (std::map<std::string, std::vector<std::pair<double, double> > >::
				const_iterator it = spectrums.begin() ; it != spectrums.end() ; ++it)
			{
				if (currInd < it->second.size())
				{
					stream << it->second[currInd].first << "\t" 
					       << it->second[currInd].second << "\t";
					oneWritten |= true;
				}
				else
					stream << "NaN\t";
			}
			stream << std::endl;
			++currInd;
		} while (oneWritten);

		allSaved &= stream.good();
	}
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void FourierTransform::Initialize()
{
	Metric::Initialize();
	spectrums.clear();
	spectrograms.clear();
	dominantFrequencies.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * FourierTransform::BuildCopy() const
{
	return new FourierTransform(); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> FourierTransform::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;

	for (std::map<std::string, std::vector<double> >::const_iterator it = 
		dominantFrequencies.begin() ; it != dominantFrequencies.end() ; ++it)
	{
		std::set<double> tmpFreq(it->second.begin(), it->second.end());
		for (double freq = minFreqToSave ; freq <= maxFreqToSave ; freq += stepFreqToSave)
		{
			// CumFreq
			double cmFreq = 0;
			for (std::set<double>::const_iterator it2 = tmpFreq.begin() ; 
				it2 != tmpFreq.end() ; ++it2)
			{
				if (*it2 >= freq)
					cmFreq += std::count(it->second.begin(), it->second.end(), *it2);
			}
			res.insert(std::make_pair(
				std::string("DomFreq_CumFreqOver_")+it->first+"_"+Stringify(freq)+"_Hz", 
				cmFreq / ((double) it->second.size())));
		}
	}

	return res;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool FourierTransform::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();

	stream >> concatenateSigs;
	stream >> useShortTimeFFT;
	stream >> stfftWinSize;
	stream >> stfftWinStep;
	stream >> domFreqThreshRat;

	stream >> minFreqToSave;
	stream >> maxFreqToSave;
	stream >> stepFreqToSave;

	elems.clear();
	unsigned int tmpNb, tmpNb2;
	std::string tmpStr, tmpStr2;
	stream >> tmpNb;
	for (unsigned int i = 0 ; i < tmpNb ; ++i)
	{
		stream >> tmpStr;
		stream >> tmpNb2;
		for (unsigned int j = 0 ; j < tmpNb2 ; ++j)
		{
			stream >> tmpStr2;
			elems[tmpStr].push_back(tmpStr2);
		}
	}

	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool FourierTransform::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< concatenateSigs << std::endl
		<< useShortTimeFFT << std::endl
		<< stfftWinSize << std::endl
		<< stfftWinStep << std::endl
		<< domFreqThreshRat << std::endl;

	stream
		<< minFreqToSave << std::endl
		<< maxFreqToSave << std::endl
		<< stepFreqToSave << std::endl;

	stream << elems.size() << std::endl;
	for (std::map<std::string, std::vector<std::string> >::const_iterator
		it = elems.begin() ; it != elems.end() ; ++it)
	{
		stream << it->first << std::endl;
		stream << it->second.size() << std::endl;
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			stream << it->second[i] << std::endl;
	}

	return stream.good();
}

//********************************************************************//
//********* W A V E L E T   T R A N S F O R M   M E T R I C **********//
//********************************************************************//

//**********************************************************************
// Dependencies
//**********************************************************************
ADD_METRIC_DEPENDENCY(WAVELET_TRANSFORM_METRIC, CONCENTRATIONS_METRIC)

//**********************************************************************
// Default Constructor
//**********************************************************************
WaveletTransform::WaveletTransform(ParamHandler & h) :
	wavelet(0), workspace(0), sigData(0), n(0), newN(0), nbLevels(0)
{
	std::vector<std::string> names = 
		h.getParam<std::vector<std::string> >("-WaveletTransformElems", 0);
	std::vector<std::string> totElems =
		h.getParam<std::vector<std::string> >("-WaveletTransformElems", 1);
	assert(names.size() == totElems.size());
	minFreqToSave = h.getParam<double>("-WaveletTransform", 0);
	maxFreqToSave = h.getParam<double>("-WaveletTransform", 1);
	stepFreqToSave = h.getParam<double>("-WaveletTransform", 2);
	domFreqThreshRat = h.getParam<double>("-WaveletTransform", 3);
	powrThreshFract = h.getParam<double>("-WaveletTransform", 4);

	for (unsigned int i = 0 ; i < names.size() ; ++i)
		elems[names[i]].push_back(totElems[i]);
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
WaveletTransform::WaveletTransform(std::ifstream & stream) :
	wavelet(0), workspace(0), sigData(0), n(0), newN(0), nbLevels(0)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
WaveletTransform::~WaveletTransform()
{
	if (sigData)
		delete[] sigData;
	if (wavelet)
		gsl_wavelet_free(wavelet);
	if (workspace)
		gsl_wavelet_workspace_free(workspace);
}

//**********************************************************************
// Compute the activated cells
//**********************************************************************
bool WaveletTransform::ComputeMetric(const ODE::ODEProblem<double, double> & prob)
{
	// Obtaining needed metrics
	ConcentrationMetrics *concentrMetr = 
		GetSpecificMetric<Metric, ConcentrationMetrics>(prob.GetAllMetrics());

	assert(concentrMetr);
	const std::vector<std::pair<double, std::vector<double> > > & concentrRaw =
		concentrMetr->GetConcentrations();
	const std::vector<std::string> & valNames = concentrMetr->GetValNames();

	double samplingFreq = 1.0 / concentrMetr->GetSavingStep();
	// First, tranform concentr to a more usable format
	unsigned int timeLen = concentrRaw.size();
	std::vector<std::vector<double> > concentr(valNames.size(), 
		std::vector<double>(concentrRaw.size(), 0));
	for (unsigned int i = 0 ; i < valNames.size() ; ++i)
		for (unsigned int t = 0 ; t < concentrRaw.size() ; ++t)
			concentr[i][t] = concentrRaw[t].second[i];

	// For each signal to be analyzed
	for (std::map<std::string, std::vector<std::string> >::iterator it = elems.begin() ; 
		it != elems.end() ; ++it)
	{
		// Construct signal mask
		std::vector<bool> sigMask(valNames.size(), true);
		for (unsigned int e = 0 ; e < elems[it->first].size() ; ++e)
			for (unsigned int i = 0 ; i < valNames.size() ; ++i)
				sigMask[i] = (sigMask[i] and (valNames[i].find(elems[it->first][e]) != std::string::npos));

		// Compute nbSubSigs
		std::vector<unsigned int> sigInds;
		for (unsigned int i = 0 ; i < sigMask.size() ; ++i)
			if (sigMask[i])
				sigInds.push_back(i);

		totNbSignals[it->first] = sigInds.size();
		double tmpMin;
		double tmpMax;
		unsigned int maxAmplInd = 0;
		std::vector<double> sigAmp(sigInds.size(), 0);
		// First, center signals
		for (unsigned int i = 0 ; i < sigInds.size() ; ++i)
		{
			double meanVal = 0;
			tmpMin = DEFAULT_MAX_VAL;
			tmpMax = -DEFAULT_MAX_VAL;
			// Compute mean values
			for (unsigned int t = 0 ; t < timeLen ; ++t)
			{
				meanVal += concentr[sigInds[i]][t];
				tmpMin = std::min(tmpMin, concentr[sigInds[i]][t]);
				tmpMax = std::max(tmpMax, concentr[sigInds[i]][t]);
			}
			sigAmp[i] = tmpMax - tmpMin;
			maxAmplInd = (sigAmp[i] > sigAmp[maxAmplInd]) ? i : maxAmplInd;
			meanVal /= (double)timeLen;
			// Center 
			for (unsigned int t = 0 ; t < timeLen ; ++t)
				concentr[sigInds[i]][t] = concentr[sigInds[i]][t] - (tmpMax + tmpMin)/2.0;
		}
		double totAmplAcrossSigs = sigAmp[maxAmplInd];

		// Compute wavelet transforms
		for (unsigned int i = 0 ; i < sigInds.size() ; ++i)
		{
			// Don't compute wavelet transform if the signal amplitude is too small
			if (sigAmp[i] > domFreqThreshRat*totAmplAcrossSigs)
			{
				computeAndAddSpectrogram(concentr[sigInds[i]], 0, concentr[sigInds[i]].size(),
					it->first, samplingFreq);
				computeAndAddDominantFrequencies(it->first, spectrograms[it->first].size() - 1);
			}
		}
	}

	return true;
}

//**********************************************************************
// Compute spectrum on the given signal delimited by iStart and iEnd, 
// add the spectrum to spectrums
//**********************************************************************
void WaveletTransform::computeAndAddSpectrogram(const std::vector<double> & sig, 
	unsigned int iStart, unsigned int iEnd, const std::string & specName, double samplingFreq)
{
	// Check that subsig length is correct
	if ((iEnd - iStart) != n)
	{
		n = iEnd - iStart;

		if (sigData)
			delete[] sigData;
		if (wavelet)
			gsl_wavelet_free(wavelet);
		if (workspace)
			gsl_wavelet_workspace_free(workspace);

		double expo = gsl_sf_log((double)n) / gsl_sf_log(2.0);
		nbLevels = floor(expo);
		newN = pow(2, nbLevels);
		sigData = new double[newN];

		wavelet = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
		workspace = gsl_wavelet_workspace_alloc(n);
	}

	// Compute the interpolation for upsampling
	double *xa = new double[n];
	double *ya = new double[n];
	for (unsigned int i = iStart ; i < iEnd ; ++i)
	{
		xa[i] = i - iStart;
		ya[i] = sig[i];
	}

	gsl_interp *interp = gsl_interp_alloc(gsl_interp_linear, n);
	gsl_interp_init(interp, xa, ya, n);
	gsl_interp_accel *tmpAcc = gsl_interp_accel_alloc();


	// Upsampling of signal to make it a power of 2
	for (unsigned int i = 0 ; i < newN ; ++i)
	{
		sigData[i] = gsl_interp_eval(interp, xa, ya, ((double)i)/((double)newN - 1.0)*((double)n - 1.0), tmpAcc);
	}

	delete[] xa;
	delete[] ya;
	gsl_interp_accel_free(tmpAcc);
	gsl_interp_free(interp);

	// Compute actual wavelet transform
	gsl_wavelet_transform_forward(wavelet, sigData, 1, newN, workspace);

	// Copy wavelet transform to spectrogram
	unsigned int ind = 1;
	double freq;
	double corrSamplFreq = samplingFreq * ((double)newN / (double)n);
	spectrograms[specName].push_back(std::vector<std::pair<double, std::vector<double> > >());
	for (unsigned int l = 0 ; l < nbLevels ; ++l)
	{
		freq = (corrSamplFreq / pow(2, nbLevels - l + 1) + 
			corrSamplFreq / pow(2, nbLevels - l)) / 2.0;
		spectrograms[specName].back().push_back(make_pair(freq, 
			std::vector<double>(pow(2,l), 0)));
		for (unsigned int i = 0 ; i < pow(2,l) ; ++i)
			spectrograms[specName].back().back().second[i] = sigData[ind++];
	}
}

//**********************************************************************
// Compute dominant frequencies on the given signal and indice in spectrogram
// (to be called after computeAndAddSpectrogram)
//**********************************************************************
void WaveletTransform::computeAndAddDominantFrequencies(
	const std::string & specName, unsigned int ind)
{
	unsigned int tmpT;
	std::vector<double> powerThreshMax(spectrograms[specName][ind].size(), 0);
	std::vector<double> powerThreshMin(spectrograms[specName][ind].size(), 0);
	double allTimeMax = -DEFAULT_MAX_VAL;
	double allTimeMin = DEFAULT_MAX_VAL;
	// Compute power threshold per level
	for (unsigned int l = 0 ; l < spectrograms[specName][ind].size() ; ++l)
	{
		powerThreshMax[l] = fabs(spectrograms[specName][ind][l].second[0]);
		powerThreshMin[l] = fabs(spectrograms[specName][ind][l].second[0]);
		// Get maximum power for this level
		for (unsigned int t = 0 ; t < spectrograms[specName][ind][l].second.size() ; ++t)
		{
			powerThreshMax[l] = (powerThreshMax[l] < fabs(spectrograms[specName][ind][l].second[t])) ?
				fabs(spectrograms[specName][ind][l].second[t]) : powerThreshMax[l];
			powerThreshMin[l] = (powerThreshMin[l] > fabs(spectrograms[specName][ind][l].second[t])) ?
				fabs(spectrograms[specName][ind][l].second[t]) : powerThreshMin[l];
		}
		if ((spectrograms[specName][ind][l].first >= minFreqToSave) and 
			(spectrograms[specName][ind][l].first <= maxFreqToSave))
		{
			allTimeMax = (allTimeMax < powerThreshMax[l]) ? powerThreshMax[l] : allTimeMax;
			allTimeMin = (allTimeMin > powerThreshMin[l]) ? powerThreshMin[l] : allTimeMin;
		}
	}

	bool oneIsHigherThanThresh = false;
	double tmpMean, tmpWeight;
	// for each lowest (high time definition) level
	for (unsigned t = 0 ; t < spectrograms[specName][ind].back().second.size() ; ++t)
	{
		tmpT = t;
		oneIsHigherThanThresh = false;
		tmpMean = 0;
		tmpWeight = 0;
		// Get the level for which the wavelet coefficient is maximal
		for (int l =  spectrograms[specName][ind].size() - 1 ; l >= 0 ; --l)
		{
			if (((fabs(spectrograms[specName][ind][l].second[tmpT]) - allTimeMin) >= 
					powrThreshFract * (allTimeMax - allTimeMin)))
				oneIsHigherThanThresh = true;

			if ((spectrograms[specName][ind][l].first >= minFreqToSave) and 
				(spectrograms[specName][ind][l].first <= maxFreqToSave))
			{
				tmpMean += fabs(spectrograms[specName][ind][l].second[tmpT]) * 
					spectrograms[specName][ind][l].first;
				tmpWeight += fabs(spectrograms[specName][ind][l].second[tmpT]);
			}
			
			tmpT /= 2;
		}
		// Add the pseudo frequency corresponding to this level
		if (oneIsHigherThanThresh)
		{
			dominantFrequencies[specName].push_back(tmpMean / tmpWeight);
		}
	}
}

//**********************************************************************
// Save the activated cells
//**********************************************************************
bool WaveletTransform::SaveMetric(ResultSaver saver) const
{
	bool allSaved = true;
	std::string WaveletSpectroName("WaveletSpectroGrams");
	if (saver.isSaving(WaveletSpectroName))
	{
		for (std::map<std::string, std::vector<std::vector<
			std::pair<double, std::vector<double> > > > >::const_iterator it = spectrograms.begin();
			it != spectrograms.end() ; ++it)
		{
			for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			{
				ResultSaver tmpSav = saver(it->first + "/" + StringifyFixed(i));
				ofstream & stream = tmpSav.getStream();
				this->AddSavedFile(WaveletSpectroName, tmpSav.getCurrFile());
				assert(stream.good());
				// Format : columns are time and rows are levels
				//          first column contains pseudo frequencies
				for (unsigned int l = 0 ; l < it->second[i].size() ; ++l)
				{
					stream << it->second[i][l].first << " ";
					for (unsigned int t = 0 ; t < it->second[i].back().second.size() ; ++t)
					{
						stream << it->second[i][l].second[t / pow(2, it->second[i].size() - 1 - l)] << " ";
					}
					stream << std::endl;
					allSaved &= stream.good();
				}
			}
		}
	}
	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void WaveletTransform::Initialize()
{
	Metric::Initialize();
	spectrograms.clear();
	dominantFrequencies.clear();
	totNbSignals.clear();
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * WaveletTransform::BuildCopy() const
{
	return new WaveletTransform(); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> WaveletTransform::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;

	for (std::map<std::string, std::vector<double> >::const_iterator it = 
		dominantFrequencies.begin() ; it != dominantFrequencies.end() ; ++it)
	{
		std::set<double> tmpFreq(it->second.begin(), it->second.end());

		std::map<std::string, unsigned int>::const_iterator it2 = totNbSignals.find(it->first);
		assert(it2 != totNbSignals.end());
		for (double freq = minFreqToSave ; freq <= maxFreqToSave ; freq += stepFreqToSave)
		{
			double cmFreq = 0;
			for (std::set<double>::const_iterator it3 = tmpFreq.begin() ; it3 != tmpFreq.end() ; ++it3)
			{
				if (*it3 >= freq)
					cmFreq += std::count(it->second.begin(), it->second.end(), *it3);
			}
			res.insert(std::make_pair(
				std::string("Wavelet_DomFreq_CumFreqOver_")+it->first+"_"+Stringify(freq)+"_Hz", 
					cmFreq / ((double) (it2->second*newN) / 2.0)));
		}
	}

	return res;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool WaveletTransform::LoadFromStream(std::ifstream & stream)
{
	this->Initialize();

	stream >> domFreqThreshRat;

	stream >> minFreqToSave;
	stream >> maxFreqToSave;
	stream >> stepFreqToSave;

	stream >> powrThreshFract;

	elems.clear();
	unsigned int tmpNb, tmpNb2;
	std::string tmpStr, tmpStr2;
	stream >> tmpNb;
	for (unsigned int i = 0 ; i < tmpNb ; ++i)
	{
		stream >> tmpStr;
		stream >> tmpNb2;
		for (unsigned int j = 0 ; j < tmpNb2 ; ++j)
		{
			stream >> tmpStr2;
			elems[tmpStr].push_back(tmpStr2);
		}
	}

	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool WaveletTransform::SaveToStream(std::ofstream & stream) const
{
	stream 
		<< domFreqThreshRat << std::endl;

	stream
		<< minFreqToSave << std::endl
		<< maxFreqToSave << std::endl
		<< stepFreqToSave << std::endl;

	stream
		<< powrThreshFract << std::endl;

	stream << elems.size() << std::endl;
	for (std::map<std::string, std::vector<std::string> >::const_iterator
		it = elems.begin() ; it != elems.end() ; ++it)
	{
		stream << it->first << std::endl;
		stream << it->second.size() << std::endl;
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			stream << it->second[i] << std::endl;
	}

	return stream.good();
}

//********************************************************************//
//************** S T O C H A S T I C   R E S O N A N C E *************//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
StochResMetric::StochResMetric(ParamHandler & )
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
StochResMetric::StochResMetric(std::ifstream & stream)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
StochResMetric::~StochResMetric()
{

}

//**********************************************************************
// Compute the activated cells
//**********************************************************************
bool StochResMetric::ComputeMetric(const ChIModel & model)
{
	ActivatedCells *actCellsMetr = 
		GetSpecificMetric<Metric, ActivatedCells>(model.GetAllMetrics());
	std::map<std::string, double> res = actCellsMetr->GetScalarStatsToSave();

	const ChIModelStochasticRes *tmp = dynamic_cast<const ChIModelStochasticRes *>(&model);
	assert(tmp);
	std::string sigStrat = tmp->getSignalStimStartName();
	std::string noiseStrat = tmp->getNoiseStimStartName();

	if (model.IsStimStratActivated(sigStrat))
	{
		if (model.IsStimStratActivated(noiseStrat))
		{
			sigAndNoise_act = res["TotalNbCellsActivated"];
			sigAndNoise_freq = res["MeanActivationFrequency"];
		}
		else
		{
			sigOnly_act = res["TotalNbCellsActivated"];
			sigOnly_freq = res["MeanActivationFrequency"];
		}
	}
	else if (model.IsStimStratActivated(noiseStrat))
	{
		noiseOnly_act = res["TotalNbCellsActivated"];
		noiseOnly_freq = res["MeanActivationFrequency"];
	}
	TRACE("sigAndNoise " << sigAndNoise_act << " // " << sigAndNoise_freq)
	TRACE("sigOnly " << sigOnly_act << " // " << sigOnly_freq)
	TRACE("NoiseOnly " << noiseOnly_act << " // " << noiseOnly_freq)
	
	return true;
}

//**********************************************************************
// Save the activated cells
//**********************************************************************
bool StochResMetric::SaveMetric(ResultSaver ) const
{
	bool allSaved = true;

	return allSaved;
}

//**********************************************************************
// Initialize the metric as if it were just created
//**********************************************************************
void StochResMetric::Initialize()
{
	Metric::Initialize();
	sigAndNoise_freq = 0;
	sigAndNoise_act = 0;
	noiseOnly_freq = 0;
	noiseOnly_act = 0;
	sigOnly_freq = 0;
	sigOnly_act = 0;
}

//**********************************************************************
// Builds a copy of the object (only the parameters are identical)
//**********************************************************************
Metric * StochResMetric::BuildCopy() const
{
	return new StochResMetric(); 
}

//**********************************************************************
// Returns the scalar values to be saved by the scalar stats metric
//  (cf SimulationMetrics.h)
//**********************************************************************
std::map<std::string, double> StochResMetric::GetScalarStatsToSave() const
{
	std::map<std::string, double> res;
	res["SigAndNoise_Act"] = sigAndNoise_act;
	res["NoiseOnly_Act"] = noiseOnly_act;
	res["SigOnly_Act"] = sigOnly_act;

	res["SigAndNoise_Freq"] = sigAndNoise_freq;
	res["NoiseOnly_Freq"] = noiseOnly_freq;
	res["SigOnly_Freq"] = sigOnly_freq;

	res["Chi_Act"] = (sigOnly_act != 0) ? 
		(sigAndNoise_act - noiseOnly_act) / sigOnly_act : DEFAULT_MAX_VAL;
	res["Chi_Freq"] = (sigOnly_freq != 0) ? 
		(sigAndNoise_freq - noiseOnly_freq) / sigOnly_freq : DEFAULT_MAX_VAL;
	TRACE("Real Values !!! : " << res["Chi_Act"] << " // " << res["Chi_Freq"])
	return res;
}

//**********************************************************************
// Loads the metric from a stream
//**********************************************************************
bool StochResMetric::LoadFromStream(std::ifstream & stream)
{
	Initialize();
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the metric to a stream
//**********************************************************************
bool StochResMetric::SaveToStream(std::ofstream & stream) const
{
	return stream.good();
}

