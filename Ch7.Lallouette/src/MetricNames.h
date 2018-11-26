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

#ifndef METRICNAMES_H
#define METRICNAMES_H

// PropagationMetrics
	#define BOOLWRP_BOOLNK_CELLSTATESAVERMETRIC "BoolWrp_BoolLnk_CellStateSaverMetric"
	#define SERSTATES_BOOLNK_CELLSTATESAVERMETRIC "SERStates_BoolLnk_CellStateSaverMetric"
// ChIModelMetrics
	#define ACTIVATED_CELLS_METRICS "ActivatedCellsMetric"
	#define CONCENTRATIONS_METRIC "ConcentrationsMetric"
	#define PROPAGATION_DISTANCE_METRIC "PropagationDistanceMetric"
	#define CORRELATIONS_METRIC "CorrelationsMetric"
	#define TRANSFER_ENTROPY_METRIC "TransferEntropyMetric"
	#define FUNCT_TOPO_BY_NET_STRAT_METRIC "FunctTopoByNetStratMetric"
	#define THRESHOLD_DETERMINATION_METRIC "ThresholdDeterminationMetric"
	#define WAVE_FRONT_DETECT_METRIC "WaveFrontDetectMetric"
	#define FOURIER_TRANSFORM_METRIC "FourierTransformMetric"
	#define WAVELET_TRANSFORM_METRIC "WaveletTransformMetric"
	#define STOCH_RES_METRIC "StochasticResonanceMetric"
// Network Metrics
	#define DEGREE_DISTRIBUTION_NETWORK_METRIC "DegreeDistributionNetworkMetric"
	#define CLUSTERING_COEFFICIENT_DISTRIBUTION_NETWORK_METRIC "ClusteringCoefficientDistributionNetworkMetric"
	#define ALL_PAIR_DISTANCES_NETWORK_METRIC "AllPairDistancesNetworkMetric"
	#define ADJACENCY_MATRIX_NETWORK_METRIC "AdjacencyMatrixNetworkMetric"
	#define SPATIAL_POSITION_NETWORK_METRIC "SpatialPositionsNetworkMetric"
	#define DIMENSION_NETWORK_METRIC "DimensionsNetworkMetric"
// Simulation Metrics
	#define PROPAGATION_DISTANCE_STAT_METRIC "PropagationDistanceStatMetric"
	#define NETWORK_TOPOLOGY_AND_SPATIAL_STAT_METRIC "NetworkTopologyAndSpatialStatMetric"
	#define SCALAR_STATISTICS_METRIC "ScalarStatisticsMetric"
	#define GRID_SEARCH_METRIC "GridSearchMetric"
// Stimulation Metrics
	#define STIMULATED_CELLS_METRIC "StimulatedCellsMetric"
	#define MEAN_PATH_STIM_METRIC "MeanPathStimulationMetric"

#endif
