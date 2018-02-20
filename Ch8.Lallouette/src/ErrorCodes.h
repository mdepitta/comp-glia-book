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

#ifndef ERRORCODES_H
#define ERRORCODES_H

#define SUBGRID_SAVING_FAILED                 1
#define GRIDSEARCH_LOADING_FAILED             2
#define GRIDSEARCH_SAVING_FAILED              3
#define SIMULATION_SAVING_FAILED              4
#define MODEL_LOADING_FAILED                  5
#define PARAMETER_PARSING_FAILED              6
#define PARAMETER_LOADING_FROM_FILE_FAILED    7

#define SIMULATION_METRIC_COMPUTATION_PROBLEM 8
#define SIMULATION_METRIC_SAVING_PROBLEM      9
#define GRIDSEARCH_METRIC_COMPUTATION_PROBLEM 10
#define GRIDSEARCH_METRIC_SAVING_PROBLEM      11
#define MODEL_METRIC_SAVING_PROBLEM           12

#endif
