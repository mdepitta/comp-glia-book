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

#include "ChISimulationManager.h"
#include "ChIModel.h"
#include "PropagationModels.h"
#include "NeuronNetModels.h"
#include "AstroNeuroModel.h"

#include "SimulationMetrics.h"
#include "ODESolvers.h"

using namespace AstroModel;
using namespace std;

//********************************************************************//
//************ C H I   S I M U L A T I O N   M A N A G E R ***********//
//********************************************************************//

template <> string RepeatSimulation<ChIModel>::ClassName("RepeatSimulationChIModel");
template <> string RepeatSimulation<KChIModel>::ClassName("RepeatSimulationKChIModel");
template <> string RepeatSimulation<ThresholdModel>::ClassName("RepeatSimulationThresholdModel");
template <> string RepeatSimulation<ModifiedThresholdModel>::ClassName("RepeatSimulationModifiedThresholdModel");
template <> string RepeatSimulation<SecondOrdNeighbSERSModel>::ClassName("RepeatSimulationSecondOrdNeighbSERSModel");
template <> string RepeatSimulation<SimpleThresholdSERSModel>::ClassName("RepeatSimulationSimpleThresholdSERSModel");
template <> string RepeatSimulation<ChIModelThresholdDetermination >::ClassName("RepeatSimulationChIModelThresholdDetermination");
template <> string RepeatSimulation<ChIModelStochasticRes>::ClassName("RepeatSimulationChIModelStochasticResonance");
template <> string RepeatSimulation<FireDiffuseModel>::ClassName("RepeatSimulationFireDiffuseModel");

template <> string RepeatSimulation<NeuronNetModel>::ClassName("RepeatSimulationNeuronNetModel");
template <> string RepeatSimulation<AstroNeuroNetModel>::ClassName("RepeatSimulationAstroNeuronNetModel");
