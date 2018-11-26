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

#include "PropagationMetrics.h"

#include "CouplingFunction.h"
#include "PropagationModels.h"
#include "MetricNames.h"

using namespace std;
using namespace AstroModel;

template <> string CellStateSaver<BooleanWrapper, BooleanLink>::ClassName(BOOLWRP_BOOLNK_CELLSTATESAVERMETRIC);
template <> string CellStateSaver<SERStates, BooleanLink>::ClassName(SERSTATES_BOOLNK_CELLSTATESAVERMETRIC);

