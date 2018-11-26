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

#include "MetricComputeStrat.h"
#include "AbstractFactory.h"
#include "utility.h"

using namespace AstroModel;
using namespace std;

//********************************************************************//
//**************************** M E T R I C ***************************//
//********************************************************************//

//**********************************************************************
// Add a dependency between two metric classes
//**********************************************************************
void Metric::AddMetricDependency(std::string m1, std::string m2)
{
	// Check that the new dependency doesn't create a cyclic dependency
	// Abort if cyclic dependency is detected
	assert(not CompareMetricsNames(m1, m2));

	// Otherwise, add the dependency
	Dependencies()[m1].insert(m2);
}

//**********************************************************************
//**********************************************************************
ParamHandler Metric::BuildModelParamHandler()
{
	ParamHandler params;
	return params;
}

//**********************************************************************
//**********************************************************************
void Metric::AddOptStatParamName(std::string str)
{
	optStatParamAddName += str;
}

//**********************************************************************
// Compare two given metrics names and return wether m2 depends on m1
//**********************************************************************
bool Metric::CompareMetricsNames(std::string m1, std::string m2)
{
	set<string> m2TotDeps, tempDeps, tempDeps2;
	tempDeps.insert(m2);
	do
	{
		tempDeps2.clear();
		for (set<string>::iterator it = tempDeps.begin() ; it != tempDeps.end() ; ++it)
			tempDeps2.insert(Dependencies()[*it].begin(), Dependencies()[*it].end());
		tempDeps = tempDeps2;
		m2TotDeps.insert(tempDeps.begin(), tempDeps.end());
	} while(tempDeps2.size() > 0);

	return (m2TotDeps.find(m1) != m2TotDeps.end());
}

const std::set<std::string> & Metric::GetDependencies(std::string mName)
{
	return Dependencies()[mName];
}

std::map<std::string, std::set<std::string> > & Metric::Dependencies()
{
	static std::map<std::string, std::set<std::string> > ans;
	return ans;
}

void Metric::AddSavedFile(std::string fileName, std::string path) const
{ 
	TRACE("Saving file : " << path) 
	savedFilePaths.push_back(std::make_pair(fileName, path)); 
}

bool Metric::HasFileBeenSaved(std::string fileName) const
{
	bool saved = false;
	for (std::vector<std::pair<std::string, std::string> >::const_iterator it = savedFilePaths.begin() ; it != savedFilePaths.end() ; ++it)
		saved |= (it->first == fileName);
	return saved;
}

std::string Metric::ModifStatName(std::string str) const
{
	return optStatParamAddName + str;
}

