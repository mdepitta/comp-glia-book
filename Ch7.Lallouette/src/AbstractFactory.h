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

#ifndef ABSTRACTFACTORY_H
#define ABSTRACTFACTORY_H

#include <fstream>
#include <map>
#include <vector>
#include <string>

#include "ParamHandler.h"

namespace AstroModel
{
	// Forward declarations
	template <typename T>
	class AbstractFactory
	{
	public:
		virtual T * Create() = 0;
		virtual T * CreateFromStream(std::ifstream & stream) = 0;
		static std::map<std::string, AbstractFactory<T>*> Factories;
		static std::vector<std::string> GetFactoriesNames() 
		{
			std::vector<std::string> result;
			for (typename std::map<std::string, AbstractFactory<T>*>::iterator it = Factories.begin() ; it != Factories.end() ; ++it)
				result.push_back(it->first);
			return result;
		}
	};

	template <typename T, typename S>
	class DerivedFactory : public AbstractFactory<T>
	{
	public:
		virtual S * Create() { return new S(); }
		virtual S * CreateFromStream(std::ifstream & stream) { return new S(stream); }
	};

}


#endif
