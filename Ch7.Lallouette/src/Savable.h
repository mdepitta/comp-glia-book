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

#ifndef SAVABLE_H
#define SAVABLE_H

#include <string>
#include <fstream>

class ResultSaver;

class Savable
{
protected:
	std::string name;
	
public:
	Savable() : name("") {}
	Savable(std::string str) : name(str) {}
	virtual std::string getName() const { return name; }

	virtual void Save(std::string, ResultSaver) const;
};

class SaveAndLoadFromStream
{
public:
	virtual ~SaveAndLoadFromStream() {}
	virtual bool LoadFromStream(std::ifstream & stream) = 0;
	virtual bool SaveToStream(std::ofstream & stream) const = 0;
};

#endif
