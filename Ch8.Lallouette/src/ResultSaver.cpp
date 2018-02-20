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

#include "ResultSaver.h"

#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <boost/filesystem/convenience.hpp>

using namespace std;

//********************************************************************//
//******************* R E S U L T   S A V E R ************************//
//********************************************************************//

ResultSaver::ParamTablesMap ResultSaver::paramsToSave;
std::map<std::string, bool> ResultSaver::headersWriten;
int ResultSaver::instCount = 0;

ofstream ResultSaver::voidStream;
ResultSaver ResultSaver::NullSaver;

//**********************************************************************
//**********************************************************************
bool ResultSaver::createDir(std::string path)
{
	return boost::filesystem::exists(path) ? boost::filesystem::is_directory(path) : boost::filesystem::create_directories(path);
}

//**********************************************************************
//**********************************************************************
ResultSaver::ResultSaver(bool _saving) : saving(_saving), path(""), name(""), ext("dat")
{
	++instCount;
}

//**********************************************************************
//**********************************************************************
ResultSaver::ResultSaver(string _path, string _name, string _ext) : saving(true), path(_path), name(_name), ext(_ext)
{
	++instCount;
	if (not createDir(path))
		cerr << "Couldn't create directory : " << path << endl;
}

//**********************************************************************
//**********************************************************************
ResultSaver::ResultSaver(const ResultSaver & rs) : 
	saving(rs.saving), path(rs.path), name(rs.name), ext(rs.ext), toSave(rs.toSave)
{
	++instCount;
}

//**********************************************************************
//**********************************************************************
ResultSaver & ResultSaver::operator=(const ResultSaver & rs)
{
	saving = rs.saving;
	path = rs.path;
	name = rs.name;
	ext = rs.ext;
	toSave = rs.toSave;
	++instCount;
	return *this;
}

//**********************************************************************
//**********************************************************************
ResultSaver::~ResultSaver()
{
	--instCount;
	if (instCount == 0)
		flushTableFiles();
	if (stream.is_open())
		stream.close();
}

//**********************************************************************
//**********************************************************************
std::ofstream & ResultSaver::getStream() 
{
	if ((not stream.is_open()) and (path != ""))
	{
		if (not createDir(path))
		{
			cerr << "Couldn't create directory : " << path << endl;
			return voidStream;
		}
		else
			stream.open((path + "/" + name + "." + ext).c_str(), ios_base::app);
	}
	return saving ? stream : voidStream;
}

//**********************************************************************
//**********************************************************************
ResultSaver ResultSaver::operator()(string subDir)
{
	ResultSaver temp(*this);
	temp.path += "/" + subDir;
	if (temp.stream.is_open())
		temp.stream.close();
	return temp;
}

//**********************************************************************
//**********************************************************************
bool ResultSaver::isSaving(string objName)
{
	if (toSave.find(objName) != toSave.end())
	{
		if ((objName != name) and stream.is_open())
			stream.close();
		name = objName;
		return true;
	}
	else
	{
		bool found = false;
		for (ParamTablesMap::iterator it = paramsToSave.begin() ; it != paramsToSave.end() and (not found) ; ++it)
			found = (it->second.second.find(objName) != it->second.second.end());
		return found;
	}
}

//**********************************************************************
//**********************************************************************
ResultSaver & ResultSaver::operator,(string objToSave)
{
	toSave.insert(objToSave);
	lastAddedName = objToSave;
	return *this;
}

//**********************************************************************
//**********************************************************************
ResultSaver & ResultSaver::operator|(string paramToSave)
{
	ParamTablesMap::iterator it;
	if ((it = paramsToSave.find(lastAddedName)) != paramsToSave.end())
		it->second.second.insert(pair<string, ValueHolder*>(paramToSave, 0));
	else
	{
		paramsToSave.insert(pair<string, ParamTable>(lastAddedName,ParamTable(path, map<string, ValueHolder*>())));
		paramsToSave[lastAddedName].second.insert(pair<string, ValueHolder*>(paramToSave, 0));
		headersWriten.insert(make_pair(lastAddedName, false));
	}
	return *this;
}

//**********************************************************************
//**********************************************************************
void ResultSaver::flushTableFiles()
{
	for (ParamTablesMap::iterator i = paramsToSave.begin() ; i != paramsToSave.end() ; ++i)
	{
		if (stream.is_open() and (name != i->first))
			stream.close();
		name = i->first;
		std::string oldPath(path);
		path = i->second.first;

		flushLine(i->first, true);

		path = oldPath;
	}
}

//**********************************************************************
//**********************************************************************
void ResultSaver::flushLine(std::string tableName, bool flushAll)
{
	map<string, ValueHolder*> &map = paramsToSave[tableName].second;
	// Write headers
	if (not headersWriten[tableName] and (path != ""))
	{
		for (std::map<string, ValueHolder*>::iterator j = map.begin() ; j != map.end() ; ++j)
		{
			if (j != map.begin())
				getStream() << "\t";
			getStream() << j->first;
		}
		getStream() << endl;
		headersWriten[tableName] = true;
	}
	// Write data
	for (std::map<string, ValueHolder*>::iterator j = map.begin() ; j != map.end() ; ++j)
	{
		if (j->second)
		{
			if (path != "")
			{
				if (j != map.begin())
					getStream() << "\t";
				j->second->ToStream(getStream());
			}
			if (flushAll or not(j->second->keepValue))
			{
				delete j->second;
				j->second = 0;
			}
		}
		else if (path != "")
		{
			if (j != map.begin())
				getStream() << "\t";
			getStream() << "NA";
		}
	}
	if (path != "")
		getStream() << std::endl;
}

//**********************************************************************
// Loads the cell from a stream
//**********************************************************************
bool ResultSaver::LoadFromStream(std::ifstream & stream)
{
	unsigned int nbObj = 0;
	stream >> nbObj;
	std::string tempStr;
	for (unsigned int i = 0 ; i < nbObj ; ++i)
	{
		stream >> tempStr;
		toSave.insert(tempStr);
	}
	return stream.good() and not stream.eof();
}

//**********************************************************************
// Saves the cell to a stream
//**********************************************************************
bool ResultSaver::SaveToStream(std::ofstream & stream) const
{
	stream << toSave.size() << std::endl;
	for (std::set<std::string>::const_iterator it = toSave.begin() ; it != toSave.end() ; ++it)
		stream << *it << std::endl;
	return stream.good();
}

