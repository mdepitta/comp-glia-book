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

#ifndef RESULTSAVER_H
#define RESULTSAVER_H

#include "Savable.h"
#include "utility.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>

#define SYS_OK_CODE 0

// Utility class for table saving
class ValueHolder
{
public:
	bool keepValue;

	ValueHolder(bool _keepVal = false) : keepValue(_keepVal) {}
	virtual void ToStream(std::ofstream &) {};

	virtual ~ValueHolder() {}
};

// Utility class for table saving
template <typename T> class SpecialValueHolder : public ValueHolder
{
protected:
	T value;
public:
	SpecialValueHolder(const T & _val, bool _keepVal = false) : 
		ValueHolder::ValueHolder(_keepVal), value(_val) {}
	virtual void ToStream(std::ofstream & stream)
	{
		stream << value;
	}
	virtual T GetVal() const
	{
		return value;
	}
};

// Utility class to facilitate the use of indexed parameters
template<typename T> class IndexedParam
{
protected:
	std::string name;
	std::vector<T> indexes;
	unsigned int w;
public:
	IndexedParam(std::string _name, T start, T end, T step = 0, unsigned int _w = 4) : 
		name(_name), w(_w)
	{
		for (T i = start ; i <= end ; i += step)
			indexes.push_back(i);
	}
	std::string GetParamName()
	{
		if (indexes.empty())
			return "";
		std::stringstream temp;
		temp << name << StringifyFixed(indexes.back(), w, '0');
		indexes.pop_back();
		return temp.str();
	}
};

class ResultSaver : public SaveAndLoadFromStream
{
public:
	typedef std::pair<std::string /*table file path*/, std::map<std::string /*column name*/, ValueHolder*> > ParamTable;
	typedef std::map<std::string /*table name*/, ParamTable> ParamTablesMap;

protected:
	bool saving;
	std::string path;
	std::string name;
	std::string ext;
	std::ofstream stream;

	std::string lastAddedName;

	std::set<std::string> toSave;

	static ParamTablesMap paramsToSave;
	static std::map<std::string, bool> headersWriten;
	static int instCount;

	static std::ofstream voidStream;

	void flushTableFiles();
	void flushLine(std::string tableName, bool flushAll = false);

	static bool createDir(std::string path);
public:
	static ResultSaver NullSaver;

	ResultSaver(bool _saving = false);
	ResultSaver(std::string _path, std::string _name = "default_file", std::string _ext = "dat");
	ResultSaver(const ResultSaver & rs);

	ResultSaver & operator=(const ResultSaver & rs);

	virtual ~ResultSaver();

	inline operator bool() const { return saving; }

	std::ofstream & getStream();
	std::string getCurrFile() const { return path + "/" + name + "." + ext; }
	inline const std::string & GetCurrPath() const { return path; }
	ResultSaver operator()(std::string subDir);
	inline ResultSaver operator()(const Savable *sav) { return operator()(sav->getName()); };

	bool isSaving(std::string objName);

	ResultSaver & operator,(std::string objToSave);
	ResultSaver & operator|(std::string);
	template <typename T> ResultSaver & operator|(IndexedParam<T> paramToSave)
	{
		std::string temp;
		while ((temp = paramToSave.GetParamName()) != "")
			this->operator|(temp);
		return *this;
	}
	inline ResultSaver & operator<=(std::string objToSave) { return operator,(objToSave); }

	void AddAttributesToTable(std::string destTable, std::string srcTable)
	{
		ParamTablesMap::const_iterator it;
		if ((it = paramsToSave.find(srcTable)) != paramsToSave.end())
		{
			for (std::map<std::string , ValueHolder*>::const_iterator it2 = it->second.second.begin() ; 
				it2 != it->second.second.end() ; ++it2)
				(*this <= destTable) | it2->first;
		}
	}

	void KeepValue(std::string tableName, std::string paramName, bool keepVal)
	{
		ParamTablesMap::iterator i;
		if ((i = paramsToSave.find(tableName)) != paramsToSave.end())
		{
			std::map<std::string, ValueHolder*>::iterator it;
			if ((it = i->second.second.find(paramName)) != i->second.second.end())
				it->second->keepValue = keepVal;
		}
	}

	template <typename T> void Param(std::string paramName, const T &val)
	{
		for (ParamTablesMap::iterator i = paramsToSave.begin() ; i != paramsToSave.end() ; ++i)
		{
			std::map<std::string, ValueHolder*>::iterator it;
			if ((it = i->second.second.find(paramName)) != i->second.second.end())
			{
				if (stream.is_open() and (name != i->first))
					stream.close();
				name = i->first;
				std::string oldPath(path);
				path = i->second.first;

				if (((path != "") and saving) or ((path == "") and (not saving)))
				{
					if (not it->second)
						it->second = new SpecialValueHolder<T>(val);
					else // If the value was already writen
					{
						bool allWriten = true;
						for (std::map<std::string, ValueHolder*>::iterator j = i->second.second.begin() ; j != i->second.second.end() ; ++j)
						{
							allWriten = allWriten and j->second;
						}
						if (allWriten) // If all values have already been writen
						{
							flushLine(i->first);
							it->second = new SpecialValueHolder<T>(val);
						}
					}
				}
				path = oldPath;
			}
		}
	}

	double GetParamVal(std::string tableName, std::string paramName) const
	{
		ParamTablesMap::const_iterator it = paramsToSave.find(tableName);
		if (it != paramsToSave.end())
		{
			std::map<std::string, ValueHolder*>::const_iterator paramIt = it->second.second.find(paramName);
			SpecialValueHolder<double>* spVh = 0;
			if ((paramIt != it->second.second.end()) and (spVh = dynamic_cast<SpecialValueHolder<double>*>(paramIt->second)))
				return spVh->GetVal();
		}
		return 0;
	}

	//===========================================================||
	// Standard Save and Load methods                            ||
	//===========================================================||
	// Loads the cell from a stream
	virtual bool LoadFromStream(std::ifstream & stream);
	// Saves the cell to a stream
	virtual bool SaveToStream(std::ofstream & stream) const;

};

#endif

