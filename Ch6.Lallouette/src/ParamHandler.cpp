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

#include "ParamHandler.h"
#include "utility.h"
#include <string.h>
#include <fstream>

using namespace std;

ParamHandler ParamHandler::GlobalParams;

template <> void AttributeValue<vector<int> >(stringstream & str, vector<int> & ref)
{
	int val;
	str >> val;
	ref.push_back(val);
}

template <> void AttributeValue<vector<double> >(stringstream & str, vector<double> & ref)
{
	double val;
	str >> val;
	ref.push_back(val);
}

template <> void AttributeValue<vector<string> >(stringstream & str, vector<string> & ref)
{
	string val;
	str >> val;
	ref.push_back(val);
}

template <> void AttributeValue<vector<bool> >(stringstream & str, vector<bool> & ref)
{
	bool val;
	str >> val;
	ref.push_back(val);
}

template <> std::string ValueToString<std::vector<int> >(std::vector<int> & ref)
{
	std::stringstream result;
	for (unsigned int i = 0 ; i < ref.size() ; ++i)
		result << (i ? "|" : "") << ref[i];
	return result.str();
}

template <> std::string ValueToString<std::vector<double> >(std::vector<double> & ref)
{
	std::stringstream result;
	for (unsigned int i = 0 ; i < ref.size() ; ++i)
		result << (i ? "|" : "") << ref[i];
	return result.str();
}
template <> std::string ValueToString<std::vector<std::string> >(std::vector<std::string> & ref)
{
	std::stringstream result;
	for (unsigned int i = 0 ; i < ref.size() ; ++i)
		result << (i ? "|" : "") << ref[i];
	return result.str();
}
template <> std::string ValueToString<std::vector<bool> >(std::vector<bool> & ref)
{
	std::stringstream result;
	for (unsigned int i = 0 ; i < ref.size() ; ++i)
		result << (i ? "|" : "") << ref[i];
	return result.str();
}

//**********************************************************************
// Default copy constructor
//**********************************************************************
ParamHandler::ParamHandler(const ParamHandler &h) : currParamName(h.currParamName)
{
	vector<ReferenceHolder*> tempVect;
	for (ParamType::const_iterator it = h.parameters.begin() ; it != h.parameters.end() ; ++it)
	{
		tempVect.clear();
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			tempVect.push_back(it->second[i]->BuildCopy());
		parameters.insert(make_pair(it->first, tempVect));
	}
}

//**********************************************************************
//**********************************************************************
ParamHandler::~ParamHandler()
{
	for (ParamType::iterator i = parameters.begin() ; i != parameters.end() ; ++i)
		for (vector<ReferenceHolder*>::iterator j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			delete *j;
		}
}

//**********************************************************************
//**********************************************************************
void ParamHandler::RemoveParam(std::string name)
{
	parameters.erase(name);
}

//**********************************************************************
//**********************************************************************
void ParamHandler::ChangeRefHolder(std::string name, unsigned int ind, 
	ReferenceHolder* r, bool delOld)
{
	if (parameters[name].size() > ind)
	{
		if (delOld and parameters[name][ind])
			delete parameters[name][ind];
		parameters[name][ind] = r;
	}
}

//**********************************************************************
//**********************************************************************
bool ParamHandler::Parse(int argc, char *argv[])
{
	ParamType::iterator it;
	bool ok = true;
	for (int i = 1 ; i < argc ; ++i)
	{
		if (string(argv[i]) == "-help")
		{
			PrintParams();
			return false;
		}
		if ((it = parameters.find(string(argv[i]))) != parameters.end())
		{
			for (vector<ReferenceHolder*>::iterator j = it->second.begin() ; j != it->second.end() ; ++j)
			{
				if ((j == it->second.begin()) and (dynamic_cast<SpecialReferenceHolder<bool>*>(*j) or 
					dynamic_cast<SpecialReferenceHolder<std::vector<bool> >*>(*j)))
				{
					it->second.front()->GiveVal("1");
					continue;
				}
				if (++i >= argc)
					break;
				else if (parameters.find(string(argv[i])) != parameters.end())
				{
					--i;
					break;
				}
				else
					(*j)->GiveVal(argv[i]);
			}
		}
		else
		{
			cerr << "Unknown parameter : " << argv[i] << endl;
			ok = false;
		}
	}
	return ok;
}

//**********************************************************************
//**********************************************************************
bool ParamHandler::SetVal(std::string name, const char *val, unsigned int ind)
{
	ParamType::iterator it;
	if (((it = parameters.find(name)) != parameters.end()) and (ind < it->second.size()))
	{
		it->second[ind]->GiveVal(val);
		return true;
	}
	else
		return false;
}

//**********************************************************************
// Loads parameters from a file
//**********************************************************************
bool ParamHandler::LoadParams(std::string path)
{
	ifstream stream(path.c_str());
	string tmpStr;
	vector<string> newParams;

	newParams.push_back("exeName");
	stream >> tmpStr;
	do
	{
		newParams.push_back(tmpStr);
		stream >> tmpStr;
	} while (not stream.eof());

	int argc = newParams.size();
	char **argv = new char*[argc];
	for (unsigned int i = 0 ; i < newParams.size() ; ++i)
	{
		char *newStr = new char[newParams[i].size() + 1];
		strcpy(newStr, newParams[i].c_str());
		argv[i] = newStr;
	}

	bool ok = Parse(argc, argv);

	for (unsigned int i = 0 ; i < newParams.size() ; ++i)
		delete[] argv[i];
	delete[] argv;

	return ok;
}

//**********************************************************************
//**********************************************************************
void ParamHandler::PrintParams(std::ostream & stream) const
{
	stream << "Parameter : type [Allowed values] (Default value), ..." << endl;
	for (ParamType::const_iterator i = parameters.begin() ; i != parameters.end() ; ++i)
	{
		stream << i->first << " : ";
		bool first = true;
		for (vector<ReferenceHolder*>::const_iterator j = i->second.begin() ; j != i->second.end() ; ++j)
		{
			if (*j)
			{
				stream << (first ? "" : ", ") 
					<< (*j)->TypeName();
				if ((*j)->HasAllowedVals()) // if allowed values specified
					stream << " [" << (*j)->AllowedValsToString() << "]";
				stream << " (" << (*j)->ValToString() << ")";
				first = false;
			}
		}
		stream << endl;
	}
}

//**********************************************************************
//**********************************************************************
std::string ParamHandler::getStringParam(std::string name, unsigned int ind) const
{
	ParamType::const_iterator it;
	if ((it = parameters.find(name)) != parameters.end())
		if(ind < it->second.size())
			return it->second[ind]->ValToString();
	return "";
}

//**********************************************************************
//**********************************************************************
ParamHandler & ParamHandler::operator<=(std::string name)
{
	currParamName = name;
	parameters.insert(make_pair(name, vector<ReferenceHolder*>()));
	return *this;
}

//**********************************************************************
//**********************************************************************
ParamHandler & ParamHandler::operator+=(const ParamHandler &h)
{
	for (ParamType::const_iterator it = h.parameters.begin() ; it != h.parameters.end() ; ++it)
	{
		if (parameters.find(it->first) == parameters.end())
			parameters.insert(make_pair(it->first, vector<ReferenceHolder*>()));
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			parameters[it->first].push_back(it->second[i]->BuildCopy());
	}
	return *this;
}

//**********************************************************************
//**********************************************************************
ParamHandler & ParamHandler::AddAllWithPrefix(const ParamHandler &h, std::string prefix)
{
	for (ParamType::const_iterator it = h.parameters.begin() ; it != h.parameters.end() ; ++it)
	{
		if (parameters.find(prefix + it->first) == parameters.end())
			parameters.insert(make_pair(prefix + it->first, vector<ReferenceHolder*>()));
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			parameters[prefix + it->first].push_back(it->second[i]->BuildCopy());
	}
	return *this;
}

//**********************************************************************
//**********************************************************************
std::vector<std::string> ParamHandler::GetParamNames() const
{
	std::vector<std::string> names;
	for (ParamType::const_iterator it = parameters.begin() ; 
			it != parameters.end() ; ++it)
		names.push_back(it->first);
	return names;
}

//**********************************************************************
// Adds a prefix to the parameter names
//**********************************************************************
ParamHandler ParamHandler::AddPrefix(std::string pre)
{
	ParamHandler tmpParams;
	for (ParamType::iterator it = parameters.begin() ; 
		it != parameters.end() ; ++it)
	{
		for (unsigned int i = 0 ; i < it->second.size() ; ++i)
			tmpParams.parameters[pre + it->first].push_back(it->second[i]->BuildCopy());
	}
	return tmpParams;
}

//**********************************************************************
// Returns the number of parameters associated to a given name
//**********************************************************************
unsigned int ParamHandler::GetNbParamsForName(const std::string & name) const
{
	ParamType::const_iterator it;
	if (((it = parameters.find(name)) != parameters.end()))
		return it->second.size();
	else
		return 0;
}
