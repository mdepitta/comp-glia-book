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

#ifndef PARAMHANDLER_H
#define PARAMHANDLER_H

#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <typeinfo>

// Forward declarations

template <typename T> void AttributeValue(std::stringstream & str, T & ref)
{
	str >> ref;
}
template <> void AttributeValue<std::vector<int> >(std::stringstream & str, std::vector<int> & ref);
template <> void AttributeValue<std::vector<double> >(std::stringstream & str, std::vector<double> & ref);
template <> void AttributeValue<std::vector<const double> >(std::stringstream & str, std::vector<const double> & ref);
template <> void AttributeValue<std::vector<std::string> >(std::stringstream & str, std::vector<std::string> & ref);
template <> void AttributeValue<std::vector<bool> >(std::stringstream & str, std::vector<bool> & ref);

template <typename T> std::string ValueToString(T & ref)
{
	std::stringstream result;
	result << ref;
	return result.str();
}
template <> std::string ValueToString<std::vector<int> >(std::vector<int> & ref);
template <> std::string ValueToString<std::vector<double> >(std::vector<double> & ref);
template <> std::string ValueToString<std::vector<std::string> >(std::vector<std::string> & ref);
template <> std::string ValueToString<std::vector<bool> >(std::vector<bool> & ref);

// Utility base class for ParamHandler
class ReferenceHolder
{
public:
	virtual ~ReferenceHolder() {}
	virtual void GiveVal(const char*) {}
	virtual void AddAllowedVal(const char*) {}
	virtual bool HasAllowedVals() const { return false; }
	virtual std::string TypeName() const { return ""; }
	virtual std::string ValToString() const { return ""; }
	virtual std::string AllowedValsToString() { return ""; }
	virtual ReferenceHolder * BuildCopy() const { return new ReferenceHolder(); }
};

// Utility class for ParamHandler
template <typename T> class SpecialReferenceHolder : public ReferenceHolder
{
protected:
	T & ref;
	std::vector<T> allowedVals;
public:
	SpecialReferenceHolder(const SpecialReferenceHolder<T> &srh) : ref(srh.ref), allowedVals(srh.allowedVals) {}
	SpecialReferenceHolder(T & _ref) : ref(_ref) {}
	virtual ~SpecialReferenceHolder()
	{
		allowedVals.clear();
	}

	virtual ReferenceHolder * BuildCopy() const { return new SpecialReferenceHolder<T>(*this); }
	virtual void GiveVal(const char* argv)
	{
		std::stringstream temp;
		temp << argv;
		AttributeValue<T>(temp, ref);
	}
	virtual T getVal() const
	{
		return ref;
	}
	virtual void AddAllowedVal(const char* argv)
	{
		T valTemp;
		std::stringstream temp;
		temp << argv;
		AttributeValue<T>(temp, valTemp);
		allowedVals.push_back(valTemp);
	}
	virtual bool HasAllowedVals() const { return not allowedVals.empty(); }
	virtual std::string TypeName() const { return typeid(ref).name(); }
	virtual std::string ValToString() const { return ValueToString<T>(ref); }
	virtual std::string AllowedValsToString() 
	{
		std::string result;
		for (unsigned int i = 0 ; i < allowedVals.size() ; ++i)
		{
			T valTemp = allowedVals[i];
			result += std::string((i ? " " : "")) + ValueToString<T>(valTemp);
		}
		return result;
	}
};

class ParamHandler
{
protected:
	typedef std::map<std::string, std::vector<ReferenceHolder*> > ParamType;
	ParamType parameters;

	std::string currParamName;
public:
	static ParamHandler GlobalParams;

	ParamHandler() : currParamName("") {}
	ParamHandler(const ParamHandler &h);
	virtual ~ParamHandler();

	void RemoveParam(std::string name);
	void ChangeRefHolder(std::string name, unsigned int ind, ReferenceHolder* r, bool delOld = false);

	bool Parse(int argc, char *argv[]);
	bool SetVal(std::string name, const char *val, unsigned int ind = 0);

	bool LoadParams(std::string path);
	void PrintParams(std::ostream & = std::cout) const;

	template <typename T> T getParam(std::string name, unsigned int ind = 0) const
	{
		ParamType::const_iterator it;
		SpecialReferenceHolder<T>* srh;
		if ((it = parameters.find(name)) != parameters.end())
			if(ind < it->second.size())
				if ((srh = dynamic_cast<SpecialReferenceHolder<T>*>(it->second[ind])))
					return srh->getVal();
		return T();
	}

	std::string getStringParam(std::string name, unsigned int ind = 0) const;

	ParamHandler & operator<=(std::string name);
	ParamHandler & operator+=(const ParamHandler &h);

	ParamHandler & AddAllWithPrefix(const ParamHandler &h, std::string prefix);

	template <typename T> ParamHandler & operator,(T & var)
	{
		ParamType::iterator it;
		if ((it = parameters.find(currParamName)) != parameters.end())
			it->second.push_back(new SpecialReferenceHolder<T>(var));
		return *this;
	}

	template <typename T> void AddAllowedVal(std::string name, unsigned int pos, T val)
	{
		ParamType::iterator it;
		if ((it = parameters.find(name)) != parameters.end() && (it->second.size() > pos))
			it->second[pos]->AddAllowedVal(ValueToString(val).c_str());
	}
	template <typename T> void AddAllowedValsList(std::string name, unsigned int pos, std::vector<T> vals)
	{
		for (unsigned int i = 0 ; i < vals.size() ; ++i)
			AddAllowedVal(name, pos, vals[i]);
	}

	std::vector<std::string> GetParamNames() const;

	// Adds a prefix to the parameter names
	ParamHandler AddPrefix(std::string pre);
	// Returns the number of parameters associated to a given name
	unsigned int GetNbParamsForName(const std::string & name) const;
};

#endif

