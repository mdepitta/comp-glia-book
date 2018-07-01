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

#ifndef METRICCOMPUTESTRAT_H
#define METRICCOMPUTESTRAT_H

#include "ResultSaver.h"
#include "Savable.h"
#include "AbstractFactory.h"
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

// Dependency declaration macro
#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)

#define UNIQUE_ID_NAME(name, count, file) TOKENPASTE2(TOKENPASTE2(name, count), file)

#define ADD_METRIC_DEPENDENCY2(count, line, str1, str2) class UNIQUE_ID_NAME(UniqueMetrClass,count,line) {public: UNIQUE_ID_NAME(UniqueMetrClass,count,line)() {Metric::AddMetricDependency(str1, str2);}  }; UNIQUE_ID_NAME(UniqueMetrClass,count,line) UNIQUE_ID_NAME(UniqueMetrObj,count,line);

#define ADD_METRIC_DEPENDENCY(str1, str2) ADD_METRIC_DEPENDENCY2(__COUNTER__, __LINE__, str1, str2)

namespace AstroModel
{
/**********************************************************************/
/* Abstract class                                                     */
/**********************************************************************/
	class Metric : public SaveAndLoadFromStream
	{
	public:
		virtual ~Metric() {}
		virtual std::string GetClassName() const = 0;
		virtual std::vector<std::pair<std::string, std::string> > PullSavedFilePaths() 
		{
			std::vector<std::pair<std::string, std::string> > temp = savedFilePaths;
			savedFilePaths.clear();
			return temp; 
		}
		virtual void Initialize() { }
		virtual Metric * BuildCopy() const = 0;
		virtual ParamHandler BuildModelParamHandler();
		virtual std::map<std::string, double> GetScalarStatsToSave() const 
			{ return std::map<std::string, double>(); }
		virtual void AddOptStatParamName(std::string str);
		virtual std::string ModifStatName(std::string str) const;

		static void AddMetricDependency(std::string m1, std::string m2);
		static bool CompareMetricsNames(std::string m1, std::string m2);
		static const std::set<std::string> & GetDependencies(std::string mName);
	
	protected:
		mutable std::vector<std::pair<std::string, std::string> > savedFilePaths;
		std::string optStatParamAddName;

		static std::map<std::string, std::set<std::string> > & Dependencies();
		virtual void AddSavedFile(std::string fileName, std::string path) const;
		virtual bool HasFileBeenSaved(std::string fileName) const;
	};

/**********************************************************************/
/* Empty classes for multiple inheritance                             */
/**********************************************************************/
	class NeedFrequentUpdateMetric {};
	class StaticMetric {};
	class AfterSimMetric {};
	class DynMetric : public NeedFrequentUpdateMetric {};

/**********************************************************************/
/* Abstract class                                                     */
/**********************************************************************/
	template <typename ObjectT>
	class SpecificMetric : public Metric
	{
	public:
		virtual ~SpecificMetric() {}
		virtual bool ComputeMetric(const ObjectT & obj) = 0;
		virtual bool SaveMetric(ResultSaver saver) const = 0;
	};

/**********************************************************************/
/* Sorted Metric vector                                               */
/**********************************************************************/
	template <typename ObjectT>
	class SortedMetrics : public std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >, SaveAndLoadFromStream
	{
	public:
		// Destructor
		virtual ~SortedMetrics()
		{
			FreeAndClean();
		}

		// Metrics order comparison functor
		class CompareMetrics
		{
		public:
			bool operator()(const std::pair<SpecificMetric<ObjectT> *, bool> & m1, const std::pair<SpecificMetric<ObjectT> *, bool> & m2)
			{
				return Metric::CompareMetricsNames(m1.first->GetClassName(), m2.first->GetClassName());
			}
		};

		// Adds the given metric to the sorted vector and adds all its dependencies
		// If the given metric is already inside it is freed (given that _f is true)
		virtual bool AddMetricAndDependencies(Metric * _m, bool _f, ObjectT * caller)
		{
			SpecificMetric<ObjectT> * metr = dynamic_cast<SpecificMetric<ObjectT> *>(_m);
			if (metr)
			{
				bool ok = true;
				if (not IsAlreadyThere(metr))
				{
					this->push_back(std::make_pair(metr, _f));
TRACE(_m->GetClassName() << " successfully added to " << caller->GetClassName())
					const std::set<std::string> & deps = Metric::GetDependencies(metr->GetClassName());
					for (std::set<std::string>::const_iterator it = deps.begin() ; it != deps.end() ; ++it)
					{
						assert(AbstractFactory<Metric>::Factories[*it]);
						ok &= caller->AddMetric(AbstractFactory<Metric>::Factories[*it]->Create(), true);
					}
				}
				else if (_f)
					delete metr;
				return ok;
			}
			else
				return false;
		}

		// Adds all the dependencies from the metrics currently added
		virtual void AddDependenciesToCaller(ObjectT *caller)
		{
			for (typename std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >::const_iterator it = this->begin() ; 
					(it != this->end()) ; ++it)
			{
				for (std::set<std::string>::const_iterator it2 =
					Metric::GetDependencies(it->first->GetClassName()).begin() ; 
					it2 != Metric::GetDependencies(it->first->GetClassName()).end() ; ++it2)
				{
					assert(AbstractFactory<Metric>::Factories[*it2]);
					caller->AddMetric(AbstractFactory<Metric>::Factories[*it2]->Create(), true);
				}
			}
		}

		virtual const std::vector<Metric*> & GetMetricsRaw() const
		{
			return metricsRaw;
		}

		inline bool ComputeMetricsDefault(const ObjectT & caller) const
		{ return ComputeMetrics<SpecificMetric<ObjectT> >(caller); }

		// Compute all metrics of given type
		template <typename MetrType>
		bool ComputeMetrics(const ObjectT & caller) const
		{
			MetrType *ptr = 0;
			bool ok = true;
			for (typename std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >::const_iterator it = this->begin() ; 
					(it != this->end()) ; ++it)
			{
				if ((ptr = dynamic_cast<MetrType *>(it->first)))
				{
					ok &= it->first->ComputeMetric(caller);
				}
			}
			return ok;
		}

		bool SaveMetricsDefault(ResultSaver saver) const
		{ return SaveMetrics<SpecificMetric<ObjectT> >(saver); }

		// Save all metrics of given type
		template <typename MetrType>
		bool SaveMetrics(ResultSaver saver) const
		{
			MetrType *ptr = 0;
			bool ok = true;
			for (typename std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >::const_iterator it = this->begin() ; 
					(it != this->end()) ; ++it)
			{
				if ((ptr = dynamic_cast<MetrType *>(it->first)))
				{
					ok &= it->first->SaveMetric(saver);
				}
			}
			return ok;
		}

		void InitializeMetricsDefault() const
		{ InitializeMetrics<SpecificMetric<ObjectT> >(); }

		// Compute all metrics of given type
		template <typename MetrType>
		void InitializeMetrics() const
		{
			MetrType *ptr = 0;
			for (typename std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >::const_iterator it = this->begin() ; 
					(it != this->end()) ; ++it)
			{
				if ((ptr = dynamic_cast<MetrType *>(it->first)))
					ptr->Initialize();
			}
		}

		virtual void push_back(const std::pair<SpecificMetric<ObjectT> *, bool> obj)
		{
			if (not IsAlreadyThere(obj.first))
			{
				std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >::push_back(obj);
				metricsRaw.push_back(obj.first);
				std::sort(this->begin(), this->end(), CompareMetrics());
			}
			else if (obj.second)
				delete obj.first;
		}

		virtual bool IsAlreadyThere(const SpecificMetric<ObjectT> * metr) const
		{
			bool alreadyThere = false;
			for (typename std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >::const_iterator it = this->begin() ; 
					(it != this->end()) and not alreadyThere ; ++it)
				alreadyThere |= (it->first->GetClassName() == metr->GetClassName());
			return alreadyThere;
		}

		virtual bool LoadFromStream(std::ifstream & stream)
		{
			bool ok = true;
			unsigned int vectSize;
			std::string metrName;
			FreeAndClean();
			stream >> vectSize;
			for (unsigned int i = 0 ; i < vectSize ; ++i)
			{
				stream >> metrName;
TRACE(metrName)
				assert(AbstractFactory<Metric>::Factories[metrName]);
				Metric *metrTemp = AbstractFactory<Metric>::Factories[metrName]->CreateFromStream(stream);
				SpecificMetric<ObjectT> *metrSpe = dynamic_cast<SpecificMetric<ObjectT> *>(metrTemp);
				if (metrSpe)
					this->push_back(std::make_pair(metrSpe, true));
				else
				{
					ok = false;
					delete metrTemp;
				}
			}
			std::sort(this->begin(), this->end(), CompareMetrics());
				
			return ok and (stream.good() or stream.eof());
		}

		virtual bool SaveToStream(std::ofstream & stream) const
		{
			bool ok = true;
			stream << this->size() << std::endl;
			for (typename std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >::const_iterator it = this->begin() ; 
					(it != this->end()) ; ++it)
			{
				assert(it->first);
				stream << it->first->GetClassName() << std::endl;
				ok &= it->first->SaveToStream(stream);
			}
			return ok and stream.good();
		}

		virtual void FreeAndClean()
		{
			for (typename std::vector<std::pair<SpecificMetric<ObjectT> *, bool> >::const_iterator it = this->begin() ; 
					(it != this->end()) ; ++it)
			{
				if (it->second)
				{
					delete it->first;
				}
			}
			this->clear();
			metricsRaw.clear();
		}

		virtual ParamHandler BuildModelParamHandler()
		{
			ParamHandler params;
			for (typename SortedMetrics::iterator it = this->begin() ; 
					it != this->end() ; ++it)
				params += it->first->BuildModelParamHandler();
			return params;
		}

	protected:
		std::vector<Metric*> metricsRaw;
	};
}

#endif
