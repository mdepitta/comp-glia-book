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

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <math.h>
#include <string>

#define HILL1(x,K) (x / (x + K))
#define HILL2(x,K) (x * x / (x * x + K * K))
#define HILLn(x, K, n) (pow(x, n) / (pow(x, n) + pow(K, n)))

#define INTPOW2(x) (x*x)
#define INTPOW3(x) (x*x*x)
#define INTPOW4(x) (x*x*x*x)

namespace AstroModel
{
	class ChICell;
	class KChICell;
	class ChIModel;
	class KChIModel;
	class Neuron;
	class Synapse;
	class TMSynapse;
	class TMSynapseOptim;
	class SFALIFNeuron;
	class DummyNeuron;
	class NeuronNetModel;
	class AstroNeuroNetModel;
	class FireDiffuseModel;
	class FireDiffuseCell;
}

namespace ODE
{
	/******************************************************************/
	/* Abstract class                                                 */
	/******************************************************************/
	template <typename Val = double, typename TimeT = double> class Function
	{
	public:
		typedef Val TVal;

		virtual void CompFunc(const TimeT & t, const Val *v, Val *f) const = 0;
		// Return class name
		virtual std::string GetClassName() const = 0;

		virtual ~Function() {}
	};
		
	/******************************************************************/
	/* Abstract Subclasses                                            */
	/******************************************************************/
	class SynapseFunct : public Function<double> {};
	class NeuronFunct : public Function<double> {};
	class AbstrChIFunct : public Function<double> {};
	class AbstrNetwChIFunct : public Function<double> {};

	/******************************************************************/
	/* ChI Cell Function                                              */
	/******************************************************************/
	class ChICellFunct : public AbstrChIFunct
	{
	protected:
		AstroModel::ChICell & cell;
		const AstroModel::ChIModel & model;

	public:
		static std::string ClassName;

		ChICellFunct(AstroModel::ChICell & _cell);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double &, const double *v, double *f) const;
	};

	/******************************************************************/
	/* KChI Cell Function                                             */
	/******************************************************************/
	class KChICellFunct : public AbstrChIFunct
	{
	protected:
		AstroModel::KChICell & cell;
		const AstroModel::KChIModel & model;

	public:
		static std::string ClassName;

		KChICellFunct(AstroModel::KChICell & _cell);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double &, const double *v, double *f) const;
	};

	/******************************************************************/
	/* Tsodyks Markram synapse Function                               */
	/******************************************************************/
	class TMSynapseFunct : public SynapseFunct
	{
	protected:
		AstroModel::TMSynapse & synapse;
		const AstroModel::NeuronNetModel & model;

	public:
		static std::string ClassName;

		TMSynapseFunct(AstroModel::TMSynapse & _syn);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double &, const double *v, double *f) const;
	};

	/******************************************************************/
	/* Tsodyks Markram synapse Optimized Function                     */
	/******************************************************************/
	class TMSynapseFunctOptim : public TMSynapseFunct
	{
	public:
		static std::string ClassName;

		TMSynapseFunctOptim(AstroModel::TMSynapseOptim & _syn);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double &, const double *v, double *f) const;
	};

	/******************************************************************/
	/* Spike Frequency Adaptation Leaky Integrate and Fire Neuron Func*/
	/******************************************************************/
	class SFALIFNeuronFunct : public NeuronFunct
	{
	protected:
		AstroModel::SFALIFNeuron & neuron;
		const AstroModel::NeuronNetModel & model;

	public:
		static std::string ClassName;

		SFALIFNeuronFunct(AstroModel::SFALIFNeuron & _neur);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double &, const double *v, double *f) const;
	};

	/******************************************************************/
	/* Dummy neuron firing spike trains loaded from a file            */
	/******************************************************************/
	class DummyNeuronFunct : public NeuronFunct
	{
	public:
		static std::string ClassName;

		DummyNeuronFunct(AstroModel::DummyNeuron & _neur);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double &, const double *v, double *f) const;
	};

	/******************************************************************/
	/* ChI Network Function                                           */
	/******************************************************************/
	class ChINetworkFunct : public AbstrNetwChIFunct
	{
	protected:
		AstroModel::ChIModel & model;

	public:
		static std::string ClassName;

		ChINetworkFunct(AstroModel::ChIModel & _model);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double & t, const double *v, double *f) const;
	};

	/******************************************************************/
	/* KChI Network Function                                          */
	/******************************************************************/
	class KChINetworkFunct : public AbstrNetwChIFunct
	{
	protected:
		AstroModel::KChIModel & model;

	public:
		static std::string ClassName;

		KChINetworkFunct(AstroModel::KChIModel & _model);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double & t, const double *v, double *f) const;
	};

	/******************************************************************/
	/* SFALIF Network Func                                            */
	/******************************************************************/
	class NeuronNetworkFunc : public Function<double>
	{
	protected:
		AstroModel::NeuronNetModel & model;

	public:
		static std::string ClassName;

		NeuronNetworkFunc(AstroModel::NeuronNetModel & _model);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double & t, const double *v, double *f) const;
	};

	/******************************************************************/
	/* Mixed astrocyte / neuron networks                              */
	/******************************************************************/
	class AstroNeuronNetFunc : public Function<double>
	{
	protected:
		AstroModel::AstroNeuroNetModel & model;

	public:
		static std::string ClassName;

		AstroNeuronNetFunc(AstroModel::AstroNeuroNetModel & _model);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double & t, const double *v, double *f) const;
	};

	/******************************************************************/
	/* Fire Diffuse model cell function                               */
	/******************************************************************/
	class FireDiffuseCellFunct : public Function<double>
	{
	protected:
		AstroModel::FireDiffuseCell & cell;
		const AstroModel::FireDiffuseModel & model;

	public:
		static std::string ClassName;

		FireDiffuseCellFunct(AstroModel::FireDiffuseCell & _cell);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double &, const double *v, double *f) const;
	};

	/******************************************************************/
	/* Fire Diffuse Network Function                                  */
	/******************************************************************/
	class FireDiffuseNetFunct : public Function<double>
	{
	protected:
		AstroModel::FireDiffuseModel & model;

	public:
		static std::string ClassName;

		FireDiffuseNetFunct(AstroModel::FireDiffuseModel & _model);

		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		virtual void CompFunc(const double & t, const double *v, double *f) const;
	};
}

#endif

