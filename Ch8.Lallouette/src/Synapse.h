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

#ifndef SYNAPSE_H
#define SYNAPSE_H

#include "ODEProblems.h"
#include "ODEFunctions.h"
#include "Savable.h"
#include "utility.h"
#include "ParamHandler.h"
#include "CouplingFunction.h"

#define TMSYNAPSE_NBVALS_PER_SYN 3 // u, x and gamma
#define TMSYNAPSE_OPTIM_NBVALS_PER_SYN 1 // gammaOpt


namespace AstroModel
{
	// Forward declarations

/**********************************************************************/
/* Abstract Synapse base class                                        */
/**********************************************************************/
	class Synapse : public NetworkEdge
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class Neuron;
		friend class SynapseFunct;
		friend class ODE::NeuronNetworkFunc;
		friend class ODE::AstroNeuronNetFunc;

	public:
		// Default Constructor
		Synapse(const NeuronNetModel * _model = 0, double *_dv = 0, 
			bool _fv = false, unsigned int _nbdv = 0);
		// Constructor from stream
		Synapse(std::ifstream & stream);
		// Copy constructor
		Synapse(const Synapse & s);
		// Destructor
		~Synapse();

		virtual std::string GetClassName() const = 0;

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the neuron to default values
		virtual void Initialize() = 0;
		// Treats a presynaptic spike
		virtual void PresynSpike(double t) = 0;
		// Update value names in ODEProblem
		virtual void SetValNamesPostfix(std::string pf) = 0;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler() = 0;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual unsigned int GetNbDynVal() const = 0;
		virtual double GetDynVal(unsigned int ind)
		{ return dynVals[ind]; }
		virtual void SetDynVals(double *_dv, bool _f);
		virtual void SetModel(const NeuronNetModel *_m);
		virtual bool IsCorrectlyInitialized() const;
		virtual void SetPreSynNeur(Neuron *_n);
		virtual void SetPostSynNeur(Neuron *_n);
		virtual const NeuronNetModel * GetModel() const;

	protected:
		//===========================================================||
		// Links to other objects                                    ||
		//===========================================================||
		ODE::SynapseFunct * funct;    // Derivative computation function
		const NeuronNetModel * model; // Pointer to parent model

		// Pointer to the pre-synaptic neuron
		Neuron *preSynNeur;
		// Pointer to the post-synaptic neuron
		Neuron *postSynNeur;

		//===========================================================||
		// Dynamic values                                            ||
		//===========================================================||
		double *dynVals; // Dynamic values
		bool freeDynVals;
		unsigned int nbDynVals;
	};

/**********************************************************************/
/* Glutamatergic Synapse class                                        */
/**********************************************************************/
	class GlutamatergicSynapse : public Synapse
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::AstroNeuronNetFunc;

	public:
		// Default Constructor
		GlutamatergicSynapse(const NeuronNetModel * _model = 0, double *_dv = 0, 
			bool _fv = false, unsigned int _nbdv = 0);
		// Constructor from stream
		GlutamatergicSynapse(std::ifstream & stream);

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual double GetGluVal() const = 0;
	};

/**********************************************************************/
/* Tsodyks Markram synapse model                                      */
/**********************************************************************/
	class TMSynapse : public GlutamatergicSynapse
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::TMSynapseFunct;
		friend class ODE::TMSynapseFunctOptim;
		// TEMP
		friend class SFALIFNeuron;
		friend class ODE::AstroNeuronNetFunc;

	public:
		static std::string ClassName;
		//===========================================================||
		// Enumeration binding indices to values names               ||
		//===========================================================||
		enum DynValNames {
			x = 0, // fraction of available ressources
			u,     // fraction of released available ressources
			gamma  // Extra-synaptic Glutamate concentration
		};

		//===========================================================||
		// Static constant equilibrium values for the 3 variables    ||
		//===========================================================||
		static double Defaultx;
		static double Defaultu;
		static double Defaultgamma;

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		TMSynapse(const NeuronNetModel * _model = 0, double *_dv = 0, 
			bool _fv = false);
		// Constructor from stream
		TMSynapse(std::ifstream & stream);
		// Copy constructor
		TMSynapse(const TMSynapse & s);
		// Destructor
		virtual ~TMSynapse();

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the synapse to default values
		virtual void Initialize();
		// Set the synapse to equilibrium (rest)
		virtual void SetToEquilibrium();
		// Treats a presynaptic spike
		virtual void PresynSpike(double t);
		// Update value names in ODEProblem
		virtual void SetValNamesPostfix(std::string pf);
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }
		// Update edge values in case parameters were changed
		virtual void UpdateEdge();

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the synapse from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the synapse to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual double GetGluVal() const;
		virtual unsigned int GetNbDynVal() const;
		// Change dynamic values to given pointer
		//virtual void SetDynVals(double *_dv, bool _f);

	protected:
		static double DefaultOd;
		static double DefaultOf;
		static double DefaultU0;
		static double DefaultrhoA;
		static double Defaultnv;
		static double DefaultGv;
		static double DefaultOc;
		static double DefaultSpillOvFract;

		//===========================================================||
		// Cell biochemical parameters                               ||
		//===========================================================||
		bool defaultBiophysParams; // Does the cell uses default 
		                           //   biophysical parameters ?

		double Od;   // Recovery rate of released synaptic 
		             //    glutamate vesicles
		double Of;   // Rate of synaptic facilitation
		double U0;   // Basal synaptic release probability
		double rhoA; // Ratio between the average volume of synaptic
		             //    vesicles and the ESS volume
		double nv;   // Number of releasable synaptic vesicles
		double Gv;   // Glutamate concentration within synaptic vesicles
		double Oc;   // Glutamate clearance rate in the ESS
		double spillOvFract; // fraction of glutamate that gets spilled over

	};

/**********************************************************************/
/* Optimized Tsodyks Markram synapse model                            */
/**********************************************************************/
	class TMSynapseOptim : public TMSynapse
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::TMSynapseFunctOptim;
		friend class SFALIFNeuron;
		friend class ODE::AstroNeuronNetFunc;

	public:
		static std::string ClassName;
		//===========================================================||
		// Enumeration binding indices to values names               ||
		//===========================================================||
		enum DynValNamesOptim {
			gammaOpt = 0 // Extra-synaptic Glutamate concentration
		};

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		TMSynapseOptim(const NeuronNetModel * _model = 0, double *_dv = 0, 
			bool _fv = false);
		// Constructor from stream
		TMSynapseOptim(std::ifstream & stream);
		// Copy constructor
		TMSynapseOptim(const TMSynapseOptim & s);
		// Destructor
		virtual ~TMSynapseOptim();

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the synapse to default values
		virtual void Initialize();
		// Set the synapse to equilibrium (rest)
		virtual void SetToEquilibrium();
		// Treats a presynaptic spike
		virtual void PresynSpike(double t);
		// Update value names in ODEProblem
		virtual void SetValNamesPostfix(std::string pf);
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual double GetGluVal() const;
		virtual unsigned int GetNbDynVal() const;

	protected:
		double lastSpikeTime;
		double lastU;
		double lastX;
	};
}

#endif
