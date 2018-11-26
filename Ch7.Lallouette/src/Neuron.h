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

#ifndef NEURON_H
#define NEURON_H

#include "ODEProblems.h"
#include "ODEFunctions.h"
#include "Savable.h"
#include "utility.h"
#include "ParamHandler.h"

#define SFALIFNEURON_NBVALS_PER_NEURON 2 // V and w
#define DUMMYNEURON_NBVALS_PER_NEURON 1  // neuron state

namespace AstroModel
{
/**********************************************************************/
/* Abstract Neuron base class                                         */
/**********************************************************************/
	class Neuron : public SaveAndLoadFromStream
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class Synapse;
		friend class ODE::NeuronNetworkFunc;

	public:
		// Default constructor
		Neuron(const NeuronNetModel * _model = 0);
		// Constructor from stream
		Neuron(std::ifstream & stream);
		// Copy constructor
		Neuron(const Neuron & n);
		// Destructor
		~Neuron();

		virtual std::string GetClassName() const = 0;

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the neuron to default values
		virtual void Initialize() = 0;
		virtual bool CheckForSpiking(double t) = 0;
		virtual void SetValNamesPostfix(std::string pf) = 0;

		virtual void ForceSpiking(double t) = 0;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler() = 0;

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual unsigned int GetNbDynVal() const = 0;
		virtual void SetDynVals(double *_dv, bool _f) = 0;
		virtual void SetModel(const NeuronNetModel *_m);
		virtual bool IsCorrectlyInitialized() const;
		virtual void AddDendrSyn(Synapse *_s);
		virtual void AddAxonSyn(Synapse *_s);
		virtual void ClearSynapses();
		virtual const NeuronNetModel * GetModel() const;

	protected:
		//===========================================================||
		// Links to other objects                                    ||
		//===========================================================||
		ODE::NeuronFunct * funct;     // Derivative computation function
		const NeuronNetModel * model; // Pointer to parent model

		// Pointers to synapses impiging on the dendritic tree
		std::vector<Synapse *> dendrSyn;
		// Pointers to synapses formed by the axon
		std::vector<Synapse *> axonSyn;
	};

/**********************************************************************/
/* Dummy Neuron (fires spikes loaded from a file)                     */
/**********************************************************************/
	class DummyNeuron : public Neuron
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::SFALIFNeuronFunct;
		friend class NeuronNetModel;
		friend class TMSynapse;

	public:
		static std::string ClassName;
		//===========================================================||
		// Enumeration binding indices to values names               ||
		//===========================================================||
		enum DynValNames {
			V = 0 // Membrane state (-1.0 when pol and 0 when depol)
		};

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		DummyNeuron(const NeuronNetModel * _model = 0, double *_dv = 0, 
			bool _fv = false);
		// Constructor from stream
		DummyNeuron(std::ifstream & stream);
		// Copy constructor
		DummyNeuron(const DummyNeuron & s);
		// Destructor
		virtual ~DummyNeuron();

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the neuron to default values
		virtual void Initialize();
		// Check for potential spiking
		virtual bool CheckForSpiking(double t);
		// Update value names in ODEProblem
		virtual void SetValNamesPostfix(std::string pf);
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }

		// Forces the neuron into emitting a spike, regardless of its state
		virtual void ForceSpiking(double t);

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the neuron from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the neuron to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual unsigned int GetNbDynVal() const;
		// Change dynamic values to given pointer
		virtual void SetDynVals(double *_dv, bool _f);

	protected:
		// Times at which the neuron is supposed to fire
		std::vector<double> spikingTimes; 
		unsigned int nextSpike;

		//===========================================================||
		// Dynamic values                                            ||
		//===========================================================||
		double *dynVals; // Dynamic values
		bool freeDynVals;
		unsigned int nbDynVals;
	};
/**********************************************************************/
/* Spike Frequency Adaptation Leaky Integrate and Fire Neuron         */
/**********************************************************************/
	class SFALIFNeuron : public Neuron
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::SFALIFNeuronFunct;
		friend class NeuronNetModel;
		// temp
		friend class TMSynapse;

	public:
		static std::string ClassName;
		//===========================================================||
		// Enumeration binding indices to values names               ||
		//===========================================================||
		enum DynValNames {
			V = 0, // Membrane potential
			w      // Hyperpolarizing current
		};

		//===========================================================||
		// Static constant equilibrium values for the 3 variables    ||
		//===========================================================||
		static double DefaultV;
		static double Defaultw;

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		SFALIFNeuron(const NeuronNetModel * _model = 0, double *_dv = 0, 
			bool _fv = false);
		// Constructor from stream
		SFALIFNeuron(std::ifstream & stream);
		// Copy constructor
		SFALIFNeuron(const SFALIFNeuron & s);
		// Destructor
		virtual ~SFALIFNeuron();

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the neuron to default values
		virtual void Initialize();
		// Set the neuron to equilibrium (rest)
		virtual void SetToEquilibrium();
		// Check for potential spiking
		virtual bool CheckForSpiking(double t);
		// Update value names in ODEProblem
		virtual void SetValNamesPostfix(std::string pf);
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }

		// Forces the neuron into emitting a spike, regardless of its state
		virtual void ForceSpiking(double t);

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the neuron from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the neuron to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual unsigned int GetNbDynVal() const;
		// Change dynamic values to given pointer
		virtual void SetDynVals(double *_dv, bool _f);

	protected:
		static double DefaultC;   
		static double DefaultgL;  
		static double DefaultE0;  
		static double DefaultVt;  
		static double DefaultVr;  
		static double Defaulttauw;
		static double Defaultb;   

		//===========================================================||
		// Cell biochemical parameters                               ||
		//===========================================================||
		bool defaultBiophysParams; // Does the cell uses default 
		                           //   biophysical parameters ?

		double C;    // Capacitance of the membrane
		double gL;   // Conductance of the membrane
		double E0;   // Resting potential
		double Vt;   // Spike threshold potential
		double Vr;   // Reset potential
		double tauw; // Time constant for w decay
		double b;    // w increases

		double inputCurr;

		//===========================================================||
		// Dynamic values                                            ||
		//===========================================================||
		double *dynVals; // Dynamic values
		bool freeDynVals;
		unsigned int nbDynVals;
	};
}

#endif

