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

#ifndef CHICELL_H
#define CHICELL_H

#include "ODEProblems.h"
#include "ODEFunctions.h"
#include "Savable.h"
#include "utility.h"
#include "ParamHandler.h"

#define CHIMODEL_NBVALS_PER_CELL 3

namespace ODE
{
	// Forward declarations
	class ChICellFunct;
	class ChINetworkFunct;
}

namespace AstroModel
{
	// Forward declarations
	class ChIModel;
	template <typename NE, typename NT> class ODENetworkDynamicsModel;
	class CouplingFunction;

/**********************************************************************/
/* ChI Cell                                                           */
/**********************************************************************/
	class ChICell : public SaveAndLoadFromStream
	{
		//===========================================================||
		// Friend declarations                                       ||
		//===========================================================||
		friend class ODE::ChICellFunct;
		friend class ODE::ChINetworkFunct;
		friend class ODE::KChICellFunct;
		friend class ODE::KChINetworkFunct;
		friend class ChIModel;
		friend class KChIModel;
		friend class ExtraCellGluStim;
		friend class PoissonianStimStrat;
		friend class ODE::AstroNeuronNetFunc;

	public:
		static std::string ClassName;

		static unsigned int NbValsPerCell;
		//===========================================================||
		// Enumeration binding indices to values names               ||
		//===========================================================||
		enum DynValNames {
			Ca = 0, // Cell-averaged Ca2+ concentration
			h,      // Fraction of open IP3R channels on the ER membrane
			IP3     // Cell-averaged concentration of IP3 second messenger
		};

		//===========================================================||
		// Static constant equilibrium values for the 3 variables    ||
		//===========================================================||
		static double DefaultCa;
		static double Defaulth;
		static double DefaultIP3;

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		ChICell(const ODENetworkDynamicsModel<CouplingFunction, ChICell> * _model, 
			double *_dv = 0, bool _fv = false);
		// Copy constructor
		ChICell(const ChICell & c);
		// Destructor
		virtual ~ChICell();

		//===========================================================||
		// Model type methods                                        ||
		//===========================================================||
		// Initializes the cell to default values
		virtual void Initialize();
		// Set the cell to equilibrium
		virtual void SetToEquilibrium();
		// Update value names in ODEProblem
		virtual void SetValNamesPostfix(std::string pf);
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the cell from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the cell to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		//===========================================================||
		// Special model parameters handling method                  ||
		//===========================================================||
		virtual ParamHandler BuildModelParamHandler();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		// Change dynamic values to given pointer
		virtual void SetDynVals(double *_dv, bool _f);
		// Return the number of dyn vals per cell
		virtual unsigned int GetNbDynVals();
		// Return the value of a dyn val
		inline double GetDynVal(int name) const
			{ return dynVals[name]; }
		// Return the number of dyn vals
		inline unsigned int GetNbDynVals() const
			{ return nbDynVals;}

	protected:
		static double InitCaVarRatio;
		static double InithVarRatio;
		static double InitIP3VarRatio;
		static double Defaultc1;
		static double DefaultrC;
		static double DefaultrL;
		static double DefaultvER;
		static double Defaulta2;
		static double DefaultKer;  
		static double Defaultvd;   
		static double DefaultKplcd;
		static double Defaultkd;   
		static double Defaultk3;   
		static double Defaultvbeta;
		static double DefaultkR;   
		static double DefaultkP;   
		static double Defaultkpi;  

		//===========================================================||
		// Links to other objects                                    ||
		//===========================================================||
		ODE::AbstrChIFunct * funct; // Derivative computation function
		const ChIModel * model;     // Pointer to parent model

		//===========================================================||
		// Cell biochemical parameters                               ||
		//===========================================================||
		bool defaultBiophysParams; // Does the cell uses default biophysical parameters ?
		double c1;    // Ratio between ER and cytosol volumes
		double rC;    // Maximal CICR rate
		double rL;    // Maximal rate of Ca2+ leak from the ER
		double vER;   // Maximal SERCA uptake rate
		double a2;    // IP3R binding rate constant for Ca2+ inhibition
		double Ker;   // SERCA Ca2+ affinity
		double vd;    // Maximal rate of IP3 synthesis by PLCd
		double Kplcd; // Ca2+ affinity of PLCd
		double kd;    // Inhibition constant of PLCd activity
		double k3;    // Half-saturation constant for Ca2+-dependent IP3-3K activation
		// G-ChI
		double vbeta; // Maximal rate of IP3 production by PLCbeta
		double kR;    // Glutamate affinity of the receptor
		double kP;    // Ca2+/PKC-dependent inhibition factor
		double kpi;   // Ca2+ affinity of PKC

		//===========================================================||
		// Dynamic values                                            ||
		//===========================================================||
		double *dynVals; // Dynamic values
		bool freeDynVals;
		unsigned int nbDynVals;
		double totFlux;  // Total fluxes of IP3 from the cell
		bool caSpontLeak; // Spontaneous leak of Ca2+ from the ER
		double gluIP3Prod;
	};
}

#endif
