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

#ifndef COUPLINGJUNCTION_H
#define COUPLINGJUNCTION_H

#include "Savable.h"
#include "ParamHandler.h"
#include <string>

namespace AstroModel
{
	class AbstractNetwork;
	template <typename LinkType> class Network;

/**********************************************************************/
/* Abstract network edge class                                        */
/**********************************************************************/
	class NetworkEdge : public SaveAndLoadFromStream
	{
	public:
		virtual std::string GetClassName() const = 0;
		// Update edge values in case parameters were changed
		virtual void UpdateEdge() = 0;
	};

/**********************************************************************/
/* Abstract class                                                     */
/**********************************************************************/
	class CouplingFunction : public NetworkEdge
	{
	protected:
		// Allowed values :
		//     0 Truncated gaussian (negative are put to 0)
		//     1 Symetrically truncated gaussian
		//     2 Same as 1 but values under 0 or over 2F are 
		//       respectively sticked to 0 and 2*F.
		enum DefCouplMeth {
			TruncGauss = 0,
			SymTruncGauss = 1,
			StickSymTruncGauss = 2,
			GammaLaw = 3
		};

	public:
		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		CouplingFunction(double _F = DefaultCouplingStrength,
			double _sigF = DefaultCouplStrengthStdDev, 
			int _meth = DefaultCouplDistrMethod);

		//===========================================================||
		// Standard Coupling Function methods                        ||
		//===========================================================||
		// Coupling Function operator (returns the flux)
		virtual double operator()(double DIP3) const = 0;

		// Build parameter handling object, handles static default params
		virtual ParamHandler BuildModelParamHandler();

		// Normalizes link strengths according to node degrees and 
		// current strength values
		static void normalizeLinkStrength(AbstractNetwork & network);

		// Update edge values in case parameters were changed
		virtual void UpdateEdge();

		//===========================================================||
		// Accessors                                                 ||
		//===========================================================||
		virtual void ChangeStrength(double _F) { F = _F; }
		virtual double GetStrength() const { return F; }
		virtual double & GetRefOnStrength() { return F; }
		virtual void SetConstStrength(bool cs) {constStrength = cs;}

	protected:
		//===========================================================||
		// Parameters                                                ||
		//===========================================================||
		double F; // Coupling Strength
		bool constStrength; // If true, F cannot be modified during updateEdge

		static double DefaultCouplingStrength;
		static double DefaultCouplStrengthStdDev;
		static int DefaultCouplDistrMethod;

		// Set Strength according to specified distribution
		virtual void redrawStrength(double _F, double _sigF, int _meth);
	};

/**********************************************************************/
/* Sigmoid Coupling Function                                          */
/**********************************************************************/
	class SigmoidCoupling : public CouplingFunction
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		SigmoidCoupling(
			double _F = CouplingFunction::DefaultCouplingStrength, 
			double _sigF = CouplingFunction::DefaultCouplStrengthStdDev, 
			double _IP3Scale = DefaultLinkScale, 
			double _IP3ScaleSD = DefaultLnkScaleStdDev,
			double _IP3Thresh = DefaultLinkThreshold,
			double _IP3ThreshSD = DefaultLnkThreshStdDev);
		// Constructor from stream
		SigmoidCoupling(std::ifstream & stream);

		//===========================================================||
		// Standard Coupling Function methods                        ||
		//===========================================================||
		// Coupling Function operator (returns the flux)
		virtual double operator()(double DIP3) const;
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;

		// Build parameter handling object, handles static default params
		virtual ParamHandler BuildModelParamHandler();

		// Update edge values in case parameters were changed
		virtual void UpdateEdge();

	protected:
		//===========================================================||
		// Parameters                                                ||
		//===========================================================||
		double IP3Scale;  // Slope factor
		double IP3Thresh; // Half-maximal diffusion IP3 Threshold

		static double DefaultLinkThreshold;
		static double DefaultLnkThreshStdDev;
		static double DefaultLinkScale;
		static double DefaultLnkScaleStdDev;
	};

/**********************************************************************/
/* Linear Coupling Function                                           */
/**********************************************************************/
	class LinearCoupling : public CouplingFunction
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		// Default Constructor
		LinearCoupling(double _F = 
			CouplingFunction::DefaultCouplingStrength);
		// Constructor from stream
		LinearCoupling(std::ifstream & stream);

		// Update edge values in case parameters were changed
		virtual void UpdateEdge();

		//===========================================================||
		// Standard Coupling Function methods                        ||
		//===========================================================||
		// Coupling Function operator (returns the flux)
		virtual double operator()(double DIP3) const;
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream);
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const;
	};

/**********************************************************************/
/* Boolean link                                                       */
/**********************************************************************/
	class BooleanLink : public NetworkEdge
	{
	public:
		static std::string ClassName;
		// Returns the class name
		virtual std::string GetClassName() const { return ClassName; }

		//===========================================================||
		// Constructors / Destructor                                 ||
		//===========================================================||
		BooleanLink() : activated(true) {}
		BooleanLink(std::ifstream & stream) : activated(true) 
			{ LoadFromStream(stream); } 

		void SetState(bool s) { activated = s; }
		static ParamHandler BuildModelParamHandler() 
			{ return ParamHandler(); }
		// Loads the coupling function from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			stream >> activated;
			return stream.good() and not stream.eof();
		}
		// Saves the coupling function to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			stream << activated << std::endl;
			return stream.good();
		}
		// Update edge values in case parameters were changed
		virtual void UpdateEdge() {}

	protected:
		bool activated;
	};

}

#endif
