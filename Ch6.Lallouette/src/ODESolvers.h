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

#ifndef ODESOLVERS_H
#define ODESOLVERS_H

#include <iostream>
#include <vector>
#include <algorithm>

#include "ODEFunctions.h"
#include "ODEProblems.h"
#include "ResultSaver.h"

namespace ODE
{
/**********************************************************************/
/* Abstract class                                                     */
/**********************************************************************/
	template <typename Val = double, typename StepT = double> class ODESolver : 
		public SaveAndLoadFromStream
	{
	public:
		ODESolver() : stepSize(0), currTime(0), needReset(false) {}
		ODESolver(std::ifstream & stream) : stepSize(0), currTime(0), needReset(false) 
		{
			LoadFromStream(stream);
			TRACE(stepSize)
		}
		ODESolver(StepT _step) : stepSize(_step), currTime(0), needReset(false) {}

		virtual void Solve(ODEProblem<Val, StepT> & prob, StepT start, StepT end)
		{
TRACE(stepSize)
			for (currTime = start ; currTime <= end ; currTime += stepSize)
			{
				DoStep(prob);
				/*******/
				if (needReset)
				{
					prob.SetToInitVals();
					needReset = false;
				}
				/*******/
				prob.UpdateVals(currTime);
			}
		}

		virtual void SetStepSize(StepT _sz) 
		{
			assert(_sz > 0);
			stepSize = _sz;
		}

		virtual void ResetVals()
		{
			needReset = true;
		}

		// Return class name
		virtual std::string GetClassName() const = 0;

		//===========================================================||
		// Standard Save and Load methods                            ||
		//===========================================================||
		// Loads the neuron from a stream
		virtual bool LoadFromStream(std::ifstream & stream)
		{
			stream >> stepSize;
			return stream.good();
		}
		// Saves the neuron to a stream
		virtual bool SaveToStream(std::ofstream & stream) const
		{
			stream << stepSize << std::endl;
			return stream.good();
		}

	protected:
		StepT stepSize;
		StepT currTime;
		std::vector<Val> tempVals;
		bool needReset;

		virtual void DoStep(ODEProblem<Val, StepT> & prob) = 0;
	};

/**********************************************************************/
/* Euler Solver                                                       */
/**********************************************************************/
	template <typename Val = double, typename StepT = double> class EulerSolver : public ODESolver<Val, StepT>
	{
	public:
		static std::string ClassName;

		EulerSolver() : ODESolver<Val, StepT>::ODESolver(), 
			deriv(0), nbVals(0) {}
		EulerSolver(std::ifstream & stream) : ODESolver<Val, StepT>::ODESolver(stream),
			deriv(0), nbVals(0) {}
		EulerSolver(StepT _step) : ODESolver<Val, StepT>::ODESolver(_step),
			deriv(0), nbVals(0) {}

		virtual void Solve(ODEProblem<Val, StepT> & prob, StepT start, StepT end)
		{
			nbVals = prob.GetNbVals();
			deriv = new Val[nbVals];

			ODESolver<Val, StepT>::Solve(prob, start, end);

			delete[] deriv;
			deriv = 0;
			nbVals = 0;
		}
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }

	protected:
		Val *deriv;
		unsigned int nbVals;

		virtual void DoStep(ODEProblem<Val, StepT> & prob)
		{
			prob.function->CompFunc(this->currTime, prob.vals, deriv);
			for (unsigned int i = 0 ; i < nbVals ; ++i)
				prob.vals[i] = prob.vals[i] + this->stepSize * deriv[i];
		}
	};


/**********************************************************************/
/* Runge Kutta Solver                                                 */
/**********************************************************************/
	template <typename Val = double, typename StepT = double> class RungeKuttaSolver : public ODESolver<Val, StepT>
	{
	public:
		static std::string ClassName;

		RungeKuttaSolver() : ODESolver<Val, StepT>::ODESolver(), 
			k1(0), k2(0), k3(0), k4(0), tempVal(0), nbVals(0) {}
		RungeKuttaSolver(std::ifstream & stream) : ODESolver<Val, StepT>::ODESolver(stream),
			k1(0), k2(0), k3(0), k4(0), tempVal(0), nbVals(0) {}
		RungeKuttaSolver(StepT _step) : ODESolver<Val, StepT>::ODESolver(_step),
			k1(0), k2(0), k3(0), k4(0), tempVal(0), nbVals(0) {}

		virtual void Solve(ODEProblem<Val, StepT> & prob, StepT start, StepT end)
		{
			nbVals = prob.GetNbVals();
			k1 = new Val[nbVals];
			k2 = new Val[nbVals];
			k3 = new Val[nbVals];
			k4 = new Val[nbVals];
			tempVal = new Val[nbVals];

			ODESolver<Val, StepT>::Solve(prob, start, end);

			delete[] k1;
			delete[] k2;
			delete[] k3;
			delete[] k4;
			delete[] tempVal;
			k1 = 0;
			k2 = 0;
			k3 = 0;
			k4 = 0;
			tempVal = 0;
			nbVals = 0;
		}
		// Return class name
		virtual std::string GetClassName() const { return ClassName; }

	protected:
		Val *k1;
		Val *k2;
		Val *k3;
		Val *k4;
		Val *tempVal;
		unsigned int nbVals;

		virtual void DoStep(ODEProblem<Val, StepT> & prob)
		{
			for (unsigned int i = 0 ; i < nbVals ; ++i)
				tempVal[i] = prob.vals[i];

			prob.function->CompFunc(this->currTime, prob.vals, k1);
			for (unsigned int i = 0 ; i < nbVals ; ++i)
			{
				k1[i] = this->stepSize * k1[i];
				prob.vals[i] = tempVal[i] + k1[i] / 2.0;
			}

			prob.function->CompFunc(this->currTime + this->stepSize / 2.0, prob.vals, k2);
			for (unsigned int i = 0 ; i < nbVals ; ++i)
			{
				k2[i] = this->stepSize * k2[i];
				prob.vals[i] = tempVal[i] + k2[i] / 2.0;
			}

			prob.function->CompFunc(this->currTime + this->stepSize / 2.0, prob.vals, k3);
			for (unsigned int i = 0 ; i < nbVals ; ++i)
			{
				k3[i] = this->stepSize * k3[i];
				prob.vals[i] = tempVal[i] + k3[i];
			}

			prob.function->CompFunc(this->currTime + this->stepSize, prob.vals, k4);
			for (unsigned int i = 0 ; i < nbVals ; ++i)
				k4[i] = this->stepSize * k4[i];

			for (unsigned int i = 0 ; i < nbVals ; ++i)
				prob.vals[i] = (tempVal[i] + k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0);
		}
	};
}

#endif
