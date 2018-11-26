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

#include "StimulationStrat.h"
#include "ChIModel.h"
#include "AbstractFactory.h"
#include "MetricComputeStrat.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include <gsl/gsl_sf_exp.h>

using namespace AstroModel;
using namespace std;

static const double NB_POISS_START_DELAY = 5.0;

//********************************************************************//
//**************** S T I M U L A T I O N   S T R A T *****************//
//********************************************************************//

string DefaultStimStrat::ClassName("DefaultStimulationStrategy");
string PoissonianStimStrat::ClassName("PoissonianStimulationStrategy");
string RandomMonoStim::ClassName("RandomMonoStimulationStrategy");
string ThresholdDeterminationStimStrat::ClassName(
	"ThresholdDeterminationStimulationStrategy");
string MeanPathStimStrat::ClassName("MeanPathLengthStimulationStrategy");
string CorrelDetermStimStrat::ClassName("CorrelationDeterminationStimulationStrategy");
string ExtraCellGluStim::ClassName("ExtraCellularGlutamateStimulation");

//********************************************************************//
//************ S T I M U L A B L E   B A S E   C L A S S *************//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
Stimulable::Stimulable(ParamHandler & h)
{
	std::vector<std::string> stimClassNames = h.getParam<std::vector<std::string> >("-StimType");
	for (unsigned int i = 0 ; i < stimClassNames.size() ; ++i)
	{
		if (AbstractFactory<StimulationStrat>::Factories[stimClassNames[i]])
			stimStrats.push_back(AbstractFactory<StimulationStrat>::Factories[stimClassNames[i]]->Create());
		else
			std::cerr << "Couldn't create the following Stimulation strategy : " << stimClassNames[i] << std::endl;
	}
	stimStratActMask = std::vector<int>(stimStrats.size(), true);
}

//**********************************************************************
// Loading Constructor
//**********************************************************************
Stimulable::Stimulable(std::ifstream & stream, ParamHandler & )
{
	if (not LoadFromStream(stream))
		cerr << "Failed to load the model !" << endl;
}

//**********************************************************************
// Destructor
//**********************************************************************
Stimulable::~Stimulable()
{
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
		delete stimStrats[i];
}

//**********************************************************************
// Initializes the model
//**********************************************************************
void Stimulable::Initialize()
{
	// Initialize stimulation strategies
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
		stimStrats[i]->Initialize();
	// Initialize stimulation strategies mask
	if (stimStratActMask.size() != stimStrats.size())
		stimStratActMask = std::vector<int>(stimStrats.size(), true);
}

//**********************************************************************
// Stimulates cells using the stimulation strategy
//**********************************************************************
void Stimulable::Stimulate(double t)
{
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
	{
		if (stimStratActMask[i] != 0)
			stimStrats[i]->Stimulate(*this, t);
	}
}

//**********************************************************************
// Is the given cell currently stimulated ?
//**********************************************************************
bool Stimulable::IsStimulated(unsigned int ind) const 
{
	bool isStim = false;
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
		isStim |= stimStrats[i]->IsStimulated(ind);
	return isStim; 
}

//********************************************************************//
// Returns true if the given stimulation strat is added to the model 
// and activated
//********************************************************************//
bool Stimulable::IsStimStratActivated(std::string name) const
{
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
		if ((stimStrats[i]->GetClassName() == name) and (stimStratActMask[i] != 0))
			return true;
	return false;
}

//********************************************************************//
// Save Stim strat metrics                                            //
//********************************************************************//
bool Stimulable::SaveStimStratMetrics(ResultSaver saver) const
{
	bool ok = true;
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
		ok &= stimStrats[i]->SaveMetrics(saver);
TRACE("Stim strats metrics saved.")
	return ok;
}

//**********************************************************************
// Loads the model from a stream
//**********************************************************************
bool Stimulable::LoadFromStream(std::ifstream & stream)
{
	bool ok = true;
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
		delete stimStrats[i];
	stimStrats.clear();

	unsigned int nbStimStrats = 0;
	stream >> nbStimStrats;
	string tempName;
	for (unsigned int i = 0 ; i < nbStimStrats ; ++i)
	{
		stream >> tempName;
		assert(AbstractFactory<StimulationStrat>::Factories[tempName]);
		stimStrats.push_back(
			AbstractFactory<StimulationStrat>::Factories[tempName]->
			CreateFromStream(stream)
		);
	}
	stimStratActMask = std::vector<int>(stimStrats.size(), true);
	return ok and (stream.good() or stream.eof());
}

//**********************************************************************
// Saves the model to a stream
//**********************************************************************
bool Stimulable::SaveToStream(std::ofstream & stream) const
{
	bool ok = true;
	stream << stimStrats.size() << endl;
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
	{
		stream << stimStrats[i]->GetClassName() << endl;
		ok &= stimStrats[i]->SaveToStream(stream);
	}
	return ok;
}

//**********************************************************************
// Returns a ParamHandler object with references to internal parameters
//**********************************************************************
ParamHandler Stimulable::BuildModelParamHandler()
{
	ParamHandler params;

	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
		params <= (stimStrats[i]->GetClassName() + "_Activ"), stimStratActMask[i];

	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
		params += stimStrats[i]->BuildModelParamHandler();

	return params;
}

//**********************************************************************
// Dynamically dispatch a metric according to its type
//**********************************************************************
bool Stimulable::AddMetric(Metric *_m, bool _f)
{
	bool workedForOneStimStrat = false;
	bool worked = false;
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
	{
		// Adding copies of the metric to stimulation strategies, 
		// only the first one can delete the metric
		worked = stimStrats[i]->AddMetric(_m, true);
		if (worked)
			_m = _m->BuildCopy();
		workedForOneStimStrat |= worked;
	}
	if (workedForOneStimStrat and _f)
		delete _m;
	return workedForOneStimStrat;
}

//**********************************************************************
// Returns all metrics and submetrics
//**********************************************************************
std::vector<Metric *> Stimulable::GetAllMetrics() const
{
	vector<Metric *> mTot, mTemp;
	for (unsigned int i = 0 ; i < stimStrats.size() ; ++i)
	{
		mTemp = stimStrats[i]->GetAllMetrics();
		for (unsigned int i = 0 ; i < mTemp.size() ; ++i)
			mTot.push_back(mTemp[i]);
	}
	return mTot;
}
//********************************************************************//
//************ S T I M U L A B L E   N E T S   C L A S S *************//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
StimulableCellNetwork::StimulableCellNetwork(ParamHandler & h) : Stimulable::Stimulable(h)
{
}

//********************************************************************//
//**************** S T I M U L A T I O N   S T R A T *****************//
//********************************************************************//

//**********************************************************************
// Destructor
//**********************************************************************
StimulationStrat::~StimulationStrat()
{
}

//**********************************************************************
// Strategy main method, calls the specialized protected stimulate method
//**********************************************************************
//void StimulationStrat::Stimulate(ChIModel & model, double t)
void StimulationStrat::Stimulate(Stimulable & model, double t)
{
	tCurr = t;
	stimulate(model);

	metrics.ComputeMetricsDefault(*this);
}

//**********************************************************************
// Save all metrics
//**********************************************************************
bool StimulationStrat::SaveMetrics(ResultSaver saver) const
{
	bool ok = true;
	ok &= metrics.SaveMetricsDefault(saver);
	return ok;
}

//**********************************************************************
// Initialize the stimulation strategy as if it were just created
//**********************************************************************
void StimulationStrat::Initialize()
{
	stimulated.clear();
	metrics.InitializeMetricsDefault();
}

//**********************************************************************
// Loads the cell from a stream
//**********************************************************************
bool StimulationStrat::LoadFromStream(std::ifstream & stream)
{
	return metrics.LoadFromStream(stream);
}

//**********************************************************************
// Saves the cell to a stream
//**********************************************************************
bool StimulationStrat::SaveToStream(std::ofstream & stream) const
{
	return metrics.SaveToStream(stream);
}

//**********************************************************************
// Builds parameter handling object
//**********************************************************************
ParamHandler StimulationStrat::BuildModelParamHandler()
{
	return metrics.BuildModelParamHandler();
}

//**********************************************************************
// Returns all metrics and submetrics
//**********************************************************************
vector<Metric *> StimulationStrat::GetAllMetrics() const
{
	vector<Metric *> mTot = metrics.GetMetricsRaw();
	return mTot;
}

//********************************************************************//
//*************** N E T W O R K   S T I M   S T R A T ****************//
//********************************************************************//

//**********************************************************************
// Stimulates the given model 
//**********************************************************************
void NetworkStimulationStrat::stimulate(Stimulable & model)
{
	StimulableCellNetwork *scn = dynamic_cast<StimulableCellNetwork*>(&model);
	if (scn)
	{
		if (scn->GetNbCells() != stimulated.size())
			stimulated = vector<bool>(scn->GetNbCells(), false);
		stimulateNet(*scn);
	}
}

//********************************************************************//
//*************** D E F A U L T   S T I M   S T R A T ****************//
//********************************************************************//

//**********************************************************************
// Default Constructor
//**********************************************************************
DefaultStimStrat::DefaultStimStrat(ParamHandler & h) : unusedRemoved(false)
{
	std::vector<std::string> names = h.getParam<vector<std::string> >("-defStim", 0);
	std::vector<int> indices       = h.getParam<vector<int> >("-defStim", 1);
	std::vector<double> starts     = h.getParam<vector<double> >("-defStim", 2);
	std::vector<double> ends       = h.getParam<vector<double> >("-defStim", 3);
	std::vector<int> expand        = h.getParam<vector<int> >("-defStim", 4);
	IP3Bias = h.getParam<double>("-IPBias", 0);
	std::string couplFunctName = h.getParam<std::string>("-defStimCoupl", 0);
	double couplingStrength = h.getParam<double>("-defStimCoupl", 1);
	
	assert(indices.size() == names.size());
	assert(indices.size() == starts.size());
	assert(indices.size() == ends.size());
	assert(indices.size() == expand.size());

	for (unsigned int i = 0 ; i < indices.size() ; ++i)
		AddCellStim(names[i], indices[i], starts[i], ends[i], std::max(0, expand[i]));

	assert(AbstractFactory<NetworkEdge>::Factories[couplFunctName]);
	NetworkEdge *netEdge = AbstractFactory<NetworkEdge>::Factories[couplFunctName]->Create();
	funct = dynamic_cast<CouplingFunction *>(netEdge);
	assert(funct);
	funct->ChangeStrength(couplingStrength);
	freeFunct = true;
}

//**********************************************************************
// Special Constructor
//**********************************************************************
DefaultStimStrat::DefaultStimStrat(unsigned int ind, double _s, double _e,
	double bias, CouplingFunction *_f) : funct(_f), freeFunct(not _f), 
	IP3Bias(bias), unusedRemoved(false)
{
	cellsToStim.push_back(StimulationInd(ClassName, ind, _s, _e));
	if (not funct)
		funct = new SigmoidCoupling();
}

//**********************************************************************
//**********************************************************************
DefaultStimStrat::DefaultStimStrat(double bias, CouplingFunction *_f) 
	: funct(_f), freeFunct(not _f), IP3Bias(bias), unusedRemoved(false)
{
	if (not funct)
		funct = new SigmoidCoupling();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
DefaultStimStrat::DefaultStimStrat(std::ifstream & stream) : funct(0), 
	freeFunct(false), unusedRemoved(false)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
DefaultStimStrat::~DefaultStimStrat()
{
	if (freeFunct and funct)
		delete funct;
}

//**********************************************************************
// Stimulates the cells according to planned stimulations
//**********************************************************************
//void DefaultStimStrat::stimulate(ChIModel & model)
void DefaultStimStrat::stimulateNet(StimulableCellNetwork & model)
{
	if (not unusedRemoved)
	{
		// Remove stimulation indices that were meant for another stim strat
		std::vector<StimulationInd> CTScopy = cellsToStim;
		cellsToStim.clear();
		for (unsigned int i = 0 ; i < CTScopy.size() ; ++i)
			if (CTScopy[i].stratName.find(GetClassName()) != std::string::npos)
				cellsToStim.push_back(CTScopy[i]);
	}

	set<unsigned int> expInds, tempInds;
	for (unsigned int i = 0 ; i < cellsToStim.size() ; ++i)
	{
		expInds.insert(cellsToStim[i].ind);
		for (; cellsToStim[i].expandToDo > 0 ; --cellsToStim[i].expandToDo)
		{
			for (set<unsigned int>::iterator it = expInds.begin() ; it != expInds.end() ; ++it)
			{
				const vector<unsigned int> & neighbs = model.GetNetwork().GetNeighbors(*it);
				tempInds.insert(neighbs.begin(), neighbs.end());
			}
			expInds.insert(tempInds.begin(), tempInds.end());
			tempInds.clear();
		}
		expInds.erase(cellsToStim[i].ind);
		for (set<unsigned int>::iterator it = expInds.begin() ; it != expInds.end() ; ++it)
			AddCellStim(GetClassName() + "_Expand", *it, cellsToStim[i].start, cellsToStim[i].end);
		expInds.clear();
	}

	vector<bool> tempTreated(stimulated.size(), false);
	for (unsigned int i = 0 ; i < cellsToStim.size() ; ++i)
	{
		if ((this->tCurr >= cellsToStim[i].start) and (this->tCurr <= cellsToStim[i].end))
		{
			assert(model.GetNbCells() > cellsToStim[i].ind);

			this->stimulated[cellsToStim[i].ind] = true;

			this->stimulateSpecific(&model, cellsToStim[i].ind);

			tempTreated[cellsToStim[i].ind] = true;
		}
		else if (not tempTreated[cellsToStim[i].ind])
			this->stimulated[cellsToStim[i].ind] = false;
	}
}

//**********************************************************************
// Actually stimulate with appropriate method
//**********************************************************************
void DefaultStimStrat::stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const
{
	ChIModel *chimod;
	if ((chimod = dynamic_cast<ChIModel*>(model)))
	{
		double diffIP3 = (chimod->GetDynVal(ind, ChICell::IP3) < IP3Bias) ? 
			(chimod->GetDynVal(ind, ChICell::IP3) - IP3Bias) : 0.0 ;
		chimod->ModifFluxes(ind, (*funct)(diffIP3));
	}
	else
		throw "Unkwnown model type";
}

//**********************************************************************
// Initialize the stimulation strategy as if it were just created
//**********************************************************************
void DefaultStimStrat::Initialize()
{
	StimulationStrat::Initialize();
	// Clear expands
	std::vector<StimulationInd> CTScopy = cellsToStim;
	cellsToStim.clear();
	for (unsigned int i = 0 ; i < CTScopy.size() ; ++i)
		if (CTScopy[i].stratName.find("_Expand") == std::string::npos)
		{
			cellsToStim.push_back(CTScopy[i]);
			cellsToStim.back().expandToDo = cellsToStim.back().expand;
		}
}

//**********************************************************************
// Adds a stimulation order on cell 'CellNbr' which will start at 
// 'tStart' and end at 'tEnd'
//**********************************************************************
void DefaultStimStrat::AddCellStim(std::string name, unsigned int cellNbr, 
	double tStart, double tEnd, unsigned int exp)
{
	cellsToStim.push_back(StimulationInd(name, cellNbr, tStart, tEnd, exp));
}

//**********************************************************************
// Returns a param handler of the stimulation strat parameters
//**********************************************************************
ParamHandler DefaultStimStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params += StimulationStrat::BuildModelParamHandler();
	params <= "IP3Bias",  IP3Bias;
	unsigned int count = 0;
	for (unsigned int i = 0 ; i < cellsToStim.size() ; ++i)
		if (cellsToStim[i].stratName == GetClassName())
		{
			params <= (GetClassName() + "_Ind_" + StringifyFixed(count)), 
				cellsToStim[i].ind;
			params <= (GetClassName() + "_Start_" + StringifyFixed(count)), 
				cellsToStim[i].start;
			params <= (GetClassName() + "_End_" + StringifyFixed(count)), 
				cellsToStim[i].end;
			params <= (GetClassName() + "_Radius_" + StringifyFixed(count)), 
				cellsToStim[i].expand; 
			++count;
		}
	return params;
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool DefaultStimStrat::LoadFromStream(std::ifstream & stream)
{
	bool ok = StimulationStrat::LoadFromStream(stream);

	if (funct and freeFunct)
		delete funct;

	string tempName;
	stream >> tempName;
	funct = dynamic_cast<CouplingFunction*>(AbstractFactory<NetworkEdge>::Factories[tempName]
		->CreateFromStream(stream));
	freeFunct = true;

	stream >> IP3Bias;

	unsigned int size;
	stream >> size;

	cellsToStim.clear();
	for (unsigned int i = 0 ; i < size ; ++i)
	{
		StimulationInd temp(this->GetClassName(), 0, 0, 0, 0);
		stream >> temp.ind;
		stream >> temp.start;
		stream >> temp.end;
		stream >> temp.expand;
		cellsToStim.push_back(temp);
	}

	return ok and (not stream.eof());
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool DefaultStimStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = StimulationStrat::SaveToStream(stream);

	if (funct)
	{
		stream << funct->GetClassName() << endl;
		ok &= funct->SaveToStream(stream);
	}
	else
		return false;

	stream << IP3Bias << endl;

	unsigned int count = 0;
	for (unsigned int i = 0 ; i < cellsToStim.size() ; ++i)
		if (cellsToStim[i].stratName == GetClassName())
			++count;

	stream << count << endl;
	for (unsigned int i = 0 ; i < cellsToStim.size() ; ++i)
	{
		if (cellsToStim[i].stratName == GetClassName())
			stream
				<< cellsToStim[i].ind    << endl
				<< cellsToStim[i].start  << endl
				<< cellsToStim[i].end    << endl
				<< cellsToStim[i].expand << endl;
	}

	return ok;
}

//********************************************************************//
//************ P O I S S O N I A N   S T I M   S T R A T *************//
//********************************************************************//

//**********************************************************************
//Default Constructor from paramHandler
//**********************************************************************
PoissonianStimStrat::PoissonianStimStrat(ParamHandler & h)
{
	period = h.getParam<double>("-poissStim", 0);
	stimTime = h.getParam<double>("-poissStim", 1);
	IP3Bias = h.getParam<double>("-IPBias", 0);
	std::string couplFunctName = h.getParam<std::string>("-poissStim", 2);
	onlyDriveToSpike = h.getParam<bool>("-poissOnlySpike", 0);
	useCaSpontRelease = h.getParam<bool>("-poissCaSpontRel", 0);
	caThresh = h.getParam<double>("-activParams", 0);
	isolateNodes = h.getParam<bool>("-poissIsolateNode", 0);
	useGluStim = h.getParam<bool>("-poissGluStim", 0);
	gluQuantalRelease = h.getParam<double>("-poissGluStim", 1);
	poissOmegaC = h.getParam<double>("-poissGluStim", 2);
	
	assert(AbstractFactory<NetworkEdge>::Factories[couplFunctName]);
	NetworkEdge *netEdge = AbstractFactory<NetworkEdge>::Factories[
		couplFunctName]->Create();
	funct = dynamic_cast<CouplingFunction *>(netEdge);
	assert(funct);
	freeFunct = true;
}

//**********************************************************************
// Default Constructor
//**********************************************************************
PoissonianStimStrat::PoissonianStimStrat(double _T, double _s, 
	double bias, CouplingFunction *_f, bool _oS, bool _cSR, double _ct, 
	bool _in, bool _ugs, double _gqr, double _poc) : 
	period(_T), stimTime(_s), IP3Bias(bias), onlyDriveToSpike(_oS), 
	useCaSpontRelease(_cSR), caThresh(_ct), isolateNodes(_in), useGluStim(_ugs),
	gluQuantalRelease(_gqr), poissOmegaC(_poc), funct(_f), freeFunct(not _f)
{
	if (not funct)
		funct = new SigmoidCoupling();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
PoissonianStimStrat::PoissonianStimStrat(std::ifstream & stream) :
	funct(0)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
PoissonianStimStrat::~PoissonianStimStrat()
{
	if (freeFunct and funct)
		delete funct;
}

//**********************************************************************
// Stimulates the cells according to planned stimulations
//**********************************************************************
void PoissonianStimStrat::stimulateNet(StimulableCellNetwork & model)
{
	unsigned int nbCells = model.GetNbCells();

	if (cellTags.size() != nbCells)
		parseCellTags(model);

	if (spikeTrain.size() != nbCells)
		generateSpikeTrain(nbCells, model.GetTEnd());
	if (spikeInd.size() != nbCells)
		spikeInd = std::vector<unsigned int>(nbCells, 0);
	if (isolateNodes and gjcShutDown.size() != nbCells)
		gjcShutDown = std::vector<std::vector<int> >(nbCells, std::vector<int>(nbCells, 0));

	for (unsigned int i = 0 ; i < this->stimulated.size() ; ++i)
	{
		if (spikeInd[i] < spikeTrain[i].size())
		{
			// If the currend cell is being stimulated and should stop being stimulated
			if (this->stimulated[i] and ((this->tCurr > spikeTrain[i][spikeInd[i]].second)))
			{
				stopStimAndSetNextSpike(i);
				if (isolateNodes)
					unIsolateANode(i, model);
			}
			// If the current cell needs to be stimulated
			else if (not this->stimulated[i] and (this->tCurr > spikeTrain[i][spikeInd[i]].first))
			{
				this->stimulated[i] = true;
				if (isolateNodes)
					isolateANode(i, model);
			}
		}

		// Stimulate
		bool toStim = !(cellTags[i] & 0x1000000);
		if (toStim and this->stimulated[i] and not isolateNodes)
		{
			stimulateSpecific(&model, i);

			ChIModel *chimod;
			if ((chimod = dynamic_cast<ChIModel*>(&model)) and onlyDriveToSpike and (chimod->GetDynVal(i, ChICell::Ca) >= caThresh))
				stopStimAndSetNextSpike(i);
		}
	}
}

//**********************************************************************
// Actually stimulate with appropriate method
//**********************************************************************
void PoissonianStimStrat::stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const
{
	ChIModel *chimod;
	if ((chimod = dynamic_cast<ChIModel*>(model)))
	{
		if (useCaSpontRelease)
			chimod->StimSpontaneousCa(ind);
		else if (useGluStim)
		{
			const ChICell & cell = chimod->GetCell(ind);
			double Calc = chimod->GetDynVal(ind, ChICell::Ca);
			double glu = gluQuantalRelease * gsl_sf_exp(-poissOmegaC*(this->tCurr - spikeTrain[ind][spikeInd[ind]].first));
			chimod->ModifFluxes(ind, - cell.vbeta * (pow(glu, 0.7) / (pow(glu, 0.7) + 
				pow(cell.kR + cell.kP*(Calc / (Calc + cell.kpi)), 0.7))));
		}
		else
		{
			double diffIP3 = chimod->GetDynVal(ind, ChICell::IP3) - IP3Bias ;
			chimod->ModifFluxes(ind, (*funct)(diffIP3));
		}
	}
	else
		throw "Unkwnown model type";
}

//**********************************************************************
// Initialize the stimulation strategy as if it were just created
//**********************************************************************
void PoissonianStimStrat::Initialize()
{
	StimulationStrat::Initialize();
	spikeInd.clear();
	spikeTrain.clear();
	gjcShutDown.clear();
	cellTags.clear();
}

//**********************************************************************
// Initialize stimulation strategy but don't delete spikeTrain
//**********************************************************************
void PoissonianStimStrat::InitializeWithSameRandomSequence()
{
	StimulationStrat::Initialize();
	spikeInd.clear();
}

//**********************************************************************
// Get to next spike for a given cell
//**********************************************************************
void PoissonianStimStrat::stopStimAndSetNextSpike(unsigned int i)
{
	this->stimulated[i] = false;
	// Take the next spike
	while ((spikeInd[i] < spikeTrain[i].size()) and 
			(spikeTrain[i][spikeInd[i]].first < this->tCurr))
		++spikeInd[i];

}

//**********************************************************************
// Returns a param handler of the stimulation strat parameters
//**********************************************************************
ParamHandler PoissonianStimStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params += StimulationStrat::BuildModelParamHandler();
	params <= "PoissonianStimPeriod", period;
	params <= "PoissonianIP3Bias", IP3Bias;
	params <= "PoissonianStimMeanTime", stimTime;
	params <= "PoissonianCouplStrength", funct->GetRefOnStrength();
	params <= "PoissonianGluQuantalRelease", gluQuantalRelease;
	params <= "PoissonianGluOmegaC", poissOmegaC;

	params <= "ActivationThreshold", caThresh;
	return params;
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool PoissonianStimStrat::LoadFromStream(std::ifstream & stream)
{
	bool ok = StimulationStrat::LoadFromStream(stream);
	if (funct and freeFunct)
		delete funct;

	string tempName;
	stream >> tempName;
	funct = dynamic_cast<CouplingFunction*>(AbstractFactory<NetworkEdge>::Factories[tempName]
		->CreateFromStream(stream));
	freeFunct = true;

	stream >> period;
	stream >> stimTime;
	stream >> IP3Bias;
	stream >> onlyDriveToSpike;
	stream >> useCaSpontRelease;
	stream >> caThresh;
	stream >> isolateNodes;
	stream >> useGluStim;
	stream >> gluQuantalRelease;
	stream >> poissOmegaC;

	return ok and (not stream.eof());
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool PoissonianStimStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = StimulationStrat::SaveToStream(stream);

	if (funct)
	{
		stream << funct->GetClassName() << endl;
		ok &= funct->SaveToStream(stream);
	}
	else
		return false;

	stream 
		<< period   << endl
		<< stimTime << endl
		<< IP3Bias  << endl
		<< onlyDriveToSpike << endl
		<< useCaSpontRelease << endl
		<< caThresh << endl
		<< isolateNodes << endl
		<< useGluStim << endl
		<< gluQuantalRelease << endl
		<< poissOmegaC << endl;

	return ok;
}

//**********************************************************************
// Returns a random value from an exponential distribution of mean 'mean'
//**********************************************************************
double PoissonianStimStrat::expDelay(double mean)
{
	return ExpRand(mean);
}

//**********************************************************************
// Returns a time during which a cell is going to stay activated
//**********************************************************************
double PoissonianStimStrat::activationTime(double mean)
{
	return mean;
}

//**********************************************************************
// Generate a spike train
//**********************************************************************
void PoissonianStimStrat::generateSpikeTrain(unsigned int nbCells, double tMax)
{
	spikeTrain.clear();
	spikeTrain = std::vector<std::vector<std::pair<double, double> > >(
		nbCells, std::vector<std::pair<double, double> >());

	// For each cell
	for (unsigned int i = 0 ; i < nbCells ; ++i)
	{
		// If the cell needs to be stimulated

		spikeTrain[i].push_back(std::make_pair(-period * NB_POISS_START_DELAY, 0));
		// Add spikes separated by exponentially distributed delays until tMax is reached
		while (spikeTrain[i].back().first <= tMax)
		{
			spikeTrain[i].back().first += expDelay(period);
			spikeTrain[i].back().second = spikeTrain[i].back().first + activationTime(stimTime);
			if ((spikeTrain[i].back().first >= 0) and (spikeTrain[i].back().first < tMax))
				spikeTrain[i].push_back(std::make_pair(spikeTrain[i].back().first, 
					spikeTrain[i].back().first + activationTime(stimTime)));
		}
	}
}

//**********************************************************************
// Isolate a node
//**********************************************************************
void PoissonianStimStrat::isolateANode(unsigned int i, StimulableCellNetwork & model)
{
	const std::vector<unsigned int> & neighb = model.GetNeighbors(i);
	for (unsigned int j = 0 ; j < neighb.size() ; ++j)
	{
		CouplingFunction *cf = dynamic_cast<CouplingFunction*>(
			model.GetNetwork().GetAbstractEdge(i, neighb[j]));
		cf->ChangeStrength(0.0);
		gjcShutDown[i][neighb[j]]++;
		gjcShutDown[neighb[j]][i]++;
	}
}

//**********************************************************************
// UnIsolate a node
//**********************************************************************
void PoissonianStimStrat::unIsolateANode(unsigned int i, StimulableCellNetwork & model)
{
	const std::vector<unsigned int> & neighb = model.GetNeighbors(i);
	for (unsigned int j = 0 ; j < neighb.size() ; ++j)
	{
		gjcShutDown[i][neighb[j]]--;
		gjcShutDown[neighb[j]][i]--;
		if (gjcShutDown[i][neighb[j]] == 0)
		{
			model.GetNetwork().GetAbstractEdge(i, neighb[j])->UpdateEdge();
		}
	}
}

//**********************************************************************
// Parse cell tags from network
//**********************************************************************
void PoissonianStimStrat::parseCellTags(StimulableCellNetwork & model)
{
	cellTags.clear();
	const Network<CouplingFunction> & net = dynamic_cast<const Network<CouplingFunction> &>(model.GetNetwork());
	for(unsigned int i = 0 ; i < net.size() ; ++i)
		cellTags.push_back(net.GetNodeTag(i));
}

//********************************************************************//
//********** R A N D O M   M O N O   S T I M   S T R A T *************//
//********************************************************************//

//**********************************************************************
//Default Constructor from paramHandler
//**********************************************************************
RandomMonoStim::RandomMonoStim(ParamHandler & h) : currStimStart(0), 
	currPauseStart(0), currStimInd(0)
{
	stimLength = h.getParam<double>("-randStim", 0);
	pauseLength = h.getParam<double>("-randStim", 1);
	IP3Bias = h.getParam<double>("-IPBias", 0);
	std::string couplFunctName = h.getParam<std::string>("-randStim", 2);
	
	assert(AbstractFactory<NetworkEdge>::Factories[couplFunctName]);
	NetworkEdge *netEdge = AbstractFactory<NetworkEdge>::Factories[couplFunctName]->Create();
	funct = dynamic_cast<CouplingFunction *>(netEdge);
	assert(funct);
	freeFunct = true;
}

//**********************************************************************
// Default Constructor
//**********************************************************************
RandomMonoStim::RandomMonoStim(double _sl, double _pl, double bias, 
	CouplingFunction *_f) : 
	stimLength(_sl), pauseLength(_pl), currStimStart(0), 
	currPauseStart(0), currStimInd(0), IP3Bias(bias), funct(_f), 
	freeFunct(not _f)
{
	if (not funct)
		funct = new SigmoidCoupling();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
RandomMonoStim::RandomMonoStim(std::ifstream & stream) : funct(0), 
	freeFunct(false)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
RandomMonoStim::~RandomMonoStim()
{
	if (freeFunct and funct)
		delete funct;
}

//**********************************************************************
// Stimulates the cells according to planend stimulations
//**********************************************************************
void RandomMonoStim::stimulateNet(StimulableCellNetwork & model)
{
	if (stimTrain.size() != model.GetNbCells())
		generateStimTrain(model.GetNbCells(), model.GetTEnd());

	double t = this->tCurr;
	if (this->stimulated[stimTrain[currStimInd]] and (currStimStart + stimLength < t))
	{
		currPauseStart = t;
		this->stimulated[stimTrain[currStimInd]] = false;
	}
	else if (not this->stimulated[stimTrain[currStimInd]] and (currPauseStart + pauseLength < t))
	{
		currStimStart = t;
		currStimInd++;
		this->stimulated[stimTrain[currStimInd]] = true;
	}

	if (this->stimulated[stimTrain[currStimInd]])
	{
		stimulateSpecific(&model, stimTrain[currStimInd]);
	}
}

//**********************************************************************
// Actually stimulate with appropriate method
//**********************************************************************
void RandomMonoStim::stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const
{
	ChIModel *chimod;
	if ((chimod = dynamic_cast<ChIModel*>(model)))
	{
		double diffIP3 = chimod->GetDynVal(ind, ChICell::IP3) - IP3Bias;
		chimod->ModifFluxes(ind, (*funct)(diffIP3));
	}
	else
		throw "Unkwnown model type";
}

//**********************************************************************
// Generates stimulation train
//**********************************************************************
void RandomMonoStim::generateStimTrain(unsigned int nbCells, double tMax)
{
	double currT = 0;
	stimTrain.clear();
	while(currT <= tMax)
	{
		stimTrain.push_back(floor(UnifRand() * (double)nbCells));
		currT += stimLength + pauseLength;
	}
}

//**********************************************************************
// Initialize the stimulation strategy as if it were just created
//**********************************************************************
void RandomMonoStim::Initialize()
{
	StimulationStrat::Initialize();
	currStimStart = 0;
	currPauseStart = 0;
	currStimInd = 0;
	stimTrain.clear();
}

//**********************************************************************
// Initialize stimulation strategy but don't delete spikeTrain
//**********************************************************************
void RandomMonoStim::InitializeWithSameRandomSequence()
{
	StimulationStrat::Initialize();
	currStimStart = 0;
	currPauseStart = 0;
	currStimInd = 0;
}

//**********************************************************************
// Returns a param handler of the stimulation strat parameters
//**********************************************************************
ParamHandler RandomMonoStim::BuildModelParamHandler() 
{
	ParamHandler params;
	params += StimulationStrat::BuildModelParamHandler();
	params <= "RandomStimLength", stimLength;
	params <= "RandomStimPause", pauseLength;
	params <= "IP3Bias", IP3Bias;
	return params;
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool RandomMonoStim::LoadFromStream(std::ifstream & stream)
{
	bool ok = StimulationStrat::LoadFromStream(stream);

	if (funct and freeFunct)
		delete funct;

	string tempName;
	stream >> tempName;
	funct = dynamic_cast<CouplingFunction*>(AbstractFactory<NetworkEdge>::Factories[tempName]
		->CreateFromStream(stream));
	freeFunct = true;

	stream >> stimLength;
	stream >> pauseLength;
	stream >> currStimStart;
	stream >> currPauseStart;
	stream >> IP3Bias;

	return ok and (not stream.eof());
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool RandomMonoStim::SaveToStream(std::ofstream & stream) const
{
	bool ok = StimulationStrat::SaveToStream(stream);

	if (funct)
	{
		stream << funct->GetClassName() << endl;
		ok &= funct->SaveToStream(stream);
	}
	else
		return false;

	stream
		<< stimLength << endl
		<< pauseLength << endl
		<< currStimStart << endl
		<< currPauseStart << endl
		<< IP3Bias << endl;

	return ok;
}

//********************************************************************//
//******* T H R E S H O L D   D E T E R M   S T I M   S T R A T ******//
//********************************************************************//

//**********************************************************************
//Default Constructor from paramHandler
//**********************************************************************
ThresholdDeterminationStimStrat::ThresholdDeterminationStimStrat(ParamHandler & h) : sinkFunct(0)
{
	nbNeighbSinks = h.getParam<unsigned int>("-ThreshDeterm", 0);
	nbStimulated = h.getParam<unsigned int>("-ThreshDeterm", 1);
	CaThresh = h.getParam<double>("-ThreshDeterm", 2);
	IP3Bias = h.getParam<double>("-IPBias", 0);
	std::string couplFunctName = h.getParam<std::string>("-ThreshDeterm", 3);
	double couplStr = h.getParam<double>("-ThreshDeterm", 4);
	useCaSpontRelease = h.getParam<bool>("-ThreshDetermCaSpontRel", 0);
	
	assert(AbstractFactory<NetworkEdge>::Factories[couplFunctName]);
	NetworkEdge *netEdge = AbstractFactory<NetworkEdge>::Factories[couplFunctName]->Create();
	funct = dynamic_cast<CouplingFunction *>(netEdge);
	assert(funct);
	funct->ChangeStrength(couplStr);
	funct->SetConstStrength(true);
	freeFunct = true;
}
//**********************************************************************
// Default Constructor
//**********************************************************************
ThresholdDeterminationStimStrat::ThresholdDeterminationStimStrat(
	unsigned int ns, unsigned int nb, double bias, CouplingFunction *_f, 
	double _Cat, bool _spontCa) : nbStimulated(nb), nbNeighbSinks(ns), IP3Bias(bias), 
		CaThresh(_Cat), useCaSpontRelease(_spontCa), funct(_f), freeFunct(not _f),
		sinkFunct(0)
{
	if (not funct)
		funct = new SigmoidCoupling();
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ThresholdDeterminationStimStrat::ThresholdDeterminationStimStrat(
	std::ifstream & stream) : funct(0), freeFunct(false), sinkFunct(0)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
ThresholdDeterminationStimStrat::~ThresholdDeterminationStimStrat()
{
	if (freeFunct and funct)
		delete funct;
}

//**********************************************************************
// Stimulates the cells according to planend stimulations
//**********************************************************************
void ThresholdDeterminationStimStrat::stimulateNet(StimulableCellNetwork & model)
{
	ChIModelThresholdDetermination *model2 = dynamic_cast<ChIModelThresholdDetermination *>(&model);
	assert(model2);

	if (neverSpiked.empty())
		neverSpiked = std::vector<bool>(model2->GetNbCells(), true);
	if (not sinkFunct)
	{
		const Network<CouplingFunction> *cfn = dynamic_cast<const Network<CouplingFunction>*>(&(model.GetNetwork()));
		sinkFunct = cfn->GetAllocatedLinks(false)[0];
	}

	unsigned int nbDeriv = model2->GetNbDerivs();

	// For all cells
	for (unsigned int i = 1 ; i < model2->GetNbCells() ; ++i)
	{
		// If cell i needs to be stimulated
		if (i <= nbStimulated/* and neverSpiked[i]*/)
		{
			// Stimulate it until CaThresh threshold is crossed (until it fires)
			if ((model2->GetDynVal(i, ChICell::Ca) < CaThresh) and neverSpiked[i])
			{
				if (useCaSpontRelease)
					model2->StimSpontaneousCa(i);
				else
				{
					double diffIP3 = model2->GetDynVal(i, ChICell::IP3) - IP3Bias;
					model2->ModifFluxes(i, (*funct)(diffIP3));
				}
				this->stimulated[i] = true;
			}
			else
			{
				this->stimulated[i] = false;
				neverSpiked[i] = false;
			}
		}
		// Otherwise, if i is a "sink" cell
		else if ((i <= model2->GetNetwork().GetNodeDegree(0)) or (i > (3+nbDeriv)*model2->GetNetwork().GetNodeDegree(0)))
		{
			double diffIP3 = model2->GetDynVal(i, ChICell::IP3) - ChICell::DefaultIP3;
			// the function linking sinks to a virtual cell should be the same as the coupling function
			// in the network.
			model2->ModifFluxes(i, (*sinkFunct)(diffIP3) * 2.0);
		}
	}
}

//**********************************************************************
// Initialize the stimulation strategy as if it were just created
//**********************************************************************
void ThresholdDeterminationStimStrat::Initialize()
{
	StimulationStrat::Initialize();
	neverSpiked.clear();
	sinkFunct = 0;
}

//**********************************************************************
// Returns a param handler of the stimulation strat parameters
//**********************************************************************
ParamHandler ThresholdDeterminationStimStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params += StimulationStrat::BuildModelParamHandler();
	params <= "ThreshDetermNbStim", nbStimulated;
	params <= "ThreshDetermNbSinks", nbNeighbSinks;
	params <= "IP3Bias", IP3Bias;
	params <= "ThreshDetermMaxCa", CaThresh;
	return params;
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool ThresholdDeterminationStimStrat::LoadFromStream(std::ifstream & stream)
{
	bool ok = StimulationStrat::LoadFromStream(stream);

	if (funct and freeFunct)
		delete funct;

	string tempName;
	stream >> tempName;
	funct = dynamic_cast<CouplingFunction*>(AbstractFactory<NetworkEdge>::Factories[tempName]
		->CreateFromStream(stream));
	freeFunct = true;

	stream >> nbStimulated;
	stream >> nbNeighbSinks;
	stream >> CaThresh;
	stream >> IP3Bias;
	stream >> useCaSpontRelease;

	return ok and (not stream.eof());
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool ThresholdDeterminationStimStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = StimulationStrat::SaveToStream(stream);

	if (funct)
	{
		stream << funct->GetClassName() << endl;
		ok &= funct->SaveToStream(stream);
	}
	else
		return false;

	stream
		<< nbStimulated << endl
		<< nbNeighbSinks << endl
		<< CaThresh << endl
		<< IP3Bias << endl
		<< useCaSpontRelease << endl;

	return ok;
}

//********************************************************************//
//************ M E A N   P A T H   S T I M   S T R A T ***************//
//********************************************************************//

//**********************************************************************
//Default Constructor from paramHandler
//**********************************************************************
MeanPathStimStrat::MeanPathStimStrat(ParamHandler & h) : 
	stimCellsChosen(false), currEffMeanPath(-1)
{
	nbStims = h.getParam<unsigned int>("-MeanPathStimParams", 0);
	meanPathLength = h.getParam<double>("-MeanPathStimParams", 1);
	tStart = h.getParam<double>("-MeanPathStimParams", 2);
	tEnd = h.getParam<double>("-MeanPathStimParams", 3);
	stimExp = h.getParam<unsigned int>("-MeanPathStimParams", 4);
	IP3Bias = h.getParam<double>("-IPBias", 0);
	std::string couplFunctName = h.getParam<std::string>("-MeanPathStimParams", 5);

	assert(AbstractFactory<NetworkEdge>::Factories[couplFunctName]);
	NetworkEdge *netEdge = AbstractFactory<NetworkEdge>::Factories[couplFunctName]->Create();
	funct = dynamic_cast<CouplingFunction *>(netEdge);
	assert(funct);
	freeFunct = true;
}

//**********************************************************************
// Default Constructor
//**********************************************************************
MeanPathStimStrat::MeanPathStimStrat(unsigned int _nS, double _mpl, 
	double _ts,	double _te, unsigned int _se, double bias, 
	CouplingFunction *_f) :
	DefaultStimStrat::DefaultStimStrat(bias, _f), nbStims(_nS), 
	meanPathLength(_mpl), tStart(_ts), tEnd(_te), stimExp(_se), 
	stimCellsChosen(false), currEffMeanPath(-1)
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
MeanPathStimStrat::MeanPathStimStrat(std::ifstream & stream) : 
	DefaultStimStrat::DefaultStimStrat(), stimCellsChosen(false), 
	currEffMeanPath(-1)
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
MeanPathStimStrat::~MeanPathStimStrat()
{
}

//**********************************************************************
// Stimulates the cells according to planend stimulations
//**********************************************************************
void MeanPathStimStrat::stimulateNet(StimulableCellNetwork & model)
{
	if (not stimCellsChosen)
	{
		AllPairDistances *allPairDist = 
			GetSpecificMetric<Metric, AllPairDistances>(model.GetAllMetrics());
		assert(allPairDist);
		const std::vector<std::vector<double> > & distances = 
			allPairDist->GetDistances();

		std::vector<unsigned int> stimCells;
		std::vector<double> predDistToMeanPath(model.GetNbCells(), 0);
		// Add randomly the first stimulated cells
		stimCells.push_back(floor(UnifRand() * model.GetNbCells()));
		for (unsigned int i = 1 ; i < nbStims ; ++i)
		{
			// Compute the current total path inside the selected cells
			double currMeanTotPath = 0;
			for (unsigned int m = 0 ; m < stimCells.size() ; ++m)
				for (unsigned int p = m+1 ; p < stimCells.size() ; ++p)
					currMeanTotPath += distances[stimCells[m]][stimCells[p]];

			// For each candidate cell
			for (unsigned int j = 0 ; j < model.GetNbCells() ; ++j)
			{
				// Compute the mean path that would result from 
				// adding the current cell to the selected cells
				double tempMeanPath = currMeanTotPath;
				for (unsigned int m = 0 ; m < stimCells.size() ; ++m)
					tempMeanPath += distances[stimCells[m]][j];

				tempMeanPath /= (double) stimCells.size() * 
					((double) stimCells.size() + 1.0) / 2.0;
				predDistToMeanPath[j] = 
					fabs(meanPathLength - tempMeanPath) + UnifRand() * MEAN_PATH_PRECISION;
			}
			// Don't take into account the cells that have already been added
			for (unsigned int k = 0 ; k < stimCells.size() ; ++k)
				predDistToMeanPath[stimCells[k]] = DEFAULT_MAX_VAL;

			// Select the cell that will affect the least the mean path length
			stimCells.push_back(GetIndiceMin(predDistToMeanPath));
		}
		// Compute the effective mean path between the selected cells
		double currMeanTotPath = 0;
		for (unsigned int m = 0 ; m < stimCells.size() ; ++m)
			for (unsigned int p = m+1 ; p < stimCells.size() ; ++p)
				currMeanTotPath += distances[stimCells[m]][stimCells[p]];
		currEffMeanPath = currMeanTotPath / ((double) stimCells.size() * 
					((double) stimCells.size() - 1.0) / 2.0);

		// Add all selected cells to DefaultStimStrat for effective stimulations
		for (unsigned int i = 0 ; i < stimCells.size() ; ++i)
			AddCellStim(GetClassName(), stimCells[i], tStart, tEnd, stimExp);

		stimCellsChosen = true;
	}

	DefaultStimStrat::stimulate(model);
}

//**********************************************************************
// Initialize the stimulation strategy as if it were just created
//**********************************************************************
void MeanPathStimStrat::Initialize()
{
	this->cellsToStim.clear();
	DefaultStimStrat::Initialize();
	stimCellsChosen = false;
	currEffMeanPath = -1;
}

//**********************************************************************
// Returns a param handler of the stimulation strat parameters
//**********************************************************************
ParamHandler MeanPathStimStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params += DefaultStimStrat::BuildModelParamHandler();

	params <= "MeanPathStimNbStims", nbStims;
	params <= "MeanPathStimMeanPathLength", meanPathLength;
	params <= "MeanPathStimTStart", tStart;
	params <= "MeanPathStimTEnd", tEnd;
	params <= "MeanPathStimExp", stimExp;

	return params;
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool MeanPathStimStrat::LoadFromStream(std::ifstream & stream)
{
	bool ok = DefaultStimStrat::LoadFromStream(stream);

	stream >> nbStims;
	stream >> meanPathLength;
	stream >> tStart;
	stream >> tEnd;
	stream >> stimExp;

	return ok and (not stream.eof());
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool MeanPathStimStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = DefaultStimStrat::SaveToStream(stream);

	stream 
		<< nbStims << std::endl
		<< meanPathLength << std::endl
		<< tStart << std::endl
		<< tEnd << std::endl
		<< stimExp << std::endl;

	return ok and stream.good();
}

//********************************************************************//
//********* C O R R E L   D E T E R M   S T I M   S T R A T **********//
//********************************************************************//

//**********************************************************************
		//Default Constructor from paramHandler
//**********************************************************************
CorrelDetermStimStrat::CorrelDetermStimStrat(ParamHandler & h)
{
	pauseTime = ParamHandler::GlobalParams.getParam<double>("-CorrelStim", 0);
	IP3Bias = h.getParam<double>("-IPBias", 0);
	std::string couplFunctName = h.getParam<std::string>("-CorrelStim", 1);

	assert(AbstractFactory<NetworkEdge>::Factories[couplFunctName]);
	NetworkEdge *netEdge = AbstractFactory<NetworkEdge>::Factories[couplFunctName]->Create();
	funct = dynamic_cast<CouplingFunction *>(netEdge);
	assert(funct);
	freeFunct = true;
}

//**********************************************************************
// Default Constructor
//**********************************************************************
CorrelDetermStimStrat::CorrelDetermStimStrat(double _ps, double bias, 
	CouplingFunction *_f) :
	DefaultStimStrat::DefaultStimStrat(bias, _f), pauseTime(_ps)
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
CorrelDetermStimStrat::CorrelDetermStimStrat(std::ifstream & stream) : 
	DefaultStimStrat::DefaultStimStrat()
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
CorrelDetermStimStrat::~CorrelDetermStimStrat()
{
}

//**********************************************************************
// Stimulates the cells according to planend stimulations
//**********************************************************************
void CorrelDetermStimStrat::stimulateNet(StimulableCellNetwork & model)
{
	// If this is the first call of the simulation,
	// initialize the stimulations.
	if (resetTimes.empty())
	{
		this->cellsToStim.clear();
		double ts = model.GetTStart();
		double te = model.GetTEnd();
		unsigned int nbCells = model.GetNbCells();
		unsigned int nbCellsToStim = 0;
		// Counting how many nodes are going to be stimulated (exclude
		// nodes that might introduce boundary effects)
		for (unsigned int i = 0 ; i < nbCells ; ++i)
			if (not model.GetNetwork().IsNodeOnEdge(i))
				++nbCellsToStim;

		double stimTime = (te - ts) / ((double)nbCellsToStim);
		double stimStart, stimEnd;
		unsigned int nbStimCurr = 0;
		for (unsigned int i = 0 ; i < nbCells ; ++i)
		{
			if (not model.GetNetwork().IsNodeOnEdge(i))
			{
				stimStart = ts + ((double)nbStimCurr)*stimTime;
				stimEnd = ts + ((double)nbStimCurr + 1.0)*stimTime;
				++nbStimCurr;

				assert(stimEnd - pauseTime > stimStart);
				AddCellStim(GetClassName(), i, stimStart, stimEnd - pauseTime, 0);
				resetTimes.push_back(stimEnd - pauseTime);
			}
		}
		nextResetTimeInd = 0;
	}

	DefaultStimStrat::stimulate(model);

	// If we need to reset all the cells (end of a stimulation)
	if (this->tCurr > resetTimes[nextResetTimeInd])
	{
		model.SetAllCellsToEquilibrium();
		++nextResetTimeInd;
	}
}

//**********************************************************************
// Initialize the stimulation strategy as if it were just created
//**********************************************************************
void CorrelDetermStimStrat::Initialize()
{
	this->cellsToStim.clear();
	DefaultStimStrat::Initialize();
	resetTimes.clear();
	nextResetTimeInd = 0;
}

//**********************************************************************
// Returns a param handler of the stimulation strat parameters
//**********************************************************************
ParamHandler CorrelDetermStimStrat::BuildModelParamHandler() 
{
	ParamHandler params;
	params += DefaultStimStrat::BuildModelParamHandler();
	params <= "CorrelDetermPauseTime", pauseTime;
	return params;
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool CorrelDetermStimStrat::LoadFromStream(std::ifstream & stream)
{
	bool ok = DefaultStimStrat::LoadFromStream(stream);

	stream >> pauseTime;

	return ok and (not stream.eof());
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool CorrelDetermStimStrat::SaveToStream(std::ofstream & stream) const
{
	bool ok = DefaultStimStrat::SaveToStream(stream);

	stream << pauseTime << std::endl;

	return ok and stream.good();
}

//********************************************************************//
//******* E X T R A   C E L L   G L U T A M A T E   S T I M **********//
//********************************************************************//

//**********************************************************************
//Default Constructor from paramHandler
//**********************************************************************
ExtraCellGluStim::ExtraCellGluStim(ParamHandler & h) : 
	DefaultStimStrat(h)
{
	glu = h.getParam<double>("-XtraCellGluStim");
}

//**********************************************************************
// Default Constructor
//**********************************************************************
ExtraCellGluStim::ExtraCellGluStim(double _glu, CouplingFunction *_f) :
	DefaultStimStrat::DefaultStimStrat(0, _f), glu(_glu)
{
}

//**********************************************************************
// Constructor from stream
//**********************************************************************
ExtraCellGluStim::ExtraCellGluStim(std::ifstream & stream) : 
	DefaultStimStrat::DefaultStimStrat()
{
	LoadFromStream(stream);
}

//**********************************************************************
// Destructor
//**********************************************************************
ExtraCellGluStim::~ExtraCellGluStim()
{
}

//**********************************************************************
// Returns a param handler of the stimulation strat parameters
//**********************************************************************
ParamHandler ExtraCellGluStim::BuildModelParamHandler() 
{
	ParamHandler params;
	params += DefaultStimStrat::BuildModelParamHandler();

	params <= "ExtraCellGluStim", glu;

	return params;
}

//**********************************************************************
// Load the object from a stream
//**********************************************************************
bool ExtraCellGluStim::LoadFromStream(std::ifstream & stream)
{
	bool ok = DefaultStimStrat::LoadFromStream(stream);

	stream >> glu;

	return ok and (not stream.eof());
}

//**********************************************************************
// Saves the object to a stream
//**********************************************************************
bool ExtraCellGluStim::SaveToStream(std::ofstream & stream) const
{
	bool ok = DefaultStimStrat::SaveToStream(stream);

	stream 
		<< glu << std::endl;

	return ok and stream.good();
}

//**********************************************************************
// Actually stimulate with appropriate method
//**********************************************************************
void ExtraCellGluStim::stimulateSpecific(StimulableCellNetwork *model, unsigned int ind) const
{
	ChIModel *chimod;
	if ((chimod = dynamic_cast<ChIModel*>(model)))
	{
		const ChICell & cell = chimod->GetCell(ind);
		double Calc = chimod->GetDynVal(ind, ChICell::Ca);
		chimod->ModifFluxes(ind, - cell.vbeta * (pow(glu, 0.7) / (pow(glu, 0.7) + 
			pow(cell.kR + cell.kP*(Calc / (Calc + cell.kpi)), 0.7))));
	}
	else
		throw "Unkwnown model type";
}

