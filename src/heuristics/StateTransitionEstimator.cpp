//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================


#include "heuristics/StateTransitionEstimator.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{

StateTransitionEstimator::StateTransitionEstimator(IndelModel* im, Definitions::OptimizationType ot, unsigned int pc, unsigned char gc, bool useEq) :
		indelModel(im), stmSamples(pc), gapCharacter(gc), useStateEq(useEq)
{
	DEBUG("Starting State Transition Estimator");
	//indelModel = new NegativeBinomialGapModel();
	maths = new Maths();

	modelParams = new OptimizedModelParameters(NULL, indelModel,0, 0, false,
	true, false, false, maths);
	modelParams->useIndelModelInitialParameters();

	bfgs = new Optimizer(modelParams, this,ot);

	maxTime = Definitions::almostZero;
}

double StateTransitionEstimator::runIteration()
{
	double result = 0;

	//modelParams->outputParameters();
	indelModel->setParameters(modelParams->getIndelParameters());
	//go through matrices
	for(auto tm : stmSamples)
	{
		//set parameter for every sample
		result += tm->getLnL();
		//calculate and add to result
	}
	return result * -1.0;
}

void StateTransitionEstimator::addTime(double time, unsigned int triplet, unsigned int pr)
{

	if (time > maxTime)
		maxTime = time;
	DEBUG("State Transition Estimator add time for triplet " << triplet << "\t pair " << pr << "\ttime " << time);
	stmSamples[2*triplet+pr] = new StateTransitionML(indelModel, time, gapCharacter, useStateEq);
}

void StateTransitionEstimator::addPair(vector<unsigned char>* s1,
		vector<unsigned char>* s2, unsigned int triplet, unsigned int pr)
{
	DUMP("State Transition Estimator add pair for triplet " << triplet << " and pair no " << pr );
	stmSamples[2*triplet+pr]->addSample(s1,s2);
}

void StateTransitionEstimator::optimize()
{
	modelParams->boundLambdaBasedOnDivergence(maxTime);
	bfgs->optimize();
	indelModel->setParameters(modelParams->getIndelParameters());
	INFO("StateTransitionEstimator results:");
	modelParams->logParameters();
	//modelParams->outputParameters();
}

void StateTransitionEstimator::clean()
{
    for(auto tm : stmSamples)
    {
    	delete tm;
    }
}

StateTransitionEstimator::~StateTransitionEstimator()
{
	delete bfgs;
	delete modelParams;
    delete maths;
    //for(auto tm : stmSamples)
    //{
    //	delete tm;
    //}
    //delete indelModel;
}

} /* namespace EBC */
