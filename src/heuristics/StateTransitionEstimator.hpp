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


#ifndef STATETRANSITIONESTIMATOR_HPP_
#define STATETRANSITIONESTIMATOR_HPP_

#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"

#include "heuristics/StateTransitionML.hpp"

namespace EBC
{

class StateTransitionEstimator : public IOptimizable
{
protected:

	Optimizer* bfgs;

	OptimizedModelParameters* modelParams;

	vector<StateTransitionML*> stmSamples;

	IndelModel* indelModel;

	Maths* maths;

	unsigned char gapCharacter;

	double maxTime;

	bool useStateEq;


public:
	StateTransitionEstimator(IndelModel* im, Definitions::OptimizationType ot, unsigned int pairCategories, unsigned char gapchar, bool useEq = true);

	void addTime(double time, unsigned int triplet, unsigned int pr);

	void addPair(vector<unsigned char>* s1, vector<unsigned char>* s2, unsigned int triplet, unsigned int pr);

	double runIteration();

	void optimize();

	void clean();

	virtual ~StateTransitionEstimator();

	OptimizedModelParameters* getModelParams()
	{
		return modelParams;
	}

	IndelModel* getIndelModel()
	{
		return indelModel;
	}
};

} /* namespace EBC */

#endif /* STATETRANSITIONESTIMATOR_HPP_ */
