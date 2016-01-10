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


#ifndef STATETRANSITIONMATRIX_HPP_
#define STATETRANSITIONMATRIX_HPP_

#include "models/IndelModel.hpp"
#include "core/Definitions.hpp"
#include "core/TransitionProbabilities.hpp"
#include "core/SequenceElement.hpp"

namespace EBC
{


//This is a 3x3 state transition matrix
class StateTransitionML
{
protected:

	TransitionProbabilities* tpb;
	double time;

	unsigned int matrixSize;

	//transition probabilities
	double g;
	double e;

	double md[Definitions::stateCount][Definitions::stateCount];
	//match,insert, delete
	unsigned int counts[Definitions::stateCount][Definitions::stateCount];
	double pis[Definitions::stateCount];

	unsigned char isGap;

	bool useStateEq;

	Definitions::StateId firstState;

public:
	virtual ~StateTransitionML();

	//get equilibrium frequencies
	void calculatePIs();

	void calculateParameters();

	StateTransitionML(IndelModel* im, double time, unsigned char, bool);

	void addSample(vector<unsigned char>*, vector<unsigned char>* s2);

	double getLnL();
};

} /* namespace EBC */

#endif /* STATETRANSITIONMATRIX_HPP_ */
