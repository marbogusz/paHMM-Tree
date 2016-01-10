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


#include "hmm/PairwiseHmmMatchState.hpp"
#include "hmm/DpMatrixFull.hpp"

namespace EBC
{

PairwiseHmmMatchState::PairwiseHmmMatchState(unsigned int x, unsigned int y)
{
	this->rows =x;
	this->cols =y;
    this->dpMatrix = new DpMatrixFull(x,y);
    stateId = Definitions::StateId::Match;
	//initializeData();
}

PairwiseHmmMatchState::PairwiseHmmMatchState(DpMatrixBase *matrix)
{
	this->dpMatrix = matrix;
	//initializeData();
}

PairwiseHmmMatchState::~PairwiseHmmMatchState()
{
}

void PairwiseHmmMatchState::initializeData(double lnl, bool backwards)
{
	if (!backwards)
	{
		//dpMatrix->setWholeCol(0,-100000);
		//dpMatrix->setWholeRow(0,-100000);
		dpMatrix->setValue(0,0,lnl);
	}
	/*else
	{
		dpMatrix->setWholeCol(this->cols-1,-100000);
		dpMatrix->setWholeRow(this->rows-1,-100000);
		dpMatrix->setValue(rows-1,cols-1,0);
	}
	*/
}

void PairwiseHmmMatchState::setDirection(unsigned int i, unsigned int j)
{
	this->dpMatrix->setDiagonalAt(i,j);
}

} /* namespace EBC */
