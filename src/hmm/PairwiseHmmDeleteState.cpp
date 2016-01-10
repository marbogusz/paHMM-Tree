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


#include "hmm/PairwiseHmmDeleteState.hpp"
#include "hmm/DpMatrixFull.hpp"

namespace EBC
{

PairwiseHmmDeleteState::PairwiseHmmDeleteState(unsigned int x, unsigned int y)
{
	this->rows =x;
	this->cols =y;
	this->dpMatrix = new DpMatrixFull(x,y);
	stateId = Definitions::StateId::Delete;
	//initializeData();
}

PairwiseHmmDeleteState::PairwiseHmmDeleteState(DpMatrixBase *matrix)
{
	this->dpMatrix = matrix;
	//initializeData();
}

void PairwiseHmmDeleteState::initializeData(double lnl, bool backwards)
{

	if(!backwards)
		dpMatrix->setValue(0,0,lnl);
	/*
	else
	{
		dpMatrix->setWholeCol(this->cols-1,-100000);
		dpMatrix->setWholeRow(this->rows-1,-100000);
		dpMatrix->setValue(rows-1,cols-1,0);
	}
	*/
}

void PairwiseHmmDeleteState::setDirection(unsigned int i, unsigned int j)
{
	dpMatrix->setHorizontalAt(i,j);
}

PairwiseHmmDeleteState::~PairwiseHmmDeleteState()
{
}

} /* namespace EBC */
