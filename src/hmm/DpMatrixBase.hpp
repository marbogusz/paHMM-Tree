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


#ifndef DPMATRIXBASE_H_
#define DPMATRIXBASE_H_

#include <limits>
#include <iostream>
#include <vector>

#include "core/SequenceElement.hpp"
#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{

class DpMatrixBase
{

protected:

	unsigned int xSize, ySize;

	double minVal;

	virtual void allocateData()=0;
public:

	virtual void setWholeRow(unsigned int row, double value)=0;

	virtual void setWholeCol(unsigned int col, double value)=0;

	DpMatrixBase(unsigned int xSize, unsigned int ySize)
	{
		this->xSize = xSize;
		this->ySize = ySize;
		this->minVal = Definitions::minMatrixLikelihood;
	}

	virtual ~DpMatrixBase() {}

	virtual void setValue(unsigned int x,unsigned int y, double value)=0;

	virtual double valueAt(unsigned int i, unsigned int j)=0;

	virtual void setSrc(unsigned int i, unsigned int j, DpMatrixBase*)=0;

	virtual void setDiagonalAt(unsigned int i, unsigned int j)=0;

	virtual void setHorizontalAt(unsigned int i, unsigned int j)=0;

	virtual void setVerticalAt(unsigned int i, unsigned int j)=0;

	virtual void traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment)=0;

	virtual void tracebackRaw(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, vector<std::pair<unsigned int, unsigned int> >&)=0;
};

} /* namespace EBC */
#endif /* DPMATRIXBASE_H_ */
