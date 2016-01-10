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


#ifndef DPLOMEMMATRIX_H_
#define DPLOMEMMATRIX_H_

#include <limits>
#include <iostream>
#include "hmm/DpMatrixBase.hpp"

using namespace std;

namespace EBC
{

class DpMatrixLoMem : public DpMatrixBase
{

protected:

	void allocateData();
	//Two lines of data!
	double* buffer[2];

	//First and second row pointers
	double* previousRow;
	double* currentRow;

	unsigned int currentRowIndex;
	unsigned int nextRowIndex;
	int previousRowIndex;

	void setValue(unsigned int col, double value);

	double valueAtColumn(unsigned int col);

	double valueAtLeft(unsigned int col);

	double valueAtTop(unsigned int col);

	double valueAtDiagonal(unsigned int col);

	inline void nextRow();

	void clear();

public:

	void setValue(unsigned int x,unsigned int y, double value);

	double valueAt(unsigned int i, unsigned int j);

	void setSrc(unsigned int i, unsigned int j, DpMatrixBase*);

	void setDiagonalAt(unsigned int i, unsigned int j);

	void setHorizontalAt(unsigned int i, unsigned int j);

	void setVerticalAt(unsigned int i, unsigned int j);

	void setWholeRow(unsigned int row, double value);

	void setWholeCol(unsigned int col, double value);

	void traceback(string& seq_a, string& seq_b, std::pair<string,string>* alignment) {}

	void tracebackRaw(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, vector<std::pair<unsigned int, unsigned int> >&) {}


	DpMatrixLoMem(unsigned int xSize, unsigned int ySize);

	virtual ~DpMatrixLoMem();




};

} /* namespace EBC */
#endif /* DPREDUCEDMATRIX_H_ */
