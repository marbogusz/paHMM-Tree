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

#ifndef DPMATRIXFULL_H_
#define DPMATRIXFULL_H_


#include "hmm/DpMatrixBase.hpp"
#include "core/Definitions.hpp"

#include <vector>
#include <limits>
#include <iostream>

using namespace std;

namespace EBC
{

class DpMatrixFull : public DpMatrixBase
{
	struct TraceStep
	{
		//3 traceback elements
		bool vert;
		bool hor;
		bool diag;
		//score/probability
		double score;
		DpMatrixFull* src;

		TraceStep() : vert(false), hor(false),
			diag(false), score(-100000.0), src(NULL) {}
	};

protected:

	void allocateData();

public:

	void setWholeRow(unsigned int row, double value);

	void setWholeCol(unsigned int col, double value);

	double  ** matrixData;
	//TraceStep ** matrixData;

	void traceback(string& seqA, string& seqB, std::pair<string,string>* alignment);

	void tracebackRaw(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, vector<std::pair<unsigned int, unsigned int> >&);

	DpMatrixFull(unsigned int xSize, unsigned int ySize);

	virtual ~DpMatrixFull();

	void setValue(unsigned int x,unsigned int y, double value);

	double valueAt(unsigned int i, unsigned int j);

	void setDiagonalAt(unsigned int i, unsigned int j);

	void setHorizontalAt(unsigned int i, unsigned int j);

	void setVerticalAt(unsigned int i, unsigned int j);

	void setSrc(unsigned int i, unsigned int j, DpMatrixBase*);

	void outputTrace(unsigned int);

	void outputValues(unsigned int);

	void outputValuesWithBands(const vector<pair<int, int> >& band, const vector<pair<int, int> >& oband1, const vector<pair<int, int> >& oband2, char os1, char os2);


};

} /* namespace EBC */
#endif /* EVALUATIONMATRIX_H_ */
