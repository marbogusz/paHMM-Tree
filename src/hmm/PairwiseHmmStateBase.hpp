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

#ifndef PAIRWISEHMMSTATEBASE_H_
#define PAIRWISEHMMSTATEBASE_H_

#include "core/Dictionary.hpp"
#include "hmm/DpMatrixBase.hpp"

#include <map>
#include <algorithm>
#include <cmath>


namespace EBC
{

class PairwiseHmmStateBase
{
protected:

	double transFromMatch;
	double transFromInsert;
	double transFromDelete;

	unsigned int rows;
	unsigned int cols;

	DpMatrixBase* dpMatrix;



public:
	Definitions::StateId stateId;

	virtual void initializeData(double lnl, bool backwards=false)=0;

	inline double getValueAt(unsigned int row, unsigned int column)
	{
		return this->dpMatrix->valueAt(row,column);
	}

	inline void setValueAt(unsigned int row, unsigned int column, double data)
	{
		this->dpMatrix->setValue(row,column,data);
	}

	virtual ~PairwiseHmmStateBase()
	{
		delete this->dpMatrix;
	}

	void setSourceMatrixPtr(unsigned int row, unsigned int column, PairwiseHmmStateBase* ptr)
	{
		dpMatrix->setSrc(row,column, ptr->getDpMatrix());
	}

	DpMatrixBase* getDpMatrix()
	{
		return this->dpMatrix;
	}

	virtual void setDirection(unsigned int, unsigned int)=0;

	//void addTransitionProbabilityFrom(PairHmmStateBase* state, double value);

	inline void setTransitionProbabilityFromMatch(double value)
	{
		transFromMatch = value;
	}

	inline void setTransitionProbabilityFromInsert(double value)
	{
		transFromInsert = value;
	}

	inline void setTransitionProbabilityFromDelete(double value)
	{
		transFromDelete = value;
	}

	inline double getTransitionProbabilityFromMatch()
	{
		return transFromMatch;
	}

	inline double getTransitionProbabilityFromInsert()
	{
		return transFromInsert;
	}
	inline double getTransitionProbabilityFromDelete()
	{
		return transFromDelete;
	}

	void traceback(string& seqA, string& seqB, std::pair<string,string>* alignment)
	{
		this->dpMatrix->traceback(seqA,seqB,alignment);
	}

	void tracebackRaw(vector<SequenceElement> s1, vector<SequenceElement> s2, Dictionary* dict, vector<std::pair<unsigned int, unsigned int> >& al)
	{
		this->dpMatrix->tracebackRaw(s1,s2,dict,al);
	}

	unsigned int getCols() const
	{
		return cols;
	}

	unsigned int getRows() const
	{
		return rows;
	}
};

} /* namespace EBC */
#endif /* PAIRWISEHMMSTATEBASE_H_ */
