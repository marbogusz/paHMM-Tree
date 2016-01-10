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



#ifndef PMATRIX_HPP_
#define PMATRIX_HPP_

#include "models/SubstitutionModelBase.hpp"
#include "core/HmmException.hpp"
#include <vector>
#include <array>

using namespace std;

namespace EBC
{

class PMatrix
{
protected:

	SubstitutionModelBase* model;

	unsigned int matrixSize;

	double time;

	unsigned int rateCategories;

	vector<double*> ptMatrices;

	unsigned int matrixFullSize;

public:
	PMatrix(SubstitutionModelBase* m);
	virtual ~PMatrix();

	void setTime(double t);

	virtual void calculate()=0;

	inline double getEquilibriumFreq(unsigned int xi)
	{
		return model->getEquilibriumFrequencies(xi);
	}

	inline double getLogEquilibriumFreq(unsigned int xi)
	{
		return model->getLogEquilibriumFrequencies(xi);
	}

	void summarize();

};

} /* namespace EBC */

#endif /* PMATRIX_HPP_ */
