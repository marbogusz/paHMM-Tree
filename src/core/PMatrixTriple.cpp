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


#include <core/PMatrixTriple.hpp>

namespace EBC
{

PMatrixTriple::PMatrixTriple(SubstitutionModelBase* m) : PMatrix(m)
{

}

PMatrixTriple::~PMatrixTriple()
{
}

void PMatrixTriple::calculate()
{
	if (time != 0)
	{
		for(unsigned int i = 0; i< rateCategories; i++)
		{
			if (ptMatrices[i] != NULL)
				delete [] ptMatrices[i];
			ptMatrices[i] = this->model->calculatePt(time, i);

		}
	}
	else
		throw HmmException("PMatrixTriple : attempting to calculate p(t) with t set to 0");
}


double PMatrixTriple::getTripleSitePattern(unsigned int root,
		const array<unsigned char, 3>& nodes, PMatrixTriple* pm2, PMatrixTriple* pm3)
{
	PMatrixTriple* pm1 = this;
	double prob = 0;
	for(int rt = 0; rt< rateCategories; rt++)
	{
		prob += getEquilibriumFreq(root) * pm1->getTransitionProb(root,nodes[0], rt) *
				pm2->getTransitionProb(root,nodes[1], rt) *
				pm3->getTransitionProb(root,nodes[2], rt) * model->gammaFrequencies[rt];
	}

	return prob;
}

double PMatrixTriple::getTransitionProb(unsigned int xi, unsigned int yi, unsigned int rateCat)
{
	if (yi >= matrixSize)
		return 1;
	return (ptMatrices[rateCat])[xi*matrixSize+yi];
}

void PMatrixTriple::summarize()
{
	cout << "P(t) matrix summary :" << endl;
	cout << "Divergence time : " << time << endl;

}

} /* namespace EBC */


