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


#ifndef PMATRIXTR_HPP_
#define PMATRIXTR_HPP_

#include "models/SubstitutionModelBase.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrix.hpp"
#include <vector>
#include <array>

using namespace std;

namespace EBC
{

//TODO - split into 2 regular pt and pairwise
class PMatrixTriple : public PMatrix
{

public:
	PMatrixTriple(SubstitutionModelBase* m);
	virtual ~PMatrixTriple();


	void calculate();

	double getTransitionProb(unsigned int xi, unsigned int yi, unsigned int rateCat = 0);

	double getTripleSitePattern(unsigned int root,const array<unsigned char, 3>& nodes, PMatrixTriple* pm2, PMatrixTriple* pm3);

	void summarize();

};

} /* namespace EBC */

#endif /* PMATRIX_HPP_ */
