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


#ifndef PMATRIXDB_HPP_
#define PMATRIXDB_HPP_

#include "models/SubstitutionModelBase.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrix.hpp"
#include "core/SequenceElement.hpp"
#include <vector>
#include <array>

using namespace std;

namespace EBC
{

//TODO - split into 2 regular pt and pairwise
class PMatrixDouble : public PMatrix
{
protected:

	double* fastPairGammaPt;
	double* fastLogPairGammaPt;

	double ** sitePatterns;

	void calculatePairSitePatterns();

public:
	PMatrixDouble(SubstitutionModelBase* m);
	virtual ~PMatrixDouble();

	void calculate();

	inline double getPairSitePattern(unsigned int xi, unsigned int yi)
	{
		return sitePatterns[xi][yi];
	}

	double getPairTransition(array<unsigned int, 2>& nodes);

	double getPairTransition(unsigned int xi, unsigned int yi);

	double getLogPairTransition(unsigned int xi, unsigned int yi);

	double getLogEquilibriumFreqClass(SequenceElement* se);

	double getLogPairTransitionClass(SequenceElement* se1, SequenceElement* se2);



	void summarize();

};

} /* namespace EBC */

#endif /* PMATRIX_HPP_ */
