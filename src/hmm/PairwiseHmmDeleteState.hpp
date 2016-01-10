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


#ifndef PAIRWISEHMMDELETIONSTATE_HPP_
#define PAIRWISEHMMDELETIONSTATE_HPP_

#include "hmm/PairwiseHmmStateBase.hpp"

namespace EBC
{

class PairwiseHmmDeleteState: public EBC::PairwiseHmmStateBase
{
public:
	void initializeData(double lnl, bool backwards=false);

	PairwiseHmmDeleteState(unsigned int x, unsigned int y);

	PairwiseHmmDeleteState(DpMatrixBase*);

	void setDirection(unsigned int i, unsigned int j);

	virtual ~PairwiseHmmDeleteState();
};

} /* namespace EBC */
#endif /* PAIRHMMDELETIONSTATE_HPP_ */
