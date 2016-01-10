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


#ifndef TRANSITIONPROBABILITIES_HPP_
#define TRANSITIONPROBABILITIES_HPP_

#include "models/IndelModel.hpp"

namespace EBC
{

class TransitionProbabilities
{
protected:

	double gapOpening;
	double gapExtension;

	double time;

	IndelModel* indelModel;

public:
	TransitionProbabilities(IndelModel* im);

	virtual ~TransitionProbabilities();

	void calculate();

	double getGapExtension() const
	{
		return gapExtension;
	}

	double getGapOpening() const
	{
		return gapOpening;
	}

	void setTime(double time)
	{
		this->time = time;
	}
};

} /* namespace EBC */

#endif /* TRANSITIONPROBABILITIES_HPP_ */
