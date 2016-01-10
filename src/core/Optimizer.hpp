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


#ifndef OPTIMIZER_HPP_
#define OPTIMIZER_HPP_

#include "core/OptimizedModelParameters.hpp"
#include "core/IOptimizable.hpp"

namespace EBC
{

class Optimizer
	{
	protected:
		column_vector initParams;
		column_vector lowerBounds;
		column_vector upperBounds;

		unsigned int paramsCount;

		OptimizedModelParameters* omp;
		IOptimizable* target;
		double accuracy;

		Definitions::OptimizationType optimizationType;

	public:
		Optimizer(OptimizedModelParameters* mp, IOptimizable* opt, Definitions::OptimizationType ot, double accuracy=Definitions::accuracyBFGS);
		virtual ~Optimizer();

		double optimize();

		void setTarget(IOptimizable* opt);

		double objectiveFunction(const column_vector& m);

		const column_vector objectiveFunctionDerivative(const column_vector& m);
	};

} /* namespace EBC */

#endif /* OPTIMIZER_HPP_ */
