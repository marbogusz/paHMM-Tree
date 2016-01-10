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


#include <core/Optimizer.hpp>

namespace EBC
{

Optimizer::Optimizer(OptimizedModelParameters* mp, IOptimizable* opt, Definitions::OptimizationType ot, double accuracy) :
		accuracy(accuracy), optimizationType(ot), omp(mp), target(opt)
{


	paramsCount = omp->optParamCount();
	this->initParams.set_size(paramsCount);
	this->lowerBounds.set_size(paramsCount);
	this->upperBounds.set_size(paramsCount);

	//omp->toDlibVector(initParams,lowerBounds,upperBounds);

	DEBUG("Numeric optimizer with " << paramsCount << " parameter(s) created");
}

void Optimizer::setTarget(IOptimizable* opt)
{
	target = opt;
}

Optimizer::~Optimizer()
{
}

double Optimizer::objectiveFunction(const column_vector& bfgsParameters)
{
	omp->fromDlibVector(bfgsParameters);
	return target->runIteration();
}


const column_vector Optimizer::objectiveFunctionDerivative(const column_vector& bfgsParameters)
{
	column_vector results(this->paramsCount);
	return results;
}


double Optimizer::optimize()
{
	omp->toDlibVector(initParams,lowerBounds,upperBounds);

	using std::placeholders::_1;
	std::function<double(const column_vector&)> f_objective= std::bind( &Optimizer::objectiveFunction, this, _1 );
	double likelihood;

	switch(optimizationType)
	{
		case Definitions::OptimizationType::BFGS:
		{
			likelihood = dlib::find_min_box_constrained(dlib::bfgs_search_strategy(),
					dlib::objective_delta_stop_strategy(accuracy),  //changed the delta drastically
					f_objective,
					derivative(f_objective),
					initParams,
					lowerBounds,
					upperBounds);
			break;
		}
		case Definitions::OptimizationType::BOBYQA:
		{
			likelihood = dlib::find_min_bobyqa(f_objective, initParams, paramsCount+4,
					lowerBounds,upperBounds, 0.05, 1e-7, 20000 );
			break;
		}
	}
	omp->fromDlibVector(initParams);
	return likelihood;
	//omp->outputParameters();
	//cout  << likelihood << "\n";

}


} /* namespace EBC */
