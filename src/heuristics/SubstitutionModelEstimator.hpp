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


#ifndef SUBSTMODELESTIMATOR_HPP_
#define SUBSTMODELESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrixTriple.hpp"
#include "core/IOptimizable.hpp"
#include "core/Optimizer.hpp"
#include "core/Sequences.hpp"

#include "models/SubstitutionModelBase.hpp"

#include <vector>
#include <array>

using namespace std;


namespace EBC
{

class SubstitutionModelEstimator : public IOptimizable
{

protected:

	Optimizer* bfgs;

	Dictionary* dict;

	Sequences* inputSequences;

	Maths* maths;

	SubstitutionModelBase* substModel;

	vector<array<PMatrixTriple* ,3> > ptMatrices;

	vector<map<array<unsigned char, 3>, unsigned int> > patterns;

	vector<double> distances;

	unsigned int gammaRateCategories;

	OptimizedModelParameters* modelParams;

	double alpha;

	bool estimateSubstitutionParams;
	bool estimateAlpha;

	unsigned int currentTriplet;

public:
	SubstitutionModelEstimator(Sequences* inputSeqs, SubstitutionModelBase* model,
			Definitions::OptimizationType ot,unsigned int rateCategories, double alpha,
			bool estimateAlpha, unsigned int matCount);

	virtual ~SubstitutionModelEstimator();

	void addTriplet(array<vector<unsigned char>*, 3> tripleAlignment, unsigned int tiplet,
			double d1 = 0.5, double d2 = 0.5, double d3 = 0.5);

	double runIteration();

	void optimize();

	void clean();

	double getTripletDivergence(unsigned int triplet, unsigned int branch)
	{
		//no bound checks - beware
		return modelParams->getDivergenceTime(((triplet*3)+branch));
	}

	OptimizedModelParameters* getModelParams()
	{
		return modelParams;
	}

	SubstitutionModelBase* getSubstModel()
	{
		return substModel;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
