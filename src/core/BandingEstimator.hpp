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


#ifndef BANDINGESTIMATOR_HPP_
#define BANDINGESTIMATOR_HPP_



#include "core/IOptimizable.hpp"
#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/BandingEstimator.hpp"
#include "core/DistanceMatrix.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/Optimizer.hpp"
#include "core/BrentOptimizer.hpp"
#include "core/PairHmmCalculationWrapper.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/IndelModel.hpp"

#include "heuristics/GuideTree.hpp"
#include "heuristics/BandCalculator.hpp"
#include "heuristics/Band.hpp"

#include "hmm/ForwardPairHMM.hpp"
#include "hmm/ViterbiPairHMM.hpp"

#include <vector>
#include <sstream>

using namespace std;


namespace EBC
{

//Estimate pairwise distances using bands
class BandingEstimator : public IOptimizable
{

private:
	class ProgressBar
	{
	private:
		//bar width
		unsigned int bw;
		//max no iter
		unsigned int n;
		//current iter;
		unsigned int curr;
	public:
		ProgressBar(unsigned int width);
		void tick();
		void setIter(unsigned int);
		void done();
	};

protected:

	BrentOptimizer* numopt;
	Dictionary* dict;
	SubstitutionModelBase* substModel;
	IndelModel* indelModel;
	Sequences* inputSequences;
	Maths* maths;
	GuideTree* gt;

	Definitions::AlgorithmType algorithm;

	unsigned int bandFactor;
	unsigned int bandSpan;
	unsigned int gammaRateCategories;

	bool bandingEnabled;

	bool estimateSubstitutionParams;
	bool estimateIndelParams;
	bool estimateDivergence;
	bool estimateAlpha;

	unsigned int pairCount;

	//vector<EvolutionaryPairHMM*> hmms;
	//delete bands in the destructor
	//vector<Band*> bands;
	vector<double> divergenceTimes;

	OptimizedModelParameters* modelParams;

public:
	BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model,std::vector<double> indel_params,
			std::vector<double> subst_params, Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, GuideTree* gt);

	virtual ~BandingEstimator();

	void runIteration(const column_vector& bfgsParameters);

	double runIteration();

	void outputDistanceMatrix(stringstream&);

	void optimizePairByPair();

	vector<double> getOptimizedTimes()
	{
		return this->divergenceTimes;
	}

	//ModelParameters getMlParameters()
	//{
	//	return this->modelParameters;
	//}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
