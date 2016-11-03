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


#ifndef MODELESTIMATOR_HPP_
#define MODELESTIMATOR_HPP_


#include "core/OptimizedModelParameters.hpp"
#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/PMatrixTriple.hpp"

#include "models/SubstitutionModelBase.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/GTRModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"

#include "heuristics/ModelEstimator.hpp"
#include "heuristics/GuideTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"
#include "heuristics/StateTransitionEstimator.hpp"
#include "heuristics/SubstitutionModelEstimator.hpp"

#include "hmm/ViterbiPairHMM.hpp"
#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"
#include "hmm/DpMatrixFull.hpp"


#include <sstream>
#include <vector>
#include <array>

using namespace std;


namespace EBC
{

class ModelEstimator
{
protected:


	Dictionary* dict;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	Sequences* inputSequences;
	Maths* maths;

	GuideTree* gtree;
	TripletSamplingTree tst;

	StateTransitionEstimator* ste;
	SubstitutionModelEstimator* sme;

	TripletAligner* tal;
	ViterbiPairHMM* vphmm;

	double bestFwdAlpha;
	double bestFwdTm;

	bool estAlpha;
	bool estIndel;
	bool estSubst;

	vector<array<vector<unsigned char>*, 3> > tripleAlignments;
	vector<array<vector<unsigned char>*, 4> > pairAlignments;
	vector<array<vector<double>*, 2> > pairwisePosteriors;
	vector<array<unsigned int, 3> > tripletIdxs;
	vector<array<double, 3> > tripletDistances;
	vector<array<ForwardPairHMM*,2> > fwdHMMs;

	vector<double> substitutionParameters;
	vector<double> indelParameters;
	double alpha;

	unsigned int gammaRateCategories;

	unsigned int tripletIdxsSize;

	double userAlpha;

	bool estimateAlpha;

	void calculateInitialHMMs(Definitions::ModelType model);

	void estimateParameters();

	Definitions::ModelType model;

	vector<double> getInitialModelParameters();

	void doSME();

public:
	ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model,
			Definitions::OptimizationType ot,
			unsigned int rateCategories, double alpha, bool estimateAlpha);

	virtual ~ModelEstimator();



	vector<double> getSubstitutionParameters();
	vector<double> getIndelParameters();
	double getAlpha();

	void recalculateHMMs();

	GuideTree* getGuideTree()
	{
		return gtree;
	}
};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
