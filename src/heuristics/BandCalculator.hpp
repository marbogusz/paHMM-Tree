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

#ifndef HEURISTICS_BANDCALCULATOR_HPP_
#define HEURISTICS_BANDCALCULATOR_HPP_

#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/SubstitutionModelBase.hpp"

#include "core/PMatrixDouble.hpp"
#include "core/TransitionProbabilities.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "models/IndelModel.hpp"
#include "core/Sequences.hpp"

#include "hmm/ForwardPairHMM.hpp"
#include "hmm/BackwardPairHMM.hpp"

#include "heuristics/Band.hpp"

#include<vector>

namespace EBC
{

class BandCalculator
{
protected:

	vector<ForwardPairHMM*> fwd;
	BackwardPairHMM* bwd;

	vector<SequenceElement*>* seq1;
	vector<SequenceElement*>* seq2;

	SubstitutionModelBase* substModel;
	IndelModel* indelModel;

	double time;


	PMatrixDouble* ptMatrix;
	TransitionProbabilities* trProbs;

	Band* band;

	double posteriorLikelihoodLimit;
	double posteriorLikelihoodDelta;

	double bestTime;
	double accuracy;

	double leftBound;
	double rightBound;

	void processPosteriorProbabilities(BackwardPairHMM* hmm, Band* band);

public:
	BandCalculator(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* sm, IndelModel* im, double divergenceTime);
	virtual ~BandCalculator();

	inline Band* getBand()
	{
		return this->band;
	}

	double getClosestDistance();

	double getBrentAccuracy();

	double getLeftBound() {
		return leftBound;
	}

	double getRightBound(){
		return rightBound;
	}
};

} /* namespace EBC */

#endif /* HEURISTICS_BANDCALCULATOR_HPP_ */
