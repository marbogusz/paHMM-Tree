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


#ifndef EVOLUTIONARYPAIRHMM_HPP_
#define EVOLUTIONARYPAIRHMM_HPP_


#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/PMatrixDouble.hpp"
#include "core/TransitionProbabilities.hpp"
#include "core/Sequences.hpp"
#include "core/Maths.hpp"
#include "core/SequenceElement.hpp"

#include "hmm/PairwiseHmmStateBase.hpp"
#include "hmm/PairwiseHmmInsertState.hpp"
#include "hmm/PairwiseHmmDeleteState.hpp"
#include "hmm/PairwiseHmmMatchState.hpp"
#include "hmm/DpMatrixLoMem.hpp"

#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/IndelModel.hpp"
#include "models/SubstitutionModelBase.hpp"

#include "heuristics/Band.hpp"

namespace EBC
{

//Thats a 3 state pair-HMM with an insertion,deletion and match state
//The states are connected by the silent (blank) state
class EvolutionaryPairHMM
{
protected:
	SubstitutionModelBase* substModel;

	PMatrixDouble* ptmatrix;

	IndelModel* indelModel;
	Sequences* inputSequences;
	TransitionProbabilities* tpb;
	Maths* maths;

	double * mlParameters;
	unsigned int totalParameters;
	unsigned int substParameters;
	unsigned int indelParameters;

	unsigned int xSize, ySize;

	unsigned int gammaRateCategories;

	bool bandingEnabled;

	Band* band;

	bool equilibriumFreqs;

	double initTransM;
	double initTransX;
	double initTransY;

	vector<SequenceElement*>* seq1;
	vector<SequenceElement*>* seq2;
	//vector<SequenceElement>::iterator itS1, itS2;
	
	//cumulative likelihood for all 3 matrices
	double totalLikelihood;

	//state transition matrix
	double md[Definitions::stateCount][Definitions::stateCount];

	//gap probs;
	double e,g;

	//end transition probs
	double xi;

	//state equilibruim frequencies
	double piM, piI, piD;

	//the following assumes a fix HMM structure
	virtual void setTransitionProbabilities();

	virtual void calculateModels();

	virtual void initializeStates(Definitions::DpMatrixType mt);

	void getStateEquilibriums();


public:

	//Match state
	PairwiseHmmStateBase* M;
	//Insert state
	PairwiseHmmStateBase* X;
	//Delete state
	PairwiseHmmStateBase* Y;

	EvolutionaryPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl,
			IndelModel* imdl, Definitions::DpMatrixType, Band* bandObj, bool useEquilibriumFreqs);

	inline void setBand(Band* bnd)
	{
		this->band = bnd;
	}

	virtual ~EvolutionaryPairHMM();

	virtual double runAlgorithm()=0;

	void summarize();

	void setDivergenceTimeAndCalculateModels(double time);

	unsigned int getIndelParameterCount()
	{
		return this->indelModel->getParamsNumber();
	}

	unsigned int getSubstitutionParameterCount()
	{
		return this->substModel->getParamsNumber();
	}

	unsigned int getTotalParameters() const
	{
		return totalParameters;
	}

	PairwiseHmmStateBase* getM()
	{
		return M;
	}

	PairwiseHmmStateBase* getX()
	{
		return X;
	}

	PairwiseHmmStateBase* getY()
	{
		return Y;
	}

	double getAlignmentLikelihood(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2);

	double getAlignmentLikelihood(vector<unsigned char>* s1, vector<unsigned char>* s2, Dictionary* dict);

	double getTotalLikelihood() const {
		return totalLikelihood;
	}

	void setTotalLikelihood(double totalLikelihood) {
		this->totalLikelihood = totalLikelihood;
	}
};

} /* namespace EBC */
#endif /* EVOLUTIONARYPAIRHMM_HPP_ */
