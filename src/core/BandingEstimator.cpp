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


#include "core/BandingEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "models/NegativeBinomialGapModel.hpp"
#include "hmm/DpMatrixFull.hpp"

namespace EBC
{

BandingEstimator::ProgressBar::ProgressBar(unsigned int width) : bw(width) {
}

void BandingEstimator::ProgressBar::tick() {
	curr++;

	if ( (curr != n) && (curr % (n/100+1) != 0) ) return;

    float ratio  =  curr/(float)n;
    int   c      =  ratio * bw;

    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<bw; x++) cout << " ";
    cout << "]\r" << flush;
}

void BandingEstimator::ProgressBar::setIter(unsigned int iters){
	n = iters;
	this->curr =0;
	//cout << setw(3) <<  "0% [";
}

void BandingEstimator::ProgressBar::done(){
	cout << endl;
}

BandingEstimator::BandingEstimator(Definitions::AlgorithmType at, Sequences* inputSeqs, Definitions::ModelType model ,std::vector<double> indel_params,
		std::vector<double> subst_params, Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, GuideTree* g) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), pairCount(inputSequences->getPairCount()),
				/*hmms(pairCount), bands(pairCount),*/ divergenceTimes(pairCount), algorithm(at), gt(g)
{
	//Banding estimator means banding enabled!

	DEBUG("Starting Banding Estimator");
	maths = new Maths();
	dict = inputSequences->getDictionary();

	//Helper models
	if (model == Definitions::ModelType::GTR)
	{
		substModel = new GTRModel(dict, maths,gammaRateCategories);
	}
	else if (model == Definitions::ModelType::HKY85)
	{
		substModel = new HKY85Model(dict, maths,gammaRateCategories);
	}
	else if (model >= Definitions::ModelType::LG)
	{
		switch(model){
			case Definitions::ModelType::LG :
				substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
				DEBUG("Using LG model");
			break;
			case Definitions::ModelType::JTT :
				substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaJttModel);
				DEBUG("Using JTT model");
			break;
			case Definitions::ModelType::WAG :
				substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaWagModel);
				DEBUG("Using WAG model");
			break;
					}
	}

	indelModel = new NegativeBinomialGapModel();

	estimateSubstitutionParams = false;
	estimateIndelParams = false;
	estimateAlpha = false;
	/*
	estimateSubstitutionParams = subst_params.size() != substModel->getParamsNumber();
	estimateIndelParams = indel_params.size() == 0;

	DEBUG("Pairwise banding model estimator starting");
	DEBUG("Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams);
	DEBUG("Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha);

	FileLogger::DebugLogger() << "Estimate substitution parameters set to : " << estimateSubstitutionParams << " Estimate indel parameters set to : " << estimateIndelParams << "\n";
	FileLogger::DebugLogger() << "Estimate alpha set to : " << estimateAlpha << " , rate categories " << gammaRateCategories << " , alpha : " << alpha << "\n";
	 */
	//pairwise comparison mode
	modelParams = new OptimizedModelParameters(substModel, indelModel,2, 1, estimateSubstitutionParams,
			estimateIndelParams, estimateAlpha, true, maths);

	modelParams->boundDivergenceBasedOnLambda(indel_params[0]);

	if(!estimateIndelParams)
		modelParams->setUserIndelParams(indel_params);
	if(!estimateSubstitutionParams)
		modelParams->setUserSubstParams(subst_params);
	modelParams->setAlpha(alpha);

	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	if (estimateSubstitutionParams == false)
	{
		//set parameters and calculate the model
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters(modelParams->getSubstParameters());
		substModel->calculateModel();
	}

	if (estimateIndelParams == false)
	{
			//set parameters and calculate the model
		indelModel->setParameters(modelParams->getIndelParameters());
	}


	EvolutionaryPairHMM *hmm;

	numopt = new BrentOptimizer(modelParams, NULL);

}

BandingEstimator::~BandingEstimator()
{
	//for(auto hmm : hmms)
	//	delete hmm;
	delete numopt;
	delete modelParams;
    delete maths;
    //for (auto bnd : bands)
    //	delete bnd;
    delete indelModel;
    delete substModel;
}

void BandingEstimator::optimizePairByPair()
{
	EvolutionaryPairHMM* hmm;
	Band* band;
	DistanceMatrix* dm = gt->getDistanceMatrix();
	PairHmmCalculationWrapper* wrapper = new PairHmmCalculationWrapper();
	double result;

	ProgressBar pb(80);
	pb.setIter(pairCount);

	for(unsigned int i =0; i< pairCount; i++)
	{
		DEBUG("Optimizing distance for pair #" << i);
		std::pair<unsigned int, unsigned int> idxs = inputSequences->getPairOfSequenceIndices(i);
		INFO("Running pairwise calculator for sequence id " << idxs.first << " and " << idxs.second
				<< " ,number " << i+1 <<" out of " << pairCount << " pairs" );
		BandCalculator* bc = new BandCalculator(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
				substModel, indelModel, gt->getDistanceMatrix()->getDistance(idxs.first,idxs.second));
		band = bc->getBand();
		if (algorithm == Definitions::AlgorithmType::Viterbi)
		{
			DEBUG("Creating Viterbi algorithm to optimize the pairwise divergence time...");
			hmm = new ViterbiPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, band);
		}
		else if (algorithm == Definitions::AlgorithmType::Forward)
		{
			DEBUG("Creating forward algorithm to optimize the pairwise divergence time...");
			hmm = new ForwardPairHMM(inputSequences->getSequencesAt(idxs.first), inputSequences->getSequencesAt(idxs.second),
					substModel, indelModel, Definitions::DpMatrixType::Full, band);
		}

		//hmm->setDivergenceTimeAndCalculateModels(modelParams->getDivergenceTime(0)); //zero as there's only one pair!

		//LikelihoodSurfacePlotter lsp;
		//lsp.setTargetHMM(hmm);
		//lsp.getLikelihoodSurface();


		wrapper->setTargetHMM(hmm);
		DUMP("Set model parameter in the hmm...");
		wrapper->setModelParameters(modelParams);
		modelParams->setUserDivergenceParams({bc->getClosestDistance()});
		numopt->setTarget(wrapper);
		numopt->setAccuracy(bc->getBrentAccuracy());
		numopt->setBounds(bc->getLeftBound(), bc->getRightBound() < 0 ? modelParams->divergenceBound : bc->getRightBound());


		result = numopt->optimize() * -1.0;
		DEBUG("Likelihood after pairwise optimization: " << result);
		if (result <= (Definitions::minMatrixLikelihood /2.0))
		{
			DEBUG("Optimization failed for pair #" << i << " Zero probability FWD");
			band->output();
			dynamic_cast<DpMatrixFull*>(hmm->M->getDpMatrix())->outputValuesWithBands(band->getMatchBand() ,band->getInsertBand(),band->getDeleteBand(),'|', '-');
			dynamic_cast<DpMatrixFull*>(hmm->X->getDpMatrix())->outputValuesWithBands(band->getInsertBand(),band->getMatchBand() ,band->getDeleteBand(),'\\', '-');
			dynamic_cast<DpMatrixFull*>(hmm->Y->getDpMatrix())->outputValuesWithBands(band->getDeleteBand(),band->getMatchBand() ,band->getInsertBand(),'\\', '|');
		}
		this->divergenceTimes[i] = modelParams->getDivergenceTime(0);

		pb.tick();

		delete band;
		delete bc;
		delete hmm;
	}

	pb.done();

	INFO("Optimized divergence times:");
	INFO(this->divergenceTimes);
}


double BandingEstimator::runIteration()
{
	double result = 0;
	double tmp;
	return result;
}

void BandingEstimator::outputDistanceMatrix(stringstream& ss)
{
	unsigned int count, pairCount;
	count = this->inputSequences->getSequenceCount();
	pairCount = this->inputSequences->getPairCount();

	ss << "\t" << this->inputSequences->getSequenceCount() << endl;

	for (unsigned int i = 0; i< count; i++)
	{
		ss << "S" << i << " ";
		for(unsigned int j=0; j< count; j++)
		{
			ss << this->modelParams->getDistanceBetween(i,j) << " ";
		}
		ss << endl;
	}
}

} /* namespace EBC */


