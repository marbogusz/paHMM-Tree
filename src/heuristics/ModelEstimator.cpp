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


#include "core/Dictionary.hpp"
#include "heuristics/ModelEstimator.hpp"
#include <chrono>
#include <array>
#include <map>

using namespace std;

namespace EBC
{

ModelEstimator::ModelEstimator(Sequences* inputSeqs, Definitions::ModelType model ,
		Definitions::OptimizationType ot, unsigned int rateCategories, double alpha, bool estimateAlpha) :
				inputSequences(inputSeqs), gammaRateCategories(rateCategories), model(model),
				gtree(new GuideTree(inputSeqs)), tst(*gtree), userAlpha(alpha), estAlpha(estimateAlpha), estIndel(true), estSubst(true)
{

	DEBUG("About to sample some triplets");
	DEBUG("Sampling triplets of sequences for gamma shape parameter estimation");
	
	maths = new Maths();
	dict = inputSequences->getDictionary();
	tal = new TripletAligner (inputSequences, gtree->getDistanceMatrix(), -1.0);

	//Up to
	tripletIdxs = tst.sampleFromTree();

	tripletIdxsSize = tripletIdxs.size();

	if (tripletIdxsSize == 0)
		throw HmmException("No triplets selected for model estimation");

	tripleAlignments.resize(tripletIdxsSize);
	pairAlignments.resize(tripletIdxsSize);
	pairwisePosteriors.resize(tripletIdxsSize);
	fwdHMMs.resize(tripletIdxsSize);
	tripletDistances.resize(tripletIdxsSize);

	//FIXME - release the memory!!!! - delete pair objects and vectors (arrays) of SeqEls

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

	this->calculateInitialHMMs(model);

	sme = new SubstitutionModelEstimator(inputSequences, substModel ,ot, rateCategories, alpha, estimateAlpha, tripletIdxsSize);
	ste = new StateTransitionEstimator(indelModel, ot, 2*tripletIdxsSize, dict->getGapID(),false);

	estimateParameters();
	//recalculateHMMs();
	//estimateParameters();

	end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;

    INFO("Model Estimator elapsed time: " << elapsed_seconds.count() << " seconds");
    //cerr <<  "|||||||||||| Model Estimator elapsed time: " << elapsed_seconds.count() << "s ||||||||||||||\n";

	//substModel->summarize();

}

vector<double> ModelEstimator::getInitialModelParameters()
{
	double initKappa = 2.0;
	vector<double> params;
	if (model == Definitions::ModelType::HKY85){
		params =  {initKappa};
	}
	else if (model == Definitions::ModelType::GTR){
		params = {{1.0,1.0/initKappa,1.0/initKappa,1.0/initKappa,1.0/initKappa}};
	}
	return params;

}

void ModelEstimator::recalculateHMMs()
{//Fwd + bwd + MPD

	ForwardPairHMM *f1, *f2;
	double tb1, tb2, tb3;

	substModel->calculateModel();

	for (int i = 0; i < tripletIdxsSize; i++){

		delete tripleAlignments[i][0];
		delete tripleAlignments[i][1];
		delete tripleAlignments[i][2];

		delete pairAlignments[i][0];
		delete pairAlignments[i][1];
		delete pairAlignments[i][2];
		delete pairAlignments[i][3];

		delete pairwisePosteriors[i][0];
		delete pairwisePosteriors[i][1];

		f1 = fwdHMMs[i][0];
		f2 = fwdHMMs[i][1];

		tb1 = sme->getTripletDivergence(i,0);
		tb2 = sme->getTripletDivergence(i,1);
		tb3 = sme->getTripletDivergence(i,2);

		f1->setDivergenceTimeAndCalculateModels(tb1+tb2);
		f2->setDivergenceTimeAndCalculateModels(tb2+tb3);

		f1->runAlgorithm();
		f2->runAlgorithm();

		BackwardPairHMM b1(inputSequences->getSequencesAt(tripletIdxs[i][0]),inputSequences->getSequencesAt(tripletIdxs[i][1]),
				substModel, indelModel, Definitions::DpMatrixType::Full, nullptr);
		BackwardPairHMM b2(inputSequences->getSequencesAt(tripletIdxs[i][1]),inputSequences->getSequencesAt(tripletIdxs[i][2]),
				substModel, indelModel, Definitions::DpMatrixType::Full, nullptr);

		b1.setDivergenceTimeAndCalculateModels(tb1+tb2);
		b2.setDivergenceTimeAndCalculateModels(tb2+tb3);

		b1.runAlgorithm();
		b2.runAlgorithm();

		b1.calculatePosteriors(f1);
		b2.calculatePosteriors(f2);

		b1.calculateMaximumPosteriorMatrix();
		b2.calculateMaximumPosteriorMatrix();

		auto mp1 = b1.getMPAlignment();
		auto mp2 = b2.getMPAlignment();

		DUMP("Pair 1 MPD alignment recalc");
		DUMP(mp1.first);
		DUMP(mp1.second);
		DUMP("Pair 2 MPD alignment recalc");
		DUMP(mp2.first);
		DUMP(mp2.second);


		//delete f1;
		//delete f2;

		//store pairs, align triplets
		pair<vector<double>*, pair<vector<unsigned char>*, vector<unsigned char>*> > alP1 = b1.getMPDWithPosteriors();
		pair<vector<double>*, pair<vector<unsigned char>*, vector<unsigned char>*> > alP2 = b2.getMPDWithPosteriors();

		pairAlignments[i][0] = alP1.second.first;
		pairAlignments[i][1] = alP1.second.second;
		pairAlignments[i][2] = alP2.second.first;
		pairAlignments[i][3] = alP2.second.second;

		pairwisePosteriors[i][0] = alP1.first;
		pairwisePosteriors[i][1] = alP2.first;

		tripleAlignments[i] = tal->alignPosteriors(alP1.second, alP2.second, alP1.first, alP2.first);

	}

	sme->clean();
	ste->clean();
}

void ModelEstimator::doSME()
{
	double d1,d2,d3;
	DUMP("Model Estimator estimate parameters iteration");

	for(unsigned int trp = 0; trp < tripletIdxsSize; trp++)
	{
		d1 = (tripletDistances[trp][0] + tripletDistances[trp][2] - tripletDistances[trp][1])/2.0;
		d2 = tripletDistances[trp][0] - d1;
		d3 = tripletDistances[trp][2] - d1;

		sme->addTriplet(tripleAlignments[trp], trp, d1, d2, d3);
	}
	sme->optimize();
}

void ModelEstimator::estimateParameters()
{


	double tb1,tb2,tb3,steT1, steT2;
	int cleaned = 0;
	bool zeroBL = false;

	//Estimate subst. model
	this->doSME();

	//Check for 0 branch lengths; delete triplets in case of 0 b.l.
	//Can't remove all though!
	for(int trp = tripletIdxsSize-1; trp >= 0; trp--)
	{
		tb1 = sme->getTripletDivergence(trp,0);
		tb2 = sme->getTripletDivergence(trp,1);
		tb3 = sme->getTripletDivergence(trp,2);
		if(tb1 <= 1.1*Definitions::almostZero || tb2 <= 1.1*Definitions::almostZero || tb3 <= 1.1 * Definitions::almostZero){
			if(tripletIdxsSize == 1){
				//can't purge all the triplets!
				continue;
			}
			WARN("Triplet estimator zero length branch in triplet " << trp << " Skipping; branch lengths " << tb1 << " " << tb2 << " " << tb3);

			delete tripleAlignments[trp][0];
			delete tripleAlignments[trp][1];
			delete tripleAlignments[trp][2];

			delete pairAlignments[trp][0];
			delete pairAlignments[trp][1];
			delete pairAlignments[trp][2];
			delete pairAlignments[trp][3];

			delete pairwisePosteriors[trp][0];
			delete pairwisePosteriors[trp][1];

			delete fwdHMMs[trp][0];
			delete fwdHMMs[trp][1];

			tripleAlignments.erase(tripleAlignments.begin() + trp);
			pairAlignments.erase(pairAlignments.begin() + trp);
			pairwisePosteriors.erase(pairwisePosteriors.begin() + trp);
			fwdHMMs.erase(fwdHMMs.begin() + trp);
			tripletIdxs.erase(tripletIdxs.begin() + trp);
			tripletDistances.erase(tripletDistances.begin() + trp);
			tripletIdxsSize--;
			cleaned++;
			zeroBL = true;
		}
	}
	if(zeroBL){
		sme->clean(cleaned);
		ste->clean(cleaned);
		this->doSME();
	}

	substitutionParameters = sme->getModelParams()->getSubstParameters();
	alpha = sme->getModelParams()->getAlpha();

	for(unsigned int trp = 0; trp < tripletIdxsSize; trp++)
	{
		tb1 = sme->getTripletDivergence(trp,0);
		tb2 = sme->getTripletDivergence(trp,1);
		tb3 = sme->getTripletDivergence(trp,2);

		steT1 = tb1 + tb2;
		steT2 = tb3 + tb2;

		if(tb1 <= 1.1*Definitions::almostZero || tb2 <= 1.1*Definitions::almostZero || tb3 <= 1.1 * Definitions::almostZero){
			INFO("Triplet estimator zero length branch " << tb1 << " " << tb2 << " " << tb3);
			WARN("Problem while estimating substitution model parameters: zero length branch in the estimator triplet tree. Model parameters might be inaccurate");
			cout << "Problem while estimating substitution model parameters: zero length branch in the estimator triplet tree.\n Model parameters might be inaccurate\n";
			steT1 = tripletDistances[trp][0] * bestFwdTm;
			steT2 = tripletDistances[trp][1] * bestFwdTm;
			alpha = bestFwdAlpha;
			//substitutionParameters = getInitialModelParameters();
		}



		ste->addTime(steT1,trp,0);
		ste->addTime(steT2,trp,1);

		ste->addPair(pairAlignments[trp][0],pairAlignments[trp][1],trp,0);
		ste->addPair(pairAlignments[trp][2],pairAlignments[trp][3],trp,1);

	}
	ste->optimize();

	indelParameters = ste->getModelParams()->getIndelParameters();



	//substModel->summarize();
	//indelModel->summarize();
}

void ModelEstimator::calculateInitialHMMs(Definitions::ModelType model)
{
	DEBUG("Estimating Triple Aligments");
	//amino acid mode
	bool aaMode = false;
	double tmpd;

	double initAlpha = 0.75;
	double initKappa = 2.0;
	double initLambda = 0.05;
	double initEpsilon = 0.5;
	//k-mers tend to underestimate the distances;
	double initTimeModifier = 1.5;

	vector<double> timeModifiers = {0.75, 1.0, 1.5};
	vector<double> lambdas = {0.03, 0.07};
	vector<double> alphas;

	array<double,3> lowMultipliers = {{0.7,1,1.3}};
	array<double,3> normalMultipliers = {{0.5,1,1.5}};
	array<double,3> highMultipliers = {{1,2,4}};

	if (!estAlpha)
		alphas = {userAlpha};
	else
		alphas = {0.5,1.0,3.0};

	ForwardPairHMM *f1, *f2;

	if (model == Definitions::ModelType::HKY85){
			DEBUG("Setting HKY85");
			this->substModel = new HKY85Model(dict, maths,gammaRateCategories);
			substModel->setParameters({initKappa});
	}
	else if (model == Definitions::ModelType::GTR){
			DEBUG("Setting GTR");
			this->substModel = new GTRModel(dict, maths,gammaRateCategories);
			substModel->setParameters({1.0,1.0/initKappa,1.0/initKappa,1.0/initKappa,1.0/initKappa});
	}
	else if (model == Definitions::ModelType::LG){
			aaMode = true;
			substModel = new AminoacidSubstitutionModel(dict, maths,gammaRateCategories,Definitions::aaLgModel);
	}


	substModel->setObservedFrequencies(inputSequences->getElementFrequencies());
	//alpha setting will have no effect if we're dealing with 1 rate category
	if (estAlpha)
		substModel->setAlpha(initAlpha);
	else
		substModel->setAlpha(userAlpha);
	//if (!aaMode)
	//	substModel->setParameters({initKappa});
	//substModel->calculateModel();

	indelModel = new NegativeBinomialGapModel();
	indelModel->setParameters({initLambda,initEpsilon});

	vector<pair<Band*, Band*> > bandPairs(tripletIdxsSize);

	vector<array<vector<SequenceElement*>*,3> > seqsA(tripletIdxsSize);

	for (int i = 0; i < tripletIdxsSize; i++)
	{


		seqsA[i][0] = inputSequences->getSequencesAt(tripletIdxs[i][0]);
		DUMP("Triplet " << i << " sequence 1:");
		DUMP(inputSequences->getRawSequenceAt(tripletIdxs[i][0]));
		seqsA[i][1] = inputSequences->getSequencesAt(tripletIdxs[i][1]);
		DUMP("Triplet " << i << " sequence 2:");
		DUMP(inputSequences->getRawSequenceAt(tripletIdxs[i][1]));
		seqsA[i][2] = inputSequences->getSequencesAt(tripletIdxs[i][2]);
		DUMP("Triplet " << i << " sequence 3:");
		DUMP(inputSequences->getRawSequenceAt(tripletIdxs[i][2]));

		unsigned int len1 = seqsA[i][0]->size();
		unsigned int len2 = seqsA[i][1]->size();
		unsigned int len3 = seqsA[i][2]->size();

		//0-1
		tmpd = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][1]);
		DEBUG("Triplet " << i << " guide distance between seq 1 and 2 " << tmpd);
		tripletDistances[i][0] = tmpd;
		//1-2
		tmpd = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][1],tripletIdxs[i][2]);
		DEBUG("Triplet " << i << " guide distance between seq 2 and 3 " << tmpd);
		tripletDistances[i][1] = tmpd;
		//0-2
		tmpd = gtree->getDistanceMatrix()->getDistance(tripletIdxs[i][0],tripletIdxs[i][2]);
		DEBUG("Triplet " << i << " guide distance between seq 1 and 3 " << tmpd);
		tripletDistances[i][2] = tmpd;
		bandPairs[i] = make_pair(new Band(len1,len2,tripletDistances[i][0] < Definitions::kmerHighDivergence ? Definitions::narrowBandFactor : Definitions::initialBandFactor ),
				new Band(len2,len3,tripletDistances[i][1] < Definitions::kmerHighDivergence ? Definitions::narrowBandFactor : Definitions::initialBandFactor ));
		//bandPairs[i] = make_pair(nullptr,nullptr);

		fwdHMMs[i][0] = new ForwardPairHMM(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first,true);
		fwdHMMs[i][1] = new ForwardPairHMM(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].second,true);
	}
	double bestLnl = Definitions::minMatrixLikelihood;
	double currentLnl;
	double bestA, bestL, bestTm;


	for (auto tm : timeModifiers){
		for(auto l : lambdas){
			indelModel->setParameters({l,initEpsilon});
			for(auto a : alphas){
				substModel->setAlpha(a);
				substModel->calculateModel();
				currentLnl = 0;
				for (int i = 0; i < tripletIdxsSize; i++){
					f1 = fwdHMMs[i][0];
					f2 = fwdHMMs[i][1];

					f1->setDivergenceTimeAndCalculateModels(tripletDistances[i][0]*tm);
					f2->setDivergenceTimeAndCalculateModels(tripletDistances[i][1]*tm);

					currentLnl += (f1->runAlgorithm() + f2->runAlgorithm()) * -1.0;

				}
				if (currentLnl > bestLnl){
					bestLnl = currentLnl;
					this->bestFwdAlpha = bestA = a;
					bestL = l;
					this->bestFwdTm = bestTm = tm;

				}
			}
		}
	}
	substModel->setAlpha(bestA);
	if (estAlpha)
		substModel->calculateModel();
	indelModel->setParameters({bestL,initEpsilon});

	DUMP("Best a " << bestA << "\tbest l " << bestL << "\ttimeMult " << bestTm );

	//Fwd + bwd + MPD
	for (int i = 0; i < tripletIdxsSize; i++){
		f1 = fwdHMMs[i][0];
		f2 = fwdHMMs[i][1];

		//remove bands for more accurate estimates
		//f1->setBand(nullptr);
		//f2->setBand(nullptr);

		f1->setDivergenceTimeAndCalculateModels(tripletDistances[i][0]*bestTm);
		f2->setDivergenceTimeAndCalculateModels(tripletDistances[i][1]*bestTm);

		f1->runAlgorithm();
		f2->runAlgorithm();

		BackwardPairHMM b1(seqsA[i][0],seqsA[i][1], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].first);
		BackwardPairHMM b2(seqsA[i][1],seqsA[i][2], substModel, indelModel, Definitions::DpMatrixType::Full, bandPairs[i].second);

		b1.setDivergenceTimeAndCalculateModels(tripletDistances[i][0]*bestTm);
		b2.setDivergenceTimeAndCalculateModels(tripletDistances[i][1]*bestTm);

		b1.runAlgorithm();
		b2.runAlgorithm();

		b1.calculatePosteriors(f1);
		b2.calculatePosteriors(f2);

		b1.calculateMaximumPosteriorMatrix();
		b2.calculateMaximumPosteriorMatrix();

		//auto mp1 = b1.getMPAlignment();
		//auto mp2 = b2.getMPAlignment();

		//DUMP("Pair 1 MPD alignment");
		//DUMP(mp1.first);
		//DUMP(mp1.second);
		//DUMP("Pair 2 MPD alignment");
		//DUMP(mp2.first);
		//DUMP(mp2.second);


		//delete f1;
		//delete f2;

		//store pairs, align triplets
		pair<vector<double>*, pair<vector<unsigned char>*, vector<unsigned char>*> > alP1 = b1.getMPDWithPosteriors();
		pair<vector<double>*, pair<vector<unsigned char>*, vector<unsigned char>*> > alP2 = b2.getMPDWithPosteriors();

		pairAlignments[i][0] = alP1.second.first;
		pairAlignments[i][1] = alP1.second.second;
		pairAlignments[i][2] = alP2.second.first;
		pairAlignments[i][3] = alP2.second.second;

		pairwisePosteriors[i][0] = alP1.first;
		pairwisePosteriors[i][1] = alP2.first;

		tripleAlignments[i] = tal->alignPosteriors(alP1.second, alP2.second, alP1.first, alP2.first);

		//TODO - delete band pairs

	}
}




ModelEstimator::~ModelEstimator()
{

	for (int i =0; i < tripletIdxsSize; i++)
	{
		delete tripleAlignments[i][0];
		delete tripleAlignments[i][1];
		delete tripleAlignments[i][2];

		delete pairAlignments[i][0];
		delete pairAlignments[i][1];
		delete pairAlignments[i][2];
		delete pairAlignments[i][3];

		delete pairwisePosteriors[i][0];
		delete pairwisePosteriors[i][1];

		delete fwdHMMs[i][0];
		delete fwdHMMs[i][1];
	}

    delete maths;
    delete sme;
    delete ste;
    delete gtree;
    delete tal;

}

vector<double> ModelEstimator::getSubstitutionParameters()
{
	//return this->sme->getModelParams()->getSubstParameters();
	return this->substitutionParameters;
}

vector<double> ModelEstimator::getIndelParameters()
{
	//return this->ste->getModelParams()->getIndelParameters();
	return this->indelParameters;
}

double ModelEstimator::getAlpha()
{
	//return this->sme->getModelParams()->getAlpha();
	return this->alpha;
}

} /* namespace EBC */
