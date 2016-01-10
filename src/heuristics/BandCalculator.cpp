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


#include <heuristics/BandCalculator.hpp>

namespace EBC
{

BandCalculator::BandCalculator(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* sm, IndelModel* im, double divergenceTime) :
		fwd(4,nullptr), seq1(s1), seq2(s2), substModel(sm), indelModel(im), time(divergenceTime)
{
	DEBUG("Band estimator running...");

	posteriorLikelihoodLimit = Definitions::bandPosteriorLikelihoodLimit;
	posteriorLikelihoodDelta = Definitions::bandPosteriorLikelihoodDelta;

	this->ptMatrix =  new PMatrixDouble(substModel);
	this->trProbs = new TransitionProbabilities(indelModel);

	bool highDivergence = false;

	//standard multipliers
	array<double,4> lowMultipliers = {{0.7,0.9,1.1,1.3}};
	array<double,4> normalMultipliers = {{0.5,0.8, 1.1, 1.4}};
	array<double,4> highMultipliers = {{0.5,1.0,2.0,4.0}};

	array<double,4> &multipliers = lowMultipliers;

	accuracy = Definitions::highDivergenceAccuracyDelta;

	if(time < Definitions::kmerLowDivergence){
		band = new Band(s1->size(),s2->size(),0.075);
		INFO("LOW divergence");
		leftBound = Definitions::almostZero;
		rightBound = 2.0;
	}
	else if (time < Definitions::kmerHighDivergence){
		//multipliers = normalMultipliers;
		band = new Band(s1->size(),s2->size(),0.1);
		INFO("MEDIUM divergence");
		leftBound = Definitions::almostZero;
		rightBound = 5.0;
		//accuracy = Definitions::highDivergenceAccuracyDelta;
	}
	else{//very high divergence
		//multipliers = highMultipliers;
		band = new Band(s1->size(),s2->size(),0.25);
		INFO("HIGH divergence");
		leftBound = 0.5;
		//use value from
		rightBound = -1.0;
		//accuracy = Definitions::ultraDivergenceAccuracyDelta;
	}


	leftBound = Definitions::almostZero;
	rightBound = -1.0;

	//band = new Band(s1->size(),s2->size());

	unsigned int best = 0;
	double tmpRes = std::numeric_limits<double>::max();
	double lnl;

	DUMP("Trying several forward calculations to assess the band...");
	for(unsigned int i = 0; i < fwd.size(); i++)
	{

		fwd[i] = new ForwardPairHMM(seq1,seq2, substModel,indelModel, Definitions::DpMatrixType::Full,band);
		fwd[i]->setDivergenceTimeAndCalculateModels(time*multipliers[i]);
		lnl = fwd[i]->runAlgorithm();
		DUMP("Calculation "<< i << " with divergence time " << time*multipliers[i] << " and lnL " << lnl);
		if(lnl < tmpRes)
		{
			best = i;
			tmpRes = lnl;
		}
		//bwd[i] = new BackwardPairHMM(seq1,seq2, false, substModel,indelModel, 0, Definitions::DpMatrixType::Full);
		//bwd[i]->setDivergenceTimeAndCalculateModels(time*multipliers[i]);
		//bwd[i]->runAlgorithm();
	}

	//TODO - perhaps band it as well ???
	bwd =  new BackwardPairHMM(seq1,seq2, substModel,indelModel, Definitions::DpMatrixType::Full,band);
	bwd->setDivergenceTimeAndCalculateModels(time*multipliers[best]);
	DUMP("Backward calculation runs...");
	bwd->runAlgorithm();

	//combine fwd and bwd metrics into one!

	bwd->calculatePosteriors(fwd[best]);
	this->processPosteriorProbabilities(bwd, band);

	bestTime = time*multipliers[best];


}

BandCalculator::~BandCalculator()
{
	delete bwd;

	for(int i=0; i< fwd.size(); i++)
	{
		delete fwd[i];
	}

	delete trProbs;
	delete ptMatrix;

}

void BandCalculator::processPosteriorProbabilities(BackwardPairHMM* hmm, Band* band)
{
	//Match state
	PairwiseHmmStateBase* M;
	//Insert state
	PairwiseHmmStateBase* X;
	//Delete state
	PairwiseHmmStateBase* Y;

	M = hmm->getM();
	X = hmm->getX();
	Y = hmm->getY();

	//cumulative posterior likelihood
	double cpl = posteriorLikelihoodLimit + posteriorLikelihoodDelta;

	int xHi, xLo, yHi, yLo, mHi,mLo, rowCount;

	unsigned int tmpRow;

	double tmpX, tmpY, tmpM;

	rowCount = M->getRows();

	for(unsigned int col = 0; col < M->getCols(); col++)
	{
		xHi=mHi=yHi=xLo=mLo=yLo=-1;
		tmpRow = 0;
		while(tmpRow < rowCount && X->getValueAt(tmpRow, col) < cpl )
			tmpRow++;
		if (tmpRow != rowCount)
		{
			//found a value
			xLo = tmpRow;
			tmpRow = rowCount-1;
			while(tmpRow >= 0 && X->getValueAt(tmpRow, col) < cpl)
				tmpRow--;
			if (tmpRow > 0)
				xHi = tmpRow;
		}

		tmpRow = 0;
		while(tmpRow < rowCount && Y->getValueAt(tmpRow, col) < cpl)
			tmpRow++;
		if (tmpRow != rowCount)
		{
			//found a value
			yLo = tmpRow;
			tmpRow = rowCount-1;
			while(tmpRow >= 0 && Y->getValueAt(tmpRow, col) < cpl)
				tmpRow--;
			if (tmpRow > 0)
				yHi = tmpRow;
		}

		tmpRow = 0;
		while(tmpRow < rowCount && M->getValueAt(tmpRow, col) < cpl)
			tmpRow++;
		if (tmpRow != rowCount)
		{
			//found a value
			mLo = tmpRow;
			tmpRow = rowCount-1;
			while(tmpRow >= 0 && M->getValueAt(tmpRow, col) < cpl)
				tmpRow--;
			if (tmpRow > 0)
				mHi = tmpRow;
		}



		band->setInsertRangeAt(col, xLo,xHi);
		band->setDeleteRangeAt(col, yLo,yHi);
		band->setMatchRangeAt(col, mLo,mHi);

		DUMP("Match/Ins/Del bands for column " << col << "\t" << band->getMatchRangeAt(col).first <<"\t" << band->getMatchRangeAt(col).second
				<< "\t" << band->getInsertRangeAt(col).first <<"\t" << band->getInsertRangeAt(col).second
				<< "\t" << band->getDeleteRangeAt(col).first <<"\t" << band->getDeleteRangeAt(col).second);
	}
}

double BandCalculator::getClosestDistance() {
	return this->bestTime;
}

double BandCalculator::getBrentAccuracy() {
	return this->accuracy;
}

} /* namespace EBC */


