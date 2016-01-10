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

#include "hmm/EvolutionaryPairHMM.hpp"
#include "hmm/DpMatrixFull.hpp"
#include "models/NegativeBinomialGapModel.hpp"

namespace EBC
{

EvolutionaryPairHMM::EvolutionaryPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2,
			SubstitutionModelBase* smdl, IndelModel* imdl,
			Definitions::DpMatrixType mt, Band* bandObj, bool useEquilibriumFreqs) :
					substModel(smdl), indelModel(imdl), band(bandObj), equilibriumFreqs(useEquilibriumFreqs)
{
	//TODO - unfix the xi terminal probability  ????
	//length distribution fixed
	xi = 0.001;

	M = X = Y = NULL;
	this->seq1 = s1;
	this->seq2 = s2;

	this->xSize = seq1->size() +1;
	this->ySize = seq2->size() +1;

	DEBUG("#######Evolutionary Pair HMM constructor for seqence 1 with size " << xSize << " and sequence 2 with size " << ySize);

	ptmatrix = new PMatrixDouble(substModel);

	this->tpb = new TransitionProbabilities(indelModel);

	piM = 0;
	piI = Definitions::minMatrixLikelihood;
	piD = Definitions::minMatrixLikelihood;

	//piM = Definitions::minMatrixLikelihood;
	//piI = 0;
	//piD = Definitions::minMatrixLikelihood;

	initTransX = initTransY = initTransM = 0;

	initializeStates(mt);
}

void EvolutionaryPairHMM::setDivergenceTimeAndCalculateModels(double time)
{
	ptmatrix->setTime(time);
	tpb->setTime(time);

	calculateModels();
	setTransitionProbabilities();
	if (this->equilibriumFreqs)
		this->getStateEquilibriums();


}

void EvolutionaryPairHMM::getStateEquilibriums()
{
	double minPi = exp(Definitions::minMatrixLikelihood);

	//TODO - change indices to state ids from the enum!!!
	md[0][0] = (1.0-2*g);
	md[1][1] = e+((1.0-e)*g);
	md[2][2] = e+((1.0-e)*g);
	md[0][1] = g;
	md[0][2] = g;

	md[1][0] = (1.0-e)*(1-2*g);
	md[2][0] = (1.0-e)*(1-2*g);

	md[2][1] = (1.0-e)*g;
	md[1][2] = (1.0-e)*g;

	piD = ((1.0-md[0][0])+(md[0][1]*(1.0-md[0][0]+md[1][0])/(md[1][1]-1.0-md[0][1])))/(((md[0][1]-md[2][1])*(1.0-md[0][0]+md[1][0])/(md[1][1]-1.0-md[0][1]))+md[2][0]-md[0][0]+1);
	piI = ((piD*(md[0][1]-md[2][1]))-md[0][1])/(md[1][1]-1.0-md[0][1]);
	piM = 1.0 -piI - piD;

	DUMP("Decimal equilibriums : PiM\t" << piM << "\tPiI\t" << piI << "\tPiD\t" << piD);

	//substract xi/3 prob
	//double xx, yy;
	//xx = yy = log((e * piI) - (xi/3.0));
	//double mm = log((((1.0-xi)*(1.0-2*g))*piM) - (xi/3.0));

	//xx = yy = log((e * piI));
	//double mm = log((((1.0-xi)*(1.0-2*g))*piM));

	piD = (piD- (xi/3.0)) < minPi ? Definitions::minMatrixLikelihood : log(piD- (xi/3.0));
	piI = (piI- (xi/3.0)) < minPi ? Definitions::minMatrixLikelihood : log(piI- (xi/3.0));
	piM = (piM- (xi/3.0)) < minPi ? Definitions::minMatrixLikelihood : log(piM- (xi/3.0));

	initTransX = maths->logSum(X->getTransitionProbabilityFromInsert() + piI, X->getTransitionProbabilityFromDelete() + piD, X->getTransitionProbabilityFromMatch() + piM);
	initTransY = maths->logSum(Y->getTransitionProbabilityFromInsert() + piI, Y->getTransitionProbabilityFromDelete() + piD, Y->getTransitionProbabilityFromMatch() + piM);
	initTransM = maths->logSum(M->getTransitionProbabilityFromInsert() + piI, M->getTransitionProbabilityFromDelete() + piD, M->getTransitionProbabilityFromMatch() + piM);



	md[0][0] = log(md[0][0]);
	md[1][1] = log(md[1][1]);
	md[2][2] = log(md[2][2]);
	md[0][1] = log(md[0][1]);
	md[0][2] = log(md[0][2]);

	md[1][0] = log(md[1][0]);
	md[2][0] = log(md[2][0]);

	md[2][1] = log(md[2][1]);
	md[1][2] = log(md[1][2]);

	//DUMP("Initial transition likelihood component : M\t" << initTransM << "\tI\t" << initTransX << "\tD\t" << initTransY);
}

void EvolutionaryPairHMM::setTransitionProbabilities()
{
	e = tpb->getGapExtension();
	g = tpb->getGapOpening();

	M->setTransitionProbabilityFromMatch(log((1-2*g)*(1-xi)));
	M->setTransitionProbabilityFromInsert(log((1-e-xi)*(1-2*g)));
	M->setTransitionProbabilityFromDelete(log((1-e-xi)*(1-2*g)));

	X->setTransitionProbabilityFromInsert(log(e+((1-e-xi)*g)));
	Y->setTransitionProbabilityFromDelete(log(e+((1-e-xi)*g)));

	X->setTransitionProbabilityFromDelete(log((1-e-xi)*g));
	Y->setTransitionProbabilityFromInsert(log((1-e-xi)*g));

	X->setTransitionProbabilityFromMatch(log(g*(1-xi)));
	Y->setTransitionProbabilityFromMatch(log(g*(1-xi)));

	/*
	DUMP(" Transition probabilities: ");
	DUMP("M->M : " << log(1-2*g));
	DUMP("I->I : " << log(e+((1-e)*g)));
	DUMP("M->I : " << log(g));
	DUMP("I->M : " << log((1-2*g)*(1-e)));
	DUMP("I->D : " << log((1-e)*g));
*/

}

void EvolutionaryPairHMM::summarize()
{
	e = tpb->getGapExtension();
	g = tpb->getGapOpening();

	DUMP(" Transition probabilities: ");
	DUMP("M->M : " << M->getTransitionProbabilityFromMatch());
	DUMP("I->I : " << X->getTransitionProbabilityFromInsert());
	DUMP("M->I : " << X->getTransitionProbabilityFromMatch());
	DUMP("I->M : " << M->getTransitionProbabilityFromInsert());
	DUMP("I->D : " << Y->getTransitionProbabilityFromInsert());

	indelModel->summarize();
	substModel->summarize();

}

void EvolutionaryPairHMM::initializeStates(Definitions::DpMatrixType mt)
{

	if (M != NULL)
		delete M;
	if (X != NULL)
		delete X;
	if (Y != NULL)
		delete Y;

	switch (mt)
	{
	case Definitions::DpMatrixType::Full :
		M = new PairwiseHmmMatchState(xSize,ySize);
		X = new PairwiseHmmInsertState(xSize,ySize);
		Y = new PairwiseHmmDeleteState(xSize,ySize);
		break;
	case Definitions::DpMatrixType::Limited :
		M = new PairwiseHmmMatchState(new DpMatrixLoMem(xSize,ySize));
		X = new PairwiseHmmInsertState(new DpMatrixLoMem(xSize,ySize));
		Y = new PairwiseHmmDeleteState(new DpMatrixLoMem(xSize,ySize));
		break;
	default :
		M = new PairwiseHmmMatchState(xSize,ySize);
		X = new PairwiseHmmInsertState(xSize,ySize);
		Y = new PairwiseHmmDeleteState(xSize,ySize);
	}
}

void EvolutionaryPairHMM::calculateModels()
{
	ptmatrix->calculate();
	tpb->calculate();
}

EvolutionaryPairHMM::~EvolutionaryPairHMM()
{

	DUMP("~~~~~~~Evolutionary pair HMM destructor");
	delete Y;
	delete X;
	delete M;
    delete ptmatrix;
    delete tpb;
}



//OBSOLETE
double EvolutionaryPairHMM::getAlignmentLikelihood(vector<unsigned char>* s1,
		vector<unsigned char>* s2, Dictionary* dict)
{
	double lnl = 0;
	return lnl;

}

} /* namespace EBC */


