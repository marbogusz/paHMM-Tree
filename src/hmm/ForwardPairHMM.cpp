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


#include "core/Definitions.hpp"
#include "hmm/ForwardPairHMM.hpp"

namespace EBC
{


ForwardPairHMM::ForwardPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl,
		IndelModel* imdl, Definitions::DpMatrixType mt, Band* bandObj, bool useEquilibriumFreqs) :
		EvolutionaryPairHMM(s1,s2, smdl, imdl, mt, bandObj, true)
{
}

ForwardPairHMM::~ForwardPairHMM()
{
}

double ForwardPairHMM::runAlgorithm()
{

	int i;
	int j;
	int k;
	int l;

	double sX,sY,sM, sS;

	double xx,xy,xm,yx,yy,ym,mx,my,mm;

	double emissionM;
	double emissionX;
	double emissionY;

	//TODO - multiple runs using the same hmm object do not require dp matrix zeroing as long as the band stays the same!

	//DUMP("Forward equilibriums : PiM\t" << piM << "\tPiI\t" << piI << "\tPiD\t" << piD);

	M->initializeData(this->piM);
	X->initializeData(this->piI);
	Y->initializeData(this->piD);

	if(this->band == NULL)
	{
		//handle 0-index rows and columns separately!
		//1st col
		X->setValueAt(1,0, ptmatrix->getLogEquilibriumFreqClass((*seq1)[0]) + initTransX);

		for(i=2,j=0; i< xSize; i++)
		{
			k = i-1;
			emissionX = ptmatrix->getLogEquilibriumFreqClass((*seq1)[i-1]);
			xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
			X->setValueAt(i,j, emissionX + xx);
		}
		//1st row

		Y->setValueAt(0,1, ptmatrix->getLogEquilibriumFreqClass((*seq2)[0]) + initTransY);
		for(j=2,i=0; j< ySize; j++)
		{
			k = j-1;
			emissionY = ptmatrix->getLogEquilibriumFreqClass((*seq2)[j-1]);
			yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
			Y->setValueAt(i,j, emissionY + yy);
		}


		for (i = 1; i<xSize; i++)
		{
			for (j = 1; j<ySize; j++)
			{

				k = i-1;
				emissionX = ptmatrix->getLogEquilibriumFreqClass((*seq1)[i-1]);
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
				X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));

				k = j-1;
				emissionY = ptmatrix->getLogEquilibriumFreqClass((*seq2)[j-1]);
				ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
				yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
				yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
				Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));

				k = i-1;
				l = j-1;
				emissionM = ptmatrix->getLogPairTransitionClass((*seq1)[i-1], (*seq2)[j-1]);
				mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
				mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
				my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
				M->setValueAt(i,j, emissionM + maths->logSum(mm,mx,my));
			}
		}
	}
	else
	{
		//banding column by column!
		//calculate 1st column for X
		int loI, hiI, loD, hiD, loM, hiM;
		auto bracket = band->getInsertRangeAt(0);
		loI = bracket.first;
		if (loI > 0)
		{
			hiI = bracket.second;
			for(i=loI,j=0; i<= hiI; i++)
			{
				k = i-1;
				emissionX = ptmatrix->getLogEquilibriumFreqClass((*seq1)[i-1]);
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
				X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));
			}
		}
		for(j=1; j<ySize; j++)
		{
			//TODO - range should use a reference perhaps
			auto bracketI = band->getInsertRangeAt(j);
			auto bracketD = band->getDeleteRangeAt(j);
			auto bracketM = band->getMatchRangeAt(j);
			loI = bracketI.first;
			loM = bracketM.first;
			loD = bracketD.first;

			if (loD > -1)
			{
				hiD = bracketD.second;
				for(i = loD; i <= hiD; i++)
				{
					k = j-1;
					emissionY = ptmatrix->getLogEquilibriumFreqClass((*seq2)[j-1]);
					ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
					yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
					yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
					Y->setValueAt(i,j, emissionY + maths->logSum(ym,yx,yy));
				}
			}
			if (loM > 0)
			{
				hiM = bracketM.second;
				for(i = loM; i <= hiM; i++)
				{
					k = i-1;
					l = j-1;
					emissionM = ptmatrix->getLogPairTransitionClass((*seq1)[i-1], (*seq2)[j-1]);
					mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
					mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
					my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
					M->setValueAt(i,j, emissionM + maths->logSum(mm,mx,my));
				}
			}

			if (loI > 0)
			{
				hiI = bracketI.second;
				for(i = loI; i <= hiI; i++)
				{
					k = i-1;
					emissionX = ptmatrix->getLogEquilibriumFreqClass((*seq1)[i-1]);
					xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
					xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
					xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();
					X->setValueAt(i,j, emissionX + maths->logSum(xm,xx,xy));
				}
			}
		}
	}

	sM = M->getValueAt(xSize-1, ySize-1);
	sX = X->getValueAt(xSize-1, ySize-1);
	sY = Y->getValueAt(xSize-1, ySize-1);

	sS = maths->logSum(sM,sX,sY) + log(xi);

	this->setTotalLikelihood(sS);

	//cerr << "\t" << sX << "\t" << sY << "\t"<< sM << "\t" << sS << endl;

	DUMP ("Forward lnls I, D, M, Total " << sX << "\t" << sY << "\t" << sM << "\t" << sS);

	return sS* -1.0;
}



} /* namespace EBC */
