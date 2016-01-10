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
#include "hmm/ViterbiPairHMM.hpp"
#include <algorithm>

namespace EBC
{


ViterbiPairHMM::ViterbiPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl,
		IndelModel* imdl,Definitions::DpMatrixType mt, Band* bandObj, bool useEquilibriumFreqs ) :
		EvolutionaryPairHMM(s1,s2, smdl, imdl,mt, bandObj, useEquilibriumFreqs)
{
	this->alignment.reserve(xSize);
}

ViterbiPairHMM::~ViterbiPairHMM()
{
}

double ViterbiPairHMM::getViterbiSubstitutionLikelihood()
{
	double lnl = 0;
	calculateModels();
	//we have the site patterns now
	//get alignment!

	for (auto it = alignment.begin(); it != alignment.end(); it ++)
	{
		lnl += this->ptmatrix->getPairSitePattern(it->first, it->second);
	}

	return lnl * -1.0;
}

double ViterbiPairHMM::getMax(double m, double x, double y, unsigned int i, unsigned int j, PairwiseHmmStateBase* state)
{
	if(m >x && m >y)
	{
		//state->setDiagonalAt(i,j);

		//cout << i << " " << j << " coming from M" << endl;
		state->setDirection(i,j);
		state->setSourceMatrixPtr(i,j,M);
		return m;
	}
	else if(x > y)
	{

		//state->setVerticalAt(i,j);
		//cout << i << " " << j << " coming from X" << endl;
		state->setDirection(i,j);
		state->setSourceMatrixPtr(i,j,X);
		return x;
	}
	else
	{
		//state->setHorizontalAt(i,j);
		//cout << i << " " << j << " coming from Y" << endl;
		state->setDirection(i,j);
		state->setSourceMatrixPtr(i,j,Y);
		return y;
	}

}

double ViterbiPairHMM::runAlgorithm()
{
	unsigned int i,j,k,l;

	double xx,xy,xm,yx,yy,ym,mx,my,mm, sS;

	double emissionM;
	double emissionX;
	double emissionY;

	DUMP("Viterbi equilibriums : PiM\t" << piM << "\tPiI\t" << piI << "\tPiD\t" << piD);

	M->initializeData(this->piM);
	X->initializeData(this->piI);
	Y->initializeData(this->piD);

		//while (i != xSize && j != ySize)

	for (i = 0; i<xSize; i++)
	{
		for (j = 0; j<ySize; j++)
		{
			if(i!=0)
			{

				k = i-1;
				emissionX = ptmatrix->getLogEquilibriumFreq((*seq1)[i-1]->getMatrixIndex());
				xm = M->getValueAt(k,j) + X->getTransitionProbabilityFromMatch();
				xx = X->getValueAt(k,j) + X->getTransitionProbabilityFromInsert();
				xy = Y->getValueAt(k,j) + X->getTransitionProbabilityFromDelete();

				X->setValueAt(i,j,getMax(xm,xx,xy,i,j,X) + emissionX);
			}
			if(j!=0)
			{
				k = j-1;
				emissionY = ptmatrix->getLogEquilibriumFreq((*seq2)[j-1]->getMatrixIndex());
				ym = M->getValueAt(i,k) + Y->getTransitionProbabilityFromMatch();
				yx = X->getValueAt(i,k) + Y->getTransitionProbabilityFromInsert();
				yy = Y->getValueAt(i,k) + Y->getTransitionProbabilityFromDelete();
				Y->setValueAt(i,j,getMax(ym,yx,yy,i,j,Y) + emissionY);
			}

			if(i!=0 && j!=0)
			{
				k = i-1;
				l = j-1;
				emissionM = ptmatrix->getLogPairTransition((*seq1)[i-1]->getMatrixIndex(), (*seq2)[j-1]->getMatrixIndex());
				mm = M->getValueAt(k,l) + M->getTransitionProbabilityFromMatch();
				mx = X->getValueAt(k,l) + M->getTransitionProbabilityFromInsert();
				my = Y->getValueAt(k,l) + M->getTransitionProbabilityFromDelete();
				M->setValueAt(i,j,getMax(mm,mx,my,i,j,M) + emissionM);
			}
		}
	}
	mx = X->getValueAt(xSize-1,ySize-1);
	my = Y->getValueAt(xSize-1,ySize-1);
	mm = M->getValueAt(xSize-1,ySize-1);

/*
	if(mm >=mx && mm >=my)
	{
		M->tracebackRaw(this->(*seq1),this->(*seq2), this->dict, this->alignment);
	}
	else if(mx >= my)
	{
		X->tracebackRaw(this->(*seq1),this->(*seq2), this->dict, this->alignment);
	}
	else
	{
		Y->tracebackRaw(this->(*seq1),this->(*seq2), this->dict, this->alignment);
	}
*/
	DUMP("Final Viterbi M  " << mm);
	DUMP("Final Viterbi X  " << mx );
	DUMP("Final Viterbi Y  " << my );

return (std::max(mm,std::max(mx,my)))*-1.0;
}


} /* namespace EBC */
