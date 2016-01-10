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


#ifndef BACKWARDPAIRHMM_HPP_
#define BACKWARDPAIRHMM_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"
#include "hmm/ForwardPairHMM.hpp"


namespace EBC
{

class BackwardPairHMM: public EBC::EvolutionaryPairHMM
{

protected:

	//maximum posterior state;
	PairwiseHmmStateBase* MPstate;

	inline bool withinBand(unsigned int line, int position, unsigned int width)
	{
		int low = line - width;
		int high = line + width;
		bool result = ((position >= low) && (position <= high));
		return result;
	}


public:
	BackwardPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2, SubstitutionModelBase* smdl, IndelModel* imdl,
			Definitions::DpMatrixType mt, Band* bandObj = nullptr);

	virtual ~BackwardPairHMM();

	double runAlgorithm();

	void calculatePosteriors(ForwardPairHMM* fwd);

	void calculateMaximumPosteriorMatrix();

	//maximum posteriori alignment
	pair<string, string> getMPAlignment();

	pair<vector<double>*, pair<vector<unsigned char>*, vector<unsigned char>*> >
	getMPDWithPosteriors();

};

} /* namespace EBC */
#endif /* BACKWARDPAIRHMM_HPP_ */
