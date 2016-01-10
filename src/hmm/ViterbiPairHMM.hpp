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


#ifndef VITERBIFULLPAIRHMM_HPP_
#define VITERBIFULLPAIRHMM_HPP_

#include "hmm/EvolutionaryPairHMM.hpp"

namespace EBC
{

class ViterbiPairHMM: public EBC::EvolutionaryPairHMM
{

protected:

	vector<double> userIndelParameters;
	vector<double> userSubstParameters;

	vector<std::pair<unsigned int, unsigned int> > alignment;

	double getMax(double m, double x, double y, unsigned int i, unsigned int j, PairwiseHmmStateBase* state);

public:
	ViterbiPairHMM(vector<SequenceElement*>* s1, vector<SequenceElement*>* s2,
			SubstitutionModelBase* smdl, IndelModel* imdl,
			Definitions::DpMatrixType mt = Definitions::DpMatrixType::Full, Band* bandObj = nullptr, bool useEquilibriumFreqs = false);

	virtual ~ViterbiPairHMM();

	double runAlgorithm();

	double getViterbiSubstitutionLikelihood();


};

} /* namespace EBC */
#endif /* FORWARDPAIRHMM_HPP_ */
