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


#include "models/HKY85Model.hpp"

namespace EBC
{

HKY85Model::HKY85Model(Dictionary* dict, Maths* alg, unsigned int rates)
	: NucleotideSubstitutionModel(dict, alg, rates, Definitions::HKY85ParamCount)
{
	this->parameterHiBounds[0] = Definitions::gtrMaxRate;
	this->parameterLoBounds[0] = Definitions::almostZero;
}


void HKY85Model::setParameters(const vector<double>& par)
{
	for (unsigned int i = 0; i< paramsNumber; i++)
	{
		this->parameters[i] = par[i];
	}
}

void HKY85Model::buildSmatrix() {
	this->k = &parameters[0];

	int s = this->matrixSize;

	for (int i=0; i<s; i++)
		for (int j=0; j<s; j++)
		{
			if (i!=j)
			{
				qMatrix[(i*s)+j] = 1;
			}
		}
	qMatrix[1] = qMatrix[4]= qMatrix[11] = qMatrix[14] = parameters[0];

}

void HKY85Model::summarize()
{
	INFO("HKY85 model summary:");
	INFO("kappa " << *k );

	cout << "Alpha " << this->alpha << endl;
 	cout << "Kappa " << *k <<endl;


	INFO("Frequencies (T C A G)");
	INFO(this->piFreqs[0] << '\t' << this->piFreqs[1] << '\t' << this->piFreqs[2] << '\t' << this->piFreqs[3] << '\t');
	//summarizeRates();
}

} /* namespace EBC */
