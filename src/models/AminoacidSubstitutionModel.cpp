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


#include "models/AminoacidSubstitutionModel.hpp"

namespace EBC
{

AminoacidSubstitutionModel::AminoacidSubstitutionModel(Dictionary* dict, Maths* alg, unsigned int alpha, Definitions::aaModelDefinition modelDef) :
	SubstitutionModelBase(dict,alg,alpha, Definitions::AAParamCount), eigenDecomposed(false)
{
	//FIXME - implement +F model with custom frequencies

	this->maxRate = modelDef.maxRate;
	this->piFreqs = new double[matrixSize];
	this->piLogFreqs = new double[matrixSize];

	std::copy(modelDef.aaRates,modelDef.aaRates + this->matrixFullSize, this->qMatrix);

/*
	std::copy(modelDef.aaFreqs,modelDef.aaFreqs + this->matrixSize, this->piFreqs);

	for(unsigned int i = 0; i< this->matrixSize; i++)
	{
			piLogFreqs[i] = log(piFreqs[i]);
	}

	this->setDiagonalMeans();
	this->doEigenDecomposition();
	eigenDecomposed=true;
*/
}


void AminoacidSubstitutionModel::calculateModel()
{
	if(this->eigenDecomposed == false)
	{
		this->setDiagonalMeans();
		this->doEigenDecomposition();
		eigenDecomposed=true;
	}
	//calculateGammaPtMatrices();
}

AminoacidSubstitutionModel::~AminoacidSubstitutionModel()
{
	//delete[] piFreqs;
	//delete[] piLogFreqs;
}

void AminoacidSubstitutionModel::summarize()
{
	//summarizeRates();
}

} /* namespace EBC */


