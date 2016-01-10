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


#include "models/NucleotideSubstitutionModel.hpp"

namespace EBC
{

NucleotideSubstitutionModel::NucleotideSubstitutionModel(Dictionary* dict, Maths* alg, unsigned int alpha, unsigned int pcnt) :
	SubstitutionModelBase(dict,alg,alpha,pcnt)
{
	this->parameters = new double[this->paramsNumber];
}


void NucleotideSubstitutionModel::calculateModel()
{
	this->buildSmatrix();
	this->setDiagonalMeans();
	this->doEigenDecomposition();
}

NucleotideSubstitutionModel::~NucleotideSubstitutionModel()
{
	if (this->parameters)
		delete[] this->parameters;
}

} /* namespace EBC */

