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



#include <core/PMatrix.hpp>

namespace EBC
{

PMatrix::PMatrix(SubstitutionModelBase* m) : model(m),  matrixSize(m->getMatrixSize()) ,time(0),
		rateCategories(model->getRateCategories()), ptMatrices(rateCategories, nullptr)
{
	this->matrixFullSize = matrixSize*matrixSize;
}

PMatrix::~PMatrix()
{
	for (int i = 0; i < ptMatrices.size(); i++)
	{
		delete [] ptMatrices[i];
	}
}

void PMatrix::setTime(double t)
{
		this->time = t;
}

void PMatrix::summarize()
{
	cout << "P(t) matrix summary :" << endl;
	cout << "Divergence time : " << time << endl;

	cout << "Pairwise Site patterns " << endl;
}

} /* namespace EBC */


