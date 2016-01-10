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


#include "models/SubstitutionModelBase.hpp"
#include <sstream>

namespace EBC
{

SubstitutionModelBase::SubstitutionModelBase(Dictionary* dict, Maths* alg, unsigned int rateCategories, unsigned int parameter_count)
	: dictionary(dict), maths(alg), rateCategories(rateCategories), paramsNumber(parameter_count),
	  parameterHiBounds(parameter_count), parameterLoBounds(parameter_count)
{
	this->matrixSize = dict->getAlphabetSize();
	this->matrixFullSize = matrixSize*matrixSize;
	this->allocateMatrices();
	this->meanRate=0.0;
	this-> roots = new double[matrixSize];
	//no alpha provided
	this->alpha = 0;
	this->piFreqs = NULL;
	this->piLogFreqs = NULL;
	this->parameters = NULL;
}


void SubstitutionModelBase::allocateMatrices()
{
	this->qMatrix = new double[matrixFullSize];
	this->uMatrix = new double[matrixFullSize];
	this->vMatrix = new double[matrixFullSize];
	this->squareRoots = new double[matrixFullSize];

	if (rateCategories > 1)
	{
		gammaFrequencies = new double[rateCategories];
		gammaRates = new double[rateCategories];
	}
	else
	{
		gammaFrequencies = new double[1];
		gammaRates = new double[1];
		gammaFrequencies[0] = gammaRates[0] = 1;
		meanRate = 1.0;
	}

	/*
	this->sitePatterns = new double*[matrixSize+1];
	this->siteProbabilities = new double*[matrixSize+1];
	for (int i =0; i<= matrixSize; i++ )
	{
		sitePatterns[i] = new double[matrixSize+1];
		siteProbabilities[i] = new double[matrixSize+1];;
	}
	*/
}


void SubstitutionModelBase::destroyMatrices()
{
	delete [] qMatrix;
	delete [] uMatrix;
	delete [] vMatrix;
	delete [] squareRoots;
	if (this->gammaFrequencies !=NULL )
		delete[] gammaFrequencies;
	if (this->gammaRates !=NULL )
		delete[] gammaRates;

}

void SubstitutionModelBase::doEigenDecomposition()
{
	std::fill(roots, roots+matrixSize, 0);
	std::fill(uMatrix, uMatrix+matrixFullSize, 0);
	std::fill(vMatrix, vMatrix+matrixFullSize, 0);
	std::fill(squareRoots, squareRoots+matrixFullSize, 0);
	this->maths->eigenQREV(qMatrix, piFreqs, matrixSize, roots, uMatrix, vMatrix, squareRoots);
}

double* SubstitutionModelBase::calculatePt(double t, unsigned int rateCategory)
{
	double *tmpRoots, *tmpUroots, *matrix;
	tmpRoots = maths->expLambdaT(roots, t*gammaRates[rateCategory], matrixSize);
	tmpUroots = maths->matrixByDiagonalMultiply(uMatrix, tmpRoots, matrixSize);
	matrix = this->maths->matrixMultiply(tmpUroots, vMatrix, matrixSize);
	delete[] tmpRoots;
	delete[] tmpUroots;
	return matrix;
}

void SubstitutionModelBase::setDiagonalMeans()
{
		unsigned int i,j;
		double sum;

		maths->matrixByDiagonalMultiplyMutable(qMatrix,piFreqs,matrixSize);
		meanRate = 0;

		for (i=0; i< this->matrixSize; i++)
		{
			qMatrix[i*matrixSize+i] = 0;
			sum = 0.0;
			for (j=0; j < this->matrixSize; j++)
			{
				sum -= qMatrix[i*matrixSize + j];
			}
			qMatrix[i*matrixSize + i] = sum;
			meanRate -= sum*piFreqs[i];
		}
		for (i=0; i< this->matrixSize; i++)
		{
				for (j=0; j < this->matrixSize; j++)
				{
					qMatrix[i*matrixSize + j] /= meanRate;
				}
		}

}

SubstitutionModelBase::~SubstitutionModelBase()
{
	delete[] this-> roots;

	if (this->piLogFreqs != NULL)
		delete[] this->piLogFreqs;
	destroyMatrices();

	//site patterns
	/*
	for (int i =0; i<= matrixSize; i++ )
	{
			delete[] sitePatterns[i];
			delete[] siteProbabilities[i];
	}
	delete[] sitePatterns;
	delete[] siteProbabilities;
	 */
}

void SubstitutionModelBase::calculateGamma()
{
	int useMedian = 0;
	if(!(alpha <=0 || rateCategories == 1))
	{
		this->maths->DiscreteGamma(gammaFrequencies, gammaRates, alpha, alpha, rateCategories, useMedian);
	}
}

void SubstitutionModelBase::setObservedFrequencies(double* observedFrequencies)
{
	//Size and order is clear
	this->piFreqs = observedFrequencies;
	if(this->piLogFreqs == NULL)
		this->piLogFreqs = new double[this->matrixSize];
	for(unsigned int i = 0; i< this->matrixSize; i++)
	{
		piLogFreqs[i] = log(piFreqs[i]);
	}
}

double SubstitutionModelBase::getEquilibriumFrequencies(unsigned int xi)
{
	return piFreqs[xi];
}

double SubstitutionModelBase::getLogEquilibriumFrequencies(unsigned int xi)
{
	return piLogFreqs[xi];
}

void SubstitutionModelBase::summarizeRates()
{
	if(this->rateCategories != 1)
	{
		cout << "Rate categories: " << endl;
		for(unsigned int r=0; r< rateCategories; r++)
		{
			cout << gammaRates[r] << "\t\t";
		}
		cout << endl;
		cout << "Rate frequencies : " << endl;
		for(unsigned int f=0; f< rateCategories; f++)
		{
			cout << this->gammaFrequencies[f] << "\t\t";
		}
		cout << endl << "Alpha : " << alpha;
		cout << endl;
	}

	stringstream matrixstr;
	matrixstr << "Rate Matrix (T C A G): " << endl;

		for (unsigned int i =0; i< matrixSize; i++ )
		{
				for (unsigned int j =0; j<  matrixSize; j++ )
				{
					matrixstr << qMatrix[(i*matrixSize)+j] << "\t\t";
				}
				matrixstr << endl;
		}
	INFO(matrixstr.str());
}

} /* namespace EBC */


