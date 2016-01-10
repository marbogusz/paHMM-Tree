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


#ifndef S_MODEL_BASE_H_
#define S_MODEL_BASE_H_

#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/HmmException.hpp"
#include <cmath>
#include <vector>

namespace EBC
{

class SubstitutionModelBase
{

protected:

	Dictionary* dictionary;

	Maths* maths;

	//number of discrete gamma rate categories
	unsigned int rateCategories;

	unsigned int paramsNumber;

	unsigned int matrixSize;

	unsigned int matrixFullSize;

	//gamma distribution alpha parameter
	double alpha;

	//equilibrium frequencies
	double* piFreqs;

	double* piLogFreqs;

	//q-rates
	double* qMatrix;

	double* vMatrix;

	double* uMatrix;

	//square roots matrix for eigen decomposition
	double* squareRoots;

	//roots vector for eigen decomposition
	double* roots;

	//ML parameters - model parameters + time!
	double* parameters;

	//patterns for sequence elements + missing data (gaps)
	//double ** sitePatterns;

	//like site patterns but not logarithmic!
	//double ** siteProbabilities;

	double meanRate;

	vector<double> parameterHiBounds;
	vector<double> parameterLoBounds;

	//Allocate the memory;
	void allocateMatrices();

	void destroyMatrices();

	void doEigenDecomposition();

	void setDiagonalMeans();

	void calculateGamma();

	//void calculateGammaPtMatrices();

	//void summarizePatterns();

public:

	double* gammaFrequencies;

	double* gammaRates;

	SubstitutionModelBase(Dictionary*, Maths*i, unsigned int, unsigned int);

	virtual ~SubstitutionModelBase();

	virtual void summarize()=0;

	virtual void calculateModel()=0;

	double* calculatePt(double time, unsigned int rateCategory = 0);

	virtual void setObservedFrequencies(double* observedFrequencies);

	//double getPiXiPXiYi(unsigned int xi, unsigned int yi);

	//double getPXiYi(unsigned int xi, unsigned int yi);

	double getEquilibriumFrequencies(unsigned int xi);

	double getLogEquilibriumFrequencies(unsigned int xi);

	//double getSitePattern(unsigned int xi, unsigned int yi);

	//double getSiteProbability(unsigned int xi, unsigned int yi);

	virtual void setParameters(const vector<double>&)=0;

	void summarizeRates();

	//virtual void calculateSitePatterns();

	/*

	void setTime(double t)
	{
		this->time = t;
	}

	*/
//getters and setters

	inline unsigned int getRateCategories()
	{
		return this->rateCategories;
	}

	inline double getHiBound(unsigned int pos)
	{
		return parameterHiBounds[pos];
	}

	inline double getLoBound(unsigned int pos)
	{
		return parameterLoBounds[pos];
	}

	void setAlpha(double a)
	{
		if (this->alpha != a)
		{
			this->alpha = a;
			//cerr << "Alpha: " << a << endl;
			calculateGamma();
		}
	}

	double getAlpha()
	{
		return alpha;
	}


	void setParameters(const double* params)
	{
		for (unsigned int i=0;i<paramsNumber;i++)
			this->parameters[i] = params[i];
			//this->parameters[i] = fabs(params[i]);
	}


	unsigned int getParamsNumber() const 
	{
		return paramsNumber;
	}

	double* getParameters() 
	{
		return parameters;
	}

	unsigned int getMatrixSize() const
	{
		return matrixSize;
	}
};

} /* namespace EBC */
#endif /* MODEL_H_ */
