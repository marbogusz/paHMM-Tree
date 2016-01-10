//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// Most of the methods of discerte gamma and eigen decomposition taken from
// Ziheng Yang's PAML package, Copyright (c) Ziheng Yang
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


#ifndef ALGEBRA_H_
#define ALGEBRA_H_

#include <ctime>
#include <cstdlib>
#include <cmath>


namespace EBC
{

//Various FUNCTIONS

//FIXME - this class should be static!
class Maths
{
private:

	unsigned int z_rndu;

public:
	Maths();

	//Get a random within bounds
	double getRandom(double, double);

	//Exponentiate a reversible model matrix Q, store the result in P
	void revMatrixExponentiation(double* Q, double* P);

	//Ziheng Yang's PAML
	int eigenRealSym(double A[], int n, double Root[], double work[]);

	int eigenQREV (double Q[], double pi[], int n,
			double Root[], double U[], double V[], double spacesqrtpi[]);

	//From PAML
	void HouseholderRealSym(double a[], int n, double d[], double e[]);

	//Ziheng Yang's PAML fast random
	double rndu (void);

	//From PAML
	void EigenSort(double d[], double U[], int n);

	//From PAML
	int EigenTridagQLImplicit(double d[], double e[], int n, double z[]);

	//returns new matrix
	double* matrixMultiply(double *matA, double *matB, int size);

	//modifies an existing one
	double* matrixByDiagonalMultiply(double *matA, double *matDiag, int size);

	void matrixByDiagonalMultiplyMutable(double *matA, double *matDiag, int size);

	//modifies existing one
	void vectorMultiply(double* vector, double factor, int size);

	//append B to A (add matrix components
	void matrixAppend(double* matA, double* matB, int size);

	void matrixScale(double* matrix, double factor, int size);

	double* expLambdaT(double* lambda, double t, int size);

	double logSum(double, double, double);

	double logSum(double, double);

	//from PAML,
	double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);

	//PAML
	double QuantileNormal (double prob);
	//PAML
	double QuantileChi2 (double prob, double v);
	//PAML
	long factorial (int n);
	//PAML
	double LnGamma (double x);
	//PAML
	int  DiscreteGamma (double freqK[], double rK[], double alpha, double beta, int K, int UseMedian);


	inline static double logistic(double x)
	{
		return 1.0/(1+exp(x*-1.0));
	}

	inline static double logit(double x)
	{
		return log(x/1.0-x);
	}

};

} /* namespace EBC */
#endif /* ALGEBRA_H_ */
