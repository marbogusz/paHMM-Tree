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


#include "models/NegativeBinomialGapModel.hpp"
#include <cmath>
#include <iostream>
#include "core/Maths.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{

NegativeBinomialGapModel::NegativeBinomialGapModel() : IndelModel(Definitions::NBIndelParamCount)
{
	//logMode = true;
	this->parameterLoBounds[0] = this->parameterLoBounds[1] = Definitions::almostZero;
	this->parameterHiBounds[0] = Definitions::lambdaHiBound;
	this->parameterHiBounds[1] = Definitions::epsilonHiBound;
}

NegativeBinomialGapModel::~NegativeBinomialGapModel()
{
}

double NegativeBinomialGapModel::calculateGapOpening(double time)
{
	double NegLambdaT = -1.0*lambda*time;
	return 1.0-exp(NegLambdaT);
}

double NegativeBinomialGapModel::calculateGapExtension(double time)
{
	return this->gapExtensionProbability;
}

void NegativeBinomialGapModel::setParameters(vector<double> vc)
{
	params[0] = vc[0];
	params[1] = vc[1];
	this->lambda = vc[0];
	this->gapExtensionProbability = vc[1];
}

void NegativeBinomialGapModel::calculate()
{
	calculateGeometricProbability(this->lambda, this->time);
}


void NegativeBinomialGapModel::calculateGeometricProbability(double lambda, double t)
{
	double exponent = 1.0-exp(-1.0*lambda*t);
	this->gapOpeningProbability = exponent;
}

double* NegativeBinomialGapModel::getParameters()
{
	return this->params;
}

void NegativeBinomialGapModel::setParameters(double* par)
{

	params[0] = par[0];
	params[1] = par[1];
	this->lambda = par[0];

	this->gapExtensionProbability = par[1];

	//std::cerr << " Opening " << this->gapOpeningProbability << "\t" << "Extension " << this->gapExtensionProbability << std::endl;

}

void NegativeBinomialGapModel::summarize()
{
	//cout << endl << "Affine geometric gap model parameters:" << endl;
	cout << "lambda:\t" << lambda << "\t" << "ext.prob\t" <<  gapExtensionProbability << endl;
}

} /* namespace EBC */
