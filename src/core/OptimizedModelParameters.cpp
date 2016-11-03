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


#include "core/OptimizedModelParameters.hpp"

namespace EBC
{

OptimizedModelParameters::OptimizedModelParameters(SubstitutionModelBase* sm, IndelModel* im, unsigned int sCount, unsigned int dCount,
		bool se, bool ie, bool ae, bool de, Maths* m) : maths(m), sm(sm), im(im), indelParameters(im != NULL ? im->getParamsNumber() : 0),
		substParameters(sm != NULL ? sm->getParamsNumber() : 0), divergenceTimes(dCount), indelHiBounds(im != NULL ? im->getParamsNumber() : 0),
		estimateIndelParams(ie), estimateSubstParams(se), estimateAlpha(ae), estimateDivergence(de)
{
	indelCount = indelParameters.size();
	substCount = substParameters.size();
	seqCount = sCount;
	distCount = dCount;
	optCount = (estimateSubstParams ? substCount : 0) +
			(estimateIndelParams ? indelCount : 0) + (estimateDivergence ? distCount : 0) + (estimateAlpha ? 1 : 0);

	this->divergenceBound  = Definitions::divergenceBound;

	if(im != NULL){
		for (int i =0; i< indelCount; i++){
			indelHiBounds[i] = im->getHiBound(i);
		}
	}
	//if(estimateIndelParams)
	//	generateInitialIndelParameters();
	//if(estimateSubstParams)
	//	generateInitialSubstitutionParameters();
	if (distCount > 0)
		generateInitialDistanceParameters();
}



void OptimizedModelParameters::useSubstitutionModelInitialParameters()
{
	if (sm != NULL){
		double* params = sm->getParameters();
		for(unsigned int i=0; i< substCount; i++)
		{
			substParameters[i] = params[i];
		}
	}

}
void OptimizedModelParameters::useIndelModelInitialParameters()
{
	if (im != NULL){
		double* params = im->getParameters();
		for(unsigned int i=0; i< indelCount; i++)
		{
			indelParameters[i] = params[i];
		}
	}
}

void OptimizedModelParameters::boundDivergenceBasedOnLambda(double lambda){
	//lambda * t must be smaller than negative ln(0.5)
	double db = ((log(0.5)*-1.0)/lambda) - 0.01; //substract 0.01 just to be sure

	this->divergenceBound = db <  Definitions::divergenceBound? db : Definitions::divergenceBound;
	DUMP("Optimised Model Parameters divergence bound : " << divergenceBound);
}

void OptimizedModelParameters::boundLambdaBasedOnDivergence(double time){
	//lambda * t must be smaller than negative ln(0.5)
	double ln05 = log(0.5)*-1.0;
	//subtract 0.00000011 to avoid the error with optimizer overrunning the bound
	indelHiBounds[0] = min((ln05/time) - (11*(Definitions::almostZero)), im->getHiBound(0));
	if (indelParameters[0] > indelHiBounds[0])
		indelParameters[0] = indelHiBounds[0];
	DUMP("Optimised Model Parameters lambda Hi bound : " << indelHiBounds[0] );
}

unsigned int OptimizedModelParameters::optParamCount()
{
	return optCount;
}

void OptimizedModelParameters::toDlibVector(column_vector& vals, column_vector& lbounds, column_vector& hbounds)
{
	unsigned int i;
	unsigned int ptr = 0;

	if(estimateSubstParams)
	{
		for (i=0; i < substCount; i++)
		{
			vals(i) = substParameters[i];
			//default probs bounds
			lbounds(i) = sm->getLoBound(i);
			hbounds(i) = sm->getHiBound(i);
		}
		ptr += substCount;
	}
	if(estimateIndelParams)
	{
		for (i=0; i < indelCount; i++)
		{
			vals(i+ptr) = indelParameters[i];
			//default probs bounds
			lbounds(i+ptr) = im->getLoBound(i);
			hbounds(i+ptr) = indelHiBounds[i];
		}
		ptr += indelCount;
	}
	if(estimateAlpha)
	{
		vals(ptr) = alpha;
		lbounds(ptr) =  Definitions::almostZero;
		hbounds(ptr) = Definitions::maxAlpha;
		ptr++;
	}
	if(estimateDivergence)
	{
		for (i=0; i < distCount; i++)
		{
			vals(i+ptr) = divergenceTimes[i];
			lbounds(i+ptr) = Definitions::almostZero;
			hbounds(i+ptr) = divergenceBound;
		}
	}
}

void OptimizedModelParameters::generateInitialSubstitutionParameters()
{
	for(unsigned int i=0; i< substCount; i++)
	{
		substParameters[i] = 0.2 + 0.1*maths->rndu();
	}
	DUMP("Model estimator initial substitution parameters:");
	DUMP(substParameters);
}

void OptimizedModelParameters::generateInitialIndelParameters()
{
	for(unsigned int i=0; i< indelCount; i++)
	{
		indelParameters[i] = 0.05 + 0.1*maths->rndu();
	}
	DUMP("Model estimator initial indel parameters:");
	DUMP(indelParameters);
}

void OptimizedModelParameters::generateInitialDistanceParameters()
{
	for(unsigned int i=0; i< distCount; i++)
	{
		divergenceTimes[i] = 0.2 + 0.1*maths->rndu();
	}

	DUMP("Model Estimator initial divergence times:");
	DUMP(divergenceTimes);
}

void OptimizedModelParameters::fromDlibVector(const column_vector& vals)
{
	unsigned int i;
	unsigned int ptr = 0;
	if(estimateSubstParams)
	{
		for (i=0; i < substCount; i++)
		{
			substParameters[i] = vals(i);
		}
		ptr += substCount;
	}
	if(estimateIndelParams)
	{
		for (i=0; i < indelCount; i++)
		{
			indelParameters[i] = vals(i+ptr);
		}
		ptr += indelCount;
	}
	if(estimateAlpha)
	{
		alpha = vals(ptr);
		ptr++;
	}
	if(estimateDivergence)
	{
		for (i=0; i < distCount; i++)
		{
			divergenceTimes[i] = vals(i+ptr);
		}
	}
}

void OptimizedModelParameters::setAlpha(double a)
{
	alpha = a;
}

void OptimizedModelParameters::setUserIndelParams(
		vector<double> allocator)
{
	indelParameters = allocator;
}

void OptimizedModelParameters::setUserDivergenceParams(
		vector<double> allocator)
{
	divergenceTimes = allocator;
}

void OptimizedModelParameters::setUserSubstParams(
		vector<double> allocator)
{
	substParameters = allocator;
}

void OptimizedModelParameters::logParameters()
{
	//for (auto p : substParameters)
	//		std::cerr << p  << '\t';
	//for (auto p : indelParameters)
	//		std::cerr << p  << '\t';
	//for (auto p : divergenceTimes)
	//	std::cerr << p  << '\t';
	//if(this->estimateAlpha)
	//	std::cerr << alpha;
	//std::cerr << std::endl;
	INFO(substParameters);
	INFO(indelParameters);
	INFO(divergenceTimes);
	if(this->estimateAlpha)
		INFO("Alpha: " << alpha);

}

void OptimizedModelParameters::outputToConsole()
{
	for (auto p : substParameters)
			std::cout << p  << '\t';
	for (auto p : indelParameters)
			std::cout << p  << '\t';
	for (auto p : divergenceTimes)
		std::cout << p  << '\t';
	if(this->estimateAlpha)
		std::cout << alpha;
	std::cout << std::endl;
}


double OptimizedModelParameters::getDistanceBetween(unsigned int i,
		unsigned int j)
{
	unsigned int total  = seqCount;
	unsigned int a,b;
	a=i;
	b=j;
	if (i>j)
	{
		b=i;
		a=j;
	}
	unsigned int index;
	if (i==j) return 0;
	else
	{
		index = (b - a - 1) + (a*total) - (((1+a)/2.0)*(a*1.0));
		return divergenceTimes[index];
	}
}

} /* namespace EBC */
