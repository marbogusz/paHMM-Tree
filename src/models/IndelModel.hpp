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


#ifndef INDELMODEL_HPP_
#define INDELMODEL_HPP_

#include <vector>

using namespace std;

namespace EBC
{

class IndelModel
{
protected:
	double gapExtensionProbability;

	double gapOpeningProbability;

	//divergence time
	double time;

	unsigned int paramsNumber;

	vector<double> parameterHiBounds;
	vector<double> parameterLoBounds;

	//bool logMode;

public:
	IndelModel(unsigned int);

	virtual double calculateGapOpening(double time) = 0;
	virtual double calculateGapExtension(double time) = 0;


	//set parameters - time + the rest of parameters
	virtual void setParameters(double*) = 0;

	virtual void setParameters(vector<double>)=0;

	void setTime(double t)
	{
		this->time = t;
	}

	virtual void summarize()=0;

	virtual void calculate()=0;

	double getGapExtensionProbability() const
	{
		return gapExtensionProbability;
	}

	double getGapOpeningProbability() const
	{
		return gapOpeningProbability;
	}

	virtual double* getParameters()=0;

	unsigned int getParamsNumber() const
	{
		return paramsNumber;
	}

	inline double getHiBound(unsigned int pos)
	{
		return parameterHiBounds[pos];
	}

	inline double getLoBound(unsigned int pos)
	{
		return parameterLoBounds[pos];
	}
};

} /* namespace EBC */
#endif /* INDELMODEL_HPP_ */
