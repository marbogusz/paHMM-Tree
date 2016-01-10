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

#ifndef OPTIMIZEDMODELPARAMETERS_H_
#define OPTIMIZEDMODELPARAMETERS_H_

#include <vector>
#include <dlib/optimization.h>
#include "models/SubstitutionModelBase.hpp"
#include "core/Maths.hpp"
#include "models/IndelModel.hpp"

typedef dlib::matrix<double,0,1> column_vector;

using namespace std;

namespace EBC
{

class OptimizedModelParameters
{
protected:

	Maths* maths;

	SubstitutionModelBase* sm;
	IndelModel* im;

	vector<double> indelParameters;
	vector<double> substParameters;
	vector<double> divergenceTimes;
	double alpha;

	vector<double> indelHiBounds;

	bool estimateIndelParams;
	bool estimateSubstParams;
	bool estimateAlpha;
	bool estimateDivergence;

	unsigned int indelCount;
	unsigned int substCount;
	unsigned int distCount;
	unsigned int seqCount;
	unsigned int optCount;


public:
	double divergenceBound;

	OptimizedModelParameters(SubstitutionModelBase*, IndelModel*, unsigned int, unsigned int, bool, bool, bool, bool, Maths*);

	unsigned int optParamCount();

	void useSubstitutionModelInitialParameters();

	void useIndelModelInitialParameters();

	void boundDivergenceBasedOnLambda(double lambda);

	void boundLambdaBasedOnDivergence(double time);

	void toDlibVector(column_vector&,column_vector&,column_vector&);

	void fromDlibVector(const column_vector&);

	void setAlpha(double);

	void setUserIndelParams(vector<double>);

	void setUserSubstParams(vector<double>);

	void setUserDivergenceParams(vector<double>);

	void setSingleDivergenceParam(unsigned int pos, double val){
		divergenceTimes[pos] = val;
	}

	void logParameters();

	void outputToConsole();

	double getDistanceBetween(unsigned int i, unsigned int j);

	void generateInitialIndelParameters();
	void generateInitialSubstitutionParameters();
	void generateInitialDistanceParameters();

	double getAlpha() const
	{
		return alpha;
	}

	double getDivergenceTime(unsigned int index)
	{
		return divergenceTimes[index];
	}

	vector<double> getIndelParameters()
	{
		return indelParameters;
	}

	vector<double> getSubstParameters()
	{
		return substParameters;
	}

	vector<double> getDivergenceTimes()
	{
		return divergenceTimes;
	}
};

} /* namespace EBC */

#endif /* OPTIMIZEDMODELPARAMETERS_H_ */
