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


#include "heuristics/SubstitutionModelEstimator.hpp"
#include "models/GTRModel.hpp"
#include "models/HKY85Model.hpp"
#include "models/AminoacidSubstitutionModel.hpp"
#include "heuristics/TripletSamplingTree.hpp"
#include "heuristics/TripletAligner.hpp"
#include "heuristics/StateTransitionEstimator.hpp"

namespace EBC
{

SubstitutionModelEstimator::SubstitutionModelEstimator(Sequences* inputSeqs, SubstitutionModelBase* model,
		Definitions::OptimizationType ot,unsigned int rateCategories, double alpha,
		bool estimateAlpha, unsigned int matCount) :
				inputSequences(inputSeqs), substModel(model), gammaRateCategories(rateCategories),
				patterns(matCount), ptMatrices(matCount), distances(3*matCount)


{

	DEBUG("Starting Substitution Model Estimator (SME)");
	DUMP("SME estimate alpha : " << estimateAlpha << " alpha value " << alpha);
	DUMP("SME rate categories : " << rateCategories << " triplets number " << matCount);

	this->estimateSubstitutionParams = true;
	this->estimateAlpha = estimateAlpha;
	this->alpha = alpha;

	currentTriplet = 0;

	maths = new Maths();
	dict = inputSequences->getDictionary();

	modelParams = new OptimizedModelParameters(substModel, NULL,3, 3*patterns.size(), estimateSubstitutionParams,
			false, estimateAlpha, true, maths);
	modelParams->useSubstitutionModelInitialParameters();

	modelParams->setAlpha(alpha);

	for(int i = 0; i < ptMatrices.size(); i++){
		DUMP("SME: creating ptMatrix");
		ptMatrices[i][0]  = new PMatrixTriple(substModel);
		ptMatrices[i][1]  = new PMatrixTriple(substModel);
		ptMatrices[i][2]  = new PMatrixTriple(substModel);
	}

	bfgs = new Optimizer(modelParams, this,ot);
}

void SubstitutionModelEstimator::clean(int nelems)
{
	//for(auto entry : ptMatrices)
	//{
	//	delete entry[0];
	//	delete entry[1];
	//	delete entry[2];
	//}
	if(nelems > 0){
		for(int i=0; i<nelems; i++){
			delete ptMatrices[i][0];
			delete ptMatrices[i][1];
			delete ptMatrices[i][2];
		}
		patterns.erase(patterns.begin(),patterns.begin() + nelems);
		ptMatrices.erase(ptMatrices.begin(),ptMatrices.begin() + nelems);
		distances.erase(distances.begin(), distances.begin() + (nelems*3));

	}

	for (auto itp = patterns.begin(); itp != patterns.end(); itp++)
	{
		for (auto itMap = itp->begin(); itMap != itp->end(); itMap++)
		{
			itMap->second = 0.0;
		}
	}
}

SubstitutionModelEstimator::~SubstitutionModelEstimator()
{
	delete bfgs;
	delete modelParams;
	//delete substModel;
	delete maths;

	for(auto entry : ptMatrices)
	{
		delete entry[0];
		delete entry[1];
		delete entry[2];
	}
}

void SubstitutionModelEstimator::addTriplet(array<vector<unsigned char>*, 3> tripleAlignment, unsigned int trp,
		double d1, double d2, double d3)
{
	DUMP("SME : adding patterns for triplet " << trp);
	for(int pos = 0; pos < tripleAlignment[0]->size(); pos++)
	{
		patterns[trp][{{(*tripleAlignment[0])[pos], (*tripleAlignment[1])[pos],(*tripleAlignment[2])[pos]}}] += 1;
	}
	//for (auto pat : patterns[trp])
	//{
	//	DUMP("SME " << pat.first[0] << " " << pat.first[1] << " " << pat.first[2] << " : " << pat.second);
	//}
	distances[trp*3] = d1 > 0 ? d1 : 0.5;
	distances[trp*3 + 1] = d2 > 0 ? d2 : 0.5;
	distances[trp*3 + 2] = d3 > 0 ? d3 : 0.5;
}

void SubstitutionModelEstimator::optimize()
{

	modelParams->setUserDivergenceParams(distances);
	bfgs->optimize();
	INFO("SubstitutionModelEstimator results:");

	substModel->setAlpha(modelParams->getAlpha());
	substModel->setParameters(modelParams->getSubstParameters());

	modelParams->logParameters();
/*
	cerr << "kappa\tlnL" << endl;

	double kappa = 0.5;
	double result;
	double partial1;
	while (kappa < 5.0)
	{
		result = 0;
		substModel->setAlpha(modelParams->getAlpha());
		substModel->setParameters({kappa});
		substModel->calculateModel();

		for (unsigned int i = 0; i< ptMatrices.size(); i++)
		{
			for(unsigned int j=0;j<Definitions::heuristicsTreeSize;j++)
			{
				ptMatrices[i][j]->setTime(modelParams->getDivergenceTime(Definitions::heuristicsTreeSize*i +j));
				ptMatrices[i][j]->calculate();
			}
		}

		for(int al = 0; al < patterns.size(); al++)
		{
			for(auto it : patterns[al])
			{
				partial1 = 0;
				for(int rt = 0; rt < dict->getAlphabetSize(); rt++)
				{
					partial1 += ptMatrices[al][0]->getTripleSitePattern(rt,it.first, ptMatrices[al][1],ptMatrices[al][2]);
				}
				result += log(partial1)* it.second;

			}
		}
		cerr << kappa << "\t" << result << endl;
		kappa += 0.025;
	}
*/

}

double SubstitutionModelEstimator::runIteration()
{
	double result = 0;

	double partial1;

	substModel->setAlpha(modelParams->getAlpha());
	substModel->setParameters(modelParams->getSubstParameters());
	substModel->calculateModel();

	for (unsigned int i = 0; i< ptMatrices.size(); i++)
	{
		for(unsigned int j=0;j<Definitions::heuristicsTreeSize;j++)
		{
			ptMatrices[i][j]->setTime(modelParams->getDivergenceTime(Definitions::heuristicsTreeSize*i +j));
			ptMatrices[i][j]->calculate();
		}
	}

	for(int al = 0; al < patterns.size(); al++)
	{
		for(auto it : patterns[al])
		{
			partial1 = 0;
			for(int rt = 0; rt < dict->getAlphabetSize(); rt++)
			{
				partial1 += ptMatrices[al][0]->getTripleSitePattern(rt,it.first, ptMatrices[al][1],ptMatrices[al][2]);
			}
			result += log(partial1)* it.second;

		}
	}
	//DEBUG("lnl result:" << result);
	//cerr << modelParams->getSubstParameters()[0] << "\t" << result << endl;

	return result * -1.0;



}

} /* namespace EBC */
