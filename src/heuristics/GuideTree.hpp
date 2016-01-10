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

#ifndef GUIDETREE_H_
#define GUIDETREE_H_


#include "core/Definitions.hpp"
#include "core/Maths.hpp"
#include "core/Dictionary.hpp"
#include "core/Sequences.hpp"
#include "core/HmmException.hpp"
#include "core/BioNJ.hpp"
#include <unordered_map>
#include <vector>
#include <string>

using namespace std;

namespace EBC
{

class GuideTree
{
protected:
	Dictionary* dict;
	Sequences* inputSequences;
	unsigned int kmerSize;
	unsigned int sequenceCount;
	vector<unordered_map<string,short>*>* kmers;
	vector<double> distances;

	vector<array<unsigned int, 3> > sampledTriplets;

	DistanceMatrix* distMat;

	string newickTree;


public:
	GuideTree(Sequences*);

	~GuideTree();

	void constructTree();

	DistanceMatrix* getDistanceMatrix()
	{
		return distMat;
	}

	const string& getNewickTree()
	{
		return newickTree;
	}

	vector<double> getDistances()
	{
		return distances;
	}

private:

	void extractKmers(string& seq, unordered_map<string,short>* umap);

	unsigned int commonKmerCount(unsigned int i, unsigned int j);

	double kimuraDist(double);

	inline double aaFunction(double x)
	{
		//fitted to a function.
		return (pow(100, x-1.04)+0.01)/0.6;

	}

	inline double nucFunction(double x)
	{
		return pow(700, x-0.95)+0.02;
	}
};

} /* namespace EBC */

#endif /* GUIDETREE_H_ */
