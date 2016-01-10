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


#ifndef TRIPLETSAMPLINGTREE_HPP_
#define TRIPLETSAMPLINGTREE_HPP_

#include <vector>
#include <array>
#include <map>
#include <unordered_set>
#include "heuristics/Node.hpp"
#include "core/DistanceMatrix.hpp"
#include "heuristics/GuideTree.hpp"

using namespace std;

namespace EBC
{

class TripletSamplingTree
{
private:
	//all nodes
	vector<Node*> nodes;

	map<unsigned int, Node*> leafNodes;

	map<unsigned int, Node*> availableNodes;

	unordered_set<Node*> usedNodes;

	Node* mostRecentAncestor(Node* n1, Node* n2);

	double distanceToParent(Node* n1, Node* par);

	inline bool isWithinRange(double val, double ran)
	{
		if(val <= ran*1.2 || val >= ran*0.8)
			return true;
		return false;
	}

	inline bool isWithinRange(double val, pair<double, double> ran)
	{
		if(val >= ran.first && val <= ran.second)
			return true;
		return false;
	}

	DistanceMatrix* distMat;

	double idealTreeSize;

	double averageLeafbranch;
	double leafBranchSd;

	void fromNewick(const string& nString);

public:
	TripletSamplingTree(GuideTree& gt);

	~TripletSamplingTree();

	//sample tripplets on a tree
	vector<array<unsigned int, 3> > sampleFromTree();

	//sample only based on the distance matrix
	vector<array<unsigned int, 3> > sampleFromDM();

	// node operator for comparision
};

} /* namespace EBC */

#endif /* TRIPLETSAMPLINGTREE_HPP_ */
