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

#ifndef TREE_HPP_
#define TREE_HPP_

#include <vector>
#include <array>
#include <map>
#include <unordered_set>
#include "heuristics/Node.hpp"
#include "core/Sequences.hpp"

using namespace std;

namespace EBC
{

class PhylogeneticTree
{
private:
	//all nodes
	vector<Node*> nodes;

	Sequences* inputSeqs;

	map<unsigned int, Node*> leafNodes;

	map<unsigned int, Node*> availableNodes;

	unordered_set<Node*> usedNodes;

	Node* mostRecentAncestor(Node* n1, Node* n2);

	double distanceToParent(Node* n1, Node* par);

public:
	PhylogeneticTree(Sequences* seqs);

	~PhylogeneticTree();

	void fromNewick(const string& nString);

	double distanceById(unsigned int n1, unsigned int n2);
	double distanceByName(string& n1, string& n2);

};

} /* namespace EBC */

#endif /* TRIPLETSAMPLINGTREE_HPP_ */
