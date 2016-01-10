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


#include <vector>
#include <string>

#ifndef NODE_HPP_
#define NODE_HPP_

using namespace std;

namespace EBC
{

struct Node
{
public:

	unsigned int nodeId;

	Node* parent;

	bool leafNode;

	bool rootNode;

	unsigned int sequenceNo;

	string nodeName;

	double distanceToParent;

	friend bool operator== (Node &cP1, Node &cP2)
	{
		return cP1.nodeId == cP2.nodeId;
	}

	vector<Node*> children;

	Node(unsigned int id);

	void setName(string name);

	inline void setSequenceNumber(unsigned int no)
	{
		this->sequenceNo = no;
	}

	void setDistance(double distance);

	void setLeaf();

	void setParent(Node* pn);

	void setChild(Node* cn);

};

} /* namespace EBC */

#endif /* NODE_HPP_ */
