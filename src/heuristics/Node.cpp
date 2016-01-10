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


#include "heuristics/Node.hpp"
#include "core/Definitions.hpp"

namespace EBC
{

Node::Node(unsigned int id) : nodeId(id), parent(NULL), leafNode(false), rootNode(false), distanceToParent(0)
{
	DUMP("New node with id " << id);
}

void Node::setName(string name)
{
	this->nodeName = name;
}

void Node::setDistance(double distance)
{
	this->distanceToParent = distance;
}

void Node::setLeaf()
{
	this->leafNode = true;
}

void Node::setParent(Node* pn)
{
	this->parent = pn;
}

void Node::setChild(Node* cn)
{
	this->children.insert(children.end(), cn->children.begin(), cn->children.end());
	this->children.push_back(cn);
}

} /* namespace EBC */
