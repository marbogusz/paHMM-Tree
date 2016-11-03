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


#include "heuristics/TripletSamplingTree.hpp"
//#include <regex>
#include <stack>
#include <random>
#include <iostream>
#include <cmath>
#include "core/Definitions.hpp"

namespace EBC
{

TripletSamplingTree::TripletSamplingTree(GuideTree& gt) : distMat(gt.getDistanceMatrix())
{
	this->averageLeafbranch = 0;
	this->leafBranchSd = 0;
	DEBUG("Creating TripletSamplingTree");

	idealTreeSize = 1.0;
	this->fromNewick(gt.getNewickTree());

}

TripletSamplingTree::~TripletSamplingTree()
{
	for (auto nod : nodes)
		delete nod;
}

Node* TripletSamplingTree::mostRecentAncestor(Node* n1, Node* n2)
{
	Node* tmpNode1 = n1;
	Node* tmpNode2 = n2;

	while(tmpNode1 != NULL)
	{
		tmpNode2 = n2;
		while(tmpNode2 != NULL)
		{
			if (*tmpNode1 == *tmpNode2)
				return tmpNode1;
			tmpNode2 = tmpNode2->parent;
		}
	tmpNode1 = tmpNode1->parent;
	}
	return NULL;
}

double TripletSamplingTree::distanceToParent(Node* n1, Node* par)
{
	double distance = 0;
	Node* tmp = n1;

	if (*n1 == *par)
		return 0;

	while (tmp != par)
	{
		distance += tmp->distanceToParent;
		tmp = tmp->parent;
	}
	return distance;
}


//TODO - this is no longer in use really.
void TripletSamplingTree::fromNewick(const string& newick)
{

}

vector<array<unsigned int, 3> > TripletSamplingTree::sampleFromDM()
{
	vector<array<unsigned int, 3> > result;

	unsigned int s1, s2, s3;
	double remainingDistance = 2.0 * this->idealTreeSize;
	auto pair = distMat->getPairWithinDistance(0.5*this->idealTreeSize, 0.8*this->idealTreeSize);
	s1 = pair.first;
	s2 = pair.second;

	remainingDistance -= distMat->getDistance(s1, s2);

	s3 = distMat->getThirdLeafWithinDistance(remainingDistance, s1, s2);

	result.push_back({{s1,s2,s3}});
	DEBUG("Triplet tree DM sampled values : " << s1 << ", " << s2 << ", " << s3);
	cout << "Triplet tree DM sampled values : " << s1 << ", " << s2 << ", " << s3 << endl;
	return result;
}

vector<array<unsigned int, 3> > TripletSamplingTree::sampleFromTree()
{
	vector<array<unsigned int, 3> > result;

	//TODO - remove magic numbers
	pair<double,double> idealRange = make_pair(0.35,0.75);
	pair<double,double> secondaryRange = make_pair(0.15,0.85);

	unsigned int treeNo = 1;
	unsigned int maxTrees = Definitions::maxSampledTriplets;


	for (unsigned int i = 0;  i < distMat->getSize(); i++)
	{
		leafNodes[i] = nullptr;
	}
	availableNodes = leafNodes;

	bool found = false;
	unsigned int s1,s2;

	double tmpd1, tmpd2;

	auto vecPairsDesired = distMat->getPairsWithinDistance(idealRange.first, idealRange.second);
	auto vecPairsSecondary = distMat->getPairsWithinDistance(secondaryRange.first, secondaryRange.second);

	if (vecPairsDesired.size() == 0 && vecPairsSecondary.size() == 0){
		vecPairsSecondary.push_back(distMat->getPairWithinDistance(idealRange.first, idealRange.second));
		DUMP("Triplet sampling tree : no pairs within desired distance found, geting the closest one");
	}
	//Check ideal range
	for(auto pr : vecPairsDesired){
		if (treeNo > maxTrees)
			break;
		s1 = pr.first;
		s2 = pr.second;

		if (availableNodes.find(s1) == availableNodes.end()
				|| availableNodes.find(s2) == availableNodes.end())
			continue;

		for (auto nd : availableNodes)
		{

			if (nd.first == s1 || nd.first == s2)
				continue;
			tmpd1 =  distMat->getDistance(nd.first, s1);
			tmpd2 =  distMat->getDistance(nd.first, s2);

			if (isWithinRange(tmpd1, idealRange) && isWithinRange(tmpd2, idealRange)) {
				result.push_back({{s2,s1,nd.first}});
				DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << nd.first);
				availableNodes.erase(s1);
				availableNodes.erase(s2);
				availableNodes.erase(nd.first);
				distMat->invalidate(pr);
				found = true;
				treeNo++;
				break;
			}
			else{
				continue;
			}
		}
	}
	//not found
	if(!found || treeNo < maxTrees){
		for(auto pr : vecPairsSecondary){
			if (treeNo > maxTrees)
				break;
			s1 = pr.first;
			s2 = pr.second;

			if (availableNodes.find(s1) == availableNodes.end()
					|| availableNodes.find(s2) == availableNodes.end())
				continue;

			for (auto nd : availableNodes)
			{

				if (nd.first == s1 || nd.first == s2)
					continue;
				tmpd1 =  distMat->getDistance(nd.first, s1);
				tmpd2 =  distMat->getDistance(nd.first, s2);

				if (isWithinRange(tmpd1, secondaryRange) && isWithinRange(tmpd2, secondaryRange)) {
					result.push_back({{s2,s1,nd.first}});
					DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << nd.first);
					availableNodes.erase(s1);
					availableNodes.erase(s2);
					availableNodes.erase(nd.first);
					distMat->invalidate(pr);
					found = true;
					treeNo++;
					break;
				}
				else{
					continue;
				}
			}
		}
	}
	//Still not found?
	if(!found){
		for(auto pr : vecPairsDesired){
			if (treeNo >maxTrees)
				break;
			s1 = pr.first;
			s2 = pr.second;

			if (availableNodes.find(s1) == availableNodes.end()
					|| availableNodes.find(s2) == availableNodes.end())
				continue;

			for (auto nd : availableNodes)
			{

				if (nd.first == s1 || nd.first == s2)
					continue;
				tmpd1 =  distMat->getDistance(nd.first, s1);
				tmpd2 =  distMat->getDistance(nd.first, s2);

				if (isWithinRange(tmpd1, idealRange)) {
					result.push_back({{s2,s1,nd.first}});
					DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << nd.first);
					availableNodes.erase(s1);
					availableNodes.erase(s2);
					availableNodes.erase(nd.first);
					distMat->invalidate(pr);
					found = true;
					treeNo++;
					break;
				}
				else if(isWithinRange(tmpd2, idealRange)){
					result.push_back({{s1,s2,nd.first}});
					DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << nd.first);
					availableNodes.erase(s1);
					availableNodes.erase(s2);
					availableNodes.erase(nd.first);
					distMat->invalidate(pr);
					found = true;
					treeNo++;
					break;
				}
				else{
					continue;
				}
			}
		}
	}
	//still nothing
	if(!found){
		for(auto pr : vecPairsSecondary){
			if (treeNo >maxTrees)
				break;
			s1 = pr.first;
			s2 = pr.second;

			if (availableNodes.find(s1) == availableNodes.end()
					|| availableNodes.find(s2) == availableNodes.end())
				continue;


			for (auto nd : availableNodes)
			{

				if (nd.first == s1 || nd.first == s2)
					continue;
				tmpd1 =  distMat->getDistance(nd.first, s1);
				tmpd2 =  distMat->getDistance(nd.first, s2);

				if (isWithinRange(tmpd1, secondaryRange)) {
					result.push_back({{s2,s1,nd.first}});
					DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << nd.first);
					availableNodes.erase(s1);
					availableNodes.erase(s2);
					availableNodes.erase(nd.first);
					distMat->invalidate(pr);
					found = true;
					treeNo++;
					break;
				}
				else if(isWithinRange(tmpd2, secondaryRange)){
					result.push_back({{s1,s2,nd.first}});
					DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << nd.first);
					availableNodes.erase(s1);
					availableNodes.erase(s2);
					availableNodes.erase(nd.first);
					distMat->invalidate(pr);
					found = true;
					treeNo++;
					break;
				}
				else{
					continue;
				}
			}
		}
	}
	//Nothing within those ranges - just get something
	if(!found){
		//check secondary range
		auto pr = vecPairsSecondary[0]; //get the best one
		s1 = pr.first;
		s2 = pr.second;

		for (auto nd : availableNodes)
		{

			if (nd.first == s1 || nd.first == s2)
				continue;
			tmpd1 =  distMat->getDistance(nd.first, s1);
			tmpd2 =  distMat->getDistance(nd.first, s2);


			if (isWithinRange(tmpd1, secondaryRange)) {
				result.push_back({{s2,s1,nd.first}});
				DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << nd.first);
				availableNodes.erase(s1);
				availableNodes.erase(s2);
				availableNodes.erase(nd.first);
				distMat->invalidate(pr);
				found = true;
				treeNo++;
				break;
			}
			else if(isWithinRange(tmpd2, secondaryRange)){
				result.push_back({{s1,s2,nd.first}});
				DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << nd.first);
				availableNodes.erase(s1);
				availableNodes.erase(s2);
				availableNodes.erase(nd.first);
				distMat->invalidate(pr);
				found = true;
				treeNo++;
				break;
			}
			else{
				continue;
			}
		}
		if(!found){
			//nothing was found. this means that either all the distances are very small or big
			//small = get the biggest, big - get the smallest
			unsigned int id = 0;
			double tdist = 0;
			//bool largeDistances = false;

			//if (distMat->getDistance(s1, s2) > secondaryRange.second){
			//	largeDistances = true;
			//}
			double currentDist = 100;

			availableNodes.erase(s1);
			availableNodes.erase(s2);

			for(auto nd : availableNodes)
			{
				//deviation from the sweetspot
				//TODO - magicnumber for the sweetspot
				tdist = abs(0.5 - distMat->getDistance(nd.first, s1)) + abs(0.5 - distMat->getDistance(nd.first, s2));

					if (tdist < currentDist){
						currentDist = tdist;
						id = nd.first;
					}

			}

			tmpd1 = distMat->getDistance(s1, id);
			tmpd2 = distMat->getDistance(s2, id);

			if (tmpd1 > distMat->getDistance(s1,s2)){
				if(tmpd1 < tmpd2){
					result.push_back({{s2,s1,id}});
					DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << id);
				}
				else{
					result.push_back({{s1,s2,id}});
					DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << id);
				}
			}
			else{
				if(tmpd1 < tmpd2){
					result.push_back({{s1,s2,id}});
					DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << id);
				}
				else{
					result.push_back({{s2,s1,id}});
					DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << id);
				}
			}

//			result.push_back({{s1,s2,id}});
//			DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << id);
			treeNo++;
		}

	}
	INFO(treeNo << " triplets sampled");
	return result;
}


} /* namespace EBC */
