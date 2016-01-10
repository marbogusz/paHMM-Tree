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

#include "core/DistanceMatrix.hpp"
#include "core/Definitions.hpp"
#include <random>

namespace EBC
{


EBC::DistanceMatrix::DistanceMatrix(int size) :taxas(size)
{

}

double EBC::DistanceMatrix::getDistance(unsigned int i, unsigned int j)
{
	double dst = 0;

	dst = this->distances[make_pair(i,j)];

	//DEBUG("Distance matrix getting distance");

	return dst;
}

void DistanceMatrix::buildMap()
{
}

void DistanceMatrix::addDistance(unsigned int s1, unsigned int s2,
		double distance)
{
	distances[make_pair(s1,s2)] = distance;
	distances[make_pair(s2,s1)] = distance;

	revdistances.insert(make_pair(distance,(make_pair(s1,s2))));
}

void  DistanceMatrix::invalidate(std::pair<unsigned int, unsigned int>& pr)
{
	for(auto el = revdistances.begin(); el != revdistances.end(); el++)
	{
		if (el->second.first == pr.first || el->second.first == pr.second ||el->second.second == pr.first || el->second.second == pr.second)
			revdistances.erase(el);
	}
}

vector<pair<unsigned int, unsigned int> > EBC::DistanceMatrix::getPairsWithinDistance(
		double lo, double hi)
{
	vector<pair<unsigned int, unsigned int> > retvec;
	auto itlow = revdistances.lower_bound(lo);
	auto ithi = revdistances.upper_bound(hi);
	auto end = revdistances.end();
	auto begin = revdistances.begin();

	if (itlow != end && ithi != begin){
		//we're in the bracket
		//unsigned int dist = std::distance(itlow,ithi);
		//uniform_int_distribution<int> dist2(0,dist);
		//std::advance(itlow, dist2(generator));
		while(itlow != ithi){
			retvec.push_back(make_pair((*itlow).second.first,(*itlow).second.second));
			itlow++;
		}
		//invalidate(ret);
	}
	return retvec;
}

pair<unsigned int, unsigned int> EBC::DistanceMatrix::getPairWithinDistance(
		double lo, double hi)
{
	//default_random_engine generator;
	//uniform_int_distribution<int> distribution(0,revdistances.size()-1);
	//iterate over keys
	auto itlow = revdistances.lower_bound(lo);
	auto ithi = revdistances.upper_bound(hi);
	auto end = revdistances.end();
	auto begin = revdistances.begin();
	pair<unsigned int, unsigned int> ret;

	if (itlow != end && ithi != begin){
		//we're in the bracket
		//unsigned int dist = std::distance(itlow,ithi);
		//uniform_int_distribution<int> dist2(0,dist);
		//std::advance(itlow, dist2(generator));
		ret = make_pair((*itlow).second.first,(*itlow).second.second);
		invalidate(ret);
	}
	//return the lowest anyway ?
	else if(ithi == begin){
		//std::advance(begin, distribution(generator) );
		ret = make_pair((*begin).second.first,(*begin).second.second);
		invalidate(ret);
	}
	else {
		ret = make_pair((*(--end)).second.first,(*(--end)).second.second);
		invalidate(ret);
	}
	return ret;
}

unsigned int DistanceMatrix::getThirdLeafWithinDistance(double targetLen, unsigned int l1, unsigned int l2)
{
	unsigned int leaves = 0;
	double tempSum;
	map<double, unsigned int> mappings;
	for (leaves =0; leaves < this->taxas; leaves++)
	{
		if(leaves == l1 || leaves ==l2)
			continue;
		tempSum = getDistance(l1,leaves) + getDistance(l2,leaves);
		mappings[tempSum] = leaves;
	}
	//range iterator over mappings!
	//find the closest to target len
	auto itlow= mappings.lower_bound(0.7*targetLen);
	auto ithi= mappings.upper_bound(1.3*targetLen);
	auto end = mappings.end();

	if (itlow == end)
		return (*(end--)).second;
	else if (ithi == end)
		return (*(mappings.begin())).second;
	else return (*itlow).second;
}

} /* namespace EBC */
