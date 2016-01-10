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

#ifndef DISTANCEMATRIX_HPP_
#define DISTANCEMATRIX_HPP_

#include <map>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

namespace EBC
{

class DistanceMatrix
{
private:

	//map the length to the pair of sequences
	multimap <double, pair<unsigned int, unsigned int> > revdistances;

	//pairwise distances ordered, double the size of dictionary!
	map<pair<unsigned int, unsigned int>,double> distances;

	unsigned int taxas;

	void buildMap();

public:
	DistanceMatrix(int size);

	unsigned int getSize()
	{
		return taxas;
	}

	void addDistance(unsigned int s1, unsigned int s2, double distance);

	double getDistance(unsigned int s1, unsigned int s2);

	pair<unsigned int, unsigned int> getPairWithinDistance(double lo, double hi);

	unsigned int getThirdLeafWithinDistance(double targetLen, unsigned int l1, unsigned int l2);

	vector<pair<unsigned int, unsigned int> > getPairsWithinDistance(double lo, double hi);

	void invalidate(std::pair<unsigned int, unsigned int>&);

};

} /* namespace EBC */

#endif /* DISTANCEMATRIX_HPP_ */
