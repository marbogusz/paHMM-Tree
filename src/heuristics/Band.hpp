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


#ifndef HEURISTICS_BAND_HPP_
#define HEURISTICS_BAND_HPP_

#include <vector>
#include "core/FileLogger.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{

class Band
{
protected:

	vector<pair<int, int> > matchBand;
	vector<pair<int, int> > insertBand;
	vector<pair<int, int> > deleteBand;


public:
	Band(unsigned int size);

	//Creates a default band that covers specified fraction of the column
	Band(unsigned int len1, unsigned int len2, double coverage=Definitions::initialBandFactor);

	virtual ~Band();

	inline void setMatchRangeAt(unsigned int pos, int start, int end)
	{
		matchBand[pos] = std::make_pair(start,end);
	}
	inline void setInsertRangeAt(unsigned int pos, int start, int end)
	{
		insertBand[pos] = std::make_pair(start,end);
	}
	inline void setDeleteRangeAt(unsigned int pos, int start, int end)
	{
		deleteBand[pos] = std::make_pair(start,end);
	}

	inline std::pair<int, int> getMatchRangeAt(unsigned int pos)
	{
		return matchBand[pos];
	}
	inline std::pair<int, int> getInsertRangeAt(unsigned int pos)
	{
		return insertBand[pos];
	}
	inline std::pair<int, int> getDeleteRangeAt(unsigned int pos)
	{
		return deleteBand[pos];
	}

	inline void output()
	{
		for(unsigned int i =0; i< matchBand.size(); i++)
		{
			DUMP("M/X/Y bands " << i << "\t" << getMatchRangeAt(i).first <<"\t" << getMatchRangeAt(i).second
							<< "\t" << getInsertRangeAt(i).first <<"\t" << getInsertRangeAt(i).second
							<< "\t" << getDeleteRangeAt(i).first <<"\t" << getDeleteRangeAt(i).second);

		}
	}

	const vector<pair<int, int> >& getDeleteBand() const {
		return deleteBand;
	}

	const vector<pair<int, int> >& getInsertBand() const {
		return insertBand;
	}

	const vector<pair<int, int> >& getMatchBand() const {
		return matchBand;
	}
};

} /* namespace EBC */

#endif /* HEURISTICS_BAND_HPP_ */
