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


#include <heuristics/Band.hpp>

namespace EBC
{

Band::Band(unsigned int size) : matchBand(size), insertBand(size), deleteBand(size)
{
}

Band::Band(unsigned int len1, unsigned int len2, double coverage) :  matchBand(len2+1), insertBand(len2+1), deleteBand(len2+1)
{
	double ratio = static_cast<double>(len1+1)/(len2+1);
	double idealWidth = coverage * (len1+1);
	int calcHalfSize = static_cast<unsigned int>(idealWidth/2);
	int halfSize = calcHalfSize < Definitions::minBandDelta ? Definitions::minBandDelta : calcHalfSize;
	int estDiagonal;
	int min, max;

	//deal with the first column separately
	this->setMatchRangeAt(0,-1,-1);
	this->setInsertRangeAt(0,0,halfSize);
	this->setDeleteRangeAt(0,-1,-1);

	for (int col = 1; col <= len2; col++)
	{
		estDiagonal = static_cast<unsigned int>(col*ratio);

		min = (estDiagonal-halfSize) < 0 ? 0 : estDiagonal-halfSize;
		max = (estDiagonal+halfSize) > len1?  len1 : estDiagonal+halfSize;

		this->setMatchRangeAt(col,min+1,max);
		this->setInsertRangeAt(col,min+1,max);
		this->setDeleteRangeAt(col,min,max);
	}
}

Band::~Band()
{
}

} /* namespace EBC */


