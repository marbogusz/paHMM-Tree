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


#ifndef SEQUENCEELEMENT_H_
#define SEQUENCEELEMENT_H_

#include <vector>
#include <string>
using namespace std;


namespace EBC
{

class SequenceElement
{

protected:
	bool isGap;
	unsigned char matrixIndex;
	char symbol;
	unsigned short clSize;
	bool fastaClass;
	unsigned char* idcs;
public:
	SequenceElement(bool isGap, unsigned char substMatrixIndex, unsigned char* altIdcs, char smbl, unsigned short classSize=1);

	inline bool isIsGap() const
	{
		return isGap;
	}

	inline bool isFastaClass() const
	{
		return fastaClass;
	}

	inline unsigned short getClassSize()
	{
		return clSize;
	}

	inline unsigned char getMatrixIndex()
	{
		return matrixIndex;
	}

	inline char getSymbol()
	{
		return symbol;
	}

	inline unsigned char* getClassIndices()
	{
		return idcs;
	}
};

} /* namespace EBC */
#endif /* SEQUENCEELEMENT_H_ */
