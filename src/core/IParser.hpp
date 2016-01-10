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


#ifndef IPARSER_H_
#define IPARSER_H_

#include "core/Definitions.hpp"
#include <string>
#include <vector>

using namespace std;

namespace EBC
{

class IParser
{
public:
	virtual string getNextName() = 0;
	virtual string getNextSequence() = 0;
	virtual unsigned int getSequenceCount() = 0;
	virtual string getSequenceAt(unsigned int) = 0;
	virtual vector<string>* getNames() =0;
	virtual vector<string>* getSequences() =0;

};

} /* namespace EBC */
#endif /* IPARSER_H_ */
