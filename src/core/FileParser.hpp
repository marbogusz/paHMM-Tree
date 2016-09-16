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


#ifndef FILEPARSER_H_
#define FILEPARSER_H_

#include <fstream>
#include "core/HmmException.hpp"
#include "core/IParser.hpp"
#include <vector>
#include <map>

using namespace std;

namespace EBC
{

class FileParser : public IParser
{
private:

	string filename;
	ifstream infile;
	vector<string>* sequences;
	vector<string>* names;
	map<string,string> mappedSeqs;
	vector<string>::iterator it;
	vector<string>::iterator itN;

public:

	FileParser(const char* filename);

	string getNextSequence();

	string getNextName();

	unsigned int getSequenceCount();

	string getSequenceAt(unsigned int);

	string getSequenceNameAt(unsigned int);

	bool isDefinitionLine(string&);

	string getSequenceName(string&);

	void trimWsChars(string&);

	virtual ~FileParser();

	inline vector<string>* getNames() {
		return names;
	}

	inline vector<string>* getSequences() {
		return sequences;
	}
};

} /* namespace EBC */
#endif /* FILEPARSER_H_ */
