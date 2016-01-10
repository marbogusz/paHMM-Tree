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


#ifndef SEQUENCES_H_
#define SEQUENCES_H_

#include <vector>
#include "core/Definitions.hpp"
#include "core/IParser.hpp"
#include "core/Dictionary.hpp"
#include "core/HmmException.hpp"
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"
#include <array>
#include <vector>


using namespace std;

namespace EBC
{

class Sequences
{
private:

	vector<string>* rawSequences;

	vector<string>* sequenceNames;

	vector<vector<SequenceElement*>* > translatedSequences;
	vector<std::pair<unsigned int, unsigned int> > pairs;
	vector<std::pair<unsigned int, unsigned int> >::iterator pairIterator;

	Dictionary* dict;
	unsigned int sequenceCount;
	double* observedFrequencies;

	bool removeGaps;

public:

	//Input from file or console
	Sequences(IParser*, Definitions::SequenceType, bool fixedAlignment=false) throw (HmmException&);

	~Sequences();

	//Return the dictionary for the input set
	Dictionary* getDictionary();

	//Get then in a dictionary order eg. T C A G, check the definitions in constants
	//TODO - change to A T C G
	double* getElementFrequencies();

	double* getElementFrequencies(array<unsigned int, 3>& triplet);

	//void getSequencePair(vector<SequenceElement> s1, vector<SequenceElement> s2 );
	vector<SequenceElement*>* getSequencesAt(unsigned int pos);

	inline unsigned int getPairCount()
	{
		unsigned int ct = translatedSequences.size();
		return (ct*(ct-1))/2;
	}

	inline unsigned int getSequenceCount()
	{
		return translatedSequences.size();
	}

	string& getSequenceName(unsigned int pos);

	unsigned int getSequenceId(string& seqname);

	string& getRawSequenceAt(unsigned int pos);

	std::pair<unsigned int, unsigned int> getPairOfSequenceIndices(unsigned int idx)
	{
		return pairs[idx];
	}
private:

	void calculateObservedFrequencies();

	void buildDictionary(Definitions::SequenceType);

};

} /* namespace EBC */



#endif /* SEQUENCES_H_ */
