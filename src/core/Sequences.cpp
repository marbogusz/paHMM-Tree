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


#include "core/Sequences.hpp"
#include <algorithm>

namespace EBC
{

Sequences::Sequences(IParser* iParser,Definitions::SequenceType st, bool rg) throw (HmmException&)
{
	//use the file parser to get sequences and build the dictionary
	removeGaps = rg;
	observedFrequencies = NULL;

	unsigned int size = iParser->getSequenceCount();
	if (size <= 0){
		throw HmmException("No FASTA sequences found in the input file. Quitting...\n");
	}
	else if (size < 3){
		throw HmmException("paHMM-Tree requires at least 3 sequences to run. Quitting...\n");
	}

	this->buildDictionary(st);

	this->sequenceCount = size;

	pairs.reserve(this->getPairCount());

	for(unsigned int i=0; i< size;i++)
		for(unsigned int j=i+1; j<size;j++)
			pairs.push_back(std::make_pair(i,j));

		pairIterator = pairs.begin();

	this->rawSequences = iParser->getSequences();
	this->sequenceNames = iParser->getNames();

	for (auto it = rawSequences->begin(); it != rawSequences->end(); it++)
		this->translatedSequences.push_back(dict->translate(*it,removeGaps));

}

vector<SequenceElement*>* Sequences::getSequencesAt(unsigned int pos){
		return translatedSequences[pos];
}

Sequences::~Sequences()
{
	delete dict;
	delete[] observedFrequencies;
}

Dictionary* Sequences::getDictionary()
{
	return dict;
}


string& Sequences::getRawSequenceAt(unsigned int pos)
{
	return (*rawSequences)[pos];
}

string& Sequences::getSequenceName(unsigned int pos)
{
	return (*sequenceNames)[pos];
}

unsigned int Sequences::getSequenceId(string& seqname)
{
	unsigned int pos = std::find(sequenceNames->begin(), sequenceNames->end(), seqname) - sequenceNames->begin();
	if(pos >= sequenceNames->size()) {
	    throw HmmException(seqname + " not found");
	}
	return pos;
}

void Sequences::calculateObservedFrequencies()
{
	this->observedFrequencies = new double[dict->getAlphabetSize()];

	unsigned int i;

	for (i=0; i<dict->getAlphabetSize(); i++)
		this->observedFrequencies[i] = 0;

	unsigned int count =0;

	for (auto it1 = translatedSequences.begin() ; it1 != translatedSequences.end(); ++it1)
	{
		for(auto it2 = (*it1)->begin(); it2 != (*it1)->end(); ++it2)
		{
			if (!((*it2)->isIsGap()))
			{
				auto elcount = (*it2)->getClassSize();
				count += elcount;
				if(elcount > 1){
					auto ids = (*it2)->getClassIndices();
					while(elcount > 0){
						observedFrequencies[ids[elcount-1]]++;
						elcount--;
					}
				}
				else
					observedFrequencies[(*it2)->getMatrixIndex()]++;
			}
		}
	}

	for (i=0; i < dict->getAlphabetSize(); i++)
		this->observedFrequencies[i] /= count;
}

double* Sequences::getElementFrequencies()
{
	if(observedFrequencies == NULL)
		calculateObservedFrequencies();
	//DEBUGV(observedFrequencies,4);
	return observedFrequencies;
}

double* Sequences::getElementFrequencies(array<unsigned int, 3>& triplet)
{

	if(observedFrequencies == NULL)
		this->observedFrequencies = new double[dict->getAlphabetSize()];
	int i;

	for (i=0; i<dict->getAlphabetSize(); i++)
		this->observedFrequencies[i] = 0;

	unsigned int count =0;

	for (auto it1 : triplet)
	{
		for(auto it2 = (translatedSequences[it1])->begin(); it2 != (translatedSequences[it1])->end(); ++it2)
		{
			if (!((*it2)->isIsGap()))
			{
				count++;
				observedFrequencies[(*it2)->getMatrixIndex()]++;
			}
		}
	}

	for (i=0; i < dict->getAlphabetSize(); i++)
		this->observedFrequencies[i] /= count;

	//DEBUGV(observedFrequencies,4);
	return observedFrequencies;
}

void Sequences::buildDictionary(Definitions::SequenceType st)
{
	//TODO differentiate based on the contents of the file!!!
	switch(st)
	{
	case (Definitions::SequenceType::Aminoacid):
		dict = new AminoacidDictionary();
		break;

	case (Definitions::SequenceType::Nucleotide):
		dict = new NucleotideDictionary();
		break;
	}
}



} /* namespace EBC */


