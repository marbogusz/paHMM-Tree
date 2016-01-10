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


#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{
	class Dictionary
	{
	protected:
		unsigned short alphabetSize;
		unsigned char gapId;
		string alphabet;
		map<char,SequenceElement*> translator;

	public:

		static const char nucleotides[5];
		static const char aminoacids[21];
		static const char gapChar;
		static const map<char,vector<char> > nucFastaClasses;
		static const map<char,vector<char> > aaFastaClasses;

		//virtual short getSymbolIndex(string &symbol);
		virtual unsigned char getSymbolIndex(char symbol);
		virtual vector<SequenceElement*>* translate(string &sequence, bool removeGaps = false);
		virtual unsigned short getAlphabetSize();

		SequenceElement* getSequenceElement(char symbol);

		virtual char getSymbolAt(unsigned char i);

		virtual void outputAlphabet();

		inline unsigned char getGapID()
		{
			return gapId;
		}

	protected:
		virtual void setAlphabet(char alphabet[], unsigned short size);
		void addFastaClasses(const map<char,vector<char> >& classmap);

	};

	class NucleotideDictionary : public Dictionary
	{
	public:
		NucleotideDictionary();
		void handleTUequivalence();
	};

	class AminoacidDictionary : public Dictionary
	{
	public:
		AminoacidDictionary();
	};
}


#endif /* DICTIONARY_H_ */
