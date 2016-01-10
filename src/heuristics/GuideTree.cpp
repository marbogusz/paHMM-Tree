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


#include "heuristics/GuideTree.hpp"
#include "heuristics/TripletSamplingTree.hpp"

namespace EBC
{

GuideTree::GuideTree(Sequences* is) : inputSequences(is)
{
	distMat = new DistanceMatrix(inputSequences->getSequenceCount());
	this->dict = inputSequences->getDictionary();
	if(dict->getAlphabetSize() == Definitions::nucleotideCount)
		this->kmerSize = 7;
	else
		this->kmerSize = 4;
	this->sequenceCount = inputSequences->getSequenceCount();
	this->kmers = new vector<unordered_map<string,short>*>(sequenceCount);
	DEBUG("Creating guide tree");
	this->constructTree();
}

GuideTree::~GuideTree()
{
	delete kmers;
}

double GuideTree::kimuraDist(double id)
{
	//TODO delete this method
	return 0;
}

void GuideTree::constructTree()
{
	unsigned int i,j;
	string currSeq;
	double identity, estIdentity;



	DEBUG("Extracting k-mers");
	for(i = 0; i< sequenceCount; i++)
	{
		(*kmers)[i] = new unordered_map<string,short>();
		currSeq = inputSequences->getRawSequenceAt(i);
		extractKmers(currSeq, (*kmers)[i]);
	}
	for(i = 0; i< sequenceCount; i++)
		for(j = i+1; j< sequenceCount; j++)
		{
			string s1 = inputSequences->getRawSequenceAt(i);
			string s2 = inputSequences->getRawSequenceAt(j);
			identity = 1.0 - commonKmerCount(i,j)/((double)(min(s1.size(),s2.size())));

			if(dict->getAlphabetSize() == Definitions::nucleotideCount)
				estIdentity = this->nucFunction(identity);
			else if(dict->getAlphabetSize() == Definitions::aminoacidCount)
				estIdentity = aaFunction(identity);

			DEBUG("k-mer distance between seq. " << i << " and " << j << " is " << identity << " adjusted distance " << estIdentity );

			distMat->addDistance(i,j,estIdentity);
			distances.push_back(estIdentity);
		}

	DEBUG("Initialized distance matrix");
	//BioNJ nj(sequenceCount, distMat);
	//newickTree = nj.calculate();
	newickTree = "";

}

void GuideTree::extractKmers(string& seq, unordered_map<string, short>* umap)
{
	string kmer;

	if(seq.size() < kmerSize)
		return;

	for(unsigned int i = 0; i< (seq.size() - kmerSize); i++)
	{
		kmer = seq.substr(i, kmerSize);
		++((*umap)[kmer]);
	}
}

unsigned int GuideTree::commonKmerCount(unsigned int i, unsigned int j)
{
	unordered_map<string, short>* m1 = (*kmers)[i];
	unordered_map<string, short>* m2 = (*kmers)[j];
	unsigned int commonCount = 0;

	for(auto it = m1->begin(); it != m1->end(); it++)
	{
		commonCount += std::min((short)((*m2)[it->first]), (short)(it->second));
	}
	//DEBUG("Common k-mer count between seq. " << i << " and " << j << " is " << commonCount);
	return commonCount;


}

} /* namespace EBC */
