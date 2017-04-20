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



#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <iostream>
#include "core/FileLogger.hpp"

//#define DEBUG_BUILD 1

#define DUMP(x) do { FileLogger::DumpLogger() << "   [DUMP]\t" << x << "\n"; } while (0)
//#define DUMP(x) do {} while (0)
#define DEBUG(x) do { FileLogger::DebugLogger() << "  [DEBUG]\t" << x << "\n"; } while (0)
#define INFO(x) do { FileLogger::InfoLogger() << " [INFO]\t"  << x<< "\n"; } while (0)
#define WARN(x) do { FileLogger::WarningLogger() << "! [WARNING]\t"  << x << "\n"; } while (0)
#define ERROR(x) do { FileLogger::ErrorLogger() << "!!! [ERROR]\t"  << x << "\n"; } while (0)

#  define DEBUGN(x) do {} while (0)
#  define DEBUGV(x,n) do {} while (0)
#  define DEBUGM(x,n) do {} while (0)

/*
#ifdef DEBUG_BUILD
#  define DEBUG(x) do { std::cerr<< "DBG! " << x << std::endl; } while (0)
#  define DEBUGN(x) do { std::cerr << x ; } while (0)
#  define DEBUGV(x,n) do { for (int i=0;i<n;i++) std::cerr << x[i] << "\t"; std::cerr << std::endl;} while (0)
#  define DEBUGM(x,h,w) do { for (int i=0;i<h;i++){std::cerr << endl; for(int j=0;j<w;j++) std::cerr << x[i][j] << "\t";} std::cerr <<std::endl;} while (0)
#else
#  define DEBUG(x) do {} while (0)
#  define DEBUGN(x) do {} while (0)
#  define DEBUGV(x,n) do {} while (0)
#  define DEBUGM(x,n) do {} while (0)
#endif
*/
namespace EBC
{



class Definitions
{
public:


	//this value coming from k-mer estimation is considered high and means that
	//this value coming from k-mer estimation is considered low and means that
	//the actual distance is likely to be < 1
	constexpr static const double kmerLowDivergence = 0.6;
	constexpr static const double kmerHighDivergence = 0.8;

	constexpr static const int maxSampledTriplets = 5;

	//for band calculations
	constexpr static const double bandPosteriorLikelihoodLimit = -3;
	constexpr static const double bandPosteriorLikelihoodDelta = -9;

	constexpr static const double normalDivergenceAccuracyDelta = 1e-3;

	constexpr static const double highDivergenceAccuracyDelta = 1e-2;

	constexpr static const double ultraDivergenceAccuracyDelta = 1e-1;


	//model param estimation accuracy
	constexpr static const double accuracyBFGS = 1e-8;

	constexpr static const int BrentMaxIter = 100;

	//band factor default for intial fwd likelihood calculations
	constexpr static const double narrowBandFactor = 0.1;
	constexpr static const double initialBandFactor = 0.33;

	//small number close to zero for param estimation
	constexpr static const double almostZero = 1e-8;

	constexpr static const double minModelBound = 0.000001;

	constexpr static const double maxAlpha = 99.999999;

	//3 states - M I D
	constexpr static const unsigned int stateCount = 3;

	//triplet trees
	constexpr static const unsigned int heuristicsTreeSize = 3;

	//obsolete
	constexpr static const unsigned int pathSampleCount = 10000;
	constexpr static const unsigned int pathInformativeCount = 50;


	//NEW sampling
	constexpr static const unsigned int samplingPathMaxCount = 1000000;
	constexpr static const unsigned int samplingPathMinCount = 50000;
	constexpr static const unsigned int samplingPathCount = 1000;
	constexpr static const double samplingPathLnLDelta = 10.0;

	//max divergence
	constexpr static const double divergenceBound = 50;

	//initial max lambda, however the max depends on divergence
	constexpr static const double lambdaHiBound = 0.3;

	constexpr static const double gtrMaxRate = 5;

	constexpr static const double epsilonHiBound = 0.95;


	constexpr static const double initialLambda = 0.05;
	constexpr static const double initialEpsilon = 0.5;

	//This makes the min band width of 15 characters
	constexpr static const unsigned int minBandDelta = 7;

	constexpr static const double minMatrixLikelihood = -1000000.0;


	constexpr static const unsigned int HKY85ParamCount = 1;
	constexpr static const unsigned int GTRParamCount = 5;
	constexpr static const unsigned int AAParamCount = 0;
	constexpr static const unsigned int NBIndelParamCount = 2;

	constexpr static const unsigned int nucleotideCount = 4;
	constexpr static const unsigned int aminoacidCount = 20;

	constexpr static auto distMatExt = ".paHMM-Tree.distmat";
	constexpr static auto paramsExt = ".paHMM-Tree.params";
	constexpr static auto treeExt = ".paHMM-Tree.tree";
	constexpr static auto logExt = ".paHMM-Tree.log";


	struct aaModelDefinition
	{
		double maxRate;
		double aaRates[400];
		double aaFreqs[20];
	};

	enum SequenceType {Nucleotide, Aminoacid, Codon};
	//TODO - CODON lookup table
	enum ModelType {GTR, HKY85, LG};

	enum OptimizationType {BFGS, BOBYQA};

	enum AlgorithmType {Forward, Viterbi, MLE};

	enum DpMatrixType {Full, Limited};

	enum StateId {Match, Insert , Delete};

	static aaModelDefinition aaLgModel;


constexpr static const auto notice ="\npaHMM-Tree - distance-based statistical phylogenetic tree estimation version 0.161101\n\
Copyright (C) 2016  Marcin Bogusz\n\
This program is free software: you can redistribute it and/or modify\n\
it under the terms of the GNU General Public License as published by\n\
the Free Software Foundation, either version 3 of the License, or\n\
(at your option) any later version.\n\
\n\
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;\n\
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\
See the GNU General Public License for more details.\n\
\n\
You should have received a copy of the GNU General Public License\n\
along with this program.  If not, see <http://www.gnu.org/licenses/>.\n\n";

};

} /* namespace EBC */
#endif /* DEFINITIONS_H_ */
