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


#ifndef COMMANDREADER_H_
#define COMMANDREADER_H_

#include "core/IParser.hpp"
#include "core/HmmException.hpp"
#include "core/Definitions.hpp"
#include <dlib/cmd_line_parser.h>

namespace EBC
{

class CommandReader
{
public:

	dlib::command_line_parser parser;

	CommandReader(int argc, char** argv);
	IParser* getParser() throw (HmmException&);


	vector<double> getIndelParams();

	vector<double> getSubstParams();

	Definitions::ModelType getModelType()
	{
		if (parser.option("GTR"))
		{
			return Definitions::ModelType::GTR;
		}
		if (parser.option("HKY"))
		{
			return Definitions::ModelType::HKY85;
		}
		if (parser.option("LG"))
		{
			return Definitions::ModelType::LG;
		}
		if (parser.option("JTT"))
		{
			return Definitions::ModelType::JTT;
		}
		if (parser.option("WAG"))
		{
			return Definitions::ModelType::WAG;
		}
		//default;
		return Definitions::ModelType::GTR;
	}

	string getInputFileName()
	{
		return parser.option("in").argument();
	}

	Definitions::OptimizationType getOptimizationType()
	{
			return (Definitions::OptimizationType::BFGS);
	}

	FileLogger::logType getLoggingLevel()
	{
		if (parser.option("lDD"))
			return FileLogger::L_DMP;
		if (parser.option("lE"))
			return FileLogger::L_ERR;
		if (parser.option("lW"))
			return FileLogger::L_WARN;
		if (parser.option("lI"))
			return FileLogger::L_INF;
		if (parser.option("lD"))
			return FileLogger::L_DBG;
		//info by default!
		return
			FileLogger::L_INF;
	}

	double getAlpha()
	{
		return get_option(parser,"initAlpha",0.5);
	}

	double getCategories()
	{
		return get_option(parser,"rateCat",4);
	}

	bool estimateAlpha()
	{
		int res = get_option(parser,"estimateAlpha",1);
		return res == 1;
	}
	Definitions::SequenceType getSequenceType()
	{
		if (parser.option("GTR") || parser.option("HKY"))
			return Definitions::SequenceType::Nucleotide;
		else return Definitions::SequenceType::Aminoacid;
	}
};

} /* namespace EBC */
#endif /* COMMANDREADER_H_ */
