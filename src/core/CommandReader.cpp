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


#include "core/CommandReader.hpp"
#include "core/FileParser.hpp"
#include <sstream>
#include <cstring>
#include <cstdlib>

using namespace std;

namespace EBC
{

CommandReader::CommandReader(int argc, char** argv)
{
	try
	{
		parser.add_option("in","This option takes one argument which specifies the name of the file we want to analyze",1);
		parser.add_option("GTR", "General time reversible substitution model");
		parser.add_option("HKY", "HKY85 substitution model");
		parser.add_option("LG", "Le & Gasquel AA substitution sodel");
		parser.add_option("indel_params","indel parameters (NB probability and rate)",2);
		parser.set_group_name("Miscellaneous Options");
		parser.add_option("h","Display this help message.");
		parser.add_option("gtr_params","GTR model parameters",5);
		parser.add_option("hky_params","HKY85 model parameters",1);
		parser.add_option("rateCat", "Specify gamma rate categories, default is 4",1);

		parser.add_option("initAlpha", "Specify initial discrete Gamma shape parameter alpha , default is 0.5",1 );

		parser.add_option("estimateAlpha", "Specify to estimate discrete Gamma shape parameter alpha 0|1, default is 1",1 );

		parser.add_option("lE", "log error");
		parser.add_option("lW", "log warning");
		parser.add_option("lI", "log info");
		parser.add_option("lD", "log debug");
		parser.add_option("lDD", "log dump");

		parser.parse(argc,argv);

		//const char* one_time_opts[] = {"in", "indel_params" ,"h","b"};
		//parser.check_one_time_options(one_time_opts);

		parser.check_incompatible_options("GTR", "HKY");
		//parser.check_incompatible_options("d", "F");

		const char* rev_sub_opts[] = {"gtr_params"};
		const char* hky_sub_opts[] = {"hky_params"};
		parser.check_sub_options("GTR", rev_sub_opts);
		parser.check_sub_options("HKY", hky_sub_opts);

		parser.check_option_arg_range("hky_params", 0.00000001, 20.0);
		parser.check_option_arg_range("gtr_params", 0.0, 10.0);
		parser.check_option_arg_range("indel_params", 0.0, 0.99);
		parser.check_option_arg_range("initAlpha", 0.00000001, 100.0);

		if (parser.option("h"))
		{
			// display all the command line options

			cout << Definitions::notice;


			cout << "Usage: HMM --in input_file --(GTR|HKY|LG) [--gtr_params ... | --hky_params ... --indel_params ...]\n";
			parser.print_options();
			throw HmmException("Exiting...\n");
		}

		if (!parser.option("in"))
		{
			cout << "paHMM-Tree - distance-based statistical phylogenetic tree estimation version 0.1512 \n\n";
		    cout << "Usage: HMM --in input_file --(GTR|HKY|LG) [--gtr_params ... | --hky_params ... --indel_params ...]\n";
		    parser.print_options();
			throw HmmException("Please specify the input file\n");
		}

		if (!(parser.option("GTR") || parser.option("HKY") || parser.option("LG")))
		{
					cout << "paHMM-Tree - distance-based statistical phylogenetic tree estimation version 0.1512 \n\n";
				    cout << "Usage: HMM --in input_file --(GTR|HKY|LG) [--gtr_params ... | --hky_params ... --indel_params ...]\n";
				    parser.print_options();
					throw HmmException("Please specify a valid substitution model \n");
		}

		parser.check_option_arg_range("estimateAlpha", 0, 1);
		parser.check_option_arg_range("rateCat", 0, 1000);


	}
	catch (exception& e)
	{
	        throw HmmException(e.what());
	}
}

vector<double> CommandReader::getSubstParams()
{
	int i;
	vector<double> vec;
	if (parser.option("HKY"))
	{
		if(parser.option("hky_params"))
		{
			for (i=0; i< 1; i++)
			{
				DEBUG("HKY parameter " << i <<  ": " << parser.option("hky_params").argument(i));
				vec.push_back(atof(parser.option("hky_params").argument(i).c_str()));
			}
		}
		else throw HmmException("Model parameters not specified");
	}
	else if (parser.option("GTR"))
	{
		if (parser.option("gtr_params"))
		{
			for (i=0; i< 5; i++)
			{
				DEBUG("GTR parameter " << i <<  ": " << parser.option("gtr_params").argument(i));
				vec.push_back(atof(parser.option("gtr_params").argument(i).c_str()));
			}
		}
		else throw HmmException("Model parameters not specified");
	}
	else if (parser.option("LG")){}
	else throw HmmException("Model not specified");

	return vec;
}

vector<double> CommandReader::getIndelParams()
{
	int i;
	vector<double> vec;
	if (parser.option("indel_params"))
	{
		for (i=0; i< 2; i++)
		{
			DEBUG("indel parameters: " << i << ": " << parser.option("indel_params").argument(i));
			vec.push_back(atof(parser.option("i").argument(i).c_str()));
		}
	}
	else throw HmmException("Indel model not specified");
	return vec;
}

IParser* CommandReader::getParser() throw (HmmException&)
{
	if (parser.option("in"))
	{
		return new FileParser((string(parser.option("in").argument())).c_str());
	}
	else
		throw HmmException("input file not specified");
}

} /* namespace EBC */
