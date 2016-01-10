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


#include "core/HmmException.hpp"
#include "core/Definitions.hpp"

namespace EBC {

HmmException::HmmException()
{
	// General message
	msg = "Parse Exception";
}

HmmException::~HmmException() throw()
{
}

HmmException::HmmException(string message) : msg(message)
{
}

const char* HmmException::what() const throw ()
{
	return this->msg.c_str();
}

} /* namespace EBC */
