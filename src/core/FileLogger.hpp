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

#ifndef FILE_LOGGER_HPP
#define FILE_LOGGER_HPP


#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace EBC
{


class FileLogger
{


public:
	enum logType { L_ERR, L_WARN, L_INF, L_DBG, L_DMP };

	static FileLogger& DumpLogger();
	static FileLogger& DebugLogger();
	static FileLogger& ErrorLogger();
	static FileLogger& InfoLogger();
	static FileLogger& WarningLogger();
	//Defualt logger without anny prefixes etc. Logs essential info like
	//estimated parameters and the tree
	static FileLogger& Logger();

	void activate()
	{
		this->active = true;
	}
	void setCerr()
	{
		this->stderrout = true;
	}

	static void start(logType priority, const std::string& fname)
	{
		switch(priority)
		{
		case L_DMP:
			dmpL.activate();
		case L_DBG:
			dbgL.activate();
		case L_INF:
			infL.activate();
		case L_WARN:
			wrnL.activate();
		default:
			errL.activate();
			clnL.activate();
			errL.setCerr();
		}
		logFile.open (fname.c_str());
	}
	static void stop()
	{
	    if (logFile.is_open())
	    {
	        logFile.close();
	    }
	}
/*
	void Debug(std::string& text)
	{
		if (this->priority > L_INF)
			this->logFile << "[DEBUG] " << text;
	}

	void Info(std::string& text)
	{
		if (this->priority > L_WARN)
				this->logFile << "[INFO] " << text;
	}

	void Warn(std::string& text)
	{
		if (this->priority > L_ERR)
			this->logFile << "[WARN] " << text;
	}

	void Error(std::string& text)
	{
		this->logFile << "[ERROR] " << text;
		std::cerr << "[ERROR] " << text;
	}

	template <typename T>
	void Debug (const std::vector<T>& v)
	{
		if (this->priority < L_DBG)
			return;
		this->logFile << "[DEBUG] ";
		for(unsigned int i = 0; i < v.size(); i++)
		{
			this->logFile << v[i] << "\t";
		}
		if (v.size() != 0)
			this->logFile << std::endl;
	}

	template <typename T>
	void Info (const std::vector<T>& v)
	{
		if (this->priority < L_INF)
			return;
		this->logFile << "[INFO] ";
		for(unsigned int i = 0; i < v.size(); i++)
		{
			this->logFile << v[i] << "\t";
		}
		if (v.size() != 0)
			this->logFile << std::endl;
	}
*/
	template <typename T>
	friend FileLogger &operator << (FileLogger &logger, const T& param) {

		if(logger.active)
		{
			logFile << param;
			if (logger.stderrout)
				std::cerr << param;
			logFile.flush();
		}
		return logger;
	}

	template <typename T>
	friend FileLogger &operator << (FileLogger &logger, const std::vector<T>& v)
	{
		if (v.size() != 0 && logger.active)
		{
			for(unsigned int i = 0; i < v.size(); i++)
			{
				logFile << v[i] << "\t\t";
				if (logger.stderrout)
				{
					std::cerr << v[i] << "\t\t";
				}

			}
			logFile.flush();
		}
		return logger;
	}
public:
	bool active;
	bool stderrout;

private:
    FileLogger() : active(false), stderrout(false) {}
    FileLogger(const FileLogger& logger) {}
    //FileLogger& operator = (const FileLogger& logger) {}
	static FileLogger errL;
	static FileLogger clnL;
	static FileLogger wrnL;
	static FileLogger dbgL;
	static FileLogger dmpL;
	static FileLogger infL;
	static std::ofstream logFile;
};

}


#endif // FILELOGGER_HPP
