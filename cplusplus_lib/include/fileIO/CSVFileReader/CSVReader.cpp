#ifndef __CSVREADER_CPP_
#define __CSVREADER_CPP_
#include <iostream>
#include <sstream> //for stringstream
#include <iterator>
#include <boost/tokenizer.hpp>
#include <boost/smart_ptr.hpp>
using std::cout;
using std::endl;

#include "CSVReader.h"

std::vector<std::string> CSVRow::getNextLineAndSplitIntoTokens(std::istream& str)
{
	std::vector<std::string>   result;
	std::string                line;
	std::getline(str, line);

	//parse line with boost::tokenizer library
	boost::tokenizer<> tok(line);
	for (boost::tokenizer<>::iterator beg = tok.begin(); beg != tok.end(); ++beg) {
		result.push_back(*beg);
	}

	return result;
}

void CSVRow::readNextRow(std::istream& str)
{
	std::string         line;
	std::getline(str, line);

	m_data.clear();

	//parse line with boost::tokenizer library
	//boost::tokenizer<> tok(line);
	//for(boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end(); ++beg){
	//seperator	
	std::string sep1("");//dont let quoted arguments escape themselves
	//std::string sep2(" ");//split on spaces
	std::string sep2(",");//split on comma
	std::string sep3("\"\'");//let it have quoted arguments

	boost::escaped_list_separator<char> els(sep1, sep2, sep3);
	boost::tokenizer< boost::escaped_list_separator<char> > tok(line, els);

	for (boost::tokenizer< boost::escaped_list_separator<char> >::iterator beg = tok.begin(); beg != tok.end(); ++beg) {
		m_data.push_back(*beg);
	}
}


#endif
