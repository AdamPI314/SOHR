#ifndef __CLCF_PARSER_H_
#define __CLCF_PARSER_H_
//clcf represents command line and config file

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <fstream>
#include <iterator>

//A helper function to simplify the main part.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
	copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
	return os;
}

class clcf_parser
{
public:
	//define some groups
	clcf_parser();
	~clcf_parser();

private:
	po::variables_map vm;

private:
	//Declare a group of options that will be 
	//allowed only on command line
	po::options_description generic;

	//Declare a group of options that will be 
	//allowed both on command line and in
	//config file
	po::options_description config;

	//Hidden options, will be allowed both on command line and
	//in config file, but will not be shown to the user.
	po::options_description hidden;
public:
	//read command line
	int parse(int ac, char** av);
	//print out
	void print();
public:
	//copy vm to vm_out
	bool get_vm(po::variables_map &vm_out) {
		vm_out = this->vm;
		return true;
	}

};


class cf_parser
{
public:
	typedef double my_time_t;

public:
	//define some groups
	cf_parser();
	~cf_parser() { ; };

private:
	po::variables_map vm;

private:
	//Declare a group of options that will be in config file
	po::options_description config;

public:
	//read config file
	int parse(std::string filename);
	//print out
	void print();
public:
	//copy vm to vm_out
	bool get_vm(po::variables_map &vm_out) {
		vm_out = this->vm;
		return true;
	}

};

#endif
