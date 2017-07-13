#ifndef __STATISTICS_GET_CONC_CPP_
#define __STATISTICS_GET_CONC_CPP_
#include <iostream>
using std::cout;
using std::endl;

#include "statistics_get_conc.h"
statistics_get_conc::statistics_get_conc(std::string str_in, std::string str_out) {
	read_in_file(str_in);
	sort_print_to_file_stat(str_out);
}//statistics_get_conc

statistics_get_conc::~statistics_get_conc() {
	;
}//~statistics_get_conc

//Read in file
void statistics_get_conc::read_in_file(std::string str_in) {
	std::string str_tmp; int num_tmp;
	try {
		std::ifstream in_file(str_in.c_str());

		while (!in_file.eof()) {
			if ((in_file >> str_tmp) && (in_file >> num_tmp)) {
				int index = str_tmp.find_last_of('S');
				//cout<<index<<"\t";
				//cout<<str_tmp.substr(index)<<"\t"<<num_tmp<<endl;
				pathway_unordered_stat[str_tmp.substr(index)] += num_tmp;

			}
		}

		in_file.close();
	}
	catch (std::ifstream::failure e) {
		std::cerr << "Exception opening/reading/closing file\n";
	}
}

//Insert pathway
void statistics_get_conc::insert_pathway_stat(std::string in_pathway) {
	pathway_unordered_stat[in_pathway] += 1;
}

//Sort by counting number and print to file
void statistics_get_conc::sort_print_to_file_stat(std::string str_out) {
	std::ofstream out_file;
	std::copy(pathway_unordered_stat.begin(), pathway_unordered_stat.end(), std::back_inserter(pathway_ordered_stat));
	std::sort(pathway_ordered_stat.begin(), pathway_ordered_stat.end(), Compare_pair());
	try
	{
		//open to use
		out_file.open(str_out.c_str());
		for (str_int_v::iterator iter = pathway_ordered_stat.begin(); iter != pathway_ordered_stat.end(); ++iter)
		{
			out_file << (*iter).first << " " << (*iter).second << std::endl;
		}
	}//try
	catch (std::ofstream::failure e) {
		std::cerr << "Exception opening/reading/closing file\n";
	}//catch

	out_file.close();
}

#endif
