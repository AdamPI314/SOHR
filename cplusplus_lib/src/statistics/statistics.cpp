#ifndef __STATISTICS_CPP_
#define __STATISTICS_CPP_
#include <iostream>
using std::cout;
using std::endl;

#include "../../include/statistics/statistics.h"

statistics::statistics() {
	;
}//statistics

statistics::~statistics() {
	;
}//~statistics

//Insert pathway
void statistics::insert_pathway_stat(std::string in_pathway) {
	pathway_unordered_stat[in_pathway] += 1;
}

//Sort by counting number and print to file, only print the pathway with counter >= path_len
void statistics::sort_print_to_file_stat(std::string str_out, int path_len) {
	std::ofstream out_file;
	std::copy(pathway_unordered_stat.begin(), pathway_unordered_stat.end(), std::back_inserter(pathway_ordered_stat));
	std::sort(pathway_ordered_stat.begin(), pathway_ordered_stat.end(), Compare_pair());
	try
	{
		//open to use
		out_file.open(str_out.c_str());
		for (str_int_v::iterator iter = pathway_ordered_stat.begin(); iter != pathway_ordered_stat.end(); ++iter)
		{
			if ((*iter).second >= path_len) {
				out_file << (*iter).first << "," << (*iter).second << std::endl;
			}
		}
	}//try
	catch (std::ofstream::failure e) {
		std::cerr << "Exception opening/reading/closing file\n";
	}//catch

	out_file.close();
}

statistics::str_int_v statistics::return_sorted_path_count()
{
	std::copy(pathway_unordered_stat.begin(), pathway_unordered_stat.end(), std::back_inserter(pathway_ordered_stat));
	std::sort(pathway_ordered_stat.begin(), pathway_ordered_stat.end(), Compare_pair());
	return pathway_ordered_stat;
}

#endif
