#ifndef __STATISTICS_GET_CONCENTRATION_H_
#define __STATISTICS_GET_CONCENTRATION_H_
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility> //for std::pair
#include <iterator> //for std::vector<T>::iterator
#include <vector>
#include <algorithm> //for std:copy


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This class is written to deal with statistics_get_concentration stuff.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class statistics_get_concentration {
private:
	typedef std::map<std::string, int> stat_map;
	typedef std::pair<std::string, int> str_int_p;
private:
	struct Compare_pair {
		bool operator()(const str_int_p &p1, const str_int_p &p2) const{
			return p1.second > p2.second;
		}
	};
private:
	typedef std::vector<str_int_p> str_int_v;


public:
	stat_map pathway_unordered_stat;
	str_int_v pathway_ordered_stat;
	//std::ofstream out_file;
public:
	statistics_get_concentration(std::string str_in, std::string str_out, int path_len = 1);
	~statistics_get_concentration();
public:
	//Read in file
	void read_in_file(std::string str_in);
	//Insert pathway
	void insert_pathway_stat(std::string in_pathway);
	//Sort by counting number and print to file, print only if pathway count >= path_len
	void sort_print_to_file_stat(std::string str_out = "./output/pathway_stat.dat", int path_len = 1);
};



#endif
