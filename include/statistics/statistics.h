#ifndef __STATISTCIS_H_
#define __STATISTCIS_H_
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility> //for std::pair
#include <iterator> //for std::vector<T>::iterator
#include <vector>
#include <algorithm> //for std:copy


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This class is written to deal with statistics stuff.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class statistics {
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
public:
	const stat_map& get_pathway_unordered_map() const {
		return pathway_unordered_stat;
	}
public:
	statistics();
	~statistics();
public:
	//Insert pathway
	void insert_pathway_stat(std::string in_pathway);
	//Insert unordered map
	void insert_unordered_map(stat_map& pathway_unordered_stat_in) {
		pathway_unordered_stat.clear();
		pathway_unordered_stat = pathway_unordered_stat_in;
	}
	//Sort by counting number and print to file, only print the pathway with counter >= path_len
	void sort_print_to_file_stat(std::string str_out = "./output/pathway_stat.csv", int path_len = 1);
	// sort and return a vector of pair<str, int>
	str_int_v return_sorted_path_count();
};



#endif
