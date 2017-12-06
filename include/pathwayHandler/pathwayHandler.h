#ifndef __PATHWAYHANDLER_
#define __PATHWAYHANDLER_
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_set>
#include <utility> //for std::pair
#include <iterator> //for std::vector<T>::iterator
#include <vector>
#include <algorithm> //for std:copy
#include <sstream>


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This class is written to deal with pathway analysis stuff.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//deal with pathway stuff
class pathwayHandler {
private:
public:
	//parse file named "pathway_stat.dat" to get a vector of topN pathways
	static void get_pathway(std::string filename, std::vector<std::string> &pathway_vec, size_t topN);
	//parse the pathway ends with spe
	static void pathway_ends_with(std::string spe, const std::vector<std::string>& str_vec_in, std::vector<std::string>& str_vec_out, size_t topN = 10);
	static void pathway_starts_with(std::string spe, const std::vector<std::string>& str_vec_in, std::vector<std::string>& str_vec_out, size_t topN = 10);

	//parse the pathway ends with spe from a multimap
	static std::vector<std::string> pathway_ends_with(const std::multimap<double, std::string, std::greater<double> > &p_map, std::size_t n_spe, std::size_t topN);
	static std::vector<std::string> pathway_ends_with(const std::multimap<double, std::string> &p_map, std::size_t n_spe, std::size_t topN);

	static void merge_pathway(std::vector<std::string> &v1, std::vector<std::string> v2);
	

public:
	pathwayHandler();
	~pathwayHandler();
};





#endif
