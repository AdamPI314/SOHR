#ifndef __PATHWAYHANDLER_CPP_
#define __PATHWAYHANDLER_CPP_
#include <iostream>
#include <sstream> //for stringstream
#include <iterator>
#include <boost/tokenizer.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/algorithm/string/predicate.hpp> //for boost::algorithm::ends_with
#include <boost/lexical_cast.hpp>
#include <numeric>
using std::cout;
using std::endl;

#include "pathwayHandler.h"
pathwayHandler::pathwayHandler() {
	;
}//pathwayHandler

pathwayHandler::~pathwayHandler() {
	;
}//~pathwayHandler


//parse file named "pathway_stat.dat" to get a vector of topN pathways
void pathwayHandler::get_pathway(std::string filename, std::vector<std::string> &pathway_vec, size_t topN) {
	std::ifstream fin(filename.c_str());
	pathway_vec.resize(0);

	std::string s; size_t icounter = 0;
	//while((icounter++<topN)&&(fin.good())){
	while ((icounter++ < topN) && (fin.good())) {
		try {
			std::getline(fin, s);
		}
		catch (std::ifstream::failure &e) {
			std::cerr << e.what() << "Exception opening/reading/closing file\n";
		}
		// there may be blank lines
		if (s.size() > 0) {
			boost::tokenizer<> tok(s);
			//just push back pathway name
			pathway_vec.push_back(std::string(*(tok.begin())));
			//for(boost::tokenizer<>::iterator beg=tok.begin(); beg!=tok.end(); ++beg){
				//cout << *beg << "\n";
			//}
		}
	}

	fin.clear(); fin.close();
}

void pathwayHandler::pathway_ends_with(std::string spe, const std::vector<std::string>& str_vec_in, std::vector<std::string>& str_vec_out, size_t topN) {
	size_t icount = 0;
	for (size_t i = 0; (icount < topN) && (i < str_vec_in.size()); ++i)
	{
		if (boost::algorithm::ends_with(str_vec_in[i], spe)) {
			++icount;
			str_vec_out.push_back(str_vec_in[i]);
		}
	}
}

void pathwayHandler::pathway_starts_with(std::string spe, const std::vector<std::string>& str_vec_in, std::vector<std::string>& str_vec_out, size_t topN)
{
	size_t icount = 0;
	for (size_t i = 0; (icount < topN) && (i < str_vec_in.size()); ++i)
	{
		if (boost::algorithm::starts_with(str_vec_in[i], spe)) {
			++icount;
			str_vec_out.push_back(str_vec_in[i]);
		}
	}
}

std::vector<std::string> pathwayHandler::pathway_ends_with(const std::multimap<double, std::string, std::greater<double>>& p_map, std::size_t n_spe, std::size_t topN)
{
	std::vector<std::string> v;
	std::vector<std::size_t> v_counter(n_spe, 0);
	for (auto x : p_map) {
		// if all species has topN path endwith itself, break
		if (std::accumulate(v_counter.begin(), v_counter.end(), 0) >= (int)(n_spe*topN))
			break;
		// iterate over species
		for (std::size_t i = 0; i < n_spe; ++i) {
			if (v_counter[i] >= topN)
				continue;
			auto spe_name = std::string("S") + boost::lexical_cast<std::string>(i);
			if (boost::algorithm::ends_with(x.second, spe_name)) {
				v.push_back(x.second);
				v_counter[i] += 1;
			}
		}
	}

	return v;
}

std::vector<std::string> pathwayHandler::pathway_ends_with(const std::multimap<double, std::string> &p_map, std::size_t n_spe, std::size_t topN) {
	std::vector<std::string> v;
	std::vector<std::size_t> v_counter(n_spe, 0);
	for (auto x : p_map) {
		// if all species has topN path endwith itself, break
		if (std::accumulate(v_counter.begin(), v_counter.end(), 0) >= (int)(n_spe*topN))
			break;
		// iterate over species
		for (std::size_t i = 0; i < n_spe; ++i) {
			if (v_counter[i] >= topN)
				continue;
			auto spe_name = std::string("S") + boost::lexical_cast<std::string>(i);
			if (boost::algorithm::ends_with(x.second, spe_name)) {
				v.push_back(x.second);
				v_counter[i] += 1;
			}
		}
	}

	return v;
}

void pathwayHandler::merge_pathway(std::vector<std::string>& v1, std::vector<std::string> v2)
{
	std::unordered_set<std::string> us;
	for (auto s1 : v1)
		us.insert(s1);
	for (auto s2 : v2)
		us.insert(s2);

	v1.clear();
	for (auto e : us)
		v1.push_back(e);

}


#endif
