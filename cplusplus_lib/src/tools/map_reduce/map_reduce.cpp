#ifndef __MAP_REDUCE_CPP_
#define __MAP_REDUCE_CPP_

#include "../../../include/tools/map_reduce/map_reduce.h"

std::map<std::string, int> merge_maps::operator() (std::map<std::string,int> &l, const std::map<std::string,int> &r) {
	//l.insert(r.begin(), r.end());
	for(std::map<std::string, int>::const_iterator itr=r.begin(); itr!=r.end(); ++itr)
	{
		l[itr->first]+=itr->second;
	}
	return l;
}

#endif
