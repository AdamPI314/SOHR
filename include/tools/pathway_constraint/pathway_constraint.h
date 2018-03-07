#ifndef __PATHWAY_CONSTRAINT_H_
#define __PATHWAY_CONSTRAINT_H_
#include <vector>
#include <set>
#include <map>
#include <tuple>
#include "../../relationshipParser/relationshipParser.h"

namespace pathway_constraint_sr {
	namespace rsp = relationshipParser_sr;

	class pathway_constraint {
	public:
		bool reaction_constraint = false;
		bool species_constraint = false;

	public:
		// the ith reaction, have to be one reaction from this vector (from this list)
		std::map<rsp::index_int_t, std::vector< rsp::index_int_t > > ith_reaction_vec;
		// the ith species, have to be one species from this vector (from this list)
		std::map<rsp::index_int_t, std::vector<rsp::index_int_t> > ith_species_vec;

	public:
		pathway_constraint();
		~pathway_constraint();

	};

}

#endif
