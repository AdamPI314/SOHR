#ifndef __PATHWAY_CONSTRAINT_H_
#define __PATHWAY_CONSTRAINT_H_
#include <vector>
#include <unordered_set>
#include <map>
#include <tuple>
#include "../../relationshipParser/relationshipParser.h"

namespace pathway_constraint_sr {
	namespace rsp = relationshipParser_sr;

	class pathway_constraint {
	public:
		bool species_sink_through_reaction_constraint = false;
		bool reaction_out_species_constraint = false;

	public:
		// species sink reaction, have to be one reaction from this vector (from this list)
		std::map<rsp::index_int_t, std::unordered_set< rsp::index_int_t > > species_sink_reaction_set_map;
		// reaction out species, have to be one species from this vector (from this list)
		std::map<rsp::index_int_t, std::unordered_set<rsp::index_int_t> > reaction_out_species_set_map;

	public:
		pathway_constraint();
		~pathway_constraint();

	};

}

#endif
