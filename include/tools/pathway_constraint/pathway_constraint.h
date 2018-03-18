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
		bool not_allowed_out_species_constraint = false;

	public:
		// species sink reaction, have to be one reaction from this vector (from this list)
		std::map<rsp::index_int_t, std::unordered_set< rsp::index_int_t > > species_sink_reaction_set_map;
		// reaction out species, have to be one species from this vector (from this list)
		std::map<rsp::index_int_t, std::unordered_set<rsp::index_int_t> > reaction_out_species_set_map;
		// a more general out species constraint, species in the set can not be out species of any reactions
		std::unordered_set<rsp::index_int_t> not_allowed_out_species_set;

		// must reaction species set, set the survial probability to 0.0, react probability to be 1.0
		std::unordered_set<rsp::index_int_t> must_react_species_set;

	public:
		pathway_constraint();
		~pathway_constraint();

	};

}

#endif
