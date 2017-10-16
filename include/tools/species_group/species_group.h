#ifndef __SPECIES_GROUP_H_
#define __SPECIES_GROUP_H_
#include <vector>
#include <set>
#include <map>
#include <tuple>

namespace species_group_sr {
	struct rxn_c1_c2 {
		//reaction index
		std::size_t r_idx;
		//coefficient of reactant
		double c1;
		//coeffcieant of product
		double c2;

		bool operator < (const rxn_c1_c2 &rhs) const {
			return std::tie(this->r_idx, this->c1, this->c2) < std::tie(rhs.r_idx, rhs.c1, rhs.c2);
		}
	};

	typedef std::map<std::pair<std::size_t, std::size_t>, std::vector<rxn_c1_c2> > species_group_pairs_rxns_t;
	//out species of a single species, suppose following a single atom
	typedef std::map<std::size_t, std::map<std::size_t, std::vector<rxn_c1_c2> > > out_species_rxns_t;

	//all species group, only one group, from a mass conservation point of view, from example to 3H2 => 2H3 and H2 => 2H reactions
	//there is no H2-->H3 transition probability thing, node to node transition, what we define here is presume we are following
	//a single atom, the branching ratio H2 becomes either H3 or H, suppose the reaction rates are the same, then p(H2, H3) vs. p(H2, H)
	//is defined as (3/4*6/6) vs. (1/4*2/2), in others words, it is reaction branching ratio * species branching ratio
	class species_group_base {
	public:
		//this contains all the reaction between two species, for example s1, and s2, all reactions from s1 to s2
		//and reaction from s2 to s1 will be calculated
		species_group_pairs_rxns_t species_group_pairs_rxns;
		out_species_rxns_t out_species_rxns;

	public:
		species_group_base();
		~species_group_base();
	};

	typedef std::map<std::pair<std::size_t, std::size_t>, std::set<rxn_c1_c2> > species_chattering_group_pairs_rxns_t;

	class chattering {
	public:
		//chattrering includes chattering reaction and chattering species from input file directly
		std::vector<std::vector<std::size_t> > chattering_rxn_idx_from_file;
		//chattering species in pair, read direcltly from file
		std::vector<std::vector<std::size_t> > chattering_spe_idx_from_file;

		//unique chattering species
		std::set<std::size_t> unique_chattering_species;
		//unique chattering reactions
		std::set<std::size_t> unique_chattering_reactions;

		//unique chattering reaction index calculated from chattering group
		//this contains all the reaction between two species, for example s1, and s2, all reactions from s1 to s2
		//and reaction from s2 to s1 will be calculated
		std::vector<species_chattering_group_pairs_rxns_t> species_chattering_group_pairs_rxns;

	public:
		//vector of species chattering groups
		std::vector<std::vector<std::size_t> > species_chattering_group_mat;
		//map species index -> its chattering group id, and index in that group 
		std::map<std::size_t, std::pair<std::size_t, std::size_t> > spe_idx_2_chattering_group_id_idx;
		//map species index -> a index in a flatten chattering super group
		std::map<std::size_t, std::size_t> spe_idx_2_super_group_idx;

	public:
		chattering();
		~chattering();

	};
}

#endif