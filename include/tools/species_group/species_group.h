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

	class spe_species_group {
	public:
		//this contains all the reaction between two species, for example s1, and s2, all reactions from s1 to s2
		//and reaction from s2 to s1 will be calculated
		std::vector< std::map<std::pair<std::size_t, std::size_t>, std::set<rxn_c1_c2> > > species_group_pairs_rxns;

	public:
		spe_species_group();
		~spe_species_group();
	};

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
		std::vector< std::map<std::pair<std::size_t, std::size_t>, std::set<rxn_c1_c2> > > species_chattering_group_pairs_rxns;

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
