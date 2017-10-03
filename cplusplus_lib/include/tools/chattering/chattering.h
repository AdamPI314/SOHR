#ifndef __CHATTERING_H_
#define __CHATTERING_H_
#include <vector>
#include <map>

namespace chattering_sr {
	class chattering {
	public:
		//chattrering includes chattering reaction and chattering species
		std::vector<std::size_t> chattering_reaction_index;
		//chattering species in pair, read direcltly from file
		std::vector<std::vector<std::size_t> > chattering_spe;

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
