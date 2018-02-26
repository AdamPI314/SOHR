#ifndef __SPECIES_GROUP_CPP_
#define __SPECIES_GROUP_CPP_

#include "../../../include/tools/species_group/species_group.h"

namespace species_group_sr {
	chattering::chattering()
	{
	}
	chattering::~chattering()
	{
	}
	bool chattering::is_in_same_chattering_group(rsp::index_int_t s1, rsp::index_int_t s2)
	{
		if (this->unique_chattering_species.find(s1) == this->unique_chattering_species.end() ||
			this->unique_chattering_species.find(s2) == this->unique_chattering_species.end())
			return false;

		return this->spe_idx_2_chattering_group_id_idx.at(s1) == this->spe_idx_2_chattering_group_id_idx.at(s2);
	}
	species_group_base::species_group_base()
	{
	}
	species_group_base::~species_group_base()
	{
	}
}

#endif
