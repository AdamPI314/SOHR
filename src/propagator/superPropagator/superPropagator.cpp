#ifndef __SUPERPROPAGATOR_CPP_
#define __SUPERPROPAGATOR_CPP_

#define PRINT_PRECISION 15

#include <numeric>
#include <limits>
#include <unordered_set>
#include <set>
#include <unordered_map>

#include <boost/property_tree/json_parser.hpp> //for json_reader
#include <boost/optional/optional.hpp> //for optional

#include "../../../include/tools/misc/fortran_routine_block_alias.h"
#include "../../../include/tools/misc/global_extern_vars.h"
#include "../../../include/mechanism/mechanism.h"
#include "../../../include/tools/union_find/unionFind.h"

//put it in the end, because it has "cubicSpline/interp_1d.h", which has to be after all boost libraries
#include "../../../include/tools/matrix/matrix_sr.h"
#include "../../../include/propagator/superPropagator/superPropagator.h"


namespace propagator_sr {
	superPropagator::superPropagator(std::vector<double> uncertainties_in, std::string cwd_in)
	{
		//current working directory
		this->cwd_pgt = cwd_in;

		std::copy(uncertainties_in.begin(), uncertainties_in.end(), std::back_inserter(this->uncertainties));

		//defer this to actually use time
		initiate_propagator();

	}//superPropagator


	superPropagator::~superPropagator()
	{
		for (std::size_t i = 0; i < spe_drc_pgt.size(); ++i) {
			if (spe_drc_pgt[i] != 0)
				delete spe_drc_pgt[i];
		}
		for (std::size_t i = 0; i < spe_drc_int_pgt.size(); ++i) {
			if (spe_drc_int_pgt[i] != 0)
				delete spe_drc_int_pgt[i];
		}
		for (std::size_t i = 0; i < spe_drc_int_time_pgt.size(); ++i) {
			if (spe_drc_int_time_pgt[i] != 0)
				delete spe_drc_int_time_pgt[i];
		}
		for (std::size_t i = 0; i < reaction_rate_pgt.size(); ++i) {
			if (reaction_rate_pgt[i] != 0)
				delete reaction_rate_pgt[i];
		}
		for (std::size_t i = 0; i < concentration_pgt.size(); ++i) {
			if (concentration_pgt[i] != 0)
				delete concentration_pgt[i];
		}

		if (time_temperature_pgt != 0)
			delete time_temperature_pgt;
		if (time_pressure_pgt != 0)
			delete time_pressure_pgt;
		if (temperature_time_pgt != 0)
			delete temperature_time_pgt;

		//chattering group
		for (std::size_t i = 0; i < chattering_group_k_data_pgt.size(); ++i) {
			if (chattering_group_k_pgt.size() == 0)
				break;
			if (chattering_group_k_pgt[i] != 0)
				delete chattering_group_k_pgt[i];
		}

		for (std::size_t i = 0; i < chattering_group_k_int_data_pgt.size(); ++i) {
			if (chattering_group_k_int_pgt.size() == 0)
				break;
			if (chattering_group_k_int_pgt[i] != 0)
				delete chattering_group_k_int_pgt[i];
		}

		for (std::size_t i = 0; i < chattering_group_k_int_data_pgt.size(); ++i) {
			if (chattering_group_k_int_time_pgt.size() == 0)
				break;
			if (chattering_group_k_int_time_pgt[i] != 0)
				delete chattering_group_k_int_time_pgt[i];
		}

		for (std::size_t i = 0; i < chattering_group_ss_prob_data_pgt.size(); ++i) {
			if (chattering_group_ss_prob_pgt.size() == 0)
				break;
			if (chattering_group_ss_prob_pgt[i] != 0)
				delete chattering_group_ss_prob_pgt[i];
		}

	}

	void superPropagator::initialize_cubic_spline_pointer()
	{
		std::fill(spe_drc_pgt.begin(), spe_drc_pgt.end(), (Linear_interp*)(0));
		std::fill(spe_drc_int_pgt.begin(), spe_drc_int_pgt.end(), (Linear_interp*)(0));
		std::fill(spe_drc_int_time_pgt.begin(), spe_drc_int_time_pgt.end(), (Linear_interp*)(0));
		std::fill(reaction_rate_pgt.begin(), reaction_rate_pgt.end(), (Linear_interp*)(0));
		std::fill(concentration_pgt.begin(), concentration_pgt.end(), (Linear_interp*)(0));

		time_temperature_pgt = (Linear_interp*)0;
		time_pressure_pgt = (Linear_interp*)0;
		temperature_time_pgt = (Linear_interp*)0;

		//chattering group
		std::fill(chattering_group_k_pgt.begin(), chattering_group_k_pgt.end(), (Linear_interp*)(0));
	}

#if defined(__CHEMKIN_AVAILABLE_)

	void superPropagator::update_temporary_data_pgt(const int nkk, const int neq, const double ti, double * const c_t, const double * const CDOT_t, const double * const DDOT_t, const double * const FWDR_t, const double * const REVR_t, const double * const xgst)
	{
		double reaction_rate_tmp = 0.0;
		//use the default dt to print out stuff
		time_data_pgt.push_back(ti);

		//destruction rate const
		for (int i = 0; i < nkk; ++i)//[for
		{
			if (c_t[i] <= 0.0)
			{
				c_t[i] = 0.0;
				spe_drc_data_pgt[i].push_back(c_t[i]);
			}
			else
			{
				//just need the destruction rate const of species
				spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
			}
		}//for]

		 //Rates of each reactions
		for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
			reaction_rate_tmp = 0.0;
			//check reactionNetwork_chemkin_index_map
			for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
				if (reactionNetwork_chemkin_index_map[i][j] > 0)
					//Fortran style index to C/C++ style index
					reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1];
				else if (reactionNetwork_chemkin_index_map[i][j] < 0)
					//Fortran style index to C/C++ style index
					reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1];
			}//for2
			reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

		}//for1

		 //print concentration and temperature
		for (int i = 0; i < nkk; ++i) {
			if (c_t[i] < 0)
				c_t[i] = 0.0;
			concentration_data_pgt[i].push_back(c_t[i]);
		}

		temperature_data_pgt.push_back(xgst[neq - 1]);
		pressure_data_pgt.push_back(ckstore.pressure);
	}

#endif // __CHEMKIN_AVAILABLE_

	void superPropagator::update_all_species_group_pairs_reactions(const std::vector<rsp::spe_info_base>& species_network_v, const std::vector<rsp::reaction_info_base>& reaction_network_v, std::string atom_followed)
	{
		species_group_sr::species_group_pairs_rxns_t pairs_rxns_map_tmp;
		//species vector
		for (std::size_t s_idx1 = 0; s_idx1 < species_network_v.size(); ++s_idx1) {//s_idx1
			//don't contain atom_followed
			if (species_network_v[s_idx1].spe_component.at(atom_followed) < 1)
				continue;

			for (auto rxn_coef : species_network_v[s_idx1].reaction_k_index_s_coef_v) {
				auto r_idx = rxn_coef.first;
				auto c1 = rxn_coef.second;

				for (auto s_idx_coef : reaction_network_v[r_idx].out_spe_index_weight_v_map.at(atom_followed)) {
					auto s_idx2 = s_idx_coef.first;
					auto c2 = s_idx_coef.second;

					//don't contain atom_followed
					if (species_network_v[s_idx2].spe_component.at(atom_followed) < 1 || c2 <= 0)
						continue;

					std::pair<std::size_t, std::size_t> s1_s2_p(s_idx1, s_idx2);

					species_group_sr::rxn_c1_c2 r_c1_c2;
					r_c1_c2.r_idx = r_idx;
					r_c1_c2.c1 = c1;
					r_c1_c2.c2 = c2 / species_network_v[s_idx2].spe_component.at(atom_followed);

					pairs_rxns_map_tmp[s1_s2_p].push_back(r_c1_c2);


				}

			}

		}//s_idx1
		this->sp_all_species_group_pgt->species_group_pairs_rxns = pairs_rxns_map_tmp;

		//write to file
		std::ofstream fout1((this->cwd_pgt + std::string("/output/species_pairs_reactions_coefs.csv")).c_str());
		for (auto x : this->sp_all_species_group_pgt->species_group_pairs_rxns)
		{
			auto s1_s2 = x.first;
			fout1 << "(" << s1_s2.first << "," << s1_s2.second << ")" << "-->"
				<< "(" << species_network_v[s1_s2.first].spe_name << "," << species_network_v[s1_s2.second].spe_name << ")" << std::endl;

			auto vec = x.second;
			for (auto rxn_c1_c2 : vec)
			{
				fout1 << "\t" << "(" << rxn_c1_c2.r_idx << "," << rxn_c1_c2.c1 << "," << rxn_c1_c2.c2 << ")" << std::endl;
				fout1 << "\t" << reaction_network_v[rxn_c1_c2.r_idx].reaction_name << std::endl;
			}

		}

		fout1.clear(); fout1.close();

		//write to json file
		boost::property_tree::ptree pt_root1;

		int counter = 0;
		for (auto x : this->sp_all_species_group_pgt->species_group_pairs_rxns) {

			auto s1_s2 = x.first;

			auto vec = x.second;
			for (auto rxn_c1_c2 : vec)
			{
				boost::property_tree::ptree pt_child1;

				pt_child1.put("from", s1_s2.first);
				pt_child1.put("to", s1_s2.second);

				pt_child1.put("r_idx", rxn_c1_c2.r_idx);
				pt_child1.put("c1", rxn_c1_c2.c1);
				pt_child1.put("c2", rxn_c1_c2.c2);
				pt_child1.put("r_name", reaction_network_v[rxn_c1_c2.r_idx].reaction_name);

				pt_root1.put_child(boost::lexical_cast<std::string>(counter + 1), pt_child1);
				++counter;
			}
		}

		std::ofstream fout_json;
		try
		{
			//open to use
			fout_json.open((this->cwd_pgt + std::string("/input/species_pairs_reactions_coefs.json")).c_str());
			boost::property_tree::write_json(fout_json, pt_root1);
		}//try
		catch (std::ofstream::failure e) {
			std::cerr << "Exception opening/reading/closing file\n";
		}//catch

		fout_json.close();

	}

	void superPropagator::update_all_species_out_species_reactions()
	{
		species_group_sr::out_species_rxns_t out_species_rxns_map_tmp;
		for (auto x : this->sp_all_species_group_pgt->species_group_pairs_rxns) {
			auto s1_idx = x.first.first;
			auto s2_idx = x.first.second;
			auto vec = x.second;

			out_species_rxns_map_tmp[s1_idx][s2_idx] = vec;
		}

		this->sp_all_species_group_pgt->out_species_rxns = out_species_rxns_map_tmp;
	}


	std::shared_ptr<species_group_sr::species_group_base> superPropagator::get_sp_of_all_species_group()
	{
		return this->sp_all_species_group_pgt;
	}

	void superPropagator::find_chattering_group_using_union_find(const std::vector<rsp::spe_info_base>& species_network_v)
	{
		//since there might be groups of chattering species, such as A=B, B=C, C=D, A,B,C,D belongs to the same fast transition group
		//use union find here
		for (auto x : this->sp_chattering_pgt->chattering_spe_pair_idx_from_file) {
			for (auto y : x) {
				this->sp_chattering_pgt->unique_chattering_species.insert(y);
			}
		}

		std::unordered_map<int, int> label_2_idx;
		std::unordered_map<int, int> idx_2_label;
		int counter = 0;
		for (auto x : this->sp_chattering_pgt->unique_chattering_species) {
			label_2_idx.emplace(counter, x);
			idx_2_label.emplace(x, counter++);
		}

		UnionFind uf((int)label_2_idx.size());
		for (std::size_t i = 0; i < this->sp_chattering_pgt->chattering_spe_pair_idx_from_file[0].size(); ++i) {
			uf.unite(idx_2_label[this->sp_chattering_pgt->chattering_spe_pair_idx_from_file[0][i]],
				idx_2_label[this->sp_chattering_pgt->chattering_spe_pair_idx_from_file[1][i]]);
		}

		//find fast transition groups
		std::unordered_set<int> unique_root_species;
		for (auto x : this->sp_chattering_pgt->chattering_spe_pair_idx_from_file) {
			for (auto y : x) {
				//find root
				auto s = label_2_idx[uf.root(idx_2_label[y])];
				unique_root_species.insert(s);
			}
		}
		std::map<int, std::set<int> > fast_transition_group_spe;
		for (auto x : unique_root_species) {
			fast_transition_group_spe.emplace(x, std::set<int>({ x }));
		}
		for (auto x : this->sp_chattering_pgt->chattering_spe_pair_idx_from_file) {
			for (auto y : x) {
				//find root
				auto s = label_2_idx[uf.root(idx_2_label[y])];
				if (fast_transition_group_spe[s].find(y) == fast_transition_group_spe[s].end()) {
					fast_transition_group_spe[s].insert(y);
				}
			}
		}

		//group_counter, spe_counter
		int group_counter = 0;
		for (auto s_c_g_pair : fast_transition_group_spe) {
			std::vector<rsp::index_int_t> s_c_g;
			int spe_counter = 0;
			for (auto s_idx : s_c_g_pair.second) {
				s_c_g.push_back(s_idx);
				this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx[s_idx] = std::make_pair(group_counter, spe_counter);
				++spe_counter;
			}
			this->sp_chattering_pgt->species_chattering_group_mat.push_back(s_c_g);
			++group_counter;
		}

		//super_group_counter
		int super_group_counter = 0;
		for (auto g1 : this->sp_chattering_pgt->species_chattering_group_mat) {
			for (auto s_idx : g1) {
				this->sp_chattering_pgt->spe_idx_2_super_group_idx[s_idx] = super_group_counter;
				++super_group_counter;
			}
		}

		//write to json file
		boost::property_tree::ptree pt_root_tmp;
		for (std::size_t i = 0; i < this->sp_chattering_pgt->species_chattering_group_mat.size(); ++i) {
			boost::property_tree::ptree pt_child_tmp;
			for (std::size_t j = 0; j < this->sp_chattering_pgt->species_chattering_group_mat[i].size(); ++j) {
				pt_child_tmp.put(boost::lexical_cast<std::string>(j), this->sp_chattering_pgt->species_chattering_group_mat[i][j]);
			}
			pt_root_tmp.put_child(boost::lexical_cast<std::string>(i + rsp::INDICATOR), pt_child_tmp);
		}

		std::ofstream out_stream;
		try
		{
			//open to use
			//this->cwd + std::string("/input/setting.json"),  std::locale()
			out_stream.open(this->cwd_pgt + std::string("/output/chattering_group_info.json"));

			boost::property_tree::write_json(out_stream, pt_root_tmp);
		}//try
		catch (std::ofstream::failure e) {
			std::cerr << "Exception opening/reading/closing file\n";
		}//catch

		out_stream.close();

	}

	void superPropagator::update_chattering_group_pairs_reactions(const std::vector<rsp::spe_info_base>& species_network_v, const std::vector<rsp::reaction_info_base>& reaction_network_v, std::string atom_followed)
	{
		this->sp_chattering_pgt->species_chattering_group_pairs_rxns.resize(0);
		this->sp_chattering_pgt->species_chattering_group_pairs_rxns.clear();

		this->sp_chattering_pgt->species_chattering_group_pairs_rxns.resize(this->sp_chattering_pgt->species_chattering_group_mat.size());
		for (std::size_t group_i = 0; group_i < this->sp_chattering_pgt->species_chattering_group_mat.size(); ++group_i) {
			species_group_sr::species_chattering_group_pairs_rxns_t pairs_rxns_map_tmp;

			//species vector
			auto s_vec_tmp = this->sp_chattering_pgt->species_chattering_group_mat[group_i];
			for (auto s_idx1 : s_vec_tmp) {//s_idx1
				//don't contain atom_followed
				if (species_network_v[s_idx1].spe_component.at(atom_followed) < 1)
					continue;

				for (auto rxn_coef : species_network_v[s_idx1].reaction_k_index_s_coef_v) {
					auto r_idx = rxn_coef.first;
					auto c1 = rxn_coef.second;
					for (auto s_idx_coef : reaction_network_v[r_idx].out_spe_index_weight_v_map.at(atom_followed)) {
						auto s_idx2 = s_idx_coef.first;
						auto c2 = s_idx_coef.second;

						//don't contain atom_followed
						if (species_network_v[s_idx2].spe_component.at(atom_followed) < 1 || c2 <= 0)
							continue;

						//if these two species are both in current chattering group
						if (std::find(s_vec_tmp.begin(), s_vec_tmp.end(), s_idx2) != s_vec_tmp.end()) {
							std::pair<std::size_t, std::size_t> s1_s2_p(s_idx1, s_idx2);

							species_group_sr::rxn_c1_c2 r_c1_c2;
							r_c1_c2.r_idx = r_idx;
							r_c1_c2.c1 = c1;
							r_c1_c2.c2 = c2 / species_network_v[s_idx2].spe_component.at(atom_followed);

							pairs_rxns_map_tmp[s1_s2_p].insert(r_c1_c2);
						}//if these two species are both in current chattering group

					}

				}

			}//s_idx1

			this->sp_chattering_pgt->species_chattering_group_pairs_rxns[group_i] = pairs_rxns_map_tmp;
		}
	}

	std::shared_ptr<species_group_sr::chattering> superPropagator::get_sp_of_chattering()
	{
		return this->sp_chattering_pgt;
	}

	void superPropagator::update_chattering_spe_pair_from_file_pgt()
	{
		std::vector<std::vector<rsp::index_int_t> > Matrix(2, std::vector<rsp::index_int_t>());

		int counter = 0;

		for (auto key1 : this->pgt_pt.get_child("pathway.chattering_species")) {
			for (auto key2 : key1.second) {
				rsp::index_int_t x = key2.second.get_value<rsp::index_int_t>();
				//std::cout << x;
				if (counter % 2 == 0)
					Matrix[0].push_back(x);
				else
					Matrix[1].push_back(x);

				++counter;
			}
		}

		this->sp_chattering_pgt->chattering_spe_pair_idx_from_file = Matrix;
	}


	void superPropagator::subtract_chattering_reaction_contribution_from_species_drc_pgt()
	{
		//for all species, either cancel fast transitions in the same chattering group
		//or cancel specified reactions
		for (std::size_t group_i = 0; group_i < this->sp_chattering_pgt->species_chattering_group_pairs_rxns.size(); ++group_i) {
			for (auto x : this->sp_chattering_pgt->species_chattering_group_pairs_rxns[group_i]) {
				auto s1_s2_p = x.first;
				auto rxn_c1_c2_set = x.second;

				for (auto rxn_c1_c2 : rxn_c1_c2_set) {
					//see if current reaction in the "expected list"
					//if (this->sp_chattering_pgt->unique_chattering_reactions.count(rxn_c1_c2.r_idx) < 1)
					//	continue;
					for (std::size_t time_j = 0; time_j < this->time_data_pgt.size(); ++time_j) {
						if (this->concentration_data_pgt[s1_s2_p.first][time_j] > 0) {
							auto drc_tmp = rxn_c1_c2.c1 *this->reaction_rate_data_pgt[rxn_c1_c2.r_idx][time_j] / this->concentration_data_pgt[s1_s2_p.first][time_j];
							if (this->spe_drc_data_pgt[s1_s2_p.first][time_j] > drc_tmp)
								this->spe_drc_data_pgt[s1_s2_p.first][time_j] -= drc_tmp;
							else
								this->spe_drc_data_pgt[s1_s2_p.first][time_j] = 0.0;
						}
						else
							this->spe_drc_data_pgt[s1_s2_p.first][time_j] = 0.0;
					}//time
				}//s1 s2 reactions
			}//s1_s2_pairs
		}//chattering group
	}

	void superPropagator::update_chattering_group_K_prob_of_groupID_at_timeIDX_using_reference_X(std::size_t group_i, std::size_t time_j)
	{
		//sum over all chattering species
		this->chattering_group_k_data_pgt[group_i][time_j] = 0.0;

		//calculate concentration at this time
		double sum_conc = 0.0;
		for (auto spe_idx : this->sp_chattering_pgt->species_chattering_group_mat[group_i]) {
			auto s_g_idx = this->sp_chattering_pgt->spe_idx_2_super_group_idx[spe_idx];
			//fast equilibrium concentration, concentration at this time
			this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j] =
				this->evaluate_concentration_at_time(this->time_data_pgt[time_j], spe_idx);
			sum_conc += this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j];
		}

		//normalize within each chattering group
		for (auto spe_idx : this->sp_chattering_pgt->species_chattering_group_mat[group_i]) {
			auto s_g_idx = this->sp_chattering_pgt->spe_idx_2_super_group_idx[spe_idx];
			//check zero
			if (sum_conc > 0) {
				this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j] /= sum_conc;

				//as soon as fast equilibrium concentration is normalized
				//this->spe_drc_data_pgt[spe_idx][time_j] *= this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j];
				this->chattering_group_k_data_pgt[group_i][time_j] +=
					this->spe_drc_data_pgt[spe_idx][time_j] * this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j];
			}
			else {
				this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j] = 1.0 / this->sp_chattering_pgt->species_chattering_group_mat[group_i].size();

				//as soon as fast equilibrium concentration is normalized
				//this->spe_drc_data_pgt[spe_idx][time_j] *= this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j];
				this->chattering_group_k_data_pgt[group_i][time_j] +=
					this->spe_drc_data_pgt[spe_idx][time_j] * this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j];
			}
		}
	}

	void superPropagator::update_chattering_group_K_prob_of_groupID_at_timeIDX_using_SSA(std::size_t group_i, std::size_t time_j)
	{

		auto spe_vec = this->sp_chattering_pgt->species_chattering_group_mat[group_i];
		// transition matrix, to each group of species, calculate a transition matrix
		std::vector<std::vector<double> > transition_mat(spe_vec.size(), std::vector<double>(spe_vec.size(), 0.0));
		// calculate equilibrium concentration based on transition matrix
		// initially represent b vector in equation K X = b
		std::vector<double> b_vector(transition_mat.size(), 0);

		// inside chattering group
		for (auto c_g : this->sp_chattering_pgt->species_chattering_group_pairs_rxns) {
			for (auto p_r_m : c_g) {
				auto s1_s2_p = p_r_m.first;

				auto s_idx_1 = s1_s2_p.first;
				auto s_idx_2 = s1_s2_p.second;

				// both in group_i
				if (this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_1).first == (rsp::index_int_t)group_i &&
					this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_2).first == (rsp::index_int_t)group_i) {

					auto rxn_c1_c2_set = p_r_m.second;
					for (auto rxn_c1_c2 : rxn_c1_c2_set) {
						auto rxn_idx = rxn_c1_c2.r_idx;
						auto s_coef_1 = rxn_c1_c2.c1;
						auto s_coef_2 = rxn_c1_c2.c2;

						if (this->concentration_data_pgt[s_idx_1][time_j] != 0) {
							// auto drc_tmp = s_coef_1 * this->reaction_rate_data_pgt[rxn_idx][time_j] / this->concentration_data_pgt[s_idx_1][time_j];
							auto drc_tmp = this->reaction_rate_data_pgt[rxn_idx][time_j] / this->concentration_data_pgt[s_idx_1][time_j];
							// assume first order transition, if it is 5A==10B, too hard to deal with
							// speaking of transition matrix, to a A==B system, 
							// k_{AB} should be the matrix element on the left bottom corner
							transition_mat[this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_2).second]
								[this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_1).second] = drc_tmp * s_coef_2;

							// A sink term for itself, only for linear algebra
							transition_mat[this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_1).second]
								[this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_1).second] -= drc_tmp * s_coef_1;

						}//if

					}//rxn_c1_c2_set
				}
			}//s1_s2_pair
		}//inside chattrering_group

		//// outside chattering group
		//for (auto s_p_rxn_c1_c2 : this->sp_all_species_group_pgt->species_group_pairs_rxns) {
		//	auto s1_s2_p = s_p_rxn_c1_c2.first;

		//	auto s_idx_1 = s1_s2_p.first;
		//	auto s_idx_2 = s1_s2_p.second;

		//	// sink terms, s1 in group, s2 not in group
		//	auto s1_c_g_id = this->sp_chattering_pgt->get_chattering_group_id(s_idx_1);
		//	auto s2_c_g_id = this->sp_chattering_pgt->get_chattering_group_id(s_idx_2);
		//	if (s1_c_g_id == (rsp::index_int_t)group_i && s2_c_g_id == -1) {
		//		// means this is a sink term for current species
		//		auto rxn_c1_c2_vec = s_p_rxn_c1_c2.second;
		//		for (auto rxn_c1_c2 : rxn_c1_c2_vec) {
		//			auto rxn_idx = rxn_c1_c2.r_idx;
		//			auto s_coef_1 = rxn_c1_c2.c1;
		//			//auto s_coef_2 = rxn_c1_c2.c2;
		//			if (this->concentration_data_pgt[s_idx_1][time_j] != 0) {
		//				// auto drc_tmp = s_coef_1 * this->reaction_rate_data_pgt[rxn_idx][time_j] / this->concentration_data_pgt[s_idx_1][time_j];
		//				auto drc_tmp = this->reaction_rate_data_pgt[rxn_idx][time_j] / this->concentration_data_pgt[s_idx_1][time_j];
		//				// A sink term for itself, only for linear algebra
		//				transition_mat[this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_1).second]
		//					[this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_1).second] -= drc_tmp * s_coef_1;
		//			}//if
		//		}//rxn_c1_c2_vector
		//	}

		//	// source term, s1 not in group, s2 in group, calculate b vector
		//	if (s1_c_g_id == -1 && s2_c_g_id == (rsp::index_int_t)group_i) {
		//		// means this is a sink term for current species
		//		auto rxn_c1_c2_vec = s_p_rxn_c1_c2.second;
		//		for (auto rxn_c1_c2 : rxn_c1_c2_vec) {
		//			auto rxn_idx = rxn_c1_c2.r_idx;
		//			//auto s_coef_1 = rxn_c1_c2.c1;
		//			auto s_coef_2 = rxn_c1_c2.c2;
		//			if (this->concentration_data_pgt[s_idx_1][time_j] != 0) {
		//				// auto drc_tmp = s_coef_2 * this->reaction_rate_data_pgt[rxn_idx][time_j];
		//				// should be negative
		//				auto drc_tmp = this->reaction_rate_data_pgt[rxn_idx][time_j];
		//				b_vector[this->sp_chattering_pgt->spe_idx_2_chattering_group_id_idx.at(s_idx_2).second] -= drc_tmp * s_coef_2;
		//			}//if
		//		}//rxn_c1_c2_vector
		//	}

		//}

		//// using transition matrix without sink terms
		//double first_real_positive_eigenvalue;
		//auto ok = matrix_sr::cal_equilibrium_ratio_from_transition_matrix(transition_mat, first_real_positive_eigenvalue, b_vector);

		// solve linear equation, with sink and source terms
		if (time_j == 10) {
			for (std::size_t i = 0; i < transition_mat.size(); ++i) {
				for (std::size_t j = 0; j < transition_mat[0].size(); ++j) {
					std::cout << transition_mat[i][j] << ",";
				}
				std::cout << std::endl;
			}
			for (std::size_t i = 0; i < b_vector.size(); ++i)
				std::cout << b_vector[i] << ",";
			std::cout << "\n";
		}
		auto ok = matrix_sr::gaussian_jordan(transition_mat, b_vector);

		if (ok == false)
			return;

		for (size_t i = 0; i < this->sp_chattering_pgt->species_chattering_group_mat[group_i].size(); ++i) {
			auto spe_idx = this->sp_chattering_pgt->species_chattering_group_mat[group_i][i];
			// super group index
			auto s_g_idx = this->sp_chattering_pgt->spe_idx_2_super_group_idx[spe_idx];
			chattering_group_ss_prob_data_pgt[s_g_idx][time_j] = b_vector[i];

			// as soon as fast equilibrium concentration is normalized
			this->chattering_group_k_data_pgt[group_i][time_j] +=
				this->spe_drc_data_pgt[spe_idx][time_j] * this->chattering_group_ss_prob_data_pgt[s_g_idx][time_j];
		}

	}

	void superPropagator::update_drc_and_equilibrium_probability_of_chattering_group()
	{
		this->chattering_group_k_data_pgt.clear();
		this->chattering_group_k_data_pgt.resize(this->sp_chattering_pgt->species_chattering_group_mat.size());
		for (auto &x : this->chattering_group_k_data_pgt)
			x.assign(this->time_data_pgt.size(), 0.0);

		this->chattering_group_ss_prob_data_pgt.clear();
		this->chattering_group_ss_prob_data_pgt.resize(this->sp_chattering_pgt->spe_idx_2_super_group_idx.size());
		for (auto &x : this->chattering_group_ss_prob_data_pgt)
			x.assign(this->time_data_pgt.size(), 0.0);


		//first order approximation, instead of zero-th order, for example, A(k1)=B(k2)->(k3)C, 
		//if there is a fast transition between A and B, there we modify the k3, k3=k2/(k1+k2) * k3
		for (std::size_t group_i = 0; group_i < this->sp_chattering_pgt->species_chattering_group_mat.size(); ++group_i) {
			for (std::size_t time_j = 0; time_j < this->time_data_pgt.size(); ++time_j) {

				//update_chattering_group_K_prob_of_groupID_at_timeIDX_using_reference_X(group_i, time_j);
				update_chattering_group_K_prob_of_groupID_at_timeIDX_using_SSA(group_i, time_j);

			}//time j
		}//chattering group

		//chattering group destruction rate constant,pseudo first order destructive rate constant
		std::string tag = "";
		if (this->pgt_pt.get<std::string>("propagator.type") == std::string("dlsode")) {
			if (this->pgt_pt.get<std::string>("propagator.convert_molar_concentration_to_mole_fraction") == std::string("no"))
				tag += "dlsode_M";
			else if (this->pgt_pt.get<std::string>("propagator.convert_molar_concentration_to_mole_fraction") == std::string("yes"))
				tag += "dlsode_fraction";
		}
		else if (this->pgt_pt.get<std::string>("propagator.type") == std::string("ssa")) {
			tag += "ssa_number";
		}

		std::ofstream fout((this->cwd_pgt + std::string("/output/chattering_group_drc_") + tag + std::string(".csv")).c_str());
		for (size_t time_i = 0; time_i < chattering_group_k_data_pgt[0].size(); ++time_i) {
			for (size_t j = 0; j < chattering_group_k_data_pgt.size(); ++j) {
				fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << chattering_group_k_data_pgt[j][time_i];
				if (j < chattering_group_k_data_pgt.size() - 1)
					fout << ",";
			}
			fout << endl;
		}
		fout.clear(); fout.close();

		//print out a specific group ss prob
		std::size_t group_i_tmp = 3;
		fout.open((this->cwd_pgt + std::string("/output/chattering_group_ss_prob_") + tag + std::string(".csv")).c_str());
		for (size_t time_i = 0; time_i < chattering_group_ss_prob_data_pgt[0].size(); ++time_i) {
			for (std::size_t j = 0; j < this->sp_chattering_pgt->species_chattering_group_mat[group_i_tmp].size(); ++j) {
				auto x = this->sp_chattering_pgt->species_chattering_group_mat[group_i_tmp][j];
				auto s_g_id = this->sp_chattering_pgt->spe_idx_2_super_group_idx[x];
				fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << chattering_group_ss_prob_data_pgt[s_g_id][time_i];
				if (j < this->sp_chattering_pgt->species_chattering_group_mat[group_i_tmp].size() - 1)
					fout << ",";
			}
			fout << endl;
		}
		fout.clear(); fout.close();

	}

	void superPropagator::set_chattering_reaction_rates_to_zero_pgt()
	{
		for (auto c_g : this->sp_chattering_pgt->species_chattering_group_pairs_rxns) {
			for (auto p_r_m : c_g) {
				//auto s1_s2_p = p_r_m.first;
				auto rxn_c1_c2_set = p_r_m.second;
				for (auto rxn_c1_c2 : rxn_c1_c2_set) {
					auto rxn_idx = rxn_c1_c2.r_idx;
					std::fill(reaction_rate_data_pgt[rxn_idx].begin(),
						reaction_rate_data_pgt[rxn_idx].end(), 0.0);
				}//rxn_c1_c2_set
			}//s1_s2_pair
		}//chattrering_group

	}


	rsp::temperature_t superPropagator::return_target_temperature() const
	{
		return this->pgt_pt.get<rsp::temperature_t>("T.target_temperature");
	}

	void superPropagator::print_temperature_targeted_concentration_pgt(std::string outfile) const
	{
		std::ofstream out_file((this->cwd_pgt + outfile).c_str());
		double time = evaluate_time_at_temperature(return_target_temperature());
		//double time = time_data_pgt.back();
		for (size_t i = 0; i < concentration_pgt.size(); ++i) {
			out_file << evaluate_concentration_at_time(time, i);
			if (i != concentration_pgt.size() - 1)
				out_file << " ";
		}

		out_file.close();
	}


#if defined(__CHEMKIN_AVAILABLE_) && defined(__USE_LOSDE_)

	void superPropagator::print_pgt()
	{
		cout << "Print:" << endl;

		std::vector<int> duplicate_reaction_ind;
		rsp::relationshipParser::read_chem_out_duplicated_reaction(duplicate_reaction_ind, "./input/chem.out");
		cout << "duplicated reaction:(C/C++ style index) " << endl;
		std::copy(duplicate_reaction_ind.begin(), duplicate_reaction_ind.end(), std::ostream_iterator<int>(cout, " "));
		cout << endl;

		//print network reaction index and ChemKin style reaction index
		for (rsp::reactionNetwork_chemkin_index_map_t::iterator itr = reactionNetwork_chemkin_index_map.begin(); itr != reactionNetwork_chemkin_index_map.end(); ++itr) {//for
			std::cout << itr->first << "\t:";
			std::copy((*itr).second.begin(), (*itr).second.end(), std::ostream_iterator<int>(std::cout, "\t"));
			std::cout << "\n";
		}//for

		std::cout << "number of reactions:\t" << this->reaction_rate_data_pgt.size() << std::endl;


		cout << "Final temperature:\t" << temperature_data_pgt.back() << endl;

		int MM_t, KK_t, II_t, NFIT_t;
		mechanism::kinetics::indx(&MM_t, &KK_t, &II_t, &NFIT_t);
		cout << "mechanism total element count: " << MM_t << endl;
		cout << "mechanism total species count: " << KK_t << endl;
		cout << "mechanism total reaction count: " << II_t << endl;
		cout << "number of coefficients in fits to thermodynamic data for a temperature range;\nNFIT=number of coefficients in polynomial fits to CP/R plus 2: " << NFIT_t << endl << endl;

		const int neq = mechanism::kinetics::nkk() + 1, nkk = mechanism::kinetics::nkk(), nii = mechanism::kinetics::nii();

		//molar fraction
		double* x_t = new double[neq - 1]; double* y_t = new double[neq - 1];
		for (int i = 0; i < neq - 1; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }
		//Read 'Temperature(K)' and 'Pressure(atm)'.
		//Read	MOLAR FRACTION for species.
		double Temp, Pressure;
		read_configuration(Temp, Pressure, neq - 1, x_t);
		//	chem_init(Temp, Pressure, neq-1, x_t, this->cwd_pgt+std::string("/input/setting.cfg"));


		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore_.rhomass' is constant in the simulation.
		mechanism::kinetics::rhoy(&Pressure, &Temp, y_t, &ckstore.rhomass);

		//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
		mechanism::kinetics::py(&ckstore.rhomass, &Temp, y_t, &Pressure);

		//molar concentration.
		double* c_t = new double[nkk];
		mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);
		std::cout << "molar concentration:" << std::endl;
		std::copy(c_t, c_t + nkk, std::ostream_iterator<double>(cout, " ")); cout << endl;

		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//Returns the molar creation and destruction rates of the species given mass density, temperature(s) and mass fractions
		mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);
		std::cout << "molar creation and destruction rates, subroutine ckcdyr:" << std::endl;
		std::copy(CDOT_t, CDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(DDOT_t, DDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;

		std::fill(CDOT_t, CDOT_t + nkk, 0.0);
		std::fill(DDOT_t, DDOT_t + nkk, 0.0);
		std::copy(CDOT_t, CDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(DDOT_t, DDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;
		//The difference is how to treat the auto-catylytic reactions
		mechanism::kinetics::cdc(&Temp, c_t, CDOT_t, DDOT_t);
		//this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);
		std::cout << "molar creation and destruction rates, subroutine ckcdc:" << std::endl;
		std::copy(CDOT_t, CDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(DDOT_t, DDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;

		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];
		//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.

		std::cout << "Pressure:\t" << Pressure << std::endl;
		double init = 0.0;
		std::cout << "molar concentration:\t" << std::accumulate(c_t, c_t + nkk, init) << std::endl;
		mechanism::kinetics::kfkr(&Pressure, &Temp, x_t, FWDR_t, REVR_t);
		std::cout << "forward and reverse reaction rates:" << std::endl;
		std::copy(FWDR_t, FWDR_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(REVR_t, REVR_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;

		//forward reaction rates constant, reverse reaction rates constant.
		double* RKFT_t = new double[nii]; double* RKRT_t = new double[nii];
		//Returns the forward and reverse reaction rates Constant for reactions given pressure and temperature(s).
		mechanism::kinetics::kfrt(&Pressure, &Temp, RKFT_t, RKRT_t);
		std::cout << "forward and reverse reaction rate constant:" << std::endl;
		std::copy(RKFT_t, RKFT_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(RKRT_t, RKRT_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;

		cout << endl << RKFT_t[5 - 1] << " " << FWDR_t[5 - 1] << " " << c_t[3 - 1] << " " << RKFT_t[5 - 1] * c_t[3 - 1] * (c_t[3 - 1] * 2.5 + c_t[1 - 1]) << endl;
		cout << endl << RKRT_t[14 - 1] << " " << REVR_t[14 - 1] << " " << c_t[3 - 1] << " " << RKRT_t[14 - 1] * c_t[1 - 1] * c_t[3 - 1]/**(c_t[3-1]+c_t[1-1])*/ << endl;

		cout << "nii: " << nii << endl;
		int I_t = 25;	double R_A = 0.0;
		mechanism::kinetics::raex(&I_t, &R_A);
		cout << R_A << endl;
		I_t = -25; R_A *= 5.0;
		mechanism::kinetics::raex(&I_t, &R_A);
		cout << R_A << endl;
		//Returns the forward and reverse reaction rates Constant for reactions given pressure and temperature(s).
		mechanism::kinetics::kfrt(&Pressure, &Temp, RKFT_t, RKRT_t);
		std::copy(RKFT_t, RKFT_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(RKRT_t, RKRT_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;

		//print pre-coefficient, activation energy, Temperature factor
		double *RA2_t = new double[nii];
		double *RB2_t = new double[nii];
		double *RE2_t = new double[nii];

		mechanism::kinetics::abe(RA2_t, RB2_t, RE2_t);
		cout << "RA: " << endl;
		std::copy(RA2_t, RA2_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		cout << "RB: " << endl;
		std::copy(RB2_t, RB2_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		cout << "RE: " << endl;
		std::copy(RE2_t, RE2_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;

		delete[] RA2_t; delete[] RB2_t; delete RE2_t;
		delete[] x_t; delete[] y_t;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
		delete[] RKFT_t; delete[] RKRT_t;

	}

#endif // __CHEMKIN_AVAILABLE_ && __USE_LOSDE_


	void superPropagator::w2f_pgt(std::string tag)
	{
		//time
		std::ofstream fout((this->cwd_pgt + std::string("/output/time_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < time_data_pgt.size(); ++i) {
			fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << time_data_pgt[i] << std::endl;
		}
		fout.clear(); fout.close();

		//temperature
		fout.open((this->cwd_pgt + std::string("/output/temperature_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < temperature_data_pgt.size(); ++i) {
			fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << temperature_data_pgt[i] << std::endl;
		}
		fout.clear(); fout.close();

		//pressure
		fout.open((this->cwd_pgt + std::string("/output/pressure_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < pressure_data_pgt.size(); ++i) {
			//fout<<pressure_data_pgt[time_j]/1013250<<"\t"<<pressure_data_pgt[time_j]<<std::endl;
			fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << pressure_data_pgt[i] << "," << pressure_data_pgt[i] / 1013250 << "," << pressure_data_pgt[i] / 1000000 << std::endl;
		}
		fout.clear(); fout.close();

		//concentration, 8 species
		fout.open((this->cwd_pgt + std::string("/output/concentration_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < concentration_data_pgt[0].size(); ++i) {
			for (size_t j = 0; j < concentration_data_pgt.size(); ++j) {
				fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << concentration_data_pgt[j][i];
				if (j < (concentration_data_pgt.size() - 1))
					fout << ",";
			}
			fout << endl;
		}
		fout.clear(); fout.close();

		////int_drc, cumulative destructive rate constant
		//fout.open((this->cwd_pgt + std::string("/output/int_drc_") + tag + std::string(".csv")).c_str());
		//for (size_t time_j = 0; time_j < spe_drc_int_data_pgt[0].size(); ++time_j) {
			//for (size_t j = 0; j < spe_drc_int_data_pgt.size(); ++j) {
				//fout << std::setprecision(std::numeric_limits<double>::max_digits10+1) << spe_drc_int_data_pgt[j][time_j];
				//if (j < spe_drc_int_data_pgt.size() - 1)
					//fout << ",";
			//}
			//fout << endl;
		//}
		//fout.clear(); fout.close();

		//destruction rate constant,pseudo first order destructive rate constant
		fout.open((this->cwd_pgt + std::string("/output/drc_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < spe_drc_data_pgt[0].size(); ++i) {
			for (size_t j = 0; j < spe_drc_data_pgt.size(); ++j) {
				fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << spe_drc_data_pgt[j][i];
				if (j < spe_drc_data_pgt.size() - 1)
					fout << ",";
			}
			fout << endl;
		}
		fout.clear(); fout.close();

		////species production rate
		//fout.open((this->cwd_pgt + std::string("/output/spe_production_rate_") + tag + std::string(".csv")).c_str());
		//std::cout << spe_production_rate_data_pgt.size() << std::endl;
		//std::cout << spe_production_rate_data_pgt[0].size() << std::endl;
		//for (size_t time_j = 0; time_j < spe_production_rate_data_pgt[0].size(); ++time_j) {
		//	for (size_t j = 0; j < spe_production_rate_data_pgt.size(); ++j) {
		//		fout << std::setprecision(std::numeric_limits<double>::max_digits10+1) << spe_production_rate_data_pgt[j][time_j];
		//		if (j < spe_production_rate_data_pgt.size() - 1)
		//			fout << ",";
		//	}
		//	fout << endl;
		//}
		//fout.clear(); fout.close();

		//reaction rate
		fout.open((this->cwd_pgt + std::string("/output/reaction_rate_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < reaction_rate_data_pgt[0].size(); ++i) {
			for (size_t j = 0; j < reaction_rate_data_pgt.size(); ++j) {
				fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << reaction_rate_data_pgt[j][i];
				if (j < reaction_rate_data_pgt.size() - 1)
					fout << ",";
			}
			fout << std::endl;
		}
		fout.clear(); fout.close();

		//	fout.open((this->cwd_pgt+std::string("/output/rate_const_pre_factor_")+tag+std::string(".csv")).c_str());
		//	int I_t=1;	double R_A=0.0;
		//	//cout<<"\nPre-exponential Constant:"<<endl;
		//	for(;static_cast<int>(I_t)<=chemkin_cpp::chemkin::nii(); ++I_t){
		//		chemkin::ckraex(&I_t, &R_A);
		//		//cout<<I_t<<","<<R_A<<"\n";
		//		fout<<I_t<<","<<R_A<<"\n";
		//	}
		//	fout.clear(); fout.close();
	}


	void superPropagator::spe_concentration_w2f_pgt(double in_time, std::string tag) const
	{
		std::cout << in_time << std::endl;
		std::ofstream fout((this->cwd_pgt + std::string("/output/spe_concentration_") + tag + std::string(".csv")).c_str());
		for (std::size_t i = 0; i < concentration_data_pgt.size(); ++i) {
			fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << this->evaluate_concentration_at_time(in_time, i);
			if (i != concentration_data_pgt.size() - 1) {
				fout << ",";
			}
		}
		fout << std::endl;

		fout.close(); fout.clear();
	}

	void superPropagator::spe_concentration_w2f_pgt(std::vector<double> time_in, std::string tag) const
	{
		std::ofstream f_out((this->cwd_pgt + std::string("/output/spe_concentration_rate_") + tag + std::string(".csv")).c_str());

		for (size_t i = 0; i < time_in.size(); ++i) {
			for (size_t j = 0; j < concentration_data_pgt.size(); ++j) {
				f_out << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << evaluate_concentration_at_time(time_in[i], j);
				if (j < (concentration_data_pgt.size() - 1)) {
					f_out << ",";
				}
			}
			f_out << std::endl;
		}

		f_out.close();
	}

	void superPropagator::read_configuration(double & Temp, double & Pressure, std::size_t nkk, double x[])
	{
		//Initialize molar fractions.
		for (std::size_t i = 0; i < nkk; ++i) { x[i] = 0.0; }
		double pressure_atm = this->pgt_pt.get<double>("chem_init.pressure_atm");

		//The pressure has to be in cgs unit, convert pascal(Pa) to dyne/square centimeter [dyn/cm**2]
		//1atm= 101325 Pa
		//1Pa= 10 dyn/cm**2
		Pressure = pressure_atm * 1013250;
		Temp = this->pgt_pt.get<double>("chem_init.init_temperature");

		//read with json_parser as property_tree
		//nice and easy
		for (auto &key1 : this->pgt_pt.get_child("chem_init.species_index_concentration")) {
			x[boost::lexical_cast<std::size_t>(key1.first)] = key1.second.get_value<double>();
		}

		//normalized initial concentration
		if (this->pgt_pt.get<std::string>("propagator.normalize_initial_concentration") == "yes") {
			//renormalization
			double total_conc = 0.0;
			for (std::size_t i = 0; i < nkk; ++i) {
				total_conc += x[i];
			}//for
			if (total_conc > 0) {
				for (std::size_t i = 0; i < nkk; ++i) {
					x[i] /= total_conc;
				}//for
			}
		}//if

	}//read_configuration


#ifdef __LSODE_AVAILABLE_

	void superPropagator::initialize_lsode(rsp::my_time_t & dt, std::string outfile)
	{
		dt = this->pgt_pt.get<double>("lsode_init.dt");

		//Read 'lsode' parameters.
		lsodestore.atol = this->pgt_pt.get<double>("lsode_init.atol");
		lsodestore.rtol = this->pgt_pt.get<double>("lsode_init.rtol");
		lsodestore.mf = this->pgt_pt.get<int>("lsode_init.mf");
		lsodestore.jt = this->pgt_pt.get<int>("lsode_init.jt");
		lsodestore.itask = this->pgt_pt.get<int>("lsode_init.itask");
		lsodestore.iopt = this->pgt_pt.get<int>("lsode_init.iopt");
		lsodestore.itol = this->pgt_pt.get<int>("lsode_init.itol");
		lsodestore.istate = this->pgt_pt.get<int>("lsode_init.istate");

		lsodestore.deltaN1 = this->pgt_pt.get<int>("lsode_init.deltaN1");
		lsodestore.deltaN2 = this->pgt_pt.get<int>("lsode_init.deltaN2");

		int lrw_lsode, liw_lsode;
		std::ofstream fout((this->cwd_pgt + outfile).c_str(), std::ios::app);

		fout << std::endl;
		const int neq = mechanism::kinetics::nkk() + 1;
		if (lsodestore.mf == 10)
		{
			lrw_lsode = 20 + 16 * neq + 1;
			liw_lsode = 20 + 1;
			fout << "  lrw_lsode =   " << lrw_lsode << std::endl;
			fout << "  liw_lsode =   " << liw_lsode << std::endl;
		}
		else if ((lsodestore.mf == 21) || (lsodestore.mf == 22))
		{
			lrw_lsode = 22 + 9 * neq + neq * neq + 1;
			liw_lsode = 20 + neq + 1;
			fout << "  lrw_lsode =   " << lrw_lsode << std::endl;
			fout << "  liw_lsode =   " << liw_lsode << std::endl;
		}
		else
		{
			fout << "  mf = " << lsodestore.mf << std::endl;
			fout << "  The program should not use this 'mf' value. Please check 'lsode.f' for more information.\n" << std::endl;
		}
		fout << std::endl << "  Please make sure that the values of constants 'lrw' and 'liw' in the Fortran subroutine 'cpplsode'"
			"\n  are larger than 'lrw_lsode' and 'liw_lsode' respectively." << std::endl;

		auto Ischild = this->pgt_pt.get_child_optional("lsode_init.jt");
		if (Ischild) {
			fout << "\n" << "  you are using the 12 November 2003 version of DLSODA: Livermore Solver for Ordinary Differential Equations,\n" <<
				"  with Automatic method switching for stiff and nonstiff problems\n";
			/*
			 *  4223 C          If the RWORK length is to be fixed, it should be at least
			 *  4224 C               MAX (LRN, LRS),
			 *  4225 C          where LRN and LRS are the RWORK lengths required when the
			 *  4226 C          current method is nonstiff or stiff, respectively.
			 *  4227 C
			 *  4228 C          The separate RWORK length requirements LRN and LRS are
			 *  4229 C          as follows:
			 *  4230 C          IF NEQ is constant and the maximum method orders have
			 *  4231 C          their default values, then
			 *  4232 C             LRN = 20 + 16*NEQ,
			 *  4233 C             LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2,
			 *  4234 C             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5.
			 */
			double lrn_pgtsode, lrs_pgtsode;
			lrn_pgtsode = 20 + 16 * neq;
			fout << "  lrn_pgtsode =   " << lrn_pgtsode << std::endl;

			if ((lsodestore.jt == 1) || (lsodestore.jt == 2))
			{
				lrs_pgtsode = 22 + 9 * neq + neq * neq;
				fout << "  lrs_pgtsode =   " << lrs_pgtsode << std::endl;
				fout << std::endl << "  Please make sure that the values of constants 'lrw' is bigger than " << std::max(lrn_pgtsode, lrs_pgtsode) << std::endl;
			}
			else if ((lsodestore.jt == 4) || (lsodestore.jt == 5))
			{
				//lrs_pgtsode= 22+10*neq+neq*neq;
				fout << "  lrs_pgtsode =   " << "check your jt value" << std::endl;
			}
			else
			{
				fout << "  jt = " << lsodestore.jt << std::endl;
				fout << "  The program should not use this 'jt' value. Please check \"include/fortran_lib/dlsode/opkdmain.f\" for more information.\n" << std::endl;
			}
			fout << std::endl;
		}

		fout.close();

	}	//initialize_lsode

#endif // __LSODE_AVAILABLE_

	void superPropagator::initiate_propagator()
	{
#ifdef __LSODE_AVAILABLE_
		//initiate chemkin
		mechanism::kinetics::chemkin_init(
			(this->cwd_pgt + std::string("/input/chem.bin")).c_str(),
			(this->cwd_pgt + std::string("/output/general_output.csv")).c_str());
#endif // __LSODE_AVAILABLE_

#ifdef __USE_CANTERA_
		mechanism::kinetics::cantera_init();
#endif //__USE_CANTERA_


		//Read the file named "chem.out", read in the chemical reactions index
		rsp::relationshipParser::read_reactionNetwork_chemkin_index_map(reactionNetwork_chemkin_index_map, this->cwd_pgt + std::string("/input/chem.out"));

		//initialize the pointer
		initialize_cubic_spline_pointer();

		//read configuration file "setting.json"
		boost::property_tree::read_json(this->cwd_pgt + std::string("/input/setting.json"), pgt_pt, std::locale());

		//set fast reactions, read fast inter-conversion reaction pairs from "setting.json"
		update_chattering_spe_pair_from_file_pgt();

		//set the reaction rate of fast reactions to be zero
		//set_chattering_reaction_rates_to_zero_pgt();
	}


#ifdef __CHEMKIN_AVAILABLE_

	void propagator_sr::superPropagator::convert_molar_concentration_to_mole_fraction()
	{
		const int nkk = mechanism::kinetics::nkk();
		double* x_t = new double[nkk]; double* c_t = new double[nkk];
		for (std::size_t k = 0; k < time_data_pgt.size(); ++k) {
			for (int i = 0; i < nkk; ++i) {
				c_t[i] = concentration_data_pgt[i][k];
			}
			mechanism::kinetics::ctx(c_t, x_t);
			//convert molar concentration to mole fraction
			for (int i = 0; i < nkk; ++i) {
				concentration_data_pgt[i][k] = x_t[i];
			}
		}
	}

	void propagator_sr::superPropagator::convert_mole_fraction_to_molar_concentration()
	{
		const int nkk = mechanism::kinetics::nkk();
		double* x_t = new double[nkk]; double* c_t = new double[nkk];
		double pressure, temperature;

		for (std::size_t k = 0; k < time_data_pgt.size(); ++k) {
			pressure = this->pressure_data_pgt[k];
			temperature = this->temperature_data_pgt[k];
			double total = 0.0;
			for (int i = 0; i < nkk; ++i) {
				x_t[i] = concentration_data_pgt[i][k];
				total += x_t[i];
			}
			//normalization of mole fraction
			for (int i = 0; i < nkk; ++i)
				x_t[i] /= total;
			//convert mole fraction to molar concentration
			mechanism::kinetics::xtcp(&pressure, &temperature, x_t, c_t);
			for (int i = 0; i < nkk; ++i) {
				concentration_data_pgt[i][k] = c_t[i];
			}
}
	}

#endif // __CHEMKIN_AVAILABLE_


	void superPropagator::integrate_propensity_function_pgt()
	{

		spe_drc_int_data_pgt.resize(spe_drc_data_pgt.size());
		std::copy(spe_drc_data_pgt.begin(), spe_drc_data_pgt.end(), spe_drc_int_data_pgt.begin());


		for (size_t i = 0; i < spe_drc_int_data_pgt.size(); ++i)
		{
			//The first time interval
			//spe_drc_int_data_pgt[time_j][0] = spe_drc_data_pgt[time_j][0] * (time_data_pgt[1] - time_data_pgt[0]);
			spe_drc_int_data_pgt[i][0] = 0.0;
		}
		//The other time interval
		for (size_t i = 0; i < spe_drc_int_data_pgt.size(); ++i) {//[for
			for (size_t j = 1; j < spe_drc_int_data_pgt[0].size(); ++j) {
				////rectangle rule
				//spe_drc_int_data_pgt[time_j][j] = spe_drc_data_pgt[time_j][j] * (time_data_pgt[j] - time_data_pgt[j - 1]) + spe_drc_int_data_pgt[time_j][j - 1];

				//trapezoidal rule
				spe_drc_int_data_pgt[i][j] = 0.5 * (spe_drc_data_pgt[i][j] + spe_drc_data_pgt[i][j - 1]) * (time_data_pgt[j] - time_data_pgt[j - 1]) + spe_drc_int_data_pgt[i][j - 1];

			}
		}//for]

	}//int_propensity_function_pgt()


	bool superPropagator::init_spe_drc_pgt()
	{
		for (std::size_t i = 0; i < spe_drc_pgt.size(); ++i) {
			if (spe_drc_pgt[i] != 0)
				delete spe_drc_pgt[i];
		}
		//Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		spe_drc_pgt.clear();
		//create the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < spe_drc_data_pgt.size(); ++i)
		{
			spe_drc_pgt.push_back(new Linear_interp(time_data_pgt, spe_drc_data_pgt[i]));
		}
		return true;
	}//init_spe_drc_pgt()

	bool superPropagator::init_spe_drc_int_pgt()
	{
		for (std::size_t i = 0; i < spe_drc_int_pgt.size(); ++i) {
			if (spe_drc_int_pgt[i] != 0)
				delete spe_drc_int_pgt[i];
		}
		//Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		spe_drc_int_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < spe_drc_int_data_pgt.size(); ++i)
		{
			spe_drc_int_pgt.push_back(new Linear_interp(time_data_pgt, spe_drc_int_data_pgt[i]));
		}
		return true;
	}

	bool superPropagator::init_spe_drc_int_time_pgt()
	{
		for (std::size_t i = 0; i < spe_drc_int_time_pgt.size(); ++i) {
			if (spe_drc_int_time_pgt[i] != 0)
				delete spe_drc_int_time_pgt[i];
		}
		//Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		spe_drc_int_time_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < spe_drc_int_data_pgt.size(); ++i)
		{
			spe_drc_int_time_pgt.push_back(new Linear_interp(spe_drc_int_data_pgt[i], time_data_pgt));
		}
		return true;
	}

	bool superPropagator::init_reaction_rate_pgt()
	{
		for (std::size_t i = 0; i < reaction_rate_pgt.size(); ++i) {
			if (reaction_rate_pgt[i] != 0)
				delete reaction_rate_pgt[i];
		}
		reaction_rate_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i)
		{
			reaction_rate_pgt.push_back(new Linear_interp(time_data_pgt, reaction_rate_data_pgt[i]));
		}
		return true;
	}

	bool superPropagator::init_concentration_pgt()
	{
		for (size_t i = 0; i < concentration_pgt.size(); ++i) {
			if (concentration_pgt[i] != 0)
				delete concentration_pgt[i];
		}
		concentration_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> concentraton_pgt
		for (size_t i = 0; i < concentration_data_pgt.size(); ++i) {
			concentration_pgt.push_back(new Linear_interp(time_data_pgt, concentration_data_pgt[i]));
		}

		return true;
	}


	bool superPropagator::init_time_temperature_pgt()
	{
		if (time_temperature_pgt != 0)
			delete time_temperature_pgt;
		time_temperature_pgt = new Linear_interp(time_data_pgt, temperature_data_pgt);

		return true;
	}

	bool superPropagator::init_temperature_time_pgt()
	{
		if (temperature_time_pgt != 0)
			delete temperature_time_pgt;
		temperature_time_pgt = new Linear_interp(temperature_data_pgt, time_data_pgt);

		return true;
	}

	bool superPropagator::init_time_pressure_pgt()
	{
		if (time_pressure_pgt != 0)
			delete time_pressure_pgt;
		time_pressure_pgt = new Linear_interp(time_data_pgt, pressure_data_pgt);

		return true;
	}

	bool superPropagator::init_time_chattering_group_k_pgt()
	{
		for (size_t i = 0; i < chattering_group_k_data_pgt.size(); ++i) {
			if (chattering_group_k_pgt.size() == 0)
				break;
			if (chattering_group_k_pgt[i] != 0)
				delete chattering_group_k_pgt[i];
		}
		chattering_group_k_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> chattering_group_k_pgt
		for (size_t i = 0; i < chattering_group_k_data_pgt.size(); ++i) {
			chattering_group_k_pgt.push_back(new Linear_interp(time_data_pgt, chattering_group_k_data_pgt[i]));
		}

		return true;
	}

	void superPropagator::integrate_chattering_group_propensity_function_pgt()
	{
		chattering_group_k_int_data_pgt.resize(chattering_group_k_data_pgt.size());
		std::copy(chattering_group_k_data_pgt.begin(), chattering_group_k_data_pgt.end(), chattering_group_k_int_data_pgt.begin());


		for (size_t i = 0; i < chattering_group_k_int_data_pgt.size(); ++i)
		{
			//The first time interval
			//chattering_group_k_int_data_pgt[time_j][0] = chattering_group_k_data_pgt[time_j][0] * (time_data_pgt[1] - time_data_pgt[0]);
			chattering_group_k_int_data_pgt[i][0] = 0.0;
		}
		//The other time interval
		for (size_t i = 0; i < chattering_group_k_int_data_pgt.size(); ++i) {//[for
			for (size_t j = 1; j < chattering_group_k_int_data_pgt[0].size(); ++j) {
				////rectangle rule
				//chattering_group_k_int_data_pgt[time_j][j] = chattering_group_k_data_pgt[time_j][j] * (time_data_pgt[j] - time_data_pgt[j - 1]) + chattering_group_k_int_data_pgt[time_j][j - 1];

				//trapezoidal rule
				chattering_group_k_int_data_pgt[i][j] = 0.5 * (chattering_group_k_data_pgt[i][j] + chattering_group_k_data_pgt[i][j - 1]) * (time_data_pgt[j] - time_data_pgt[j - 1]) + chattering_group_k_int_data_pgt[i][j - 1];

			}
		}//for]

	}

	bool superPropagator::init_time_chattering_group_ss_prob_pgt()
	{
		for (size_t i = 0; i < chattering_group_ss_prob_data_pgt.size(); ++i) {
			if (chattering_group_ss_prob_pgt.size() == 0)
				break;
			if (chattering_group_ss_prob_pgt[i] != 0)
				delete chattering_group_ss_prob_pgt[i];
		}
		chattering_group_ss_prob_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> chattering_group_ss_prob_pgt
		for (size_t i = 0; i < chattering_group_ss_prob_data_pgt.size(); ++i) {
			chattering_group_ss_prob_pgt.push_back(new Linear_interp(time_data_pgt, chattering_group_ss_prob_data_pgt[i]));
		}

		return true;
	}

	bool superPropagator::init_chattering_group_k_int_pgt()
	{
		for (std::size_t i = 0; i < chattering_group_k_int_pgt.size(); ++i) {
			if (chattering_group_k_int_pgt[i] != 0)
				delete chattering_group_k_int_pgt[i];
		}
		//removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		chattering_group_k_int_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < chattering_group_k_int_data_pgt.size(); ++i)
		{
			chattering_group_k_int_pgt.push_back(new Linear_interp(time_data_pgt, chattering_group_k_int_data_pgt[i]));
		}
		return true;
	}

	bool superPropagator::init_chattering_group_k_int_time_pgt()
	{
		for (std::size_t i = 0; i < chattering_group_k_int_time_pgt.size(); ++i) {
			if (chattering_group_k_int_time_pgt[i] != 0)
				delete chattering_group_k_int_time_pgt[i];
		}
		//Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		chattering_group_k_int_time_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < chattering_group_k_int_data_pgt.size(); ++i)
		{
			chattering_group_k_int_time_pgt.push_back(new Linear_interp(chattering_group_k_int_data_pgt[i], time_data_pgt));
		}
		return true;
	}

	double superPropagator::evaluate_temperature_at_time(double in_time) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		return  time_temperature_pgt->interp(in_time);
	}

	double superPropagator::evaluate_time_at_temperature(double in_temperature) const
	{
		if (in_temperature >= temperature_data_pgt.back())
			in_temperature = temperature_data_pgt.back();
		return temperature_time_pgt->interp(in_temperature);
	}

	double superPropagator::evaluate_pressure_at_time(double in_time) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		return  time_pressure_pgt->interp(in_time);
	}

	double superPropagator::evaluate_concentration_at_time(double in_time, size_t index) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((index < 0) || (index >= concentration_pgt.size()))
			index = concentration_pgt.size() - 1;
		return concentration_pgt[index]->interp(in_time);
	}

	double superPropagator::evaluate_spe_drc_at_time(double in_time, size_t index) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((index < 0) || (index >= spe_drc_pgt.size()))
			index = spe_drc_pgt.size() - 1;
		return spe_drc_pgt[index]->interp(in_time);
	}

	double superPropagator::evaluate_spe_drc_int_at_time(double in_time, size_t index) const
	{
		if ((index < 0) || (index >= spe_drc_int_pgt.size()))
			index = spe_drc_int_pgt.size() - 1;

		double diff = in_time - time_data_pgt.back();
		if (diff > 0) {
			in_time = time_data_pgt.back();
			return spe_drc_int_pgt[index]->interp(in_time) + diff * spe_drc_data_pgt[index].back();
		}
		else {
			return spe_drc_int_pgt[index]->interp(in_time);
		}
	}

	double superPropagator::evaluate_time_at_spe_drc_int(double integral, size_t index) const
	{
		//the last 3 species's propensity function is zero in the time region for H2-O2 system
		//can not solve for time
		//if((index<0)||(index+3>= spe_drc_int_time_pgt.size()))
		//index=spe_drc_int_time_pgt.size()-1-3;
		if ((index < 0) || (index >= spe_drc_int_time_pgt.size())) {
			//std::cout<<"index out of bound!"<<std::endl;
			std::cerr << "index out of bound!" << std::endl;
		}
		double diff = integral - spe_drc_int_data_pgt[index].back();
		if (diff > 0) {
			integral = spe_drc_int_data_pgt[index].back();
			//return spe_drc_int_time_pgt[index]->interp(integral)+diff/spe_drc_data_pgt[index].back();
			return time_data_pgt.back() + diff / spe_drc_data_pgt[index].back();
		}
		else {
			return spe_drc_int_time_pgt[index]->interp(integral);
		}
	}

	double superPropagator::evaluate_reaction_rate_at_time(double in_time, size_t index) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((index < 0) || (index >= reaction_rate_pgt.size()))
			index = reaction_rate_pgt.size() - 1;
		return reaction_rate_pgt[index]->interp(in_time);
	}

	double superPropagator::evaluate_chattering_group_ss_prob_at_time(double in_time, size_t index) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((index < 0) || (index >= chattering_group_ss_prob_pgt.size()))
			index = chattering_group_ss_prob_pgt.size() - 1;
		return chattering_group_ss_prob_pgt[index]->interp(in_time);
	}

	double superPropagator::evaluate_chattering_group_k_at_time(double in_time, size_t chattering_group_id) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((chattering_group_id < 0) || (chattering_group_id >= chattering_group_k_pgt.size()))
			chattering_group_id = chattering_group_k_pgt.size() - 1;
		return chattering_group_k_pgt[chattering_group_id]->interp(in_time);
	}

	double superPropagator::evaluate_chattering_group_k_int_at_time(double in_time, size_t index) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((index < 0) || (index >= chattering_group_k_int_pgt.size()))
			index = chattering_group_k_int_pgt.size() - 1;
		return chattering_group_k_int_pgt[index]->interp(in_time);
	}

	double superPropagator::evaluate_time_at_chattering_group_k_int(double integral, size_t index) const
	{
		if ((index < 0) || (index >= chattering_group_k_int_time_pgt.size())) {
			std::cerr << "index out of bound!" << std::endl;
		}
		double diff = integral - chattering_group_k_int_data_pgt[index].back();
		if (diff > 0) {
			integral = chattering_group_k_int_data_pgt[index].back();
			return time_data_pgt.back() + diff / chattering_group_k_data_pgt[index].back();
		}
		else {
			return chattering_group_k_int_time_pgt[index]->interp(integral);
		}
	}



	}//namespace propagator_sr


#endif


