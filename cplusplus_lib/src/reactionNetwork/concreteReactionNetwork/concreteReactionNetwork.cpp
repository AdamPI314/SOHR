#ifndef __CONCRETEREACTIONNETWORK_CPP_
#define __CONCRETEREACTIONNETWORK_CPP_

#include <boost/property_tree/json_parser.hpp> //for json_reader
#include <boost/property_tree/xml_parser.hpp> //for write_xml

#include "../../../include/reactionNetwork/concreteReactionNetwork/concreteReactionNetwork.h"
#include "../../../include/tools/union_find/unionFind.h"

//infinitesimal dt
#define INFINITESIMAL_DT 1.0E-14

namespace reactionNetwork_sr {

	concreteReactionNetwork::concreteReactionNetwork(std::vector<double> uncertainties, std::size_t random_seed_for_this_core, std::string cwd_in)
		:superReactionNetwork(uncertainties, random_seed_for_this_core, cwd_in)
	{
		//can use polymorphism here, suppose we have multiple propagator
		//defer selection of propagator to meta data, dependency injection
		//actually we can use factory pattern here, we are not going to try it
		if (this->rnk_pt.get<std::string>("propagator.type") == std::string("dlsode") || this->rnk_pt.get<std::string>("propagator.type") == std::string("DLSODE")) {
			this->propagator = new pgt::dlsodePropagator(uncertainties, this->cwd);
			//set species initial concentration
			set_init_spe_concentration(sys_min_time);
		}
		else if (this->rnk_pt.get<std::string>("propagator.type") == std::string("sohr") || this->rnk_pt.get<std::string>("propagator.type") == std::string("SOHR")) {
			this->propagator = new pgt::SOHRPropagator(uncertainties, this->cwd);
		}
		else if (this->rnk_pt.get<std::string>("propagator.type") == std::string("ssa") || this->rnk_pt.get<std::string>("propagator.type") == std::string("SSA")) {
			this->propagator = new pgt::ssaPropagator(uncertainties, random_seed_for_this_core, this->cwd);
		}
		else {
			this->propagator = new pgt::superPropagator(uncertainties, this->cwd);
		}

		if (this->rnk_pt.get<std::string>("network.merge_chatterings") == std::string("yes")) {
			this->sp_chattering_rnk = this->propagator->get_sp_of_chattering();
			//std::cout << this->sp_chattering_rnk.use_count();
			this->merge_chatterings();
		}
		else {
			//not a good idea to leave smart pointer nullptr
			this->sp_chattering_rnk = std::make_shared<chattering_sr::chattering>();
		}

	}

	concreteReactionNetwork::~concreteReactionNetwork()
	{
		delete this->propagator;
	}

	void concreteReactionNetwork::merge_chatterings()
	{
		this->propagator->find_chattering_group(this->species_network_v);

		this->propagator->update_chattering_group_pairs_reactions(this->species_network_v, this->reaction_network_v, this->rnk_pt.get<std::string>("pathway.atom_followed"));

		this->propagator->subtract_chattering_reaction_contribution_from_species_drc_pgt();

		this->propagator->update_drc_and_equilibrium_probability_of_chattering_group();

		//set the reaction rate of fast reactions to be zero
		this->propagator->set_chattering_reaction_rates_to_zero_pgt();

		//re construct the cubic spline
		this->propagator->initiate_cubic_spline();

		//update_species_chattering_group_id
		this->update_species_chattering_group_id();

	}

	void concreteReactionNetwork::update_species_chattering_group_id()
	{
		for (auto x : this->sp_chattering_rnk->spe_idx_2_chattering_group_id_idx) {
			auto spe_idx = x.first;
			auto group_id_idx = x.second;
			this->species_network_v[spe_idx].chattering_group_id = group_id_idx.first;
		}
	}

	bool concreteReactionNetwork::set_init_spe_concentration(rsp::my_time_t in_time)
	{
		for (size_t i = 0; i < this->species_network_v.size(); ++i) {
			this->species_network_v[i].spe_conc = propagator->evaluate_concentration_at_time(in_time, i);
		}
		return true;
	}

	bool concreteReactionNetwork::update_reaction_rate(double in_time, vertex_t curr_vertex)
	{
		for (std::size_t i = 0; i < this->species_network_v[curr_vertex].reaction_k_index_s_coef_v.size(); ++i) {
			this->reaction_network_v[this->species_network_v[curr_vertex].reaction_k_index_s_coef_v[i].first].reaction_rate =
				propagator->evaluate_reaction_rate_at_time(in_time, this->species_network_v[curr_vertex].reaction_k_index_s_coef_v[i].first);
		}

		return true;
	}

	bool concreteReactionNetwork::update_reaction_rate(double in_time)
	{
		for (size_t i = 0; i < this->reaction_network_v.size(); ++i) {
			this->reaction_network_v[i].reaction_rate = propagator->evaluate_reaction_rate_at_time(in_time, i);
		}
		return true;
	}

	double concreteReactionNetwork::reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time, vertex_t curr_spe, double Y)
	{
		//if current species is a dead species, found
		if (this->dead_species.count(curr_spe) >= 1) {
			return std::numeric_limits<rsp::my_time_t>::max();
		}
		else {//not found
			  //the ln of 1.0/(1.0-Y)
			double ln_Y_reciprocal = log(1.0 / (1.0 - Y));
			//use the initial integral value at this time
			double init_spe_drc_int = propagator->evaluate_spe_drc_int_at_time(curr_time, curr_spe);
			//exact integral
			double exact_integral = ln_Y_reciprocal + init_spe_drc_int;
			//Solve for the first reaction time
			double reaction_time = propagator->evaluate_time_at_spe_drc_int(exact_integral, curr_spe);
			if (reaction_time < sys_min_time)
				return sys_min_time;
			else
				return reaction_time;
		}
	}

	double concreteReactionNetwork::reaction_time_from_importance_sampling(rsp::my_time_t curr_time, vertex_t curr_spe, double Y)
	{
		//if current species is a dead species, found
		if (this->dead_species.count(curr_spe) >= 1) {
			return std::numeric_limits<rsp::my_time_t>::max();
		}

		else {//not found
			  //the ln of 1.0/(1.0-Y)
			double ln_Y_reciprocal = log(1.0 / (1.0 - Y));
			//use the initial integral value at this time

			double init_spe_drc_int = propagator->evaluate_spe_drc_int_at_time(curr_time, curr_spe);

			//exact integral
			double exact_integral = ln_Y_reciprocal + init_spe_drc_int;
			//Solve for the first reaction time
			double reaction_time = propagator->evaluate_time_at_spe_drc_int(exact_integral, curr_spe);

			if (reaction_time < sys_min_time)
				return sys_min_time;
			else if (reaction_time > tau)
				return tau;
			else
				return reaction_time;
		}
	}

	double concreteReactionNetwork::chattering_group_reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time, vertex_t curr_group, double Y)
	{
		if (curr_time >= this->tau) {
			return curr_time;
		}
		else {
			//the ln of 1.0/(1.0-Y)
			double ln_Y_reciprocal = log(1.0 / (1.0 - Y));
			//use the initial integral value at this time
			double init_chattering_group_k_int = propagator->evaluate_chattering_group_k_int_at_time(curr_time, curr_group);
			//exact integral
			double exact_integral = ln_Y_reciprocal + init_chattering_group_k_int;
			//Solve for the first reaction time
			double reaction_time = propagator->evaluate_time_at_chattering_group_k_int(exact_integral, curr_group);

			if (reaction_time < sys_min_time)
				return sys_min_time;
			else
				return reaction_time;
		}
	}

	bool concreteReactionNetwork::set_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t)
	{
		//pro_max= 1-prob_min= 1-exp[-integrate_{init_time}^{end_time}{propensity function}];
		if (init_time <= end_time) {
			this->species_network_v[index_t].prob_max = 1.0 - exp(-(propagator->evaluate_spe_drc_int_at_time(end_time, index_t) - propagator->evaluate_spe_drc_int_at_time(init_time, index_t)));
		}
		else
			this->species_network_v[index_t].prob_max = 0.0;
		return true;
	}

	bool concreteReactionNetwork::set_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t)
	{
		//prob_min= exp[-integrate_{init_time}^{end_time}{propensity function}];
		if (init_time <= end_time)
			this->species_network_v[index_t].prob_min = exp(-(propagator->evaluate_spe_drc_int_at_time(end_time, index_t) - propagator->evaluate_spe_drc_int_at_time(init_time, index_t)));
		else
			this->species_network_v[index_t].prob_min = 1.0;
		return true;
	}

	double concreteReactionNetwork::get_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t) const
	{
		if (init_time >= end_time)
			return 0.0;
		else
			return 1.0 - exp(-(propagator->evaluate_spe_drc_int_at_time(end_time, index_t) - propagator->evaluate_spe_drc_int_at_time(init_time, index_t)));
	}

	double concreteReactionNetwork::get_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t) const
	{
		if (init_time >= end_time)
			return 1.0;
		else
			return exp(-(propagator->evaluate_spe_drc_int_at_time(end_time, index_t) - propagator->evaluate_spe_drc_int_at_time(init_time, index_t)));
	}

	double concreteReactionNetwork::evaluate_spe_concentration_at_time(double time, std::size_t index) const
	{
		return propagator->evaluate_concentration_at_time(time, index);
	}

	double concreteReactionNetwork::evaluate_spe_drc_at_time(double time, std::size_t index) const
	{
		return propagator->evaluate_spe_drc_at_time(time, index);
	}

	double concreteReactionNetwork::evaluate_chattering_group_ss_prob_at_time(double in_time, size_t index) const
	{
		return propagator->evaluate_chattering_group_ss_prob_at_time(in_time, index);
	}

	double concreteReactionNetwork::prob_chattering_group_will_react_in_a_time_range(double init_time, double pathway_end_time, size_t curr_chattering_group)
	{
		//pro_max= 1-prob_min= 1-exp[-integrate_{init_time}^{end_time}{propensity function}];
		double pro_max = 0.0;
		if (init_time <= pathway_end_time) {
			pro_max = 1.0 - exp(-(propagator->evaluate_chattering_group_k_int_at_time(pathway_end_time, curr_chattering_group) - propagator->evaluate_chattering_group_k_int_at_time(init_time, curr_chattering_group)));
		}
		else
			pro_max = 0.0;

		return pro_max;
	}

	double concreteReactionNetwork::pathway_AT_sim_move_one_step(double when_time, size_t curr_spe, size_t next_reaction, size_t next_spe)
	{
		//update rate in the reaction network
		update_reaction_rate(when_time, curr_spe);

		//generage the random number u_1
		double u_1 = 0.0;
		do {
			u_1 = rand->random01();
		} while (u_1 == 1.0);

		//when_time= reaction_time_from_importance_sampling(when_time, curr_chattering_group, u_1);
		when_time = reaction_time_from_importance_sampling_without_cutoff(when_time, curr_spe, u_1);
		return when_time;
	}

	double concreteReactionNetwork::pathway_AT_input_pathway_sim_once(double init_time, double pathway_end_time, std::string pathway_in)
	{
		//parse pathway to spe and reaction vector
		std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;
		this->parse_pathway_to_vector(pathway_in, spe_vec, reaction_vec);

		//set pathway end time
		set_tau(pathway_end_time);

		double when_time = init_time;

		//start from the first reaction
		for (size_t i = 0; i < reaction_vec.size(); ++i)
		{
			when_time = pathway_AT_sim_move_one_step(when_time, spe_vec[i], reaction_vec[i], spe_vec[i + 1]);
		}

		return when_time;
	}


	double concreteReactionNetwork::pathway_prob_input_pathway_recursive_relation(double tau_j_minus_1, std::vector<std::size_t> spe_vec, std::vector<std::size_t> reaction_vec, std::vector<std::size_t> N_subvolume, std::size_t j_th)
	{
		/*
		* assume spe_vec.size()= reaction_vec.size()+1
		* {{I^{(M)}}({\tau ^{(M)}},{S^{(M - 1)}},{R^{(M)}},{S^{(M)}},{N^{(M)}})
		* = exp\left[ { - \int_{\tau _i^{(M)}}^{{\tau _f}} {a(\tau '|t + \tau _i^{(M)};{S^{(M)}})} d\tau '} \right]}
		* if j_th==M, the last integral, or says no reaction after current species
		* our index start from j_th=0
		*/
		if (j_th == reaction_vec.size()) {
			//std::cout<<"j_th:\t"<<j_th<<"\t"<<get_spe_prob_min_at_a_time(tau_j_minus_1, pathway_end_time, spe_vec.back())<<std::endl;
			return get_spe_prob_min_at_a_time(tau_j_minus_1, tau, spe_vec.back());
		}
		else {
			/*
			* rectangle rule, trapezoidal rule are all Newton–Cotes formulas, use simpson's rule here
			* http://nbviewer.ipython.org/github/ResearchComputing/HPSC-Fall-2013/blob/master/lab/lab-08-integrate/simpsons.ipynb
			* # From Wikipeia
			* def simpson_iter(f, a, b, n):
			*     """Approximates the definite integral of f from a to b by
			*      the composite Simpson's rule, using n subintervals"""
			*     h = (b - a) / float(n)
			*     s = f(a) + f(b)
			*     for i in range(1, n, 2):
			*     s += 4 * f(a + i * h)
			*     for i in range(2, n-1, 2):
			*     s += 2 * f(a + i * h)

			*     return s * h / 3
			*/

			//need to do the transformation essential to importance sampling
			//be careful, the prob_spe_will_react_in_a_time_range is taken cared in this integral range!!!
			double a = 0.0, b = get_spe_prob_max_at_a_time(tau_j_minus_1, tau, spe_vec[j_th]);
			double h = (b - a) / N_subvolume[j_th];

			//basically says reaction occur at time_a and time_b immediately, so the next recursive relation start from this time point a
			double tau_j = reaction_time_from_importance_sampling(tau_j_minus_1, spe_vec[j_th], a);
			if (tau_j <= sys_min_time)
				tau_j = sys_min_time;

			//double s= prob_spe_will_react_in_a_time_range(tau_j_minus_1, pathway_end_time, spe_vec[j_th])*
			double s = reaction_spe_branching_ratio(tau_j, spe_vec[j_th], reaction_vec[j_th], spe_vec[j_th + 1])*
				pathway_prob_input_pathway_recursive_relation(tau_j, spe_vec, reaction_vec, N_subvolume, j_th + 1);

			//b
			tau_j = reaction_time_from_importance_sampling(tau_j_minus_1, spe_vec[j_th], b);
			//s+= prob_spe_will_react_in_a_time_range(tau_j_minus_1, pathway_end_time, spe_vec[j_th])*
			s += reaction_spe_branching_ratio(tau_j, spe_vec[j_th], reaction_vec[j_th], spe_vec[j_th + 1])*
				pathway_prob_input_pathway_recursive_relation(tau_j, spe_vec, reaction_vec, N_subvolume, j_th + 1);

			for (size_t i = 1; i < N_subvolume[j_th]; i += 2) {
				tau_j = reaction_time_from_importance_sampling(tau_j_minus_1, spe_vec[j_th], a + i*h);
				s += 4 *
					//prob_spe_will_react_in_a_time_range(tau_j_minus_1, pathway_end_time, spe_vec[j_th])*
					reaction_spe_branching_ratio(tau_j, spe_vec[j_th], reaction_vec[j_th], spe_vec[j_th + 1])*
					pathway_prob_input_pathway_recursive_relation(tau_j, spe_vec, reaction_vec, N_subvolume, j_th + 1);
			}
			for (size_t i = 2; i < N_subvolume[j_th] - 1; i += 2) {
				tau_j = reaction_time_from_importance_sampling(tau_j_minus_1, spe_vec[j_th], a + i*h);
				s += 2 *
					//prob_spe_will_react_in_a_time_range(tau_j_minus_1, pathway_end_time, spe_vec[j_th])*
					reaction_spe_branching_ratio(tau_j, spe_vec[j_th], reaction_vec[j_th], spe_vec[j_th + 1])*
					pathway_prob_input_pathway_recursive_relation(tau_j, spe_vec, reaction_vec, N_subvolume, j_th + 1);
			}
			return s*h / 3;
		}
	}


	/*
	* out flux of a species at a specific time
	* return a vector<species index, flux>
	* Note that a local class/struct is not allowed as predicate to find_if, got to define it outside the local function or class
	*/
	struct comp {
		std::size_t spe_index_;
		comp(std::size_t const &spe_index_in) :spe_index_(spe_index_in) {};

		bool operator()(spe_index_flux_t const &p) {
			return p.first == spe_index_;
		}
	};

	std::vector<spe_index_flux_t> concreteReactionNetwork::spe_out_flux_at_a_time(rsp::my_time_t in_time, std::size_t curr_spe, std::string atom_followed)
	{
		this->update_reaction_rate(in_time, curr_spe);

		std::vector<spe_index_flux_t> out_flux;
		//search all out reaction of a species
		for (std::size_t i = 0; i < this->species_network_v[curr_spe].reaction_k_index_s_coef_v.size(); ++i) {//for
																											//search all out species of a reaction
			for (std::size_t j = 0; j < this->reaction_network_v[this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first].out_spe_index_weight_v_map[atom_followed].size(); ++j) {
				//if current pointing species is not in "out_flux"
				std::vector<spe_index_flux_t>::iterator ite = std::find_if(out_flux.begin(), out_flux.end(),
					comp(this->reaction_network_v[this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first].out_spe_index_weight_v_map[atom_followed][j].first));
				if (ite != out_flux.end()) {
					//find it
					//std::cout<<"find it!"<<std::endl;
					ite->second +=/*this->species_network_v[curr_chattering_group].reaction_k_index_s_coef_v[i].second* //s_coef_product*/
						this->reaction_network_v[this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first].reaction_rate; //reaction rate
				}
				else {
					//not find it
					out_flux.push_back(std::make_pair(this->reaction_network_v[this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first].out_spe_index_weight_v_map[atom_followed][j].first,
						this->reaction_network_v[this->species_network_v[curr_spe].reaction_k_index_s_coef_v[i].first].reaction_rate));
				}
			}
		}

		////normalize the flux to get flux ratio
		//double total_flux= std::accumulate(out_flux.begin(), out_flux.end(), 0.0, myplus());
		////std::cout<<total_flux<<std::endl;
		//for(std::size_t i=0; i<out_flux.size(); ++i){
		//	out_flux[i].second/=total_flux;
		//}

		return out_flux;
	}

	/*
	* calculate the flux of all species
	* normalize them
	* print to file
	*/
	typedef std::pair<std::size_t, std::size_t> spe_index_pair_t;
	typedef std::pair<spe_index_pair_t, double> spe_index_pair_flux_t;

	struct comp2 {
		spe_index_pair_t lhs;
		comp2(spe_index_pair_t const &in_lhs) :lhs(in_lhs) {};
		bool operator()(spe_index_pair_flux_t const &rhs) const {
			return (lhs == rhs.first);
		}
	};

	//Almost. The functor must be a binary operator taking the return value type as the first and the range type as the second argument
	//http://stackoverflow.com/questions/6935118/how-to-apply-stdaccumulate-algorithm-for-associative-containers
	struct myplus2 {
		double operator()(double const &lhs, spe_index_pair_flux_t const &rhs) const {
			return lhs + rhs.second;
		}
	};

	void concreteReactionNetwork::flux_analysis_w2f(rsp::my_time_t in_time)
	{
		//get forward and backward flux separately
		//A->B, B->A
		std::vector<spe_index_pair_flux_t> spe_pair_flux_v;

		std::vector<spe_index_flux_t> out_flux_tmp;
		for (std::size_t k = 0; k < get_num_vertices(); ++k) {
			out_flux_tmp = spe_out_flux_at_a_time(in_time, k);
			for (std::size_t i = 0; i < out_flux_tmp.size(); ++i) {
				//both contain "H" atom
				if ((this->species_network_v[k].spe_component["H"] != 0) && (this->species_network_v[out_flux_tmp[i].first].spe_component["H"] != 0)) {
					spe_pair_flux_v.push_back(std::make_pair(std::make_pair(k, out_flux_tmp[i].first), out_flux_tmp[i].second));
				}
			}
		}

		//combine forward and backward flux
		//if flux(A->B)> flux(B->A), keep A->B, let flux_new(A->B)= flux(A->B)- flux(B->A)
		//and vice-versa
		std::vector<spe_index_pair_flux_t> spe_pair_flux_v2;
		for (std::size_t i = 0; i < get_num_vertices(); ++i) {
			for (std::size_t j = 0; j < i; ++j) {
				std::vector<spe_index_pair_flux_t>::iterator ite1 = std::find_if(spe_pair_flux_v.begin(), spe_pair_flux_v.end(), comp2(std::make_pair(i, j)));
				std::vector<spe_index_pair_flux_t>::iterator ite2 = std::find_if(spe_pair_flux_v.begin(), spe_pair_flux_v.end(), comp2(std::make_pair(j, i)));
				if (ite1 != spe_pair_flux_v.end()) {
					if (ite2 != spe_pair_flux_v.end()) {
						if (ite1->second >= ite2->second) {
							spe_pair_flux_v2.push_back(std::make_pair(std::make_pair(i, j), ite1->second - ite2->second));
						}
						else {
							spe_pair_flux_v2.push_back(std::make_pair(std::make_pair(j, i), ite2->second - ite1->second));
						}
					}
				}

			}
		}

		std::cout << "num of edges:\t" << spe_pair_flux_v2.size() << std::endl;
		for (std::size_t i = 0; i < spe_pair_flux_v2.size(); ++i) {
			std::cout << spe_pair_flux_v2[i].first.first << "\t" << spe_pair_flux_v2[i].first.second << "\t" << spe_pair_flux_v2[i].second << std::endl;
		}

		//re-normalize the out flux of a species, basically calculate the out flux ratio
		std::vector<std::size_t> first_ele;
		std::vector< std::vector<spe_index_pair_flux_t> > spe_index_pair_flux_v_v;
		for (std::size_t i = 0; i < spe_pair_flux_v2.size(); ++i) {
			std::vector<std::size_t>::iterator ite = find(first_ele.begin(), first_ele.end(), spe_pair_flux_v2[i].first.first);

			if (ite == first_ele.end()) {
				//not in the tag vector, need create a new vector of spe_index_pair_flux_t
				std::vector<spe_index_pair_flux_t> spe_index_pair_flux_v;
				spe_index_pair_flux_v.push_back(spe_pair_flux_v2[i]);
				spe_index_pair_flux_v_v.push_back(spe_index_pair_flux_v);
				first_ele.push_back(spe_pair_flux_v2[i].first.first);
			}
			else {
				//found in the tag vector
				std::size_t index = std::distance(first_ele.begin(), ite);
				spe_index_pair_flux_v_v[index].push_back(spe_pair_flux_v2[i]);
			}

		}

		std::cout << "num of classes:\t" << first_ele.size() << std::endl;

		for (std::size_t i = 0; i < spe_index_pair_flux_v_v.size(); ++i) {
			for (std::size_t j = 0; j < spe_index_pair_flux_v_v[i].size(); ++j) {
				std::cout << spe_index_pair_flux_v_v[i][j].first.first << "\t" << spe_index_pair_flux_v_v[i][j].first.second << "\t" << spe_index_pair_flux_v_v[i][j].second << endl;
			}
		}

		//re-normalize
		for (std::size_t i = 0; i < spe_index_pair_flux_v_v.size(); ++i) {
			double total_flux = std::accumulate(spe_index_pair_flux_v_v[i].begin(), spe_index_pair_flux_v_v[i].end(), 0.0, myplus2());
			for (std::size_t j = 0; j < spe_index_pair_flux_v_v[i].size(); ++j) {
				spe_index_pair_flux_v_v[i][j].second /= total_flux;
			}
		}

		std::cout << "normalized result:\n";

		for (std::size_t i = 0; i < spe_index_pair_flux_v_v.size(); ++i) {
			for (std::size_t j = 0; j < spe_index_pair_flux_v_v[i].size(); ++j) {
				std::cout << spe_index_pair_flux_v_v[i][j].first.first << "\t" << spe_index_pair_flux_v_v[i][j].first.second << "\t" << spe_index_pair_flux_v_v[i][j].second << endl;
			}
		}

		/*
		* generate Mathematica Format graph input file
		* DirectedEdge["A", "B"]
		*/

		std::ofstream fout((this->cwd + std::string("/output/H2_O2_graph.csv")).c_str());

		for (std::size_t i = 0; i < spe_index_pair_flux_v_v.size(); ++i) {
			for (std::size_t j = 0; j < spe_index_pair_flux_v_v[i].size(); ++j) {
				fout << "DirectedEdge[\"" << species_network_v[spe_index_pair_flux_v_v[i][j].first.first].spe_name
					<< "\", \"" << species_network_v[spe_index_pair_flux_v_v[i][j].first.second].spe_name << "\"]" << std::endl;
			}
		}
		fout.clear(); fout.close();

		std::ofstream fout2((this->cwd + std::string("/output/H2_O2_graph_edge_weight.csv")).c_str());

		for (std::size_t i = 0; i < spe_index_pair_flux_v_v.size(); ++i) {
			for (std::size_t j = 0; j < spe_index_pair_flux_v_v[i].size(); ++j) {
				fout2 << spe_index_pair_flux_v_v[i][j].second << std::endl;
			}
		}
		fout2.clear(); fout2.close();
	}

	void concreteReactionNetwork::spe_concentration_w2f_rnk(double in_time, std::string str) const
	{
		propagator->spe_concentration_w2f_pgt(in_time, str);
	}

}/*namespace reactionNetwork_sr*/

#endif
