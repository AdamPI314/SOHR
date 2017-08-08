#ifndef __REACTIONNETWORKODESOLVER_CPP_
#define __REACTIONNETWORKODESOLVER_CPP_

#include "reactionNetworkODESolver.h"
#include "../../pathwayHandler/pathwayHandler.h"
#include "../../../include/tools/my_debug/my_debug.h"
//infinitesimal dt
#define INFINITESIMAL_DT 1.0E-14
#define PRINT_PRECISION 15


namespace reactionNetworkODESolver_sr {
	reactionNetworkODESolver::reactionNetworkODESolver(std::vector<double> uncertainties,
		std::size_t random_seed_for_this_core, std::string cwd_in)
		: superReactionNetwork(uncertainties, random_seed_for_this_core, cwd_in),
		SOHRPropagator(uncertainties, cwd_in) {

		this->update_number_of_ways_making_species();

		this->single_source = this->pgt_pt.get<int>("SOHR_init.single_source_species");
		// if there is a single source, change the time-dependent single correction factor
		if (this->single_source >= 0) {
			this->single_source_spe_dr_rnos.assign(this->time_data_pgt.size(), 0.0);
			this->single_source_additional_concentration_data_rnos.assign(this->time_data_pgt.size(), 0.0);

			this->update_single_source_additional_concentration_data_rnos_s_ct_np();
			this->init_time_single_source_additional_concentration_pointer_rnos();
		}

	}

	reactionNetworkODESolver::~reactionNetworkODESolver() {
		if (this->time_single_source_additional_concentration_pointer_rnos != nullptr)
			delete this->time_single_source_additional_concentration_pointer_rnos;
	}

	//return number of species and number of time points
	std::pair<std::size_t, std::size_t> reactionNetworkODESolver::get_size_of_concentration_data() const {
		return std::make_pair(this->concentration_data_pgt.size(), this->concentration_data_pgt[0].size());
	}

	std::vector<double> reactionNetworkODESolver::get_initial_concentration() const
	{
		//probability
		std::vector<double> prob(this->species_network_v.size());
		for (std::size_t i = 0; i < prob.size(); ++i) {
			prob[i] = species_network_v[i].spe_conc;
		}
		return prob;
	}

	bool reactionNetworkODESolver::set_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t) {
		//pro_max= 1-prob_min= 1-exp[-integrate_{init_time}^{end_time}{propensity function}];
		if (init_time <= end_time) {
			this->species_network_v[index_t].prob_max = 1.0 -
				exp(-(this->evaluate_spe_drc_int_at_time(end_time, index_t) -
					this->evaluate_spe_drc_int_at_time(init_time, index_t)));
		}
		else
			this->species_network_v[index_t].prob_max = 0.0;
		return true;
	}

	bool reactionNetworkODESolver::set_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t) {
		//prob_min= exp[-integrate_{init_time}^{end_time}{propensity function}];
		if (init_time <= end_time)
			this->species_network_v[index_t].prob_min = exp(-(this->evaluate_spe_drc_int_at_time(end_time, index_t) -
				this->evaluate_spe_drc_int_at_time(init_time, index_t)));
		else
			this->species_network_v[index_t].prob_min = 1.0;
		return true;
	}

	double reactionNetworkODESolver::get_spe_prob_max_at_a_time(double init_time, double end_time,
		size_t index_t) const {
		if (init_time >= end_time)
			return 0.0;
		else
			return 1.0 - exp(-(this->evaluate_spe_drc_int_at_time(end_time, index_t) -
				this->evaluate_spe_drc_int_at_time(init_time, index_t)));
	}

	double reactionNetworkODESolver::get_spe_prob_min_at_a_time(double init_time, double end_time,
		size_t index_t) const {
		if (init_time >= end_time)
			return 1.0;
		else
			return exp(-(this->evaluate_spe_drc_int_at_time(end_time, index_t) -
				this->evaluate_spe_drc_int_at_time(init_time, index_t)));
	}

	bool reactionNetworkODESolver::update_reaction_rate(double in_time, rnk::vertex_t curr_vertex) {
		for (std::size_t i = 0; i < this->species_network_v[curr_vertex].reaction_k_index_s_coef_v.size(); ++i) {
			this->reaction_network_v[this->species_network_v[curr_vertex].reaction_k_index_s_coef_v[i].first].reaction_rate =
				this->evaluate_reaction_rate_at_time(in_time,
					this->species_network_v[curr_vertex].reaction_k_index_s_coef_v[i].first);
		}

		return true;
	}

	bool reactionNetworkODESolver::update_reaction_rate(double in_time) {
		for (size_t i = 0; i < this->reaction_network_v.size(); ++i) {
			this->reaction_network_v[i].reaction_rate = this->evaluate_reaction_rate_at_time(in_time, i);
		}
		return true;
	}

	double reactionNetworkODESolver::reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time,
		rnk::vertex_t curr_spe,
		double Y) {
		//if current species is a dead species, found
		if (std::find(this->dead_species.begin(), this->dead_species.end(), curr_spe) != this->dead_species.end()) {
			return std::numeric_limits<rsp::my_time_t>::max();
		}
		else {//not found
			//the ln of 1.0/(1.0-Y)
			double ln_Y_reciprocal = log(1.0 / (1.0 - Y));
			//use the initial integral value at this time
			double init_spe_drc_int = this->evaluate_spe_drc_int_at_time(curr_time, curr_spe);
			//exact integral
			double exact_integral = ln_Y_reciprocal + init_spe_drc_int;
			//Solve for the first reaction time
			double reaction_time = this->evaluate_time_at_spe_drc_int(exact_integral, curr_spe);
			if (reaction_time < sys_min_time)
				return sys_min_time;
			else
				return reaction_time;
		}
	}

	double reactionNetworkODESolver::reaction_time_from_importance_sampling(rsp::my_time_t curr_time,
		rnk::vertex_t curr_spe, double Y) {
		//if current species is a dead species, found
		if (std::find(this->dead_species.begin(), this->dead_species.end(), curr_spe) != this->dead_species.end()) {
			return std::numeric_limits<rsp::my_time_t>::max();
		}

		else {
			//not found
			//the ln of 1.0/(1.0-Y)
			double ln_Y_reciprocal = log(1.0 / (1.0 - Y));
			//use the initial integral value at this time

			double init_spe_drc_int = this->evaluate_spe_drc_int_at_time(curr_time, curr_spe);

			//exact integral
			double exact_integral = ln_Y_reciprocal + init_spe_drc_int;
			//Solve for the first reaction time
			double reaction_time = this->evaluate_time_at_spe_drc_int(exact_integral, curr_spe);

			if (reaction_time < sys_min_time)
				return sys_min_time;
			else if (reaction_time > path_end_time)
				return path_end_time;
			else
				return reaction_time;
		}
	}

	void reactionNetworkODESolver::set_concentration_at_time_zero_to_initial_fraction_or_concentration() {
		for (std::size_t i = 0; i < concentration_data_pgt.size(); ++i) {
			concentration_data_pgt[i][0] = this->species_network_v[i].spe_conc;
		}
	}

	void reactionNetworkODESolver::set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(std::vector<std::vector<double> > & prob_Mat) {
		for (std::size_t i = 0; i < this->species_network_v.size(); ++i)
			prob_Mat[i][0] = this->species_network_v[i].spe_conc;
	}

	void reactionNetworkODESolver::set_probability_matrix_of_a_species_constant(std::vector<std::vector<double>>& prob_Mat, std::size_t ind)
	{
		std::fill(prob_Mat[ind].begin(), prob_Mat[ind].end(), this->species_network_v[ind].spe_conc);
	}

	void reactionNetworkODESolver::set_probability_matrix_of_single_source_constant(std::vector<std::vector<double>>& prob_Mat)
	{
		if (this->single_source < 0)
			return;

		set_probability_matrix_of_a_species_constant(prob_Mat, this->single_source);
	}

	void reactionNetworkODESolver::update_single_source_additional_concentration_data_rnos_s_ct_np()
	{
		if (this->single_source < 0)
			return;

		update_single_source_species_dr_based_on_spe_concentration_s_ct_np();

		// don't include initial concentration of source, renormalize among other non-source species
		single_source_additional_concentration_data_rnos[0] = 0.0;

		double Npoints = this->time_data_pgt.size();
		double delta_t = this->time_data_pgt.back() / (double)(Npoints - 1);

		//// Rectangle rule based numerical integral
		//for (std::size_t i = 1; i < single_source_additional_concentration_data_rnos.size(); ++i) {
		//single_source_additional_concentration_data_rnos[i] = this->single_source_spe_dr_rnos[i] * delta_t;
		//}

		// Trapezoidal rule
		for (std::size_t ti = 1; ti < single_source_additional_concentration_data_rnos.size(); ++ti) {
			single_source_additional_concentration_data_rnos[ti] = 0.5 * (this->single_source_spe_dr_rnos[ti - 1] + this->single_source_spe_dr_rnos[ti]) * delta_t;
		}
	}

	void reactionNetworkODESolver::init_time_single_source_additional_concentration_pointer_rnos()
	{
		if (this->single_source < 0)
			return;

		if (time_single_source_additional_concentration_pointer_rnos != nullptr)
			delete time_single_source_additional_concentration_pointer_rnos;
		time_single_source_additional_concentration_pointer_rnos = new Linear_interp(time_data_pgt, single_source_additional_concentration_data_rnos);
	}

	double reactionNetworkODESolver::evaluate_single_source_additional_concentration_at_time(double in_time) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		return  time_single_source_additional_concentration_pointer_rnos->interp(in_time);
	}

	void reactionNetworkODESolver::update_single_source_species_dr_based_on_spe_concentration_s_ct_np()
	{
		if (this->single_source < 0)
			return;

		//in which nkk is number of species
		const int nkk = this->species_network_v.size();

		//molar concentration.
		double *c_t = new double[nkk];
		for (int si = 0; si < nkk; ++si) { c_t[si] = 0.0; }

		//species creation rates and destruction rates.
		double *CDOT_t = new double[nkk];
		double *DDOT_t = new double[nkk];

		//temperature, dummy variable in Lotka-Volterra model
		double Temp = 298.0;

		for (std::size_t tk = 0; tk < this->time_data_pgt.size(); ++tk) {
			for (int si = 0; si < nkk; ++si) {
				c_t[si] = concentration_data_pgt[si][tk];//ith species, kth point
			}
			//The difference is how to treat the auto-catylytic reactions
			//chemkincpp_sr::chemkin::ckcdc(&Temp, c_t, CDOT_t, DDOT_t);
			this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);
			this->single_source_spe_dr_rnos[tk] = DDOT_t[this->single_source];
		}

		delete[] c_t;
		delete[] CDOT_t;
		delete[] DDOT_t;
	}


	void reactionNetworkODESolver::check_zero_concentration(double deltaConcentration)
	{
		for (std::size_t i = 0; i < concentration_data_pgt.size(); ++i) {
			// not the initial concentration
			for (std::size_t j = 1; j < concentration_data_pgt[i].size(); ++j) {
				if (concentration_data_pgt[i][j] == 0)
					concentration_data_pgt[i][j] = deltaConcentration;
			}
		}
	}

	void reactionNetworkODESolver::divide_concentration_by_number_of_ways_making_it()
	{
		// find available elements, based on initial concentrations
		std::vector<bool> element_available(this->element_v.size(), false);
		for (auto x : this->species_network_v) {
			// initial concentration greater than zero
			if (x.spe_conc > 0) {
				for (std::size_t i = 0; i < this->element_v.size(); ++i)
					// find a element
					if (x.spe_component.at(this->element_v[i].ele_name) > 0) {
						element_available[i] = true;
					}

			}//if
		}//for

		std::vector<double> num_ways(this->species_network_v.size(), 0.0);
		// for each species, if conatains one available element, #ways make it +1
		for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
			for (std::size_t j = 0; j < this->element_v.size(); ++j) {
				if (element_available[j] == false)
					continue;
				else if (element_available[j] == true) {
					if (this->species_network_v[i].spe_component.at(this->element_v[j].ele_name) > 0)
						num_ways[i] += 1.0;
				}
			}//for
		}//for

		for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
			if (num_ways[i] > 0)
				num_ways[i] = 1.0 / num_ways[i];
		}

		this->rescale_concentration_data(num_ways);
	}

	void reactionNetworkODESolver::divide_prob_matrix_by_number_of_ways_making_species(std::vector<std::vector<double>>& prob_Mat)
	{
		this->rescale_prob_matrix_data(prob_Mat, this->num_ways_making_spe);
	}

	void reactionNetworkODESolver::update_number_of_ways_making_species()
	{
		this->num_ways_making_spe.assign(this->species_network_v.size(), 0.0);


		// find available elements, based on initial concentrations
		std::vector<bool> element_available(this->element_v.size(), false);
		for (auto x : this->species_network_v) {
			// initial concentration greater than zero
			if (x.spe_conc > 0) {
				for (std::size_t i = 0; i < this->element_v.size(); ++i)
					// find a element
					if (x.spe_component.at(this->element_v[i].ele_name) > 0) {
						element_available[i] = true;
					}

			}//if
		}//for

		// for each species, if conatains one available element, #ways make it +1
		for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
			for (std::size_t j = 0; j < this->element_v.size(); ++j) {
				if (element_available[j] == false)
					continue;
				else if (element_available[j] == true) {
					if (this->species_network_v[i].spe_component.at(this->element_v[j].ele_name) > 0)
						this->num_ways_making_spe[i] += 1.0;
				}
			}//for
		}//for

		for (std::size_t i = 0; i < this->species_network_v.size(); ++i) {
			if (this->num_ways_making_spe[i] > 0)
				this->num_ways_making_spe[i] = 1.0 / this->num_ways_making_spe[i];
		}
	}

	void reactionNetworkODESolver::set_single_source_concentration_constant()
	{
		if (this->single_source < 0)
			return;

		std::fill(this->concentration_data_pgt[this->single_source].begin(), this->concentration_data_pgt[this->single_source].end(),
			this->species_network_v[this->single_source].spe_conc);
	}

	void reactionNetworkODESolver::update_concentration_oriented_prob_matrix_from_single_source_path_based_prob_matrix(const std::vector<std::string>& path_v,
		const std::vector<double> &P2C, const std::vector<std::vector<double> >& path_prob_Mat, std::vector<std::vector<double> > &conc_prob_Mat)
	{

		if (this->single_source < 0)
			return;

		std::size_t Ntime = this->time_data_pgt.size();
		std::size_t Npath = path_v.size();

		for (std::size_t pi = 0; pi < Npath; ++pi) {
			std::string atom_followed; std::string real_path;
			this->split_atom_followed_and_pathway(path_v[pi], atom_followed, real_path);

			std::vector<size_t> spe_vec; std::vector<size_t> reaction_vec;
			this->parse_pathway_to_vector(real_path, spe_vec, reaction_vec);

			for (std::size_t tj = 1; tj < Ntime; ++tj) {
				conc_prob_Mat[spe_vec.back()][tj] +=
					path_prob_Mat[pi][tj]
					* P2C[spe_vec.back()]
					* this->species_network_v[spe_vec[0]].spe_component[atom_followed]
					/ this->species_network_v[spe_vec.back()].spe_component[atom_followed]
					/ this->num_ways_making_spe[spe_vec.back()];

			}//time tj
		}//path pi
	}

	rnk::when_where_t reactionNetworkODESolver::ODE_pathway_move_one_step(double time_in, rnk::vertex_t curr_vertex) {
		//Monte-Carlo simulation
		//generate the random number u_1 between 0 and 1.0
		double u_1 = rand->random01();

		rsp::my_time_t curr_time = reaction_time_from_importance_sampling_without_cutoff(time_in, curr_vertex, u_1);
		update_reaction_rate(curr_time, curr_vertex);
		rnk::when_where_t when_where(curr_time, curr_vertex);

		//if current species is not a dead species, not found
		if (std::find(this->dead_species.begin(), this->dead_species.end(), curr_vertex) == this->dead_species.end()) {//if1
			rsp::index_int_t next_reaction_index = random_pick_next_reaction(curr_vertex);
			//random pick next spe
			rnk::vertex_t next_vertex = random_pick_next_spe(next_reaction_index);

			when_where.first = curr_time;
			when_where.second = next_vertex;

		}//if1

		return when_where;
	}

	void reactionNetworkODESolver::ODE_pathway_sim_once(double init_time, double end_time, rnk::vertex_t init_vertex) {
		//set the pathway end time
		set_path_end_time(end_time);

		rnk::when_where_t when_where(init_time, init_vertex);
		std::size_t preceding_spe_index = 0;
		std::size_t preceding_time_index = 0;
		std::size_t current_spe_index = 0;
		std::size_t current_time_index = 0;

		while (when_where.first < this->path_end_time) {
			preceding_time_index = this->evaluate_index_at_time(when_where.first);
			preceding_spe_index = when_where.second;

			when_where = ODE_pathway_move_one_step(when_where.first, when_where.second);

			if (when_where.first < this->return_path_end_time()) {
				current_time_index = this->evaluate_index_at_time(when_where.first);
				current_spe_index = when_where.second;
				//preceding species in [preceding time, current_time], still preceding species
				this->update_spe_concentration_at_time_range(preceding_time_index, current_time_index,
					preceding_spe_index, this->number2Concentration);
				//preceding species at preceding time, the concentration should plus one (+1)
				this->update_spe_concentration_at_time(current_time_index, current_spe_index,
					this->number2Concentration);
			}
			else {
				// nothing happens before the end time point, include the very last point
				// because function "update_spe_conc_at_time_range" exclude the very last point, add it manually here
				this->update_spe_concentration_at_time_range(preceding_time_index,
					static_cast<std::size_t>(index_data_pgt.back() + 1),
					preceding_spe_index, this->number2Concentration);

			}

		}
	}

	void reactionNetworkODESolver::ODE_pathway_sim_N(double init_time, double end_time, std::size_t numberofTrajectory) {
		//probability
		std::vector<double> prob(this->species_network_v.size());
		for (std::size_t i = 0; i < prob.size(); ++i) {
			prob[i] = species_network_v[i].spe_conc;
		}

		std::size_t init_spe;

		for (std::size_t i = 0; i < this->iterationNumber; ++i) {
			//not first iteration
			if (i != 0) {
				this->update_spe_drc_based_on_spe_concentration_s_ct_np();
				this->integrate_propensity_function_pgt();
				this->init_spe_drc_int_pgt();

				//have to evaluate time every hopping
				this->init_spe_drc_int_time_pgt();

				this->update_reaction_rate_based_on_spe_concentration_s_ct_np();
				this->init_reaction_rate_pgt();
			}

			//reset conc_data_sr
			set_concentration_data_zero();
			set_concentration_at_time_zero_to_initial_fraction_or_concentration();

			for (size_t j = 0; j < numberofTrajectory; ++j) {
				//pick a initial spe randomly
				init_spe = rand->return_index_randomly_given_probability_vector(prob);
				ODE_pathway_sim_once(init_time, end_time, init_spe);

			}//for 2

			////re-scale conc_data_pgt, basically, convert from pathway probability to concentration
			//rescale_concentration_data(1.0);

		}//for 1
	}

	void reactionNetworkODESolver::get_probability_matrix(std::vector<std::vector<double>>& prob_Mat) const
	{
		std::copy(this->concentration_data_pgt.begin(), this->concentration_data_pgt.end(), prob_Mat.begin());
	}

	void reactionNetworkODESolver::ODEdirectlyEvaluatePathwayProbability_si_tj(const std::size_t si, const std::size_t tj, const std::vector<std::string>& pathway_vec, const double P2C, const std::size_t Nlocal, std::vector<std::vector<double>>& prob_Mat)
	{
		// species name
		std::string S = std::string("S") + boost::lexical_cast<std::string>(si);
		double conc_S = 0.0;//local concentration for S
		for (std::size_t k = 0; k < pathway_vec.size(); ++k) {
			// notice here that letters before first "S" represents the followed atom
			std::string atom_followed; std::string pathway;
			this->split_atom_followed_and_pathway(pathway_vec[k], atom_followed, pathway);

			std::vector<size_t> spe_vec; std::vector<size_t> reaction_vec;
			this->parse_pathway_to_vector(pathway, spe_vec, reaction_vec);

			// if the concentration of initial species is zero, no need to calculate
			if (this->species_network_v[spe_vec[0]].spe_conc > 0.0) {
				for (std::size_t l = 0; l < Nlocal; ++l) {
					// the initial concentration of starting species* path_prob* p2c of ending species,
					// notice p2c is calculated from species components
					// if the concentration of initial species is zero, no need to calculate
					conc_S += this->species_network_v[spe_vec[0]].spe_conc*
						pathway_prob_input_pathway_sim_once(0.0, this->time_data_pgt[tj], spe_vec, reaction_vec, atom_followed)
						* P2C
						* this->species_network_v[spe_vec[0]].spe_component[atom_followed]
						/ this->species_network_v[si].spe_component[atom_followed];

				}//l-->trajectory

			}//if
		}//k-->pathway
		 //same time points, same end speices
		prob_Mat[si][tj] = conc_S;
	}

	/*
	 * pathway based ODE solver, call it a ODE_solver, evaluate the pathway probability directly
	 * input-->pathway vector, pathways that are used to evaluate pathway probability so that to evaluate concentrations
	 *		-->Ntrajectory, number of trajectories to run to evaluate the pathway probability
	 *		-->list, factors that convert pathway probability to concentration,namely P2C
	 * just for a single core, got to renormalize the concentration at the end, basically divided by the total number of trajectories
	 */
	void reactionNetworkODESolver::ODEdirectlyEvaluatePathwayProbability(const std::vector<std::string> &pathway_vec, const std::vector<double> &P2C, const std::size_t Nlocal, const std::size_t topN, std::vector<std::vector<double> > &prob_Mat) {
		//set all the concentration to zero
		set_concentration_data_zero();

		//std::cout << "Number of trajectory on this core:\t" << Nlocal << std::endl;

		const std::size_t Nspecies = this->concentration_data_pgt.size();
		const std::size_t Ntimepoints = this->time_data_pgt.size();
		assert(P2C.size() == Nspecies); //make sure that P2C has the same number of element as conc_data_sr

		for (std::size_t si = 0; si < Nspecies; ++si) {
			//for every species, got to evaluate the concentration at every time point
			std::string S = std::string("S") + boost::lexical_cast<std::string>(si);
			std::vector<std::string> pathway_vec_t;
			//topN ends with S
			pathwayHandler::pathway_ends_with(S, pathway_vec, pathway_vec_t, topN); //topN pathways

			//not the first time point
			for (std::size_t tj = 1; tj < Ntimepoints; ++tj) {
				ODEdirectlyEvaluatePathwayProbability_si_tj(si, tj, pathway_vec_t, P2C[si], Nlocal, prob_Mat);
			}//j-->Ntimepoints
		}//i-->Nspecies

	}

	double reactionNetworkODESolver::single_source_ODEdirectlyEvaluatePathwayProbability_pi_tj_integrate_over_t0_MC(std::string path, std::size_t tj, const std::size_t Nlocal)
	{
		//local probability for S
		double prob_S = 0.0;
		// notice here that letters before first "S" represents the followed atom
		std::string atom_followed; std::string real_path;
		this->split_atom_followed_and_pathway(path, atom_followed, real_path);

		std::vector<size_t> spe_vec; std::vector<size_t> reaction_vec;
		this->parse_pathway_to_vector(real_path, spe_vec, reaction_vec);

		for (std::size_t l = 0; l < Nlocal; ++l) {

			// integrate over points before, Monte Carlo integral, uniformaly sample initial time
			for (std::size_t k = 0; k < tj; ++k) {
				double t0 = rand->random_min_max(0, time_data_pgt[tj]);

				prob_S += pathway_prob_input_pathway_sim_once(t0, time_data_pgt[tj], spe_vec, reaction_vec, atom_followed)
					* this->evaluate_single_source_additional_concentration_at_time(t0);
			}//k-->time points before

		}//l-->trajectory

		return prob_S;
	}

	double reactionNetworkODESolver::single_source_ODEdirectlyEvaluatePathwayProbability_pi_tj_integrate_over_t0_rectangle(std::string path, std::size_t tj, const std::size_t Nlocal)
	{
		//local probability for S
		double prob_S = 0.0;
		// notice here that letters before first "S" represents the followed atom
		std::string atom_followed; std::string real_path;
		this->split_atom_followed_and_pathway(path, atom_followed, real_path);

		std::vector<size_t> spe_vec; std::vector<size_t> reaction_vec;
		this->parse_pathway_to_vector(real_path, spe_vec, reaction_vec);

		for (std::size_t l = 0; l < Nlocal; ++l) {

			// integrate over points before, Rectangle rule
			for (std::size_t k = 0; k < tj; ++k) {
				prob_S += pathway_prob_input_pathway_sim_once(time_data_pgt[k], time_data_pgt[tj], spe_vec, reaction_vec, atom_followed)
					* this->evaluate_single_source_additional_concentration_at_time(time_data_pgt[k]);
			}//k-->time points before

		}//l-->trajectory

		return prob_S;
	}

	void reactionNetworkODESolver::single_source_ODEdirectlyEvaluatePathwayProbability_pi(std::string path, const std::size_t Nlocal, std::vector<double>& prob_v)
	{
		std::size_t Npoints = this->time_data_pgt.size();
		prob_v.assign(Npoints, 0.0);

		for (std::size_t tj = 1; tj < Npoints; ++tj)
			prob_v[tj] = this->single_source_ODEdirectlyEvaluatePathwayProbability_pi_tj_integrate_over_t0_MC(path, tj, Nlocal);
	}

	void reactionNetworkODESolver::single_source_ODEdirectlyEvaluatePathwayProbability_pathVector(std::vector<std::string> path_v, const std::size_t Nlocal, std::vector<std::vector<double>>& path_prob_Mat)
	{

		if (this->single_source < 0)
			return;

		std::size_t Npath = path_v.size();

		path_prob_Mat.resize(Npath);

		for (std::size_t pi = 0; pi < Npath; ++pi) {
			this->single_source_ODEdirectlyEvaluatePathwayProbability_pi(path_v[pi], Nlocal, path_prob_Mat[pi]);
		}
	}



}/*namespace reactionNetworkODESolver_sr*/

#endif
