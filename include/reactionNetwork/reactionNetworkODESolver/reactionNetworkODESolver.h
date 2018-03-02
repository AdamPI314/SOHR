#ifndef __REACTIONNETWORKODESOLVER_H_
#define __REACTIONNETWORKODESOLVER_H_

#include "../superReactionNetwork/superReactionNetwork.h"
#include "../../propagator/SOHRPropagator/SOHRPropagator.h"

namespace reactionNetworkODESolver_sr {
	namespace rnk = reactionNetwork_sr;
	namespace rsp = relationshipParser_sr;
	namespace pgt = propagator_sr;

	class reactionNetworkODESolver :public rnk::superReactionNetwork, public pgt::SOHRPropagator {
		//private:
		//number of ways making species
		std::vector<double> num_ways_making_spe;

		//Whether there is a single source
		int single_source = -1;

		//single source destruction rate, rnos-->reaction network ODE solver
		std::vector<rsp::reaction_rate_t> single_source_spe_dr_rnos;

		//single source additional concentration over time, rnos-->reaction network ODE solver
		std::vector<rsp::concentration_t> single_source_additional_concentration_data_rnos;
		/*
		* 	cubic spline
		* 	species destruction rate constant
		* 	got to delete pointer in destructor!!!
		*/
		Linear_interp* time_single_source_additional_concentration_pointer_rnos = nullptr;

	public:
		reactionNetworkODESolver(std::vector<double> uncertainties, std::size_t random_seed_for_this_core, std::string cwd_in);
		~reactionNetworkODESolver();

	public:
		//return number of species and number of time points
		std::pair<std::size_t, std::size_t> get_size_of_concentration_data() const;
		//return initial concentration
		std::vector<double> get_initial_concentration() const;

	public:
		//set Prob_max(tau^{j}|t+tau^{j-1};S^{j-1}), with pathway_end_time fixed
		bool set_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t) override;
		//set the  prob_min of a species at a given time
		bool set_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t) override;
		double get_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t) const override;
		double get_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t) const override;

	public:
		//update out reaction rate of current vertex at a specific time point
		bool update_reaction_rate(double in_time, rnk::vertex_t curr_vertex) override;
		//update rate at specific time point
		bool update_reaction_rate(double in_time) override;

	public:
		//got to instantiate this function otherwise this sub-class is still pure-virtual class can not be instantiated
		double evaluate_spe_concentration_at_time(double time, std::size_t index = 0) const override { return 0.0; };
		double evaluate_spe_drc_at_time(double time, std::size_t index = 0) const override { return 0.0; };		
		double evaluate_chattering_group_ss_prob_at_time(double in_time, size_t index = 0) const override { return 0.0; };

	public:
		double prob_chattering_group_will_react_in_a_time_range(double init_time, double pathway_end_time, size_t curr_chattering_group) override { return 0.0;  };

	public:
		/*
		* reaction time from importance sampling, exact time
		* if reaction_time> tau, let it be, don't cut it off
		*/
		double reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time, rnk::vertex_t curr_spe, double Y) override;
		/*
		* reaction time from importance sampling, exact time
		* if reaction_time> tau, cut it off
		*/
		double reaction_time_from_importance_sampling(rsp::my_time_t curr_time, rnk::vertex_t curr_spe, double Y) override;

		/*
		* chattering group reaction time from importance sampling, exact time
		* if reaction_time> tau, let it be, don't cut it off
		*/
		double chattering_group_reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time, rnk::vertex_t curr_group, double Y) override;
		/*
		* inside chattering group, randomly pick next species, for example, A<=>B chattering set,
		* either pick A or B based on their SSA(or equilibrium) concentration
		*/
		std::vector<double> chattering_group_probability_vector(rsp::index_int_t chattering_group_id, double time) override;
		rnk::vertex_t inside_chattering_group_random_pick_next_spe(rsp::index_int_t chattering_group_id, double time) override;

	public:
		/*
		* set the conc_data_sr[i][0] to the initial probability
		*/
		void set_concentration_at_time_zero_to_initial_fraction_or_concentration();
		void set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(std::vector<std::vector<double> > & prob_Mat);

		void set_probability_matrix_of_a_species_constant(std::vector<std::vector<double> > & prob_Mat, std::size_t ind);
		void set_probability_matrix_of_single_source_constant(std::vector<std::vector<double> > & prob_Mat);

		//update single_source_additional_concentration_data_rnos, from previous iteration
		void update_single_source_additional_concentration_data_rnos_s_ct_np();
		//initiate single source additional concentration pointer
		void init_time_single_source_additional_concentration_pointer_rnos();
		double evaluate_single_source_additional_concentration_at_time(double in_time) const;

		/*update source species dr based on spe concentration s ct np*/
		void update_single_source_species_dr_based_on_spe_concentration_s_ct_np();

		/*
		* check zero concentration, if it is zero, set it to be "SOHR_init.deltaConcentration"
		*/
		void check_zero_concentration(double deltaConcentration);

		//divide the concentration of a species, like OH, by #ways that can make it, here for OH, we can either follow
		//O or H to make OH, the two ways should give the same result so divide the concentration by a factor of 2
		void divide_concentration_by_number_of_ways_making_it();
		//divide the concentration of a species, like OH, by #ways that can make it, here for OH, we can either follow
		//O or H to make OH, the two ways should give the same result so divide the concentration by a factor of 2
		void divide_prob_matrix_by_number_of_ways_making_species(std::vector<std::vector<double> > & prob_Mat);
		void update_number_of_ways_making_species();

		//set single source concentration constant
		void set_single_source_concentration_constant();
		//update concentration data from pathway based probability matrix
		void update_concentration_oriented_prob_matrix_from_single_source_path_based_prob_matrix(const std::vector<std::string> &path_v,
			const std::vector<double> &P2C, const std::vector<std::vector<double> > &path_prob_Mat, std::vector<std::vector<double> > &conc_prob_Mat);

	public:
		/*
		* return when it is, where we are, for MPI
		*/
		rnk::when_where_t ODE_pathway_move_one_step(double time, rnk::vertex_t curr_vertex);
		/*
		* pathway based ODE solver, call it a propagator, every trajectory corresponds to a pathway
		* along every trajectory, update the concentration of related species
		*/
		void ODE_pathway_sim_once(double init_time, double end_time, rnk::vertex_t init_vertex);
		/*
		* run trajectory multiple times, which include inner loop and outer loop
		* outer loop update destruction rate constant (drc)
		* inner loop run lots of trajectories with the same drc
		* randomly select the initial spe based the inital concentration
		*/
		void ODE_pathway_sim_N(double init_time, double end_time, std::size_t numberofTrajectory);
		/*
		* pathway based ODE solver, call it a ODE_solver, evaluate the pathway probability directly
		* input-->pathway vector, pathways that are used to evaluate pathway probability so that to evaluate concentrations
		* 		-->Ntrajectory, number of trajectories to run to evaluate the pathway probability
		* 		-->list, factors that convert pathway probability to concentration,namely P2C
		* just for a single core, got to renormalize the concentration at the end, basically divided by the total number of trajectories
		*/
		/*
		* return probability matrix, which is concentration vs. time actually
		*/
		void get_probability_matrix(std::vector<std::vector<double> > &prob_Mat) const;
		//for species si, time point tj 
		void ODEdirectlyEvaluatePathwayProbability_si_tj(const std::size_t si, const std::size_t tj, const std::vector<std::string> &pathway_vec,
			const double P2C, const std::size_t Nlocal, std::vector<std::vector<double> > &prob_Mat);
		void ODEdirectlyEvaluatePathwayProbability(const std::vector<std::string> &pathway_vec,
			const std::vector<double> &P2C, const std::size_t Nlocal, const std::size_t topN, std::vector<std::vector<double> > &prob_Mat);

		//path starting with source species, input path pi evaluate pathway probability at time point tj, integrate over initial time t0, Monte Carlo integration
		double single_source_ODEdirectlyEvaluatePathwayProbability_pi_tj_integrate_over_t0_MC(std::string path, std::size_t tj, const std::size_t Nlocal);
		//Rectangle integration rule
		double single_source_ODEdirectlyEvaluatePathwayProbability_pi_tj_integrate_over_t0_rectangle(std::string path, std::size_t tj, const std::size_t Nlocal);

		//path starting with source species, input path pi, evaluate pathway probability each time point
		void single_source_ODEdirectlyEvaluatePathwayProbability_pi(std::string path, const std::size_t Nlocal, std::vector<double> &prob_v);
		//path starting with source species, input path vector, evaluate pathway probability at each time point
		void single_source_ODEdirectlyEvaluatePathwayProbability_pathVector(std::vector<std::string> path_v, const std::size_t Nlocal, std::vector<std::vector<double> > &path_prob_Mat);

	};//class reactionNetworkODESolver




}/*namespace reactionNetworkODESolver_sr*/

#endif
