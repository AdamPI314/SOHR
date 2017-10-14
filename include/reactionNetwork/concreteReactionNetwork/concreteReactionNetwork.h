#ifndef __CONCRETEREACTIONNETWORK_H_
#define __CONCRETEREACTIONNETWORK_H_

#include "../superReactionNetwork/superReactionNetwork.h"


//Try including nr3.h after all boost and system headers!!!
//#include "../cubicSpline/interp_1d.h"

namespace reactionNetwork_sr {
	namespace mt = misc_template;
	namespace rsp = relationshipParser_sr;
	namespace pgt = propagator_sr;

	class concreteReactionNetwork :public superReactionNetwork {

	protected:
		//propagator
		//polymorphism
		pgt::superPropagator* propagator;

	public:
		concreteReactionNetwork(std::vector<double> uncertainties, std::size_t random_seed_for_this_core, std::string cwd_in);
		~concreteReactionNetwork();
	public:
		//merge fast transitions, for example A=B first order reaction, assume the forward and backward reactions are both
		//very fast, if we keep this fast transition in the network, this will result in long path in a short time range
		//Here is how we do it, the idea is to treat these two species as a single one
		//One philosophy is to minimize the change
		//1) set fast reaction rates to be zeros
		//2) chattering with a time scale
		//re-call initiate_cubic_spline()
		void merge_chatterings();
		void update_species_chattering_group_id();

	public:
		/*
		* set the initial concentration
		*/
		bool set_init_spe_concentration(rsp::my_time_t in_time);

	public:
		//update out reaction rate of current vertex at a specific time point
		bool update_reaction_rate(double in_time, vertex_t curr_vertex) override;
		//update rate at specific time point
		bool update_reaction_rate(double in_time) override;

	public:
		/*
		* reaction time from importance sampling, exact time
		* if reaction_time> tau, let it be, don't cut it off
		*/
		double reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time, vertex_t curr_spe, double Y) override;
		/*
		* reaction time from importance sampling, exact time
		* if reaction_time> tau, cut it off
		*/
		double reaction_time_from_importance_sampling(rsp::my_time_t curr_time, vertex_t curr_spe, double Y) override;
		/*
		* chattering group reaction time from importance sampling, exact time
		* if reaction_time> tau, let it be, don't cut it off
		*/
		double chattering_group_reaction_time_from_importance_sampling_without_cutoff(rsp::my_time_t curr_time, vertex_t curr_group, double Y) override;

	public:
		//set Prob_max(tau^{j}|t+tau^{j-1};S^{j-1}), with pathway_end_time fixed
		bool set_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t) override;
		//set the  prob_min of a species at a given time
		bool set_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t) override;
		double get_spe_prob_max_at_a_time(double init_time, double end_time, size_t index_t) const override;
		double get_spe_prob_min_at_a_time(double init_time, double end_time, size_t index_t) const override;

	public:
		double evaluate_spe_concentration_at_time(double time, std::size_t index = 0) const override;
		double evaluate_spe_drc_at_time(double time, std::size_t index = 0) const override;
		double evaluate_chattering_group_ss_prob_at_time(double in_time, size_t index = 0) const override;

	public:
		double prob_chattering_group_will_react_in_a_time_range(double init_time, double pathway_end_time, size_t curr_chattering_group) override;

	public:
		//return temperature target time
		rsp::my_time_t return_temperature_target_time() const {
			return propagator->evaluate_time_at_temperature(propagator->return_target_temperature());
		}

	public:
		/*
		* return where we are, when it is, for MPI ,for arrival time
		* force the next reaction to be  next_reaction, next species to be next_spe, because what we want is just reaction time
		* directly set the next reaction and next species the ones we want
		*/
		double pathway_AT_sim_move_one_step(double when_time, size_t curr_spe, size_t next_reaction, size_t next_spe);
		//input a pathway, return its arrival time
		double pathway_AT_input_pathway_sim_once(double init_time, double pathway_end_time, std::string pathway_in);

	public:
		/*
		* recursive relation
		* {I^{(j)}}({\tau ^{(j)}},{S^{(j - 1)}},{R^{(j)}},{S^{(j)}},{N^{(j)}})
		* = \frac{1}{{{N^{(j)}}}}\sum\limits_{i = 1}^{{N^{(j)}}} {{f^{(j)}}(\tau _i^{(j)},{S^{(j - 1)}},{R^{(j)}},{S^{(j)}}) \times }
		* {I^{(j + 1)}}({\tau ^{(j + 1)}},{S^{(j)}},{R^{(j + 1)}},{S^{(j + 1)}},{N^{(j + 1)}})
		* neglect the the first i=0 term here, we assume reaction must occur
		* tau_j_minus_1- current time, lower limit of the integral, or says the initial time of current species
		* spe_vec- species vector
		* reaction_vec- reaction vector
		* N_subvolume- how may sub-interval every interval will be divided into
		* j_th- the j_th integral
		*/
		double pathway_prob_input_pathway_recursive_relation(double tau_j_minus_1, std::vector<std::size_t> spe_vec, std::vector<std::size_t> reaction_vec,
			std::vector<std::size_t> N_subvolume, std::size_t j_th);

	public:
		/*
		* Flux analysis
		* Flux based pathway analysis
		*/
		/*
		* out flux of a species at a specific time
		* return a vector<species index, flux>
		*/
		std::vector<spe_index_flux_t> spe_out_flux_at_a_time(rsp::my_time_t in_time, std::size_t curr_spe, std::string atom_followed = "H");
		/*
		* calculate the flux of all species
		* normalize them
		* print to file
		*/
		void flux_analysis_w2f(rsp::my_time_t in_time);

	public:
		void spe_concentration_w2f_rnk(double in_time, std::string str) const;

	};//class reactionNetwork




}/*namespace reactionNetwork_sr*/

#endif
