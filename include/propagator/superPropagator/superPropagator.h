#ifndef __SUPERPROPAGATOR_H_
#define __SUPERPROPAGATOR_H_

#include <memory>
#include "../../tools/misc/misc_template.h"
#include "../../odeSolver/odesolver.h"
#include "../../relationshipParser/relationshipParser.h"
#include <boost/property_tree/ptree.hpp> //for property_tree
#include "../../tools/species_group/species_group.h"

#include "../../cubicSpline/interp_1d.h"
#include "../../cubicSpline/interp_linear.h"

namespace propagator_sr {
	namespace mt = misc_template;
	namespace rsp = relationshipParser_sr;

	//species destruction rate constant(drc) type
	typedef double spe_drc_t;
	//species production rate type
	typedef double spe_production_rate_t;

	class superPropagator {
	public:
		//current working directory
		std::string cwd_pgt;
	public:
		//configurations, read configuration file named "setting.json"
		//super propagator pt
		boost::property_tree::ptree pgt_pt;
	public:
		std::vector<double> uncertainties;
	public:
		rsp::reactionNetwork_chemkin_index_map_t reactionNetwork_chemkin_index_map;
	public:
		//pgt represents propagator
		std::vector<rsp::my_time_t> time_data_pgt;
		std::vector<rsp::temperature_t> temperature_data_pgt;
		std::vector<rsp::pressure_t> pressure_data_pgt;

		////total concentration, if there is not external source, the total concentration is constant over time
		////if there is a source, the total concentration is a integral over rate*time
		////didn't define pointer up to now, not necessary, since this is only used in SOHR ODE at this moment
		//std::vector<rsp::concentration_t> total_concentration_data_pgt;

		std::vector<std::vector<rsp::concentration_t> > concentration_data_pgt;
		std::vector<std::vector<rsp::reaction_rate_t> > reaction_rate_data_pgt;
		std::vector<std::vector<spe_drc_t> > spe_drc_data_pgt;
		//integral with respect to the initial drc, which means this is the integral minus the initial value
		std::vector<std::vector<spe_drc_t> > spe_drc_int_data_pgt;

		//species production rate, all the elementary source terms of species production rate
		//only created in cv_s2m_dl, constant volumn propagator
		//std::vector<std::vector<spe_production_rate_t> > spe_production_rate_data_pgt;

	public:
		//shared pointer
		std::shared_ptr<species_group_sr::chattering> sp_chattering_pgt = std::make_shared<species_group_sr::chattering>();

		//chattering group data
		//(1) chattering group drc, the time scale of this chattering mode, there might be multiple modes, only consider
		//the first mode
		//(2) project the mode on species direction, probability being that species, or says steady state probability
		std::vector<std::vector<double> > chattering_group_k_data_pgt;
		std::vector<std::vector<double> > chattering_group_k_int_data_pgt;
		std::vector<std::vector<double> > chattering_group_ss_prob_data_pgt;

	public:
		/*
		* 	cubic spline
		* 	species destruction rate constant
		* 	got to delete pointer in destructor!!!
		*/
		std::vector<Linear_interp*> spe_drc_pgt;
		//integral
		std::vector<Linear_interp*> spe_drc_int_pgt;
		//integral reverse
		std::vector<Linear_interp*> spe_drc_int_time_pgt;
		//reaction rate
		std::vector<Linear_interp*> reaction_rate_pgt;
		//concentration
		std::vector<Linear_interp*>	concentration_pgt;
		////total concentration
		//Linear_interp* total_concentration_pgt;

		//temperature versus time
		Linear_interp* time_temperature_pgt;
		//time vs. temperature
		Linear_interp* temperature_time_pgt;

		//Pressure versus time
		Linear_interp* time_pressure_pgt;

	public:
		//chattering group rate constant k, time scale pointer
		std::vector<Linear_interp*> chattering_group_k_pgt;
		//integral
		std::vector<Linear_interp*> chattering_group_k_int_pgt;
		//integral reverse
		std::vector<Linear_interp*> chattering_group_k_int_time_pgt;
		//steady state probability
		std::vector<Linear_interp*> chattering_group_ss_prob_pgt;

	public:
		superPropagator(std::vector<double> uncertainties_in, std::string cwd_in);
		/*
		//http://stackoverflow.com/questions/461203/when-to-use-virtual-destructors
		//Virtual destructors are useful when you can delete an instance of a derived class
		//through a pointer to base class:

		//class Base
		//{
		////some virtual methods
		//};

		//class Derived : public Base
		//{
		//~Derived()
		//{
		////Do some important cleanup
		//}
		//}
		//Here, you'll notice that I didn't declare Base's destructor to be virtual.
		//Now, let's have a look at the following snippet:

		//Base *b = new Derived();
		////use b
		//delete b; //Here's the problem!
		//Since Base's destructor is not virtual and b is a Base* pointing to a Derived object,
		//delete b has undefined behaviour. In most implementations, the call to the destructor
		//will be resolved like any non-virtual code, meaning that the destructor of the base
		//class will be called but not the one of the derived class, resulting in resources leak.
		*/
		virtual ~superPropagator();
	public:
		//initialize the pointer, kind of template method
		void initialize_cubic_spline_pointer();
	public:
		void update_temporary_data_pgt(const int nkk, const int neq,
			const double ti,
			const double * const c_t,
			const double * const CDOT_t,
			const double * const DDOT_t,
			const double * const FWDR_t,
			const double * const REVR_t,
			const double * const xgst);

	public:
		void find_chattering_group(const std::vector<rsp::spe_info_base> &species_network_v);
		//update chattering group-pairs-reactions
		void update_chattering_group_pairs_reactions(const std::vector<rsp::spe_info_base> &species_network_v, const std::vector<rsp::reaction_info_base> &reaction_network_v, std::string atom_followed = "H");

		//return shared pointer of chattering
		std::shared_ptr<species_group_sr::chattering> get_sp_of_chattering();

		//chattering species and reaction, local reaction with fast inter-conversion rate
		void set_chattering_spe_from_file_pgt();
		void set_chattering_reactions_from_file_pgt();

		//cancel fast transition within each chattering group
		void subtract_chattering_reaction_contribution_from_species_drc_pgt();
		//set fast transition A=B's pseudo-first order rate constant given list of trapped species pair
		void update_drc_and_equilibrium_probability_of_chattering_group();
		//set the reaction rate of fast reactions to be zero
		void set_chattering_reaction_rates_to_zero_pgt();

	public:
		//return target temperature
		rsp::temperature_t return_target_temperature() const;

	public:
		//print temperature-target concentration into file name "temperature_targeted_conc.csv"
		void print_temperature_targeted_concentration_pgt(std::string outfile = "/output/temperature_targeted_conc.csv") const;
		//Print some interesting stuff to study...
		void print_pgt();
		//Write some interesting stuff to file
		void w2f_pgt(std::string tag = "");
		//write species concentration to file
		void spe_concentration_w2f_pgt(double in_time, std::string tag = "") const;
		//print conc at specific time into file name "conc.csv", mostly for speciation
		void spe_concentration_w2f_pgt(std::vector<double> time_in, std::string tag = "") const;

		//just for test at this moment
		void print_final_concentration_of_a_spe_pgt(std::size_t index = 0) {
			std::cout << std::setprecision(16) << concentration_data_pgt[index].back() << std::endl;
		}

		//return temperature target time
		rsp::my_time_t return_temperature_target_time() const {
			return evaluate_time_at_temperature(return_target_temperature());
		}


	public:
		//Read in the initial Temperature, Pressure, molar fractions. from file "setting.json"
		void read_configuration(double & Temp, double & Pressure, std::size_t nkk, double x[]);
		//initialize lsode, read lsode settings
		void initialize_lsode(rsp::my_time_t& dt, std::string outfile = "/output/general_output.out");

	public:
		//template function, use template design pattern here
		void initiate_propagator();
		//template function, use template design pattern here
		void propagate_and_initiate_cublic_spline_pgt() {
			//propagator
			propagate_pgt();
			//initiate cubic spline
			initiate_cubic_spline();

		}
		//template function, use template design pattern here 
		//which is a virtual/abstract method, and a pure virtual method
		virtual void propagate_pgt() {};
		//initiate cubic spline and integral of drc, etc.
		void initiate_cubic_spline() {
			integrate_propensity_function_pgt();

			init_spe_drc_pgt();
			init_spe_drc_int_pgt();
			init_spe_drc_int_time_pgt();
			init_reaction_rate_pgt();
			init_concentration_pgt();
			//init_time_total_concentration_pgt();
			init_time_temperature_pgt();
			init_temperature_time_pgt();
			init_time_pressure_pgt();

			//chatterings
			if (this->chattering_group_k_data_pgt.size() > 0 && this->chattering_group_k_data_pgt[0].size() > 0) {
				init_time_chattering_group_k_pgt();
			}
			this->integrate_chattering_group_propensity_function_pgt();
			if (this->chattering_group_k_int_data_pgt.size() > 0 && this->chattering_group_k_int_data_pgt[0].size() > 0) {
				init_chattering_group_k_int_pgt();
			}
			if (this->chattering_group_k_int_data_pgt.size() > 0 && this->chattering_group_k_int_data_pgt[0].size() > 0) {
				init_chattering_group_k_int_time_pgt();
			}

			if (this->chattering_group_ss_prob_data_pgt.size() > 0 && this->chattering_group_ss_prob_data_pgt[0].size() > 0) {
				init_time_chattering_group_ss_prob_pgt();
			}
		}//initiate_cubic_spline


#if defined(__CHEMKIN_AVAILABLE_) && defined(__LSODE_AVAILABLE_)
		
		void convert_molar_concentration_to_mole_fraction();
		void convert_mole_fraction_to_molar_concentration();

#endif // defined(__CHEMKIN_AVAILABLE_) && defined(__LSODE_AVAILABLE_)

		virtual void propagator_from_file(std::string tag = "") {};

#if defined(__CHEMKIN_AVAILABLE_) && defined(__LSODE_AVAILABLE_)

		/*
		* constant volume-->cv
		*/
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*/
		virtual void time_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) {};
		/* constant volume-->cv
		* constant temperature
		*/
		virtual void time_propagator_cv_ct_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) {};
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_temperature as stopping criteria
		* s2m means storing to memory
		*/
		virtual void temperature_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t, double end_temperature) {};
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, stop when the temperature doesn't change,
		* or says the system reaches equilibrium
		* s2m means storing to memory
		*/
		virtual void equilibrium_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t) {};

		/*
		* constant pressure-->cp
		*/
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*/
		virtual void time_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) {};
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_temperature as stopping criteria
		* s2m means storing to memory
		*/
		virtual void temperature_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t, double end_temperature) {};
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, stop when the temperature doesn't change,
		* or says the system reaches equilibrium
		* s2m means storing to memory
		*/
		virtual void equilibrium_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t) {};
		/*
		* constant pressure-->cp
		* constant temperature-->ct
		*/
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*/
		virtual void time_propagator_cp_ct_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) {};

		/*
		* surface reaction -->s
		* constant temperature-->ct
		* not dependent of pressure, independent of pressure
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*
		* This solver can be used to solve general kinetic equations, which maybe have no relation with temperature, pressure or surface
		* Like Lotka-Volterra model, which turns out to be a set of ODE. This model comprises a simple mechanism
		* http://www.idea.wsu.edu/OscilChem/
		*
		*/
		virtual void time_propagator_s_ct_np_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) {};
		//cc-->constant concentration, hold the concentration of the first species to be constant
		virtual void time_propagator_s_ct_np_cc1_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) {};
		//cc-->constant concentration, hold the concentration of the first two species to be constant
		virtual void time_propagator_s_ct_np_cc2_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) {};

#endif // defined(__CHEMKIN_AVAILABLE_) && defined(__LSODE_AVAILABLE_)

	public:
		//integral of the propensity function or says the destruction rate constant function
		void integrate_propensity_function_pgt();
		//initialize species destruction rate constant
		bool init_spe_drc_pgt();
		//initialize cumulative species destruction rate constant
		bool init_spe_drc_int_pgt();
		//initialize time vs. cumulative species destruction rate constant
		bool init_spe_drc_int_time_pgt();
		//initialize reaction rates
		bool init_reaction_rate_pgt();
		//initialize concentration
		bool init_concentration_pgt();
		////initialize total concentration vs. time
		//bool init_time_total_concentration_pgt();
		//initialize temperature vs. time
		bool init_time_temperature_pgt();
		//initialize time vs. temperature
		bool init_temperature_time_pgt();
		//initialize pressure vs. time
		bool init_time_pressure_pgt();

		//chattering group k
		bool init_time_chattering_group_k_pgt();
		//integral of the propensity function or says the destruction rate constant function
		void integrate_chattering_group_propensity_function_pgt();
		bool init_time_chattering_group_ss_prob_pgt();
		//initialize cumulative species destruction rate constant
		bool init_chattering_group_k_int_pgt();
		//initialize time vs. cumulative species destruction rate constant
		bool init_chattering_group_k_int_time_pgt();

	public:
		double evaluate_temperature_at_time(double in_time) const;
		double evaluate_time_at_temperature(double in_temperature) const;
		double evaluate_pressure_at_time(double in_time) const;

		double evaluate_concentration_at_time(double in_time, size_t index = 0) const;
		double evaluate_spe_drc_at_time(double in_time, size_t index = 0) const;
		double evaluate_spe_drc_int_at_time(double in_time, size_t index = 0) const;
		double evaluate_time_at_spe_drc_int(double integral, size_t index = 0) const;
		double evaluate_reaction_rate_at_time(double in_time, size_t index = 0) const;

	public:
		double evaluate_chattering_group_ss_prob_at_time(double in_time, size_t index = 0) const;
		//evaluate chattering group's k or says time scale
		double evaluate_chattering_group_k_at_time(double in_time, size_t chattering_group_id = 0) const;
		double evaluate_chattering_group_k_int_at_time(double in_time, size_t index = 0) const;
		double evaluate_time_at_chattering_group_k_int(double integral, size_t index = 0) const;

	};


}//namespace propagator_sr

#endif
