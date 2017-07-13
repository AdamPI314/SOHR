#ifndef __SSAPROPAGATOR_H_
#define __SSAPROPAGATOR_H_

#include "../../tools/misc/misc_template.h"
#include "../../relationshipParser/relationshipParser.h"

#include "../superPropagator/superPropagator.h"
#include "../../srkin/srkin.h"
#include "../../random/random.h"

#include "../../cubicSpline/interp_1d.h"

namespace propagator_sr {
	class ssaPropagator : public superPropagator, public srkin_sr::srkin {
	protected:
		boost::uint32_t random_seed_for_this_core;
		//random number generator
		random_sr::random* rand;

	public:
		ssaPropagator(std::vector<double> uncertainties, std::size_t random_seed_for_this_core, std::string cwd_in);
		~ssaPropagator();

	public:
		//template method
		virtual void propagate_pgt() override;
		/*
		* constant volume-->cv
		*/
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*/
		void time_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;
		/* constant volume-->cv
		* constant temperature-->ct
		*/
		void time_propagator_cv_ct_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_temperature as stopping criteria
		* s2m means storing to memory
		*/
		void temperature_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t, double end_temperature) override;
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, stop when the temperature doesn't change,
		* or says the system reaches equilibrium
		* s2m means storing to memory
		*/
		void equilibrium_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t) override;
		/*
		* constant ckstore.pressure-->cp
		*/
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*/
		void time_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_temperature as stopping criteria
		* s2m means storing to memory
		*/
		void temperature_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t, double end_temperature) override;
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, stop when the temperature doesn't change,
		* or says the system reaches equilibrium
		* s2m means storing to memory
		*/
		void equilibrium_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t) override;
		/*
		* constant ckstore.pressure-->cp
		* constant temperature-->ct
		*/
		/*
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*/
		void time_propagator_cp_ct_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;
		/*
		* surface reaction -->s
		* constant temperature-->ct
		* not dependent of ckstore.pressure, independent of ckstore.pressure
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*
		* This solver can be used to solve general kinetic equations, which maybe have no relation with temperature, ckstore.pressure or surface
		* Like Lotka-Volterra model, which turns out to be a set of ODE. This model comprises of a simple mechanism
		* http://www.idea.wsu.edu/OscilChem/
		*/
		void time_propagator_s_ct_np_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;
		void time_propagator_s_ct_np_s2m_find_one_transition_pgt(std::vector<double> uncertainties, double critical_time, double end_time);

		// cc-->constant concentration, hold the concentration of the first species to be constant
		void time_propagator_s_ct_np_cc1_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;
		// cc-->constant concentration, hold the concentration of the first two species to be constant
		void time_propagator_s_ct_np_cc2_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;

		//Gillespie Stochastic Simulation Algorithm update step, here decouple c_t and reaction_rates so that this method is applicable to
		//different systems. move one reaction ahead, namelly NEXT REACTION METHOD
		void ssa_update_next_reaction_pgt(double* c_t, const std::vector<double>& reaction_rates, double& tau);
		//Move time dt, For this specific system-->const temperature, no pressure, caculate reaction rates from c_t
		void ssa_update_dt_s_ct_np_pgt(const double Temp, const double dt, double* c_t, std::vector<double>& reaction_rate_v_tmp, double& tau);
	};

}//namespace propagator_sr

#endif
