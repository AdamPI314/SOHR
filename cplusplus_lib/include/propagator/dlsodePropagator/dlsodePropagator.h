#ifndef __DLSODEPROPAGATOR_H_
#define __DLSODEPROPAGATOR_H_

#include "../../tools/misc/misc_template.h"
#include "../../relationshipParser/relationshipParser.h"

#include "../superPropagator/superPropagator.h"

#include "../../cubicSpline/interp_1d.h"

namespace propagator_sr {
	class dlsodePropagator : public superPropagator {

	public:
		dlsodePropagator(std::vector<double> uncertainties, std::string cwd_in);
		~dlsodePropagator();

	public:
		//template method
		virtual void propagate_pgt() override;
		/*
		* read time, pressure, temperature, concentration, drc from file
		*/
		void propagator_from_file(std::string tag = "") override;
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
		* constant pressure-->cp
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
		* constant pressure-->cp
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
		* not dependent of pressure, independent of pressure
		* Read in the initial condition and kinetic law, solve the first order equations, with end_time as stopping criteria
		* s2m means storing to memory
		*
		* This solver can be used to solve general kinetic equations, which maybe have no relation with temperature, pressure or surface
		* Like Lotka-Volterra model, which turns out to be a set of ODE. This model comprises of a simple mechanism
		* http://www.idea.wsu.edu/OscilChem/
		*
		*/
		void time_propagator_s_ct_np_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;
		//cc-->constant concentration, hold the concentration of the first species to be constant
		void time_propagator_s_ct_np_cc1_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;
		//cc-->constant concentration, hold the concentration of the first two species to be constant
		void time_propagator_s_ct_np_cc2_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time) override;

	};

}//namespace propagator_sr

#endif
