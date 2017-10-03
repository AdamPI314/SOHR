#ifndef __SOHRPROPAGATOR_H_
#define __SOHRPROPAGATOR_H_

#include "../../tools/misc/misc_template.h"
#include "../../relationshipParser/relationshipParser.h"

#include "../superPropagator/superPropagator.h"
#include "../../srkin/srkin.h"


#include "../../cubicSpline/interp_1d.h"


//Notice in this ODESlover problem, we are saving molar concentration instead of molar fraction
//The reason, we use molar concentration to calculate first order rate coefficient directly instead of molar fraction
//We need distinguish molar concentration and molar fraction
//if is a surface problem, constant pressure, constant temperature, then concentration_data_pgt store the molar concentration
//if it is not a surface problem, like reaction coefficient depends on pressure, then concentration_data_pgt store the molar fraction
//should make them to be consistent though
namespace propagator_sr {
	class SOHRPropagator : public superPropagator, public srkin_sr::srkin {
	protected:
		//number of time points
		//before critical time
		std::size_t timeN1;
		//after critical time
		std::size_t timeN2;
		//conversion factor from number of species to concentration
		double number2Concentration;
		//negative conversion factor
		double negative_number2Concentration;
		//number of trajectory
		std::size_t trajectoryNumber;
		//number of iterations
		std::size_t iterationNumber;
		//conversion factor from pathway prob to concentration
		std::vector<double> P2C;

		/*
		* only for sohrPropagator use
		*/
	public:
		//index vector, Spline_interp doesn't support integer, use double instead here
		std::vector<double> index_data_pgt;
		//index of element in vector of time versus time
		//given time, gives index
		Spline_interp* time_index_pgt;

		//private:
		//	//index vector, Spline_interp doesn't support integer, use double instead here
		//	std::vector<double> index_data_pgt;
		//	//index of element in vector of time versus time
		//	//given time, gives index
		//	Spline_interp* time_index_pgt;

	public:
		SOHRPropagator(std::vector<double> uncertainties, std::string cwd_in);
		~SOHRPropagator();

	public:
		//initialize index vs. time
		bool initiate_time_index();
		std::size_t evaluate_index_at_time(rsp::my_time_t in_time);
		/*
		* initialize concentration based initial guess, given critical time and end_time
		* initialize vector->time, conc, temperature, pressure, drc, rxn_rate
		* notice here that if initial concentration is zero, then the pseudo-first order rate coefficient
		* or says propensity function might be zero, when invert the integral_of_K, might throw error, to prevent this from
		* happening, we set the initial guesses to be a small number if they are zeros
		*/
		void time_initializer_based_on_guess(rsp::my_time_t critical_time, rsp::my_time_t end_time);
		/*at const volume, the input are mole fraction*/
		void time_initializer_based_on_guess_cv(rsp::my_time_t critical_time, rsp::my_time_t end_time);
		/*
		* initialize concentration from file, for which concentration may be calculated from lsode
		* initialize vector->time, conc, temperature, pressure, drc, rxn_rate
		* input-->
		* 	filename of time
		* 	filename of concentration
		*/
		void time_initializer_from_file(std::string filename_time, std::string filename_conc);
		/*
		* lsode initializer, solve for trajectory using losode
		* const volume
		*/
		void time_initializer_using_dlsode_cv(rsp::my_time_t critical_time, rsp::my_time_t end_time);
		/*
		* lsode initializer, solve for trajectory using losode
		* surface reaction, constant temperature, np pressure, keep the concentration of the first species constant
		*/
		/* 
		 * constant volume and constant temperature
		 */
		void time_initializer_using_dlsode_cv_ct(rsp::my_time_t critical_time, rsp::my_time_t end_time);
		void time_initializer_using_dlsode_s_ct_np(rsp::my_time_t critical_time, rsp::my_time_t end_time);
		void time_initializer_using_dlsode_s_ct_np_cc1(rsp::my_time_t critical_time, rsp::my_time_t end_time);

		/*
		* for pathway based ODE solver
		*/
		void update_spe_concentration_at_time(std::size_t time_index, std::size_t spe_index, double delta_conc);
		//doesn't include the initial time point and end_time point
		void update_spe_concentration_at_time_range(std::size_t init_time_index, std::size_t end_time_index, std::size_t spe_index, double delta_conc);

		/*
		* pressure is not a related factor here, everything is independent of pressure
		* surface reaction -->s
		* constant temperature-->ct
		*/
		void update_spe_drc_based_on_spe_concentration_s_ct_np();
		void update_reaction_rate_based_on_spe_concentration_s_ct_np();

		/*
		* constant volume-->cv
		*/
		void update_pressure_based_on_spe_concentration_cv();
		void update_temperature_pressure_based_on_spe_concentration_cv();
		void update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

		/*
		* set conc_data_sr vectors to zero
		*/
		void set_concentration_data_zero();

		/*
		* multiply every conc_data_sr element by a factor
		*/
		void rescale_concentration_data(double factor);
		void rescale_concentration_data(const std::vector<double> &factor_v);
		void rescale_prob_matrix_data(std::vector<std::vector<double> > &prob_Mat, const double factor);
		void rescale_prob_matrix_data(std::vector<std::vector<double> > &prob_Mat, const std::vector<double> &factor_v);
		//normalize prob_matrix_data
		void normalize_prob_matrix_data(std::vector<std::vector<double> > &prob_Mat);

		/*
		* set conc data from species based probability matrix
		*/
		void update_concentration_data_from_spe_based_probability_matrix(const std::vector<std::vector<double> > &prob_Mat);

	public:
		//initialize lsode, read lsode settings
		void initialize_sohr();


	public:
		//template method
		virtual void propagate_pgt() override;
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
