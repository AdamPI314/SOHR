#ifndef __SSAPROPAGATOR_CPP_
#define __SSAPROPAGATOR_CPP_

#include "ssaPropagator.h"
#include "../../chemkinCpp/chemkincpp.h"

#include "../../tools/misc/fortran_routine_block_alias.h"
#include "../../tools/misc/global_extern_vars.h"

#include <numeric>
#include <climits>

namespace propagator_sr {

	ssaPropagator::ssaPropagator(std::vector<double> uncertainties, std::size_t random_seed_for_this_core, std::string cwd_in)
		:superPropagator(uncertainties, cwd_in), srkin(cwd_in, 0.0, 0.0, 0.0)
	{
		//random seed for this core
		this->random_seed_for_this_core = static_cast<boost::uint32_t>(random_seed_for_this_core);
		rand = new random_sr::random(this->random_seed_for_this_core);

		//defer this to actual use time
		propagate_and_initiate_cublic_spline_pgt();

	}

	ssaPropagator::~ssaPropagator()
	{
		delete rand;
	}

	void ssaPropagator::propagate_pgt()
	{
		std::string sub_job_type = pgt_pt.get<std::string>("propagator.sub_type");
		if (sub_job_type == std::string("temperature_propagator_cv_s2m_pgt"))
			temperature_propagator_cv_s2m_pgt(uncertainties, this->pgt_pt.get<double>("T.critical_temperature"), this->pgt_pt.get<double>("T.end_temperature"));
		else if (sub_job_type == std::string("time_propagator_cv_s2m_pgt"))
			time_propagator_cv_s2m_pgt(uncertainties, this->pgt_pt.get<rsp::my_time_t>("time.critical_time"), this->pgt_pt.get<rsp::my_time_t>("time.max_time"));
		else if (sub_job_type == std::string("time_propagator_cv_ct_s2m_pgt"))
			time_propagator_cv_ct_s2m_pgt(uncertainties, this->pgt_pt.get<rsp::my_time_t>("time.critical_time"), this->pgt_pt.get<rsp::my_time_t>("time.max_time"));
		else if (sub_job_type == std::string("equilibrium_propagator_cv_s2m_pgt"))
			equilibrium_propagator_cv_s2m_pgt(uncertainties, this->pgt_pt.get<double>("T.critical_temperature"));
		else if (sub_job_type == std::string("temperature_propagator_cp_s2m_pgt"))
			temperature_propagator_cp_s2m_pgt(uncertainties, this->pgt_pt.get<double>("T.critical_temperature"), this->pgt_pt.get<double>("T.end_temperature"));
		else if (sub_job_type == std::string("time_propagator_cp_s2m_pgt"))
			time_propagator_cp_s2m_pgt(uncertainties, this->pgt_pt.get<rsp::my_time_t>("time.critical_time"), this->pgt_pt.get<rsp::my_time_t>("time.max_time"));
		else if (sub_job_type == std::string("equilibrium_propagator_cp_s2m_pgt"))
			equilibrium_propagator_cp_s2m_pgt(uncertainties, this->pgt_pt.get<double>("T.critical_temperature"));
		else if (sub_job_type == std::string("time_propagator_cp_ct_s2m_pgt"))
			time_propagator_cp_ct_s2m_pgt(uncertainties, this->pgt_pt.get<double>("time.critical_time"), this->pgt_pt.get<double>("time.max_time"));
		else if (sub_job_type == std::string("time_propagator_s_ct_np_s2m_pgt"))
			time_propagator_s_ct_np_s2m_pgt(uncertainties, this->pgt_pt.get<double>("time.critical_time"), this->pgt_pt.get<double>("time.max_time"));
		else if (sub_job_type == std::string("time_propagator_s_ct_np_s2m_find_one_transition_pgt"))
			time_propagator_s_ct_np_s2m_find_one_transition_pgt(uncertainties, this->pgt_pt.get<double>("time.critical_time"), this->pgt_pt.get<double>("time.max_time"));

		else if (sub_job_type == std::string("time_propagator_s_ct_np_cc1_s2m_pgt"))
			time_propagator_s_ct_np_cc1_s2m_pgt(uncertainties, this->pgt_pt.get<double>("time.critical_time"), this->pgt_pt.get<double>("time.max_time"));
		else if (sub_job_type == std::string("time_propagator_s_ct_np_cc2_s2m_pgt"))
			time_propagator_s_ct_np_cc2_s2m_pgt(uncertainties, this->pgt_pt.get<double>("time.critical_time"), this->pgt_pt.get<double>("time.max_time"));

		if (pgt_pt.get<std::string>("propagator.convert_molar_concentration_to_mole_fraction") == std::string("yes"))
			convert_molar_concentration_to_mole_fraction();

	}

	void ssaPropagator::time_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

	void ssaPropagator::time_propagator_cv_ct_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

	void ssaPropagator::temperature_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t, double end_temperature)
	{
	}

	void ssaPropagator::equilibrium_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t)
	{
	}

	void ssaPropagator::time_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

	void ssaPropagator::temperature_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t, double end_temperature)
	{
	}

	void ssaPropagator::equilibrium_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t)
	{
	}

	void ssaPropagator::time_propagator_cp_ct_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

	void ssaPropagator::time_propagator_s_ct_np_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = this->reaction_v_sk.size();

		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ style index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = this->species_v_sk.size();

		//molar concentration.
		double* c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)', and MOLAR FRACTION for species.
		double Temp, Pressure;
		read_configuration(Temp, Pressure, nkk, c_t);
		//time step
		double dt = this->pgt_pt.get<double>("ssa_init.dt");
		std::size_t deltaN1 = this->pgt_pt.get<std::size_t>("ssa_init.deltaN1");
		std::size_t deltaN2 = this->pgt_pt.get<std::size_t>("ssa_init.deltaN2");

		double ti = 0.0, tau = 0.0;
		int print_count = 0;

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction); std::vector<double> reaction_rate_v_tmp(num_reaction, 0.0);
		spe_drc_data_pgt.resize(nkk);

		//initialization of srkin
		this->set_temperature(Temp);
		this->set_species_conc_v(c_t);
		this->cal_reaction_rate_discrete();

		for (size_t i = 0; i < this->reaction_v_sk.size(); i++)
		{
			reaction_rate_v_tmp[i] = this->get_reaction_rate(i);
		}

		do
		{
			//[ print out
			if (((ti >= critical_time) && (print_count%deltaN2 == 0)) || ((end_time - ti) < 0.001*dt)) {
				//use the default dt to print out stuff
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i) {
					if (c_t[i] == 0)
						spe_drc_data_pgt[i].push_back(0.0);
					else
						spe_drc_data_pgt[i].push_back(this->cal_spe_destruction_rate(i) / c_t[i]);
				}

				for (size_t i = 0; i < num_reaction; i++)
					reaction_rate_data_pgt[i].push_back(reaction_rate_v_tmp[i]);

				//print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(Temp);
				pressure_data_pgt.push_back(Pressure);
			}
			else if ((print_count%deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				//print out every * dt
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i) {
					if (c_t[i] == 0)
						spe_drc_data_pgt[i].push_back(0.0);
					else
						spe_drc_data_pgt[i].push_back(this->cal_spe_destruction_rate(i) / c_t[i]);
				}

				for (size_t i = 0; i < num_reaction; i++)
					reaction_rate_data_pgt[i].push_back(reaction_rate_v_tmp[i]);

				//print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(Temp);
				pressure_data_pgt.push_back(Pressure);
			}
			//] print out

			//update step, run until delta_t >= dt, avoid saving too many data points
			this->ssa_update_dt_s_ct_np_pgt(Temp, dt, c_t, reaction_rate_v_tmp, tau);

			ti += tau;
			//in case print_Count is too big
			print_count = (print_count + 1) % INT_MAX;
		} while ((end_time - time_data_pgt.back()) > 0.001*dt);

		delete[] c_t;
	}

	void ssaPropagator::time_propagator_s_ct_np_s2m_find_one_transition_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = this->reaction_v_sk.size();

		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ style index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = this->species_v_sk.size();

		//molar concentration.
		double* c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)', and MOLAR FRACTION for species.
		double Temp, Pressure;
		read_configuration(Temp, Pressure, nkk, c_t);

		// linear combination prefactors
		std::vector<double> pre_factor = { 0, 0, 1, -1, 2, -2, 2, -2 };
		double initial_state = 0.0, final_state = 0.0;
		for (int i = 0; i < nkk; ++i)
			initial_state += pre_factor[i] * c_t[i];
		final_state = initial_state;

		//time step
		double dt = this->pgt_pt.get<double>("ssa_init.dt");
		std::size_t deltaN1 = this->pgt_pt.get<std::size_t>("ssa_init.deltaN1");
		std::size_t deltaN2 = this->pgt_pt.get<std::size_t>("ssa_init.deltaN2");

		double ti = 0.0, tau = 0.0;
		int print_count = 0;

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction); std::vector<double> reaction_rate_v_tmp(num_reaction, 0.0);
		spe_drc_data_pgt.resize(nkk);

		//initialization of srkin
		this->set_temperature(Temp);
		this->set_species_conc_v(c_t);
		this->cal_reaction_rate_discrete();

		for (size_t i = 0; i < this->reaction_v_sk.size(); i++)
		{
			reaction_rate_v_tmp[i] = this->get_reaction_rate(i);
		}

		do
		{
			//[ print out
			if (((ti >= critical_time) && (print_count%deltaN2 == 0)) || ((end_time - ti) < 0.001*dt)) {
				//use the default dt to print out stuff
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i) {
					if (c_t[i] == 0)
						spe_drc_data_pgt[i].push_back(0.0);
					else
						spe_drc_data_pgt[i].push_back(this->cal_spe_destruction_rate(i) / c_t[i]);
				}

				for (size_t i = 0; i < num_reaction; i++)
					reaction_rate_data_pgt[i].push_back(reaction_rate_v_tmp[i]);

				//print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(Temp);
				pressure_data_pgt.push_back(Pressure);
			}
			else if ((print_count%deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				//print out every * dt
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i) {
					if (c_t[i] == 0)
						spe_drc_data_pgt[i].push_back(0.0);
					else
						spe_drc_data_pgt[i].push_back(this->cal_spe_destruction_rate(i) / c_t[i]);
				}

				for (size_t i = 0; i < num_reaction; i++)
					reaction_rate_data_pgt[i].push_back(reaction_rate_v_tmp[i]);

				//print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(Temp);
				pressure_data_pgt.push_back(Pressure);
			}
			//] print out

			//update step, run until delta_t >= dt, avoid saving too many data points
			this->ssa_update_dt_s_ct_np_pgt(Temp, dt, c_t, reaction_rate_v_tmp, tau);

			ti += tau;
			//in case print_Count is too big
			print_count = (print_count + 1) % INT_MAX;

			// update final state
			final_state = 0.0;
			for (int i = 0; i < nkk; ++i)
				final_state += pre_factor[i] * c_t[i];

		} while ((end_time - time_data_pgt.back()) > 0.001*dt && final_state >= -1 * initial_state*0.6);

		delete[] c_t;
	}

	void ssaPropagator::time_propagator_s_ct_np_cc1_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

	void ssaPropagator::time_propagator_s_ct_np_cc2_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

	void ssaPropagator::ssa_update_next_reaction_pgt(double* c_t, const std::vector<double>& reaction_rates, double & tau)
	{
		double sum_propensity = std::accumulate(reaction_rates.begin(), reaction_rates.end(), 0.0);
		//generage the random number u_1
		double u_1 = 0.0;
		do {
			u_1 = rand->random01();
		} while (u_1 == 1.0 || u_1 == 0.0);

		tau = (-1.0 / sum_propensity) * log(u_1);

		//select raction
		std::size_t rxn_ind = this->rand->return_index_randomly_given_probability_vector(reaction_rates);

		//update species number/concentration
		for (auto spe : this->reaction_v_sk[rxn_ind].net_reactant)
			c_t[spe.first] -= spe.second;

		for (auto spe : this->reaction_v_sk[rxn_ind].net_product)
			c_t[spe.first] += spe.second;
	}

	void ssaPropagator::ssa_update_dt_s_ct_np_pgt(const double Temp, const double dt, double* c_t, std::vector<double>& reaction_rate_v_tmp, double& tau)
	{
		double tau_tmp = 0.0;
		//have to reset tau here
		tau = 0.0;
		while (tau < dt) {
			//discrete version of calculating reaction rates
			this->set_temperature(Temp);
			this->set_species_conc_v(c_t);
			this->cal_reaction_rate_discrete();

			for (std::size_t i = 0; i < this->reaction_v_sk.size(); ++i) {
				reaction_rate_v_tmp[i] = this->get_reaction_rate(i);
			}

			this->ssa_update_next_reaction_pgt(c_t, reaction_rate_v_tmp, tau_tmp);
			tau += tau_tmp;
		}

	}

}


#endif