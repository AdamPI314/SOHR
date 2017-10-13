#ifndef __SSAPROPAGATOR_CPP_
#define __SSAPROPAGATOR_CPP_

#include "../../../include/propagator/ssaPropagator/ssaPropagator.h"
#include "../../../include/mechanism/mechanism.h"

#include "../../../include/tools/misc/fortran_routine_block_alias.h"
#include "../../../include/tools/misc/global_extern_vars.h"

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

#if defined(__CHEMKIN_AVAILABLE_) && defined(__LSODE_AVAILABLE_)
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

#endif // defined(__CHEMKIN_AVAILABLE_) && defined(__LSODE_AVAILABLE_)


	}


#if defined(__CHEMKIN_AVAILABLE_) && defined(__LSODE_AVAILABLE_)

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
		int num_reaction = this->reaction_v_sk.size();

		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ style index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
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
				update_temporary_data_pgt(nkk, num_reaction,
					ti, Temp, Pressure,
					c_t, reaction_rate_v_tmp);
			}
			else if ((print_count%deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, num_reaction,
					ti, Temp, Pressure,
					c_t, reaction_rate_v_tmp);
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

	void ssaPropagator::update_temporary_data_pgt(const int nkk, const int num_reaction, const double ti, const double Temp, const double Pressure, const double * const c_t, const std::vector<double>& reaction_rate_v_tmp)
	{
		//use the default dt to print out stuff
		time_data_pgt.push_back(ti);

		//destruction rate const
		for (int i = 0; i < nkk; ++i) {
			if (c_t[i] == 0)
				spe_drc_data_pgt[i].push_back(0.0);
			else
				spe_drc_data_pgt[i].push_back(this->cal_spe_destruction_rate(i) / c_t[i]);
		}

		for (int i = 0; i < num_reaction; i++)
			reaction_rate_data_pgt[i].push_back(reaction_rate_v_tmp[i]);

		//print concentration and temperature
		for (int i = 0; i < nkk; ++i)
			concentration_data_pgt[i].push_back(c_t[i]);

		temperature_data_pgt.push_back(Temp);
		pressure_data_pgt.push_back(Pressure);
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
		int num_reaction = this->reaction_v_sk.size();

		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ style index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
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

		int fixed_size = this->pgt_pt.get<int>("ssa_init.fixed_size");
		double order_parameter_ratio = this->pgt_pt.get<double>("ssa_init.order_parameter_ratio");
		std::list<double> time_data_list_pgt(fixed_size, 0.0);
		std::list<double> temperature_data_list_pgt(fixed_size, 0.0);
		std::list<double> pressure_data_list_pgt(fixed_size, 0.0);

		std::vector<std::list<double> > concentration_data_list_pgt(nkk, std::list<double>(fixed_size, 0.0));
		std::vector<std::list<double> > reaction_rate_data_list_pgt(num_reaction, std::list<double>(fixed_size, 0.0));
		std::vector<std::list<double> > spe_drc_data_list_pgt(nkk, std::list<double>(fixed_size, 0.0));

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
				update_temporary_data_pgt(
					time_data_list_pgt,
					temperature_data_list_pgt,
					pressure_data_list_pgt,
					concentration_data_list_pgt,
					reaction_rate_data_list_pgt,
					spe_drc_data_list_pgt,
					nkk, num_reaction,
					ti, Temp, Pressure,
					c_t, reaction_rate_v_tmp);
			}
			else if ((print_count%deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(
					time_data_list_pgt,
					temperature_data_list_pgt,
					pressure_data_list_pgt,
					concentration_data_list_pgt,
					reaction_rate_data_list_pgt,
					spe_drc_data_list_pgt,
					nkk, num_reaction,
					ti, Temp, Pressure,
					c_t, reaction_rate_v_tmp);
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

		} while ((end_time - time_data_list_pgt.back()) > 0.001*dt && final_state > -1 * initial_state*order_parameter_ratio);
		// data of the last time step is not saved, save it here
		update_temporary_data_pgt(
			time_data_list_pgt,
			temperature_data_list_pgt,
			pressure_data_list_pgt,
			concentration_data_list_pgt,
			reaction_rate_data_list_pgt,
			spe_drc_data_list_pgt,
			nkk, num_reaction, ti, Temp, Pressure,
			c_t, reaction_rate_v_tmp);

		//time, temperature, pressure, concentration, reaction rate and drc
		this->time_data_pgt.assign(time_data_list_pgt.begin(), time_data_list_pgt.end());
		for (auto &x : this->time_data_pgt)
			x -= time_data_list_pgt.front();

		this->temperature_data_pgt.assign(temperature_data_list_pgt.begin(), temperature_data_list_pgt.end());
		this->pressure_data_pgt.assign(pressure_data_list_pgt.begin(), pressure_data_list_pgt.end());
		for (int i = 0; i < nkk; ++i)
			this->concentration_data_pgt[i].assign(concentration_data_list_pgt[i].begin(), concentration_data_list_pgt[i].end());
		for (int i = 0; i < num_reaction; ++i)
			this->reaction_rate_data_pgt[i].assign(reaction_rate_data_list_pgt[i].begin(), reaction_rate_data_list_pgt[i].end());
		for (int i = 0; i < nkk; ++i)
			this->spe_drc_data_pgt[i].assign(spe_drc_data_list_pgt[i].begin(), spe_drc_data_list_pgt[i].end());

		delete[] c_t;
	}

	void ssaPropagator::update_temporary_data_pgt(std::list<double>& time_list, std::list<double>& temperature_list, std::list<double>& pressure_list, std::vector<std::list<double>>& concentration_list, std::vector<std::list<double>>& reaction_rate_list, std::vector<std::list<double>>& spe_drc_list, const int nkk, const int num_reaction, const double ti, const double Temp, const double Pressure, const double * const c_t, const std::vector<double>& reaction_rate_v_tmp)
	{
		//print out every * dt
		time_list.pop_front();
		time_list.push_back(ti);

		//destruction rate const
		for (int i = 0; i < nkk; ++i) {
			if (c_t[i] == 0)
			{
				spe_drc_list[i].pop_front();
				spe_drc_list[i].push_back(0.0);
			}
			else
			{
				spe_drc_list[i].pop_front();
				spe_drc_list[i].push_back(this->cal_spe_destruction_rate(i) / c_t[i]);
			}
		}

		for (int i = 0; i < num_reaction; i++)
		{
			reaction_rate_list[i].pop_front();
			reaction_rate_list[i].push_back(reaction_rate_v_tmp[i]);
		}

		//print concentration and temperature
		for (int i = 0; i < nkk; ++i) {
			concentration_list[i].pop_front();
			concentration_list[i].push_back(c_t[i]);
		}

		temperature_list.pop_front();
		temperature_list.push_back(Temp);
		pressure_list.pop_front();
		pressure_list.push_back(Pressure);
	}

	void ssaPropagator::time_propagator_s_ct_np_cc1_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

	void ssaPropagator::time_propagator_s_ct_np_cc2_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

#endif // defined(__CHEMKIN_AVAILABLE_) && defined(__LSODE_AVAILABLE_)


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
