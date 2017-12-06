#ifndef __DLSODEPROPAGATOR_CPP_
#define __DLSODEPROPAGATOR_CPP_

#include <boost/filesystem.hpp>

#include "../../../include/propagator/dlsodePropagator/dlsodePropagator.h"
#include "../../../include/mechanism/mechanism.h"
#include "../../../include/fileIO/CSVFileReader/CSVReader.h"
#include "../../../include/fileIO/fileIO/fileIO.h"

#include "../../../include/tools/misc/fortran_routine_block_alias.h"
#include "../../../include/tools/misc/global_extern_vars.h"


namespace propagator_sr {

	dlsodePropagator::dlsodePropagator(std::vector<double> uncertainties, std::string cwd_in) :superPropagator(uncertainties, cwd_in)
	{
		//defer this to actual use time
		propagate_and_initiate_cublic_spline_pgt();

	}

	dlsodePropagator::~dlsodePropagator()
	{
	}

	void dlsodePropagator::propagate_pgt()
	{
		if (pgt_pt.get<std::string>("propagator.primary_type") == std::string("from_file")) {
			// search file, if found, read and return, can not find, keep going
			if (boost::filesystem::exists(this->cwd_pgt + "/output/time_dlsode_M.csv")) {
				propagator_from_file(std::string("_M"));
				integrate_propensity_function_pgt();
				return;
			}
			else if (boost::filesystem::exists(this->cwd_pgt + "/output/time_dlsode_fraction.csv")) {
				propagator_from_file(std::string("_fraction"));
				integrate_propensity_function_pgt();
				convert_mole_fraction_to_molar_concentration();
				return;
			}
		}

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
		else if (sub_job_type == std::string("time_propagator_s_ct_np_cc1_s2m_pgt"))
			time_propagator_s_ct_np_cc1_s2m_pgt(uncertainties, this->pgt_pt.get<double>("time.critical_time"), this->pgt_pt.get<double>("time.max_time"));
		else if (sub_job_type == std::string("time_propagator_s_ct_np_cc2_s2m_pgt"))
			time_propagator_s_ct_np_cc2_s2m_pgt(uncertainties, this->pgt_pt.get<double>("time.critical_time"), this->pgt_pt.get<double>("time.max_time"));

		if (pgt_pt.get<std::string>("propagator.convert_molar_concentration_to_mole_fraction") == std::string("yes"))
			convert_molar_concentration_to_mole_fraction();

	}

	void dlsodePropagator::propagator_from_file(std::string tag)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		std::size_t topN = std::numeric_limits<std::size_t>::max() - 100;
		// time
		std::vector<std::vector<double> > time_Mat;
		time_Mat = fileIO::fileIO::read_topN_line_csv_matrix(this->cwd_pgt + "/output/time_dlsode" + tag + ".csv", topN);
		this->time_data_pgt.resize(time_Mat.size());
		for (std::size_t i = 0; i < time_Mat.size(); ++i)
			this->time_data_pgt[i] = time_Mat[i][0];
		time_Mat.resize(0); time_Mat.clear();
		// pressure
		std::vector<std::vector<double> > pressure_Mat;
		pressure_Mat = fileIO::fileIO::read_topN_line_csv_matrix(this->cwd_pgt + "/output/pressure_dlsode" + tag + ".csv", topN);
		this->pressure_data_pgt.resize(pressure_Mat.size());
		for (std::size_t i = 0; i < pressure_Mat.size(); ++i)
			this->pressure_data_pgt[i] = pressure_Mat[i][0];
		pressure_Mat.resize(0); pressure_Mat.clear();
		// temperature
		std::vector<std::vector<double> > temperature_Mat;
		temperature_Mat = fileIO::fileIO::read_topN_line_csv_matrix(this->cwd_pgt + "/output/temperature_dlsode" + tag + ".csv", topN);
		this->temperature_data_pgt.resize(temperature_Mat.size());
		for (std::size_t i = 0; i < temperature_Mat.size(); ++i)
			this->temperature_data_pgt[i] = temperature_Mat[i][0];
		temperature_Mat.resize(0); temperature_Mat.clear();
		// concentration
		std::vector<std::vector<double> > concentration_Mat;
		concentration_Mat = fileIO::fileIO::read_topN_line_csv_matrix(this->cwd_pgt + "/output/concentration_dlsode" + tag + ".csv", topN);
		this->concentration_data_pgt.resize(concentration_Mat[0].size());
		for (std::size_t i = 0; i < concentration_Mat[0].size(); ++i) {
			this->concentration_data_pgt[i].resize(concentration_Mat.size());
			for (std::size_t j = 0; j < concentration_Mat.size(); ++j) {
				this->concentration_data_pgt[i][j] = concentration_Mat[j][i];
			}
		}
		concentration_Mat.resize(0); concentration_Mat.clear();
		// drc
		std::vector<std::vector<double> > drc_Mat;
		drc_Mat = fileIO::fileIO::read_topN_line_csv_matrix(this->cwd_pgt + "/output/drc_dlsode" + tag + ".csv", topN);
		this->spe_drc_data_pgt.resize(drc_Mat[0].size());
		for (std::size_t i = 0; i < drc_Mat[0].size(); ++i) {
			this->spe_drc_data_pgt[i].resize(drc_Mat.size());
			for (std::size_t j = 0; j < drc_Mat.size(); ++j) {
				this->spe_drc_data_pgt[i][j] = drc_Mat[j][i];
			}
		}
		drc_Mat.resize(0); drc_Mat.clear();
		// reaction rate
		std::vector<std::vector<double> > reaction_rate_Mat;
		reaction_rate_Mat = fileIO::fileIO::read_topN_line_csv_matrix(this->cwd_pgt + "/output/reaction_rate_dlsode" + tag + ".csv", topN);
		this->reaction_rate_data_pgt.resize(reaction_rate_Mat[0].size());
		for (std::size_t i = 0; i < reaction_rate_Mat[0].size(); ++i) {
			this->reaction_rate_data_pgt[i].resize(reaction_rate_Mat.size());
			for (std::size_t j = 0; j < reaction_rate_Mat.size(); ++j) {
				this->reaction_rate_data_pgt[i][j] = reaction_rate_Mat[j][i];
			}
		}
		reaction_rate_Mat.resize(0); reaction_rate_Mat.clear();


	}


	void dlsodePropagator::time_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();
		//spe_production_rate_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());


		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		// xgst[0]~xgst[nkk-1] are the MASS FRACTION of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();

		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double* x_t = new double[nkk]; double* y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));
		this->read_configuration(Temp, ckstore.pressure, neq - 1, x_t);


		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		double dt = 0;
		this->initialize_lsode(dt);

		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = y_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;
		//molar concentration.
		double* c_t = new double[nkk];
		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);
		//spe_production_rate_data_pgt.resize(nkk);

		//while (tout < end_time)
		do
		{
			//convert mass fractions to molar fractions
			mechanism::kinetics::ytx(y_t, x_t);
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::py(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);
			//molar concentration
			mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::kfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//destruction relative rate Constant of species
			//[ print out
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN2 == 0)) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if ((print_Count%lsodestore.deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				//print out every * dt
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodev(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodav(&ti, &tout, &neq, xgst);
			//update mass fractions
			for (int i = 0; i < nkk; ++i)
				y_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		} while ((end_time - time_data_pgt.back()) > 0.001*dt);

		delete[] x_t; delete[] y_t; delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::time_propagator_cv_ct_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();
		//spe_production_rate_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());

		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		// xgst[0]~xgst[nkk-1] are the MASS FRACTION of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();

		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double* x_t = new double[nkk]; double* y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));
		this->read_configuration(Temp, ckstore.pressure, neq - 1, x_t);


		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		double dt = 0;
		this->initialize_lsode(dt);

		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = y_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;
		//molar concentration.
		double* c_t = new double[nkk];
		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);
		//spe_production_rate_data_pgt.resize(nkk);

		//while (tout < end_time)
		do
		{
			//convert mass fractions to molar fractions
			mechanism::kinetics::ytx(y_t, x_t);
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::py(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);
			//molar concentration
			mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::kfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//[ print out
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN2 == 0)) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if ((print_Count%lsodestore.deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodevt(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodavt(&ti, &tout, &neq, xgst);
			//update mass fractions
			for (int i = 0; i < nkk; ++i)
				y_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		} while ((end_time - time_data_pgt.back()) > 0.001*dt);

		delete[] x_t; delete[] y_t; delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::temperature_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t, double end_temperature)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();
		//spe_production_rate_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());


		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		// xgst[0]~xgst[nkk-1] are the MASS FRACTION of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();

		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double* x_t = new double[nkk]; double* y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));
		this->read_configuration(Temp, ckstore.pressure, neq - 1, x_t);


		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		double dt = 0;
		this->initialize_lsode(dt);

		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = y_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;
		//molar concentration.
		double* c_t = new double[nkk];
		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);
		//spe_production_rate_data_pgt.resize(nkk);

		//while (Temp<end_temperature)
		while ((Temp < end_temperature) || (temperature_data_pgt.back() < end_temperature))
		{
			//convert mass fractions to molar fractions
			mechanism::kinetics::ytx(y_t, x_t);
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::py(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);
			mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::kfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//[ print out
			if ((Temp >= critical_temperature_t) && (print_Count%lsodestore.deltaN2 == 0)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if (print_Count%lsodestore.deltaN1 == 0) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodev(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodav(&ti, &tout, &neq, xgst);
			//update mass fractions
			for (int i = 0; i < nkk; ++i)
				y_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		}//while

		delete[] x_t; delete[] y_t; delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::equilibrium_propagator_cv_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());

		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		// xgst[0]~xgst[nkk-1] are the MASS FRACTION of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double* x_t = new double[nkk]; double* y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp, Temp_t = 0.0;
		read_configuration(Temp, ckstore.pressure, neq - 1, x_t);
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));

		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		//time step
		double dt = 0;
		initialize_lsode(dt);
		//	lsode_init(dt, this->cwd_dl+std::string("/input/setting.cfg"));

		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = y_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;
		//molar concentration.
		double* c_t = new double[nkk];
		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);

		//while (Temp<end_temperature)
		//while ((Temp<end_temperature) || (temperature_data_pgt.back()<end_temperature))
		while (Temp != Temp_t)
		{
			Temp_t = Temp;
			//convert mass fractions to molar fractions
			mechanism::kinetics::ytx(y_t, x_t);
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::py(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);
			mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::kfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//[ print out
			if ((Temp >= critical_temperature_t) && (print_Count%lsodestore.deltaN2 == 0)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if (print_Count%lsodestore.deltaN1 == 0) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodev(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodav(&ti, &tout, &neq, xgst);
			//update mass fractions
			for (int i = 0; i < nkk; ++i)
				y_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;

		}//while

		delete[] x_t; delete[] y_t; delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::time_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());

		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		// xgst[0]~xgst[nkk-1] are the MASS FRACTION of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double* x_t = new double[nkk]; double* y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		read_configuration(Temp, ckstore.pressure, neq - 1, x_t);
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));

		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		//time step
		double dt = 0;
		initialize_lsode(dt);
		//	lsode_init(dt, this->cwd_dl+std::string("/input/setting.cfg"));

		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant pressure condition.
		//So 'ckstore.pressure' is consant, 'ckstore.rhomass' is not in the simulation.
		mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = y_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;
		//molar concentration.
		double* c_t = new double[nkk];
		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);

		while (tout < end_time)
		{
			//convert mass fractions to molar fractions
			mechanism::kinetics::ytx(y_t, x_t);
			//Returns the 'ckstore.rhomass' of the gas mixture given mass density, temperature(s) and pressure
			mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);
			mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::kfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//destruction relative rate Constant of species
			//[ print out
			//int NPrecision=16;
			if ((tout >= critical_time) && (print_Count%lsodestore.deltaN2 == 0)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if (print_Count%lsodestore.deltaN1 == 0) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodep(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodap(&ti, &tout, &neq, xgst);
			//update mass fractions
			for (int i = 0; i < nkk; ++i)
				y_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		}//while

		delete[] x_t; delete[] y_t; delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::temperature_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t, double end_temperature)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());


		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		// xgst[0]~xgst[nkk-1] are the MASS FRACTION of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double* x_t = new double[nkk]; double* y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		read_configuration(Temp, ckstore.pressure, neq - 1, x_t);

		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		//time step
		double dt = 0;
		initialize_lsode(dt);
		//	lsode_init(dt, this->cwd_dl+std::string("/input/setting.cfg"));

		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant pressure condition.
		//So 'ckstore.pressure' is consant, 'ckstore.rhomass' is not in the simulation.
		mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = y_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;
		//molar concentration.
		double* c_t = new double[nkk];
		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);

		//while (Temp<end_temperature)
		while ((Temp < end_temperature) || (temperature_data_pgt.back() < end_temperature))
		{
			//convert mass fractions to molar fractions
			mechanism::kinetics::ytx(y_t, x_t);
			//Returns the 'ckstore.rhomass' of the gas mixture given mass density, temperature(s) and pressure
			mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);
			mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::kfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//destruction relative rate Constant of species
			//[ print out
			//int NPrecision=16;
			//std::fill(vec_for_rr_tmp.begin(), vec_for_rr_tmp.end(), 0);
			if ((Temp >= critical_temperature_t) && (print_Count%lsodestore.deltaN2 == 0)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if (print_Count%lsodestore.deltaN1 == 0) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodep(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodap(&ti, &tout, &neq, xgst);
			//update mass fractions
			for (int i = 0; i < nkk; ++i)
				y_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		}//while

		delete[] x_t; delete[] y_t; delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::equilibrium_propagator_cp_s2m_pgt(std::vector<double> uncertainties, double critical_temperature_t)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());

		//set new rate constant
		int I_t = 1;	double R_A = 0.0;
		//cout<<"\nPre-exponential Constant, old and new: "<<endl;
		for (; static_cast<size_t>(I_t) <= uncertainties.size(); ++I_t) {
			mechanism::kinetics::raex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::raex(&I_t, &R_A);
			////To see whether it changed or not
			I_t = abs(I_t);
			//chemkin::ckraex(&I_t, &R_A);
			//cout<<R_A<<endl;
		}

		// xgst[0]~xgst[nkk-1] are the MASS FRACTION of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double* x_t = new double[nkk]; double* y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp, Temp_t = 0.0;
		read_configuration(Temp, ckstore.pressure, neq - 1, x_t);
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));

		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		//time step
		double dt = 0;
		initialize_lsode(dt);
		//	lsode_init(dt, this->cwd_dl+std::string("/input/setting.cfg"));

		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant pressure condition.
		//So 'ckstore.pressure' is consant, 'ckstore.rhomass' is not in the simulation.
		mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = y_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;
		//molar concentration.
		double* c_t = new double[nkk];
		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);

		//while (Temp<end_temperature)
		//while ((Temp<end_temperature) || (temperature_data_pgt.back()<end_temperature))
		while (Temp != Temp_t)
		{
			Temp_t = Temp;
			//convert mass fractions to molar fractions
			mechanism::kinetics::ytx(y_t, x_t);
			//Returns the 'ckstore.rhomass' of the gas mixture given mass density, temperature(s) and pressure
			mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);
			mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::kfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//[ print out
			if ((Temp >= critical_temperature_t) && (print_Count%lsodestore.deltaN2 == 0)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if (print_Count%lsodestore.deltaN1 == 0) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodep(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodap(&ti, &tout, &neq, xgst);
			//update mass fractions
			for (int i = 0; i < nkk; ++i)
				y_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;

		}//while

		delete[] x_t; delete[] y_t; delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::time_propagator_cp_ct_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());

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

		// xgst[0]~xgst[nkk-1] are the MASS FRACTION of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double* x_t = new double[nkk]; double* y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		read_configuration(Temp, ckstore.pressure, neq - 1, x_t);
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));

		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		//time step
		double dt = 0;
		initialize_lsode(dt);
		//	lsode_init(dt, this->cwd_dl+std::string("/input/setting.cfg"));

		//convert mole fractions to mass fractions
		mechanism::kinetics::xty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant pressure condition.
		//So 'ckstore.pressure' is consant, 'ckstore.rhomass' is not in the simulation.
		mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = y_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;
		//molar concentration.
		double* c_t = new double[nkk];
		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);

		//while (tout < end_time)
		do
		{
			//convert mass fractions to molar fractions
			mechanism::kinetics::ytx(y_t, x_t);
			//Returns the 'ckstore.rhomass' of the gas mixture given mass density, temperature(s) and pressure
			mechanism::kinetics::rhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);
			mechanism::kinetics::ytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::cdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::kfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//destruction relative rate Constant of species
			//[ print out
			//int NPrecision=16;
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN2 == 0)) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if ((print_Count%lsodestore.deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodept(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodapt(&ti, &tout, &neq, xgst);
			//update mass fractions
			for (int i = 0; i < nkk; ++i)
				y_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		} while ((end_time - time_data_pgt.back()) > 0.001*dt);

		delete[] x_t; delete[] y_t; delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::time_propagator_s_ct_np_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());

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

		// xgst[0]~xgst[nkk-1] are the MASS molar concentration of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar concentration.
		double* c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));
		read_configuration(Temp, ckstore.pressure, neq - 1, c_t);
		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		//time step
		double dt = 0;
		//	lsode_init(dt, this->cwd_dl+std::string("/input/setting.cfg"));
		initialize_lsode(dt);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = c_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;

		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);

		//"imaginary pressure"
		//double i_pressure=0.0;

		//while (tout < end_time)
		do
		{
			//Returns the molar creation and destruction rates of the species given temperature(s) and molar concentration
			mechanism::kinetics::cdc(&Temp, c_t, CDOT_t, DDOT_t);

			// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
			// Applicable for reactions with rate constant independent of pressure
			// where sr stands for Shirong
			mechanism::kinetics::kfkrsr(&Temp, c_t, FWDR_t, REVR_t);

			//[ print out
			//int NPrecision=16;
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN2 == 0)) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if ((print_Count%lsodestore.deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodest(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodast(&ti, &tout, &neq, xgst);
			//update  molar concentration
			for (int i = 0; i < nkk; ++i)
				c_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		} while ((end_time - time_data_pgt.back()) > 0.001*dt);

		//delete[] x_t;
		delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::time_propagator_s_ct_np_cc1_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());

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

		// xgst[0]~xgst[nkk-1] are the MASS molar concentration of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar concentration.
		double* c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));
		read_configuration(Temp, ckstore.pressure, neq - 1, c_t);
		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		//time step
		double dt = 0;
		//	lsode_init(dt, this->cwd_dl+std::string("/input/setting.cfg"));
		initialize_lsode(dt);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = c_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;

		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);

		//while (tout < end_time)
		do
		{
			//Returns the molar creation and destruction rates of the species given temperature(s) and molar concentration
			//The difference is how to treat the auto-catylytic reactions
			mechanism::kinetics::cdc(&Temp, c_t, CDOT_t, DDOT_t);
			//this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);

			// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
			// Applicable for reactions with rate constant independent of pressure
			// where sr stands for Shirong
			mechanism::kinetics::kfkrsr(&Temp, c_t, FWDR_t, REVR_t);

			//[ print out
			//int NPrecision=16;
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN2 == 0)) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if ((print_Count%lsodestore.deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodest(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodastcc1(&ti, &tout, &neq, xgst);
			//update  molar concentration
			for (int i = 0; i < nkk; ++i)
				c_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		} while ((end_time - time_data_pgt.back()) > 0.001*dt);

		//delete[] x_t;
		delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}

	void dlsodePropagator::time_propagator_s_ct_np_cc2_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(), reactionNetwork_chemkin_index_map.end());

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

		// xgst[0]~xgst[nkk-1] are the MASS molar concentration of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq = mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double* xgst = new double[neq];	for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar concentration.
		double* c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp;
		//	chem_init(Temp, ckstore.pressure, neq-1, x_t, this->cwd_dl+std::string("/input/setting.cfg"));
		read_configuration(Temp, ckstore.pressure, neq - 1, c_t);
		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		// Read time step 'dt(s)'.
		// Read 'lsode' parameters.
		//time step
		double dt = 0;
		//	lsode_init(dt, this->cwd_dl+std::string("/input/setting.cfg"));
		initialize_lsode(dt);

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = c_t[i];

		double ti = 0.0, tout = ti + dt;
		int print_Count = 0;

		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		spe_drc_data_pgt.resize(nkk);

		//while (tout < end_time)
		do
		{
			//Returns the molar creation and destruction rates of the species given temperature(s) and molar concentration
			//The difference is how to treat the auto-catylytic reactions
			mechanism::kinetics::cdc(&Temp, c_t, CDOT_t, DDOT_t);
			//this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);

			// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
			// Applicable for reactions with rate constant independent of pressure
			// where sr stands for Shirong
			mechanism::kinetics::kfkrsr(&Temp, c_t, FWDR_t, REVR_t);

			//[ print out
			//int NPrecision=16;
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN2 == 0)) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			else if ((print_Count%lsodestore.deltaN1 == 0) || ((end_time - ti) < 0.001*dt)) {
				update_temporary_data_pgt(nkk, neq, ti, c_t, CDOT_t, DDOT_t, FWDR_t, REVR_t, xgst);
			}
			//] print out

			//ODE::solver::cppdlsodest(&ti, &tout, &neq, xgst);
			ODE::solver::cppdlsodastcc2(&ti, &tout, &neq, xgst);
			//update  molar concentration
			for (int i = 0; i < nkk; ++i)
				c_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];
			//tout=ti+dt;
			ti = tout;
			tout += dt;
			//in case print_Count is too big
			print_Count = (print_Count + 1) % 1000000;
		} while ((end_time - time_data_pgt.back()) > 0.001*dt);

		//delete[] x_t;
		delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
	}


}


#endif
