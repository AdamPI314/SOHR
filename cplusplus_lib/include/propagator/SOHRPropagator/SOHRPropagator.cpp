#ifndef __SOHRPROPAGATOR_CPP_
#define __SOHRPROPAGATOR_CPP_

#include <numeric>
#include "SOHRPropagator.h"
#include "../../mechanism/mechanism.h"

#include "../../tools/misc/fortran_routine_block_alias.h"
#include "../../tools/misc/global_extern_vars.h"


namespace propagator_sr {

	SOHRPropagator::SOHRPropagator(std::vector<double> uncertainties, std::string cwd_in) : superPropagator(
		uncertainties, cwd_in), srkin(cwd_in, 0.0, 0.0, 0.0) {
		time_index_pgt = static_cast<Spline_interp *>(0);

		////defer this to actual use time
		//propagate_and_initiate_cublic_spline_pgt();
		initialize_sohr();

		// only constant volume
		if (this->pgt_pt.get<std::string>("system.condition") == std::string("cv")) {
			// dlsode
			if (this->pgt_pt.get<std::string>("system.initializer") == std::string("dlsode")) {
				time_initializer_using_dlsode_cv(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
			// initial guess
			else if (this->pgt_pt.get<std::string>("system.initializer") == std::string("guess")) {
				time_initializer_based_on_guess_cv(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
			// others
			else {
				time_initializer_based_on_guess_cv(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}

		}
		// constant volume and constant temperature
		else if (this->pgt_pt.get<std::string>("system.condition") == std::string("cv_ct")) {
			if (this->pgt_pt.get<std::string>("system.initializer") == std::string("dlsode")) {
				time_initializer_using_dlsode_cv_ct(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
			else {
				time_initializer_based_on_guess_cv(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));

			}

		}
		// not constant volume
		// surface reaction, constant temperature, np pressure
		else if (this->pgt_pt.get<std::string>("system.condition") == std::string("s_ct_np")) {
			if (this->pgt_pt.get<std::string>("system.initializer") == std::string("dlsode")) {
				time_initializer_using_dlsode_s_ct_np(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
			else if (this->pgt_pt.get<std::string>("system.initializer") == std::string("guess")) {
				time_initializer_based_on_guess(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
			else {
				time_initializer_based_on_guess(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
		}
		// surface reaction, constant temperature, np pressure, cc1
		else if (this->pgt_pt.get<std::string>("system.condition") == std::string("s_ct_np_cc1")) {
			if (this->pgt_pt.get<std::string>("system.initializer") == std::string("dlsode")) {
				time_initializer_using_dlsode_s_ct_np_cc1(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
			else if (this->pgt_pt.get<std::string>("system.initializer") == std::string("guess")) {
				time_initializer_based_on_guess(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
			else {
				time_initializer_based_on_guess(this->pgt_pt.get<double>("time.critical_time"),
					this->pgt_pt.get<double>("time.max_time"));
			}
		}
		// others
		else {
			time_initializer_based_on_guess(this->pgt_pt.get<double>("time.critical_time"),
				this->pgt_pt.get<double>("time.max_time"));
		}
		initiate_cubic_spline();
		initiate_time_index();
	}

	SOHRPropagator::~SOHRPropagator() {
	}


	bool SOHRPropagator::initiate_time_index() {
		//re-initialize index vector
		index_data_pgt.clear();
		for (std::size_t i = 0; i < time_data_pgt.size(); ++i) {
			index_data_pgt.push_back(i);
		}
		//create Spline_interp pointer
		if (time_index_pgt != 0) {
			std::cout << time_index_pgt << std::endl;
			delete time_index_pgt;
		}
		time_index_pgt = new Spline_interp(time_data_pgt, index_data_pgt);

		return true;
	}

	std::size_t SOHRPropagator::evaluate_index_at_time(rsp::my_time_t in_time) {
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		return static_cast<std::size_t>(round(time_index_pgt->interp(in_time)));
	}

	void SOHRPropagator::time_initializer_based_on_guess(rsp::my_time_t critical_time, rsp::my_time_t end_time) {
		//std::cout <<"deltaConcentration\t" <<this->pt.get<double>("SOHR_init.deltaConcentration") << std::endl;
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();
		//total_concentration_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(),
			reactionNetwork_chemkin_index_map.end());

		// xgst[0]~xgst[nkk-1] are the MASS molar concentration of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq =
			mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double *xgst = new double[neq];
		for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double *x_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; }

		//molar concentration.
		double *c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR CONCENTRATION for species.
		double Temp;
		read_configuration(Temp, ckstore.pressure, neq - 1, c_t);

		//recalculate the pressure using equation P= nRT
		//The pressure has to be in cgs unit, convert pascal(Pa) to dyne/square centimeter [dyn/cm**2]
		// 1atm= 101325 Pa
		// 1Pa= 10 dyn/cm**2
		// P = nRT
		double ru, ruc, pa;
		mechanism::kinetics::ckrp(&ru, &ruc, &pa);
		//ckstore.pressure = pressure_data_pgt[0];
		ckstore.pressure = std::accumulate(c_t, c_t + nkk, 0.0)*ru*Temp;

		/*
		 * notice here that if initial concentration is zero, then the pseudo-first order rate coefficient
		 * or says propensity function might be zero, when invert the integral_of_K, might throw error, to prevent this from
		 * happening, we set the initial guesses to be a small number if they are zeros
		 */
		for (int i = 0; i < nkk; ++i)
		{
			if (c_t[i] == 0.0)
				c_t[i] = this->pgt_pt.get<double>("SOHR_init.deltaConcentration");
		}

		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		//Read time step 'dt(s)'.
		//Read 'sohr' parameters.
		//initialize_sohr();

		//time step
		double dt1, dt2;

		if (critical_time == 0.0) {
			dt1 = dt2 = end_time / this->timeN2;
		}
		else if (critical_time == end_time) {
			dt1 = dt2 = end_time / (this->timeN1);
		}
		else {
			dt1 = critical_time / (this->timeN1);
			dt2 = (end_time - critical_time) / this->timeN2;
		}

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = c_t[i];

		double ti = 0.0, tout = 0.0;

		//species creation rates and destruction rates.
		double *CDOT_t = new double[nkk];
		double *DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double *FWDR_t = new double[nii];
		double *REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		double reaction_rate_tmp = 0.0;
		spe_drc_data_pgt.resize(nkk);

		do {
			//////////////////////////////////////////////////////////////////////////
			//Calculate useful stuff
			//////////////////////////////////////////////////////////////////////////

			//Returns the molar creation and destruction rates of the species given temperature(s) and molar concentration
			this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);
			this->cal_reaction_rate(&Temp, c_t, FWDR_t, REVR_t);

			if (tout >= critical_time) {
				tout += dt2;
				//use the default dt to print out stuff
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0) {
						spe_drc_data_pgt[i].push_back(c_t[i]);
					}
					else {
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
					}
				}//for]

				//Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					reaction_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] -
							1]; //Fortran style index to C/C++ style index
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) -
							1]; //Fortran style index to C/C++ style index
					}//for2
					reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

				}//for1

				//print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);

				//double tmp_conc = 0.0;
				//for (int i = 0; i < nkk; ++i)
				//	tmp_conc += concentration_data_pgt[i].back();
				//total_concentration_data_pgt.push_back(tmp_conc);
			}
			else if (tout < critical_time) {
				tout += dt1;
				//print out every * dt1
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0) {
						spe_drc_data_pgt[i].push_back(c_t[i]);
					}
					else {
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
					}
				}//for]

				//Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					reaction_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0;
						j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] -
							1]; //Fortran style index to C/C++ style index
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) -
							1]; //Fortran style index to C/C++ style index
					}//for2
					reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

				}//for1

				//print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);

				//double tmp_conc = 0.0;
				//for (int i = 0; i < nkk; ++i)
				//	tmp_conc += concentration_data_pgt[i].back();
				//total_concentration_data_pgt.push_back(tmp_conc);
			}
			//] print out

			//ckcppdlsodest_(&ti, &tout, &neq, xgst);
			//ckcppdlsodast_(&ti, &tout, &neq, xgst); /*currently used*/
			//update  molar concentration
			for (int i = 0; i < nkk; ++i)
				c_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];

			ti = tout;
			if ((end_time - ti) < 0.0001*(dt1 + dt2))
				ti = end_time;
			//tout+=dt;
		//} while (tout <= end_time);
		} while (time_data_pgt.back() < end_time);
		//while

		delete[] x_t;
		delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t;
		delete[] DDOT_t;
		delete[] FWDR_t;
		delete[] REVR_t;
	}

	void SOHRPropagator::time_initializer_based_on_guess_cv(rsp::my_time_t critical_time, rsp::my_time_t end_time)
	{
		//std::cout <<"deltaConcentration\t" <<this->pt.get<double>("SOHR_init.deltaConcentration") << std::endl;
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();
		//total_concentration_data_pgt.clear();

		concentration_data_pgt.clear();
		reaction_rate_data_pgt.clear();
		spe_drc_data_pgt.clear();

		//number of reaction in reaction network space
		std::size_t num_reaction = std::distance(reactionNetwork_chemkin_index_map.begin(),
			reactionNetwork_chemkin_index_map.end());

		// xgst[0]~xgst[nkk-1] are the MASS molar concentration of each species. xgst[nkk] is the temperature.
		//in which nkk is number of species
		//number of equations neq=nkk+1
		const int nkk = mechanism::kinetics::nkk(), neq =
			mechanism::kinetics::nkk() + 1, nii = mechanism::kinetics::nii();
		//initial conditions for lsode.
		double *xgst = new double[neq];
		for (int i = 0; i < neq; ++i) xgst[i] = 0.0;

		//molar fraction
		double *x_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; }

		//molar concentration.
		double *c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		//mass fraction
		double *y_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { y_t[i] = 0.0; }

		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLE FRACTION for species.
		double Temp;
		read_configuration(Temp, ckstore.pressure, neq - 1, x_t);
		//normalization
		double sum = 0.0;
		for (int i = 0; i < nkk; ++i) { sum += x_t[i]; }
		for (int i = 0; i < nkk; ++i) { x_t[i] /= sum; }

		// Convert mole fraction to molar concentration given pressure and temperature
		mechanism::kinetics::ckxtcp(&ckstore.pressure, &Temp, x_t, c_t);

		/*
		 * notice here that if initial concentration is zero, then the pseudo-first order rate coefficient
		 * or says propensity function might be zero, when invert the integral_of_K, might throw error, to prevent this from
		 * happening, we set the initial guesses to be a small number if they are zeros
		 */
		for (int i = 0; i < nkk; ++i)
		{
			if (c_t[i] == 0.0)
				c_t[i] = this->pgt_pt.get<double>("SOHR_init.deltaConcentration");
		}

		//convert molar concentration to mass fractions
		mechanism::kinetics::ckcty(c_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::ckrhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		mechanism::kinetics::ckytx(y_t, x_t);

		//add temperature to initial conditions for lsode.
		xgst[neq - 1] = Temp;
		//Read time step 'dt(s)'.
		//Read 'sohr' parameters.
		//initialize_sohr();

		//time step
		double dt1, dt2;

		if (critical_time == 0.0) {
			dt1 = dt2 = end_time / this->timeN2;
		}
		else if (critical_time == end_time) {
			//like 100 points, only 99 intervals
			dt1 = dt2 = end_time / (this->timeN1);
		}
		else {
			dt1 = critical_time / (this->timeN1);
			dt2 = (end_time - critical_time) / this->timeN2;
		}

		//initialize the 1st order ode.
		for (int i = 0; i < nkk; ++i)
			xgst[i] = c_t[i];

		double ti = 0.0, tout = 0.0;

		//species creation rates and destruction rates.
		double *CDOT_t = new double[nkk];
		double *DDOT_t = new double[nkk];
		//forward reaction rates, reverse reaction rates.
		double *FWDR_t = new double[nii];
		double *REVR_t = new double[nii];

		concentration_data_pgt.resize(nkk);
		reaction_rate_data_pgt.resize(num_reaction);
		double reaction_rate_tmp = 0.0;
		spe_drc_data_pgt.resize(nkk);

		do {
			//////////////////////////////////////////////////////////////////////////
			//Calculate useful stuff
			//////////////////////////////////////////////////////////////////////////
			////convert molar concentration to mass fractions
			//chemkincpp_sr::chemkin::ckcty(c_t, y_t);
			////Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			//chemkincpp_sr::chemkin::ckpy(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::ckcdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			////convert molar concentration to mole fraction
			//chemkincpp_sr::chemkin::ckctx(c_t, x_t);
			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::ckkfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			if (tout >= critical_time) {
				tout += dt2;
				//use the default dt to print out stuff
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0) {
						spe_drc_data_pgt[i].push_back(c_t[i]);
					}
					else {
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
					}
				}//for]

				//Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					reaction_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] -
							1]; //Fortran style index to C/C++ style index
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) -
							1]; //Fortran style index to C/C++ style index
					}//for2
					reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

				}//for1

				//print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);

				//double tmp_conc = 0.0;
				//for (int i = 0; i < nkk; ++i)
				//	tmp_conc += concentration_data_pgt[i].back();
				//total_concentration_data_pgt.push_back(tmp_conc);
			}
			else if (tout < critical_time) {
				tout += dt1;
				//print out every * dt1
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0) {
						spe_drc_data_pgt[i].push_back(c_t[i]);
					}
					else {
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
					}
				}//for]

				//Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					reaction_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0;
						j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] -
							1]; //Fortran style index to C/C++ style index
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) -
							1]; //Fortran style index to C/C++ style index
					}//for2
					reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

				}//for1

				//print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);

				//double tmp_conc = 0.0;
				//for (int i = 0; i < nkk; ++i)
				//	tmp_conc += concentration_data_pgt[i].back();
				//total_concentration_data_pgt.push_back(tmp_conc);
			}
			//] print out

			//ckcppdlsodest_(&ti, &tout, &neq, xgst);
			//ckcppdlsodast_(&ti, &tout, &neq, xgst); /*currently used*/
			//update  molar concentration
			for (int i = 0; i < nkk; ++i)
				c_t[i] = xgst[i];
			//update Temperature of system
			Temp = xgst[neq - 1];

			ti = tout;
			if ((end_time - ti) < 0.0001*(dt1 + dt2))
				ti = end_time;
			//tout+=dt;
		//} while (tout <= end_time);
		} while (time_data_pgt.back() < end_time);
		//while

		delete[] x_t;
		delete[] xgst;
		delete[] c_t;
		delete[] CDOT_t;
		delete[] DDOT_t;
		delete[] FWDR_t;
		delete[] REVR_t;
	}


	void SOHRPropagator::time_initializer_from_file(std::string filename_time, std::string filename_conc) {
		std::cout << filename_time << std::endl;
		std::cout << filename_conc << std::endl;
	}

	void SOHRPropagator::time_initializer_using_dlsode_cv(rsp::my_time_t critical_time, rsp::my_time_t end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();
		//total_concentration_data_pgt.clear();

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
			mechanism::kinetics::ckraex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::ckraex(&I_t, &R_A);
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
		mechanism::kinetics::ckxty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::ckrhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

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
		reaction_rate_data_pgt.resize(num_reaction); double reaction_rate_tmp = 0.0;
		spe_drc_data_pgt.resize(nkk);
		//spe_production_rate_data_pgt.resize(nkk);

		//while (tout < end_time)
		do
		{
			//////////////////////////////////////////////////////////////////////////
			//Calculate useful stuff
			//////////////////////////////////////////////////////////////////////////
			//convert mass fractions to molar fractions
			mechanism::kinetics::ckytx(y_t, x_t);
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::ckpy(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);
			//molar concentration
			mechanism::kinetics::ckytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s) 
			//and mass fractions
			mechanism::kinetics::ckcdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::ckkfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//////////////////////////////////////////////////////////////////////////
			//print out
			//////////////////////////////////////////////////////////////////////////
			//destruction relative rate Constant of species
			//[ print out
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN1 == 0)) || ((end_time - ti) < 0.001*dt)) {
				//use the default dt to print out stuff
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0)
					{
						spe_drc_data_pgt[i].push_back(c_t[i]);
						//spe_production_rate_data_pgt[i].push_back(CDOT_t[i]);
					}
					else
					{
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
						//spe_production_rate_data_pgt[i].push_back(CDOT_t[i]);
					}
				}//for]

				 //Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					reaction_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1]; //Fortran style index to C/C++ style index
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1]; //Fortran style index to C/C++ style index
					}//for2
					reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

				}//for1

				 //print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);

			}
			else if ((print_Count%lsodestore.deltaN2 == 0) || ((end_time - ti) < 0.001*dt)) {
				//print out every * dt		
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0)
					{
						spe_drc_data_pgt[i].push_back(c_t[i]);
						//spe_production_rate_data_pgt[i].push_back(CDOT_t[i]);
					}
					else
					{
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
						//spe_production_rate_data_pgt[i].push_back(CDOT_t[i]);
					}
				}//for]

				 //Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					reaction_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1]; //Fortran style index to C/C++ style index
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1]; //Fortran style index to C/C++ style index
					}//for2
					reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

				}//for1

				 //print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);

			}
			//] print out

			//ckcppdlsodev(&ti, &tout, &neq, xgst);
			ckcppdlsodav(&ti, &tout, &neq, xgst);
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

	void SOHRPropagator::time_initializer_using_dlsode_cv_ct(rsp::my_time_t critical_time, rsp::my_time_t end_time)
	{
		time_data_pgt.clear();
		temperature_data_pgt.clear();
		pressure_data_pgt.clear();
		//total_concentration_data_pgt.clear();

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
			mechanism::kinetics::ckraex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ sytle index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::ckraex(&I_t, &R_A);
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
		mechanism::kinetics::ckxty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::ckrhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

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
		reaction_rate_data_pgt.resize(num_reaction); double reaction_rate_tmp = 0.0;
		spe_drc_data_pgt.resize(nkk);
		//spe_production_rate_data_pgt.resize(nkk);

		//while (tout < end_time)
		do
		{
			//////////////////////////////////////////////////////////////////////////
			//Calculate useful stuff
			//////////////////////////////////////////////////////////////////////////
			//convert mass fractions to molar fractions
			mechanism::kinetics::ckytx(y_t, x_t);
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::ckpy(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);
			//molar concentration
			mechanism::kinetics::ckytcr(&ckstore.rhomass, &Temp, y_t, c_t);

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s) 
			//and mass fractions
			mechanism::kinetics::ckcdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::ckkfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//////////////////////////////////////////////////////////////////////////
			//print out
			//////////////////////////////////////////////////////////////////////////
			//destruction relative rate Constant of species
			//[ print out
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN1 == 0)) || ((end_time - ti) < 0.001*dt)) {
				//use the default dt to print out stuff
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0)
					{
						spe_drc_data_pgt[i].push_back(c_t[i]);
						//spe_production_rate_data_pgt[i].push_back(CDOT_t[i]);
					}
					else
					{
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
						//spe_production_rate_data_pgt[i].push_back(CDOT_t[i]);
					}
				}//for]

				 //Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					reaction_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1]; //Fortran style index to C/C++ style index
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1]; //Fortran style index to C/C++ style index
					}//for2
					reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

				}//for1

				 //print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);

			}
			else if ((print_Count%lsodestore.deltaN2 == 0) || ((end_time - ti) < 0.001*dt)) {
				//print out every * dt		
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0)
					{
						spe_drc_data_pgt[i].push_back(c_t[i]);
						//spe_production_rate_data_pgt[i].push_back(CDOT_t[i]);
					}
					else
					{
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
						//spe_production_rate_data_pgt[i].push_back(CDOT_t[i]);
					}
				}//for]

				 //Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					reaction_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1]; //Fortran style index to C/C++ style index
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1]; //Fortran style index to C/C++ style index
					}//for2
					reaction_rate_data_pgt[i].push_back(reaction_rate_tmp);

				}//for1

				 //print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);

			}
			//] print out

			//ckcppdlsodev(&ti, &tout, &neq, xgst);
			ckcppdlsodavt(&ti, &tout, &neq, xgst);
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

	void SOHRPropagator::time_initializer_using_dlsode_s_ct_np(rsp::my_time_t critical_time, rsp::my_time_t end_time)
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
			mechanism::kinetics::ckraex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ style index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::ckraex(&I_t, &R_A);
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

		//	//molar fraction
		//	double* x_t=new double[nkk]; double* y_t=new double[nkk];
		//	for (int i=0; i<nkk; ++i){ x_t[i]=0.0; y_t[i]=0.0;}

		////molar fraction
		//double* x_t = new double[nkk];
		//for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; }

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

		//	//convert mole fractions to mass fractions
		//	chemkincpp_sr::chemkin::ckxty(x_t, y_t);

		//convert molar concentration to molar fraction
		//chemkincpp_sr::chemkin::ckctx(c_t, x_t);


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
		reaction_rate_data_pgt.resize(num_reaction); double rxn_rate_tmp = 0.0;
		spe_drc_data_pgt.resize(nkk);

		//"imaginary pressure"
		//double i_pressure=0.0;

		//while (tout < end_time)
		do
		{
			//Returns the molar creation and destruction rates of the species given temperature(s) and molar concentration
			//The difference is how to treat the auto-catylytic reactions
			//chemkincpp_sr::chemkin::ckcdc(&Temp, c_t, CDOT_t, DDOT_t);
			this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			//chemkincpp_sr::chemkin::ckkfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
			// Applicable for reactions with rate constant independent of pressure
			// where sr stands for Shirong
			mechanism::kinetics::ckkfkrsr(&Temp, c_t, FWDR_t, REVR_t);

			//[ print out
			//int NPrecision=16;
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN1 == 0)) || ((end_time - ti) < 0.001*dt)) {
				//use the default dt to print out stuff
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0)
					{
						spe_drc_data_pgt[i].push_back(c_t[i]);
					}
					else
					{
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
					}
				}//for]

				 //Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					rxn_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							//Fortran style index to C/C++ style index
							rxn_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1];
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							//Fortran style index to C/C++ style index
							rxn_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1];
					}//for2
					reaction_rate_data_pgt[i].push_back(rxn_rate_tmp);

				}//for1

				 //print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);
			}
			else if ((print_Count%lsodestore.deltaN2 == 0) || ((end_time - ti) < 0.001*dt)) {
				//print out every * dt
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0)
					{
						spe_drc_data_pgt[i].push_back(c_t[i]);
					}
					else
					{
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
					}
				}//for]

				 //Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					rxn_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							//Fortran style index to C/C++ style index
							rxn_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1];
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							//Fortran style index to C/C++ style index
							rxn_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1];
					}//for2
					reaction_rate_data_pgt[i].push_back(rxn_rate_tmp);

				}//for1

				 //print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);
			}
			//] print out

			//ckcppdlsodest(&ti, &tout, &neq, xgst);
			ckcppdlsodast(&ti, &tout, &neq, xgst);
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


	void SOHRPropagator::time_initializer_using_dlsode_s_ct_np_cc1(rsp::my_time_t critical_time, rsp::my_time_t end_time)
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
			mechanism::kinetics::ckraex(&I_t, &R_A);
			//cout<<I_t<<" "<<R_A<<"\t";
			//C/C++ style index to Fortran style index
			I_t = -I_t; R_A *= uncertainties[abs(I_t) - 1];
			mechanism::kinetics::ckraex(&I_t, &R_A);
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

		//	//molar fraction
		//	double* x_t=new double[nkk]; double* y_t=new double[nkk];
		//	for (int i=0; i<nkk; ++i){ x_t[i]=0.0; y_t[i]=0.0;}

		////molar fraction
		//double* x_t = new double[nkk];
		//for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; }

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

		//	//convert mole fractions to mass fractions
		//	chemkincpp_sr::chemkin::ckxty(x_t, y_t);

		//convert molar concentration to molar fraction
		//chemkincpp_sr::chemkin::ckctx(c_t, x_t);


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
		reaction_rate_data_pgt.resize(num_reaction); double rxn_rate_tmp = 0.0;
		spe_drc_data_pgt.resize(nkk);

		//"imaginary pressure"
		//double i_pressure=0.0;

		//while (tout < end_time)
		do
		{
			//Returns the molar creation and destruction rates of the species given temperature(s) and molar concentration
			//The difference is how to treat the auto-catylytic reactions
			//chemkincpp_sr::chemkin::ckcdc(&Temp, c_t, CDOT_t, DDOT_t);
			this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);

			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			//chemkincpp_sr::chemkin::ckkfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
			// Applicable for reactions with rate constant independent of pressure
			// where sr stands for Shirong
			mechanism::kinetics::ckkfkrsr(&Temp, c_t, FWDR_t, REVR_t);

			//[ print out
			//int NPrecision=16;
			if (((tout >= critical_time) && (print_Count%lsodestore.deltaN1 == 0)) || ((end_time - ti) < 0.001*dt)) {
				//use the default dt to print out stuff
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0)
					{
						spe_drc_data_pgt[i].push_back(c_t[i]);
					}
					else
					{
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
					}
				}//for]

				 //Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					rxn_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							//Fortran style index to C/C++ style index
							rxn_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1];
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							//Fortran style index to C/C++ style index
							rxn_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1];
					}//for2
					reaction_rate_data_pgt[i].push_back(rxn_rate_tmp);

				}//for1

				 //print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);
			}
			else if ((print_Count%lsodestore.deltaN2 == 0) || ((end_time - ti) < 0.001*dt)) {
				//print out every * dt
				time_data_pgt.push_back(ti);

				//destruction rate const
				for (int i = 0; i < nkk; ++i)//[for
				{
					if (c_t[i] <= 0.0)
					{
						spe_drc_data_pgt[i].push_back(c_t[i]);
					}
					else
					{
						//just need the destruction rate const of species
						spe_drc_data_pgt[i].push_back(DDOT_t[i] / c_t[i]);
					}
				}//for]

				 //Rates of each reactions
				for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
					rxn_rate_tmp = 0.0;
					//check reactionNetwork_chemkin_index_map
					for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
						if (reactionNetwork_chemkin_index_map[i][j] > 0)
							//Fortran style index to C/C++ style index
							rxn_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] - 1];
						else if (reactionNetwork_chemkin_index_map[i][j] < 0)
							//Fortran style index to C/C++ style index
							rxn_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) - 1];
					}//for2
					reaction_rate_data_pgt[i].push_back(rxn_rate_tmp);

				}//for1

				 //print concentration and temperature
				for (int i = 0; i < nkk; ++i)
					concentration_data_pgt[i].push_back(c_t[i]);

				temperature_data_pgt.push_back(xgst[neq - 1]);
				pressure_data_pgt.push_back(ckstore.pressure);
			}
			//] print out

			//ckcppdlsodest(&ti, &tout, &neq, xgst);
			ckcppdlsodastcc1(&ti, &tout, &neq, xgst);
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

	void SOHRPropagator::update_spe_concentration_at_time(std::size_t time_index, std::size_t spe_index,
		double delta_conc) {
		this->concentration_data_pgt[spe_index][time_index] += delta_conc*this->P2C[spe_index];
	}

	void SOHRPropagator::update_spe_concentration_at_time_range(std::size_t init_time_index,
		std::size_t end_time_index, std::size_t spe_index,
		double delta_conc) {
		for (std::size_t i = init_time_index + 1; i < end_time_index; ++i) {
			this->concentration_data_pgt[spe_index][i] += delta_conc*this->P2C[spe_index];
		}
	}


	void SOHRPropagator::update_spe_drc_based_on_spe_concentration_s_ct_np() {
		//in which nkk is number of species
		const int nkk = mechanism::kinetics::nkk();

		//molar concentration.
		double *c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		//species creation rates and destruction rates.
		double *CDOT_t = new double[nkk];
		double *DDOT_t = new double[nkk];

		//temperature, dummy variable in Lotka-Volterra model
		double Temp = 298.0;

		for (std::size_t k = 0; k < this->time_data_pgt.size(); ++k) {
			//Calculate drc
			for (int i = 0; i < nkk; ++i) {
				c_t[i] = concentration_data_pgt[i][k];//ith species, kth point
			}
			//The difference is how to treat the auto-catylytic reactions
			//chemkincpp_sr::chemkin::ckcdc(&Temp, c_t, CDOT_t, DDOT_t);
			this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);

			//std::copy(DDOT_t, DDOT_t+nkk, std::ostream_iterator<double>(std::cout, "\t")); std::cout<<std::endl;

			//update drc at every time point
			//destruction relative rate Constant of species
			for (int i = 0; i < nkk; ++i)//[for
			{
				if (c_t[i] <= 0.0) {
					spe_drc_data_pgt[i][k] = c_t[i];
				}
				else {
					//just need the destruction rate const of species
					//ith species, kth point
					spe_drc_data_pgt[i][k] = DDOT_t[i] / c_t[i];
				}
			}//for]

		};

		delete[] c_t;
		delete[] CDOT_t;
		delete[] DDOT_t;
	}

	void SOHRPropagator::update_reaction_rate_based_on_spe_concentration_s_ct_np() {
		//in which nkk is number of species
		const int nkk = mechanism::kinetics::nkk(), nii = mechanism::kinetics::nii();

		////molar fraction
		//double *x_t = new double[nkk];
		//for (int i = 0; i < nkk; ++i) { x_t[i] = 0.0; }

		//molar concentration.
		double *c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = 0.0; }

		//temperature, dummy variable in Lotka-Volterra model
		double Temp = temperature_data_pgt[0];

		//forward reaction rates, reverse reaction rates.
		double *FWDR_t = new double[nii];
		double *REVR_t = new double[nii];

		double reaction_rate_tmp = 0.0;
		for (std::size_t k = 0; k < time_data_pgt.size(); ++k) {
			//Calculate reaction rates
			for (int i = 0; i < nkk; ++i) {
				c_t[i] = concentration_data_pgt[i][k];//ith species, kth point
			}
			////convert molar concentration to molar fractions
			//chemkincpp_sr::chemkin::ckctx(c_t, x_t);

			this->cal_reaction_rate(&Temp, c_t, FWDR_t, REVR_t);

			//update reaction rate at every time point
			//reaction rate of reactions
			for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
				reaction_rate_tmp = 0.0;
				//check reactionNetwork_chemkin_index_map
				for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
					if (reactionNetwork_chemkin_index_map[i][j] > 0)
						reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] -
						1]; //Fortran style index to C/C++ style index
					else if (reactionNetwork_chemkin_index_map[i][j] < 0)
						reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) -
						1]; //Fortran style index to C/C++ style index
				}//for2
				reaction_rate_data_pgt[i][k] = reaction_rate_tmp;

			}//for1

		};

		//delete[] x_t;
		delete[] c_t;
		delete[] FWDR_t;
		delete[] REVR_t;
	}

	void SOHRPropagator::update_pressure_based_on_spe_concentration_cv()
	{
		//in which nkk is number of species
		const int nkk = mechanism::kinetics::nkk();

		//mass fractions
		double* y_t = new double[nkk];
		//molar concentration.
		double* c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = concentration_data_pgt[i][0]; }

		double Temp = temperature_data_pgt[0];
		//temperature derivative with respect to time
		ckstore.pressure = pressure_data_pgt[0];

		//convert molar concentration to mass fractions
		mechanism::kinetics::ckcty(c_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::ckrhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//update pressure, not temperature
		for (std::size_t k = 0; k < this->time_data_pgt.size(); ++k) {
			//retrieve molar concentration
			for (int i = 0; i < nkk; ++i) {
				c_t[i] = concentration_data_pgt[i][k];//ith species, kth point
			}

			//convert molar concentration to mass fractions
			mechanism::kinetics::ckcty(c_t, y_t);

			//get temperature
			Temp = temperature_data_pgt[k];
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::ckpy(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);

			//update pressures
			pressure_data_pgt[k] = ckstore.pressure;
		}
	}

	void SOHRPropagator::update_temperature_pressure_based_on_spe_concentration_cv()
	{
		//in which nkk is number of species
		const int nkk = mechanism::kinetics::nkk();

		//mass fractions
		double* y_t = new double[nkk];
		//molar concentration.
		double* c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = concentration_data_pgt[i][0]; }

		double Temp = temperature_data_pgt[0];
		//temperature derivative with respect to time
		double Temp_dot = 0.0;
		ckstore.pressure = pressure_data_pgt[0];

		//convert molar concentration to mass fractions
		mechanism::kinetics::ckcty(c_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::ckrhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		//update temperature and pressure
		for (std::size_t k = 1; k < this->time_data_pgt.size(); ++k) {
			//retrieve molar concentration
			for (int i = 0; i < nkk; ++i) {
				c_t[i] = concentration_data_pgt[i][k];//ith species, kth point
			}
			//previous temperature
			Temp = temperature_data_pgt[k - 1];
			//convert molar concentration to mass fractions
			mechanism::kinetics::ckcty(c_t, y_t);
			mechanism::kinetics::calculate_t_dot_cv(y_t, &Temp, &Temp_dot);

			//new temperature
			Temp += Temp_dot*(time_data_pgt[k] - time_data_pgt[k - 1]);
			temperature_data_pgt[k] = Temp;

			//convert molar concentration to mass fractions
			mechanism::kinetics::ckcty(c_t, y_t);
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::ckpy(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);

			//update pressures
			pressure_data_pgt[k] = ckstore.pressure;
		}
	}

	//x_t molar concentration
	//y_t mass fraction
	//c_t molar concentration
	void SOHRPropagator::update_spe_drc_reaction_rate_based_on_spe_concentration_cv()
	{
		//in which nkk is number of species
		const int nkk = mechanism::kinetics::nkk(), nii = mechanism::kinetics::nii();

		//molar fraction.
		double *x_t = new double[nkk];
		//mass fractions
		double* y_t = new double[nkk];
		//molar concentration.
		double* c_t = new double[nkk];
		for (int i = 0; i < nkk; ++i) { c_t[i] = concentration_data_pgt[i][0]; }

		//species creation rates and destruction rates.
		double *CDOT_t = new double[nkk];
		double *DDOT_t = new double[nkk];

		//forward reaction rates, reverse reaction rates.
		double *FWDR_t = new double[nii];
		double *REVR_t = new double[nii];

		//temperature, dummy variable in Lotka-Volterra model
		double Temp = temperature_data_pgt[0];

		////The pressure has to be in cgs unit, convert pascal(Pa) to dyne/square centimeter [dyn/cm**2]
		//// 1atm= 101325 Pa
		//// 1Pa= 10 dyn/cm**2
		//// P = nRT
		//double ru, ruc, pa;
		//chemkincpp_sr::chemkin::ckrp(&ru, &ruc, &pa);
		////ckstore.pressure = pressure_data_pgt[0];
		//ckstore.pressure = std::accumulate(c_t, c_t+nkk, 0.0)*ru*Temp;
		ckstore.pressure = pressure_data_pgt[0];

		//convert molar concentration to mass fractions
		mechanism::kinetics::ckcty(c_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore.rhomass' is constant in the simulation.
		mechanism::kinetics::ckrhoy(&ckstore.pressure, &Temp, y_t, &ckstore.rhomass);

		double reaction_rate_tmp = 0.0;
		for (std::size_t k = 0; k < this->time_data_pgt.size(); ++k) {
			//supposed to use different temperature and pressure every time point
			Temp = temperature_data_pgt[k];
			ckstore.pressure = pressure_data_pgt[k];
			//Calculate drc
			//retrieve molar concentration
			for (int i = 0; i < nkk; ++i) {
				c_t[i] = concentration_data_pgt[i][k];//ith species, kth point
			}

			//convert molar concentration to mass fractions
			mechanism::kinetics::ckcty(c_t, y_t);
			//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
			mechanism::kinetics::ckpy(&ckstore.rhomass, &Temp, y_t, &ckstore.pressure);

			////update pressures
			//pressure_data_pgt[k] = ckstore.pressure;

			//Returns the molar creation and destruction rates of the species given mass density, temperature(s)
			//and mass fractions
			mechanism::kinetics::ckcdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);

			//convert molar concentration to mole fraction
			mechanism::kinetics::ckctx(c_t, x_t);
			//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.
			mechanism::kinetics::ckkfkr(&ckstore.pressure, &Temp, x_t, FWDR_t, REVR_t);

			//update drc at every time point
			//destruction relative rate Constant of species
			for (int i = 0; i < nkk; ++i)//[for
			{
				if (c_t[i] <= 0.0) {
					spe_drc_data_pgt[i][k] = c_t[i];
				}
				else {
					//just need the destruction rate const of species
					//ith species, kth point
					spe_drc_data_pgt[i][k] = DDOT_t[i] / c_t[i];
				}
			}//for]

			//update reaction rate at every time point
			//reaction rate of reactions
			for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i) {//for1
				reaction_rate_tmp = 0.0;
				//check reactionNetwork_chemkin_index_map
				for (std::size_t j = 0; j < reactionNetwork_chemkin_index_map[i].size(); ++j) {//for2
					if (reactionNetwork_chemkin_index_map[i][j] > 0)
						reaction_rate_tmp += FWDR_t[reactionNetwork_chemkin_index_map[i][j] -
						1]; //Fortran style index to C/C++ style index
					else if (reactionNetwork_chemkin_index_map[i][j] < 0)
						reaction_rate_tmp += REVR_t[abs(reactionNetwork_chemkin_index_map[i][j]) -
						1]; //Fortran style index to C/C++ style index
				}//for2
				reaction_rate_data_pgt[i][k] = reaction_rate_tmp;

			}//for1

		};

		delete[] x_t;
		delete[] y_t;
		delete[] c_t;

		delete[] CDOT_t;
		delete[] DDOT_t;

		delete[] FWDR_t;
		delete[] REVR_t;
	}


	void SOHRPropagator::set_concentration_data_zero() {
		for (std::size_t i = 0; i < concentration_data_pgt.size(); ++i) {
			std::fill(concentration_data_pgt[i].begin(), concentration_data_pgt[i].end(), 0.0);
		}
	}


	void SOHRPropagator::rescale_concentration_data(double factor) {
		if (factor != 1.0) {
			for (std::size_t i = 0; i < concentration_data_pgt.size(); ++i) {
				for (std::size_t j = 0; j < concentration_data_pgt[i].size(); ++j) {
					concentration_data_pgt[i][j] *= factor;
				}
			}
		}//if
	}

	void SOHRPropagator::rescale_concentration_data(const std::vector<double> &factor_v)
	{
		assert(concentration_data_pgt.size() == factor_v.size());
		for (std::size_t i = 0; i < concentration_data_pgt.size(); ++i) {
			if (factor_v[i] != 1.0) {
				for (std::size_t j = 0; j < concentration_data_pgt[i].size(); ++j) {
					concentration_data_pgt[i][j] *= factor_v[i];
				}
			}
		}//if
	}

	void SOHRPropagator::rescale_prob_matrix_data(std::vector<std::vector<double>>& prob_Mat, const double factor)
	{
		if (factor != 1.0) {
			for (std::size_t i = 0; i < prob_Mat.size(); ++i) {
				for (std::size_t j = 0; j < prob_Mat[i].size(); ++j) {
					prob_Mat[i][j] *= factor;
				}
			}
		}//if
	}

	void SOHRPropagator::rescale_prob_matrix_data(std::vector<std::vector<double>>& prob_Mat, const std::vector<double>& factor_v)
	{
		assert(prob_Mat.size() == factor_v.size());
		for (std::size_t i = 0; i < prob_Mat.size(); ++i) {
			if (factor_v[i] != 1.0) {
				for (std::size_t j = 0; j < prob_Mat[i].size(); ++j) {
					prob_Mat[i][j] *= factor_v[i];
				}
			}
		}//if
	}

	void SOHRPropagator::normalize_prob_matrix_data(std::vector<std::vector<double>> & prob_Mat)
	{
		double total = 0.0;
		for (std::size_t i = 0; i < prob_Mat[0].size(); ++i) {
			total = 0.0;
			for (std::size_t j = 0; j < prob_Mat.size(); ++j) {
				total += prob_Mat[j][i];
			}
			for (std::size_t j = 0; j < prob_Mat.size(); ++j) {
				prob_Mat[j][i] /= total;
			}
		}
	}

	void SOHRPropagator::update_concentration_data_from_spe_based_probability_matrix(const std::vector<std::vector<double>> &prob_Mat) {
		std::copy(prob_Mat.begin(), prob_Mat.end(), this->concentration_data_pgt.begin());
	}

	void SOHRPropagator::initialize_sohr() {
		this->timeN1 = this->pgt_pt.get<int>("SOHR_init.timeN1");
		this->timeN2 = this->pgt_pt.get<int>("SOHR_init.timeN2");
		this->number2Concentration = this->pgt_pt.get<double>("SOHR_init.massConservationFactor") / this->pgt_pt.get<int>("pathway.trajectoryNumber");
		this->negative_number2Concentration = (-1.0) * this->number2Concentration;
		this->trajectoryNumber = this->pgt_pt.get<int>("pathway.trajectoryNumber");
		this->iterationNumber = this->pgt_pt.get<int>("SOHR_init.iterationNumber");
		for (auto key1 : pgt_pt.get_child("SOHR_init.P2C")) {
			P2C.push_back(key1.second.get_value<double>());
		}
	}

	void SOHRPropagator::propagate_pgt() {
		time_propagator_s_ct_np_s2m_pgt(uncertainties, this->pgt_pt.get<double>("time.critical_time"),
			this->pgt_pt.get<double>("time.max_time"));
	}

	void SOHRPropagator::time_propagator_s_ct_np_s2m_pgt(std::vector<double> uncertainties, double critical_time,
		double end_time) {

	}

	void SOHRPropagator::time_propagator_s_ct_np_cc1_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}

	void SOHRPropagator::time_propagator_s_ct_np_cc2_s2m_pgt(std::vector<double> uncertainties, double critical_time, double end_time)
	{
	}


}


#endif
