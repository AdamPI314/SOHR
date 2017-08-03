#ifndef __SUPERPROPAGATOR_CPP_
#define __SUPERPROPAGATOR_CPP_

#define PRINT_PRECISION 15

#include <numeric>
#include <limits>

#include <boost/property_tree/json_parser.hpp> //for json_reader
#include <boost/foreach.hpp> //for BOOST_FOREACH
#include <boost/optional/optional.hpp> //for optional

#include "../../tools/misc/fortran_routine_block_alias.h"
#include "../../tools/misc/global_extern_vars.h"
#include "../../chemkinCpp/chemkincpp.h"

//put it in the end, because it has "cubicSpline/interp_1d.h", which has to be after all boost libraries
#include "superPropagator.h"

namespace propagator_sr {
	superPropagator::superPropagator(std::vector<double> uncertainties_in, std::string cwd_in)
	{
		//current working directory
		this->cwd_pgt = cwd_in;

		std::copy(uncertainties_in.begin(), uncertainties_in.end(), std::back_inserter(this->uncertainties));

		//defer this to actually use time
		initiate_propagator();

	}//superPropagator


	superPropagator::~superPropagator()
	{
		for (std::size_t i = 0; i < spe_drc_pgt.size(); ++i) {
			if (spe_drc_pgt[i] != 0)
				delete spe_drc_pgt[i];
		}
		for (std::size_t i = 0; i < spe_drc_int_pgt.size(); ++i) {
			if (spe_drc_int_pgt[i] != 0)
				delete spe_drc_int_pgt[i];
		}
		for (std::size_t i = 0; i < spe_drc_int_time_pgt.size(); ++i) {
			if (spe_drc_int_time_pgt[i] != 0)
				delete spe_drc_int_time_pgt[i];
		}
		for (std::size_t i = 0; i < reaction_rate_pgt.size(); ++i) {
			if (reaction_rate_pgt[i] != 0)
				delete reaction_rate_pgt[i];
		}
		for (std::size_t i = 0; i < concentration_pgt.size(); ++i) {
			if (concentration_pgt[i] != 0)
				delete concentration_pgt[i];
		}
		//if (total_concentration_pgt != 0)
		//	delete total_concentration_pgt;
		if (time_temperature_pgt != 0)
			delete time_temperature_pgt;
		if (time_pressure_pgt != 0)
			delete time_pressure_pgt;
		if (temperature_time_pgt != 0)
			delete temperature_time_pgt;
	}

	void superPropagator::initialize_cubic_spline_pointer()
	{
		std::fill(spe_drc_pgt.begin(), spe_drc_pgt.end(), (Linear_interp*)(0));
		std::fill(spe_drc_int_pgt.begin(), spe_drc_int_pgt.end(), (Linear_interp*)(0));
		std::fill(spe_drc_int_time_pgt.begin(), spe_drc_int_time_pgt.end(), (Linear_interp*)(0));
		std::fill(reaction_rate_pgt.begin(), reaction_rate_pgt.end(), (Linear_interp*)(0));
		std::fill(concentration_pgt.begin(), concentration_pgt.end(), (Linear_interp*)(0));

		//total_concentration_pgt = (Linear_interp*)0;
		time_temperature_pgt = (Linear_interp*)0;
		time_pressure_pgt = (Linear_interp*)0;
		temperature_time_pgt = (Linear_interp*)0;
	}

	void superPropagator::set_fast_reaction_pgt()
	{
		std::vector<std::size_t> fast_reaction_index;
		//read with json_parser as property_tree
		//nice and easy
		BOOST_FOREACH(boost::property_tree::ptree::value_type const& key1, this->pgt_pt.get_child("pathway.fast_reaction")) {
			fast_reaction_index.push_back(boost::lexical_cast<std::size_t>(key1.first));
			fast_reaction_index.push_back(key1.second.get_value<std::size_t>());
		}
		//	std::copy(fast_reaction_index.begin(), fast_reaction_index.end(), std::ostream_iterator<std::size_t>(std::cout, "\t"));

		this->fast_reaction_pgt = fast_reaction_index;
	}

	void superPropagator::set_fast_reaction_rate_to_zero_pgt()
	{
		for (std::size_t i = 0; i < this->fast_reaction_pgt.size(); ++i) {
			std::fill(reaction_rate_data_pgt[fast_reaction_pgt[i]].begin(), reaction_rate_data_pgt[fast_reaction_pgt[i]].end(), 0.0);
		}
	}

	rsp::temperature_t superPropagator::return_target_temperature() const
	{
		return this->pgt_pt.get<rsp::temperature_t>("T.target_temperature");
	}

	void superPropagator::print_temperature_targeted_concentration_pgt(std::string outfile) const
	{
		std::ofstream out_file((this->cwd_pgt + outfile).c_str());
		double time = evaluate_time_at_temperature(return_target_temperature());
		//double time = time_data_pgt.back();
		for (size_t i = 0; i < concentration_pgt.size(); ++i) {
			out_file << evaluate_concentration_at_time(time, i);
			if (i != concentration_pgt.size() - 1)
				out_file << " ";
		}

		out_file.close();
	}

	void superPropagator::print_pgt()
	{
		cout << "Print:" << endl;

		std::vector<int> duplicate_reaction_ind;
		rsp::relationshipParser::read_chem_out_duplicated_reaction(duplicate_reaction_ind, "./input/chem.out");
		cout << "duplicated reaction:(C/C++ style index) " << endl;
		std::copy(duplicate_reaction_ind.begin(), duplicate_reaction_ind.end(), std::ostream_iterator<int>(cout, " "));
		cout << endl;

		//print network reaction index and ChemKin style reaction index
		for (rsp::reactionNetwork_chemkin_index_map_t::iterator itr = reactionNetwork_chemkin_index_map.begin(); itr != reactionNetwork_chemkin_index_map.end(); ++itr) {//for
			std::cout << itr->first << "\t:";
			std::copy((*itr).second.begin(), (*itr).second.end(), std::ostream_iterator<int>(std::cout, "\t"));
			std::cout << "\n";
		}//for

		std::cout << "number of reactions:\t" << this->reaction_rate_data_pgt.size() << std::endl;


		cout << "Final temperature:\t" << temperature_data_pgt.back() << endl;

		int MM_t, KK_t, II_t, NFIT_t;
		chemkincpp_sr::chemkin::ckindx(&MM_t, &KK_t, &II_t, &NFIT_t);
		cout << "mechanism total element count: " << MM_t << endl;
		cout << "mechanism total species count: " << KK_t << endl;
		cout << "mechanism total reaction count: " << II_t << endl;
		cout << "number of coefficients in fits to thermodynamic data for a temperature range;\nNFIT=number of coefficients in polynomial fits to CP/R plus 2: " << NFIT_t << endl << endl;

		const int neq = chemkincpp_sr::chemkin::nkk() + 1, nkk = chemkincpp_sr::chemkin::nkk(), nii = chemkincpp_sr::chemkin::nii();

		//molar fraction
		double* x_t = new double[neq - 1]; double* y_t = new double[neq - 1];
		for (int i = 0; i < neq - 1; ++i) { x_t[i] = 0.0; y_t[i] = 0.0; }
		// Read 'Temperature(K)' and 'Pressure(atm)'.
		// Read	MOLAR FRACTION for species.
		double Temp, Pressure;
		read_configuration(Temp, Pressure, neq - 1, x_t);
		//	chem_init(Temp, Pressure, neq-1, x_t, this->cwd_pgt+std::string("/input/setting.cfg"));


		//convert mole fractions to mass fractions
		chemkincpp_sr::chemkin::ckxty(x_t, y_t);
		//Returns the mass density of the gas mixture given pressure, temperature(s) and mass fractions
		//The simulation is done in the constant volume condition. So 'ckstore_.rhomass' is constant in the simulation.
		chemkincpp_sr::chemkin::ckrhoy(&Pressure, &Temp, y_t, &ckstore.rhomass);

		//Returns the pressure of the gas mixture given mass density, temperature(s) and mass fractions.
		chemkincpp_sr::chemkin::ckpy(&ckstore.rhomass, &Temp, y_t, &Pressure);

		//molar concentration.
		double* c_t = new double[nkk];
		chemkincpp_sr::chemkin::ckytcr(&ckstore.rhomass, &Temp, y_t, c_t);
		std::cout << "molar concentration:" << std::endl;
		std::copy(c_t, c_t + nkk, std::ostream_iterator<double>(cout, " ")); cout << endl;

		//species creation rates and destruction rates.
		double* CDOT_t = new double[nkk]; double* DDOT_t = new double[nkk];
		//Returns the molar creation and destruction rates of the species given mass density, temperature(s) and mass fractions
		chemkincpp_sr::chemkin::ckcdyr(&ckstore.rhomass, &Temp, y_t, CDOT_t, DDOT_t);
		std::cout << "molar creation and destruction rates, subroutine ckcdyr:" << std::endl;
		std::copy(CDOT_t, CDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(DDOT_t, DDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;

		std::fill(CDOT_t, CDOT_t + nkk, 0.0);
		std::fill(DDOT_t, DDOT_t + nkk, 0.0);
		std::copy(CDOT_t, CDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(DDOT_t, DDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;
		//The difference is how to treat the auto-catylytic reactions
		chemkincpp_sr::chemkin::ckcdc(&Temp, c_t, CDOT_t, DDOT_t);
		//this->cal_spe_destruction_rate(&Temp, c_t, CDOT_t, DDOT_t);
		std::cout << "molar creation and destruction rates, subroutine ckcdc:" << std::endl;
		std::copy(CDOT_t, CDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(DDOT_t, DDOT_t + neq, std::ostream_iterator<double>(cout, " ")); cout << endl;

		//forward reaction rates, reverse reaction rates.
		double* FWDR_t = new double[nii]; double* REVR_t = new double[nii];
		//Returns the forward and reverse reaction rates for reactions given pressure, temperature(s) and mole fractions.

		std::cout << "Pressure:\t" << Pressure << std::endl;
		double init = 0.0;
		std::cout << "molar concentration:\t" << std::accumulate(c_t, c_t + nkk, init) << std::endl;
		chemkincpp_sr::chemkin::ckkfkr(&Pressure, &Temp, x_t, FWDR_t, REVR_t);
		std::cout << "forward and reverse reaction rates:" << std::endl;
		std::copy(FWDR_t, FWDR_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(REVR_t, REVR_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;

		//forward reaction rates constant, reverse reaction rates constant.
		double* RKFT_t = new double[nii]; double* RKRT_t = new double[nii];
		//Returns the forward and reverse reaction rates Constant for reactions given pressure and temperature(s).
		chemkincpp_sr::chemkin::ckkfrt(&Pressure, &Temp, RKFT_t, RKRT_t);
		std::cout << "forward and reverse reaction rate constant:" << std::endl;
		std::copy(RKFT_t, RKFT_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(RKRT_t, RKRT_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;

		cout << endl << RKFT_t[5 - 1] << " " << FWDR_t[5 - 1] << " " << c_t[3 - 1] << " " << RKFT_t[5 - 1] * c_t[3 - 1] * (c_t[3 - 1] * 2.5 + c_t[1 - 1]) << endl;
		cout << endl << RKRT_t[14 - 1] << " " << REVR_t[14 - 1] << " " << c_t[3 - 1] << " " << RKRT_t[14 - 1] * c_t[1 - 1] * c_t[3 - 1]/**(c_t[3-1]+c_t[1-1])*/ << endl;

		cout << "nii: " << nii << endl;
		int I_t = 25;	double R_A = 0.0;
		chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
		cout << R_A << endl;
		I_t = -25; R_A *= 5.0;
		chemkincpp_sr::chemkin::ckraex(&I_t, &R_A);
		cout << R_A << endl;
		//Returns the forward and reverse reaction rates Constant for reactions given pressure and temperature(s).
		chemkincpp_sr::chemkin::ckkfrt(&Pressure, &Temp, RKFT_t, RKRT_t);
		std::copy(RKFT_t, RKFT_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		std::copy(RKRT_t, RKRT_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;

		//print pre-coefficient, activation energy, Temperature factor
		double *RA2_t = new double[nii];
		double *RB2_t = new double[nii];
		double *RE2_t = new double[nii];

		chemkincpp_sr::chemkin::ckabe(RA2_t, RB2_t, RE2_t);
		cout << "RA: " << endl;
		std::copy(RA2_t, RA2_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		cout << "RB: " << endl;
		std::copy(RB2_t, RB2_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;
		cout << "RE: " << endl;
		std::copy(RE2_t, RE2_t + nii, std::ostream_iterator<double>(cout, " ")); cout << endl;

		delete[] RA2_t; delete[] RB2_t; delete RE2_t;
		delete[] x_t; delete[] y_t;
		delete[] c_t;
		delete[] CDOT_t; delete[] DDOT_t;
		delete[] FWDR_t; delete[] REVR_t;
		delete[] RKFT_t; delete[] RKRT_t;

	}

	void superPropagator::w2f_pgt(std::string tag)
	{
		// time
		std::ofstream fout((this->cwd_pgt + std::string("/output/time_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < time_data_pgt.size(); ++i) {
			fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << time_data_pgt[i] << std::endl;
		}
		fout.clear(); fout.close();

		// temperature
		fout.open((this->cwd_pgt + std::string("/output/temperature_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < temperature_data_pgt.size(); ++i) {
			fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << temperature_data_pgt[i] << std::endl;
		}
		fout.clear(); fout.close();

		// pressure
		fout.open((this->cwd_pgt + std::string("/output/pressure_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < pressure_data_pgt.size(); ++i) {
			//fout<<pressure_data_pgt[i]/1013250<<"\t"<<pressure_data_pgt[i]<<std::endl;
			fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << pressure_data_pgt[i] / 1013250 << "\t" << pressure_data_pgt[i] / 1000000 << "\t" << pressure_data_pgt[i] << std::endl;
		}
		fout.clear(); fout.close();

		// concentration, 8 species
		fout.open((this->cwd_pgt + std::string("/output/concentration_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < concentration_data_pgt[0].size(); ++i) {
			for (size_t j = 0; j < concentration_data_pgt.size(); ++j) {
				fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << concentration_data_pgt[j][i];
				if (j < (concentration_data_pgt.size() - 1))
					fout << ",";
			}
			fout << endl;
		}
		fout.clear(); fout.close();

		//// int_drc, cumulative destructive rate constant
		//fout.open((this->cwd_pgt + std::string("/output/int_drc_") + tag + std::string(".csv")).c_str());
		//for (size_t i = 0; i < spe_drc_int_data_pgt[0].size(); ++i) {
			//for (size_t j = 0; j < spe_drc_int_data_pgt.size(); ++j) {
				//fout << std::setprecision(std::numeric_limits<double>::max_digits10+1) << spe_drc_int_data_pgt[j][i];
				//if (j < spe_drc_int_data_pgt.size() - 1)
					//fout << ",";
			//}
			//fout << endl;
		//}
		//fout.clear(); fout.close();

		// destruction rate constant,pseudo first order destructive rate constant
		fout.open((this->cwd_pgt + std::string("/output/drc_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < spe_drc_data_pgt[0].size(); ++i) {
			for (size_t j = 0; j < spe_drc_data_pgt.size(); ++j) {
				fout << std::setprecision(std::numeric_limits<double>::max_digits10+1) << spe_drc_data_pgt[j][i];
				if (j < spe_drc_data_pgt.size() - 1)
					fout << ",";
			}
			fout << endl;
		}
		fout.clear(); fout.close();

		////species production rate
		//fout.open((this->cwd_pgt + std::string("/output/spe_production_rate_") + tag + std::string(".csv")).c_str());
		//std::cout << spe_production_rate_data_pgt.size() << std::endl;
		//std::cout << spe_production_rate_data_pgt[0].size() << std::endl;
		//for (size_t i = 0; i < spe_production_rate_data_pgt[0].size(); ++i) {
		//	for (size_t j = 0; j < spe_production_rate_data_pgt.size(); ++j) {
		//		fout << std::setprecision(std::numeric_limits<double>::max_digits10+1) << spe_production_rate_data_pgt[j][i];
		//		if (j < spe_production_rate_data_pgt.size() - 1)
		//			fout << ",";
		//	}
		//	fout << endl;
		//}
		//fout.clear(); fout.close();

		// reaction rate
		fout.open((this->cwd_pgt + std::string("/output/reaction_rate_") + tag + std::string(".csv")).c_str());
		for (size_t i = 0; i < reaction_rate_data_pgt[0].size(); ++i) {
			for (size_t j = 0; j < reaction_rate_data_pgt.size(); ++j) {
				fout << std::setprecision(std::numeric_limits<double>::max_digits10+1) << reaction_rate_data_pgt[j][i];
				if (j < reaction_rate_data_pgt.size() - 1)
					fout << ",";
			}
			fout << std::endl;
		}
		fout.clear(); fout.close();

		//	fout.open((this->cwd_pgt+std::string("/output/rate_const_pre_factor_")+tag+std::string(".csv")).c_str());
		//	int I_t=1;	double R_A=0.0;
		//	//cout<<"\nPre-exponential Constant:"<<endl;
		//	for(;static_cast<int>(I_t)<=chemkin_cpp::chemkin::nii(); ++I_t){
		//		chemkin::ckraex(&I_t, &R_A);
		//		//cout<<I_t<<","<<R_A<<"\n";
		//		fout<<I_t<<","<<R_A<<"\n";
		//	}
		//	fout.clear(); fout.close();
	}


	void superPropagator::spe_concentration_w2f_pgt(double in_time, std::string str) const
	{
		std::ofstream fout((this->cwd_pgt + str).c_str());
		for (std::size_t i = 0; i < concentration_data_pgt.size(); ++i) {
			fout << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << this->evaluate_concentration_at_time(in_time, i);
			if (i != concentration_data_pgt.size() - 1) {
				fout << ",";
			}
		}
		fout << std::endl;

		fout.close(); fout.clear();
	}

	void superPropagator::spe_concentration_w2f_pgt(std::vector<double> time_in, std::string str) const
	{
		std::ofstream f_out((this->cwd_pgt + str).c_str());

		for (size_t i = 0; i < time_in.size(); ++i) {
			for (size_t j = 0; j < concentration_data_pgt.size(); ++j) {
				f_out << std::setprecision(std::numeric_limits<double>::max_digits10 + 1) << evaluate_concentration_at_time(time_in[i], j);
				if (j < (concentration_data_pgt.size() - 1)) {
					f_out << ",";
				}
			}
			f_out << std::endl;
		}

		f_out.close();
	}

	void superPropagator::read_configuration(double & Temp, double & Pressure, std::size_t nkk, double x[])
	{
		//Initialize molar fractions.
		for (std::size_t i = 0; i < nkk; ++i) { x[i] = 0.0; }
		double pressure_atm = this->pgt_pt.get<double>("chem_init.pressure_atm");

		//The pressure has to be in cgs unit, convert pascal(Pa) to dyne/square centimeter [dyn/cm**2]
		// 1atm= 101325 Pa
		// 1Pa= 10 dyn/cm**2
		Pressure = pressure_atm * 1013250;
		Temp = this->pgt_pt.get<double>("chem_init.init_temperature");

		//read with json_parser as property_tree
		//nice and easy
		BOOST_FOREACH(boost::property_tree::ptree::value_type const& key1, this->pgt_pt.get_child("chem_init.species_index_concentration")) {
			x[boost::lexical_cast<std::size_t>(key1.first)] = key1.second.get_value<double>();
		}

	}//read_configuration

	void superPropagator::initialize_lsode(rsp::my_time_t & dt, std::string outfile)
	{
		dt = this->pgt_pt.get<double>("lsode_init.dt");

		// Read 'lsode' parameters.
		lsodestore.atol = this->pgt_pt.get<double>("lsode_init.atol");
		lsodestore.rtol = this->pgt_pt.get<double>("lsode_init.rtol");
		lsodestore.mf = this->pgt_pt.get<int>("lsode_init.mf");
		lsodestore.jt = this->pgt_pt.get<int>("lsode_init.jt");
		lsodestore.itask = this->pgt_pt.get<int>("lsode_init.itask");
		lsodestore.iopt = this->pgt_pt.get<int>("lsode_init.iopt");
		lsodestore.itol = this->pgt_pt.get<int>("lsode_init.itol");
		lsodestore.istate = this->pgt_pt.get<int>("lsode_init.istate");

		lsodestore.deltaN1 = this->pgt_pt.get<int>("lsode_init.deltaN1");
		lsodestore.deltaN2 = this->pgt_pt.get<int>("lsode_init.deltaN2");

		int lrw_lsode, liw_lsode;
		std::ofstream fout((this->cwd_pgt + outfile).c_str(), std::ios::app);

		fout << std::endl;
		const int neq = chemkincpp_sr::chemkin::nkk() + 1;
		if (lsodestore.mf == 10)
		{
			lrw_lsode = 20 + 16 * neq + 1;
			liw_lsode = 20 + 1;
			fout << "  lrw_lsode =   " << lrw_lsode << std::endl;
			fout << "  liw_lsode =   " << liw_lsode << std::endl;
		}
		else if ((lsodestore.mf == 21) || (lsodestore.mf == 22))
		{
			lrw_lsode = 22 + 9 * neq + neq*neq + 1;
			liw_lsode = 20 + neq + 1;
			fout << "  lrw_lsode =   " << lrw_lsode << std::endl;
			fout << "  liw_lsode =   " << liw_lsode << std::endl;
		}
		else
		{
			fout << "  mf = " << lsodestore.mf << std::endl;
			fout << "  The program should not use this 'mf' value. Please check 'lsode.f' for more information.\n" << std::endl;
		}
		fout << std::endl << "  Please make sure that the values of constants 'lrw' and 'liw' in the Fortran subroutine 'cpplsode'"
			"\n  are larger than 'lrw_lsode' and 'liw_lsode' respectively." << std::endl;

		//	if(this->vm.count("lsode_init.jt")){
		boost::optional< boost::property_tree::ptree& > Ischild = this->pgt_pt.get_child_optional("lsode_init.jt");
		if (Ischild) {
			fout << "\n" << "  you are using the 12 November 2003 version of DLSODA: Livermore Solver for Ordinary Differential Equations,\n" <<
				"  with Automatic method switching for stiff and nonstiff problems\n";
			/*
			 *  4223 C          If the RWORK length is to be fixed, it should be at least
			 *  4224 C               MAX (LRN, LRS),
			 *  4225 C          where LRN and LRS are the RWORK lengths required when the
			 *  4226 C          current method is nonstiff or stiff, respectively.
			 *  4227 C
			 *  4228 C          The separate RWORK length requirements LRN and LRS are
			 *  4229 C          as follows:
			 *  4230 C          IF NEQ is constant and the maximum method orders have
			 *  4231 C          their default values, then
			 *  4232 C             LRN = 20 + 16*NEQ,
			 *  4233 C             LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2,
			 *  4234 C             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5.
			 */
			double lrn_pgtsode, lrs_pgtsode;
			lrn_pgtsode = 20 + 16 * neq;
			fout << "  lrn_pgtsode =   " << lrn_pgtsode << std::endl;

			if ((lsodestore.jt == 1) || (lsodestore.jt == 2))
			{
				lrs_pgtsode = 22 + 9 * neq + neq*neq;
				fout << "  lrs_pgtsode =   " << lrs_pgtsode << std::endl;
				fout << std::endl << "  Please make sure that the values of constants 'lrw' is bigger than " << std::max(lrn_pgtsode, lrs_pgtsode) << std::endl;
			}
			else if ((lsodestore.jt == 4) || (lsodestore.jt == 5))
			{
				//lrs_pgtsode= 22+10*neq+neq*neq;
				fout << "  lrs_pgtsode =   " << "check your jt value" << std::endl;
			}
			else
			{
				fout << "  jt = " << lsodestore.jt << std::endl;
				fout << "  The program should not use this 'jt' value. Please check \"include/fortran_lib/dlsode/opkdmain.f\" for more information.\n" << std::endl;
			}
			fout << std::endl;
		}

		fout.close();

	}	//initialize_lsode

	void superPropagator::initiate_propagator()
	{
		//initiate chemkin
		chemkincpp_sr::chemkin::chemkin_init();

		//Read the file named "chem.out", read in the chemical reactions index
		rsp::relationshipParser::read_reactionNetwork_chemkin_index_map(reactionNetwork_chemkin_index_map, this->cwd_pgt + std::string("/input/chem.out"));

		//initialize the pointer
		initialize_cubic_spline_pointer();

		//read configuration file "setting.json"
		boost::property_tree::read_json(this->cwd_pgt + std::string("/input/setting.json"), pgt_pt, std::locale());

		//set fast reactions, read fast inter-conversion reaction pairs from "setting.json"
		set_fast_reaction_pgt();
		//set the reaction rate of fast reactions to be zero
		//set_fast_reaction_rate_to_zero_pgt();
	}

	void propagator_sr::superPropagator::convert_molar_concentration_to_mole_fraction()
	{
		const int nkk = chemkincpp_sr::chemkin::nkk();
		double* x_t = new double[nkk]; double* c_t = new double[nkk];
		for (std::size_t k = 0; k < time_data_pgt.size(); ++k) {
			for (int i = 0; i < nkk; ++i) {
				c_t[i] = concentration_data_pgt[i][k];
			}
			chemkincpp_sr::chemkin::ckctx(c_t, x_t);
			//convert molar concentration to mole fraction
			for (int i = 0; i < nkk; ++i) {
				concentration_data_pgt[i][k] = x_t[i];
			}
		}
	}

	void propagator_sr::superPropagator::convert_mole_fraction_to_molar_concentration()
	{
		const int nkk = chemkincpp_sr::chemkin::nkk();
		double* x_t = new double[nkk]; double* c_t = new double[nkk];
		double pressure, temperature;

		for (std::size_t k = 0; k < time_data_pgt.size(); ++k) {
			pressure = this->pressure_data_pgt[k];
			temperature = this->temperature_data_pgt[k];
			double total = 0.0;
			for (int i = 0; i < nkk; ++i) {
				x_t[i] = concentration_data_pgt[i][k];
				total += x_t[i];
			}
			//normalization of mole fraction
			for (int i = 0; i < nkk; ++i)
				x_t[i] /= total;
			//convert mole fraction to molar concentration
			chemkincpp_sr::chemkin::ckxtcp(&pressure, &temperature, x_t, c_t);
			for (int i = 0; i < nkk; ++i) {
				concentration_data_pgt[i][k] = c_t[i];
			}
		}
	}

	void superPropagator::integrate_propensity_function_pgt()
	{

		spe_drc_int_data_pgt.resize(spe_drc_data_pgt.size());
		std::copy(spe_drc_data_pgt.begin(), spe_drc_data_pgt.end(), spe_drc_int_data_pgt.begin());


		for (size_t i = 0; i < spe_drc_int_data_pgt.size(); ++i)
		{
			//The first time interval
			//spe_drc_int_data_pgt[i][0] = spe_drc_data_pgt[i][0] * (time_data_pgt[1] - time_data_pgt[0]);
			spe_drc_int_data_pgt[i][0] = 0.0;
		}
		//The other time interval
		for (size_t i = 0; i < spe_drc_int_data_pgt.size(); ++i) {//[for
			for (size_t j = 1; j < spe_drc_int_data_pgt[0].size(); ++j) {
				////rectangle rule
				//spe_drc_int_data_pgt[i][j] = spe_drc_data_pgt[i][j] * (time_data_pgt[j] - time_data_pgt[j - 1]) + spe_drc_int_data_pgt[i][j - 1];

				//trapezoidal rule
				spe_drc_int_data_pgt[i][j] = 0.5 * (spe_drc_data_pgt[i][j] + spe_drc_data_pgt[i][j - 1]) * (time_data_pgt[j] - time_data_pgt[j - 1]) + spe_drc_int_data_pgt[i][j - 1];

			}
		}//for]

	}//int_propensity_function_pgt()


	bool superPropagator::init_spe_drc_pgt()
	{
		for (std::size_t i = 0; i < spe_drc_pgt.size(); ++i) {
			if (spe_drc_pgt[i] != 0)
				delete spe_drc_pgt[i];
		}
		//Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		spe_drc_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < spe_drc_data_pgt.size(); ++i)
		{
			spe_drc_pgt.push_back(new Linear_interp(time_data_pgt, spe_drc_data_pgt[i]));
		}
		return true;
	}//init_spe_drc_pgt()

	bool superPropagator::init_spe_drc_int_pgt()
	{
		for (std::size_t i = 0; i < spe_drc_int_pgt.size(); ++i) {
			if (spe_drc_int_pgt[i] != 0)
				delete spe_drc_int_pgt[i];
		}
		//Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		spe_drc_int_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < spe_drc_int_data_pgt.size(); ++i)
		{
			spe_drc_int_pgt.push_back(new Linear_interp(time_data_pgt, spe_drc_int_data_pgt[i]));
		}
		return true;
	}

	bool superPropagator::init_spe_drc_int_time_pgt()
	{
		for (std::size_t i = 0; i < spe_drc_int_time_pgt.size(); ++i) {
			if (spe_drc_int_time_pgt[i] != 0)
				delete spe_drc_int_time_pgt[i];
		}
		//Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
		spe_drc_int_time_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < spe_drc_int_data_pgt.size(); ++i)
		{
			spe_drc_int_time_pgt.push_back(new Linear_interp(spe_drc_int_data_pgt[i], time_data_pgt));
		}
		return true;
	}

	bool superPropagator::init_reaction_rate_pgt()
	{
		for (std::size_t i = 0; i < reaction_rate_pgt.size(); ++i) {
			if (reaction_rate_pgt[i] != 0)
				delete reaction_rate_pgt[i];
		}
		reaction_rate_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> spe_drc_pgt;
		for (std::size_t i = 0; i < reaction_rate_data_pgt.size(); ++i)
		{
			reaction_rate_pgt.push_back(new Linear_interp(time_data_pgt, reaction_rate_data_pgt[i]));
		}
		return true;
	}

	bool superPropagator::init_concentration_pgt()
	{
		for (size_t i = 0; i < concentration_pgt.size(); ++i) {
			if (concentration_pgt[i] != 0)
				delete concentration_pgt[i];
		}
		concentration_pgt.clear();
		//creat the cubic spline and link it to std::vector<Linear_interp*> concentraton_pgt
		for (size_t i = 0; i < concentration_data_pgt.size(); ++i) {
			concentration_pgt.push_back(new Linear_interp(time_data_pgt, concentration_data_pgt[i]));
		}

		return true;
	}

	//bool propagator_sr::superPropagator::init_time_total_concentration_pgt()
	//{
	//	if (total_concentration_pgt != 0)
	//		delete total_concentration_pgt;
	//	total_concentration_pgt = new Linear_interp(time_data_pgt, total_concentration_data_pgt);

	//	return true;
	//}


	bool superPropagator::init_time_temperature_pgt()
	{
		if (time_temperature_pgt != 0)
			delete time_temperature_pgt;
		time_temperature_pgt = new Linear_interp(time_data_pgt, temperature_data_pgt);

		return true;
	}

	bool superPropagator::init_temperature_time_pgt()
	{
		if (temperature_time_pgt != 0)
			delete temperature_time_pgt;
		temperature_time_pgt = new Linear_interp(temperature_data_pgt, time_data_pgt);

		return true;
	}

	bool superPropagator::init_time_pressure_pgt()
	{
		if (time_pressure_pgt != 0)
			delete time_pressure_pgt;
		time_pressure_pgt = new Linear_interp(time_data_pgt, pressure_data_pgt);

		return true;
	}

	double superPropagator::evaluate_temperature_at_time(double in_time) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		return  time_temperature_pgt->interp(in_time);
	}

	double superPropagator::evaluate_time_at_temperature(double in_temperature) const
	{
		if (in_temperature >= temperature_data_pgt.back())
			in_temperature = temperature_data_pgt.back();
		return temperature_time_pgt->interp(in_temperature);
	}

	double superPropagator::evaluate_pressure_at_time(double in_time) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		return  time_pressure_pgt->interp(in_time);
	}

	double superPropagator::evaluate_concentration_at_time(double in_time, size_t index) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((index < 0) || (index >= concentration_pgt.size()))
			index = concentration_pgt.size() - 1;
		return concentration_pgt[index]->interp(in_time);
	}

	double superPropagator::evaluate_spe_drc_at_time(double in_time, size_t index) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((index < 0) || (index >= spe_drc_pgt.size()))
			index = spe_drc_pgt.size() - 1;
		return spe_drc_pgt[index]->interp(in_time);
	}

	double superPropagator::evaluate_spe_drc_int_at_time(double in_time, size_t index) const
	{
		if ((index < 0) || (index >= spe_drc_int_pgt.size()))
			index = spe_drc_int_pgt.size() - 1;

		double diff = in_time - time_data_pgt.back();
		if (diff > 0) {
			in_time = time_data_pgt.back();
			return spe_drc_int_pgt[index]->interp(in_time) + diff*spe_drc_data_pgt[index].back();
		}
		else {
			return spe_drc_int_pgt[index]->interp(in_time);
		}
	}

	double superPropagator::evaluate_time_at_spe_drc_int(double integral, size_t index) const
	{
		//the last 3 species's propensity function is zero in the time region for H2-O2 system
		//can not solve for time
		//if((index<0)||(index+3>= spe_drc_int_time_pgt.size()))
		//index=spe_drc_int_time_pgt.size()-1-3;
		if ((index < 0) || (index >= spe_drc_int_time_pgt.size())) {
			//std::cout<<"index out of bound!"<<std::endl;
			std::cerr << "index out of bound!" << std::endl;
		}
		double diff = integral - spe_drc_int_data_pgt[index].back();
		if (diff > 0) {
			integral = spe_drc_int_data_pgt[index].back();
			//return spe_drc_int_time_pgt[index]->interp(integral)+diff/spe_drc_data_pgt[index].back();
			return time_data_pgt.back() + diff / spe_drc_data_pgt[index].back();
		}
		else {
			return spe_drc_int_time_pgt[index]->interp(integral);
		}
	}

	double superPropagator::evaluate_reaction_rate_at_time(double in_time, size_t index) const
	{
		if (in_time >= time_data_pgt.back())
			in_time = time_data_pgt.back();
		if ((index < 0) || (index >= reaction_rate_pgt.size()))
			index = reaction_rate_pgt.size() - 1;
		return reaction_rate_pgt[index]->interp(in_time);
	}




}//namespace propagator_sr


#endif


