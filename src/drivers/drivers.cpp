#ifndef __DRIVERS_CPP_
#define __DRIVERS_CPP_

#include "../../include/drivers/drivers.h"
#include "../../include/tools/my_debug/my_debug.h"

#include <iostream>
#include <limits>

#include <boost/filesystem.hpp>

#define PRINT_PRECISION 10


void driver::parse_parameters(const int argc, char ** argv, po::variables_map & vm, std::string &main_cwd, boost::property_tree::ptree & pt)
{
	// read vm
	clcf_parser clcf_parser_obj;
	clcf_parser_obj.parse(argc, argv);
	clcf_parser_obj.get_vm(vm);

	// read main_cwd
	if (vm.count("cwd")) {
		main_cwd = vm["cwd"].as<std::string>();
	}
	else {
		main_cwd = boost::filesystem::current_path().string();
	}
	//std::cout << "main_cwd:\t" << boost::filesystem::canonical(main_cwd) << std::endl;

	// read pt
	boost::property_tree::read_json(main_cwd + string("/input/setting.json"), pt, std::locale());

}

#ifdef __NO_USE_MPI_

void driver::generate_pathway_running_Monte_carlo_trajectory(const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//initialize MPI stuff
	int local_N = pt.get<int>("pathway.trajectoryNumber");


	fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
		main_cwd + std::string("/input/uncertainties.inp"));

	rnk::concreteReactionNetwork rnk_obj(uncertainties, 0, main_cwd);
	//rnk_obj.print_network();

	//double target_time_db = rnk_obj.return_temperature_target_time();
	double tau = pt.get<double>("time.tau");

	double init_time = 0.0 * tau;
	double end_time = pt.get<double>("pathway.tau")* tau;

	int trajectory_count_limit = pt.get<int>("pathway.trajectory_count_limit");

	//statistics
	statistics stat;

	// generate pathway one trajectory
	for (int i = 0; i < local_N; ++i) {
		auto str_t = rnk_obj.pathway_sim_once(init_time, end_time, rnk_obj.return_initial_spe(), pt.get<std::string>("pathway.atom_followed")); //vertex 2 is H2
		stat.insert_pathway_stat(str_t);
	}

	stat.sort_print_to_file_stat(main_cwd + std::string("/output/pathway_stat.csv"), trajectory_count_limit);

}

void driver::generate_species_pathway_running_Monte_carlo_trajectory(const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//initialize MPI stuff
	int local_N = pt.get<int>("pathway.trajectoryNumber");


	fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
		main_cwd + std::string("/input/uncertainties.inp"));

	rnk::concreteReactionNetwork rnk_obj(uncertainties, 0, main_cwd);
	//rnk_obj.print_network();

	//double target_time_db = rnk_obj.return_temperature_target_time();
	double tau = pt.get<double>("time.tau");

	double init_time = 0.0 * tau;
	double end_time = pt.get<double>("pathway.tau")* tau;

	int trajectory_count_limit = pt.get<int>("pathway.trajectory_count_limit");

	//statistics
	statistics stat;

	// generate pathway one trajectory
	for (int i = 0; i < local_N; ++i) {
		auto str_t = rnk_obj.species_pathway_sim_once(init_time, end_time, rnk_obj.return_initial_spe(), pt.get<std::string>("pathway.atom_followed")); //vertex 2 is H2
		stat.insert_pathway_stat(str_t);
	}

	stat.sort_print_to_file_stat(main_cwd + std::string("/output/species_pathway_stat.csv"), trajectory_count_limit);
}

void driver::evaluate_path_integral_over_time(const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;
	//time points, every pathway has its own set of time points
	std::vector<std::vector<double> > time_Mat;


	fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
		main_cwd + std::string("/input/uncertainties.inp"));

	//get the pathway name only on the processor 0
	std::vector<std::string> pathway_vec_t;
	pathwayHandler::get_pathway(main_cwd + std::string("/output/pathway_name_candidate.csv"), pathway_vec_t,
		std::numeric_limits<int>::max() - 1000); //all pathways


	std::vector<std::size_t> topN_vec;
	for (auto key1 : pt.get_child("pathway.topN")) {
		topN_vec.push_back(key1.second.get_value<size_t>());
	}
	std::size_t topN = topN_vec.front();
	topN = (topN <= pathway_vec_t.size()) ? topN : pathway_vec_t.size();

	if (pt.get<std::string>("pathway.pathwayEndWith") == "ALL") {
		pathway_vec.assign(pathway_vec_t.begin(), pathway_vec_t.begin() + topN);
	}
	else {
		pathwayHandler::pathway_ends_with(pt.get<std::string>("pathway.pathwayEndWith"), pathway_vec_t, pathway_vec,
			topN); //topN pathways
	}

	//initialize the time points we are going to calculate
	//read file pathway_time.csv
	time_Mat = fileIO::fileIO::read_topN_line_csv_matrix(main_cwd + std::string("/output/pathway_time_candidate.csv"), topN);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(pathway_vec.size(), std::vector<double>(time_Mat[0].size(), 0.0));

	////create rxn_network, generate pathway
	size_t trajectoryNumber_local = pt.get<std::size_t>("pathway.trajectoryNumber");

	//different seed for different core/CPU
	rnk::concreteReactionNetwork rnk_obj(uncertainties, 0, main_cwd);

	double tau = pt.get<double>("time.tau");

	// evaluate path integral on each core
	double p_p_db = 0.0;
	std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;

	for (size_t i = 0; i < prob_Mat.size(); ++i) {
		for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
			for (size_t k = 0; k < trajectoryNumber_local; ++k) {
				rnk_obj.parse_pathway_to_vector(pathway_vec[i], spe_vec, reaction_vec);
				p_p_db = rnk_obj.pathway_prob_input_pathway_sim_once(0.0, time_Mat[i][j] * tau,
					spe_vec, reaction_vec, pt.get<std::string>("pathway.atom_followed"));
				prob_Mat[i][j] += p_p_db / trajectoryNumber_local;
			}
		}
	}

	std::ofstream fout((main_cwd + std::string("/output/pathway_prob.csv")).c_str(), std::ofstream::out);
	for (size_t i = 0; i < prob_Mat.size(); ++i) {
		for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
			fout << setprecision(PRINT_PRECISION) << prob_Mat[i][j];
			if (j != (prob_Mat[0].size() - 1)) {
				fout << ",";
			}
		}
		//std::cout<<std::endl;
		fout << std::endl;
	}

	fout.clear();
	fout.close();
	fout.open((main_cwd + std::string("/output/pathway_name_selected.csv")).c_str(), std::ofstream::out);
	for (size_t i = 0; i < pathway_vec.size(); ++i) {
		fout << pathway_vec[i] << std::endl;
	}
	fout.clear();
	fout.close();

}

void driver::evaluate_species_path_integral_over_time(const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;
	//time points, every pathway has its own set of time points
	std::vector<std::vector<double> > time_Mat;


	fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
		main_cwd + std::string("/input/uncertainties.inp"));

	//get the pathway name only on the processor 0
	std::vector<std::string> pathway_vec_t;
	pathwayHandler::get_pathway(main_cwd + std::string("/output/species_pathway_name_candidate.csv"), pathway_vec_t,
		std::numeric_limits<int>::max() - 1000); //all pathways


	std::vector<std::size_t> topN_vec;
	for (auto key1 : pt.get_child("pathway.topN")) {
		topN_vec.push_back(key1.second.get_value<size_t>());
	}
	std::size_t topN = topN_vec.front();
	topN = (topN <= pathway_vec_t.size()) ? topN : pathway_vec_t.size();

	if (pt.get<std::string>("pathway.pathwayEndWith") == "ALL") {
		pathway_vec.assign(pathway_vec_t.begin(), pathway_vec_t.begin() + topN);
	}
	else {
		pathwayHandler::pathway_ends_with(pt.get<std::string>("pathway.pathwayEndWith"), pathway_vec_t, pathway_vec,
			topN); //topN pathways
	}

	//initialize the time points we are going to calculate
	//read file pathway_time.csv
	time_Mat = fileIO::fileIO::read_topN_line_csv_matrix(main_cwd + std::string("/output/species_pathway_time_candidate.csv"), topN);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(pathway_vec.size(), std::vector<double>(time_Mat[0].size(), 0.0));

	////create rxn_network, generate pathway
	size_t trajectoryNumber_local = pt.get<std::size_t>("pathway.trajectoryNumber");

	//different seed for different core/CPU
	rnk::concreteReactionNetwork rnk_obj(uncertainties, 0, main_cwd);

	double tau = pt.get<double>("time.tau");

	// evaluate path integral on each core
	double p_p_db = 0.0;
	std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;

	for (size_t i = 0; i < prob_Mat.size(); ++i) {
		for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
			for (size_t k = 0; k < trajectoryNumber_local; ++k) {
				rnk_obj.parse_pathway_to_vector(pathway_vec[i], spe_vec, reaction_vec);
				p_p_db = rnk_obj.species_pathway_prob_input_pathway_sim_once(0.0, time_Mat[i][j] * tau,
					spe_vec, reaction_vec, pt.get<std::string>("pathway.atom_followed"));
				prob_Mat[i][j] += p_p_db / trajectoryNumber_local;
			}
		}
	}

	std::ofstream fout((main_cwd + std::string("/output/species_pathway_prob.csv")).c_str(), std::ofstream::out);
	for (size_t i = 0; i < prob_Mat.size(); ++i) {
		for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
			fout << setprecision(PRINT_PRECISION) << prob_Mat[i][j];
			if (j != (prob_Mat[0].size() - 1)) {
				fout << ",";
			}
		}
		//std::cout<<std::endl;
		fout << std::endl;
	}

	fout.clear();
	fout.close();
	fout.open((main_cwd + std::string("/output/species_pathway_name_selected.csv")).c_str(), std::ofstream::out);
	for (size_t i = 0; i < pathway_vec.size(); ++i) {
		fout << pathway_vec[i] << std::endl;
	}
	fout.clear();
	fout.close();
}

void driver::evaluate_path_AT_over_time(const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;

	fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
		main_cwd + std::string("/input/uncertainties.inp"));

	//get the pathway name only on the processor 0
	std::vector<std::string> pathway_vec_t;
	pathwayHandler::get_pathway(main_cwd + std::string("/output/species_pathway_name_candidate.csv"), pathway_vec_t,
		std::numeric_limits<int>::max() - 1000); //all pathways

	std::vector<std::size_t> topN_vec;
	for (auto key1 : pt.get_child("pathway.topN")) {
		topN_vec.push_back(key1.second.get_value<size_t>());
	}
	std::size_t topN = topN_vec.front();
	topN = (topN <= pathway_vec_t.size()) ? topN : pathway_vec_t.size();

	if (pt.get<std::string>("pathway.pathwayEndWith") == "ALL") {
		pathway_vec.assign(pathway_vec_t.begin(), pathway_vec_t.begin() + topN);
	}
	else {
		pathwayHandler::pathway_ends_with(pt.get<std::string>("pathway.pathwayEndWith"), pathway_vec_t, pathway_vec,
			topN); //topN pathways
	}

	size_t trajectoryNumber_local = pt.get<std::size_t>("pathway.trajectoryNumber");
	std::vector< std::vector<double> > path_AT_vec(pathway_vec.size(), std::vector<double>(trajectoryNumber_local, 0.0));

	//different seed for different core/CPU
	rnk::concreteReactionNetwork rnk_obj(uncertainties, 0, main_cwd);

	double time = pt.get<double>("time.tau") * pt.get<double>("pathway.tau");

	// evaluate path AT on each core
	std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;

	for (std::size_t i = 0; i < pathway_vec.size(); ++i) {
		rnk_obj.parse_pathway_to_vector(pathway_vec[i], spe_vec, reaction_vec);
		for (size_t j = 0; j < trajectoryNumber_local; ++j) {
			path_AT_vec[i][j] = rnk_obj.pathway_AT_input_pathway_sim_once(0.0, time, spe_vec, reaction_vec);
		}
	}

	std::ofstream fout((main_cwd + std::string("/output/species_pathway_AT.csv")).c_str(), std::ofstream::out);
	for (size_t i = 0; i < path_AT_vec.size(); ++i) {
		for (size_t j = 0; j < path_AT_vec[0].size(); ++j) {
			fout << setprecision(PRINT_PRECISION) << path_AT_vec[i][j];
			if (j != (path_AT_vec[0].size() - 1)) {
				fout << ",";
			}
		}
		fout << std::endl;
	}

	fout.clear();
	fout.close();
}

void driver::evaluate_path_AT_with_SP_over_time(const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;

	fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
		main_cwd + std::string("/input/uncertainties.inp"));

	//get the pathway name only on the processor 0
	std::vector<std::string> pathway_vec_t;
	pathwayHandler::get_pathway(main_cwd + std::string("/output/species_pathway_name_candidate.csv"), pathway_vec_t,
		std::numeric_limits<int>::max() - 1000); //all pathways

	std::vector<std::size_t> topN_vec;
	for (auto key1 : pt.get_child("pathway.topN")) {
		topN_vec.push_back(key1.second.get_value<size_t>());
	}
	std::size_t topN = topN_vec.front();
	topN = (topN <= pathway_vec_t.size()) ? topN : pathway_vec_t.size();

	if (pt.get<std::string>("pathway.pathwayEndWith") == "ALL") {
		pathway_vec.assign(pathway_vec_t.begin(), pathway_vec_t.begin() + topN);
	}
	else {
		pathwayHandler::pathway_ends_with(pt.get<std::string>("pathway.pathwayEndWith"), pathway_vec_t, pathway_vec,
			topN); //topN pathways
	}

	size_t trajectoryNumber_local = pt.get<std::size_t>("pathway.trajectoryNumber");
	std::vector< std::vector<double> > path_AT_mat(pathway_vec.size(), std::vector<double>(trajectoryNumber_local, 0.0));
	std::vector< std::vector<double> > path_AT_prob_mat(pathway_vec.size(), std::vector<double>(trajectoryNumber_local, 0.0));

	//different seed for different core/CPU
	rnk::concreteReactionNetwork rnk_obj(uncertainties, 0, main_cwd);

	double time = pt.get<double>("time.tau") * pt.get<double>("pathway.tau");

	// evaluate path AT on each core
	std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;

	for (std::size_t i = 0; i < pathway_vec.size(); ++i) {
		rnk_obj.parse_pathway_to_vector(pathway_vec[i], spe_vec, reaction_vec);
		for (size_t j = 0; j < trajectoryNumber_local; ++j) {
			std::tie(path_AT_mat[i][j], path_AT_prob_mat[i][j]) = rnk_obj.pathway_AT_with_SP_input_pathway_sim_once(0.0, time, spe_vec, reaction_vec);
		}
	}

	std::ofstream fout1((main_cwd + std::string("/output/species_pathway_AT_with_SP.csv")).c_str(), std::ofstream::out);
	std::ofstream fout2((main_cwd + std::string("/output/species_pathway_SP.csv")).c_str(), std::ofstream::out);
	for (size_t i = 0; i < path_AT_mat.size(); ++i) {
		for (size_t j = 0; j < path_AT_mat[0].size(); ++j) {
			fout1 << setprecision(PRINT_PRECISION) << path_AT_mat[i][j];
			fout2 << setprecision(PRINT_PRECISION) << path_AT_prob_mat[i][j];
			if (j != (path_AT_mat[0].size() - 1)) {
				fout1 << ",";
				fout2 << ",";
			}
		}
		fout1 << std::endl;
		fout2 << std::endl;
	}

	fout1.clear();
	fout1.close();
	fout2.clear();
	fout2.close();
}

#endif // __NO_USE_MPI_


#if defined(__USE_MPI_) 

void driver::write_concentration_at_time_to_file(const boost::mpi::communicator & world, std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	if (world.rank() == 0) {
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));
	}//if

	 //boradcast
	broadcast(world, uncertainties, 0);

	if (world.rank() == 0) {
		pgt::dlsodePropagator pgt_obj(uncertainties, main_cwd);
		if (pt.get<std::string>("propagator.convert_molar_concentration_to_mole_fraction") == std::string("yes")) {
			pgt_obj.convert_molar_concentration_to_mole_fraction();
			pgt_obj.spe_concentration_w2f_pgt(pt.get<double>("time.tau") * pt.get<double>("pathway.tau"),
				pt.get<std::string>("pathway.tau") + std::string("_dlsode_fraction"));
		}
		else
			pgt_obj.spe_concentration_w2f_pgt(pt.get<double>("time.tau") * pt.get<double>("pathway.tau"),
				pt.get<std::string>("pathway.tau") + std::string("_dlsode_M"));
	}
}

void driver::solve_ODEs_for_concentration_using_LSODE(const boost::mpi::communicator & world, std::string & main_cwd, const boost::property_tree::ptree &pt)
{
	std::vector<double> uncertainties;

	if (world.rank() == 0) {
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));

	}//if

	 //boradcast
	broadcast(world, uncertainties, 0);

	if (world.rank() == 0) {
		pgt::dlsodePropagator pgt_obj(uncertainties, main_cwd);
		if (pt.get<std::string>("propagator.convert_molar_concentration_to_mole_fraction") == std::string("yes"))
		{
			pgt_obj.convert_molar_concentration_to_mole_fraction();
			pgt_obj.w2f_pgt("dlsode_fraction");
		}
		else
			pgt_obj.w2f_pgt("dlsode_M");
		//pgt_obj.convert_molar_concentration_to_mole_fraction();

	}
}

void driver::solve_ODEs_for_concentration_using_SSA(const boost::mpi::communicator & world, std::string & main_cwd, const boost::property_tree::ptree &pt)
{
	std::vector<double> uncertainties;

	if (world.rank() == 0) {
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));
	}//if

	 //boradcast
	broadcast(world, uncertainties, 0);

	if (world.rank() == 0) {
		pgt::ssaPropagator pgt_obj(uncertainties, world.rank(), main_cwd);
		pgt_obj.w2f_pgt("ssa_number");
	}
}

void driver::generate_pathway_running_Monte_carlo_trajectory(const boost::mpi::communicator & world, const std::string &main_cwd, const boost::property_tree::ptree &pt)
{
	std::vector<double> uncertainties;

	//initialize MPI stuff
	int trajectoryNumber_total = pt.get<std::size_t>("pathway.trajectoryNumber");
	int P = world.size();
	int local_N = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	if (world.rank() == 0) {
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));
	}//if

	 //boradcast
	broadcast(world, uncertainties, 0);

	rnk::concreteReactionNetwork rnk_obj(uncertainties, world.rank(), main_cwd);
	//rnk_obj.print_network();

	//double target_time_db = rnk_obj.return_temperature_target_time();
	double tau = pt.get<double>("time.tau");

	double init_time = 0.0 * tau;
	double end_time = pt.get<double>("pathway.tau")* tau;

	int trajectory_count_limit = pt.get<int>("pathway.trajectory_count_limit");

	//statistics
	statistics stat;
	std::string str_t;

	// generate pathway one trajectory
	for (int i = 0; i < local_N; ++i) {
		str_t = rnk_obj.pathway_sim_once(init_time, end_time, rnk_obj.return_initial_spe(), pt.get<std::string>("pathway.atom_followed")); //vertex 2 is H2
		stat.insert_pathway_stat(str_t);
	}

	//map reduce
	std::map<std::string, int> result;
	reduce(world, stat.get_pathway_unordered_map(),
		result, merge_maps(), 0);

	if (world.rank() == 0) {
		stat.insert_unordered_map(result);
		stat.sort_print_to_file_stat(main_cwd + std::string("/output/pathway_stat.csv"), trajectory_count_limit);
	}
}

void driver::generate_species_pathway_running_Monte_carlo_trajectory(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//initialize MPI stuff
	int trajectoryNumber_total = pt.get<std::size_t>("pathway.trajectoryNumber");
	int P = world.size();
	int local_N = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	if (world.rank() == 0) {
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));
	}//if

	 //boradcast
	broadcast(world, uncertainties, 0);

	rnk::concreteReactionNetwork rnk_obj(uncertainties, world.rank(), main_cwd);
	//rnk_obj.print_network();

	//double target_time_db = rnk_obj.return_temperature_target_time();
	double tau = pt.get<double>("time.tau");

	double init_time = 0.0 * tau;
	double end_time = pt.get<double>("pathway.tau")* tau;

	int trajectory_count_limit = pt.get<int>("pathway.trajectory_count_limit");

	//statistics
	statistics stat;
	std::string str_t;

	// generate pathway one trajectory
	for (int i = 0; i < local_N; ++i) {
		str_t = rnk_obj.species_pathway_sim_once(init_time, end_time, rnk_obj.return_initial_spe(), pt.get<std::string>("pathway.atom_followed")); //vertex 2 is H2
		stat.insert_pathway_stat(str_t);
	}

	//map reduce
	std::map<std::string, int> result;
	reduce(world, stat.get_pathway_unordered_map(),
		result, merge_maps(), 0);

	if (world.rank() == 0) {
		stat.insert_unordered_map(result);
		stat.sort_print_to_file_stat(main_cwd + std::string("/output/species_pathway_stat.csv"), trajectory_count_limit);
	}
}

void driver::evaluate_path_integral_over_time(const boost::mpi::communicator & world, const std::string &main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;
	//time points, every pathway has its own set of time points
	std::vector<std::vector<double> > time_Mat;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));

		//get the pathway name only on the processor 0
		std::vector<std::string> pathway_vec_t;
		pathwayHandler::get_pathway(main_cwd + std::string("/output/pathway_name_candidate.csv"), pathway_vec_t,
			std::numeric_limits<int>::max() - 1000); //all pathways


		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		std::size_t topN = topN_vec.front();
		topN = (topN <= pathway_vec_t.size()) ? topN : pathway_vec_t.size();

		if (pt.get<std::string>("pathway.pathwayEndWith") == "ALL") {
			pathway_vec.assign(pathway_vec_t.begin(), pathway_vec_t.begin() + topN);
		}
		else {
			pathwayHandler::pathway_ends_with(pt.get<std::string>("pathway.pathwayEndWith"), pathway_vec_t, pathway_vec,
				topN); //topN pathways
		}

		//initialize the time points we are going to calculate
		//read file pathway_time.csv
		time_Mat = fileIO::fileIO::read_topN_line_csv_matrix(main_cwd + std::string("/output/pathway_time_candidate.csv"), topN);
	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, pathway_vec, 0);
	broadcast(world, time_Mat, 0);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(pathway_vec.size(), std::vector<double>(time_Mat[0].size(), 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(pathway_vec.size(), std::vector<double>(time_Mat[0].size(), 0.0));

	////create rxn_network, generate pathway
	size_t trajectoryNumber_total = pt.get<std::size_t>("pathway.trajectoryNumber");
	size_t P = world.size();
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	//different seed for different core/CPU
	rnk::concreteReactionNetwork rnk_obj(uncertainties, world.rank(), main_cwd);

	double tau = pt.get<double>("time.tau");

	// evaluate path integral on each core
	double pathway_prob_db_t = 0.0;
	std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;

	for (size_t i = 0; i < prob_Mat.size(); ++i) {
		for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
			for (size_t k = 0; k < trajectoryNumber_local; ++k) {
				rnk_obj.parse_pathway_to_vector(pathway_vec[i], spe_vec, reaction_vec);
				pathway_prob_db_t = rnk_obj.pathway_prob_input_pathway_sim_once(0.0, time_Mat[i][j] * tau,
					spe_vec, reaction_vec, pt.get<std::string>("pathway.atom_followed"));
				prob_Mat[i][j] += pathway_prob_db_t / trajectoryNumber_total;
			}
		}
	}
	// map reduce
	for (size_t i = 0; i < prob_Mat.size(); ++i) {
		for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
			reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
		}
	}

	if (world.rank() == 0) {
		std::ofstream fout((main_cwd + std::string("/output/pathway_prob.csv")).c_str(), std::ofstream::out);
		for (size_t i = 0; i < prob_Mat_reduce.size(); ++i) {
			for (size_t j = 0; j < prob_Mat_reduce[0].size(); ++j) {
				fout << setprecision(PRINT_PRECISION) << prob_Mat_reduce[i][j];
				if (j != (prob_Mat_reduce[0].size() - 1)) {
					fout << ",";
				}
			}
			fout << std::endl;
		}

		fout.clear();
		fout.close();
		fout.open((main_cwd + std::string("/output/pathway_name_selected.csv")).c_str(), std::ofstream::out);
		for (size_t i = 0; i < pathway_vec.size(); ++i) {
			fout << pathway_vec[i] << std::endl;
		}
		fout.clear();
		fout.close();
	}

}

void driver::evaluate_species_path_integral_over_time(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;
	//time points, every pathway has its own set of time points
	std::vector<std::vector<double> > time_Mat;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));

		//get the pathway name only on the processor 0
		std::vector<std::string> pathway_vec_t;
		pathwayHandler::get_pathway(main_cwd + std::string("/output/species_pathway_name_candidate.csv"), pathway_vec_t,
			std::numeric_limits<int>::max() - 1000); //all pathways


		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		std::size_t topN = topN_vec.front();
		topN = (topN <= pathway_vec_t.size()) ? topN : pathway_vec_t.size();

		if (pt.get<std::string>("pathway.pathwayEndWith") == "ALL") {
			pathway_vec.assign(pathway_vec_t.begin(), pathway_vec_t.begin() + topN);
		}
		else {
			pathwayHandler::pathway_ends_with(pt.get<std::string>("pathway.pathwayEndWith"), pathway_vec_t, pathway_vec,
				topN); //topN pathways
		}

		//initialize the time points we are going to calculate
		//read file pathway_time.csv
		time_Mat = fileIO::fileIO::read_topN_line_csv_matrix(main_cwd + std::string("/output/species_pathway_time_candidate.csv"), topN);
	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, pathway_vec, 0);
	broadcast(world, time_Mat, 0);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(pathway_vec.size(), std::vector<double>(time_Mat[0].size(), 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(pathway_vec.size(), std::vector<double>(time_Mat[0].size(), 0.0));

	////create rxn_network, generate pathway
	size_t trajectoryNumber_total = pt.get<std::size_t>("pathway.trajectoryNumber");
	size_t P = world.size();
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	//different seed for different core/CPU
	rnk::concreteReactionNetwork rnk_obj(uncertainties, world.rank(), main_cwd);

	double tau = pt.get<double>("time.tau");

	// evaluate path integral on each core
	double pathway_prob_db_t = 0.0;
	std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;

	for (size_t i = 0; i < prob_Mat.size(); ++i) {
		for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
			for (size_t k = 0; k < trajectoryNumber_local; ++k) {
				rnk_obj.parse_pathway_to_vector(pathway_vec[i], spe_vec, reaction_vec);
				pathway_prob_db_t = rnk_obj.species_pathway_prob_input_pathway_sim_once(0.0, time_Mat[i][j] * tau,
					spe_vec, reaction_vec, pt.get<std::string>("pathway.atom_followed"));
				prob_Mat[i][j] += pathway_prob_db_t / trajectoryNumber_total;
			}
		}
	}
	// map reduce
	for (size_t i = 0; i < prob_Mat.size(); ++i) {
		for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
			reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
		}
	}

	if (world.rank() == 0) {
		std::ofstream fout((main_cwd + std::string("/output/species_pathway_prob.csv")).c_str(), std::ofstream::out);
		for (size_t i = 0; i < prob_Mat_reduce.size(); ++i) {
			for (size_t j = 0; j < prob_Mat_reduce[0].size(); ++j) {
				fout << setprecision(PRINT_PRECISION) << prob_Mat_reduce[i][j];
				if (j != (prob_Mat_reduce[0].size() - 1)) {
					fout << ",";
				}
			}
			fout << std::endl;
		}

		fout.clear();
		fout.close();
		fout.open((main_cwd + std::string("/output/species_pathway_name_selected.csv")).c_str(), std::ofstream::out);
		for (size_t i = 0; i < pathway_vec.size(); ++i) {
			fout << pathway_vec[i] << std::endl;
		}
		fout.clear();
		fout.close();
	}

}

void driver::evaluate_path_AT_over_time(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{

	if (world.rank() == 0) {
		std::vector<double> uncertainties;

		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;

		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));

		//get the pathway name only on the processor 0
		std::vector<std::string> pathway_vec_t;
		pathwayHandler::get_pathway(main_cwd + std::string("/output/species_pathway_name_candidate.csv"), pathway_vec_t,
			std::numeric_limits<int>::max() - 1000); //all pathways

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		std::size_t topN = topN_vec.front();
		topN = (topN <= pathway_vec_t.size()) ? topN : pathway_vec_t.size();

		if (pt.get<std::string>("pathway.pathwayEndWith") == "ALL") {
			pathway_vec.assign(pathway_vec_t.begin(), pathway_vec_t.begin() + topN);
		}
		else {
			pathwayHandler::pathway_ends_with(pt.get<std::string>("pathway.pathwayEndWith"), pathway_vec_t, pathway_vec,
				topN); //topN pathways
		}

		size_t trajectoryNumber_local = pt.get<std::size_t>("pathway.trajectoryNumber");
		std::vector< std::vector<double> > path_AT_vec(pathway_vec.size(), std::vector<double>(trajectoryNumber_local, 0.0));

		//different seed for different core/CPU
		rnk::concreteReactionNetwork rnk_obj(uncertainties, 0, main_cwd);

		double time = pt.get<double>("time.tau") * pt.get<double>("pathway.tau");

		// evaluate path AT on each core
		std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;

		for (std::size_t i = 0; i < pathway_vec.size(); ++i) {
			rnk_obj.parse_pathway_to_vector(pathway_vec[i], spe_vec, reaction_vec);
			for (size_t j = 0; j < trajectoryNumber_local; ++j) {
				path_AT_vec[i][j] = rnk_obj.pathway_AT_input_pathway_sim_once(0.0, time, spe_vec, reaction_vec);
			}
		}

		std::ofstream fout((main_cwd + std::string("/output/species_pathway_AT.csv")).c_str(), std::ofstream::out);
		for (size_t i = 0; i < path_AT_vec.size(); ++i) {
			for (size_t j = 0; j < path_AT_vec[0].size(); ++j) {
				fout << setprecision(PRINT_PRECISION) << path_AT_vec[i][j];
				if (j != (path_AT_vec[0].size() - 1)) {
					fout << ",";
				}
			}
			fout << std::endl;
		}

		fout.clear();
		fout.close();
	}//world.rank() ==0

}

void driver::evaluate_path_AT_with_SP_over_time(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	if (world.rank() == 0) {
		std::vector<double> uncertainties;

		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;

		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));

		//get the pathway name only on the processor 0
		std::vector<std::string> pathway_vec_t;
		pathwayHandler::get_pathway(main_cwd + std::string("/output/species_pathway_name_candidate.csv"), pathway_vec_t,
			std::numeric_limits<int>::max() - 1000); //all pathways

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		std::size_t topN = topN_vec.front();
		topN = (topN <= pathway_vec_t.size()) ? topN : pathway_vec_t.size();

		if (pt.get<std::string>("pathway.pathwayEndWith") == "ALL") {
			pathway_vec.assign(pathway_vec_t.begin(), pathway_vec_t.begin() + topN);
		}
		else {
			pathwayHandler::pathway_ends_with(pt.get<std::string>("pathway.pathwayEndWith"), pathway_vec_t, pathway_vec,
				topN); //topN pathways
		}

		size_t trajectoryNumber_local = pt.get<std::size_t>("pathway.trajectoryNumber");
		std::vector< std::vector<double> > path_AT_mat(pathway_vec.size(), std::vector<double>(trajectoryNumber_local, 0.0));
		std::vector< std::vector<double> > path_AT_prob_mat(pathway_vec.size(), std::vector<double>(trajectoryNumber_local, 0.0));

		//different seed for different core/CPU
		rnk::concreteReactionNetwork rnk_obj(uncertainties, 0, main_cwd);

		double time = pt.get<double>("time.tau") * pt.get<double>("pathway.tau");

		// evaluate path AT on each core
		std::vector<rsp::index_int_t> spe_vec; std::vector<rsp::index_int_t> reaction_vec;

		for (std::size_t i = 0; i < pathway_vec.size(); ++i) {
			rnk_obj.parse_pathway_to_vector(pathway_vec[i], spe_vec, reaction_vec);
			for (size_t j = 0; j < trajectoryNumber_local; ++j) {
				std::tie(path_AT_mat[i][j], path_AT_prob_mat[i][j]) = rnk_obj.pathway_AT_with_SP_input_pathway_sim_once(0.0, time, spe_vec, reaction_vec);
			}
		}

		std::ofstream fout1((main_cwd + std::string("/output/species_pathway_AT_with_SP.csv")).c_str(), std::ofstream::out);
		std::ofstream fout2((main_cwd + std::string("/output/species_pathway_SP.csv")).c_str(), std::ofstream::out);
		for (size_t i = 0; i < path_AT_mat.size(); ++i) {
			for (size_t j = 0; j < path_AT_mat[0].size(); ++j) {
				fout1 << setprecision(PRINT_PRECISION) << path_AT_mat[i][j];
				fout2 << setprecision(PRINT_PRECISION) << path_AT_prob_mat[i][j];
				if (j != (path_AT_mat[0].size() - 1)) {
					fout1 << ",";
					fout2 << ",";
				}
			}
			fout1 << std::endl;
			fout2 << std::endl;
		}

		fout1.clear();
		fout1.close();
		fout2.clear();
		fout2.close();
	}//world.rank() ==0
}

void driver::speciation_evaluate_concentrations_for_different_sets_rate_coefficients(const boost::mpi::communicator & world, const std::string &main_cwd, const boost::property_tree::ptree &pt)
{
	std::vector<double> uncertainties;

	size_t number_of_Ks = pt.get<std::size_t>("speciation.number_of_Ks");
	size_t P = world.size();
	//size_t number_of_Ks_local = get_num_block_decomposition_2(world.rank(), number_of_Ks, P);

	for (int ith_k = get_first_block_decomposition_2(world.rank(), number_of_Ks, P); ith_k <= get_last_block_decomposition_2(world.rank(), number_of_Ks, P); ++ith_k) {
		fileIO::fileIO::read_generate_uncertainties_w2f_random(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"), main_cwd + std::string("/output/uncertainties_random.csv"),
			main_cwd + std::string("/input/chem.out"), ith_k);

		std::cout << "ith_k:\t" << ith_k << std::endl;

		pgt::dlsodePropagator pgt_obj(uncertainties, main_cwd);

		std::vector<std::vector<double> > time_mat = fileIO::fileIO::read_topN_line_csv_matrix(main_cwd + std::string("/input/speciation_time_points.csv"), 1);
		std::vector<double> time_point = time_mat[0];

		double target_time_db = pgt_obj.return_temperature_target_time();

		for (std::size_t i = 0; i < time_point.size(); ++i) {
			time_point[i] *= target_time_db;
		}

		// print concentration to file "./output/conc.csv"
		pgt_obj.spe_concentration_w2f_pgt(time_point, std::string("/output/conc_") + boost::lexical_cast<std::string>(ith_k) + std::string(".csv"));

		// print target time to file
		std::ofstream tt_fout((main_cwd + std::string("/output/target_time_") + boost::lexical_cast<std::string>(ith_k) + std::string(".csv")).c_str());
		tt_fout << std::setprecision(16) << target_time_db;
		tt_fout.close();
	}

}

void driver::ODE_solver_MC_trajectory_single_core(const boost::mpi::communicator & world, const std::string &main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;
	if (world.rank() == 0) {
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties,
			main_cwd + std::string("/input/uncertainties.inp"));

	}//if

	 //boradcast
	broadcast(world, uncertainties, 0);

	if (world.rank() == 0) {
		rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);

		double tau = pt.get<double>("time.tau");
		double init_time = 0.0 * tau;
		double end_time = 1.0 * tau;

		rnkODEs_obj.ODE_pathway_sim_N(init_time, end_time, pt.get<int>("pathway.trajectoryNumber"));
		rnkODEs_obj.w2f_pgt("SOHR");
	}
}

void driver::ODE_solver_MC_trajectory_s_ct_np_parallel(const boost::mpi::communicator & world, const std::string &main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;

	//initialize MPI stuff
	std::size_t iterationNumber;
	size_t trajectoryNumber_total;
	double init_time;
	double end_time;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		trajectoryNumber_total = pt.get<std::size_t>("pathway.trajectoryNumber");
		init_time = 0.0*pt.get<double>("time.tau");
		end_time = pt.get<double>("pathway.tau")*pt.get<double>("time.tau");
	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_total, 0);
	broadcast(world, init_time, 0);
	broadcast(world, end_time, 0);


	//create reaction_network
	size_t P = world.size();
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	//different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	//probability
	std::vector<double> prob = rnkODEs_obj.get_initial_concentration();

	std::pair<std::size_t, std::size_t> conc_data_size = rnkODEs_obj.get_size_of_concentration_data();

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(conc_data_size.first, std::vector<double>(conc_data_size.second, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(conc_data_size.first,
		std::vector<double>(conc_data_size.second, 0.0));

	for (std::size_t ni = 0; ni < iterationNumber; ++ni) {
		//not first iteration
		if (ni != 0) {
			rnkODEs_obj.update_spe_drc_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();

			rnkODEs_obj.update_reaction_rate_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.init_reaction_rate_pgt();
		}

		//reset conc_data_sr
		rnkODEs_obj.set_concentration_data_zero();
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		//on every local core
		std::size_t init_spe;
		for (size_t j = 0; j < trajectoryNumber_local; ++j) {
			//pick a initial spe randomly
			init_spe = rnkODEs_obj.return_index_randomly_given_probability_vector(prob);
			rnkODEs_obj.ODE_pathway_sim_once(init_time, end_time, init_spe);
		}//for 2
		 //get probability matrix
		rnkODEs_obj.get_probability_matrix(prob_Mat);
		//map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				// think about it later
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		broadcast(world, prob_Mat_reduce, 0);

		rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

	}

	if (world.rank() == 0) {
		rnkODEs_obj.w2f_pgt("SOHR");
	}

}

void driver::ODE_solver_MC_trajectory_cv_parallel(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;
	//initialize MPI stuff
	std::size_t iterationNumber;
	size_t trajectoryNumber_total;
	double init_time;
	double end_time;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		trajectoryNumber_total = pt.get<std::size_t>("pathway.trajectoryNumber");
		init_time = 0.0*pt.get<double>("time.tau");
		end_time = pt.get<double>("pathway.tau")*pt.get<double>("time.tau");
	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_total, 0);
	broadcast(world, init_time, 0);
	broadcast(world, end_time, 0);


	//create reaction_network
	size_t P = world.size();
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	//different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs(uncertainties, world.rank(), main_cwd);
	//probability
	std::vector<double> prob = rnkODEs.get_initial_concentration();
	double prob_total = std::accumulate(prob.begin(), prob.end(), 0.0);
	for (std::size_t i = 0; i < prob.size(); ++i)
		prob[i] /= prob_total;

	std::pair<std::size_t, std::size_t> conc_data_size = rnkODEs.get_size_of_concentration_data();

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(conc_data_size.first, std::vector<double>(conc_data_size.second, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(conc_data_size.first,
		std::vector<double>(conc_data_size.second, 0.0));

	for (std::size_t ni = 0; ni < iterationNumber; ++ni) {
		//not first iteration
		if (ni != 0) {
			rnkODEs.update_temperature_pressure_based_on_spe_concentration_cv();
			rnkODEs.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();
			rnkODEs.integrate_propensity_function_pgt();
			rnkODEs.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs.init_spe_drc_int_time_pgt();

			rnkODEs.init_reaction_rate_pgt();
		}

		//reset conc_data_sr
		rnkODEs.set_concentration_data_zero();
		rnkODEs.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		//on every local core
		std::size_t trajectoryNumber_local_O2 = static_cast<std::size_t>(trajectoryNumber_local*prob[0]);
		std::size_t trajectoryNumber_local_H2 = static_cast<std::size_t>(trajectoryNumber_local*prob[2]);

		// follow O
		for (std::size_t j = 0; j < trajectoryNumber_local_O2; ++j) {
			rnkODEs.set_reaction_out_spe_info("O");
			rnkODEs.ODE_pathway_sim_once(init_time, end_time, 0);
		}

		// follow H
		for (std::size_t j = 0; j < trajectoryNumber_local_H2; ++j) {
			rnkODEs.set_reaction_out_spe_info("H");
			rnkODEs.ODE_pathway_sim_once(init_time, end_time, 2);
		}

		//get probability matrix
		rnkODEs.get_probability_matrix(prob_Mat);
		//map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				// think about it later
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		broadcast(world, prob_Mat_reduce, 0);

		rnkODEs.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
		rnkODEs.set_concentration_at_time_zero_to_initial_fraction_or_concentration();
		rnkODEs.check_zero_concentration(pt.get<double>("SOHR_init.deltaConcentration"));
	}

	if (world.rank() == 0) {
		//rnkODEs.convert_molar_concentration_to_mole_fraction();
		rnkODEs.w2f_pgt("SOHR");
	}
}

void driver::ODE_solver_path_integral_parallel_s_ct_np_v1(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;
	//initialize MPI stuff

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;

	std::size_t iterationNumber;
	std::size_t topN;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		topN = topN_vec.front();

		//get the pathway name only on the processor 0
		pathwayHandler::get_pathway(main_cwd + std::string("/output/pathway_name_candidate.csv"), pathway_vec,
			std::numeric_limits<int>::max() - 1000); //topN pathways
		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");

	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, pathway_vec, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, topN, 0);

	////create rxn_network, generate pathway
	size_t trajectoryNumber_total = pt.get<size_t>("pathway.trajectoryNumber");
	size_t P = world.size();
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	//different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::pair<std::size_t, std::size_t> conc_data_size = rnkODEs_obj.get_size_of_concentration_data();

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(conc_data_size.first, std::vector<double>(conc_data_size.second, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(conc_data_size.first,
		std::vector<double>(conc_data_size.second, 0.0));

	//conversion factor from pathway prob to concentration
	std::vector<double> P2C;
	for (auto key1 : pt.get_child("SOHR_init.P2C")) {
		P2C.push_back(key1.second.get_value<double>());
	}
	assert(P2C.size() == conc_data_size.first);

	rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();
	for (std::size_t ni = 0; ni < iterationNumber; ++ni) {
		rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability(pathway_vec, P2C, trajectoryNumber_local, topN, prob_Mat);
		//map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}
		broadcast(world, prob_Mat_reduce, 0);

		rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);

		if (ni != 0) {
			rnkODEs_obj.update_spe_drc_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();

			rnkODEs_obj.update_reaction_rate_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.init_reaction_rate_pgt();
		}
	}

	if (world.rank() == 0) {
		rnkODEs_obj.w2f_pgt("SOHR_0");
	}
}

void driver::ODE_solver_path_integral_parallel_s_ct_np_v2(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;
	//initialize MPI stuff

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;

	std::size_t iterationNumber;
	std::size_t topN;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		topN = topN_vec.front();

		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");

	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, topN, 0);

	////create rxn_network, generate pathway
	size_t trajectoryNumber_total = pt.get<size_t>("pathway.trajectoryNumber");
	size_t P = world.size();
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	//different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	auto conc_data_size = rnkODEs_obj.get_size_of_concentration_data();

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(conc_data_size.first, std::vector<double>(conc_data_size.second, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(conc_data_size.first,
		std::vector<double>(conc_data_size.second, 0.0));

	//conversion factor from pathway prob to concentration
	std::vector<double> P2C;
	for (auto key1 : pt.get_child("SOHR_init.P2C")) {
		P2C.push_back(key1.second.get_value<double>());
	}
	assert(P2C.size() == conc_data_size.first);

	for (std::size_t ni = 0; ni < iterationNumber; ++ni) {
		if (world.rank() == 0) {
			// generate pathname
			std::string filename = main_cwd + std::string("/output") + std::string("/pathname_")
				+ boost::lexical_cast<std::string>(ni) + std::string(".csv");
			rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_length_all_elements(ni, filename);
			//get the pathway name only on the processor 0
			pathwayHandler::get_pathway(filename, pathway_vec, topN); //topN pathways
		}

		broadcast(world, pathway_vec, 0);

		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();
		rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability(pathway_vec, P2C, trajectoryNumber_local, topN, prob_Mat);
		//map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}
		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}
		broadcast(world, prob_Mat_reduce, 0);

		rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);

		if (ni != 0) {
			rnkODEs_obj.update_spe_drc_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();

			rnkODEs_obj.update_reaction_rate_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.init_reaction_rate_pgt();
		}
	}

	if (world.rank() == 0) {
		rnkODEs_obj.w2f_pgt("SOHR_0");
	}
}

void driver::ODE_solver_path_integral_parallel_s_ct_np_v3(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;
	//initialize MPI stuff

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;

	std::size_t iterationNumber;
	std::size_t topN;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		topN = topN_vec.front();

		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");

	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, topN, 0);

	////create rxn_network, generate pathway
	size_t trajectoryNumber_total = pt.get<size_t>("pathway.trajectoryNumber");
	size_t P = world.size();
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	//different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	auto conc_data_size = rnkODEs_obj.get_size_of_concentration_data();

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(conc_data_size.first, std::vector<double>(conc_data_size.second, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(conc_data_size.first,
		std::vector<double>(conc_data_size.second, 0.0));

	//conversion factor from pathway prob to concentration
	std::vector<double> P2C;
	for (auto key1 : pt.get_child("SOHR_init.P2C")) {
		P2C.push_back(key1.second.get_value<double>());
	}
	assert(P2C.size() == conc_data_size.first);

	for (std::size_t ni = 0; ni < iterationNumber; ++ni) {
		if (ni != 0) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);

			rnkODEs_obj.update_spe_drc_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();

			rnkODEs_obj.update_reaction_rate_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.init_reaction_rate_pgt();
		}

		if (world.rank() == 0) {
			// generate pathname
			std::string filename = main_cwd + std::string("/output") + std::string("/pathname_")
				+ boost::lexical_cast<std::string>(ni) + std::string(".csv");
			rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_length_all_elements(ni, filename);
			//get the pathway name only on the processor 0
			pathwayHandler::get_pathway(filename, pathway_vec, topN); //topN pathways
		}

		broadcast(world, pathway_vec, 0);

		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();
		rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability(pathway_vec, P2C, trajectoryNumber_local, topN, prob_Mat);
		//map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}
		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}

		if (ni != iterationNumber - 1)
			broadcast(world, prob_Mat_reduce, 0);

		if (ni == iterationNumber - 1) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);

			rnkODEs_obj.update_spe_drc_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();

			rnkODEs_obj.update_reaction_rate_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.init_reaction_rate_pgt();
		}
	}

	if (world.rank() == 0) {
		rnkODEs_obj.w2f_pgt("SOHR_0");
	}
}

void driver::ODE_solver_path_integral_parallel_s_ct_np_cc1_v1(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	std::vector<double> uncertainties;
	//initialize MPI stuff

	//pathway name-pathway we are interested
	std::vector<std::string> pathway_vec;
	//pathway starting with single source
	std::vector<std::string> pathway_swss_vec;

	std::size_t iterationNumber;
	std::size_t topN;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		topN = topN_vec.front();

		//get the pathway name only on the processor 0
		pathwayHandler::get_pathway(main_cwd + std::string("/output/pathway_name_candidate.csv"), pathway_vec,
			std::numeric_limits<int>::max() - 1000); //topN pathways
													 //single source species
		if (pt.get<int>("SOHR_init.single_source_species") >= 0) {
			std::string sss = pt.get<std::string>("pathway.atom_followed") + std::string("S") + pt.get<std::string>("SOHR_init.single_source_species");
			pathwayHandler::pathway_starts_with(sss, pathway_vec, pathway_swss_vec);
			//Not the source species itself
			//You can not do that directly. You need to use std::remove algorithm to move the element to be erased
			//to the end of the vector and then use erase function.

			pathway_swss_vec.erase(std::remove(pathway_swss_vec.begin(), pathway_swss_vec.end(), sss));
		}

		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");

	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, pathway_vec, 0);
	broadcast(world, pathway_swss_vec, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, topN, 0);

	////create rxn_network, generate pathway
	size_t trajectoryNumber_total = pt.get<size_t>("pathway.trajectoryNumber");
	size_t P = world.size();
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, P);

	//different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::pair<std::size_t, std::size_t> conc_data_size = rnkODEs_obj.get_size_of_concentration_data();

	//concentration oriented prob matrix result
	std::vector<std::vector<double> > conc_prob_Mat(conc_data_size.first, std::vector<double>(conc_data_size.second, 0.0));
	std::vector<std::vector<double> > conc_prob_Mat_reduce(conc_data_size.first,
		std::vector<double>(conc_data_size.second, 0.0));

	//pathway based prob matrix result
	std::vector<std::vector<double> > path_prob_Mat(pathway_swss_vec.size(), std::vector<double>(conc_data_size.second, 0.0));

	//conversion factor from pathway prob to concentration
	std::vector<double> P2C;
	for (auto key1 : pt.get_child("SOHR_init.P2C")) {
		P2C.push_back(key1.second.get_value<double>());
	}
	assert(P2C.size() == conc_data_size.first);

	rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();
	for (std::size_t ni = 0; ni < iterationNumber; ++ni) {
		rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability(pathway_vec, P2C, trajectoryNumber_local, topN, conc_prob_Mat);

		//single source term, additional add-ons to the concentration from initial species
		rnkODEs_obj.single_source_ODEdirectlyEvaluatePathwayProbability_pathVector(pathway_swss_vec, trajectoryNumber_local, path_prob_Mat);
		rnkODEs_obj.update_concentration_oriented_prob_matrix_from_single_source_path_based_prob_matrix(pathway_swss_vec, P2C, path_prob_Mat, conc_prob_Mat);


		//map reduce
		for (size_t i = 0; i < conc_prob_Mat.size(); ++i) {
			for (size_t j = 0; j < conc_prob_Mat[0].size(); ++j) {
				reduce(world, conc_prob_Mat[i][j], conc_prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(conc_prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(conc_prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(conc_prob_Mat_reduce);
			// single source concentration to be constant
			rnkODEs_obj.set_probability_matrix_of_single_source_constant(conc_prob_Mat_reduce);

		}

		broadcast(world, conc_prob_Mat_reduce, 0);

		rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(conc_prob_Mat_reduce);

		if (ni != iterationNumber - 1) {
			rnkODEs_obj.update_spe_drc_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();

			rnkODEs_obj.update_reaction_rate_based_on_spe_concentration_s_ct_np();
			rnkODEs_obj.init_reaction_rate_pgt();

			//single source additional concentration data and pointer
			rnkODEs_obj.update_single_source_additional_concentration_data_rnos_s_ct_np();
			rnkODEs_obj.init_time_single_source_additional_concentration_pointer_rnos();
		}

		//write out every iteration
		if (world.rank() == 0) {
			//rnkODEs_obj.w2f_pgt("SOHR_0");
			rnkODEs_obj.w2f_pgt(std::string("SOHR_M_") + boost::lexical_cast<std::string>(ni + 1));
		}

	}//ni
}


void driver::ODE_solver_path_integral_parallel_cv_v9(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	//initialize MPI stuff
	std::vector<double> uncertainties;

	std::size_t iterationNumber;
	std::size_t topN;

	//different trajectory nubmers for different iterations
	std::vector<size_t> trajectoryNumber_vector;

	//conversion factor from pathway prob to concentration
	std::vector<double> P2C;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		topN = topN_vec.front();

		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		for (auto key1 : pt.get_child("pathway.trajectoryNumberVector")) {
			trajectoryNumber_vector.push_back(key1.second.get_value<size_t>());
		}

		for (auto key1 : pt.get_child("SOHR_init.P2C")) {
			P2C.push_back(key1.second.get_value<double>());
		}

	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_vector, 0);
	broadcast(world, P2C, 0);
	broadcast(world, topN, 0);

	//create reaction network, generate pathway, different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::pair<std::size_t, std::size_t> species_time_size = rnkODEs_obj.get_size_of_concentration_data();
	assert(P2C.size() == species_time_size.first);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(species_time_size.first, std::vector<double>(species_time_size.second, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(species_time_size.first, std::vector<double>(species_time_size.second, 0.0));

	for (std::size_t ni = 1; ni <= iterationNumber; ++ni) {
		if (ni != 1) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);

			rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

			rnkODEs_obj.update_temperature_pressure_based_on_spe_concentration_cv();
			rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();

			rnkODEs_obj.init_reaction_rate_pgt();
		}

		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;
		//trajectory numbers
		size_t trajectoryNumber_total;
		//number of followed atoms
		std::size_t N_followed_atoms;

		//got to generate paths after resetting the concentration, pseudo-first order rate constant and reaction rates, etc.
		if (world.rank() == 0) {
			//generate pathname, not necessay that all nodes have the same set of pathway vector, modify later
			N_followed_atoms = rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(ni - 1, pathway_vec, topN);

			if (ni <= trajectoryNumber_vector.size())
				trajectoryNumber_total = trajectoryNumber_vector[ni - 1];
			else
				trajectoryNumber_total = trajectoryNumber_vector.back();
		}
		broadcast(world, pathway_vec, 0);
		broadcast(world, N_followed_atoms, 0);
		broadcast(world, trajectoryNumber_total, 0);

		size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, world.size());

		// this step is necessary since we convert mole fraction to molar concentration previously
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		// at most number of atoms * topN
		rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability(pathway_vec, P2C, trajectoryNumber_local, N_followed_atoms*topN, prob_Mat);
		//map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			rnkODEs_obj.normalize_prob_matrix_data(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}

		if (ni != iterationNumber)
			broadcast(world, prob_Mat_reduce, 0);
		if (ni == iterationNumber) {
			//last iteration, print concentration data to file
			if (world.rank() == 0) {
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				// this is necessary since we also need to update pressre and temperature
				rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

				rnkODEs_obj.update_temperature_pressure_based_on_spe_concentration_cv();
				rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

				rnkODEs_obj.integrate_propensity_function_pgt();
				rnkODEs_obj.init_spe_drc_int_pgt();

				//have to evaluate time every hopping
				rnkODEs_obj.init_spe_drc_int_time_pgt();

				rnkODEs_obj.init_reaction_rate_pgt();

				//rnkODEs_obj.w2f_pgt(std::string("SOHR_M_") + boost::lexical_cast<std::string>(ni));
				//rnkODEs_obj.convert_molar_concentration_to_mole_fraction();

				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				rnkODEs_obj.w2f_pgt(std::string("SOHR_fraction_") + boost::lexical_cast<std::string>(ni));

			}// if
		}

	}// for ni

}

void driver::ODE_solver_path_integral_parallel_cv_v10(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	//initialize MPI stuff
	std::vector<double> uncertainties;
	std::size_t iterationNumber;
	std::size_t topN;
	//different trajectory nubmers for different iterations
	std::vector<size_t> trajectoryNumber_vector;
	//conversion factor from pathway prob to concentration
	std::vector<double> P2C;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		topN = topN_vec.front();

		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		for (auto key1 : pt.get_child("pathway.trajectoryNumberVector")) {
			trajectoryNumber_vector.push_back(key1.second.get_value<size_t>());
		}

		for (auto key1 : pt.get_child("SOHR_init.P2C")) {
			P2C.push_back(key1.second.get_value<double>());
		}

	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_vector, 0);
	broadcast(world, P2C, 0);
	broadcast(world, topN, 0);

	//create reaction network, generate pathway, different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::pair<std::size_t, std::size_t> species_time_size = rnkODEs_obj.get_size_of_concentration_data();
	assert(P2C.size() == species_time_size.first);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(species_time_size.first, std::vector<double>(species_time_size.second, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(species_time_size.first, std::vector<double>(species_time_size.second, 0.0));

	for (std::size_t ni = 1; ni <= iterationNumber; ++ni) {
		if (ni != 1) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
			rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

			rnkODEs_obj.update_temperature_pressure_based_on_spe_concentration_cv();
			rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();
			rnkODEs_obj.init_reaction_rate_pgt();
		}

		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;
		//trajectory numbers
		size_t trajectoryNumber_total;
		//number of followed atoms
		std::size_t N_followed_atoms;

		std::size_t i_stage = ni - 1;
		//got to generate paths after resetting the concentration, pseudo-first order rate constant and reaction rates, etc.
		if (world.rank() == 0) {
			//generate pathname, not necessay that all nodes have the same set of pathway vector, modify later
			N_followed_atoms = rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(i_stage, pathway_vec, topN);

			if (ni <= trajectoryNumber_vector.size())
				trajectoryNumber_total = trajectoryNumber_vector[ni - 1];
			else
				trajectoryNumber_total = trajectoryNumber_vector.back();
		}
		broadcast(world, pathway_vec, 0);
		broadcast(world, N_followed_atoms, 0);
		broadcast(world, trajectoryNumber_total, 0);

		size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, world.size());

		// this step is necessary since we convert mole fraction to molar concentration previously
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		// at most number of atoms * topN
		rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability(pathway_vec, P2C, trajectoryNumber_local, N_followed_atoms*topN, prob_Mat);
		//map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			rnkODEs_obj.normalize_prob_matrix_data(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}

		if (ni != iterationNumber)
			broadcast(world, prob_Mat_reduce, 0);
		if (ni == iterationNumber) {
			//last iteration, print concentration data to file
			if (world.rank() == 0) {
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				// this is necessary since we also need to update pressre and temperature
				rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

				rnkODEs_obj.update_temperature_pressure_based_on_spe_concentration_cv();
				rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

				rnkODEs_obj.integrate_propensity_function_pgt();
				rnkODEs_obj.init_spe_drc_int_pgt();

				//have to evaluate time every hopping
				rnkODEs_obj.init_spe_drc_int_time_pgt();
				rnkODEs_obj.init_reaction_rate_pgt();

				//rnkODEs_obj.w2f_pgt(std::string("SOHR_M_") + boost::lexical_cast<std::string>(ni));
				//rnkODEs_obj.convert_molar_concentration_to_mole_fraction();
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				rnkODEs_obj.w2f_pgt(std::string("SOHR_fraction_") + boost::lexical_cast<std::string>(ni));

			}// if
		}

	}// for ni
}

void driver::ODE_solver_path_integral_parallel_cv_v11(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	//initialize MPI stuff
	std::vector<double> uncertainties;
	std::size_t iterationNumber;
	//different topN for different species
	std::vector<std::size_t> topN_vec;
	//different trajectory nubmers for different iterations
	std::vector<size_t> trajectoryNumber_vector;
	//conversion factor from pathway prob to concentration, different factor for different species
	std::vector<double> P2C;
	//path length number for all stages
	std::vector<std::size_t> path_n_v;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		for (auto key1 : pt.get_child("pathway.trajectoryNumberVector")) {
			trajectoryNumber_vector.push_back(key1.second.get_value<size_t>());
		}
		for (auto key1 : pt.get_child("SOHR_init.P2C")) {
			P2C.push_back(key1.second.get_value<double>());
		}
		for (auto key : pt.get_child("pathway.max_path_length"))
			path_n_v.push_back(key.second.get_value<std::size_t>());
	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_vector, 0);
	broadcast(world, P2C, 0);
	broadcast(world, topN_vec, 0);
	broadcast(world, path_n_v, 0);

	//create reaction network, generate pathway, different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::size_t Nspecies, Ntimepoints;
	std::tie(Nspecies, Ntimepoints) = rnkODEs_obj.get_size_of_concentration_data();
	assert(P2C.size() == Nspecies);
	assert(topN_vec.size() == Nspecies);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(Nspecies, std::vector<double>(Ntimepoints, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(Nspecies, std::vector<double>(Ntimepoints, 0.0));

	for (std::size_t ni = 1; ni <= iterationNumber; ++ni) {
		if (ni != 1) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
			rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

			rnkODEs_obj.update_temperature_pressure_based_on_spe_concentration_cv();
			rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();
			rnkODEs_obj.init_reaction_rate_pgt();
		}

		//trajectory numbers
		std::size_t trajectoryNumber_total;
		//got to generate paths after resetting the concentration, pseudo-first order rate constant and reaction rates, etc.
		if (world.rank() == 0) {
			if (ni <= trajectoryNumber_vector.size())
				trajectoryNumber_total = trajectoryNumber_vector[ni - 1];
			else
				trajectoryNumber_total = trajectoryNumber_vector.back();
		}
		broadcast(world, trajectoryNumber_total, 0);

		size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, world.size());
		// this step is necessary since we convert mole fraction to molar concentration previously
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		std::vector<std::size_t> i_stage_vec;
		//generate pathway first, different pathway for different time points
		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;
		std::size_t N_followed_atoms;
		//not the first time point
		for (std::size_t tj = 1; tj < Ntimepoints; ++tj) {
			//total stage number is 10-ish
			std::size_t i_stage = (std::size_t)get_stage_number(tj, Ntimepoints, path_n_v.size());
			// a new stage, not exist
			if (!std::binary_search(i_stage_vec.begin(), i_stage_vec.end(), i_stage)) {
				i_stage_vec.push_back(i_stage);
				if (world.rank() == 0) {
					double end_time_ratio_tmp = (i_stage + 1)*1.0 / path_n_v.size();
					double end_time_ratio = (end_time_ratio_tmp > 1.0 ? 1.0 : end_time_ratio_tmp);
					N_followed_atoms = rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(i_stage, pathway_vec, *std::max_element(topN_vec.begin(), topN_vec.end()), end_time_ratio);
				}
				broadcast(world, pathway_vec, 0);
				broadcast(world, N_followed_atoms, 0);

			}

			for (std::size_t si = 0; si < Nspecies; ++si) {
				//for every species, got to evaluate the concentration at every time point
				std::string S = std::string("S") + boost::lexical_cast<std::string>(si);
				std::vector<std::string> pathway_vec_t;
				//topN ends with S
				pathwayHandler::pathway_ends_with(S, pathway_vec, pathway_vec_t, N_followed_atoms*topN_vec[si]); //topN pathways
				rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability_si_tj(si, tj, pathway_vec_t, P2C[si], trajectoryNumber_local, prob_Mat);

			}//Nspecies
		}//Ntimepoints

		 //map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			rnkODEs_obj.normalize_prob_matrix_data(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}

		if (ni != iterationNumber)
			broadcast(world, prob_Mat_reduce, 0);
		if (ni == iterationNumber) {
			//last iteration, print concentration data to file
			if (world.rank() == 0) {
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				// this is necessary since we also need to update pressre and temperature
				rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

				rnkODEs_obj.update_temperature_pressure_based_on_spe_concentration_cv();
				rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

				rnkODEs_obj.integrate_propensity_function_pgt();
				rnkODEs_obj.init_spe_drc_int_pgt();

				//have to evaluate time every hopping
				rnkODEs_obj.init_spe_drc_int_time_pgt();
				rnkODEs_obj.init_reaction_rate_pgt();

				//rnkODEs_obj.w2f_pgt(std::string("SOHR_M_") + boost::lexical_cast<std::string>(ni));
				//rnkODEs_obj.convert_molar_concentration_to_mole_fraction();
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				rnkODEs_obj.w2f_pgt(std::string("SOHR_fraction_") + boost::lexical_cast<std::string>(ni));

			}// if
		}

	}//Niterations

}

void driver::ODE_solver_path_integral_parallel_cv_ct_v1(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	//initialize MPI stuff
	std::vector<double> uncertainties;

	std::size_t iterationNumber;
	std::size_t topN;

	//different trajectory nubmers for different iterations
	std::vector<size_t> trajectoryNumber_vector;

	//conversion factor from pathway prob to concentration
	std::vector<double> P2C;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));

		std::vector<std::size_t> topN_vec;
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		topN = topN_vec.front();

		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		for (auto key1 : pt.get_child("pathway.trajectoryNumberVector")) {
			trajectoryNumber_vector.push_back(key1.second.get_value<size_t>());
		}

		for (auto key1 : pt.get_child("SOHR_init.P2C")) {
			P2C.push_back(key1.second.get_value<double>());
		}

	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_vector, 0);
	broadcast(world, P2C, 0);
	broadcast(world, topN, 0);

	//create reaction network, generate pathway, different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::pair<std::size_t, std::size_t> species_time_size = rnkODEs_obj.get_size_of_concentration_data();
	assert(P2C.size() == species_time_size.first);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(species_time_size.first, std::vector<double>(species_time_size.second, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(species_time_size.first, std::vector<double>(species_time_size.second, 0.0));

	for (std::size_t ni = 1; ni <= iterationNumber; ++ni) {
		if (ni != 1) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);

			rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

			rnkODEs_obj.update_pressure_based_on_spe_concentration_cv();
			rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();

			rnkODEs_obj.init_reaction_rate_pgt();
		}

		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;
		//trajectory numbers
		size_t trajectoryNumber_total;
		//number of followed atoms
		std::size_t N_followed_atoms;

		//got to generate paths after resetting the concentration, pseudo-first order rate constant and reaction rates, etc.
		if (world.rank() == 0) {
			//generate pathname, not necessay that all nodes have the same set of pathway vector, modify later
			N_followed_atoms = rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(ni - 1, pathway_vec, topN);

			if (ni <= trajectoryNumber_vector.size())
				trajectoryNumber_total = trajectoryNumber_vector[ni - 1];
			else
				trajectoryNumber_total = trajectoryNumber_vector.back();
		}
		broadcast(world, pathway_vec, 0);
		broadcast(world, N_followed_atoms, 0);
		broadcast(world, trajectoryNumber_total, 0);

		size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, world.size());

		// this step is necessary since we convert mole fraction to molar concentration previously
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		// at most number of atoms * topN
		rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability(pathway_vec, P2C, trajectoryNumber_local, N_followed_atoms*topN, prob_Mat);
		//map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			rnkODEs_obj.normalize_prob_matrix_data(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}

		if (ni != iterationNumber)
			broadcast(world, prob_Mat_reduce, 0);
		if (ni == iterationNumber) {
			//last iteration, print concentration data to file
			if (world.rank() == 0) {
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				// this is necessary since we also need to update pressre and temperature
				rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

				rnkODEs_obj.update_pressure_based_on_spe_concentration_cv();
				rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

				rnkODEs_obj.integrate_propensity_function_pgt();
				rnkODEs_obj.init_spe_drc_int_pgt();

				//have to evaluate time every hopping
				rnkODEs_obj.init_spe_drc_int_time_pgt();

				rnkODEs_obj.init_reaction_rate_pgt();

				//rnkODEs_obj.w2f_pgt(std::string("SOHR_M_") + boost::lexical_cast<std::string>(ni));
				//rnkODEs_obj.convert_molar_concentration_to_mole_fraction();

				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				rnkODEs_obj.w2f_pgt(std::string("SOHR_fraction_") + boost::lexical_cast<std::string>(ni));

			}// if
		}

	}// for ni
}

void driver::ODE_solver_path_integral_parallel_cv_ct_v2(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	//initialize MPI stuff
	std::vector<double> uncertainties;
	std::size_t iterationNumber;
	//different topN for different species
	std::vector<std::size_t> topN_vec;
	//different trajectory nubmers for different iterations
	std::vector<size_t> trajectoryNumber_vector;
	//conversion factor from pathway prob to concentration, different factor for different species
	std::vector<double> P2C;
	//path length number for all stages
	std::vector<std::size_t> path_n_v;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		for (auto key1 : pt.get_child("pathway.trajectoryNumberVector")) {
			trajectoryNumber_vector.push_back(key1.second.get_value<size_t>());
		}
		for (auto key1 : pt.get_child("SOHR_init.P2C")) {
			P2C.push_back(key1.second.get_value<double>());
		}
		for (auto key : pt.get_child("pathway.max_path_length"))
			path_n_v.push_back(key.second.get_value<std::size_t>());
	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_vector, 0);
	broadcast(world, P2C, 0);
	broadcast(world, topN_vec, 0);
	broadcast(world, path_n_v, 0);

	//create reaction network, generate pathway, different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::size_t Nspecies, Ntimepoints;
	std::tie(Nspecies, Ntimepoints) = rnkODEs_obj.get_size_of_concentration_data();
	assert(P2C.size() == Nspecies);
	assert(topN_vec.size() == Nspecies);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(Nspecies, std::vector<double>(Ntimepoints, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(Nspecies, std::vector<double>(Ntimepoints, 0.0));

	for (std::size_t ni = 1; ni <= iterationNumber; ++ni) {
		if (ni != 1) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
			rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

			rnkODEs_obj.update_pressure_based_on_spe_concentration_cv();
			rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();
			rnkODEs_obj.init_reaction_rate_pgt();
		}

		//trajectory numbers
		std::size_t trajectoryNumber_total;
		//got to generate paths after resetting the concentration, pseudo-first order rate constant and reaction rates, etc.
		if (world.rank() == 0) {
			if (ni <= trajectoryNumber_vector.size())
				trajectoryNumber_total = trajectoryNumber_vector[ni - 1];
			else
				trajectoryNumber_total = trajectoryNumber_vector.back();
		}
		broadcast(world, trajectoryNumber_total, 0);

		size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, world.size());
		// this step is necessary since we convert mole fraction to molar concentration previously
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		std::vector<std::size_t> i_stage_vec;
		//generate pathway first, different pathway for different time points
		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;
		std::size_t N_followed_atoms;
		//not the first time point
		for (std::size_t tj = 1; tj < Ntimepoints; ++tj) {
			//total stage number is 10-ish
			std::size_t i_stage = (std::size_t)get_stage_number(tj, Ntimepoints, path_n_v.size());
			// a new stage, not exist
			if (!std::binary_search(i_stage_vec.begin(), i_stage_vec.end(), i_stage)) {
				i_stage_vec.push_back(i_stage);
				if (world.rank() == 0) {
					double end_time_ratio_tmp = (i_stage + 1)*1.0 / path_n_v.size();
					double end_time_ratio = (end_time_ratio_tmp > 1.0 ? 1.0 : end_time_ratio_tmp);
					N_followed_atoms = rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(i_stage, pathway_vec, *std::max_element(topN_vec.begin(), topN_vec.end()), end_time_ratio);
				}
				broadcast(world, pathway_vec, 0);
				broadcast(world, N_followed_atoms, 0);

			}

			for (std::size_t si = 0; si < Nspecies; ++si) {
				//for every species, got to evaluate the concentration at every time point
				std::string S = std::string("S") + boost::lexical_cast<std::string>(si);
				std::vector<std::string> pathway_vec_t;
				//topN ends with S
				pathwayHandler::pathway_ends_with(S, pathway_vec, pathway_vec_t, N_followed_atoms*topN_vec[si]); //topN pathways
				rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability_si_tj(si, tj, pathway_vec_t, P2C[si], trajectoryNumber_local, prob_Mat);

			}//Nspecies
		}//Ntimepoints

		 //map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			rnkODEs_obj.normalize_prob_matrix_data(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}

		if (ni != iterationNumber)
			broadcast(world, prob_Mat_reduce, 0);
		if (ni == iterationNumber) {
			//last iteration, print concentration data to file
			if (world.rank() == 0) {
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				// this is necessary since we also need to update pressre and temperature
				rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

				rnkODEs_obj.update_pressure_based_on_spe_concentration_cv();
				rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

				rnkODEs_obj.integrate_propensity_function_pgt();
				rnkODEs_obj.init_spe_drc_int_pgt();

				//have to evaluate time every hopping
				rnkODEs_obj.init_spe_drc_int_time_pgt();
				rnkODEs_obj.init_reaction_rate_pgt();

				//rnkODEs_obj.w2f_pgt(std::string("SOHR_M_") + boost::lexical_cast<std::string>(ni));
				//rnkODEs_obj.convert_molar_concentration_to_mole_fraction();
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				rnkODEs_obj.w2f_pgt(std::string("SOHR_fraction_") + boost::lexical_cast<std::string>(ni));

			}// if
		}

	}//Niterations
}

std::vector<std::string> driver::generate_pathway_running_Monte_carlo_trajectory_s2m(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt, rnkODEs::reactionNetworkODESolver & rnkODEs_obj, double end_time_ratio)
{
	std::size_t Nspecies = rnkODEs_obj.get_num_vertices();

	// pick the top N
	std::size_t mc_topN = pt.get<std::size_t>("pathway.topN_monte_carlo_path");
	std::vector<std::string> pathway_monte_carlo_vec;

	auto element_vector = rnkODEs_obj.return_element_vecotr();

	std::size_t trajectoryNumber_total = pt.get<std::size_t>("pathway.monteCarloTrajectoryNumber");
	size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, world.size());

	for (auto e : element_vector) {
		// generate pathway by running monte carlo trajectories
		std::vector<statistics > statistics_v(Nspecies);
		rnkODEs_obj.generate_path_by_running_monte_carlo_trajectory_s2m(statistics_v, trajectoryNumber_local, e.ele_name, end_time_ratio);
		// map reduce
		std::vector<std::map<std::string, int> > result_v(Nspecies);
		for (std::size_t mi = 0; mi < Nspecies; ++mi)
			reduce(world, statistics_v[mi].get_pathway_unordered_map(), result_v[mi], merge_maps(), 0);

		// core 0 only
		if (world.rank() == 0) {
			auto spe_ind_initial_conc = rnkODEs_obj.return_species_index_and_initial_concentration();
			// end species index
			std::vector<std::multimap<double, std::string, std::greater<double> > > prob_path_map_v(Nspecies);

			for (auto ic : spe_ind_initial_conc) {
				// multiply by the initial concentration
				// endwith string
				for (auto sp : result_v[ic.first]) {
					// parse the index of last element
					auto ind = sp.first.find_last_of("S");
					auto spe_ind = boost::lexical_cast<std::size_t>(sp.first.substr(ind + 1));
					// initial concentration * path count
					prob_path_map_v[spe_ind].emplace(ic.second * sp.second, sp.first);

				} // for sp
			} // for ic

			  // save topN to vector
			for (auto ppm : prob_path_map_v) {
				// for this end species following this element
				std::size_t counter = 0;
				for (auto pp : ppm) {
					if (counter < mc_topN) {
						pathway_monte_carlo_vec.push_back(pp.second);
						++counter;
					} // if
					else
						break;

				} // for pp
			} // for ppm
		} // if
	}


	return pathway_monte_carlo_vec;
}

void driver::ODE_solver_path_integral_parallel_cv_ct_v3(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	//initialize MPI stuff
	std::vector<double> uncertainties;
	std::size_t iterationNumber;
	//different topN for different species
	std::vector<std::size_t> topN_vec;
	//different trajectory nubmers for different iterations
	std::vector<size_t> trajectoryNumber_vector;
	//conversion factor from pathway prob to concentration, different factor for different species
	std::vector<double> P2C;
	//path length number for all stages
	std::vector<std::size_t> path_n_v;

	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		for (auto key1 : pt.get_child("pathway.trajectoryNumberVector")) {
			trajectoryNumber_vector.push_back(key1.second.get_value<size_t>());
		}
		for (auto key1 : pt.get_child("SOHR_init.P2C")) {
			P2C.push_back(key1.second.get_value<double>());
		}
		for (auto key : pt.get_child("pathway.max_path_length"))
			path_n_v.push_back(key.second.get_value<std::size_t>());
	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_vector, 0);
	broadcast(world, P2C, 0);
	broadcast(world, topN_vec, 0);
	broadcast(world, path_n_v, 0);

	//create reaction network, generate pathway, different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::size_t Nspecies, Ntimepoints;
	std::tie(Nspecies, Ntimepoints) = rnkODEs_obj.get_size_of_concentration_data();
	assert(P2C.size() == Nspecies);
	assert(topN_vec.size() == Nspecies);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(Nspecies, std::vector<double>(Ntimepoints, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(Nspecies, std::vector<double>(Ntimepoints, 0.0));

	for (std::size_t ni = 1; ni <= iterationNumber; ++ni) {
		if (ni != 1) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
			rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

			rnkODEs_obj.update_pressure_based_on_spe_concentration_cv();
			rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();
			rnkODEs_obj.init_reaction_rate_pgt();
		}

		//trajectory numbers
		std::size_t trajectoryNumber_total;
		//got to generate paths after resetting the concentration, pseudo-first order rate constant and reaction rates, etc.
		if (world.rank() == 0) {
			if (ni <= trajectoryNumber_vector.size())
				trajectoryNumber_total = trajectoryNumber_vector[ni - 1];
			else
				trajectoryNumber_total = trajectoryNumber_vector.back();
		}
		broadcast(world, trajectoryNumber_total, 0);

		size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, world.size());
		// this step is necessary since we convert mole fraction to molar concentration previously
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		std::vector<std::size_t> i_stage_vec;
		//generate pathway first, different pathway for different time points
		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;
		std::vector<std::string> pathway_monte_carlo_vec;
		std::size_t N_followed_atoms;
		//not the first time point
		for (std::size_t tj = 1; tj < Ntimepoints; ++tj) {
			//total stage number is 10-ish
			std::size_t i_stage = (std::size_t)get_stage_number(tj, Ntimepoints, path_n_v.size());
			// a new stage, not exist
			if (!std::binary_search(i_stage_vec.begin(), i_stage_vec.end(), i_stage)) {
				i_stage_vec.push_back(i_stage);
				if (world.rank() == 0) {
					double end_time_ratio_tmp = (i_stage + 1)*1.0 / path_n_v.size();
					double end_time_ratio = (end_time_ratio_tmp > 1.0 ? 1.0 : end_time_ratio_tmp);
					N_followed_atoms = rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(i_stage, pathway_vec, *std::max_element(topN_vec.begin(), topN_vec.end()), end_time_ratio);
				}

				pathway_monte_carlo_vec.clear();
				pathway_monte_carlo_vec = generate_pathway_running_Monte_carlo_trajectory_s2m(world, main_cwd, pt, rnkODEs_obj);

				broadcast(world, pathway_vec, 0);
				broadcast(world, pathway_monte_carlo_vec, 0);
				broadcast(world, N_followed_atoms, 0);

			}

			for (std::size_t si = 0; si < Nspecies; ++si) {
				// for every species, got to evaluate the concentration at every time point
				std::string S = std::string("S") + boost::lexical_cast<std::string>(si);
				std::vector<std::string> pathway_vec_t; std::vector<std::string> pathway_monte_carlo_vec_t;
				// topN ends with S
				pathwayHandler::pathway_ends_with(S, pathway_vec, pathway_vec_t, N_followed_atoms*topN_vec[si]); //topN pathways
				pathwayHandler::pathway_ends_with(S, pathway_monte_carlo_vec, pathway_monte_carlo_vec_t, N_followed_atoms*pt.get<std::size_t>("pathway.topN_monte_carlo_path")); //topN pathways
																																												 // merge pathway
				pathwayHandler::merge_pathway(pathway_vec_t, pathway_monte_carlo_vec_t);

				rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability_si_tj(si, tj, pathway_vec_t, P2C[si], trajectoryNumber_local, prob_Mat);

			}// Nspecies
		}// Ntimepoints

		 // map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			rnkODEs_obj.normalize_prob_matrix_data(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}

		if (ni != iterationNumber)
			broadcast(world, prob_Mat_reduce, 0);
		if (ni == iterationNumber) {
			//last iteration, print concentration data to file
			if (world.rank() == 0) {
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				// this is necessary since we also need to update pressre and temperature
				rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

				rnkODEs_obj.update_pressure_based_on_spe_concentration_cv();
				rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

				rnkODEs_obj.integrate_propensity_function_pgt();
				rnkODEs_obj.init_spe_drc_int_pgt();

				//have to evaluate time every hopping
				rnkODEs_obj.init_spe_drc_int_time_pgt();
				rnkODEs_obj.init_reaction_rate_pgt();

				//rnkODEs_obj.w2f_pgt(std::string("SOHR_M_") + boost::lexical_cast<std::string>(ni));
				//rnkODEs_obj.convert_molar_concentration_to_mole_fraction();
				rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
				rnkODEs_obj.w2f_pgt(std::string("SOHR_fraction_") + boost::lexical_cast<std::string>(ni));

			}// if
		}

	}//Niterations

}

void driver::ODE_solver_path_integral_parallel_cv_ct_v4(const boost::mpi::communicator & world, const std::string & main_cwd, const boost::property_tree::ptree & pt)
{
	//initialize MPI stuff
	std::vector<double> uncertainties;
	std::size_t iterationNumber;
	//different topN for different species
	std::vector<std::size_t> topN_vec;
	//different trajectory nubmers for different iterations
	std::vector<size_t> trajectoryNumber_vector;
	//conversion factor from pathway prob to concentration, different factor for different species
	std::vector<double> P2C;
	//path length number for all stages
	std::vector<std::size_t> path_n_v;

	//calculate the uncertainties only on the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));
		for (auto key1 : pt.get_child("pathway.topN")) {
			topN_vec.push_back(key1.second.get_value<size_t>());
		}
		iterationNumber = pt.get<std::size_t>("SOHR_init.iterationNumber");
		for (auto key1 : pt.get_child("pathway.trajectoryNumberVector")) {
			trajectoryNumber_vector.push_back(key1.second.get_value<size_t>());
		}
		for (auto key1 : pt.get_child("SOHR_init.P2C")) {
			P2C.push_back(key1.second.get_value<double>());
		}
		for (auto key : pt.get_child("pathway.max_path_length"))
			path_n_v.push_back(key.second.get_value<std::size_t>());
	}

	//broadcast
	broadcast(world, uncertainties, 0);
	broadcast(world, iterationNumber, 0);
	broadcast(world, trajectoryNumber_vector, 0);
	broadcast(world, P2C, 0);
	broadcast(world, topN_vec, 0);
	broadcast(world, path_n_v, 0);

	//create reaction network, generate pathway, different seed for different core/CPU
	rnkODEs::reactionNetworkODESolver rnkODEs_obj(uncertainties, world.rank(), main_cwd);
	std::size_t Nspecies, Ntimepoints;
	std::tie(Nspecies, Ntimepoints) = rnkODEs_obj.get_size_of_concentration_data();
	assert(P2C.size() == Nspecies);
	assert(topN_vec.size() == Nspecies);

	//pathway prob result
	std::vector<std::vector<double> > prob_Mat(Nspecies, std::vector<double>(Ntimepoints, 0.0));
	std::vector<std::vector<double> > prob_Mat_reduce(Nspecies, std::vector<double>(Ntimepoints, 0.0));

	//print out initial guess result
	if (world.rank() == 0) {
		rnkODEs_obj.convert_molar_concentration_to_mole_fraction();
		rnkODEs_obj.w2f_pgt(std::string("SOHR_fraction_") + boost::lexical_cast<std::string>(0));
	}

	for (std::size_t ni = 1; ni <= iterationNumber; ++ni) {
		if (ni != 1) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);

			rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

			rnkODEs_obj.update_pressure_based_on_spe_concentration_cv();
			rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();
			rnkODEs_obj.init_reaction_rate_pgt();
		}

		//trajectory numbers
		std::size_t trajectoryNumber_total;
		//got to generate paths after resetting the concentration, pseudo-first order rate constant and reaction rates, etc.
		if (world.rank() == 0) {
			if (ni <= trajectoryNumber_vector.size())
				trajectoryNumber_total = trajectoryNumber_vector[ni - 1];
			else
				trajectoryNumber_total = trajectoryNumber_vector.back();
		}
		broadcast(world, trajectoryNumber_total, 0);

		size_t trajectoryNumber_local = get_num_block_decomposition_2(world.rank(), trajectoryNumber_total, world.size());
		// this step is necessary since we convert mole fraction to molar concentration previously
		rnkODEs_obj.set_concentration_at_time_zero_to_initial_fraction_or_concentration();

		std::vector<std::size_t> i_stage_vec;
		//generate pathway first, different pathway for different time points
		//pathway name-pathway we are interested
		std::vector<std::string> pathway_vec;
		std::vector<std::string> pathway_monte_carlo_vec;
		std::size_t N_followed_atoms;
		//not the first time point
		for (std::size_t tj = 1; tj < Ntimepoints; ++tj) {
			//total stage number is 10-ish
			std::size_t i_stage = (std::size_t)get_stage_number(tj, Ntimepoints, path_n_v.size());
			// a new stage, not exist
			if (!std::binary_search(i_stage_vec.begin(), i_stage_vec.end(), i_stage)) {
				i_stage_vec.push_back(i_stage);
				double end_time_ratio_tmp = (i_stage + 1)*1.0 / path_n_v.size();
				double end_time_ratio = (end_time_ratio_tmp > 1.0 ? 1.0 : end_time_ratio_tmp);

				if (world.rank() == 0)
					N_followed_atoms = rnkODEs_obj.heuristic_path_string_vector_by_stage_number_path_prob_all_elements_s2m(i_stage, pathway_vec, *std::max_element(topN_vec.begin(), topN_vec.end()), end_time_ratio);

				pathway_monte_carlo_vec.clear();
				pathway_monte_carlo_vec = generate_pathway_running_Monte_carlo_trajectory_s2m(world, main_cwd, pt, rnkODEs_obj, end_time_ratio);

				broadcast(world, pathway_vec, 0);
				broadcast(world, pathway_monte_carlo_vec, 0);
				broadcast(world, N_followed_atoms, 0);
			}

			for (std::size_t si = 0; si < Nspecies; ++si) {
				// for every species, got to evaluate the concentration at every time point
				std::string S = std::string("S") + boost::lexical_cast<std::string>(si);
				std::vector<std::string> pathway_vec_t; std::vector<std::string> pathway_monte_carlo_vec_t;
				// topN ends with S
				pathwayHandler::pathway_ends_with(S, pathway_vec, pathway_vec_t, N_followed_atoms*topN_vec[si]); //topN pathways
				pathwayHandler::pathway_ends_with(S, pathway_monte_carlo_vec, pathway_monte_carlo_vec_t,
					N_followed_atoms*pt.get<std::size_t>("pathway.topN_monte_carlo_path")); //topN pathways
																							// merge pathway
				pathwayHandler::merge_pathway(pathway_vec_t, pathway_monte_carlo_vec_t);

				rnkODEs_obj.ODEdirectlyEvaluatePathwayProbability_si_tj(si, tj, pathway_vec_t, P2C[si], trajectoryNumber_local, prob_Mat);

			}// Nspecies
		}// Ntimepoints

		 // map reduce
		for (size_t i = 0; i < prob_Mat.size(); ++i) {
			for (size_t j = 0; j < prob_Mat[0].size(); ++j) {
				reduce(world, prob_Mat[i][j], prob_Mat_reduce[i][j], std::plus<double>(), 0);
			}
		}

		// normalize prob_Mat_reduce so that it represents concentration
		if (world.rank() == 0) {
			rnkODEs_obj.rescale_prob_matrix_data(prob_Mat_reduce, 1.0 / trajectoryNumber_total);
			rnkODEs_obj.divide_prob_matrix_by_number_of_ways_making_species(prob_Mat_reduce);
			// normalize pathway probability so that it represents mole fraction
			//rnkODEs_obj.normalize_prob_matrix_data(prob_Mat_reduce);
			// set the probability of the first time point to initial probability
			rnkODEs_obj.set_probability_matrix_at_time_zero_to_initial_fraction_or_concentration(prob_Mat_reduce);
		}

		if (ni != iterationNumber)
			broadcast(world, prob_Mat_reduce, 0);

		//each iteration, print concentration data to file
		if (world.rank() == 0) {
			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);

			// this is necessary since we also need to update pressre and temperature
			rnkODEs_obj.convert_mole_fraction_to_molar_concentration();

			rnkODEs_obj.update_pressure_based_on_spe_concentration_cv();
			rnkODEs_obj.update_spe_drc_reaction_rate_based_on_spe_concentration_cv();

			rnkODEs_obj.integrate_propensity_function_pgt();
			rnkODEs_obj.init_spe_drc_int_pgt();

			//have to evaluate time every hopping
			rnkODEs_obj.init_spe_drc_int_time_pgt();
			rnkODEs_obj.init_reaction_rate_pgt();

			rnkODEs_obj.update_concentration_data_from_spe_based_probability_matrix(prob_Mat_reduce);
			rnkODEs_obj.w2f_pgt(std::string("SOHR_fraction_") + boost::lexical_cast<std::string>(ni));

		}// if

	}//Niterations

}


void driver::M_matrix_R_matrix(const boost::mpi::communicator & world, const std::string & main_cwd)
{
	std::vector<double> uncertainties;
	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));
		//fileIO::fileIO::read_generate_uncertainties_w2f_random(uncertainties, main_cwd + std::string("/input/uncertainties.inp"));
		//fileIO::fileIO::read_uncertainties_from_file(uncertainties, main_cwd + std::string("/input/uncertainties.inp"));
	}
	//boradcast
	broadcast(world, uncertainties, 0);

	if (world.rank() == 0) {
		rnk::concreteReactionNetwork rnk_concrete(uncertainties, world.rank(), main_cwd);
		//rnk_concrete.initiate_M_matrix("H");
		//rnk_concrete.initiate_M_matrix();
		//rnk_concrete.print_M_matrix("H");
		//rnk_concrete.print_M_matrix("O");

		//rnk_concrete.initiate_R_matrix("H");
		//rnk_concrete.initiate_R_matrix();
		//rnk_concrete.print_R_matrix("H");
		//rnk_concrete.print_R_matrix("O");

		//std::cout << "(1,5)\t" << rnk_concrete.get_M_matrix_element("H", 1, 5) << std::endl;

		//matrix_sr::path_R_matrix_element_t p = rnk_concrete.get_R_matrix_element("O", 0, 7);

		rnk_concrete.print();

		//for (auto x : p) {
		//	// convert to path string
		//	std::cout << rnk_concrete.R_matrix_path_representation_to_string(x);
		//	//for (auto y : x) {
		//	//	std::cout << y << "-->";
		//	//}
		//	std::cout << std::endl;
		//}

		//std::vector<std::string> vs = rnk_concrete.get_path_string_i_j_with_length_n("H", 1, 5, 2);
		//std::vector<std::string> vs = rnk_concrete.get_non_zero_path_string_i_j_with_length_n("H", 1, 5, 1);
		//rnk_concrete.path_string_vector_s2f(vs);

		//std::string filename_H = main_cwd + std::string("/output") + std::string("/heuristic_pathname_H.csv");
		//std::cout << filename_H << std::endl;
		//rnk_concrete.heuristic_path_string_vector_s2f("H", 3, filename_H);

		//std::string filename_O = main_cwd + std::string("/output") + std::string("/heuristic_pathname_O.csv");
		//std::cout << filename_O << std::endl;
		//rnk_concrete.heuristic_path_string_vector_s2f("O", 3, filename_O);

		//std::string filename = main_cwd + std::string("/output") + std::string("/pathname_")
		//	+ boost::lexical_cast<std::string>(0) + std::string(".csv");
		//std::cout << filename << std::endl;
		//rnk_concrete.heuristic_path_string_vector_by_stage_number_path_length_all_elements(0, filename);

		std::cout << "Test.\n";
	}

}

void driver::MISC(const boost::mpi::communicator & world, const std::string & main_cwd)
{
	std::vector<double> uncertainties;
	//calculate the uncertainties only in the first node
	if (world.rank() == 0) {
		//read in uncertainties
		fileIO::fileIO::read_generate_uncertainties_w2f_nominal(uncertainties, (main_cwd + "/input/uncertainties.inp"));
		//fileIO::fileIO::read_generate_uncertainties_w2f_random(uncertainties, main_cwd + std::string("/input/uncertainties.inp"));
		//fileIO::fileIO::read_uncertainties_from_file(uncertainties, main_cwd + std::string("/input/uncertainties.inp"));
	}
	//boradcast
	broadcast(world, uncertainties, 0);

	if (world.rank() == 0) {
		rnk::concreteReactionNetwork rnk_concrete(uncertainties, world.rank(), main_cwd);
		//rnk_concrete.print();
		double target_time_db = rnk_concrete.return_temperature_target_time();
		std::cout << std::setprecision(15) << "time at target temperature is:\t" << target_time_db << std::endl;

		//std::vector<std::vector<double> > transition_mat = { {0.0, 1.0}, {2.0, 0.0} };
		//double first_real_positive_eigenvalue;
		//std::vector<double> equil_ratio;
		//auto result = matrix_sr::cal_equilibrium_ratio_from_transition_matrix(transition_mat, first_real_positive_eigenvalue, equil_ratio);
		//std::cout << result << std::endl;

		std::cout << "MISC\n";
	}
}


#endif // __USE_MPI_


#endif
