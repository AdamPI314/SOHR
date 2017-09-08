#ifndef __DRIVERS_H_
#define __DRIVERS_H_

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <string>
#include <exception>
#include <cstdlib>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/property_tree/json_parser.hpp> //for json_reader

#include "../time/myclock_us.h"
#include "../pathwayHandler/pathwayHandler.h"
#include "../statistics/statistics.h"
#include "../statistics/statistics_get_concentration.h"
#include "../fileIO/CSVFileReader/CSVReader.h"
#include "../fileIO/commandLineConfigFileReader/clcf_parser.h"
#include "../tools/block_decomposition/block_decomposition.h"
#include "../tools/map_reduce/map_reduce.h"
#include "../tools/misc/misc_template.h"
#include "../tools/misc/global_extern_vars.h"
#include "../tools/misc/graph_bundled.h"
#include "../fileIO/fileIO/fileIO.h"
#include "../mechanism/mechanism.h"
#include "../relationshipParser/relationshipParser.h"
#include "../reactionNetwork/reactionNetworkODESolver/reactionNetworkODESolver.h"
#include "../tools/matrix/matrix_sr.h"


#include "../srkin/srkin.h"

//#include "../propagator/dlsodePropagator/dlsodePropagator.h" /*contains boost and cubic spline*/
#include "../reactionNetwork/concreteReactionNetwork/concreteReactionNetwork.h" /*contains dlsodePropagator so that contains boost and cubic spline*/


//#include <boost/mpi/environment.hpp>
//#include <boost/mpi/communicator.hpp>

//Try including nr3.h after all boost and system headers!!!
//#include "../cubicSpline/interp_1d.h"



namespace driver {
	namespace rsp = relationshipParser_sr;
	namespace pgt = propagator_sr;
	namespace rnk = reactionNetwork_sr;
	namespace rnkODEs = reactionNetworkODESolver_sr;

	/*0. Parse parameters*/
	void parse_parameters(const int argc, char **argv, po::variables_map &vm, std::string &main_cwd, boost::property_tree::ptree &pt);

	/*1. Solve ODEs for concentrations using LSODE*/
	void solve_ODEs_for_concentration_using_LSODE(const boost::mpi::communicator &world, std::string &main_cwd, const boost::property_tree::ptree &pt);
	//Gillespie Stochastic Simulation Algorithm
	void solve_ODEs_for_concentration_using_SSA(const boost::mpi::communicator &world, std::string &main_cwd, const boost::property_tree::ptree &pt);

	/*2. Generate pathway first, it might take some time, depends on what kinds of pathway you want*/
	void generate_pathway_running_Monte_carlo_trajectory(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);

	/*3. Evaluate pathway integral using Importance based Marte Carlo simulation*/
	void evaluate_path_integral_over_time(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);

	/*4. For speciation, print out concentration at for different sets of rate constant*/
	void speciation_evaluate_concentrations_for_different_sets_rate_coefficients(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);

	/*
	 * 5. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	 * version 2
	 * Monte Carlo simulation
	 * not parallel code
	 */
	void ODE_solver_MC_trajectory_single_core(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);

	/*
	 * 6. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	 * version 4
	 * constant temperature, surface reactions, independent of pressure
	 * Monte Carlo Simulation
	 * Parallel version
	 */
	void ODE_solver_MC_trajectory_s_ct_np_parallel(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);

	/*
	 * 7. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	 * version 5
	 * varing temperature, varing pressure, constant volume, H2-O2 system, not working because the initiation
	 * is actually a rare event, stochastic trajectory converge too slow
	 * Monte Carlo Simulation
	 * Parallel version
	 */
	void ODE_solver_MC_trajectory_cv_parallel(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	/*
	 * 8. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	 * version 3
	 * surface reaction, independent of pressure, constant temperature
	 * directly evaluate the pathway probability
	 * pathname is fixed for every iteration, read pathname from file
	 */
	 /*read pathname from file*/
	void ODE_solver_path_integral_parallel_s_ct_np_v1(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	/* pathname varies with iteration, path importance depends on path length*/
	void ODE_solver_path_integral_parallel_s_ct_np_v2(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	/* pathname varies with iteration, path importance depends on heuristic pathway probabilities*/
	void ODE_solver_path_integral_parallel_s_ct_np_v3(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);

	/*read pathname from file, hold the concentration of the frist species to be constant*/
	void ODE_solver_path_integral_parallel_s_ct_np_cc1_v1(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);


	/*
	* 9. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	* constant volume
	*/
	// get rid of intermediate files, in other words, do not write intermediate files to disks anymore, like path prob matrix
	void ODE_solver_path_integral_parallel_cv_v9(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	// stage number here is the iteration number
	void ODE_solver_path_integral_parallel_cv_v10(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	// generate pathway separately every time points, stage number here is the time stage number
	void ODE_solver_path_integral_parallel_cv_v11(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	/*
	* 10. constant volume, constant temperature
	*/
	void ODE_solver_path_integral_parallel_cv_ct_v1(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	// generate pathway separately every time points, stage number here is the time stage number
	void ODE_solver_path_integral_parallel_cv_ct_v2(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	// generate monte carlo pathways, going to generate pathway on multiple cores, that's why we write this function here
	std::vector<std::string> generate_pathway_running_Monte_carlo_trajectory_s2m(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt, rnkODEs::reactionNetworkODESolver &rnkODEs_obj, double end_time_ratio = 1.0);
	// incorpate monte carlo pathways into the pathway gengerating step, besides matrix based method
	void ODE_solver_path_integral_parallel_cv_ct_v3(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);
	// for genetic toggle switch system
	void ODE_solver_path_integral_parallel_cv_ct_v4(const boost::mpi::communicator &world, const std::string &main_cwd, const boost::property_tree::ptree &pt);

	/*
	 * 12. k-shortest path algorithm
	 */
	void k_shortest_path_algorithms(const boost::mpi::communicator &world, const std::string &main_cwd);

	/*
	 * 13. M-matrix and R-matrix
	 */
	void M_matrix_R_matrix(const boost::mpi::communicator &world, const std::string &main_cwd);

	/*
	* Miscellaneous routines, For TEST purpose.
	*/
	void MISC(const boost::mpi::communicator &world, const std::string &main_cwd);
} /* namespace driver*/

#endif
