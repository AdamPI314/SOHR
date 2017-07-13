/*
 * boostMPI template, boost+boostMPI, include boost include/lib MS-MPI include, and MS-MPI library, which is msmpi.lib
 * Author: Shirong Bai
 * Email: bunnysirah@hotmail.com or shirong.bai@colorado.edu
 * Date: 10/08/2015
 */

#include "include/drivers/drivers.h"
#include "include/tools/my_debug/my_debug.h"

int main(int argc, char **argv) {

	boost::mpi::environment env(argc, argv);
	boost::mpi::communicator world;

	MyClock_us timer;
	if (world.rank() == 0)
		timer.begin();

	/*************************************************************************************************/
	/*
	* 0. Parse parameters
	* Default, has to run
	*/
	/*************************************************************************************************/
	po::variables_map vm;
	// current working directory
	std::string main_cwd;
	// boost property tree
	boost::property_tree::ptree pt;

	driver::parse_parameters(argc, argv, vm, main_cwd, pt);


	/*************************************************************************************************/
	/*
	* 1. Solve for concentration of Lokta-Voltera system or Dimerization or Michaelis�Menten, using LSODE
	*/
	/*************************************************************************************************/
	if (pt.get<std::string>("job.job_type") == std::string("solve_ODEs_for_concentration_using_LSODE"))
		driver::solve_ODEs_for_concentration_using_LSODE(world, main_cwd);
	else if (pt.get<std::string>("job.job_type") == std::string("solve_ODEs_for_concentration_using_SSA"))
		driver::solve_ODEs_for_concentration_using_SSA(world, main_cwd);

	/*************************************************************************************************/
	/*
	* 2. Generate pathway first, it might take some time, depends on what kinds of pathway you want
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("generate_pathway_running_Monte_carlo_trajectory"))
		driver::generate_pathway_running_Monte_carlo_trajectory(world, main_cwd, pt);

	/*************************************************************************************************/
	/*
	* 3. Evaluate pathway integral using Importance based Monte Carlo simulation
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("evaluate_path_integral_over_time"))
		driver::evaluate_path_integral_over_time(world, main_cwd, pt);

	/*************************************************************************************************/
	/*
	* 4. For speciation, print out concentration at for different sets of rate constant
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") ==
		std::string("speciation_evaluate_concentrations_for_different_sets_rate_coefficients"))
		driver::speciation_evaluate_concentrations_for_different_sets_rate_coefficients(world, main_cwd, pt);

	/*************************************************************************************************/
	/*
	* 5. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	* version 2
	* Monte Carlo simulation
	* not parallel code
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_MC_trajectory_single_core"))
		driver::ODE_solver_MC_trajectory_single_core(world, main_cwd, pt);

	/*************************************************************************************************/
	/*
	* 6. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	* version 4
	* constant temperature, surface reactions, independent of pressure
	* Monte Carlo Simulation
	* Parallel version
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_MC_trajectory_s_ct_np_parallel"))
		driver::ODE_solver_MC_trajectory_s_ct_np_parallel(world, main_cwd, pt);

	/*************************************************************************************************/
	/*
	* 7. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	* version 5
	* varing temperature, varing pressure, constant volume, H2-O2 system, not working because the initiation
	* is actually a rare event, stochastic trajectory converge too slow
	* Monte Carlo Simulation
	* Parallel version
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_MC_trajectory_cv_parallel"))
		driver::ODE_solver_MC_trajectory_cv_parallel(world, main_cwd, pt);

	/*************************************************************************************************/
	/*
	* 8. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	* version 3
	* surface reaction, independent of pressure, constant temperature, or independent or temperature
	* directly evaluate the pathway probability
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_s_ct_np_v1"))
		driver::ODE_solver_path_integral_parallel_s_ct_np_v1(world, main_cwd, pt);
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_s_ct_np_v2"))
		driver::ODE_solver_path_integral_parallel_s_ct_np_v2(world, main_cwd, pt);
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_s_ct_np_v3"))
		driver::ODE_solver_path_integral_parallel_s_ct_np_v3(world, main_cwd, pt);

	/*************************************************************************************************/
	/*
	* 9. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	* surface reaction, independent of pressure, constant temperature, or independent or temperature
	* directly evaluate the pathway probability
	* hold concentration of some species to be constant
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_s_ct_np_cc1_v1"))
		driver::ODE_solver_path_integral_parallel_s_ct_np_cc1_v1(world, main_cwd, pt);


	/*************************************************************************************************/
	/*
	* 10. Pathway based equation solver, solve differential equations derived from chemical mechanisms
	* version 3
	* constant volume, varying temperature, varying pressure
	* directly evaluate the pathway probability
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_cv_v9"))
		driver::ODE_solver_path_integral_parallel_cv_v9(world, main_cwd, pt);
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_cv_v10"))
		driver::ODE_solver_path_integral_parallel_cv_v10(world, main_cwd, pt);
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_cv_v11"))
		driver::ODE_solver_path_integral_parallel_cv_v11(world, main_cwd, pt);
	/*
	* 11. constant volume, constant temperature
	*/
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_cv_ct_v1"))
		driver::ODE_solver_path_integral_parallel_cv_ct_v1(world, main_cwd, pt);
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_cv_ct_v2"))
		driver::ODE_solver_path_integral_parallel_cv_ct_v2(world, main_cwd, pt);
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_cv_ct_v3"))
		driver::ODE_solver_path_integral_parallel_cv_ct_v3(world, main_cwd, pt);
	else if (pt.get<std::string>("job.job_type") == std::string("ODE_solver_path_integral_parallel_cv_ct_v4"))
		driver::ODE_solver_path_integral_parallel_cv_ct_v4(world, main_cwd, pt);

	/*************************************************************************************************/
	/*
	* 12. Dijkstra algorithm, finding the shortest path
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("k_shortest_path_algorithms"))
		driver::k_shortest_path_algorithms(world, main_cwd);

	/*************************************************************************************************/
	/*
	* 13. M-matrix test
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("M_matrix_R_matrix"))
		driver::M_matrix_R_matrix(world, main_cwd);
	/*************************************************************************************************/
	/*
	* Test
	*/
	/*************************************************************************************************/
	else if (pt.get<std::string>("job.job_type") == std::string("Test"))
		driver::Test(world, main_cwd);


	return EXIT_SUCCESS;
}