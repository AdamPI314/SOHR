#ifndef __MATRIX_SR_H_
#define __MATRIX_SR_H_

#include<vector>

namespace matrix_sr {
	// size_t_matrix type
	typedef std::vector<std::vector<std::size_t> > size_t_matrix_t;
	size_t_matrix_t matrix_multiplication(size_t_matrix_t m_k, size_t_matrix_t k_n);
	size_t_matrix_t matrix_power(size_t_matrix_t n_n, std::size_t k);

	typedef std::size_t path_edge_t;

	// matrix element, so each matrix element is actually a vector of vector, the inside vector is used
	// to store a single path, means each matrix element is a bunch of path(s) with equal length
	typedef std::vector<path_edge_t> path_t;
	typedef std::vector<path_t > path_R_matrix_element_t;
	typedef std::vector<std::vector<path_R_matrix_element_t> > path_R_matrix_t;

	// notice that the matrix multiplication has property: Associativity, but not Commutativity
	path_R_matrix_t matrix_multiplication(path_R_matrix_t m_k, path_R_matrix_t k_n);
	path_R_matrix_t matrix_power(path_R_matrix_t n_n, std::size_t k);

	// calculate equilibrium concentration/ratios from transition matrix
	bool cal_equilibrium_ratio_from_transition_matrix(std::vector<std::vector<double> > &transition_mat, double &first_positive_eigenvalue, std::vector<double> &equil_ratio);

}


#endif //__MATRIX_SR_H_





