#ifndef __MATRIX_SR_CPP_
#define __MATRIX_SR_CPP_
#include "matrix_sr.h"
#include <iostream>
#include <cassert>

namespace matrix_sr {
	size_t_matrix_t matrix_multiplication(size_t_matrix_t m_k, size_t_matrix_t k_n) {
		std::size_t m = m_k.size();
		std::size_t m_k_k = m_k[0].size();
		std::size_t k_n_k = k_n.size();
		std::size_t n = k_n[0].size();
		// #columns of first matrix equals #rows of second matrix
		assert(m_k_k == k_n_k);

		size_t_matrix_t return_matrix(m, std::vector<std::size_t>(n, 0));

		for (std::size_t i = 0; i < m; ++i) {
			for (std::size_t j = 0; j < n; ++j) {
				return_matrix[i][j] = 0;
				for (std::size_t k = 0; k < m_k_k; ++k) {
					return_matrix[i][j] += m_k[i][k] * k_n[k][j];
				}
			}
		}

		return return_matrix;
	}

	size_t_matrix_t matrix_power(size_t_matrix_t n_n, std::size_t k)
	{
		if (k == 0)
			return size_t_matrix_t();
		assert(k >= 1);
		size_t_matrix_t A = n_n;
		for (std::size_t i = 0; i < k - 1; ++i) {
			A = matrix_multiplication(A, n_n);
		}

		return A;
	}

	path_R_matrix_t matrix_multiplication(path_R_matrix_t m_k, path_R_matrix_t k_n)
	{
		std::size_t m = m_k.size();
		std::size_t m_k_k = m_k[0].size();
		std::size_t k_n_k = k_n.size();
		std::size_t n = k_n[0].size();
		// #columns of first matrix equals #rows of second matrix
		assert(m_k_k == k_n_k);
		// each matrix element is a potential matrix
		path_R_matrix_t return_matrix(m, std::vector<path_R_matrix_element_t>(n));

		for (std::size_t i = 0; i < m; ++i) {
			for (std::size_t j = 0; j < n; ++j) {
				// define multiplication
				// combine two vectors into a single vector
				path_R_matrix_element_t p_new;
				for (std::size_t k = 0; k < m_k_k; ++k) {
					for (std::size_t l1 = 0; l1 < m_k[i][k].size(); ++l1) {
						// for all none-zero path, zero means no path, just a empty vector
						if (m_k[i][k][l1].size() != 0) {
							for (std::size_t l2 = 0; l2 < k_n[k][j].size(); ++l2) {
								if (k_n[k][j][l2].size() != 0) {
									std::vector<std::size_t > v_new;
									// combine two vectors
									v_new.reserve(m_k[i][k][l1].size() + k_n[k][j][l2].size()); // preallocate memory
									v_new.insert(v_new.end(), m_k[i][k][l1].begin(), m_k[i][k][l1].end());
									v_new.insert(v_new.end(), k_n[k][j][l2].begin(), k_n[k][j][l2].end());
									p_new.push_back(v_new);
								} // if
							} // l2

						} // if

					} // l1
				} // k

				if (p_new.size() == 0)
					p_new = {};
				return_matrix[i][j] = p_new;
			} // j
		} // i

		return return_matrix;
	}

	path_R_matrix_t matrix_power(path_R_matrix_t n_n, std::size_t k)
	{
		if (k == 0)
			return path_R_matrix_t();
		assert(k >= 1);
		path_R_matrix_t A = n_n;
		for (std::size_t i = 0; i < k - 1; ++i) {
			A = matrix_multiplication(A, n_n);
		}

		return A;

	}

	//matrix_sr::size_t_matrix_t m1;
	//m1.resize(5);
	//m1[0] = { 0,1,1,1,1 };
	//m1[1] = { 0,0,1,1,1 };
	//m1[2] = { 0,0,0,1,1 };
	//m1[3] = { 0,0,0,0,1 };
	//m1[4] = { 0,0,0,0,0 };

	//matrix_sr::size_t_matrix_t m2 = matrix_sr::matrix_multiplication(m1, m1);

	//for (auto x : m2)
	//	std::cout << x << std::endl;

	//matrix_sr::size_t_matrix_t A = matrix_sr::matrix_power(m1, 2);
	//for (auto x : A)
	//	std::cout << x << std::endl;

	//matrix_sr::path_R_matrix_t m1;
	//m1.resize(5);
	//// each matrix element is a matrix
	//m1[0] = { {},{ { 1 } },{ { 2 } },{ { 3 } },{ { 4 } } };
	//m1[1] = { {},{},{ { 5 } },{ { 6 } },{ { 7 } } };
	//m1[2] = { {},{},{},{ { 8 } },{ { 9 } } };
	//m1[3] = { {},{},{},{},{ {10}, {11} } };
	//m1[4] = { {},{},{},{},{} };

	//for (std::size_t i = 0; i < m1.size(); ++i) {
	//	for (std::size_t j = 0; j < m1[i].size(); ++j) {
	//		std::cout << "(" << i << "," << j << ")\t";
	//		if (m1[i][j].size() == 0) {
	//			continue;
	//		}
	//		for (std::size_t k = 0; k < m1[i][j].size(); ++k) {
	//			for (std::size_t l = 0; l < m1[i][j][k].size(); ++l)
	//				std::cout << m1[i][j][k][l] << "\t";
	//		}
	//	}
	//	std::cout << std::endl;
	//}

	//matrix_sr::path_R_matrix_t m2 = matrix_sr::matrix_multiplication(m1, m1);

	//for (std::size_t i = 0; i < m2.size(); ++i) {
	//	for (std::size_t j = 0; j < m2[i].size(); ++j) {
	//		std::cout << "(" << i << "," << j << ")\t";
	//		if (m2[i][j].size() == 0) {
	//			continue;
	//		}
	//		for (std::size_t k = 0; k < m2[i][j].size(); ++k) {
	//			for (std::size_t l = 0; l < m2[i][j][k].size(); ++l)
	//				std::cout << m2[i][j][k][l] << "\t";
	//		}
	//	}
	//	std::cout << std::endl;
	//}

	//matrix_sr::path_R_matrix_t m3 = matrix_sr::matrix_power(m1, 3);

	//for (std::size_t i = 0; i < m3.size(); ++i) {
	//	for (std::size_t j = 0; j < m3[i].size(); ++j) {
	//		std::cout << "(" << i << "," << j << ")\t";
	//		if (m3[i][j].size() == 0) {
	//			continue;
	//		}
	//		for (std::size_t k = 0; k < m3[i][j].size(); ++k) {
	//			for (std::size_t l = 0; l < m3[i][j][k].size(); ++l)
	//				std::cout << m3[i][j][k][l] << "\t";
	//		}
	//	}
	//	std::cout << std::endl;
	//}

}

#endif // !__MATRIX_SR_CPP_