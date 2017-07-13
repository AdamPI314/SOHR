#ifndef __BLOCK_DECOMPOSITION_CPP_
#define __BLOCK_DECOMPOSITION_CPP_

#include <cmath>
#include <algorithm>
#include "block_decomposition.h"

//ith process
//amount of task N
//number of processor p
int get_first_block_decomposition_1(int i, int N, int p) {
	int r = N%p;
	return i*(int)floor(N / p) + std::min(i, r);
}
int get_last_block_decomposition_1(int i, int N, int p) {
	int r = N%p;
	return (i + 1)*(int)floor(N / p) + std::min(i + 1, r) - 1;
}
int get_num_block_decomposition_1(int i, int N, int p) {
	return get_last_block_decomposition_1(i, N, p) - get_first_block_decomposition_1(i, N, p) + 1;
}

int get_first_block_decomposition_2(int i, int N, int p) {
	return (int)floor(i*N / p);
}
int get_last_block_decomposition_2(int i, int N, int p) {
	return (int)floor((i + 1)*N / p) - 1;
}
int get_num_block_decomposition_2(int i, int N, int p) {
	return get_last_block_decomposition_2(i, N, p) - get_first_block_decomposition_2(i, N, p) + 1;
}

int get_stage_number(int i, int N, int p)
{
	int s = (int)floor(i*p / N);
	s = ((s >= p) ? p : s);
	return s;
}

#endif
