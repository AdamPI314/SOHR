#ifndef __BLOCK_DECOMPOSITION_H_
#define __BLOCK_DECOMPOSITION_H_

#include <cmath>

//ith process
//amount of task N
//number of processor p
//method 1 
int get_first_block_decomposition_1(int i, int N, int p);
int get_last_block_decomposition_1(int i, int N, int p);
int get_num_block_decomposition_1(int i, int N, int p);
//method 2
int get_first_block_decomposition_2(int i, int N, int p);
int get_last_block_decomposition_2(int i, int N, int p);
int get_num_block_decomposition_2(int i, int N, int p);

//get which stage it will be
//amount of task N
//current index of task i
//total number of stages p
int get_stage_number(int i, int N, int p);



#endif
