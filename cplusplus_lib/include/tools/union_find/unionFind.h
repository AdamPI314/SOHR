//#pragma once
#ifndef __UNIONFIND_H_
#define __UNIONFIND_H_
#include <vector>

class UnionFind {
	//Weighted Quick Union With Path Compression
public:
	UnionFind(int n);
public:
	//id
	std::vector<int> id;
	//weight
	std::vector<int> w;

public:
	int root(int i);
	bool find(int p, int q);
	void unite(int p, int q);

};

#endif