//#pragma once
#ifndef __UNIONFIND_CPP_
#define __UNIONFIND_CPP_

#include "../../../include/tools/union_find/unionFind.h"

UnionFind::UnionFind(int n)
{
	//Weighted Quick Union With Path Compression
	this->id.assign(n, 0);
	for (int i = 0; i < n; ++i)
		this->id[i] = i;
	this->w.assign(n, -1);
}

int UnionFind::root(int i)
{
	//path compression, flat tree, avoid deep tree
	while (i != id[i]) {
		id[i] = id[id[i]];
		i = id[i];
	}
	return i;
}

bool UnionFind::find(int p, int q)
{
	return root(p) == root(q);
}

void UnionFind::unite(int p, int q)
{
	auto i = root(p);
	auto j = root(q);

	if (w[i] < w[j]) {
		id[i] = j;
		w[j] += w[i];
	}
	else {
		id[j] = i;
		w[i] += w[j];
	}
}


#endif


