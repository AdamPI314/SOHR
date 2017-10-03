//#pragma once
#ifndef __EPPSTEIN_ALGORITHM_H_
#define __EPPSTEIN_ALGORITHM_H_
#include <iostream>
#include "../../relationshipParser/relationshipParser.h"
#include "../../tools/misc/tree_bundled.h"
#include <queue>

namespace eppstein_algorithm
{
	namespace rsp = relationshipParser_sr;

	typedef rsp::index_int_t vertex_index_t;
	typedef rsp::index_int_t edge_index_t;
	typedef rsp::index_int_t path_index_t;
	typedef int path_length_t;

	struct VertexProperties_tree
	{
		//index of current node in new tree
		vertex_index_t vertex_index_in_tree;
		//source edge of sidetrack edge in the original graph
		vertex_index_t from_vertex;
		//target edge of sidetrack edge in the original graph
		vertex_index_t to_vertex;
		//index of sidetrack edge in the original graph
		edge_index_t reaction_index;
		//cost of the corresponding path
		double cost;

		//path length from the origin node to "from_vertex", right before "to_vertex"
		std::size_t path_length_before_to_vertex;
		VertexProperties_tree() :vertex_index_in_tree(0), from_vertex(0), to_vertex(0), reaction_index(0), cost(0.0), path_length_before_to_vertex(0) {};

	};

	struct EdgeProperties_tree
	{
		//doesn't matter actually
		//edge_t edge;
	};

	class sidetrack_tree : public Tree_bundled<VertexProperties_tree, EdgeProperties_tree>
	{
	public:
		sidetrack_tree();
		~sidetrack_tree();
	};

	//length of path, path priority struct
	struct len_path_t
	{
		//path cost
		double cost;
		//node index in sidetrack tree
		vertex_index_t vertex_index_in_sidetrack_tree;

		//parent path's or says reference path's index in PATH vector
		path_index_t ref_path_index;

		bool operator < (const len_path_t& rhs) const {
			return cost > rhs.cost;
		}
		len_path_t() :cost(0.0), vertex_index_in_sidetrack_tree(0), ref_path_index(0) {}
	};

	//path cost, reference path, length of path before sidetrack edge
	struct path_info_t
	{
		//path cost
		double cost;
		//reference path, index in vector
		path_index_t ref_path_index;
		//length of path before sidetrack edge
		path_length_t path_length_before_to_vertex;
		//target edge of sidetrack edge in the original graph
		vertex_index_t to_vertex;
		//index of sidetrack edge in the original graph, that edge ends with to_vertex
		edge_index_t reaction_index;

		path_info_t():cost(0.0), ref_path_index(0), path_length_before_to_vertex(0), to_vertex(0), reaction_index(0){}

	};

	typedef std::priority_queue<len_path_t, std::vector<len_path_t> > len_path_pq_t;

}//eppstein_algorithm


#endif
