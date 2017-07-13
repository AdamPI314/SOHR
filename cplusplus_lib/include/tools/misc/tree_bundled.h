///////////////////////////////////////////////////////////////////////////////////////////////////
//This is a boost graph template with bundled property
///////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************/
//http://stackoverflow.com/questions/671714/modifying-vertex-properties-in-a-boostgraph
/*   Examples:
struct VertexProperties_graph {
int i;
};

struct EdgeProperties_graph {
};

typedef Graph_bundled<VertexProperties_graph, EdgeProperties_graph> MyGraph;

MyGraph g;

VertexProperties_graph vp;
vp.i = 42;

MyGraph::Vertex v = g.AddVertex(vp);

g.properties(v).i = 23;*/
/************************************************************************/

#ifndef __TREE_BUNDLED_H_
#define __TREE_BUNDLED_H_

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/static_assert.hpp>
#include <algorithm> //for std::find_if
#include "graph_bundled.h"

using namespace boost;

/* the graph base class template */
template < typename VERTEXPROPERTIES, typename EDGEPROPERTIES >
class Tree_bundled
{
public:
	/* an adjacency_list like we need it */
	typedef adjacency_list<
		//vecS, // allow parallel edges
		listS, // doesn't allow parallel edges
		vecS, // vertex container
		bidirectionalS, // directed graph
		/*property<vertex_properties_t, VERTEXPROPERTIES>,*/
		//property<edge_properties_t, EDGEPROPERTIES>

		property<vertex_index_t, std::size_t, property<vertex_properties_t, VERTEXPROPERTIES> >,
		property<edge_index_t, std::size_t, property<edge_properties_t, EDGEPROPERTIES> > //multiple properties
		//edge_index_t is not installed for default, if you want to use get() method to get edge_index, you have to install it explicitly
	> TreeContainer;

	/* a bunch of graph-specific typedefs */
	typedef typename graph_traits<TreeContainer>::vertex_descriptor Vertex;
	typedef std::pair<Vertex, Vertex> VertexPair;
	typedef typename graph_traits<TreeContainer>::edge_descriptor Edge;
	typedef std::pair<Edge, Edge> EdgePair;

	typedef typename graph_traits<TreeContainer>::vertex_iterator vertex_iter;
	typedef typename graph_traits<TreeContainer>::edge_iterator edge_iter;
	typedef typename graph_traits<TreeContainer>::adjacency_iterator adjacency_iter;
	typedef typename graph_traits<TreeContainer>::out_edge_iterator out_edge_iter;

	typedef typename graph_traits<TreeContainer>::degree_size_type degree_t;

	typedef std::pair<adjacency_iter, adjacency_iter> adjacency_vertex_range_t;
	typedef std::pair<out_edge_iter, out_edge_iter> out_edge_range_t;
	typedef std::pair<vertex_iter, vertex_iter> vertex_range_t;
	typedef std::pair<edge_iter, edge_iter> edge_range_t;

	//typedef typename property_map<TreeContainer, edge_index_t>::type edge_index_map_t;

	//typedef for reversed graph
	typedef reverse_graph<TreeContainer> reverse_graph_t;
	typedef typename graph_traits<reverse_graph_t>::vertex_descriptor reverse_vertex_t;
	typedef typename graph_traits<reverse_graph_t>::edge_descriptor reverse_edge_t;

	/* constructors etc. */
	Tree_bundled() {}

	Tree_bundled(const Tree_bundled& t) : tree(t.tree) {}

	virtual ~Tree_bundled() {}

	/* structure modification methods */
	void Clear()
	{
		tree.clear();
	}

	Vertex AddVertex(const VERTEXPROPERTIES& prop)
	{
		Vertex v = add_vertex(tree);
		properties(v) = prop;
		return v;
	}

	void RemoveVertex(const Vertex& v)
	{
		clear_vertex(v, tree);
		remove_vertex(v, tree);
	}

	Edge AddEdge(const Vertex& v1, const Vertex& v2, const EDGEPROPERTIES& prop_12)
	{
		/* TODO: maybe one wants to check if this edge could be inserted */
		Edge addedEdge = add_edge(v1, v2, tree).first;

		properties(addedEdge) = prop_12;

		return addedEdge;
	}

	EdgePair AddEdge(const Vertex& v1, const Vertex& v2, const EDGEPROPERTIES& prop_12, const EDGEPROPERTIES& prop_21)
	{
		/* TODO: maybe one wants to check if this edge could be inserted */
		Edge addedEdge1 = add_edge(v1, v2, tree).first;
		Edge addedEdge2 = add_edge(v2, v1, tree).first;

		properties(addedEdge1) = prop_12;
		properties(addedEdge2) = prop_21;

		return EdgePair(addedEdge1, addedEdge2);
	}

	/* property access */
	VERTEXPROPERTIES& properties(const Vertex& v)
	{
		typename property_map<TreeContainer, vertex_properties_t>::type param = get(vertex_properties, tree);
		return param[v];
	}

	const VERTEXPROPERTIES& properties(const Vertex& v) const
	{
		typename property_map<TreeContainer, vertex_properties_t>::const_type param = get(vertex_properties, tree);
		return param[v];
	}

	EDGEPROPERTIES& properties(const Edge& v)
	{
		typename property_map<TreeContainer, edge_properties_t>::type param = get(edge_properties, tree);
		return param[v];
	}

	const EDGEPROPERTIES& properties(const Edge& v) const
	{
		typename property_map<TreeContainer, edge_properties_t>::const_type param = get(edge_properties, tree);
		return param[v];
	}

	/* selectors and properties */
	const TreeContainer& getTree() const
	{
		return tree;
	}

	vertex_range_t getVertices() const
	{
		return vertices(tree);
	}
	//method one
	std::size_t get_num_vertices() const {
		vertex_range_t vp;
		int count = 0;
		for (vp = vertices(tree); vp.first != vp.second; ++vp.first)
			count++;
		return count;
	}
	//method two, recommended
	std::size_t getVertexCount() const
	{
		return num_vertices(tree);
	}

	edge_range_t getEdges() const
	{
		return edges(tree);
	}

	//method one
	std::size_t get_num_edges() const {
		edge_range_t ep;
		int count = 0;
		for (ep = edges(tree); ep.first != ep.second; ++ep.first)
			count++;
		return count;
	}
	//method two, recommended
	std::size_t getEdgeCount() const
	{
		return num_edges(tree);
	}

	adjacency_vertex_range_t getAdjacentVertices(const Vertex& v) const
	{
		return adjacent_vertices(v, tree);
	}

	out_edge_range_t getOutEdges(const Vertex&v) const {
		return out_edges(v, tree);
	}

	std::size_t getVertexDegree(const Vertex& v) const
	{
		return out_degree(v, tree);
	}

	/* operators */
	Tree_bundled& operator=(const Tree_bundled &rhs)
	{
		tree = rhs.tree;
		return *this;
	}

	//find the first edge with specific edge_index
	edge_iter find_edge_with_an_edge_index(size_t x) {
		edge_iter iter_beg, iter_end;
		boost::tie(iter_beg, iter_end) = getEdges();

		//return std::find_if(iter_beg, iter_end, check_x(x));		
		for (; iter_beg != iter_end; ++iter_beg) {
			if (properties(*iter_beg).edge_index == x) {
				break;
			}
		}
		return iter_beg;
	}

protected:
	TreeContainer tree;
};

#endif
